
const use_std_lib_threads = @load_preference("use_std_lib_threads", false)
const use_fast_math = @load_preference("use_fast_math", true)
const use_inbounds = @load_preference("use_inbounds", true)


"""
Controls which multi-threading library to use.
"""
macro threads(expr)
    if use_std_lib_threads
        return esc(quote
            Threads.@threads :static $(expr)
        end)
    else
        return esc(quote
            @batch $(expr)
        end)
    end
end


macro fast(expr)
    expr = use_fast_math ? :(@fastmath $expr) : expr
    expr = use_inbounds  ? :(@inbounds $expr) : expr
    return esc(expr)
end


# TODO: check if all function calls (excluding Base, Core, CUDA, AMDGPU, KA...) in a kernel are for
#  functions annotated with @kernel_function
macro kernel_function(func)
    func = :(@inline $func)
    func = use_fast_math ? :(@fastmath $func) : func
    return esc(func)
end


function make_threaded_loop(expr::Expr; choice=:dynamic)
    with = :(@inbounds @threads $(expr))
    without = :($(expr))

    if choice == :dynamic
        return quote
            if params.use_threading
                $(with)
            else
                $(without)
            end
        end
    elseif choice == :with
        return with
    elseif choice == :without
        return without
    else
        error("Unknown 'choice' value: $choice")
    end
end


function make_simd_loop(expr::Expr; choice=:dynamic)
    with = :(@fastmath @inbounds @simd ivdep $(expr))
    without = :(@inbounds $(expr))

    if choice == :dynamic
        return quote
            if params.use_simd
                $(with)
            else
                $(without)
            end
        end
    elseif choice == :with
        return with
    elseif choice == :without
        return without
    else
        error("Unknown 'choice' value: $choice")
    end
end


function extract_range_expr(expr::Expr)
    if Meta.isexpr(expr, :block, 2)
        if expr.args[1] isa LineNumberNode && Meta.isexpr(expr.args[2], :for, 2)
            # Something of the form `block # line number # for i = r <body> end end`
            expr_copy = copy(expr)
            range_expr = expr_copy.args[2].args[1]
            loop_body = expr_copy.args[2].args[2]
            return expr_copy, range_expr, loop_body
        end
    elseif Meta.isexpr(expr, :for, 2)
        # Something of the form `for i = r <body> end`
        expr_copy = copy(expr)
        range_expr = expr_copy.args[1]
        loop_body = expr_copy.args[2]
        return expr_copy, range_expr, loop_body
    end
    error("Expected a valid for loop")
end


function make_simd_threaded_loop(expr::Expr; threading=:dynamic, simd=:dynamic, add_iteration_indexes=false)
    loop_expr, range_expr, loop_body = extract_range_expr(expr)

    if @capture(range_expr, i_ = r_)
        loop_range = r
        range_expr.args[2] = :(__ideb:__step:__ifin)
    else
        error("Expected a loop expression of the form `i = range_var` or `i = 1:2:10`")
    end

    if add_iteration_indexes
        i = range_expr.args[1]

        # Replace `i = range_var` by `__i_idx = 1:__loop_length`
        range_expr.args[1] = :(__i_idx)
        range_expr.args[2] = :(1:__loop_length)

        # Add in the loop body, before the first statement: `i = range_var[__i_idx] - 1 + __j_iter`
        first_expr_idx = loop_body.args[1] isa LineNumberNode ? 2 : 1
        insert!(loop_body.args, first_expr_idx, :($i = __loop_range[__i_idx] - 1 + __j_iter))
    end

    return quote
        let __loop_range = $loop_range, __loop_length = length(__loop_range),
                __total_iter = length(__loop_range), __num_threads = Threads.nthreads(),
                # Equivalent to __total_iter ÷ __num_threads
                __batch = convert(Int, cld(__total_iter, __num_threads))::Int,
                __first_i = first(__loop_range),
                __last_i = last(__loop_range),
                __step = step(__loop_range)
            $(make_threaded_loop(:(for __i_thread = 1:__num_threads
                __j_iter = __batch * (__i_thread - 1)
                __ideb = __first_i + __j_iter * __step
                __ifin = min(__ideb + (__batch - 1) * __step, __last_i)
                $(make_simd_loop(loop_expr, choice=simd))
            end), choice=threading))
        end
    end
end


function make_simd_threaded_iter(range, expr::Expr; threading=:dynamic, simd=:dynamic, add_iteration_indexes=false)
    loop_expr, range_expr, loop_body = extract_range_expr(expr)

    if @capture(range_expr, i_ = r_)
        loop_range = r
        range_expr.args[2] = :(__loop_range .+ (__j - 1))
    else
        error("Expected a loop expression of the form `i = range_var` or `i = 1:2:10`")
    end

    if add_iteration_indexes
        i = range_expr.args[1]

        # Replace `i = range_var` by `__i_idx = 1:__loop_length)`
        range_expr.args[1] = :(__i_idx)
        range_expr.args[2] = :(1:__loop_length)

        # Add in the loop body, before the first statement: `i = range_var[__i_idx] - 1 + __j`
        first_expr_idx = loop_body.args[1] isa LineNumberNode ? 2 : 1
        insert!(loop_body.args, first_expr_idx, :($i = __loop_range[__i_idx] - 1 + __j))

        # By adding `__j_idx` and `__j_iter`, we can know where we are in the total iterations with
        # `__j_iter + __i_idx`
        main_loop_idx = :(__j_idx)
        main_loop_expr = :(Base.OneTo(length(__main_range)))
        init_inner_loop = quote
            __j = first(__main_range) + step(__main_range) * (__j_idx - 1)
            __j_iter = (__j_idx - 1) * length(__loop_range)
        end
    else
        main_loop_idx = :(__j)
        main_loop_expr = :(__main_range)
        init_inner_loop = Expr(:block)
    end

    return quote
        let __main_range = $range, __loop_range = $loop_range, __loop_length = length(__loop_range)
            $(make_threaded_loop(:(for $main_loop_idx in $main_loop_expr
                $init_inner_loop
                $(make_simd_loop(loop_expr, choice=simd))
            end), choice=threading))
        end
    end
end


"""
    @simd_loop(expr)

Allows to enable/disable SIMD optimisations for a loop.
When SIMD is enabled, it is assumed that there is no dependencies between each iterations of the loop.

```julia
    @simd_loop for i = 1:n
        y[i] = x[i] * (x[i-1])
    end
```
"""
macro simd_loop(expr)
    return esc(make_simd_loop(expr))
end


"""
    @threaded(expr)

Allows to enable/disable multithreading of the loop depending on the parameters.

```julia
    @threaded for i = 1:n
        y[i] = log10(x[i]) + x[i]
    end
```
"""
macro threaded(expr)
    return esc(make_threaded_loop(expr))
end


"""
    @simd_threaded_loop(expr)

Allows to enable/disable multithreading and/or SIMD of the loop depending on the parameters.
When using SIMD, `@fastmath` and `@inbounds` are used.

In order to use SIMD and multithreading at the same time, the range of the loop is split in even 
batches.
Each batch has a size of `params.simd_batch` iterations, meaning that the inner `@simd` loop has a
fixed number of iterations, while the outer threaded loop will have `N ÷ params.simd_batch`
iterations.

The loop range is assumed to be increasing, i.e. this is correct: 1:2:100, this is not: 100:-2:1
The inner `@simd` loop assumes there is no dependencies between each iteration.

```julia
    @simd_threaded_loop for i = 1:n
        y[i] = log10(x[i]) + x[i]
    end
```
"""
macro simd_threaded_loop(expr)
    return esc(make_simd_threaded_loop(expr))
end


"""
    @simd_threaded_iter(range, expr)

Same as `@simd_threaded_loop(expr)`, but instead of slicing the range of the for loop in `expr`,
we slice the `range` given as the first parameter and distribute the slices evenly to the threads.

The inner `@simd` loop assumes there is no dependencies between each iteration.

```julia
    @simd_threaded_iter 4:2:100 for i in 1:100
        y[i] = log10(x[i]) + x[i]
    end
    # is equivalent to (without threading and SIMD)
    for j in 4:2:100
        for i in (1:100) .+ (j - 1)
            y[i] = log10(x[i]) + x[i]
        end
    end
```
"""
macro simd_threaded_iter(range, expr)
    return esc(make_simd_threaded_iter(range, expr))
end


kernel_macro_error() = error("Either 'kernel_body_pass!' failed to replace this macro with the proper indexing expression, or you forgot to use `@generic_kernel`")


"""
    @index_1D_lin()

Indexing macro to use in a `@generic_kernel` function. 
Returns a linear index to access the 1D arrays used by the kernel.

Cannot be used alongside `@index_2D_lin()`.
"""
macro index_1D_lin() kernel_macro_error() end


"""
    @index_2D_lin()

Indexing macro to use in a `@generic_kernel` function. 
Returns a linear index to access the 2D arrays used by the kernel.

Cannot be used alongside `@index_1D_lin()`.
"""
macro index_2D_lin() kernel_macro_error() end


"""
    @iter_idx()

Indexing macro to use in a `@generic_kernel` function. 
Returns a linear index to access the 2D arrays used by the kernel.

Equivalent to KernelAbstractions.jl's `@index(Global, Linear)`.
"""
macro iter_idx()     kernel_macro_error() end


"""
    @kernel_options(options...)

Must be used (once, and explicitly) in the definition of a `@generic_kernel` function.

Gives options for `@generic_kernel` to adjust the resulting functions.

The possible options are:
 - `debug`:
        Prints the generated functions to stdout at compile time.

```julia
@generic_kernel function add_kernel(a, b)
    @kernel_options(debug)
    i = @index_1D_lin()
    a[i] += b[i]
end
```
"""
macro kernel_options(options...) kernel_macro_error() end


"""
    @kernel_init(expr)

Allows to initialize some internal variables of the kernel before the loop.
The given expression must NOT depend on any index. 
You must not use any indexing macro (`@index_1D_lin()`, etc...) in the expression.

All paramerters of the kernel are available during the execution of the init expression.
On GPU, this expression will be executed by each thread.

This is a workaround for a limitation of Polyester.jl which prevents you from typing variables.

```julia
@generic_kernel function simple_kernel(a, b)
    @kernel_init begin
        c::Float64 = sin(0.123)
    end
    i = @index_1D_lin()
    a[i] += b[i] * c
end
```
"""
macro kernel_init(expr) kernel_macro_error() end


struct KernelWithSIMD end
struct KernelWithoutSIMD end
struct KernelWithThreading end
struct KernelWithoutThreading end


const default_kernel_options = Dict(
    :debug => false
)


function options_list_to_dict(raw_options)
    isnothing(raw_options) && return default_kernel_options
    options = Dict{Symbol, Any}()
    for raw_option in raw_options
        if raw_option isa Expr
            @capture(raw_option, opt_ = s_) || error("Expected an option of the form `name=value`, got: $(string(raw_option))")
            name = opt
            val = s
        elseif raw_option isa Symbol
            name = raw_option
            val = true
        else
            error("Unexpected option type: $raw_option")
        end

        if !(name in keys(default_kernel_options))
            error("Unknown option: $name")
        end

        options[name] = val
    end
    return merge(default_kernel_options, options)
end


function kernel_body_pass!(body::Expr, indexing_replacements::Dict{Symbol, Expr}, device::Symbol)
    indexing_type = :none
    used_indexing_types = Set{Symbol}()
    options = nothing
    init_expr = Expr(:block)

    # TODO: call macroexpand on the body to allow macros to customize kernels
    #  => yes but we need to filter out some macros (like the ones from KernelAbstractions)
    #  => maybe filter by the package which defined the macro

    new_body = MacroTools.postwalk(body) do expr
        if @capture(expr, @index_1D_lin())
            stmt_indexing_type = :lin_1D
        elseif @capture(expr, @index_2D_lin())
            stmt_indexing_type = :lin_2D
        elseif @capture(expr, @iter_idx())
            stmt_indexing_type = :iter_idx
        elseif @capture(expr, @kernel_options(opts__))
            options !== nothing && error("@kernel_options can only be used once per kernel.")
            options = opts
            return Expr(:block)
        elseif @capture(expr, @kernel_init(body_))
            init_expr = body
            return Expr(:block)
        elseif expr isa Expr && expr.head == :macrocall
            # Expressions of the form `@index_1D_lin()` can have a line number node in the Expr,
            # which prevents MacroTools.@capture to identify the macro call, this case is handled here.
            if expr.args[1] == Symbol("@index_1D_lin")
                stmt_indexing_type = :lin_1D
            elseif expr.args[1] == Symbol("@index_2D_lin")
                stmt_indexing_type = :lin_2D
            elseif expr.args[1] == Symbol("@iter_idx")
                stmt_indexing_type = :iter_idx
            elseif expr.args[1] == Symbol("@kernel_options")
                # TODO: allow to combine options
                options !== nothing && error("@kernel_options can only be used once per kernel.")
                opts = expr.args[2:end]
                filter!((opt) -> !(opt isa LineNumberNode), opts)
                options = opts
                return Expr(:block)
            elseif expr.args[1] == Symbol("@kernel_init")
                # TODO: allow to combine init blocks
                !isempty(init_expr.args) && error("@kernel_init can only be used once per kernel.")
                init_expr = expr.args[end]
                return Expr(:block)
            else
                return expr
            end
        elseif device == :CPU && @capture(expr, @print(opts__))
            return Expr(:call, :print, opts...)
        else
            return expr
        end

        if stmt_indexing_type == :iter_idx
            # @iter_idx can be mixed with the other indexing macros
        elseif indexing_type == :none
            indexing_type = stmt_indexing_type
        elseif stmt_indexing_type != indexing_type
            error("Cannot mix calls to @index_1D_lin() and @index_2D_lin() in the same kernel")
        end

        push!(used_indexing_types, stmt_indexing_type)

        # Replace the macro accordingly
        return unblock(indexing_replacements[stmt_indexing_type])
    end

    new_body = quote
        @fast $new_body
    end

    return new_body, init_expr, indexing_type, used_indexing_types, options_list_to_dict(options)
end


function pack_struct_fields(args, struct_t)
    struct_args = []
    filtered_args = []
    data_fields = fieldnames(struct_t)
    for arg in args
        arg_name = splitarg(arg)[1]
        if arg_name in data_fields
            push!(struct_args, arg_name)
        else
            push!(filtered_args, arg)
        end
    end
    return filtered_args, struct_args
end


function make_kokkos_kernel_call(func_name, cpu_kernel_def, is_V_in_where, loop_params_names)
    kokkos_def = deepcopy(cpu_kernel_def)
    kokkos_def[:name] = Symbol("kokkos_", kokkos_def[:name])

    _, params_args = pack_struct_fields(kokkos_def[:args], ArmonParameters)

    special_arguments = Dict{Symbol, Any}()
    if func_name === :acoustic_GAD!
        push!(special_arguments, :limiter_tag => (:Cint, :(params.backend_options.limiter_index)))
    elseif func_name === :init_test
        push!(special_arguments, :test_case => (:Cint, :(params.backend_options.test_case_index)))
    elseif func_name === :boundary_conditions!
        # The C++ `boundary_conditions` is iterated a bit differently from the Julia kernel, and the
        # `inner_range` replaces those two explicit arguments
        push!(special_arguments, :stride => (nothing, nothing), :i_start => (nothing, nothing))
    end

    additional_arguments = Dict(
        :init_test => (:T, :(typeof(params.test) <: Sedov ? params.test.r : zero(T)))
    )

    map_type(type) = type in (:Int, :Int64, Int, Int64) ? Expr(:$, :Idx) : type

    kernel_args = []
    kernel_call = []
    ccall_args = []
    ccall_types = []
    for arg in kokkos_def[:args]
        arg_name, arg_type, is_splat, _ = splitarg(arg)        
        is_splat && error("splat argument is incompatible with ccall")

        if !(arg_name in params_args)
            push!(kernel_args, isnothing(arg_name) ? :(::$arg_type) : :($arg_name::$arg_type))
            push!(kernel_call, isnothing(arg_name) ? :($arg_type()) : :($arg_name))
        end

        if arg_name in keys(special_arguments)
            ccall_type, ccall_arg = special_arguments[arg_name]
            isnothing(ccall_type) && continue
            push!(ccall_args, ccall_arg)
            push!(ccall_types, ccall_type)
        elseif isnothing(arg_name)
            # Skip unnamed arguments (method dispatch guides)
        elseif is_V_in_where && arg_type === :V
            push!(ccall_args, arg_name)
            push!(ccall_types, :(Ref{V}))
        elseif arg_name === :params || arg_name === :data
            # Only the useful fields of special struct should be passed to the kernel
        elseif arg_name in loop_params_names
            # Loop params are converted below
        elseif startswith(string(arg_type), r"NTuple")
            # Unpack NTuples and convert their types if needed
            tuple_size = arg_type.args[2]
            tuple_type = arg_type.args[3]
            for i in 1:tuple_size
                push!(ccall_args, :($arg_name[$i]))
                push!(ccall_types, map_type(tuple_type))
            end
        else
            push!(ccall_args, arg_name)
            push!(ccall_types, map_type(arg_type))
        end
    end

    if func_name in keys(additional_arguments)
        type, value = additional_arguments[func_name]
        push!(ccall_args, value)
        push!(ccall_types, map_type(type))
    end

    kokkos_def[:args] = kernel_args

    # Transform the `main_range` and `inner_range` into a `ArmonKokkos.Range` and
    # `ArmonKokkos.InnerRange1D` (or 2D). For convenience they are stored in the `backend_options`
    # as mutable structs, allowing to use `pointer_from_objref`.
    if length(loop_params_names) == 1
        # Iterate a 1D range
        pushfirst!(ccall_args, :(pointer_from_objref(range)), :(pointer_from_objref(inner_range_1D)))
        pushfirst!(ccall_types, :(Ptr{Cvoid}), :(Ptr{Cvoid}))

        if func_name === :boundary_conditions!
            range_start_expr = :(i_start + stride - 1)
            range_step_expr = :(stride)
        else
            range_start_expr = :(first(loop_range) - 1)
            range_step_expr = :(step(loop_range))
        end

        transform_range_expr = quote
            (; range, inner_range_1D) = params.backend_options

            range.start = 0
            range.end = length(loop_range)

            inner_range_1D.start = $range_start_expr
            inner_range_1D.step  = $range_step_expr
        end
    else
        # Iterate a 2D range
        pushfirst!(ccall_args, :(pointer_from_objref(range)), :(pointer_from_objref(inner_range_2D)))
        pushfirst!(ccall_types, :(Ptr{Cvoid}), :(Ptr{Cvoid}))

        transform_range_expr = quote
            (; range, inner_range_2D) = params.backend_options

            range.start = 0
            range.end = length(main_range) * length(inner_range)
        
            inner_range_2D.main_range_start = first(main_range) - 1
            inner_range_2D.main_range_step  = step(main_range)
            inner_range_2D.row_range_start  = first(inner_range) - 1
            inner_range_2D.row_range_length = length(inner_range)
        end
    end

    # Get the C++ kernel symbol
    if endswith(string(func_name), "!")
        func_name = string(func_name)[1:end-1] |> Symbol
    end
    func_name_quote = QuoteNode(func_name)

    kokkos_def[:body] = quote
        $transform_range_expr
        (; $(params_args...)) = params
        # `CCallTypes` is manually replaced below
        ccall(Main.Kokkos.get_symbol(params.backend_options.lib, $func_name_quote), Cvoid, CCallTypes)
    end

    body_ccall_args = kokkos_def[:body].args[end].args
    body_ccall_args[end] = Expr(:tuple, ccall_types...)  # CCallTypes => ($ccall_types...)
    append!(body_ccall_args, ccall_args)

    # The kernel is a @generated function, therefore we need to return the body expression.
    # The `Idx` type is retrieved from ArmonKokkos. Any mention of it in the quoted body interpolates it.
    kokkos_def[:body] = quote
        kokkos_options_t = params.parameters[3]
        Idx = kokkos_options_t.parameters[1]
        $(Expr(:quote, kokkos_def[:body]))
    end

    kokkos_call = Expr(:call, kokkos_def[:name], kernel_call...)

    return kokkos_def, kokkos_call
end


function transform_kernel(func::Expr)
    def = splitdef(func)
    func_name = def[:name]
    kernel_func_name = Symbol(func_name, "_kernel")
    def[:name] = kernel_func_name

    loop_index_name = gensym(:loop_index)

    # Remove all `@Const()` macros in the argument list. They are kept only for the GPU kernel.
    args = map(def[:args]) do arg
        arg isa Expr || return arg
        (arg.head !== :macrocall || arg.args[1] !== Symbol("@Const")) && return arg
        return last(arg.args)
    end

    is_T_in_where = false
    is_V_in_where = false
    for expr in def[:whereparams]
        if expr === :T
            is_T_in_where = true
        elseif expr === :V
            is_V_in_where = true
        elseif expr isa Expr && @capture(expr, V <: AbstractArray{T})
            is_V_in_where = true
        end
    end

    # -- CPU --

    cpu_def = deepcopy(def)

    cpu_body, init_expr, indexing_type, used_indexing_types, options = kernel_body_pass!(cpu_def[:body], Dict(
        :lin_1D => quote $loop_index_name end,
        :lin_2D => quote $loop_index_name end,
        :iter_idx => quote __j_iter + __i_idx end
    ), :CPU)

    uses_iteration_index = :iter_idx in used_indexing_types

    # Wrap the kernel body with a multi-threaded loop with SIMD
    if indexing_type == :lin_1D
        main_loop = quote 
            for $loop_index_name in loop_range
                $cpu_body
            end
        end

        loop_params = (:(loop_range::OrdinalRange{Int}),)
        loop_params_names = (:loop_range,)
        main_loop_arg = :(loop_range::OrdinalRange{Int})
        main_loop_arg_unpack = Expr(:block)
    elseif indexing_type == :lin_2D
        main_loop = quote 
            for $loop_index_name in inner_range
                $cpu_body
            end
        end

        loop_params = (:(main_range::OrdinalRange{Int}), :(inner_range::OrdinalRange{Int}))
        loop_params_names = (:main_range, :inner_range)
        main_loop_arg = :(range::DomainRange)
        main_loop_arg_unpack = :((main_range, inner_range) = (range.col, range.row))
    else
        error("There is no indexing macro explicitly used in the kernel, therefore the CPU loop cannot be created")
    end

    # Cleanup the LineNumberNode and the redundant block
    main_loop = main_loop |> MacroTools.rmlines |> MacroTools.flatten

    # Add the extra parameters needed for the multi-threaded loop
    cpu_def[:args] = [:(params::ArmonParameters), loop_params..., args...]
    original_cpu_def = deepcopy(cpu_def)

    cpu_def[:name] = Symbol("cpu_$func_name")  # Rename the CPU function

    # Build a function for each possible combination of enable/disable SIMD and threading
    push!(cpu_def[:args], Expr(:threading_switch), Expr(:simd_switch))  # Add two arguments to the kernel
    cpu_block = quote end
    for threading in (:with, :without), simd in (:with, :without)
        cpu_def[:args][end-1] = threading == :with ? :(::KernelWithThreading) : :(::KernelWithoutThreading)
        cpu_def[:args][end]   = simd      == :with ? :(::KernelWithSIMD)      : :(::KernelWithoutSIMD)
        if indexing_type == :lin_1D
            cpu_def[:body] = make_simd_threaded_loop(main_loop; threading, simd, add_iteration_indexes=uses_iteration_index)
        else
            cpu_def[:body] = make_simd_threaded_iter(:main_range, main_loop; threading, simd, add_iteration_indexes=uses_iteration_index)
        end

        if !isempty(init_expr.args)
            # Add the init expression before the `let` block
            cpu_def[:body] = Expr(:block, init_expr, cpu_def[:body])
        end

        push!(cpu_block.args, combinedef(cpu_def))
    end

    # Compute the switches before calling the cpu kernel
    setup_cpu_call = quote
        threading = (!no_threading && params.use_threading) ? KernelWithThreading() : KernelWithoutThreading()
        simd      = params.use_simd ? KernelWithSIMD() : KernelWithoutSIMD()
    end

    # -- Kokkos --

    kokkos_def, kokkos_call = make_kokkos_kernel_call(func_name, original_cpu_def, is_V_in_where, loop_params_names)
    kokkos_block = quote
        @generated $(combinedef(kokkos_def))
    end

    setup_kokkos_call = Expr(:block)

    # -- GPU --

    gpu_def = deepcopy(def)

    var_global_lin = gensym(:I_gl_lin)
    var_1D_lin = gensym(:I_1D_lin)
    var_2D_lin = gensym(:I_2D_lin)

    gpu_def[:body], _, _, _, _ = kernel_body_pass!(gpu_def[:body], Dict(
        :lin_1D => quote $var_1D_lin end,
        :lin_2D => quote $var_2D_lin end,
        :iter_idx => quote $var_global_lin end
    ), :GPU)

    use_1D_lin = :lin_1D in used_indexing_types
    use_2D_lin = :lin_2D in used_indexing_types
    use_global_lin = :iter_idx in used_indexing_types || use_1D_lin || use_2D_lin

    # KernelAbstractions parses only the first layer statements of the kernel body, and doesn't
    # recurse into it. Therefore in order to properly initialize the `@index` macros we store into
    # tmp variables the result of the `@index` macro used in the loop body. Those variables can then
    # be accessed from everywhere in the kernel body.
    indexing_init = quote
        $(use_global_lin ? :($var_global_lin = @index(Global, Linear)) : Expr(:block))
        $(use_1D_lin     ? :($var_1D_lin = $var_global_lin + i_0)      : Expr(:block))
        $(use_2D_lin ? :(
                $var_2D_lin = let
                    ix, iy = divrem($var_global_lin - 1, __ranges_info[4])
                    j = __ranges_info[1] + ix * __ranges_info[2] - 1  # first index in of the row
                    i = __ranges_info[3] + iy + j  # cell index
                    i
                end
            ) : Expr(:block)
        )
    end
    pushfirst!(gpu_def[:body].args, indexing_init.args..., init_expr)

    # Adjust the GPU parameters
    if indexing_type == :lin_1D
        # Note: For the initial index `i_0`, we subtract 1 because of how the main and inner ranges
        # are defined.
        gpu_loop_params = (:(i_0::Int),)
        gpu_loop_params_names = (:(first(loop_range) - 1),)
        gpu_ndrange = :(length(loop_range))
    elseif indexing_type == :lin_2D
        gpu_loop_params = (:(__ranges_info::NTuple{4, Int}),)
        gpu_loop_params_names = (:(
            (first(main_range), step(main_range), first(inner_range), length(inner_range))
        ),)
        gpu_ndrange = :(length(main_range) * length(inner_range))
    end

    pushfirst!(gpu_def[:args], gpu_loop_params...)

    # Define the GPU kernel with KernelAbstractions' @kernel
    gpu_block = quote
        @kernel $(combinedef(gpu_def))
    end

    setup_gpu_call = quote
        gpu_kernel_func = $kernel_func_name(params.device, params.block_size)
        ndrange = ($gpu_ndrange, 1, 1)
    end

    # -- Wrapping function --

    # Create the definition of the main function, which will take care of dispatching the arguments
    # between the different implementations, using the definition of the kernel.
    # The arguments of the main function will be:
    # $func_name(params::ArmonParameters, loop_params..., kernel_args...; kernel_kwargs...)
    main_def = deepcopy(def)
    main_def[:name] = func_name

    # Filter out the arguments to the main function that are present in the ArmonData or 
    # ArmonParameters structs, to then unpack them from a struct instance.
    main_args, data_args = pack_struct_fields(args, ArmonData)
    main_args, params_args = pack_struct_fields(main_args, ArmonParameters)

    params_type = is_T_in_where ? :(ArmonParameters{T}) : :(ArmonParameters)
    data_type   = is_V_in_where ? :(ArmonData{V})       : :(ArmonData)

    if isempty(data_args)
        main_def[:args] = [:(params::$params_type), main_loop_arg, main_args...]
        data_unpack = Expr(:block)
    else
        main_def[:args] = [:(params::$params_type), :(data::$data_type), main_loop_arg, main_args...]
        data_unpack = :((; $(data_args...)) = data)
    end

    params_unpack = isempty(params_args) ? Expr(:block) : :((; $(params_args...)) = params)
 
    # Add our keyword args of the main function
    push!(main_def[:kwargs],
        Expr(:kw, :no_threading, false)
    )

    # Build the parameter list needed to call the CPU or GPU kernel
    call_args = map(Iterators.flatten((args, def[:kwargs]))) do arg
        arg_def = splitarg(arg)
        arg_name = arg_def[1]
        if isnothing(arg_name)
            # Singleton type argument: `::SomeType` -> `SomeType()`
            arg_type = arg_def[2]
            return :($arg_type())
        else
            return arg_name
        end
    end

    # Profiling
    @gensym profiling_state
    profiling_start = quote
        $profiling_state = kernel_start(params, $(QuoteNode(func_name)))
    end

    profiling_end = quote
        kernel_end(params, $(QuoteNode(func_name)), $profiling_state)
    end

    # Build the kernel call expressions
    cpu_call = Expr(:call, cpu_def[:name], :params, loop_params_names..., call_args..., :threading, :simd)
    # Equivalent to: gpu_kernel_func(loop_params_names..., args...; ndrange)
    gpu_call = Expr(:call, :gpu_kernel_func, Expr(:parameters, :ndrange), gpu_loop_params_names..., call_args...)

    call_block = quote
        if params.use_kokkos
            $setup_kokkos_call
        elseif params.use_gpu
            $setup_gpu_call
        else
            $setup_cpu_call
        end

        if params.enable_profiling
            $profiling_start
        end

        try
            return if params.use_kokkos
                $kokkos_call
            elseif params.use_gpu
                $gpu_call
            else
                $cpu_call
            end
        finally
            if params.enable_profiling
                $profiling_end
            end
        end
    end

    # Define the body of the main function
    main_def[:body] = quote
        $params_unpack
        $data_unpack
        $main_loop_arg_unpack
        $call_block
    end

    main_block = quote
        $(combinedef(main_def))
    end

    kernel_def = quote
        $(cpu_block)
        $(kokkos_block)
        $(gpu_block)
        $(main_block)
    end

    if options[:debug]
        println("Generated result for function '$func_name':\n", MacroTools.prettify(kernel_def))
    end

    return esc(kernel_def)
end


# TODO : parameter: fold=<struct> to fold the parameters of the kernel into a `(; <params>) = <struct>` expression
#  + take care of the T and V in where params
"""
    @generic_kernel(function definition)

Transforms a single kernel function into six different functions:
 - 4 which run on the CPU using Polyester.jl's multi-threading or not, as well as SIMD or not.
 - one which uses KernelAbstractions.jl to make a GPU-compatible kernel
 - a main function, which will take care of calling the two others depending if we want to use the
   GPU or not.

To do this, two things are done:
 - All calls to `@index_1D_lin()`, `@index_2D_lin()` and `@iter_idx()` are replaced by their
   equivalent in their respective platforms: a simple loop index for CPUs, and a call to KA.jl's
   `@index` for GPUs.
 - Arguments to each function are edited 

A kernel function must call one of `@index_1D_lin()` or `@index_2D_lin()` at least once, since this 
will determine which type of indexing to use as well as which parameters to add.

The indexing macro `@iter_idx` gives the linear index to the current iteration (on CPU) or global 
thread (on GPU).

This means that in order to call the new main function, one needs to take into account which indexing
macro was used:
 - In all cases, `params::ArmonParameters` is the first argument
 - Then, depending on the indexing macro used:
    - `@index_1D_lin()` : `loop_range::OrdinalRange{Int}`
    - `@index_2D_lin()` : `main_range::OrdinalRange{Int}`, `inner_range::OrdinalRange{Int}`
 - An optional keyword argument, `no_threading`, allows to override the `use_threading`
   parameter, which can be useful in asynchronous contexts. It defaults to `false`.

Using KA.jl's `@Const` to annotate arguments is supported, but they will be present only in the GPU
kernel definition.

Further customisation of the kernel and main function can be obtained using `@kernel_options` and
`@kernel_init`.

```julia
@generic_kernel f_kernel(A)
    i = @index_1D_lin()
    A[i] += 1
end

params.use_gpu = false
f(params, 1:10)  # CPU call

params.use_gpu = true
f(params, 1:10)  # GPU call
```
"""
macro generic_kernel(func)
    return transform_kernel(func)
end
