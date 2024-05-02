
abstract type BackendParams end
struct EmptyParams <: BackendParams end


"""
    ArmonParameters(; options...)

The parameters and current state of the solver.

The state is reset at each call to [`armon`](@ref).

There are many options. Each backend can add their own.


# Options


## Backend and MPI

    device = :CUDA

Device to use. Supported values:
 - `:CPU_HP`: `Polyester.jl` CPU multithreading  (default if `use_gpu=false`)
 - `:CUDA`: `CUDA.jl` GPU (default if `use_gpu=true`)
 - `:ROCM`: `AMDGPU.jl` GPU
 - `:CPU`: `KernelAbstractions.jl` CPU multithreading (using the standard `Threads.jl`)


    use_MPI = true, P = (1, 1), reorder_grid = true, global_comm = nothing

MPI config. The MPI domain will be a process grid of size `P`.
`global_comm` is the global communicator to use, defaults to `MPI.COMM_WORLD`.
`reorder_grid` is passed to `MPI.Cart_create`.


    gpu_aware = true

Store MPI buffers on the device. This requires to use a GPU-aware MPI implementation. Does nothing
when using the CPU only.


## Kernels

    use_threading = true, use_simd = true

Switches for [`CPU_HP`](@ref) kernels.
`use_threading` enables [`@threaded`](@ref) for outer loops.
`use_simd` enables [`@simd_loop`](@ref) for inner loops.


    use_gpu = false

Enables the use of `KernelAbstractions.jl` kernels.


    use_kokkos = false

Use kernels for `Kokkos.jl`.


    use_cache_blocking = true

Separate the domain into semi-independant blocks, improving the cache-locality of memory accesses
and therefore memory throughput.


    async_cycle = false

Apply all steps of the solver to all blocks asynchronously, fully taking advantage of cache blocking.


    block_size = 1024

Size of blocks for cache blocking. Can be a tuple. If `use_cache_blocking == false`, this option
only controls the size of GPU blocks.


    use_two_step_reduction = false

Reduction kernels (`dtCFL_kernel` and `conservation_vars`) use some optimizations to perform the
reduction in a single step. It might cause issues on some GPU backends: a more "gentle" approach
could avoid those by doing it in two steps.


## Profiling

    profiling = Symbol[]

List of profiling callbacks to use:
 - `:TimerOutputs`: `TimerOutputs.jl` sections (added if `measure_time=true`)
 - `:NVTX_sections`: `NVTX.jl` sections
 - `:NVTX_kernels`: `NVTX.jl` sections for kernels
 - `:CUDA_kernels`: equivalent to `CUDA.@profile` in front of all kernels


    measure_time = true

`measure_time=false` can remove any overhead caused by profiling.


    time_async = true

`time_async=false` will add a barrier at the end of every section. Useful for GPU kernels.


## Scheme and CFD solver

    scheme = :GAD, riemann_limiter = :minmod

`scheme` is the Riemann solver scheme to use:
 - `:Godunov` (1st order)
 - `:GAD` (2nd order, with limiter).
 
`riemann_limiter` is the limiter to use for the Riemann solver: `:no_limiter`, `:minmod` or `:superbee`.


    projection = :euler_2nd

Scheme for the Eulerian remap step:
 - `:euler` (1st order)
 - `:euler_2nd` (2nd order, +minmod limiter)


    axis_splitting = :Sequential

Axis splitting to use:
 - `:Sequential`: X then Y
 - `:SequentialSym` (or `:Godunov`): X and Y then Y and X, alternating
 - `:Strang`: ½X, Y, ½X then ½Y, X, ½Y, alternating (½ is for halved time step)
 - `:X_only`
 - `:Y_only`


    N = (10, 10)

Number of cells of the global domain in each axes.


    nghost = 4

Number of ghost cells. Must be greater or equal to the minimum number of ghost cells (min 1,
`scheme=:GAD` adds one, `projection=:euler_2nd` adds one, `scheme=:GAD` + `projection=:euler_2nd`
adds another one)


    Dt = 0., cst_dt = false, dt_on_even_cycles = false

`Dt` is the initial time step, it is computed after initialization by default.
If `cst_dt=true` then the time step is always `Dt` and no reduction over the entire domain occurs. 
If `dt_on_even_cycles=true` then then time step is only updated at even cycles (the first cycle is
even).


    data_type = Float64

Data type for all variables. Should be an `AbstractFloat`.


## Test case and domain

    test = :Sod, domain_size = nothing, origin = nothing

`test` is the test case name to use:
 - `:Sod`: Sod shock tube test
 - `:Sod_y`: Sod shock tube test along the Y axis
 - `:Sod_circ`: Circular Sod shock tube test (centered in the domain)
 - `:Bizarrium`: Bizarrium test, similar to the Sod shock tube but with a special equation of state
 - `:Sedov`: Sedov blast-wave test (centered in the domain, reaches the border at `t=1` by default)
 - `:DebugIndexes`: Set all variables to their index in the global domain. Debug only.


    cfl = 0., maxtime = 0., maxcycle = 500_000

`cfl` defaults to the test's default value, same for `maxtime`.
The solver stops when `t` reaches `maxtime` or `maxcycle` iterations were done (`maxcycle=0` stops
after initialization).


## Output

    silent = 0

`silent=0` for maximum verbosity. `silent=3` doesn't print info at each cycle. `silent=5` doesn't
print anything.


    output_dir = ".", output_file = "output"

`joinpath(output_dir, output_file)` will be path to the output file.


    write_output = false, write_ghosts = false

`write_output=true` will write all `saved_vars()` to the output file.
If `write_ghosts=true`, ghost cells will also be included.


    write_slices = false

Will write all `saved_vars()` to 3 output files, one for the middle X row, another for the middle
Y column, and another for the diagonal. If `write_ghosts=true`, ghost cells will also be included.


    output_precision = nothing

Numbers are saved with `output_precision` digits of precision. Defaults to enough numbers for an
exact decimal representation.


    animation_step = 0

If `animation_step ≥ 1`, then every `animation_step` cycles, variables will be saved as with
`write_output=true`.


    compare = false, is_ref = false, comparison_tolerance = 1e-10

If `compare=true`, then at every sub step of each iteration of the solver all variables will:
 - (`is_ref=false`) be compared with a reference file found in `output_dir`
 - (`is_ref=true`) be saved to a reference file in `output_dir`
When comparing, a relative `comparison_tolerance` (the `rtol` kwarg of `isapprox`) is accepted
between values.


    check_result = false

Check if conservation of mass and energy is verified between initialization and the last iteration.
An error is thrown otherwise. Accepts a relative `comparison_tolerance`.


    return_data = false

If `return_data=true`, then in the [`SolverStats`](@ref) returned by [`armon`](@ref), the `data`
field will contain the [`BlockGrid`](@ref) used by the solver.
"""
mutable struct ArmonParameters{Flt_T, Device, DeviceParams}
    # Test problem type, riemann solver and solver scheme
    test::TestCase
    riemann_scheme::RiemannScheme
    riemann_limiter::Limiter
    projection_scheme::ProjectionScheme
    axis_splitting::SplittingMethod

    # Domain parameters
    nghost::Int
    N::NTuple{2, Int}
    N_origin::NTuple{2, Int}  # Position of the first cell in the global domain
    domain_size::NTuple{2, Flt_T}
    origin::NTuple{2, Flt_T}
    cfl::Flt_T
    Dt::Flt_T
    cst_dt::Bool
    dt_on_even_cycles::Bool
    steps_ranges::Vector{StepsRanges}

    # Bounds
    maxtime::Flt_T
    maxcycle::Int

    # Output
    silent::Int
    output_dir::String
    output_file::String
    write_output::Bool
    write_ghosts::Bool
    write_slices::Bool
    output_precision::Int
    animation_step::Int
    measure_time::Bool
    timer::TimerOutput
    time_async::Bool
    enable_profiling::Bool
    profiling_info::Set{Symbol}
    log_blocks::Bool
    estimated_blk_log_size::Int
    return_data::Bool

    # Performance
    use_threading::Bool
    use_simd::Bool
    use_gpu::Bool
    use_kokkos::Bool
    use_cache_blocking::Bool
    use_two_step_reduction::Bool
    async_cycle::Bool
    device::Device  # A KernelAbstractions.Backend, Kokkos.ExecutionSpace or CPU_HP
    backend_options::DeviceParams
    block_size::NTuple{2, Int}

    # MPI
    use_MPI::Bool
    is_root::Bool
    rank::Int
    root_rank::Int
    proc_size::Int
    proc_dims::NTuple{2, Int}
    global_comm::MPI.Comm
    cart_comm::MPI.Comm
    cart_coords::NTuple{2, Int}  # Coordinates of this process in the cartesian grid (0-indexed)
    neighbours::Dict{Side.T, Int}  # Ranks of the neighbours of this process
    global_grid::NTuple{2, Int}  # Dimensions of the global grid
    reorder_grid::Bool
    gpu_aware::Bool

    # Tests & Comparison
    compare::Bool
    is_ref::Bool
    comparison_tolerance::Float64
    check_result::Bool
    initial_mass::Flt_T
    initial_energy::Flt_T

    function ArmonParameters(; data_type = Float64, N = (10, 10), options...)
        device, options = get_device(; options...)

        params = new{data_type, typeof(device), Any}()
        params.N = N
        params.device = device

        # Each initialization step consumes the options it needs. At the end no option should remain.
        # This allows to easily add new options, as well as for backends to have their own custom
        # set of options.
        options = init_scheme(params; options...)
        options = init_test(params; options...)
        options = init_MPI(params; options...)
        options = init_device(params; options...)
        options = init_profiling(params; options...)
        options = init_indexing(params; options...)
        options = init_output(params; options...)
        options = init_solver_state(params; options...)
        options = init_backend(params, params.device; options...)

        if !isempty(options)
            invalid_options = join(map(k -> "'$k'", keys(options)), ", ", " and ")
            error("$(length(options)) unconsumed options:\n$invalid_options")
        end

        for field in fieldnames(typeof(params))
            isdefined(params, field) && continue
            error("Uninitialized field: $field")
        end

        # TODO: this is ugly, but allows to circumvent a circular dependency between `init_device`,
        # `ArmonParameters` and the Kokkos backend of `@generic_kernel`: this way we can access the
        # index type from the `@generated` function without relying on external functions.
        complete_params = new{data_type, typeof(device), typeof(params.backend_options)}()
        for field in fieldnames(typeof(params))
            setfield!(complete_params, field, getfield(params, field))
        end

        return complete_params
    end
end


function get_device(; device = :CUDA, options...)
    use_kokkos = get(options, :use_kokkos, false)
    use_gpu = get(options, :use_gpu, false)

    if use_kokkos
        device_tag = Val(:Kokkos)
    elseif use_gpu
        device_tag = Val(device)
    else
        device_tag = Val(:CPU_HP)
    end

    return create_device(device_tag), options
end


function init_MPI(params::ArmonParameters;
    use_MPI = true, P = (1, 1), reorder_grid = true, global_comm = nothing, gpu_aware = true,
    options...
)
    global_comm = something(global_comm, MPI.COMM_WORLD)
    params.global_comm = global_comm

    if length(P) != length(params.N)
        solver_error(:config, "Mismatched dimensions: expected a grid of $(length(N)) processes, got: $(length(P))")
    end

    params.use_MPI = use_MPI
    params.reorder_grid = reorder_grid
    params.gpu_aware = gpu_aware

    if use_MPI
        !MPI.Initialized() && solver_error(:config, "'use_MPI=true' but MPI has not yet been initialized")

        params.rank = MPI.Comm_rank(global_comm)
        params.proc_size = MPI.Comm_size(global_comm)
        params.proc_dims = P

        # Create a cartesian grid communicator of P processes. `reorder=true` can be very important
        # for performance since it will optimize the layout of the processes.
        C_COMM = MPI.Cart_create(global_comm, P; reorder=reorder_grid)
        if C_COMM == MPI.COMM_NULL
            p_str = join(P, '×')
            solver_error(:config, "`MPI_Cart_create` could not create a $p_str cartesian topology \
                                   using $(params.proc_size) processes")
        end

        params.cart_comm = C_COMM
        params.cart_coords = Tuple(MPI.Cart_coords(C_COMM))

        # TODO: dimension agnostic
        params.neighbours = Dict(
            Side.Left   => MPI.Cart_shift(C_COMM, 0, -1)[2],
            Side.Right  => MPI.Cart_shift(C_COMM, 0,  1)[2],
            Side.Bottom => MPI.Cart_shift(C_COMM, 1, -1)[2],
            Side.Top    => MPI.Cart_shift(C_COMM, 1,  1)[2]
        )
    else
        params.rank = 0
        params.proc_size = 1
        params.proc_dims = ntuple(Returns(1), length(params.N))
        params.cart_comm = global_comm
        params.cart_coords = ntuple(Returns(0), length(params.N))
        params.neighbours = Dict(
            Side.Left   => MPI.PROC_NULL,
            Side.Right  => MPI.PROC_NULL,
            Side.Bottom => MPI.PROC_NULL,
            Side.Top    => MPI.PROC_NULL
        )
    end

    params.root_rank = 0
    params.is_root = params.rank == params.root_rank

    return options
end


function init_device(params::ArmonParameters;
    use_threading = true, use_simd = true,
    use_gpu = false, use_kokkos = false,
    block_size = nothing, use_cache_blocking = true, async_cycle = false,
    use_two_step_reduction = false,
    options...
)
    params.use_threading = use_threading
    params.use_simd = use_simd
    params.use_kokkos = use_kokkos
    params.use_gpu = use_gpu
    params.use_cache_blocking = use_cache_blocking
    params.use_two_step_reduction = use_two_step_reduction
    params.async_cycle = async_cycle

    if use_cache_blocking && use_threading && params.use_MPI
        thread_level = MPI.Query_thread()
        if thread_level < MPI.THREAD_MULTIPLE
            solver_error(:config, "Using multithreading with cache blocking requires MPI to be \
                                   initialized with `threadlevel ≥ MPI.THREAD_MULTIPLE`, \
                                   got: $thread_level")
        end
    end

    if !use_cache_blocking
        if use_gpu
            # The literal block size for GPU kernels
            block_size = something(block_size, 1024)
        else
            # Disable cache blocking by using an empty block size
            block_size = (0, 0)
        end
    elseif isnothing(block_size)
        # TODO: Estimate the optimal block size, given the solver's stencils
        if !use_gpu
            block_size = (64, 64)
        else
            # TODO: GPU block size ?? 1024? but how?
            block_size = (32, 32)
        end
    end

    length(block_size) > 2 && solver_error(:config, "Expected `block_size` to contain up to 2 elements, got: $block_size")
    params.block_size = tuple(block_size..., ntuple(Returns(1), 2 - length(block_size))...)

    return options
end


function init_profiling(params::ArmonParameters;
    profiling = Symbol[], measure_time = true, time_async = true,
    log_blocks = false, estimated_blk_log_size = 0,
    options...
)
    params.profiling_info = Set{Symbol}(profiling)
    measure_time && push!(params.profiling_info, :TimerOutputs)
    params.enable_profiling = !isempty(params.profiling_info)

    missing_prof = setdiff(params.profiling_info, (cb.name for cb in Armon.SECTION_PROFILING_CALLBACKS))
    setdiff!(missing_prof, (cb.name for cb in Armon.KERNEL_PROFILING_CALLBACKS))
    if !isempty(missing_prof)
        solver_error(:config, "Unknown profiler$(length(missing_prof) > 1 ? "s" : ""): " * join(missing_prof, ", "))
    end

    params.measure_time = measure_time
    params.time_async = time_async

    params.log_blocks = log_blocks
    if estimated_blk_log_size == 0 && log_blocks
        sweep_count = length(split_axes(params.axis_splitting, data_type(params), 1))
        estimated_blk_log_size = min(params.maxcycle, 1000) * sweep_count
    end
    params.estimated_blk_log_size = log_blocks ? estimated_blk_log_size : 0

    params.timer = TimerOutput()
    disable_timer!(params.timer)

    return options
end


function init_scheme(params::ArmonParameters{T};
    scheme = :GAD, projection = :euler_2nd,
    riemann_limiter = :minmod,
    axis_splitting = :Sequential,
    nghost = 4,
    cst_dt = false, Dt = 0., dt_on_even_cycles = false,
    options...
) where {T}
    if projection isa Symbol
        projection = scheme_from_name(projection)
    elseif !(projection isa ProjectionScheme)
        solver_error(:config, "Expected a ProjectionScheme type or a Symbol, got: $projection")
    end

    if scheme isa Symbol
        scheme = scheme_from_name(scheme)
    elseif !(scheme isa RiemannScheme)
        solver_error(:config, "Expected a RiemannScheme type or a Symbol, got: $scheme")
    end

    if axis_splitting isa Symbol
        axis_splitting = splitting_from_name(axis_splitting)
    elseif !(axis_splitting isa SplittingMethod)
        solver_error(:config, "Expected a SplittingMethod type or a Symbol, got: $axis_splitting")
    end

    if riemann_limiter isa Symbol
        riemann_limiter = limiter_from_name(riemann_limiter)
    elseif !(riemann_limiter isa Limiter)
        solver_error(:config, "Expected a Limiter type or a Symbol, got: $riemann_limiter")
    end

    min_nghost = stencil_width(scheme) * stencil_width(projection)
    if nghost < min_nghost
        solver_error(:config, "Not enough ghost cells for the riemann solver and projection, \
                               at least $min_nghost are needed, got $nghost")
    end

    if cst_dt && Dt == zero(T)
        solver_error(:config, "Dt == 0 with constant step enabled")
    end

    params.nghost = nghost
    params.riemann_scheme = scheme
    params.projection_scheme = projection
    params.riemann_limiter = riemann_limiter
    params.axis_splitting = axis_splitting
    params.cst_dt = cst_dt
    params.Dt = Dt
    params.dt_on_even_cycles = dt_on_even_cycles

    return options
end


function init_test(params::ArmonParameters{T};
    test = :Sod,
    domain_size = nothing, origin = nothing,
    cfl = 0., maxtime = 0., maxcycle = 500_000,
    options...
) where {T}
    if test isa Symbol
        test_type = test_from_name(test)
        test = nothing
    elseif test isa TestCase
        test_type = typeof(test)
    else
        solver_error(:config, "Expected a TestCase type or a symbol, got: $test")
    end

    if isnothing(domain_size)
        domain_size = default_domain_size(test_type)
    end
    params.domain_size = Tuple(T.(domain_size))

    if isnothing(origin)
        origin = default_domain_origin(test_type)
    end
    params.origin = Tuple(T.(origin))

    if isnothing(test)
        Δx = params.domain_size ./ params.N
        test = create_test(Δx, test_type)
    end
    params.test = test
    params.maxcycle = maxcycle

    has_source_term(test) && error("Inhomogeneous test cases are not yet supported")

    params.cfl     = cfl     != 0 ? cfl     : default_CFL(test)
    params.maxtime = maxtime != 0 ? maxtime : default_max_time(test)

    return options
end


function init_indexing(params::ArmonParameters; options...)
    # Dimensions of the global domain
    params.global_grid = params.N

    # Dimensions of an array of the sub-domain
    params.N =
        # Spread the global domain evenly to all processes
        params.global_grid .÷ params.proc_dims .+
        # The processes at the edge of the process grid get the remaining cells
        ifelse.(params.cart_coords .== params.proc_dims .- 1, params.global_grid .% params.proc_dims, 0)

    if any(params.proc_dims .> 1 .&& params.N .< params.nghost)
        # We want more real cells than ghost cells to avoid having to depend on processes farther
        # than the direct neighbours.
        solver_error(:config, "domain $(params.global_grid) is too small to be split by \
                               $(params.proc_dims) processes while keeping more than $(params.nghost) \
                               cells along each axis")
    end

    params.N_origin = params.cart_coords .* (params.global_grid .÷ params.proc_dims) .+ 1

    compute_steps_ranges(params)

    return options
end


function init_output(params::ArmonParameters{T};
    silent = 0, output_dir = ".", output_file = "output",
    write_output = false, write_ghosts = false, write_slices = false, output_precision = nothing,
    animation_step = 0,
    compare = false, is_ref = false, comparison_tolerance = 1e-10,
    check_result = false, return_data = false,
    options...
) where {T}
    if isnothing(output_precision)
        output_precision = T == Float64 ? 17 : 9  # Exact decimal output by default
    end

    params.silent = silent
    params.output_dir = output_dir
    params.output_file = output_file
    params.write_output = write_output
    params.write_ghosts = write_ghosts
    params.write_slices = write_slices
    params.output_precision = output_precision
    params.animation_step = animation_step
    params.compare = compare
    params.is_ref = is_ref
    params.comparison_tolerance = comparison_tolerance
    params.check_result = check_result
    params.return_data = return_data

    return options
end


function init_solver_state(params::ArmonParameters{T}; options...) where {T}
    params.initial_mass = zero(T)
    params.initial_energy = zero(T)
    return options
end


"""
    create_device(::Val{:device_name})

Create a device object from its name.

Default devices:
 - `:CPU`: the CPU backend of `KernelAbstractions.jl`
 - `:CPU_HP`: `Polyester.jl` multithreading

Extensions:
 - `:Kokkos`: the default `Kokkos.jl` device
 - `:CUDA`: the `CUDA.jl` backend of `KernelAbstractions.jl`
 - `:ROCM`: the `AMDGPU.jl` backend of `KernelAbstractions.jl`
"""
function create_device end


create_device(::Val{:CPU}) = CPU()
create_device(::Val{:CPU_HP}) = CPU_HP()


"""
    init_backend(params::ArmonParameters, ::Dev; options...)

Initialize the backend corresponding to the `Dev` device returned by `create_device` using
`options`. Set the `params.backend_options` field.

It must return `options`, with the backend-specific options removed.
"""
function init_backend(params::ArmonParameters, ::Dev; options...) where {Dev}
    params.backend_options = EmptyParams()
    return options
end


function init_backend(params::ArmonParameters, ::CPU; options...)
    # The CPU backend of KernelAbstractions can be useful in some cases for debugging
    params.is_root && @warn "`use_gpu=true` but the device is set to the CPU. \
                              Therefore no kernel will run on a GPU." maxlog=1
    params.backend_options = EmptyParams()
    return options
end


function print_parameter(io::IO, pad::Int, name::String, value; nl=true, suffix="")
    print(io, "  ", rpad(name * ": ", pad), value, suffix)
    nl && println(io)
end


function print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, CPU_HP})
    print_parameter(io, pad, "multithreading", p.use_threading, nl=!p.use_threading)
    if p.use_threading
        println(io, " ($(Threads.nthreads()) $(use_std_lib_threads ? "standard " : "")thread",
            Threads.nthreads() != 1 ? "s" : "", ")")
    end
    print_parameter(io, pad, "use_simd", p.use_simd)
    print_parameter(io, pad, "use_gpu", false)
    print_parameter(io, pad, "use_kokkos", false)
end


function print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, CPU})
    print_parameter(io, pad, "GPU", true, nl=false)
    println(io, ": KA.jl's CPU backend (block size: ", join(p.block_size, '×'), ")")
end


function print_parameters(io::IO, p::ArmonParameters; pad = 20)
    println(io, "Armon parameters:")
    print_parameter(io, pad, "data_type", data_type(p))
    print_device_info(io, pad, p)
    print_parameter(io, pad, "async_cycle", p.async_cycle)
    print_parameter(io, pad, "MPI", p.use_MPI)

    println(io, " ", "─" ^ (pad*2+2))

    print_parameter(io, pad, "test", p.test)
    print_parameter(io, pad, "riemann scheme", p.riemann_scheme, nl=false)
    if uses_limiter(p.riemann_scheme)
        print(io, " (with $(p.riemann_limiter))")
    end
    println(io, ", ", p.projection_scheme)

    print_parameter(io, pad, "axis splitting", p.axis_splitting)
    print_parameter(io, pad, "time step", "", nl=false)
    print(io, p.Dt != 0 ? "starting at $(p.Dt), " :  "initialized automatically, ")
    if p.cst_dt
        println(io, "constant")
    else
        println(io, "updated ", p.dt_on_even_cycles ? "only at even cycles" : "every cycle")
    end
    print_parameter(io, pad, "CFL", p.cfl)

    println(io, " ", "─" ^ (pad*2+2))

    print_parameter(io, pad, "max time", p.maxtime, suffix=" sec")
    print_parameter(io, pad, "max cycle", p.maxcycle)
    print_parameter(io, pad, "measure step times", p.measure_time)
    if !isempty(p.profiling_info)
        profilers = copy(p.profiling_info)
        p.measure_time && delete!(profilers, :TimerOutputs)
        profilers_str = length(profilers) > 0 ? ": " * join(profilers, ", ") : ""
        print_parameter(io, pad, "profiling", (p.enable_profiling ? "ON" : "OFF") * profilers_str)
    end
    print_parameter(io, pad, "block log", p.log_blocks)
    print_parameter(io, pad, "verbosity", p.silent)
    print_parameter(io, pad, "check result", p.check_result)

    if p.write_output
        print_parameter(io, pad, "write output", p.write_output, nl=false)
        print(io, " (precision: $(p.output_precision) digits)")
        println(io, p.write_ghosts ? "with ghosts" : "")
        print_parameter(io, pad, "to", "'$(p.output_file)'")
        p.write_slices && print_parameter(io, pad, "write slices", p.write_slices)
        if p.compare
            print_parameter(io, pad, "compare", p.compare, nl=false)
            println(io, p.is_ref ? ", as reference" : "with $(p.comparison_tolerance) of tolerance")
        end
    end

    println(io, " ", "─" ^ (pad*2+2))

    if p.use_MPI
        print_parameter(io, pad, "global domain", join(p.global_grid, "×"), nl=false)
        print(io, ", spread on a ", join(p.proc_dims, "×"), " process grid")
        print(io, " (", p.reorder_grid ? "" : "not ", "reordered)")
        println(io, " ($(p.proc_size) in total)")
    end

    domain_str = p.use_MPI ? "sub-domain" : "domain"
    N_str = join(p.N, '×')
    print_parameter(io, pad, domain_str, "$N_str with $(p.nghost) ghosts", nl=false)
    println(io, " (", @sprintf("%g", prod(p.N)), " real cells)")

    if p.use_MPI
        print_parameter(io, pad, "coords", join(p.cart_coords, "×"), nl=false)
        print(io, " (rank: ", p.rank, "/", p.proc_size-1, ")")
        neighbours_list = filter(≠(MPI.PROC_NULL) ∘ last, p.neighbours) |> collect .|> first
        neighbours_str = join(neighbours_list, ", ", " and ") |> lowercase
        print(io, ", with $(neighbour_count(p)) neighbour", neighbour_count(p) != 1 ? "s" : "")
        println(io, neighbour_count(p) > 0 ? " on the " * neighbours_str : "")
        print_parameter(io, pad, "gpu aware", p.gpu_aware)
    end

    print_parameter(io, pad, "domain size", join(p.domain_size, " × "), nl=false)
    println(io, ", origin: (", join(p.origin, ", "), ")")

    grid_size, static_sized_grid, _ = grid_dimensions(p)
    print_grid_dimensions(io, grid_size, static_sized_grid, p.block_size, p.N, p.nghost; pad)
end


print_parameters(p::ArmonParameters) = print_parameters(stdout, p)
Base.show(io::IO, p::ArmonParameters) = print_parameters(io::IO, p::ArmonParameters)


"""
    memory_info(params)

The total and free memory the current process can store on the `params.device`.
"""
function memory_info(params::ArmonParameters)
    mem_info = device_memory_info(params.device)
    # TODO: MPI support
    return mem_info
end


"""
    device_memory_info(device)

The total and free memory on the device, in bytes.
"""
function device_memory_info(::Union{CPU_HP, CPU})
    return (
        total = UInt64(Sys.total_physical_memory()),
        free  = UInt64(Sys.free_physical_memory())
    )
end


"""
    data_type(::ArmonParameters{T})

Get `T`, the type used for numbers by the solver
"""
data_type(::ArmonParameters{T}) where T = T


neighbour_at(params::ArmonParameters, side::Side.T) = params.neighbours[side]
has_neighbour(params::ArmonParameters, side::Side.T) = params.neighbours[side] ≠ MPI.PROC_NULL

neighbour_count(params::ArmonParameters) = count(≠(MPI.PROC_NULL), values(params.neighbours))
neighbour_count(params::ArmonParameters, dir::Axis.T) = count(≠(MPI.PROC_NULL), neighbour_at.(params, sides_along(dir)))


# Default copy method
function Base.copy(p::ArmonParameters{T}) where T
    return ArmonParameters([getfield(p, k) for k in fieldnames(ArmonParameters{T})]...)
end


host_array_type(::D) where D = Array
device_array_type(::D) where D = Array


function alloc_host_kwargs(params::ArmonParameters)
    if params.use_kokkos
        return Base.Pairs(NamedTuple{(:use_label,), Tuple{Bool}}((true,)), (:use_label,))
    else
        return Base.Pairs((), ())
    end
end


function alloc_device_kwargs(params::ArmonParameters)
    if params.use_kokkos
        return Base.Pairs(NamedTuple{(:use_label,), Tuple{Bool}}((true,)), (:use_label,))
    else
        return Base.Pairs((), ())
    end
end


function alloc_array_kwargs(; label, kwargs...)
    if haskey(kwargs, :use_label) && kwargs[:use_label]
        return Base.Pairs(NamedTuple{(:label,), Tuple{String}}((label,)), (:label,))
    else
        return Base.Pairs((), ())
    end
end

#
# Steps indexing
#

function compute_steps_ranges(params::ArmonParameters)
    params.steps_ranges = collect(compute_steps_ranges.(instances(Axis.T), params.nghost, Ref(params.projection_scheme)))
end

function compute_steps_ranges(axis::Axis.T, ghosts::Int, projection::ProjectionScheme)
    # Extra cells to compute in each step for the projection
    extra_FLX = stencil_width(projection)
    extra_UP = stencil_width(projection)

    # Real domain
    bl_corner = (0, 0)  # Bottom-Left corner offset
    tr_corner = (0, 0)  # Top-Right   corner offset
    real_range = (bl_corner, tr_corner)
    real_domain = real_range

    # Full domain (real + ghosts), for initialization + first-touch
    full_domain = (bl_corner .- ghosts, bl_corner .+ ghosts)

    # Steps ranges, computed so that there is no need for an extra BC step before the projection
    EOS = real_range  # The BC overwrites any changes to the ghost cells right after

    if axis == Axis.X
        # Fluxes are computed between 'i-s' and 'i', we need one more cell on the right to have all fluxes
        fluxes_bl  = (extra_FLX, 0); fluxes_tr  = (extra_FLX+1, 0)
        cell_up_bl = (extra_UP,  0); cell_up_tr = (extra_UP,    0)
        advec_bl   = (0,         0); advec_tr   = (1,           0)
    else
        fluxes_bl  = (0, extra_FLX); fluxes_tr  = (0, extra_FLX+1)
        cell_up_bl = (0, extra_UP ); cell_up_tr = (0, extra_UP   )
        advec_bl   = (0, 0        ); advec_tr   = (0, 1          )
    end

    fluxes      = (bl_corner .- fluxes_bl,  tr_corner .+ fluxes_tr )
    cell_update = (bl_corner .- cell_up_bl, tr_corner .+ cell_up_tr)
    advection   = (bl_corner .- advec_bl,   tr_corner .+ advec_tr  )
    projection  = real_range

    return StepsRanges(
        axis, real_domain, full_domain,
        EOS, fluxes, cell_update, advection, projection
    )
end

#
# Synchronisation
#

function Base.wait(::ArmonParameters{<:Any, <:Union{CPU, CPU_HP}})
    # CPU backends are synchronous
end


function Base.wait(params::ArmonParameters{<:Any, <:GPU})
    KernelAbstractions.synchronize(params.device)
end
