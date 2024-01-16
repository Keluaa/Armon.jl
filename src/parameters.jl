
abstract type BackendParams end
struct EmptyParams <: BackendParams end


mutable struct ArmonParameters{Flt_T, Device, DeviceParams}
    # Test problem type, riemann solver and solver scheme
    test::TestCase
    riemann::Symbol
    scheme::Symbol
    riemann_limiter::Limiter

    # Domain parameters
    nghost::Int
    nx::Int
    ny::Int
    dx::Flt_T
    domain_size::NTuple{2, Flt_T}
    origin::NTuple{2, Flt_T}
    cfl::Flt_T
    Dt::Flt_T
    cst_dt::Bool
    dt_on_even_cycles::Bool
    axis_splitting::Symbol
    projection::Symbol

    # Indexing
    row_length::Int
    col_length::Int
    nbcell::Int
    ideb::Int
    ifin::Int
    index_start::Int
    idx_row::Int
    idx_col::Int
    current_axis::Axis
    s::Int  # Stride
    stencil_width::Int

    # Bounds
    maxtime::Flt_T
    maxcycle::Int

    # Current solver state
    cycle::Int
    time::Flt_T
    cycle_dt::Flt_T  # Time step used by kernels, scaled according to the axis splitting
    curr_cycle_dt::Flt_T  # Current unscaled time step
    next_cycle_dt::Flt_T  # Time step of the next cycle
    steps_ranges::StepsRanges

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
    return_data::Bool

    # Performance
    use_threading::Bool
    use_simd::Bool
    use_gpu::Bool
    use_kokkos::Bool
    device::Device  # A KernelAbstractions.Backend, Kokkos.ExecutionSpace or CPU_HP
    backend_options::DeviceParams
    block_size::NTuple{3, Int}
    tasks_storage::Dict{Symbol, Union{Nothing, IdDict}}

    # MPI
    use_MPI::Bool
    is_root::Bool
    rank::Int
    root_rank::Int
    proc_size::Int
    proc_dims::NTuple{2, Int}
    global_comm::MPI.Comm
    cart_comm::MPI.Comm
    cart_coords::NTuple{2, Int}  # Coordinates of this process in the cartesian grid
    neighbours::Dict{Side, Int}  # Ranks of the neighbours of this process
    global_grid::NTuple{2, Int}  # Dimensions (nx, ny) of the global grid
    reorder_grid::Bool
    comm_array_size::Int
    async_comms::Bool

    # Tests & Comparison
    compare::Bool
    is_ref::Bool
    comparison_tolerance::Float64
    debug_indexes::Bool
    check_result::Bool
    initial_mass::Flt_T
    initial_energy::Flt_T

    function ArmonParameters(; ieee_bits = 64, nx = 10, ny = 10, options...)
        flt_type = (ieee_bits == 64) ? Float64 : Float32
        device, options = get_device(; options...)

        params = new{flt_type, typeof(device), Any}()
        params.nx = nx
        params.ny = ny
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

        update_axis_parameters(params, first(split_axes(params))[1])
        update_steps_ranges(params)

        # TODO: this is ugly, but allows to circumvent a circular dependency between `init_device`,
        # `ArmonParameters` and the Kokkos backend of `@generic_kernel`: this way we can access the
        # index type from the `@generated` function without relying on external functions.
        complete_params = new{flt_type, typeof(device), typeof(params.backend_options)}()
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
    use_MPI = true, px = 1, py = 1, reorder_grid = true, global_comm = nothing,
    options...
)
    global_comm = something(global_comm, MPI.COMM_WORLD)
    params.global_comm = global_comm

    (; nx, ny) = params
    if (nx % px != 0) || (ny % py != 0)
        solver_error(:config, "The dimensions of the global domain ($nx x $ny) are not divisible by the number of processors ($px x $py)")
    end

    params.use_MPI = use_MPI
    params.reorder_grid = reorder_grid

    if use_MPI
        !MPI.Initialized() && solver_error(:config, "'use_MPI=true' but MPI has not yet been initialized")

        params.rank = MPI.Comm_rank(global_comm)
        params.proc_size = MPI.Comm_size(global_comm)
        params.proc_dims = (px, py)

        # Create a cartesian grid communicator of px × py processes. reorder=true can be very
        # important for performance since it will optimize the layout of the processes.
        C_COMM = MPI.Cart_create(global_comm, [Int32(px), Int32(py)], [Int32(0), Int32(0)], reorder_grid)
        params.cart_comm = C_COMM
        params.cart_coords = tuple(MPI.Cart_coords(C_COMM)...)

        params.neighbours = Dict(
            Left   => MPI.Cart_shift(C_COMM, 0, -1)[2],
            Right  => MPI.Cart_shift(C_COMM, 0,  1)[2],
            Bottom => MPI.Cart_shift(C_COMM, 1, -1)[2],
            Top    => MPI.Cart_shift(C_COMM, 1,  1)[2]
        )
    else
        params.rank = 0
        params.proc_size = 1
        params.proc_dims = (px, py)
        params.cart_comm = global_comm
        params.cart_coords = (0, 0)
        params.neighbours = Dict(
            Left   => MPI.PROC_NULL,
            Right  => MPI.PROC_NULL,
            Bottom => MPI.PROC_NULL,
            Top    => MPI.PROC_NULL
        )
    end

    params.root_rank = 0
    params.is_root = params.rank == params.root_rank

    return options
end


function init_device(params::ArmonParameters;
    async_comms = false,
    use_threading = true, use_simd = true,
    use_gpu = false, use_kokkos = false,
    block_size = 1024,
    options...
)
    if async_comms && !use_gpu && !use_kokkos
        @warn "Asynchronous communications only work when using a GPU, or with a multithreading \
               backend which supports tasking." maxlog=1
    end

    if async_comms && use_kokkos
        @warn "Asynchronous communications with Kokkos NYI" maxlog=1
    end

    params.async_comms = async_comms
    params.use_threading = use_threading
    params.use_simd = use_simd
    params.use_kokkos = use_kokkos
    params.use_gpu = use_gpu

    length(block_size) > 3 && solver_error(:config, "Expected `block_size` to contain up to 3 elements, got: $block_size")
    params.block_size = tuple(block_size..., ntuple(Returns(1), 3 - length(block_size))...)

    params.tasks_storage = Dict{Symbol, Union{Nothing, IdDict}}()

    return options
end


function init_profiling(params::ArmonParameters;
    profiling = Symbol[], measure_time = true, time_async = true,
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

    params.timer = TimerOutput()
    disable_timer!(params.timer)

    return options
end


function init_scheme(params::ArmonParameters{T};
    scheme = :GAD, projection = :euler,
    riemann = :acoustic, riemann_limiter = :minmod,
    axis_splitting = :Sequential,
    nghost = 2, stencil_width = nothing,
    cst_dt = false, Dt = 0., dt_on_even_cycles = false,
    options...
) where {T}
    min_nghost = 1
    min_nghost += (scheme != :Godunov)
    min_nghost += (projection == :euler_2nd)
    min_nghost += (projection == :euler_2nd) && (scheme != :Godunov)

    if nghost < min_nghost
        solver_error(:config, "Not enough ghost cells for the scheme and/or projection, \
                               at least $min_nghost are needed, got $nghost")
    end

    if projection == :none
        solver_error(:config, "Lagrangian mode unsupported")
    end

    if isnothing(stencil_width)
        stencil_width = min_nghost
    elseif stencil_width < min_nghost
        @warn "The detected minimum stencil width is $min_nghost, but $stencil_width was given. \
               The Boundary conditions might be false." maxlog=1
    elseif stencil_width > nghost
        solver_error(:config, "The stencil width given ($stencil_width) cannot be bigger than the \
                               number of ghost cells ($nghost)")
    end

    if riemann_limiter isa Symbol
        riemann_limiter = limiter_from_name(riemann_limiter)
    elseif !(riemann_limiter isa Limiter)
        solver_error(:config, "Expected a Limiter type or a symbol, got: $riemann_limiter")
    end

    if cst_dt && Dt == zero(T)
        solver_error(:config, "Dt == 0 with constant step enabled")
    end

    params.nghost = nghost
    params.stencil_width = stencil_width
    params.scheme = scheme
    params.projection = projection
    params.riemann = riemann
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
        (sx, sy) = params.domain_size
        Δx::T = sx / params.nx
        Δy::T = sy / params.ny
        test = create_test(Δx, Δy, test_type)
    end
    params.test = test
    params.maxcycle = maxcycle

    params.cfl     = cfl     != 0 ? cfl     : default_CFL(test)
    params.maxtime = maxtime != 0 ? maxtime : default_max_time(test)

    return options
end


function init_indexing(params::ArmonParameters; options...)
    (; nx, ny) = params

    # Dimensions of the global domain
    g_nx = nx
    g_ny = ny

    # Dimensions of an array of the sub-domain
    (px, py) = params.proc_dims
    nx ÷= px
    ny ÷= py
    row_length = params.nghost * 2 + nx
    col_length = params.nghost * 2 + ny

    params.nx = nx
    params.ny = ny
    params.row_length = row_length
    params.col_length = col_length
    params.global_grid = (g_nx, g_ny)

    # First and last index of the real domain of an array
    params.ideb = row_length * params.nghost + params.nghost + 1
    params.ifin = row_length * (ny - 1 + params.nghost) + params.nghost + nx
    params.index_start = params.ideb - row_length - 1  # Used only by the `@i` macro

    # Used only for indexing with the `@i` macro
    params.idx_row = row_length
    params.idx_col = 1

    params.dx = params.domain_size[1] / g_nx

    # Array allocation sizes
    params.nbcell = row_length * col_length
    params.comm_array_size = params.use_MPI ? max(nx, ny) * params.nghost * 7 : 0

    params.steps_ranges = StepsRanges()

    return options
end


function init_output(params::ArmonParameters{T};
    silent = 0, output_dir = ".", output_file = "output",
    write_output = false, write_ghosts = false, write_slices = false, output_precision = nothing,
    animation_step = 0,
    compare = false, is_ref = false, comparison_tolerance = 1e-10, debug_indexes = false,
    check_result = false, return_data = false,
    options...
) where {T}
    if compare && params.async_comms
        solver_error(:config, "Cannot compare when using asynchronous communications")
    end

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
    params.debug_indexes = debug_indexes
    params.check_result = check_result
    params.return_data = return_data

    return options
end


function init_solver_state(params::ArmonParameters{T}; options...) where {T}
    params.cycle = 0
    params.time = zero(T)
    params.cycle_dt = zero(T)
    params.curr_cycle_dt = zero(T)
    params.next_cycle_dt = zero(T)
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
    print_parameter(io, pad, "ieee_bits", sizeof(data_type(p)) * 8)
    print_device_info(io, pad, p)
    print_parameter(io, pad, "MPI", p.use_MPI)

    println(io, " ", "─" ^ (pad*2+2))

    print_parameter(io, pad, "test", p.test)
    print_parameter(io, pad, "solver", p.riemann, nl=false)
    print(io, ", ", p.scheme, " scheme")
    if p.scheme != :Godunov
        print(io, "(with $(p.riemann_limiter))")
    end
    proj_str = p.projection === :euler ? "1ˢᵗ order" : p.projection === :euler_2nd ? "2ⁿᵈ order" : "<unknown>"
    println(io, ", ", proj_str, " projection")
    print_parameter(io, pad, "axis splitting", "", nl=false)
    if p.axis_splitting === :Sequential
        println(io, "X, Y ; X, Y")
    elseif p.axis_splitting === :SequentialSym
        println(io, "X, Y ; Y, X")
    elseif p.axis_splitting === :Strang
        println(io, "½X, Y, ½X ; ½Y, X, ½Y")
    elseif p.axis_splitting === :X_only
        println(io, "X ; X")
    elseif p.axis_splitting === :Y_only
        println(io, "Y ; Y")
    else
        println(io, "<unknown>")
    end
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

    domain_str = p.use_MPI ? "sub-domain" : "domain"
    print_parameter(io, pad, domain_str, "$(p.nx)×$(p.ny) with $(p.nghost) ghosts", nl=false)
    println(io, " (", @sprintf("%g", p.nx * p.ny), " real cells, ", @sprintf("%g", p.nbcell), " in total)")
    print_parameter(io, pad, "domain size", join(p.domain_size, " × "), nl=false)
    println(io, ", origin: (", join(p.origin, ", "), ")")

    if p.use_MPI
        print_parameter(io, pad, "global domain", join(p.global_grid, "×"), nl=false)
        print(io, ", spread on a ", join(p.proc_dims, "×"), " process grid")
        print(io, " (", p.reorder_grid ? "" : "not ", "reordered)")
        println(io, " ($(p.proc_size) in total)")
        print_parameter(io, pad, "coords", join(p.cart_coords, "×"), nl=false)
        print(io, " (rank: ", p.rank, "/", p.proc_size-1, ")")
        neighbours_list = filter(≠(MPI.PROC_NULL) ∘ last, p.neighbours) |> collect .|> first
        neighbours_str = join(neighbours_list, ", ", " and ") |> lowercase
        print(io, ", with $(neighbour_count(p)) neighbour", neighbour_count(p) != 1 ? "s" : "")
        println(io, neighbour_count(p) > 0 ? " on the " * neighbours_str : "")
        print_parameter(io, pad, "async comms", p.async_comms)
    else
        print_parameter(io, pad, "async code path", p.async_comms)
    end
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


neighbour_at(params::ArmonParameters, side::Side) = params.neighbours[side]
has_neighbour(params::ArmonParameters, side::Side) = params.neighbours[side] ≠ MPI.PROC_NULL

neighbour_count(params::ArmonParameters) = count(≠(MPI.PROC_NULL), values(params.neighbours))
neighbour_count(params::ArmonParameters, dir::Axis) = count(≠(MPI.PROC_NULL), neighbour_at.(params, sides_along(dir)))


function grid_coord_along(params::ArmonParameters, dir::Axis = params.current_axis)
    dir == X_axis ? params.cart_coords[1] : params.cart_coords[2]
end


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
# Axis splitting
#

function split_axes(params::ArmonParameters{T}) where T
    axis_1, axis_2 = X_axis, Y_axis
    if iseven(params.cycle)
        axis_1, axis_2 = axis_2, axis_1
    end

    if params.axis_splitting == :Sequential
        return [
            (X_axis, T(1.0)),
            (Y_axis, T(1.0)),
        ]
    elseif params.axis_splitting == :SequentialSym
        return [
            (axis_1, T(1.0)),
            (axis_2, T(1.0)),
        ]
    elseif params.axis_splitting == :Strang
        return [
            (axis_1, T(0.5)),
            (axis_2, T(1.0)),
            (axis_1, T(0.5)),
        ]
    elseif params.axis_splitting == :X_only
        return [(X_axis, T(1.0))]
    elseif params.axis_splitting == :Y_only
        return [(Y_axis, T(1.0))]
    else
        solver_error(:config, "Unknown axis splitting method: $(params.axis_splitting)")
    end
end


function update_axis_parameters(params::ArmonParameters{T}, axis::Axis) where T
    (; row_length, global_grid, domain_size) = params
    (g_nx, g_ny) = global_grid
    (sx, sy) = domain_size

    params.current_axis = axis

    if axis == X_axis
        params.s = 1
        params.dx = sx / g_nx
    else  # axis == Y_axis
        params.s = row_length
        params.dx = sy / g_ny
    end
end

#
# Steps indexing
#

function update_steps_ranges(params::ArmonParameters)
    (; nx, ny, nghost, row_length, current_axis) = params
    @indexing_vars(params)

    ax = current_axis
    steps = params.steps_ranges

    # Extra cells to compute in each step
    extra_FLX = 1
    extra_UP = 1

    if params.projection == :euler
        # No change
    elseif params.projection == :euler_2nd
        extra_FLX += 1
        extra_UP  += 1
    else
        solver_error(:config, "Unknown scheme: $(params.projection)")
    end

    # Real domain
    col_range = @i(1,1):row_length:@i(1,ny)
    row_range = 1:nx
    real_range = DomainRange(col_range, row_range)
    steps.real_domain = real_range

    # Full domain (real + ghosts), for initialization + first-touch
    full_range = inflate_dir(real_range, X_axis, params.nghost)
    full_range = inflate_dir(full_range, Y_axis, params.nghost)
    steps.full_domain = full_range

    # Steps ranges, computed so that there is no need for an extra BC step before the projection
    steps.EOS = real_range  # The BC overwrites any changes to the ghost cells right after
    steps.fluxes = inflate_dir(real_range, ax, extra_FLX)
    steps.cell_update = inflate_dir(real_range, ax, extra_UP)
    steps.advection = expand_dir(real_range, ax, 1)
    steps.projection = real_range

    # Fluxes are computed between 'i-s' and 'i', we need one more cell on the right to have all fluxes
    steps.fluxes = expand_dir(steps.fluxes, ax, 1)

    # Inner ranges: real domain without sides
    steps.inner_EOS = inflate_dir(steps.EOS, ax, -nghost)
    steps.inner_fluxes = steps.inner_EOS

    rt_offset = direction_length(steps.inner_EOS, ax)

    # Outer ranges: sides of the real domain
    if ax == X_axis
        steps.outer_lb_EOS = DomainRange(col_range, row_range[1:min(nx, nghost)])
    else
        steps.outer_lb_EOS = DomainRange(col_range[1:min(ny, nghost)], row_range)
    end

    steps.outer_rt_EOS = shift_dir(steps.outer_lb_EOS, ax, nghost + rt_offset)

    if rt_offset == 0
        # Correction when the side regions overlap
        overlap_width = direction_length(real_range, ax) - 2*nghost
        steps.outer_rt_EOS = expand_dir(steps.outer_rt_EOS, ax, overlap_width)
    end

    steps.outer_lb_fluxes = prepend_dir(steps.outer_lb_EOS, ax, extra_FLX)
    steps.outer_rt_fluxes = expand_dir(steps.outer_rt_EOS, ax, extra_FLX + 1)
end


function boundary_conditions_indexes(params::ArmonParameters, side::Side)
    (; row_length, nx, ny) = params
    @indexing_vars(params)

    stride::Int = 1
    d::Int = 1

    if side == Left
        stride = row_length
        i_start = @i(0,1)
        loop_range = 1:ny
        d = 1
    elseif side == Right
        stride = row_length
        i_start = @i(nx+1,1)
        loop_range = 1:ny
        d = -1
    elseif side == Top
        stride = 1
        i_start = @i(1,ny+1)
        loop_range = 1:nx
        d = -row_length
    elseif side == Bottom
        stride = 1
        i_start = @i(1,0)
        loop_range = 1:nx
        d = row_length
    else
        solver_error(:config, "Unknown side: $side")
    end

    return i_start, loop_range, stride, d
end


function border_domain(params::ArmonParameters, side::Side)
    (; nghost, nx, ny, row_length) = params
    @indexing_vars(params)

    if side == Left
        main_range = @i(1, 1):row_length:@i(1, ny)
        inner_range = 1:nghost
        side_length = ny
    elseif side == Right
        main_range = @i(nx-nghost+1, 1):row_length:@i(nx-nghost+1, ny)
        inner_range = 1:nghost
        side_length = ny
    elseif side == Top
        main_range = @i(1, ny-nghost+1):row_length:@i(1, ny)
        inner_range = 1:nx
        side_length = nx
    elseif side == Bottom
        main_range = @i(1, 1):row_length:@i(1, nghost)
        inner_range = 1:nx
        side_length = nx
    end

    return DomainRange(main_range, inner_range)
end


function ghost_domain(params::ArmonParameters, side::Side)
    (; nghost, nx, ny, row_length) = params
    @indexing_vars(params)

    if side == Left
        main_range = @i(1-nghost, 1):row_length:@i(1-nghost, ny)
        inner_range = 1:nghost
    elseif side == Right
        main_range = @i(nx+1, 1):row_length:@i(nx+1, ny)
        inner_range = 1:nghost
    elseif side == Top
        main_range = @i(1, ny+1):row_length:@i(1, ny+nghost)
        inner_range = 1:nx
    elseif side == Bottom
        main_range = @i(1, 1-nghost):row_length:@i(1, 0)
        inner_range = 1:nx
    end

    return DomainRange(main_range, inner_range)
end


function Base.wait(::ArmonParameters{<:Any, <:Union{CPU, CPU_HP}})
    # CPU backends are synchronous
end


function Base.wait(params::ArmonParameters{<:Any, <:GPU})
    KernelAbstractions.synchronize(params.device)
end
