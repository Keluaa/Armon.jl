
using TOML
using ThreadPinning
using KernelAbstractions
import .Armon: ArmonDualData, memory_required, init_test, MinmodLimiter
import .Armon: @indexing_vars, @i, DomainRange, update_steps_ranges, get_device_array, linear_range


include(joinpath(@__DIR__, "device_info.jl"))


const PERF_FILE_CPU = joinpath(@__DIR__, "performance_cpu.toml")
const PERF_FILE_GPU = joinpath(@__DIR__, "performance_gpu.toml")

const SMALL_SIZE = 1000    # 1e6 cells
const BIG_SIZE   = 10000   # 1e8 cells
const KERNEL_REPEATS = 100


# (name, run on a single row or not, lambda to launch the kernel)
const KERNELS = [
    ("acoustic!",                  false, (p, d, r, _ , deps) -> Armon.acoustic!(                 p, d, "label", r, d.ustar, d.pstar, d.umat;    dependencies=deps) ),
    ("acoustic_GAD!",              false, (p, d, r, dt, deps) -> Armon.acoustic_GAD!(             p, d, "label", r, dt, d.umat, MinmodLimiter(); dependencies=deps) ),
    ("update_perfect_gas_EOS!",    false, (p, d, r, _ , deps) -> Armon.update_perfect_gas_EOS!(   p, d, "label", r, 7/5;                         dependencies=deps) ),
    ("update_bizarrium_EOS!",      false, (p, d, r, _ , deps) -> Armon.update_bizarrium_EOS!(     p, d, "label", r;                              dependencies=deps) ),
    ("cell_update!",               false, (p, d, r, dt, deps) -> Armon.cell_update!(              p, d, linear_range(r), dt, d.umat;             dependencies=deps) ),
    ("cell_update_lagrange!",      false, (p, d, r, dt, deps) -> Armon.cell_update_lagrange!(     p, d, linear_range(r), p.ifin, dt, d.x;        dependencies=deps) ),
    ("euler_projection!",          false, (p, d, r, dt, deps) -> Armon.euler_projection!(         p, d, r, dt, d.work_array_1, d.work_array_2, d.work_array_3, d.work_array_4; dependencies=deps)                ),
    ("first_order_euler_remap!",   false, (p, d, r, dt, deps) -> Armon.first_order_euler_remap!(  p, d, r, dt, d.work_array_1, d.work_array_2, d.work_array_3, d.work_array_4; dependencies=deps)                ),
    ("second_order_euler_remap!",  false, (p, d, r, dt, deps) -> Armon.second_order_euler_remap!( p, d, r, dt, d.work_array_1, d.work_array_2, d.work_array_3, d.work_array_4; dependencies=deps)                ),
    ("boundaryConditions! left",   true,  (p, d, _, _ , deps) -> Armon.boundaryConditions!(       p, d, 1:p.ny, p.row_length, p.index_start + p.idx_row #= @i(0,1) =#,            1, -1., 1.; dependencies=deps) ),
    ("boundaryConditions! bottom", true,  (p, d, _, _ , deps) -> Armon.boundaryConditions!(       p, d, 1:p.nx,            1, p.index_start + p.idx_col #= @i(1,0) =#, p.row_length,  1., 1.; dependencies=deps) ),
    ("read_border_array!",         true,  (p, d, r, _ , deps) -> Armon.read_border_array!(        p, d, r, p.nx, d.tmp_comm_array;                 dependencies=deps) ),
    ("write_border_array!",        true,  (p, d, r, _ , deps) -> Armon.write_border_array!(       p, d, r, p.nx, d.tmp_comm_array;                 dependencies=deps) ),
]


function check_julia_options(cpu_info)
    options = Base.JLOptions()

    if options.opt_level < 3 
        error("Julia's optimisation level should be O3 for performance tests, restart Julia with '-O3'")
    end

    if options.check_bounds != 2
        error("Automatic bounds checking must be disabled, restart Julia with '--check-bounds=no'")
    end

    num_cores_per_socket = cpu_info[:cores] / cpu_info[:sockets]  # We suppose this is a symetric node
    if Threads.nthreads() != num_cores_per_socket
        error("Performance tests must use all available physical cores on a single socket, restart Julia with '-t $num_cores_per_socket'")
    end
end


function get_perf_file(device_type::Symbol)
    if device_type == :CPU
        return PERF_FILE_CPU
    elseif device_type == :GPU
        return PERF_FILE_GPU
    else
        error("Wrong device type: $device_type")
    end
end


function read_performance_data(device_type::Symbol)
    return TOML.parsefile(get_perf_file(device_type))
end


function write_performance_data(device_type::Symbol, data)
    file_name = get_perf_file(device_type)
    tmp_file_name = file_name * "_TMP"
    try
        open(tmp_file_name, "w") do file
            TOML.print(file, data; sorted=true)
        end
        mv(tmp_file_name, file_name; force=true)
    catch
        rm(tmp_file_name; force=true)
        rethrow()
    end
end


function get_perf_data_for_device(device_type::Symbol, device_info)
    perf_data = read_performance_data(device_type)
    info_hash = hash(device_info)
    device_perf_data = get(perf_data, repr(info_hash), nothing)
    return isnothing(device_perf_data) ? nothing : device_perf_data["performance"]
end


function get_performance_params(test::Symbol, type::Type, device_type::Symbol, n::Int)
    if device_type == :CPU
        ArmonParameters(; 
            ieee_bits=sizeof(type)*8,
            test, scheme=:GAD, projection=:euler_2nd, riemann_limiter=:minmod,
            nghost=5, nx=n, ny=n, 
            maxcycle=50, maxtime=1,
            silent=5, write_output=false, measure_time=false,
            use_MPI=false, async_comms=false,
            use_threading=true, use_simd=true, use_gpu=false)
    elseif device_type == :GPU
        gpu_type = get_available_gpu()
        ArmonParameters(; 
            ieee_bits=sizeof(type)*8,
            test, scheme=:GAD, projection=:euler_2nd, riemann_limiter=:minmod,
            nghost=5, nx=n, ny=n, 
            maxcycle=50, maxtime=1,
            silent=5, write_output=false, measure_time=false,
            use_MPI=false, async_comms=false,
            use_gpu=true, device=gpu_type, block_size=1024)
    end
end


function setup_cpu_threads()
    pinthreads(:compact)
end


function has_enough_memory_for(params, memory_available)
    memory_needed = memory_required(params)
    not_enough_mem = memory_needed > memory_available
    if not_enough_mem
        mem_req_str = round(memory_needed / 1e9; digits=1)
        mem_cur_str = round(memory_available / 1e9; digits=1)
        @warn "There is not enough memory on this device to run this measurement: \
               $(params.nx)x$(params.ny) (needs $mem_req_str GB, has $mem_cur_str GB)"
    end
    return not_enough_mem
end


function measure_solver_performance(params)
    stats = armon(params)
    return round(stats.giga_cells_per_sec; sigdigits=4)
end


function setup_kernel_tests(params, memory_available)
    not_enough_mem = has_enough_memory_for(params, memory_available)
    not_enough_mem && return true, nothing

    data = ArmonDualData(params)
    init_test(params, data)
    device_to_host!(data)

    return false, data
end


function measure_kernel_performance(params::ArmonParameters{T}, data::ArmonDualData, 
        kernel_lambda, single_row_kernel) where T
    (; row_length, nghost, nx) = params
    @indexing_vars(params)

    if single_row_kernel
        # Same as `read_border_array!` domain for the bottom side
        main_range = @i(1, 1):row_length:@i(1, nghost)
        inner_range = 1:nx
        range = DomainRange((main_range, inner_range))
    else
        update_steps_ranges(params)
        range = params.steps_ranges.real_domain
    end

    dt = T(1e-6)
    event = NoneEvent()

    # Compile the kernel and make sure it runs as fast as possible
    # Kernels which weren't used during the solver's performance measurement can take a few
    # iterations to be as fast as possible (reason unknown).
    for _ in 1:KERNEL_REPEATS
        event = kernel_lambda(params, data, range, dt, event)
    end
    wait(event)

    kernel_time = @elapsed begin
        for _ in 1:KERNEL_REPEATS
            event = kernel_lambda(params, data, range, dt, event)
        end
        wait(event)  # Wait for the last one since on GPU kernel launches are asynchronous
    end

    return round(kernel_time / KERNEL_REPEATS; sigdigits=4)
end
