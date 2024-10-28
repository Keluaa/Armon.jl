
using Printf
using MPI
using ThreadPinning

if !MPI.Initialized()
    MPI.Init(; threadlevel=:multiple)
end

MPI.Barrier(MPI.COMM_WORLD)
global_rank = MPI.Comm_rank(MPI.COMM_WORLD)

TEST_CUDA_MPI = if parse(Bool, get(ENV, "TEST_CUDA_MPI", "false"))
    import CUDA
    CUDA.functional()
else
    false
end

TEST_ROCM_MPI = if parse(Bool, get(ENV, "TEST_ROCM_MPI", "false"))
    import AMDGPU
    AMDGPU.functional()
else
    false
end


TEST_TYPES_MPI = (Float64,)
TEST_CASES_MPI = (:Sod, :Sod_y, :Sod_circ)
# TEST_CASES_MPI = (:Sod, :Sod_y, :Sod_circ, :Sedov, :Bizarrium)
is_root && @warn "Sedov and Bizarrium test cases are broken with MPI"  # TODO: fix


TEST_KOKKOS_MPI = parse(Bool, get(ENV, "TEST_KOKKOS_MPI", "false"))
if TEST_KOKKOS_MPI
    using Kokkos

    if is_root && !Kokkos.is_initialized()
        Kokkos.load_wrapper_lib()
    end

    MPI.Barrier(MPI.COMM_WORLD)

    if !Kokkos.is_initialized()
        !is_root && Kokkos.load_wrapper_lib(; no_compilation=true, no_git=true)
        Kokkos.set_omp_vars()
        Kokkos.initialize()
    end
end


function read_sub_domain_from_global_domain_file!(params::ArmonParameters, data::BlockGrid, file::IO)
    # TODO: use HDF5 for this

    # Ranges of the global domain
    global_cols = 1:params.global_grid[2]
    global_rows = 1:params.global_grid[1]
    if params.write_ghosts
        global_cols = Armon.inflate(global_cols, params.nghost)
        global_rows = Armon.inflate(global_rows, params.nghost)
    end

    # Position of the origin and end of this sub-domain
    offset = params.write_ghosts ? params.nghost : 0
    pos_x = params.N_origin[1] - offset
    pos_y = params.N_origin[2] - offset
    end_x = pos_x + params.N[1] - 1 + offset * 2
    end_y = pos_y + params.N[2] - 1 + offset * 2
    col_offset = params.cart_coords[2] * params.N[2]

    # Ranges of the sub-domain in the global domain
    sub_domain_rows = pos_y:end_y
    cols_before = first(global_cols):(pos_x-1)
    cols_after = (end_x+1):last(global_cols)

    skip_cells(range) = for _ in range
        skipchars(!=('\n'), file)
        skip(file, 1)  # Skip the '\n'
    end

    for iy in global_rows
        if iy in sub_domain_rows
            skip_cells(cols_before)

            # We only read one row at a time
            col_idx = iy - col_offset  # Transforms `iy` to the local index of the row
            read_domain = (CartesianIndex(1, col_idx), CartesianIndex(params.N[1], col_idx))
            Armon.read_data_from_file(params, data, file, read_domain; global_ghosts=params.write_ghosts)

            skip_cells(cols_after)
        else
            # Skip the entire row, made of g_nx cells (+ ghosts if any)
            skip_cells(global_cols)
        end

        skip(file, 1)  # Skip the additional '\n' at the end of each row of cells
    end
end


function ref_data_for_sub_domain(params::ArmonParameters{T}) where T
    file_path = get_reference_data_file_name(params.test, T)
    ref_data = BlockGrid(params)
    ref_dt::T = 0
    ref_cycles = 0

    open(file_path, "r") do ref_file
        ref_dt = parse(T, readuntil(ref_file, ','))
        ref_cycles = parse(Int, readuntil(ref_file, '\n'))
        read_sub_domain_from_global_domain_file!(params, ref_data, ref_file)
    end

    return ref_dt, ref_cycles, ref_data
end


function ref_params_for_sub_domain(test::Symbol, type::Type, P; overriden_options...)
    ref_options = Dict{Symbol, Any}(
        :use_MPI => true, :P => P, :reorder_grid => true
    )
    merge!(ref_options, overriden_options)
    return get_reference_params(test, type; ref_options...)
end


function non_mpi_params(P, cart_pos; options...)
    # Utility to create an `ArmonParameters` for any sub-domain, but with MPI disabled
    params = ArmonParameters(; options..., use_MPI=false)
    params.proc_size = prod(P)
    params.proc_dims = P
    params.cart_coords = cart_pos
    Armon.init_indexing(params)
    return params
end


function set_comm_for_grid(P)
    new_grid_size = prod(P)
    global_rank = MPI.Comm_rank(MPI.COMM_WORLD)
    # Only the first `new_grid_size` ranks will be part of the new communicator
    in_grid = global_rank < new_grid_size
    color = in_grid ? 0 : MPI.API.MPI_UNDEFINED[]
    sub_comm = MPI.Comm_split(MPI.COMM_WORLD, color, global_rank)
    return sub_comm, in_grid
end


macro MPI_test(comm, expr, kws...)
    # Similar to @test in Test.jl, but the results are gathered to the root process
    skip = [kw.args[2] for kw in kws if kw.args[1] === :skip]
    kws = filter(kw -> kw.args[1] ∉ (:skip, :broken), kws)
    length(skip) > 1 && error("'skip' only allowed once")
    length(kws) > 1 && error("Cannot handle keywords other than 'skip'")
    skip = length(skip) > 0 ? first(skip) : false

    # Run `expr` only if !skip, reduce the test result, then only the root prints and calls @test
    return esc(quote
        let comm = $comm, skip::Bool = $skip, test_rank_ok::Int = if skip
                false
            else
                try
                    $expr
                catch e
                    # Print the error as a single string, to avoid interleaved messages
                    global_rank = MPI.Comm_rank(MPI.COMM_WORLD)
                    local_rank = MPI.Comm_rank(comm)
                    rank_str = "[$global_rank (local: $local_rank)] caught an error: "
                    err_str = "ERROR: " * sprint(showerror, e; context=stdout)
                    bt_str = sprint(Base.show_backtrace, catch_backtrace(); context=stdout)
                    print(rank_str * "\n" * err_str * "\n" * bt_str * "\n")
                    MPI.Abort(MPI.COMM_WORLD, 1)  # Cannot recover in an MPI app
                end
            end;
            test_result = skip ? 0 : MPI.Allreduce(test_rank_ok, MPI.PROD, comm)
            test_result = test_result > 0
            if !test_result && !skip
                # Print which ranks failed the test
                test_results = MPI.Gather(test_rank_ok, 0, comm)
                ranks = MPI.Gather(MPI.Comm_rank(comm), 0, comm)
                if is_root
                    test_results = Bool.(test_results)
                    if !any(test_results)
                        println("All ranks failed this test:")
                    else
                        failed_ranks = ranks[.!test_results]
                        println("$(length(failed_ranks)) ranks failed this test (ranks: $failed_ranks):")
                    end
                end
            end
            is_root && @test test_result skip=skip
        end
    end)
end


macro root_test(expr, kws...)
    return esc(quote
        if is_root
            @test($expr, $(kws...))
        end
    end)
end


function test_neighbour_coords(P, global_comm)
    ref_params = ref_params_for_sub_domain(:Sod, Float64, P; global_comm)
    coords = ref_params.cart_coords

    all_ok = true
    for (coord, axis) in zip(coords, instances(Armon.Axis.T)),
            side in (iseven(coord) ? Armon.sides_along(axis) : reverse(Armon.sides_along(axis)))
        Armon.has_neighbour(ref_params, side) || continue
        neighbour_rank = Armon.neighbour_at(ref_params, side)
        neighbour_coords = zeros(Int, 2)
        MPI.Sendrecv!(collect(coords), neighbour_coords, ref_params.cart_comm;
                      dest=neighbour_rank, source=neighbour_rank)
        neighbour_coords = tuple(neighbour_coords...)

        expected_coords = coords .+ Armon.offset_to(side)
        @test expected_coords == neighbour_coords
        if expected_coords != neighbour_coords
            all_ok = false
            @debug "[$(ref_params.rank)] $neighbour_rank at $side: expected $expected_coords, got $neighbour_coords"
        end
    end

    return all_ok
end


function dump_neighbours(P, proc_in_grid, global_comm)
    !proc_in_grid && return

    ref_params = ref_params_for_sub_domain(:Sod, Float64, P; N=(100, 100), global_comm)
    coords = ref_params.cart_coords
    neighbour_coords = Dict{Armon.Side.T, Tuple{Int, Int}}()

    for (coord, axis) in zip(coords, instances(Armon.Axis.T)),
            side in (iseven(coord) ? Armon.sides_along(axis) : reverse(Armon.sides_along(axis)))
        Armon.has_neighbour(ref_params, side) || continue
        neighbour_rank = Armon.neighbour_at(ref_params, side)
        n_coords = zeros(Int, 2)
        MPI.Sendrecv!(collect(coords), n_coords, ref_params.cart_comm;
                      dest=neighbour_rank, source=neighbour_rank)
        neighbour_coords[side] = tuple(n_coords...)
    end

    MPI.Barrier(ref_params.cart_comm)

    # Wait for the previous rank
    if ref_params.rank > 0
        MPI.Recv(Bool, ref_params.cart_comm; source=ref_params.rank-1)
    end

    println("[$(ref_params.rank)]: $(coords)")
    for side in instances(Armon.Side.T)
        if Armon.has_neighbour(ref_params, side)
            neighbour_rank = Armon.neighbour_at(ref_params, side)
            @printf(" - %6s: [%2d] = %6s (expected: %6s)",
                string(side), neighbour_rank, string(neighbour_coords[side]),
                string(coords .+ Armon.offset_to(side)))
        else
            @printf(" - %6s: ∅", string(side))
        end
        println()
    end

    # Notify the next rank
    if ref_params.rank < prod(P)-1
        MPI.Send(true, ref_params.cart_comm; dest=ref_params.rank+1)
    end

    MPI.Barrier(ref_params.cart_comm)
end


function fill_domain_idx(array, domain, val)
    for (iy, j) in enumerate(domain.col), (ix, i) in enumerate(domain.row)
        idx = i + (j - 1)
        array[idx] = val + iy * 1000 + ix
    end
end


function check_domain_idx(array, domain, val, my_rank, neighbour_rank)
    diff_count = 0
    for (iy, j) in enumerate(domain.col), (ix, i) in enumerate(domain.row)
        idx = i + j - 1
        expected_val = val + iy * 1000 + ix
        if expected_val != array[idx]
            diff_count += 1
            @debug "[$my_rank] With $neighbour_rank: at ($i,$j) (or $ix,$iy), expected $expected_val, got $(array[idx])"
        end
    end
    return diff_count
end


function positions_along(grid::BlockGrid, side::Armon.Side.T)
    axis = Armon.axis_of(side)
    side_pos  = ifelse.(side in Armon.first_sides(), 1, grid.grid_size)
    first_pos = ifelse.(instances(Armon.Axis.T) .== axis, side_pos, 1)
    last_pos  = ifelse.(instances(Armon.Axis.T) .== axis, side_pos, grid.grid_size)
    return CartesianIndex(first_pos):CartesianIndex(last_pos)
end


function test_halo_exchange(P, global_comm; opts...)
    ref_params = ref_params_for_sub_domain(:DebugIndexes, Float64, P; N=(100, 100), global_comm, opts...)
    block_grid = BlockGrid(ref_params)
    coords = ref_params.cart_coords

    if WRITE_FAILED
        Armon.init_test(ref_params, block_grid)
        for blk in Armon.all_blocks(block_grid)
            Armon.block_host_data(blk).ρ .= 0
        end
        Armon.device_to_host!(block_grid)
    end

    total_diff = 0
    for (coord, axis) in zip(coords, instances(Armon.Axis.T)),
            side in (iseven(coord) ? Armon.sides_along(axis) : reverse(Armon.sides_along(axis)))
        Armon.has_neighbour(ref_params, side) || continue
        neighbour_rank = Armon.neighbour_at(ref_params, side)

        for pos in positions_along(block_grid, side)
            @testset let blk = Armon.block_at(block_grid, pos)
                test_array = Armon.block_host_data(blk).ρ

                # Fill the real domain we send with predictable data, with indexes encoded into the data
                domain = Armon.border_domain(blk.size, side; single_strip=false)
                fill_domain_idx(test_array, domain, (ref_params.rank + 1) * 1_000_000)
                Armon.host_to_device!(blk)

                # Halo exchange, but with one neighbour at a time
                remote_blk = blk.neighbours[Int(side)]
                send_buffer = only(Armon.Communications.unsafe_send_buffer(remote_blk.comm_data))
                @root_test length(domain) * length(Armon.comm_vars()) == length(send_buffer)

                @test Armon.start_exchange(ref_params, blk, remote_blk, side)
                Armon.Communications.wait_send_completed(remote_blk.comm_data)
                Armon.Communications.wait_recv_completed(remote_blk.comm_data)
                @test Armon.finish_exchange(ref_params, blk, remote_blk, side)

                # Check if the received array was correctly pasted into our ghost domain
                Armon.device_to_host!(blk)
                ghost_domain = Armon.ghost_domain(blk.size, side; single_strip=false)
                diff_count = check_domain_idx(test_array, ghost_domain, (neighbour_rank + 1) * 1_000_000, ref_params.rank, neighbour_rank)
                total_diff += diff_count
            end
        end
    end

    if WRITE_FAILED
        global_diff_count = MPI.Allreduce(total_diff, MPI.SUM, global_comm)
        if global_diff_count > 0
            p_str = join(P, '×')
            Armon.write_sub_domain_file(
                ref_params, block_grid, "xchg_$(p_str)";
                no_msg=true, all_ghosts=true, vars=(:x, :y, :ρ)
            )
        end
    end

    return total_diff == 0
end


function test_matched_distribution(P, global_comm, parts, grid_size)
    ref_params = ref_params_for_sub_domain(:Sod, Float64, P; global_comm)
    comm = ref_params.cart_comm

    workload = Armon.thread_workload_distribution(
        parts, grid_size;
        scotch=true, simple=false, perimeter_first=false, check=true,
        match_neighbour_domains=comm,
    )
    blk_grid = Armon.block_grid_from_workload(grid_size, workload)

    distrib_count = length.(workload)
    @MPI_test comm sum(distrib_count) == prod(grid_size)

    blk_grid = Armon.block_grid_from_workload(grid_size, workload)
    @MPI_test comm count(==(0), blk_grid) == 0  # All blocks are assigned to a thread

    distrib_ok = Armon.check_matched_distribution(blk_grid, comm; throw_error=false)
    @MPI_test comm distrib_ok

    if WRITE_FAILED && !distrib_ok
        p_str = join(P, '×')
        gs_str = join(grid_size, '×')
        file = "matched_distrib_P=$(p_str)_parts=$(parts)_gs=$(gs_str).grid"
        Armon.write_workload_distribution(file, ref_params, blk_grid)
    end
end


function test_reference(prefix, comm, test, type, P; kwargs...)
    ref_params = ref_params_for_sub_domain(test, type, P; N=(100, 100), global_comm=comm, kwargs...)

    diff_count, data, ref_data = try
        dt, cycles, data = run_armon_reference(ref_params)
        ref_dt, ref_cycles, ref_data = ref_data_for_sub_domain(ref_params)

        atol = abs_tol(type, ref_params.test)
        rtol = rel_tol(type, ref_params.test)
        @root_test dt ≈ ref_dt atol=atol rtol=rtol
        @root_test cycles == ref_cycles

        diff_count, _ = count_differences(ref_params, data, ref_data)
        diff_count += (cycles != ref_cycles) + !isapprox(dt, ref_dt; atol, rtol)

        diff_count, data, ref_data
    catch e
        # We cannot throw exceptions since it would create a deadlock
        err_str = "[$(MPI.Comm_rank(comm))] threw an exception:\n"
        err_str *= "ERROR: " * sprint(showerror, e; context=stdout) * "\n"
        err_str *= sprint(Base.show_backtrace, catch_backtrace(); context=stdout) * "\n"
        print(err_str)
        -1, nothing, nothing
    end

    global_diff_count = MPI.Allreduce(diff_count, MPI.SUM, comm)
    if WRITE_FAILED && global_diff_count > 0
        println("[$(MPI.Comm_rank(comm))]: found $diff_count")
        if global_diff_count > 0 && diff_count >= 0
            prefix *= isempty(prefix) ? "" : "_"
            p_str = join(P, '×')
            Armon.write_sub_domain_file(ref_params, data, "$(prefix)test_$(test)_$(type)_$(p_str)"; no_msg=true)
            Armon.write_sub_domain_file(ref_params, ref_data, "$(prefix)ref_$(test)_$(type)_$(p_str)"; no_msg=true)
        end
    end

    return global_diff_count == 0
end


function test_conservation(test, P, N; maxcycle=10000, maxtime=10000, kwargs...)
    ref_params = ref_params_for_sub_domain(test, Float64, P;
        maxcycle, maxtime, N, kwargs...
    )

    data = BlockGrid(ref_params)
    Armon.init_test(ref_params, data)

    init_mass, init_energy = Armon.conservation_vars(ref_params, data)
    Armon.time_loop(ref_params, data)
    end_mass, end_energy = Armon.conservation_vars(ref_params, data)

    @root_test   init_mass ≈ end_mass   atol=1e-12
    @root_test init_energy ≈ end_energy atol=1e-12

    return true
end


function test_communication_model(
    comm, P, comm_model_name, comm_model_kwargs,
    use_threading, enough_processes, proc_in_grid
)
    comm_model = try
        Armon.Communications.communication_model(comm_model_name; comm_model_kwargs...)
    catch
        # The model cannot be created (e.g. partitioned comms unsupported)
        @test true skip=true
        return
    end
    !Armon.supports_point_to_point(comm_model) && return

    if comm_model isa Armon.Communications.NoCommunicationModel
        # Since processes cannot communicate with each other, we restrict the test to only tests
        # whose result does not depend on communications.
        # Here we place ourselves in the case where we have 1 process along the X axis, with the Sod
        # tube test: the result must be constant along the Y axis, therefore with communications or
        # not, the result must be the same.
        P[1] != 1 && return
        test_cases = TEST_CASES_MPI ∩ (:Sod,)
    else
        test_cases = TEST_CASES_MPI
    end

    use_threading &= Armon.Communications.is_thread_safe(comm_model)
    use_cache_blocking = Armon.Communications.is_async(comm_model)
    async_cycle = use_cache_blocking

    opts = (; use_threading, use_cache_blocking, async_cycle, comm_model, global_comm=comm)

    @testset "Halo exchange" begin
        test_halo_exchange(P, comm; opts...)
    end

    @testset "Reference" begin
        @testset "$test with $type" for type in TEST_TYPES_MPI, test in test_cases
            @MPI_test comm begin
                test_reference("CPU", comm, test, type, P; global_comm=comm, opts...)
            end skip=!enough_processes || !proc_in_grid
        end
    end
end


mpi_precomp_done = false
function local_precompilation()
    mpi_precomp_done && return
    for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
        ref_params = get_reference_params(test, type; use_MPI=false, maxcycle=2)
        run_armon_reference(ref_params)

        ref_params = get_reference_params(test, type; use_MPI=false, maxcycle=2, async_cycle=true)
        run_armon_reference(ref_params)

        ref_params = get_reference_params(test, type; use_MPI=false, maxcycle=2, use_cache_blocking=false)
        run_armon_reference(ref_params)

        if TEST_CUDA_MPI
            ref_params = get_reference_params(test, type; use_MPI=false, maxcycle=2, use_gpu=true, device=:CUDA)
            run_armon_reference(ref_params)
        end

        if TEST_ROCM_MPI
            ref_params = get_reference_params(test, type; use_MPI=false, maxcycle=2, use_gpu=true, device=:ROCM)
            run_armon_reference(ref_params)
        end

        if TEST_KOKKOS_MPI
            ref_params = get_reference_params(test, type; use_MPI=false, maxcycle=2, use_kokkos=true)
            run_armon_reference(ref_params)
        end
    end
    global mpi_precomp_done = true
end


function mpi_local_thread_pinning()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    node_local_comm = MPI.Comm_split_type(MPI.COMM_WORLD, MPI.COMM_TYPE_SHARED, rank)
    local_size = MPI.Comm_size(node_local_comm)
    local_rank = MPI.Comm_rank(node_local_comm)
    is_local_root = local_rank == 0

    if local_size * Threads.nthreads() > ThreadPinning.ncores()
        # Not enough cores on the node
        !is_local_root && return
        println("[$rank] cannot use one Julia thread per CPU core: $local_size ranks on one node with \
                 $(ThreadPinning.ncores()) cores want to use $(Threads.nthreads()) threads each")
        MPI.Abort(MPI.COMM_WORLD, 3)
    end

    cores = (0:Threads.nthreads()-1) .+ local_rank * Threads.nthreads()
    pinthreads(cores; warn=false)

    # If two threads of different ranks are pinned to the same core, then deadlocks are near unavoidable.
    rank_cores = collect(cores)
    all_local_cores = MPI.Gather(rank_cores, node_local_comm)
    all_local_ranks = MPI.Gather(rank, node_local_comm)  # local to global rank
    MPI.free(node_local_comm)
    (!is_local_root || allunique(all_local_cores)) && return

    # Only the local root prints and aborts
    node_name = gethostname()
    core_count = maximum(all_local_cores)
    if core_count > ThreadPinning.ncores()
        # We are pinning to hyperthreads as there is too many threads per rank: same deadlock issue.
        print("[$rank] ranks $(join(all_local_ranks, ", ", " and ")) at node '$node_name' \
               use too many threads and are therefore pinned to hyperthreads.\n")
    end

    assigned_cores = zeros(Int, core_count)
    for (i, core) in enumerate(all_local_cores)
        local_rank = mod1(i, length(rank_cores))
        if assigned_cores[core] == 0
            assigned_cores[core] = local_rank
        else
            print("[$rank] core n°$core of node '$node_name' is already assigned to rank \
                   $(all_local_ranks[local_rank])\n")
        end
    end

    MPI.Abort(MPI.COMM_WORLD, 3)
end


total_proc_count = MPI.Comm_size(MPI.COMM_WORLD)

if parse(Bool, get(ENV, "CI", "false"))
    # Limit the use of CPU cores only in the CI, as it can only run on a single node, with very
    # limited resources, which forces us to use hyperthreads. However, using more CPU-threads than
    # available will cause deadlocks.
    total_core_count = Sys.CPU_THREADS
else
    mpi_local_thread_pinning()

    core_limit = parse(Int, get(ENV, "TEST_MPI_CORE_LIMIT", "0"))
    if core_limit > 0
        total_core_count = core_limit
    else
        total_core_count = typemax(Int)
    end
end


@testset "MPI" begin
    @testset "$(join(P, '×'))" for P in (
        (1, 1),
        (1, 2),
        (1, 4),
        (4, 1),
        (2, 2),
        (4, 4),
        (5, 2),
        (2, 5),
        (5, 5),
    )
        enough_processes = prod(P) ≤ total_proc_count
        if enough_processes
            is_root && @info "Testing with a $P domain"
            comm, proc_in_grid = set_comm_for_grid(P)
        else
            is_root && @info "Not enough processes to test a $P domain"
            comm, proc_in_grid = MPI.COMM_NULL, false
        end

        enough_cores = prod(P) * Threads.nthreads() ≤ total_core_count
        if enough_processes && !enough_cores && is_root
            @warn "Not enough cores to run the test with multithreading: \
                   need $(prod(P) * Threads.nthreads()), got $total_core_count"
        end
        use_threading = enough_cores

        if comm == MPI.COMM_NULL
            # This rank will test nothing for this test iteration.
            # Instead of doing nothing but waiting, we trigger compilation on the non-MPI parts of the solver.
            local_precompilation()
        end

        # dump_neighbours(P, proc_in_grid, comm)

        @testset "Neighbours" begin
            @MPI_test comm begin
                test_neighbour_coords(P, comm)
            end skip=!enough_processes || !proc_in_grid
        end

        @testset "Halo exchange" begin
            @MPI_test comm begin
                test_halo_exchange(P, comm)
            end skip=!enough_processes || !proc_in_grid
        end

        @testset "Match neighbours" begin
            @testset "$threads - $grid_size" for threads in (1, 3, 8, 37), grid_size in ((4, 4), (8, 1), (16, 15))
                (!enough_processes || !proc_in_grid) && continue
                test_matched_distribution(P, comm, threads, grid_size)
            end
        end

        @testset "CPU" begin
            @testset "Reference" begin
                @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                    @MPI_test comm begin
                        test_reference("CPU", comm, test, type, P; use_threading)
                    end skip=!enough_processes || !proc_in_grid
                end
            end

            @testset "No Blocking" begin
                @MPI_test comm begin
                    test_reference("CPU", comm, :Sod_circ, Float64, P; use_threading, use_cache_blocking=false)
                end skip=!enough_processes || !proc_in_grid
            end

            @testset "Conservation" begin
                @MPI_test comm begin
                    test_conservation(:Sod_circ, P, (100, 100); use_threading, global_comm=comm)
                end skip=!enough_processes || !proc_in_grid
            end

            @testset "Async cycle" begin
                @testset "Reference" begin
                    @MPI_test comm begin
                        test_reference("CPU", comm, :Sod_circ, Float64, P; use_threading, async_cycle=true)
                    end skip=!enough_processes || !proc_in_grid
                end

                @testset "No patience" begin
                    @MPI_test comm begin
                        test_reference("CPU", comm, :Sod_circ, Float64, P; use_threading, async_cycle=true, busy_wait_limit=0)
                    end skip=!enough_processes || !proc_in_grid
                end

                @testset "Conservation" begin
                    @MPI_test comm begin
                        test_conservation(:Sod_circ, P, (100, 100); use_threading, async_cycle=true, global_comm=comm)
                    end skip=!enough_processes || !proc_in_grid
                end

                @testset "Thread split comm" begin
                    @MPI_test comm begin
                        test_reference("CPU", comm, :Sod_circ, Float64, P;
                            use_threading, async_cycle=true, global_comm=comm,
                            thread_split_comm=true, workload_distribution=:weighted_sorted_scotch
                        )
                    end skip=!enough_processes || !proc_in_grid
                end
            end

            @testset "Uneven domain" begin
                @testset "$(join(domain, '×'))" for domain in (
                        (107, 113),
                        (20, 20),
                        (37, 241)
                    )
                    @MPI_test comm begin
                        test_conservation(:Sod_circ, P, domain; use_threading, maxcycle=100, global_comm=comm)
                    end skip=!enough_processes || !proc_in_grid
                end
            end

            @testset "Communication models" begin
                @testset "$(comm_model_name)" for (comm_model_name, comm_model_kwargs) in (
                    (:async_safe, (;)),
                    (:async, (;)),
                    (:sync, (;)),
                    (:rma, (;)),
                    (:partitioned, (; partition_size=56*4)),  # "default real cells in a block per axis" * "default number of ghost cells"
                    (:no_comms, (;)),
                )
                    test_communication_model(
                        comm, P, comm_model_name, comm_model_kwargs,
                        use_threading, enough_processes, proc_in_grid
                    )
                end
            end
        end

        @testset "CUDA" begin
            @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                @MPI_test comm begin
                    test_reference("CUDA", comm, test, type, P; use_threading, use_gpu=true, device=:CUDA)
                end skip=!TEST_CUDA_MPI ||!enough_processes || !proc_in_grid
            end
        end

        @testset "ROCm" begin
            @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                @MPI_test comm begin
                    test_reference("ROCm", comm, test, type, P; use_threading, use_gpu=true, device=:ROCM)
                end skip=!TEST_ROCM_MPI || !enough_processes || !proc_in_grid
            end
        end

        @testset "Kokkos" begin
            @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                @MPI_test comm begin
                    test_reference("kokkos", comm, test, type, P; use_threading, use_kokkos=true)
                end skip=!TEST_KOKKOS_MPI || !enough_processes || !proc_in_grid
            end
        end

        # TODO: thread pinning tests (no overlaps, no gaps)
        # TODO: GPU assignment tests (overlaps only if there is more processes than GPUs in a node)
        # TODO: add @debug statements for those tests to get a view of the structure of cores&gpus assigned to each rank

        MPI.free(comm)
    end
end
