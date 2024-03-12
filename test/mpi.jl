
using MPI
import CUDA
using AMDGPU

MPI.Init()
MPI.Barrier(MPI.COMM_WORLD)


TEST_CUDA_MPI = CUDA.functional()   && parse(Bool, get(ENV, "TEST_CUDA_MPI", "false"))
TEST_ROCM_MPI = AMDGPU.functional() && parse(Bool, get(ENV, "TEST_ROCM_MPI", "false"))

TEST_TYPES_MPI = (Float64,)
TEST_CASES_MPI = (:Sod, :Sod_y, :Sod_circ, #=:Sedov,=# :Bizarrium)


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
    error("NYI for BlockGrid")
    # TODO: use HDF5 for this

    (g_nx, g_ny) = params.global_grid
    (cx, cy) = params.cart_coords
    (; nghost) = params

    # Ranges of the global domain
    global_cols = 1:g_nx
    global_rows = 1:g_ny
    if params.write_ghosts
        global_cols = Armon.inflate(global_cols, nghost)
        global_rows = Armon.inflate(global_rows, nghost)
    end

    # Ranges of the sub-domain
    col_range = 1:ny
    row_range = 1:nx
    if params.write_ghosts
        offset = nghost
    else
        offset = 0
    end

    # Position of the origin and end of this sub-domain
    pos_x = cx * length(row_range) + 1 - offset
    pos_y = cy * length(col_range) + 1 - offset
    end_x = pos_x + length(row_range) - 1 + offset * 2
    end_y = pos_y + length(col_range) - 1 + offset * 2
    col_offset = cy * length(col_range)

    # Ranges of the sub-domain in the global domain
    sub_domain_rows = pos_y:end_y
    cols_before = first(global_cols):(pos_x-1)
    cols_after = (end_x+1):last(global_cols)

    if params.write_ghosts
        col_range = Armon.inflate(col_range, nghost)
        row_range = Armon.inflate(row_range, nghost)
        offset = nghost
    else
        offset = 0
    end

    skip_cells(range) = for _ in range
        skipchars(!=('\n'), file)
        skip(file, 1)  # Skip the '\n'
    end

    for iy in global_rows
        if iy in sub_domain_rows
            skip_cells(cols_before)

            # `col_range = iy:iy` since we can only read one row at a time
            # The offset then transforms `iy` to the local index of the row
            col_range = (iy:iy) .- col_offset
            Armon.read_data_from_file(params, data, col_range, row_range, file)

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
    # Similar to @test in Test.jl
    skip = [kw.args[2] for kw in kws if kw.args[1] === :skip]
    kws = filter(kw -> kw.args[1] ∉ (:skip, :broken), kws)
    length(skip) > 1 && error("'skip' only allowed once")
    length(kws) > 1 && error("Cannot handle keywords other than 'skip'")
    skip = length(skip) > 0 ? first(skip) : false

    # Run `expr` only if !skip, reduce the test result, then only the root prints and calls @test
    return esc(quote
        let comm = $comm, skip::Bool = $skip, test_rank_ok::Int = skip ? false : $expr;
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


function test_neighbour_coords(P, proc_in_grid, global_comm)
    !proc_in_grid && return

    ref_params = ref_params_for_sub_domain(:Sod, Float64, P; global_comm)
    coords = ref_params.cart_coords

    for (coord, sides) in ((coords[1], Armon.sides_along(Armon.Axis.X)), (coords[2], Armon.sides_along(Armon.Axis.Y))), 
            side in (coord % 2 == 0 ? sides : reverse(sides))
        Armon.has_neighbour(ref_params, side) || continue
        neighbour_rank = Armon.neighbour_at(ref_params, side)
        neighbour_coords = zeros(Int, 2)
        MPI.Sendrecv!(collect(coords), neighbour_coords, ref_params.cart_comm;
                      dest=neighbour_rank, source=neighbour_rank)
        neighbour_coords = tuple(neighbour_coords...)

        expected_coords = coords .+ Armon.offset_to(side)
        @test expected_coords == neighbour_coords
        if expected_coords != neighbour_coords
            @debug "[$(ref_params.rank)] $neighbour_rank at $side: expected $expected_coords, got $neighbour_coords"
        end
    end
end

NX = 100
NY = 100
using Printf

function dump_neighbours(P, proc_in_grid, global_comm)
    !proc_in_grid && return

    ref_params = ref_params_for_sub_domain(:Sod, Float64, P; N=(NX, NY), global_comm)

    coords = ref_params.cart_coords

    neighbour_coords = Dict{Armon.Side.T, Tuple{Int, Int}}()

    for (coord, sides) in ((coords[1], Armon.sides_along(Armon.Axis.X)), (coords[2], Armon.sides_along(Armon.Axis.Y))), 
            side in (coord % 2 == 0 ? sides : reverse(sides))
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
    for (iy, j) in enumerate(domain.col), ix in 1:length(domain.row)
        i = ix + (j - 1)
        array[i] = val + iy * 1000 + ix
    end
end


function check_domain_idx(array, domain, val, my_rank, neighbour_rank)
    diff_count = 0
    for (iy, j) in enumerate(domain.col), ix in 1:length(domain.row)
        i = ix + j - 1
        expected_val = val + iy * 1000 + ix
        if expected_val != array[i]
            diff_count += 1
            @debug "[$my_rank] With $neighbour_rank: at ($i,$j) (or $ix,$iy), expected $expected_val, got $(Int(array[i]))"
        end
    end
    return diff_count
end


function test_halo_exchange(P, proc_in_grid, global_comm)
    !proc_in_grid && return

    ref_params = ref_params_for_sub_domain(:Sod, Int64, P; N=(NX, NY), global_comm)
    data = BlockGrid(ref_params)
    coords = ref_params.cart_coords

    for (coord, sides) in ((coords[1], Armon.sides_along(Armon.Axis.X)), (coords[2], Armon.sides_along(Armon.Axis.Y))), 
            side in (coord % 2 == 0 ? sides : reverse(sides))
        Armon.has_neighbour(ref_params, side) || continue
        neighbour_rank = Armon.neighbour_at(ref_params, side)

        # TODO: redo
        # Fill the domain we send with predictable data, with indexes encoded into the data
        domain = Armon.border_domain(ref_params, side)
        fill_domain_idx(device(data).rho, domain, ref_params.rank * 1_000_000)

        # "Halo exchange", but with one neighbour at a time
        comm_array = Armon.get_send_comm_array(data, side)
        Armon.read_border_array!(ref_params, data, comm_array, side)
        wait(ref_params)

        requests = data.requests[side]
        MPI.Start(requests.send)
        MPI.Start(requests.recv)

        MPI.Wait(requests.send)
        MPI.Wait(requests.recv)

        comm_array = Armon.get_recv_comm_array(data, side)
        Armon.copy_from_recv_buffer!(data, comm_array, side)
        Armon.write_border_array!(ref_params, data, comm_array, side)
        wait(ref_params)

        # Check if the received array was correctly pasted into our ghost domain
        g_domain = Armon.ghost_domain(ref_params, side)
        diff_count = check_domain_idx(device(data).rho, g_domain, neighbour_rank * 1_000_000, ref_params.rank, neighbour_rank)

        @test diff_count == 0
    end
end


function test_reference(prefix, comm, test, type, P; kwargs...)
    ref_params = ref_params_for_sub_domain(test, type, P; N=(NX, NY), global_comm=comm, kwargs...)

    diff_count, data, ref_data = try
        dt, cycles, data = run_armon_reference(ref_params)
        ref_dt, ref_cycles, ref_data = ref_data_for_sub_domain(ref_params)

        @root_test dt ≈ ref_dt atol=abs_tol(type, ref_params.test) rtol=rel_tol(type, ref_params.test)
        @root_test cycles == ref_cycles

        diff_count, _ = count_differences(ref_params, data, ref_data)

        diff_count, data, ref_data
    catch e
        # We cannot throw exceptions since it would create a deadlock
        println("[$(MPI.Comm_rank(comm))] threw an exception:")
        Base.showerror(stdout, e, catch_backtrace(); backtrace=true)
        -1, nothing, nothing
    end

    if WRITE_FAILED
        global_diff_count = MPI.Allreduce(diff_count, MPI.SUM, comm)
        if global_diff_count > 0 && diff_count >= 0
            prefix *= isempty(prefix) ? "" : "_"
            p_str = join(P, '×')
            Armon.write_sub_domain_file(ref_params, data, "$(prefix)test_$(test)_$(type)_$(p_str)"; no_msg=true)
            Armon.write_sub_domain_file(ref_params, ref_data, "$(prefix)ref_$(test)_$(type)_$(p_str)"; no_msg=true)
        end
        println("[$(MPI.Comm_rank(comm))]: found $diff_count")
    end

    return diff_count == 0
end


# All grids must be able to perfectly divide the number of cells in each direction
# of the reference case (100×100)
domain_combinations = [
    (1, 1),
    (1, 2),
    (1, 4),
    (4, 1),
    (2, 2),
    (4, 4),
    (5, 2),
    (2, 5),
    (5, 5)
]

total_proc_count = MPI.Comm_size(MPI.COMM_WORLD)

@testset "MPI" begin
    @testset "$(join(P, '×'))" for P in domain_combinations
        enough_processes = prod(P) ≤ total_proc_count
        if enough_processes
            is_root && @info "Testing with a $P domain"
            comm, proc_in_grid = set_comm_for_grid(P)
        else
            is_root && @info "Not enough processes to test a $P domain"
            comm, proc_in_grid = MPI.COMM_NULL, false
        end

        # dump_neighbours(P, proc_in_grid, comm)

        @testset "Neighbours" begin
            test_neighbour_coords(P, proc_in_grid, comm)
        end

        @testset "Halo exchange" begin
            test_halo_exchange(P, proc_in_grid, comm)
        end


        @testset "CPU" begin
            @testset "Reference" begin
                @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                    @MPI_test comm begin
                        test_reference("CPU", comm, test, type, P)
                    end skip=!enough_processes || !proc_in_grid
                end
            end

            @testset "Conservation" begin
                @testset "$test" for test in (:Sod, :Sod_y, :Sod_circ)
                    if enough_processes && proc_in_grid
                        ref_params = ref_params_for_sub_domain(test, Float64, P;
                            maxcycle=10000, maxtime=10000,
                            N=(NX, NY), global_comm=comm
                        )

                        data = BlockGrid(ref_params)
                        init_test(ref_params, data)

                        init_mass, init_energy = Armon.conservation_vars(ref_params, data)
                        Armon.time_loop(ref_params, data)
                        end_mass, end_energy = Armon.conservation_vars(ref_params, data)
                    else
                        init_mass, init_energy = 0., 0.
                        end_mass,  end_energy  = 0., 0.
                    end

                    @root_test   init_mass ≈ end_mass    atol=1e-12  skip=!enough_processes
                    @root_test init_energy ≈ end_energy  atol=1e-12  skip=!enough_processes
                end
            end
        end


        @testset "CUDA" begin
            @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                @MPI_test comm begin
                    test_reference("CUDA", comm, test, type, P; use_gpu=true, device=:CUDA)
                end skip=!TEST_CUDA_MPI ||!enough_processes || !proc_in_grid
            end
        end


        @testset "ROCm" begin
            @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                @MPI_test comm begin
                    test_reference("ROCm", comm, test, type, P; use_gpu=true, device=:ROCM)
                end skip=!TEST_ROCM_MPI || !enough_processes || !proc_in_grid
            end
        end


        @testset "Kokkos" begin
            @testset "$test with $type" for type in TEST_TYPES_MPI, test in TEST_CASES_MPI
                @MPI_test comm begin
                    test_reference("kokkos", comm, test, type, P; use_kokkos=true)
                end skip=!TEST_KOKKOS_MPI || !enough_processes || !proc_in_grid
            end
        end

        # TODO: thread pinning tests (no overlaps, no gaps)
        # TODO: GPU assignment tests (overlaps only if there is more processes than GPUs in a node)
        # TODO: add @debug statements for those tests to get a view of the structure of cores&gpus assigned to each rank

        MPI.free(comm)
    end
end
