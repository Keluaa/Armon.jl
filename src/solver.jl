
"""
    SolverStats

Solver output.

`data` is nothing if `parameters.return_data` is `false`.

`timer` is nothing if `parameters.measure_time` is `false`.
"""
struct SolverStats
    final_time::Float64
    last_dt::Float64
    cycles::Int
    solve_time::Float64  # in seconds
    cell_count::Int
    giga_cells_per_sec::Float64
    data::Union{Nothing, BlockGrid}
    timer::Union{Nothing, TimerOutput}
end


function Base.show(io::IO, ::MIME"text/plain", stats::SolverStats)
    println(io, "Solver stats:")
    println(io, " - final time:  ", @sprintf("%.18f", stats.final_time))
    println(io, " - last Δt:     ", @sprintf("%.18f", stats.last_dt))
    println(io, " - cycles:      ", stats.cycles)
    println(io, " - performance: ",
        round(stats.giga_cells_per_sec * 1e3, digits=3), " ×10⁶ cell-cycles/sec ",
        "(", round(stats.solve_time, digits=3), " sec, ", stats.cell_count, " cells)")
    if !isnothing(stats.timer)
        println(io, "Steps time breakdown:")
        show(io, stats.timer; compact=false, allocations=true, sortby=:firstexec)
    end
end


macro checkpoint(step_label)
    esc(:(step_checkpoint(params, state, data, string($step_label))))
end


"""
    block_state_machine(params::ArmonParameters, blk::LocalTaskBlock)

Advances the [`SolverStep`](@ref) state of the `blk`, apply each step of the solver on the `blk`.
This continues until the current cycle is done, or the block needs to wait for another block to do
the ghost cells exchange ([`block_ghost_exchange`](@ref)) or compute the its time step
([`next_time_step`](@ref)).

Returns `true` if we reached the end of the current cycle for `blk`.
"""
function block_state_machine(params::ArmonParameters, blk::LocalTaskBlock)
    @label next_step
    state = blk.state
    blk_state = state.step
    new_state = blk_state
    stop_processing = false

    #=
    Roughly equivalent to:
    ```
    if cycle == 0
        update_EOS!(params, state, blk)
    end
    next_time_step(params, blk)  # Yields to other blocks until this cycle's time step is available
    for (axis, dt_factor) in split_axes(state)
        update_solver_state!(params, state, axis, dt_factor)
        update_EOS!(params, state, blk)
        block_ghost_exchange(params, state, blk)  # Yields to other blocks until all neighbours are updated
        numerical_fluxes!(params, state, blk)
        cell_update!(params, state, blk)
        projection_remap!(params, state, blk)
    end
    # Yield to other blocks until all blocks have finished processing this cycle
    ```
    =#

    # TODO: put `init_test` here too in order to maximize first-touch accuracy and ease of use

    if blk_state == SolverStep.NewCycle
        if start_cycle(state)
            if state.global_dt.cycle == 0
                update_EOS!(params, state, blk)
            end
            new_state = SolverStep.TimeStep
        else
            # Wait for the other blocks to finish the previous cycle
            stop_processing = true
            new_state = SolverStep.NewCycle
        end

    elseif blk_state in (SolverStep.TimeStep, SolverStep.InitTimeStep)
        # If not given at config-time, the time step of the first cycle will be the same as the
        # second cycle, requiring all blocks to finish computing the time step before starting the
        # first cycle, hence the `InitTimeStep` state.
        already_contributed = blk_state == SolverStep.InitTimeStep
        must_wait = next_time_step(params, state, blk; already_contributed)
        if must_wait
            stop_processing = true
            if state.dt == 0
                new_state = SolverStep.InitTimeStep
            end
        else
            new_state = SolverStep.NewSweep
        end

    elseif blk_state == SolverStep.NewSweep
        if next_axis_sweep!(params, state)
            new_state = SolverStep.EndCycle
        else
            new_state = SolverStep.EOS
        end

    elseif blk_state == SolverStep.EOS
        update_EOS!(params, state, blk)
        new_state = SolverStep.Exchange

    elseif blk_state == SolverStep.Exchange
        must_wait = block_ghost_exchange(params, state, blk)
        if must_wait
            stop_processing = true
        else
            new_state = SolverStep.Fluxes
        end

    elseif blk_state == SolverStep.Fluxes
        numerical_fluxes!(params, state, blk)
        new_state = SolverStep.CellUpdate

    elseif blk_state == SolverStep.CellUpdate
        cell_update!(params, state, blk)
        new_state = SolverStep.Remap

    elseif blk_state == SolverStep.Remap
        projection_remap!(params, state, blk)
        new_state = SolverStep.NewSweep

    elseif blk_state == SolverStep.EndCycle
        end_cycle!(state)
        stop_processing = true
        new_state = SolverStep.NewCycle

    else
        error("unknown state: $blk_state")
    end

    state.step = new_state
    !stop_processing && @goto next_step
    return new_state == SolverStep.NewCycle
end


function simple_block_distribution(tid, threads, block_count)
    # TODO: improve by taking into account the individual workload of each block
    # Spread `block_count` among `threads`, the remaining blocks are given to the first
    # `remaining_blocks` threads: there will be max ±1 blocks per thread.
    blocks_per_thread = fld(block_count, threads)
    remaining_blocks = block_count - threads * blocks_per_thread

    prev_tids_blocks = blocks_per_thread * (tid - 1)
    tid_blocks = blocks_per_thread
    if tid > remaining_blocks
        prev_tids_blocks += remaining_blocks
    else
        prev_tids_blocks += tid - 1
        tid_blocks += 1
    end

    return (1:tid_blocks) .+ prev_tids_blocks
end


function solver_cycle_async_init(params, tree::BlockTree)
    is_leaf(tree) && error("BlockTree must have sub blocks")
    length(params.block_tree_levels) != 4 && error("expected a BlockTree of depth 4, got: $(length(params.block_tree_levels))")

    threads_count = params.use_threading ? min(length(tree.sub_blocks), Threads.nthreads()) : 1

    threads_depth = length(params.block_tree_levels) - 1
    if params.block_tree_levels[threads_depth] != threads_count
        error("expected the second-to-last level to have $threads_count blocks, got: $(params.block_tree_levels[threads_depth])")
    end

    return threads_count
end


function solver_cycle_async_init_thread(params, ::BlockTree, _, _)
    return length(params.block_tree_levels) - 1  # threads_depth
end


function solver_cycle_async_step(params, tree::BlockTree, tid, threads_depth)
    all_finished_cycle = true
    visit_all_tree_block(tree; min_depth=threads_depth-1, max_depth=threads_depth-1) do _, sub_t
        tid > length(sub_t.sub_blocks) && return
        sb_all_finished_cycle = apply_until_true(sub_t.sub_blocks[tid]) do blk
            return block_state_machine(params, blk)
        end
        all_finished_cycle &= sb_all_finished_cycle
    end
    return all_finished_cycle
end


function solver_cycle_async_init(params, ::BlockGrid)
    threads_count = params.use_threading ? Threads.nthreads() : 1
    return threads_count
end


function solver_cycle_async_init_thread(_, grid::BlockGrid, tid, threads_count)
    return simple_block_distribution(tid, threads_count, prod(grid.grid_size))  # thread_blocks_idx
end


function solver_cycle_async_step(params, grid::BlockGrid, tid, thread_blocks_idx)
    all_finished_cycle = true
    for blk_idx in thread_blocks_idx
        blk_pos = CartesianIndices(grid.grid_size)[blk_idx]
        # One path for each type of block to avoid runtime dispatch
        if in_grid(blk_pos, grid.static_sized_grid)
            blk = grid.blocks[block_idx(grid, blk_pos)]
            all_finished_cycle &= block_state_machine(params, blk)
        else
            blk = grid.edge_blocks[edge_block_idx(grid, blk_pos)]
            all_finished_cycle &= block_state_machine(params, blk)
        end
    end
    return all_finished_cycle
end


function thread_workload(_, bt::BlockTree, tid, threads_depth)
    tid_level_blocks = 0
    blocks_count = 0

    visit_all_tree_block(bt; min_depth=threads_depth-1, max_depth=threads_depth-1) do _, sub_t
        tid > length(sub_t.sub_blocks) && return
        last_level_blocks += 1
        blocks_count += tree_block_count(sub_t.sub_blocks[tid])
    end

    return tid_level_blocks, blocks_count
end


function thread_workload(_, bg::BlockGrid, tid, threads_depth)
    return 1, length(thread_blocks_idx)
end


"""
    threads_workload(params::ArmonParameters, grid)

Estimate the workload of each thread when using `solver_cycle_async` (i.e. `params.async_cycle == true`).
Returns a `thread_count × 3` `Matrix`: the first column is the thread ID (sorted by ascending order),
the second is the number of batches, the third is the number of [`LocalTaskBlock`](@ref).
"""
function threads_workload(params::ArmonParameters, grid)
    threads_count = solver_cycle_async_init(params, grid)
    workload = Matrix{Int64}(undef, (threads_count, 3))
    workload .= 0

    @threaded :outside_kernel for i in 1:threads_count
        tid = Threads.threadid()
        distrib = solver_cycle_async_init_thread(params, grid, tid, threads_count)
        tid_level_blocks, blocks_count = thread_workload(params, grid, tid, distrib)
        workload[i, 1]  = tid
        workload[i, 2] += tid_level_blocks
        workload[i, 3] += blocks_count
    end

    return sortslices(workload; by=first, dims=1)  # Sort by tid
end


function solver_cycle_async(params::ArmonParameters, grid, block_grid::BlockGrid, max_step_count=typemax(Int))
    timeout = UInt(120e9)  # 120 sec  # TODO: should depend on the total workload, or be deactivatable
    threads_count = solver_cycle_async_init(params, grid)

    # TODO: is Polyester.jl still needed here? The low overhead is not important here
    @threaded :outside_kernel for _ in 1:threads_count
        tid = Threads.threadid()
        distrib = solver_cycle_async_init_thread(params, grid, tid, threads_count)

        t_start = time_ns()
        step_count = 1  # TODO: monitor the maximum `step_count` reached, if it is small, then ok, but with MPI this will not be the case
        while step_count < max_step_count
            if step_count % 100 == 0
                # A safepoint might be needed in some cases as threads waiting for other threads
                # would never allocate and therefore might prevent the GC to run.
                GC.safepoint()

                if time_ns() - t_start > timeout
                    solver_error(:timeout, "cycle took too long in thread $tid")
                end

                # TODO: we are busy waiting for MPI comms/dependencies!! Stop using Polyester in this case and `yield()`!
            end

            # Since blocks may need multiple calls to `block_state_machine` for them to reach the end
            # of the cycle, we must loop until all blocks are updated. One step is one call to
            # `block_state_machine` on every block of `grid`.
            solver_cycle_async_step(params, grid, tid, distrib) && break
            step_count += 1
        end
    end

    # We cannot use checkpoints for each individual blocks, as it would require barriers.
    # Therefore we only rely on a last one at the end of a cycle as well as for the time step.
    step_checkpoint(params, first_state(block_grid), block_grid, "time_step")        && return true
    step_checkpoint(params, first_state(block_grid), block_grid, "projection_remap") && return true
    return false
end


function solver_cycle(params::ArmonParameters, data::BlockGrid)
    state = first_state(data)

    if state.global_dt.cycle == 0
        @checkpoint("init_test") && return true
        @section "EOS_init" update_EOS!(params, state, data)
        @checkpoint("EOS_init") && return true
    end

    (@section "time_step" next_time_step(params, state, data)) && return true
    @checkpoint("time_step") && return true

    @section "$axis" for (axis, dt_factor) in split_axes(state)
        update_solver_state!(params, state, axis, dt_factor)

        @section "EOS" update_EOS!(params, state, data)
        @checkpoint("EOS") && return true

        @section "BC" block_ghost_exchange(params, state, data)
        @checkpoint("boundary_conditions") && return true

        @section "fluxes" numerical_fluxes!(params, state, data)
        @checkpoint("numerical_fluxes") && return true

        @section "update" cell_update!(params, state, data)
        @checkpoint("cell_update") && return true

        @section "remap" projection_remap!(params, state, data)
        @checkpoint("projection_remap") && return true
    end

    return false
end


function time_loop(params::ArmonParameters, grid::BlockGrid)
    (; maxtime, maxcycle, silent, animation_step, is_root, initial_mass, initial_energy) = params

    reset!(grid, params)
    (; global_dt) = grid

    total_cycles_time = 0.
    t1 = time_ns()

    use_block_tree = !isempty(params.block_tree_levels)
    if use_block_tree
        block_tree = BlockTree(grid, params.block_tree_levels)
    end

    # Main solver loop
    while global_dt.time < maxtime && global_dt.cycle < maxcycle
        cycle_start = time_ns()

        stop = @section "solver_cycle" begin
            try
                if params.async_cycle
                    if use_block_tree
                        solver_cycle_async(params, block_tree, grid)
                    else
                        solver_cycle_async(params, grid, grid)
                    end
                else
                    solver_cycle(params, grid)
                end
            catch e
                if e isa SolverException && !params.is_root
                    # Avoid exceeding large error messages by only throwing in the root process.
                    # `SolverException`s must be thrown by all processes during a cycle.
                    true
                else
                    rethrow(e)
                end
            end
        end
        stop && break

        next_cycle!(params, global_dt)

        total_cycles_time += time_ns() - cycle_start

        if is_root
            if silent <= 1
                wait(params)
                current_mass, current_energy = conservation_vars(params, grid)
                ΔM = abs(initial_mass - current_mass)     / initial_mass   * 100
                ΔE = abs(initial_energy - current_energy) / initial_energy * 100
                @printf("Cycle %4d: dt = %.18f, t = %.18f, |ΔM| = %#8.6g%%, |ΔE| = %#8.6g%%\n",
                    global_dt.cycle, global_dt.current_dt, global_dt.time, ΔM, ΔE)
            end
        elseif silent <= 1
            wait(params)
            conservation_vars(params, grid)
        end

        if animation_step != 0 && (global_dt.cycle - 1) % animation_step == 0
            wait(params)
            frame_index = (global_dt.cycle - 1) ÷ animation_step
            frame_file = joinpath("anim", params.output_file) * "_" * @sprintf("%03d", frame_index)
            write_sub_domain_file(params, grid, frame_file)
        end
    end

    @section "Last fence" wait(params)

    t2 = time_ns()

    solve_time = t2 - t1
    grind_time = solve_time / (global_dt.cycle * prod(params.N))

    if is_root
        if silent < 3
            println(" ")
            println("Total time:  ", round(solve_time / 1e9, digits=5),        " sec")
            println("Cycles time: ", round(total_cycles_time / 1e9, digits=5), " sec")
            println("Grind time:  ", round(grind_time / 1e3, digits=5),        " µs/cell/cycle")
            println("Cells/sec:   ", round(1 / grind_time * 1e3, digits=5),    " Mega cells/sec")
            println("Cycles:      ", global_dt.cycle)
            println("Last cycle:  ", 
                @sprintf("%.18f", global_dt.time), " sec, Δt=", 
                @sprintf("%.18f", global_dt.current_dt), " sec")
        end
    end

    return global_dt.time, global_dt.current_dt, global_dt.cycle, 1 / grind_time, solve_time
end


"""
    armon(::ArmonParameters)

Main entry point of the solver. Returns a [`SolverStats`](@ref).
"""
function armon(params::ArmonParameters{T}) where T
    (; silent, is_root, timer) = params

    if is_root && silent < 3 && !isinteractive()
        print_parameters(params)
    end

    if params.use_MPI && silent < 3
        (; rank, proc_size, cart_coords) = params

        # Local info
        node_local_comm = MPI.Comm_split_type(params.global_comm, MPI.COMM_TYPE_SHARED, rank)
        local_rank = MPI.Comm_rank(node_local_comm)
        local_size = MPI.Comm_size(node_local_comm)

        is_root && println("\nProcesses info:")
        rank > 0 && MPI.Recv(Bool, rank-1, 1, params.global_comm)
        @printf(" - %2d/%-2d, local: %2d/%-2d, coords: (%2d,%-2d), cores: %3d to %3d\n", 
            rank, proc_size, local_rank, local_size, cart_coords[1], cart_coords[2], 
            minimum(getcpuids()), maximum(getcpuids()))
        rank < proc_size-1 && MPI.Send(true, rank+1, 1, params.global_comm)
    end

    if is_root && params.animation_step != 0
        if isdir("anim")
            rm.("anim/" .* readdir("anim"))
        else
            mkdir("anim")
        end
    end

    if params.measure_time
        reset_timer!(timer)
        enable_timer!(timer)
    end

    @section "init" begin
        # Allocatation and initialisation are separated in order to correctly map the NUMA space
        # using the first-touch policy when working on CPU
        data = @section "alloc" BlockGrid(params)
        @section "init_test" begin
            init_test(params, data)
            wait(params)
        end
    end

    if params.check_result || params.silent <= 1
        @section "Conservation variables" begin
            params.initial_mass, params.initial_energy = conservation_vars(params, data) 
        end
    end

    final_time, dt, cycles, cells_per_sec, solve_time = time_loop(params, data)

    if params.check_result && is_conservative(params.test)
        @section "Conservation variables" begin
            final_mass, final_energy = conservation_vars(params, data) 
        end

        if params.is_root
            Δm = abs(final_mass   - params.initial_mass)   / params.initial_mass
            Δe = abs(final_energy - params.initial_energy) / params.initial_energy

            # Scale the tolerance with the progress in the default test case, therefore roughly
            # accounting for the number of cells (more cells -> slower time step -> more precision).
            rtol = 1e-2 * min(1, final_time / default_max_time(params.test))

            # 1% of relative error, or 10⁻¹¹ of absolute error, whichever is greater.
            Δm_ok = isapprox(Δm, 0; atol=1e-12, rtol)
            Δe_ok = isapprox(Δe, 0; atol=1e-12, rtol)

            if !(Δm_ok && Δe_ok)
                @warn "Mass and energy are not constant, the solution might not be valid!\n\
                    |mₑ-mᵢ|/mᵢ = $(@sprintf("%#8.6g", Δm))\n\
                    |Eₑ-Eᵢ|/Eᵢ = $(@sprintf("%#8.6g", Δe))\n"
            end
        end
    end

    if params.measure_time
        disable_timer!(timer)
    end

    stats = SolverStats(
        final_time, dt, cycles, solve_time / 1e9, prod(params.N), cells_per_sec,
        params.return_data ? data : nothing,
        params.measure_time ? flatten_sections(timer, ("Inner blocks", "Edge blocks")) : nothing
    )

    if params.return_data || params.write_output || params.write_slices
        device_to_host!(data)  # No-op if the host is the device
    end

    params.write_output && write_sub_domain_file(params, data, params.output_file)
    params.write_slices && write_slices_files(params, data, params.output_file)

    if params.measure_time && params.silent < 3 && !isinteractive()
        show(params.timer)
        println()
    end

    return stats
end
