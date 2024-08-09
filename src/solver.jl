
"""
    SolverStats

Solver output.

`data` is `nothing` if `parameters.return_data` is `false`.

`timer` is `nothing` if `parameters.measure_time` is `false`.

`grid_log` is `nothing` if `parameters.log_blocks` is `false`.
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
    grid_log::Union{Nothing, BlockGridLog}
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

Returns the new step of the block.
If `SolverStep.NewCycle` is returned, the `blk` reached the end of the current cycle and will not
progress any further until all other blocks have reached the same point.
"""
function block_state_machine(params::ArmonParameters, blk::LocalTaskBlock)
    state = blk.state
    steps_completed = 0
    steps_vars = zero(UInt16)
    steps_var_count = 0

    @label next_step
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
    steps_completed += 1
    if params.log_blocks
        axis_dependent, var_flags = SOLVER_STEPS_VARS[blk_state]
        if axis_dependent
            # TODO: dimension agnostic
            var_flags |= state.axis == Axis.X ? STEPS_VARS_FLAGS.u : STEPS_VARS_FLAGS.v
        end
        steps_vars |= var_flags
        steps_var_count += count_ones(var_flags)
    end
    !stop_processing && @goto next_step

    if params.log_blocks
        # `> 1` since `stop_processing` means that we did not finish `blk_state`, apart from `EndCycle`.
        # This allows to exclude stalls from logs.
        if steps_completed > 1
            push_log!(state, BlockLogEvent(state, new_state, steps_completed, steps_vars, steps_var_count))
        else
            state.total_stalls += 1
        end
    end

    return new_state
end


function stop_busy_waiting(params::ArmonParameters, grid::BlockGrid, first_waiting_block::CartesianIndex, stop_count)
    wait_start = time_ns()

    # A safepoint might be needed in some cases as threads waiting for other threads
    # would never allocate and therefore might prevent the GC to run.
    GC.safepoint()

    if params.use_MPI && !iszero(first_waiting_block)
        # MPI_Wait on the `first_waiting_block`'s remote neighbour
        blk = block_at(grid, first_waiting_block)
        for neighbour in blk.neighbours
            !(neighbour isa RemoteTaskBlock) && continue
            MPI.Testall(neighbour.requests) && continue
            # Only wait for a single side, expecting that once one is done, there is more work to do.
            MPI.Waitall(neighbour.requests)
            return time_ns() - wait_start, true
        end
    end

    # Yield to the OS scheduler, incase some multithreading schenanigans are preventing us to
    # continue further (e.g. another process' thread is bound to the same core as this thread).
    # Wait twice as long as the previous time, starting from 2µs and up to 8ms
    µs_to_wait = 2^clamp(stop_count, 1, 13)
    # TODO: the thread is not marked as GC_STATE_SAFE, therefore this wait may prevent GC
    Libc.systemsleep(µs_to_wait * 1e-6)  # this is `usleep` on Linux btw
    return time_ns() - wait_start, false
end


function solver_cycle_async(params::ArmonParameters, grid::BlockGrid, max_step_count=typemax(Int))
    # TODO: use meta-blocks, one for each core/thread, containing a set of `LocalTaskBlock`,
    # with a predefined repartition and device

    timeout = UInt(120e9)  # 120 sec  # TODO: should depend on the total workload, or be deactivatable
    threads_count = params.use_threading ? Threads.nthreads() : 1

    @threaded :outside_kernel for _ in 1:threads_count
        # TODO: thread block iteration should be done along the current axis

        tid = Threads.threadid()
        thread_blocks_idx = grid.threads_workload[tid]

        t_start = time_ns()
        step_count = 0
        no_progress_count = 0
        total_wait_time = 0
        total_mpi_waits = 0
        while step_count < max_step_count
            # Repeatedly parse through all blocks assigned to the current thread, each time advancing
            # them through the solver steps, until all of them are done with the cycle.
            all_finished_cycle = true
            no_progress = true
            first_waiting_block = zero(eltype(thread_blocks_idx))
            for blk_pos in thread_blocks_idx
                # One path for each type of block to avoid runtime dispatch
                if in_grid(blk_pos, grid.static_sized_grid)
                    blk = grid.blocks[block_idx(grid, blk_pos)]
                    prev_state = blk.state.step
                    new_state = block_state_machine(params, blk)
                else
                    blk = grid.edge_blocks[edge_block_idx(grid, blk_pos)]
                    prev_state = blk.state.step
                    new_state = block_state_machine(params, blk)
                end

                all_finished_cycle &= new_state == SolverStep.NewCycle
                no_progress &= prev_state == new_state
                if prev_state == new_state && iszero(first_waiting_block)
                    first_waiting_block = blk_pos
                end
            end
            step_count += 1
            all_finished_cycle && break
            no_progress_count += no_progress

            if no_progress_count % params.busy_wait_limit == 0
                # No block did any progress for more than `params.busy_wait_limit` calls to
                # `block_state_machine`, to prevent deadlocks (caused by MPI or multithreading),
                # we should stop busy waiting.
                if time_ns() - t_start > timeout
                    solver_error(:timeout, "cycle took too long in thread $tid")
                end
                stop_count = no_progress_count ÷ params.busy_wait_limit
                wait_time, waited_for_mpi = stop_busy_waiting(params, grid, first_waiting_block, stop_count)
                total_wait_time += wait_time
                total_mpi_waits += waited_for_mpi
            end
        end

        if params.log_blocks && !isempty(thread_blocks_idx)
            t_end = time_ns()
            stop_count = no_progress_count ÷ params.busy_wait_limit
            push_log!(grid, tid, ThreadLogEvent(
                grid, tid, step_count, no_progress_count, stop_count,
                total_mpi_waits, total_wait_time, t_end - t_start
            ))
        end
    end

    # We cannot use checkpoints for each individual blocks, as it would require barriers.
    # Therefore we only rely on the last one at the end of a cycle as well as the time step.
    step_checkpoint(params, first_state(grid), grid, "time_step")        && return true
    step_checkpoint(params, first_state(grid), grid, "projection_remap") && return true
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

    # Main solver loop
    while global_dt.time < maxtime && global_dt.cycle < maxcycle
        cycle_start = time_ns()

        stop = @section "solver_cycle" begin
            try
                if params.async_cycle
                    solver_cycle_async(params, grid)
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

        rank_info = @sprintf(" - %2d/%-2d, local: %2d/%-2d, coords: (%2d,%-2d), cores: %3d to %3d",
                             rank, proc_size-1, local_rank, local_size-1, cart_coords[1], cart_coords[2],
                             minimum(getcpuids()), maximum(getcpuids()))

        is_root && println("\nProcesses info:")
        rank > 0 && MPI.Recv(Bool, rank-1, 1, params.global_comm)
        println(rank_info)
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
        params.measure_time ? flatten_sections(timer, ("Inner blocks", "Edge blocks")) : nothing,
        params.log_blocks ? collect_logs(data) : nothing
    )

    if params.return_data || params.write_output || params.write_slices
        device_to_host!(data)  # No-op if the host is the device
    end

    params.write_output && write_sub_domain_file(params, data, params.output_file)
    params.write_slices && write_slices_files(params, data, params.output_file)

    if is_root && params.measure_time && params.silent < 3 && !isinteractive()
        show(params.timer)
        println()
    end

    return stats
end
