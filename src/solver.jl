
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
    println(io, " - final time:  ", @sprintf("%.18f", stats.final_time), " sec")
    println(io, " - last Δt:     ", @sprintf("%.18f", stats.last_dt), " sec")
    println(io, " - cycles:      ", stats.cycles)
    println(io, " - performance: ",
        round(stats.giga_cells_per_sec * 1e3, digits=3), " ×10⁶ cells-cycles/sec ",
        "(", round(stats.solve_time, digits=3), " sec, ", stats.cell_count, " cells)")
    if !isnothing(stats.timer)
        println(io, "Steps time breakdown:")
        show(io, stats.timer; compact=false, allocations=true, sortby=:firstexec)
    end
end


macro checkpoint(step_label)
    esc(:(step_checkpoint(params, state, data, string($step_label))))
end


function block_state_machine(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    @label next_step
    blk_state = block_state(blk)
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
    ```
    =#

    # TODO: put `init_test` here too in order to maximize first-touch accuracy and ease of use

    if blk_state == BlockSolverState.NewCycle
        if state.global_dt.cycle == 0
            update_EOS!(params, state, blk)
        end
        new_state = BlockSolverState.TimeStep

    elseif blk_state in (BlockSolverState.TimeStep, BlockSolverState.InitTimeStep)
        # If not given at config-time, the time step of the first cycle will be the same as the
        # second cycle, requiring all blocks to finish computing the time step before starting the
        # first cycle, hence the `InitTimeStep` state.
        already_contributed = blk_state == BlockSolverState.InitTimeStep
        must_wait = next_time_step(params, state, blk; already_contributed)
        if must_wait
            stop_processing = true
            if state.dt == 0
                new_state = BlockSolverState.InitTimeStep
            end
        else
            new_state = BlockSolverState.NewSweep
        end

    elseif blk_state == BlockSolverState.NewSweep
        if next_axis_sweep!(params, state)
            new_state = BlockSolverState.EndCycle
        else
            new_state = BlockSolverState.EOS
        end

    elseif blk_state == BlockSolverState.EOS
        update_EOS!(params, state, blk)
        new_state = BlockSolverState.Exchange

    elseif blk_state == BlockSolverState.Exchange
        must_wait = block_ghost_exchange(params, state, blk)
        if must_wait
            stop_processing = true
        else
            new_state = BlockSolverState.Fluxes
        end

    elseif blk_state == BlockSolverState.Fluxes
        numerical_fluxes!(params, state, blk)
        new_state = BlockSolverState.CellUpdate

    elseif blk_state == BlockSolverState.CellUpdate
        cell_update!(params, state, blk)
        new_state = BlockSolverState.Remap

    elseif blk_state == BlockSolverState.Remap
        projection_remap!(params, state, blk)
        new_state = BlockSolverState.NewSweep

    elseif blk_state == BlockSolverState.EndCycle
        stop_processing = true
        new_state = BlockSolverState.NewCycle

    else
        error("unknown state: $blk_state")
    end

    block_state!(blk, new_state)
    !stop_processing && @goto next_step
    return new_state == BlockSolverState.NewCycle
end


function simple_block_distribution(tid, threads, grid_size)
    # TODO: improve by taking into account the individual workload of each block
    block_count = prod(grid_size)
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


# TODO: remove (debug only)
function test_distrib(f, T, GS)
    blocks = Array{Int, length(N)}(undef, GS)
    blocks .= 0

    for tid in 1:T
        blocks[f(tid)] .= tid
    end

    return blocks
end


function solver_cycle_async(
    params::ArmonParameters, block_states::Vector{SState}, grid::BlockGrid, max_step_count=typemax(Int)
) where {SState <: SolverState}
    # TODO: use meta-blocks, one for each core/thread, containing a set of `LocalTaskBlock`,
    # with a predefined repartition and device

    timeout = 60e9  # 60 sec  # TODO: should depend on the total workload, or be deactivatable
    threads_count = params.use_threading ? Threads.nthreads() : 1

    @threaded for _ in 1:threads_count
        # TODO: optimize thread block iteration to make iter-block comms faster by iterating over the first-wise edges first
        # TODO: thread block iteration should be done along the current axis

        tid = Threads.threadid()
        thread_blocks = simple_block_distribution(tid, threads_count, grid.grid_size)

        states = @view block_states[thread_blocks]
        blocks = device_block.(Ref(grid), CartesianIndices(grid.grid_size)[thread_blocks])
        reached_end_of_cycle = zeros(Bool, length(thread_blocks))

        t_start = time_ns()
        step_count = 0
        while !all(reached_end_of_cycle) && step_count < max_step_count
            if time_ns() - t_start > timeout
                blocks_pos = Tuple.(CartesianIndices(grid.grid_size)[thread_blocks])
                eoc_str = join(string.(blocks_pos) .* ": " .* string.(reached_end_of_cycle), ", ")
                solver_error(:timeout, "cycle took too long in thread $tid, blocks: $eoc_str")
            end

            for (i, (state, blk)) in enumerate(zip(states, blocks))
                reached_end_of_cycle[i] && continue
                reached_end_of_cycle[i] |= block_state_machine(params, state, blk)
            end

            step_count += 1
        end
    end

    # We cannot use checkpoints for each individual blocks, as it would require barriers.
    # Therefore we only rely on the last one at the end of a cycle as well as the time step.
    step_checkpoint(params, first(block_states), grid, "time_step")        && return true
    step_checkpoint(params, first(block_states), grid, "projection_remap") && return true
    return false
end


function solver_cycle(params::ArmonParameters, state::SolverState, data::BlockGrid)
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

    global_dt = GlobalTimeStep(params, grid)

    if params.async_cycle
        states = map(_ -> SolverState(params, global_dt), 1:prod(grid.grid_size))
    else
        state = SolverState(params, global_dt)
    end

    total_cycles_time = 0.
    t1 = time_ns()

    # Main solver loop
    while global_dt.time < maxtime && global_dt.cycle < maxcycle
        cycle_start = time_ns()

        stop = @section "solver_cycle" begin
            try
                if params.async_cycle
                    solver_cycle_async(params, states, grid)
                else
                    solver_cycle(params, state, grid)
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
