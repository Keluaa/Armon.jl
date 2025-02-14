
disp_blk(blk, var; on_device=true) = reshape(getfield(block_data(blk; on_device), var), block_size(blk))'
disp_real_blk(blk, var; on_device=true) = view(disp_blk(blk, var; on_device)', (.+).(Base.oneto.(real_block_size(blk.size)), ghosts(blk.size))...)'
disp_grid_state(grid) = permutedims(reshape(solver_step.(solver_state.(all_blocks(grid))), grid.grid_size))
disp_mirror_y(A) = view(A, size(A, 1):-1:1, :)  # places the bottom-left cell at the bottom-left of the display


function debug_grid_state(io::IO, params::ArmonParameters, grid::BlockGrid; legend=true)
    (; global_dt) = grid

    println(io, "State of sub-domain ", params.cart_coords, " (rank $(params.rank)):")
    println(io, " - Cycle: ", global_dt.cycle)
    println(io, " - Time step state: ", (@atomic global_dt.state.x), " (lock: ", (@atomic global_dt.state_lock.x), ")")
    println(io, " - Time step contributions: ", (@atomic global_dt.contributions.x), " / ", global_dt.expected_count)

    println(io, "Thread assignments:")
    tw_grid = thread_workload_to_grid(grid.grid_size, grid.threads_workload)
    show(io, MIME"text/plain"(), tw_grid)
    println(io)

    println(io, "Block states of the grid (state, cycle, axis):")
    states = map(all_blocks(grid)) do blk
        blk_state = solver_state(blk)
        return Int(solver_step(blk_state)), blk_state.cycle, Int(blk_state.axis)
    end
    states = permutedims(reshape(states, grid.grid_size))
    show(io, MIME"text/plain"(), states)
    println(io)

    if legend
        println(io, "Legend:")
        for s in instances(SolverStep.T)
            s == SolverStep.ErrorState && continue
            println(io, " - $(Int(s)): $(Symbol(s))")
        end
    end
end
