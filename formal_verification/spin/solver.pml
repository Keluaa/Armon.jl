
bool has_progress_been_made[NUM_PROC];
byte waiting_threads[NUM_PROC];
byte active_threads[NUM_PROC];
inline stop_busy_waiting()
{
    // We prevent fake non-progress cycles by imposing a wait when a thread cannot make any progress.
    // Without this, SPIN will find fake non-progress cycles, e.g. one thread did not start and
    // therefore the thread depending on its blocks cannot make progress.
    // While it is not implemented like such in the Julia solver, it forces SPIN to simulate a
    // realistic execution of the solver.
    atomic {
        has_progress_been_made[RANK] = false;
        waiting_threads[RANK] = waiting_threads[RANK] + 1;
    }
    // Whenever a thread parses through its blocks and makes progress on at least one of them, all
    // threads waiting for some neighbouring block will resume work.
    // The condition on active threads is here to prevent deadlocks, which could happen if e.g. some
    // threads start waiting on each other at the same time.
    // The 'progress' label is important, as it marks this wait state as useful work and allows
    // infinite wait to happen there.
    skip;
    // TODO: wait on active MPI requests?
progress_thread_wait:  // TODO: this cannot be right...
    (has_progress_been_made[RANK] || (active_threads[RANK] == waiting_threads[RANK])) -> {
        atomic {
            waiting_threads[RANK] = waiting_threads[RANK] - 1;
        }
    }
}


inline notify_work_done()
{
    atomic {
        has_progress_been_made[RANK] = true;
    }
}


inline solver_cycle_async(current_cycle)
{
    byte i_block, idx;
    bool all_finished_cycle = false;
    bool no_progress = false;
    mtype:BlockStates prev_state;
    byte workload[MAX_WORKLOAD];
    byte num_blocks;

    get_workload(num_blocks);
    assert(0 <= num_blocks <= MAX_WORKLOAD);

    do
    :: ( all_finished_cycle) -> break;
    :: (!all_finished_cycle) -> {
        all_finished_cycle = true;
        no_progress = true;
        for (i_block : 0 .. (num_blocks-1)) {
            idx = workload[i_block];
            if
            :: (idx == NULL_BLOCK) -> skip;
            :: else -> {
                prev_state = block_grid[idx].state;
                block_state_machine(block_grid[idx], idx);
                all_finished_cycle = all_finished_cycle & (block_grid[idx].state == NewCycle   && block_grid[idx].cycle == current_cycle + 1);
                no_progress        = no_progress        & (block_grid[idx].state == prev_state && block_grid[idx].cycle == current_cycle    );
            }
            fi
        }

        // TODO: remove or not?
        // Guarantee that progress always being made, otherwise wait for other threads
        if
        :: (no_progress && num_blocks > 0) -> stop_busy_waiting();
        :: else ->
progress_thread_work:  // Mirror of 'progress_thread_wait': this ensures there is a progress path in all threads
                notify_work_done();
        fi
    }
    od

    atomic {
        active_threads[RANK] = active_threads[RANK] - 1;
    }
}


hidden byte _idx;
hidden byte _count_to_check;
inline check_rank_state(_cycle)
{
    d_step {
        assert(SOLVER_STATE.dt_state == DT_Ready);
        assert(SOLVER_STATE.dt_contributions == 0);
        assert(SOLVER_STATE.cycle == _cycle);

        get_all_elements(_count_to_check, 0);  // fills '_all_idx' with the block indexes of the current MPI rank
        for (_idx : 0 .. _count_to_check - 1) {
            assert(block_grid[_all_idx[_idx]].state == NewCycle);
            assert(block_grid[_all_idx[_idx]].sweep_num == 0);
            assert(block_grid[_all_idx[_idx]].cycle == _cycle);
        }

        get_all_elements(_count_to_check, 1);  // fills '_all_idx' with the interface indexes of the current MPI rank
        for (_idx : 0 .. _count_to_check - 1) {
            assert(block_interfaces[_all_idx[_idx]].state == XCHG_NotReady);
        }

#if USE_MPI
        MPI_Assert_Iallreduce(global_time_step_reduction);
#if MPI_SPIN_CHECK_MESSAGE
        get_all_elements(_count_to_check, 2);  // fills '_all_idx' with the remote block indexes of the current MPI rank
        for (_idx : 0 .. _count_to_check - 1) {
            assert(remote_blocks[_all_idx[_idx]].req.counter == _cycle * MAX_SWEEPS);
        }
#endif // MPI_SPIN_CHECK_MESSAGE
#endif // USE_MPI
    }
}


inline check_all_states(_cycle)
{
    d_step {
        for (_idx : 0 .. NUM_PROC - 1) {
            assert(solver_states[_idx].dt_state == DT_Ready);
            assert(solver_states[_idx].dt_contributions == 0);
            assert(solver_states[_idx].cycle == _cycle);
        }

        for (_idx : 0 .. TOTAL_BLOCKS - 1) {
            assert(block_grid[_idx].state == NewCycle);
            assert(block_grid[_idx].sweep_num == 0);
            assert(block_grid[_idx].cycle == _cycle);
        }

        for (_idx : 0 .. TOTAL_INTERFACES - 1) {
            assert(block_interfaces[_idx].state == XCHG_NotReady);
        }

#if USE_MPI && MPI_SPIN_CHECK_MESSAGE
        for (_idx : 0 .. TOTAL_REMOTE_BLOCKS - 1) {
            assert(remote_blocks[_idx].req.counter == _cycle * MAX_SWEEPS);
        }
#endif
    }
}


byte join_threads[NUM_PROC];
byte join_count[NUM_PROC];
inline solver_thread()
{
    // Instead of modelling the thread creation/destruction when entering the only parallel region
    // of the solver, threads are created only once to simplify verification.
    byte current_cycle;
    do
    :: (SOLVER_STATE.cycle < NUM_CYCLES) -> {
        current_cycle = SOLVER_STATE.cycle;

        // Start of the "parallel region"
        solver_cycle_async(current_cycle);

        // End of the "parallel" region: explicit barrier to join all threads together
        thread_barrier(join_threads[RANK], join_count[RANK], NUM_THREADS);

        if
        :: (TID == 0) -> {
            // Only the main thread concludes the current cycle
            printf("Completed cycle %d\n", SOLVER_STATE.cycle);
            check_rank_state(current_cycle + 1);
            active_threads[RANK] = NUM_THREADS;
            has_progress_been_made[RANK] = false;
            next_cycle();
        }
        :: else -> {
            // The other threads wait for the main thread to start the new cycle
            (SOLVER_STATE.cycle == current_cycle + 1);
        }
        fi
    }
    :: else -> break;
    od
}


hidden byte _init_rank;
inline init_global_vars()
{
    d_step {
        for (_init_rank : 0 .. NUM_PROC-1) {
            // initialize all global variables for the rank
            has_progress_been_made[_init_rank] = false;
            waiting_threads[_init_rank] = 0;
            active_threads[_init_rank] = NUM_THREADS;

            join_threads[_init_rank] = 0;
            join_count[_init_rank] = 0;

            solver_states[_init_rank].dt_state = DT_Ready;
            solver_states[_init_rank].dt_contributions = 0;
            solver_states[_init_rank].cycle = 0;
        }
    }
}
