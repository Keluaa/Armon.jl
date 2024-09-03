
#ifndef NUM_CYCLES
#define NUM_CYCLES 3
#endif

#include "utils.pml"
#include "blocks.pml"


bool has_progress_been_made = false;
byte waiting_threads = 0;
byte active_threads = NUM_THREADS;
inline stop_busy_waiting()
{
    // We prevent fake non-progress cycles by imposing a wait when a thread cannot make any progress.
    // Without this, SPIN will find fake non-progress cycles, e.g. one thread did not start and
    // therefore the thread depending on its blocks cannot make progress.
    // While it is not implemented like such in the Julia solver, it forces SPIN to simulate a
    // realistic execution of the solver.
    atomic {
        has_progress_been_made = false;
        waiting_threads = waiting_threads + 1;
    }
    // Whenever a thread parses through its blocks and makes progress on at least one of them, all
    // threads waiting for some neighbouring block will resume work.
    // The condition on active threads is here to prevent deadlocks, which could happen if e.g. some
    // threads start waiting on each other at the same time.
    // The 'progress' label is important, as it marks this wait state as useful work and allows
    // infinite wait to happen there.
    skip;
progress_thread_wait:
    (has_progress_been_made || (active_threads == waiting_threads)) -> {
        atomic {
            waiting_threads = waiting_threads - 1;
        }
    }
}


inline notify_work_done()
{
    atomic {
        has_progress_been_made = true;
    }
}


inline solver_cycle_async(tid, current_cycle)
{
    byte i_block, idx;
    bool all_finished_cycle = false;
    bool no_progress = false;
    mtype:BlockStates prev_state;
    byte workload[MAX_WORKLOAD];
    byte num_blocks;

    get_workload(tid, num_blocks);
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
        active_threads = active_threads - 1;
    }
}


hidden byte _idx;
hidden byte _cycle_incr;
inline check_all_blocks(is_init)
{
    d_step {
        _cycle_incr = (is_init -> 0 : 1);
        for (_idx : 0 .. (TOTAL_BLOCKS - 1)) {
            assert(block_grid[_idx].sweep_num == 0);
            assert(block_grid[_idx].cycle == global_cycle + _cycle_incr);
        }
    };
}


inline check_end_state()
{
    d_step {
        assert(global_dt_state == DT_Ready);
        assert(dt_contributions == 0);
        assert(global_cycle == NUM_CYCLES);

        for (_idx : 0 .. TOTAL_BLOCKS - 1) {
            assert(block_grid[_idx].state == NewCycle);
            assert(block_grid[_idx].sweep_num == 0);
            assert(block_grid[_idx].cycle == NUM_CYCLES);
        }

        for (_idx : 0 .. TOTAL_INTERFACES - 1) {
            assert(block_interfaces[_idx].state == XCHG_NotReady);
        }
    }
}


byte join_threads = 0;
byte join_count = 0;
proctype solver_thread(byte tid)
{
    // Instead of modelling the thread creation/destruction when entering the only parallel region
    // of the solver, threads are created only once to simplify verification.
    byte current_cycle;
    do
    :: (global_cycle < NUM_CYCLES) -> {
        current_cycle = global_cycle;

        // Start of the "parallel region"
        solver_cycle_async(tid, current_cycle);

        // End of the "parallel" region: explicit barrier to join all threads together
        thread_barrier(join_threads, join_count, NUM_THREADS);

        if
        :: (tid == 0) -> {
            // Only the main thread concludes the current cycle
            printf("Completed cycle %d\n", global_cycle);
            check_all_blocks(false);
            active_threads = NUM_THREADS;
            has_progress_been_made = false;
            next_cycle();
        }
        :: else -> {
            // The other threads wait for the main thread to start the new cycle
            (global_cycle == current_cycle + 1);
        }
        fi
    }
    :: else -> break;
    od
}


init {
    byte tid;
    pid init_pid = _nr_pr;

    printf("Solving for a "); printf(GRID_SIZE_STR); printf(" grid using %d threads\n", NUM_THREADS);

    init_grid();
    check_all_blocks(true);

    for (tid : 0 .. NUM_THREADS-1) {
        run solver_thread(tid);
    }

    (_nr_pr == init_pid);  // Wait for the solver to complete (i.e. until all solver threads have terminated)
    check_end_state();
}
