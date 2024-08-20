
#ifndef NUM_THREADS
#define NUM_THREADS  4
#endif
#ifndef NUM_CYCLES
#define NUM_CYCLES   3
#endif

#include "utils.pml"
#include "blocks.pml"
#include "threads_workload.pml"


inline solver_cycle_async(tid)
{
    byte i_block, idx;
    bool all_finished_cycle = false;
    do
    :: ( all_finished_cycle) -> break;
    :: (!all_finished_cycle) -> {
        all_finished_cycle = true;
        for (i_block : 0 .. threads_workload[tid].num_blocks-1) {
            idx = threads_workload[tid].blocks[i_block];
            block_state_machine(block_grid[idx]);
            all_finished_cycle = all_finished_cycle & (block_grid[idx].state == NewCycle -> true : false);
        }
    }
    od
}


inline check_all_blocks(is_init)
{
    byte idx;
    byte cycle_incr = (is_init -> 0 : 1);
    d_step {
        for (idx : 0 .. TOTAL_BLOCKS - 1) {
            assert(block_grid[idx].sweep_num == 0);
            assert(block_grid[idx].cycle == global_cycle + cycle_incr);
        }
    };
}


byte join_threads = 0;
byte join_count = 0;
proctype solver_thread(byte tid)
{
    byte current_cycle;
    printf("Thread %d will work on %d blocks\n", tid, threads_workload[tid].num_blocks);
    do
    :: (global_cycle < NUM_CYCLES) -> {
        current_cycle = global_cycle;

        // Start of the "parallel region"
        solver_cycle_async(tid);

        // End of the "parallel" region: explicit barrier to join all threads together
        thread_barrier(join_threads, join_count, NUM_THREADS);

        if
        :: (tid == 0) -> {
            // Only the main thread concludes the current cycle
            printf("Completed cycle %d\n", global_cycle);
            check_all_blocks(false);
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

    init_grid();
    printf("Solving for a %d x %d grid using %d threads\n", GRID_SIZE_X, GRID_SIZE_Y, NUM_THREADS);
    init_threads_workload();
    check_all_blocks(true);

    for (tid : 0 .. NUM_THREADS-1) {
        run solver_thread(tid);
    }
}
