
#include "solver.h"
#include "time_step.h"
#include "debug.h"

bool solver_cycle_async(struct ThreadWorkload* thread_workload)
{
    int not_finished_cycle = 1;
    int step_count = 0;
    while (not_finished_cycle != 0 && step_count < 1000) {  // TODO: remove the max step count???
        not_finished_cycle = 0;
        step_count++;
        for (int block_i = 0; block_i < thread_workload->num_blocks; block_i++) {
            struct Block* block = thread_workload->threads_blocks[block_i];
            block_state_machine(block);
            not_finished_cycle += (block->state == NewCycle && block->cycle == global_cycle + 1) ? 0 : 1;
        }
    }
    return not_finished_cycle != 0;
}

AtomicVar int solver_join_count = 0;
AtomicVar int solver_barrier_counter = 0;
AtomicVar int solver_abort = 0;
bool solver_thread(int tid, struct ThreadWorkload* thread_workload, struct BlockGrid* block_grid)
{
    for (int current_cycle = 0; current_cycle < NUM_CYCLES;) {
        current_cycle = global_cycle;

        // Start of the "parallel" region
        bool abort = solver_cycle_async(thread_workload);

        CIVL_atomic_add_fetch(solver_abort, abort, abort);
        // End of the "parallel" region: explicit barrier to join all threads
        thread_barrier(&solver_join_count, &solver_barrier_counter, NUM_THREADS);

        if (solver_abort > 0) {
#ifndef _CIVL
            if (tid == 0) printf("Solver abort\n");
#endif
            return true;
        } else if (tid == 0) {
            // Only the main thread concludes the current cycle
#ifndef _CIVL
            printf("Completed cycle %d\n", global_cycle);
#endif
//            pretty_print_solver_state(block_grid);
            next_cycle();
        } else {
            // The other threads wait for the main thread to start the new cycle
        }

        thread_barrier(&solver_join_count, &solver_barrier_counter, NUM_THREADS);
    }

    return false;
}
