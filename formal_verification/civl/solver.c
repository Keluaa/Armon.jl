
#include "solver.h"
#include "time_step.h"
#include "debug.h"

AtomicVar int threads_progress = 0;
AtomicVar int waiting_threads = 0;
AtomicVar int active_threads = NUM_THREADS;

void stop_busy_waiting(struct ThreadWorkload* thread_workload)
{
#ifdef _CIVL
    // TODO: but with MPI threads are unaware of the progress of communications!
    // To make sure that CIVL does not get lost in repeats of non-progressing cycle steps, threads will notify others
    // that they did some work, implying that at least one block exchange/time step reduction has been completed.
    // Threads that must wait for others will do it here.
    int old_progress = threads_progress;
    waiting_threads += 1;
    $when ((threads_progress > old_progress) || (active_threads == waiting_threads)) waiting_threads -= 1;
#else
    // TODO: MPI busy waiting
    for (int block_i = 0; block_i < thread_workload->num_blocks; block_i++) {
        struct Block* block = thread_workload->threads_blocks[block_i];
//        block.
    }
#endif
}

void notify_work_done(struct ThreadWorkload* thread_workload)
{
#ifdef _CIVL
    // Notify other threads that some work has been done
    int _;
    CIVL_atomic_incr_fetch(threads_progress, _);
#endif
}

bool solver_cycle_async(struct ThreadWorkload* thread_workload)
{
    int not_finished_cycle = 1;
    int no_progress;
    while (not_finished_cycle != 0) {  // TODO: remove the max step count???
        not_finished_cycle = 0;
        no_progress = 0;

        for (int block_i = 0; block_i < thread_workload->num_blocks; block_i++) {
            struct Block* block = thread_workload->threads_blocks[block_i];
            enum BlockState prev_state = block->state;
            block_state_machine(block);
            not_finished_cycle += (block->state == NewCycle && block->cycle == global_cycle + 1) ? 0 : 1;
            no_progress += (block->state == prev_state && block->cycle == global_cycle) ? 0 : 1;
        }

        // TODO: new way of waiting for progress, is it implementable in Julia?
        if (no_progress == thread_workload->num_blocks && thread_workload->num_blocks > 0) {
            // No threads did any progress this step: slow neighbouring thread, MPI communication in progress or deadlock
            // are possible explanations.
            stop_busy_waiting(thread_workload);  // In the Julia solver, this would be called once every Nth step without progress
        } else {
            notify_work_done(thread_workload);
        }
    }
    active_threads -= 1;
    return false;
}

AtomicVar int solver_join_count = 0;
AtomicVar int solver_barrier_counter = 0;
AtomicVar int solver_abort = 0;
bool solver_thread(int tid, struct ThreadWorkload* thread_workload, struct BlockGrid* block_grid)
{
    for (int current_cycle = 0; current_cycle < NUM_CYCLES; current_cycle = global_cycle) {
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
            CIVL_assert(active_threads == 0);
            active_threads = NUM_THREADS;
        } else {
            // The other threads wait for the main thread to start the new cycle
        }

        thread_barrier(&solver_join_count, &solver_barrier_counter, NUM_THREADS);
        CIVL_assert(global_cycle > current_cycle);
    }

    return false;
}
