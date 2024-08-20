
#include "vars.h"
#include "blocks.h"
#include "solver.h"
#include "debug.h"

#include <omp.h>


int main()
{
    int max_threads = omp_get_max_threads();
    if (max_threads != NUM_THREADS) {
        fprintf(stderr, "Expected %d OpenMP threads, got %d\n", NUM_THREADS, max_threads);
        return 1;
    }

    struct BlockGrid* block_grid = init_grid(NUM_THREADS);

    bool aborted = false;
#pragma omp parallel shared(block_grid) default(none) reduction(| : aborted)
    {
        int tid = omp_get_thread_num();
        aborted |= solver_thread(tid, &block_grid->threads_workload[tid], tid == 0 ? block_grid : NULL);
    }

    if (aborted) {
        pretty_print_solver_state(block_grid);
    }

    free_grid(block_grid);
    return 0;
}
