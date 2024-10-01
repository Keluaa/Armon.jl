
#include "vars.h"
#include "blocks.h"
#include "solver.h"
#include "debug.h"
#include "mpi_comms.h"

#if NUM_THREADS != 1
#error "Sequential run requires NUM_THREADS=1"
#endif


int main(int argc, char** argv)
{
    int cart_coords[] = { 0, 0 };
    int neighbours[] = { MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL };
    struct BlockGrid* block_grid = init_grid(MPI_COMM_WORLD, cart_coords, neighbours, NUM_THREADS);

    int tid = 0;
    bool aborted = solver_thread(tid, &block_grid->threads_workload[tid], block_grid);

    if (aborted) {
        pretty_print_solver_state(block_grid);
    }

    free_grid(block_grid);
    return 0;
}
