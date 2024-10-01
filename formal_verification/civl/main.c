
#include "vars.h"
#include "blocks.h"
#include "solver.h"
#include "debug.h"
#include "mpi_comms.h"

#include <omp.h>
#include <mpi.h>


int main(int argc, char** argv)
{
    int provided;
    MPI_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided));
    if (provided < MPI_THREAD_MULTIPLE) {
        fprintf(stderr, "MPI_THREAD_MULTIPLE is required, got: %d\n", provided);
        MPI_Finalize();
        return 1;
    }

    int max_threads = omp_get_max_threads();
    if (max_threads != NUM_THREADS) {
        fprintf(stderr, "Expected %d OpenMP threads, got %d\n", NUM_THREADS, max_threads);
        MPI_Finalize();
        return 2;
    }

    int comm_size;
    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &comm_size));

    MPI_Comm cart_comm;
    int cart_coords[] = { 0, 0 };
    int neighbours[] = { MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL };
    if (PROC_GRID_X == 1 && PROC_GRID_Y == 1) {
        // No MPI
        cart_comm = MPI_COMM_SELF;
    } else {
        const int dims[] = { PROC_GRID_X, PROC_GRID_Y };
        const int periods[] = { 0, 0 };
        if (MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &cart_comm) != MPI_SUCCESS) {
            fprintf(stderr, "Failed to create a %d x %d process grid with %d ranks\n", dims[0], dims[1], comm_size);
            MPI_Finalize();
            return 3;
        }

        if (cart_comm == MPI_COMM_NULL) {
            // This process is excluded from the global domain decomposition
            MPI_Finalize();
            return 0;
        }

        int cart_rank;
        MPI_CHECK(MPI_Comm_rank(cart_comm, &cart_rank));
        MPI_CHECK(MPI_Cart_coords(cart_comm, cart_rank, 2, cart_coords));

        int rank_source;
        MPI_CHECK(MPI_Cart_shift(cart_comm, 0, -1, &rank_source, &neighbours[0]));
        MPI_CHECK(MPI_Cart_shift(cart_comm, 0,  1, &rank_source, &neighbours[1]));
        MPI_CHECK(MPI_Cart_shift(cart_comm, 1, -1, &rank_source, &neighbours[2]));
        MPI_CHECK(MPI_Cart_shift(cart_comm, 1,  1, &rank_source, &neighbours[3]));
    }

    struct BlockGrid* block_grid = init_grid(cart_comm, cart_coords, neighbours, NUM_THREADS);

    bool aborted = false;
#pragma omp parallel shared(block_grid) default(none) reduction(| : aborted)
    {
        int tid = omp_get_thread_num();
        aborted |= solver_thread(tid, &block_grid->threads_workload[tid], tid == 0 ? block_grid : NULL);
    }

    if (aborted) {
        pretty_print_solver_state(block_grid);
    } else {
        check_end_state(block_grid);
    }

    free_grid(block_grid);
    MPI_Finalize();
    return 0;
}
