
#include "grid_definition.h"
#include "utils.pml"
#include "mpi_spin.pml"
#include "block_grid.pml"
#include "blocks.pml"
#include "solver.pml"
#include "grid_definition.pml"

init {
    pid init_pid = _nr_pr;

#if USE_MPI
    printf("Solving using "); printf(PROC_GRID_STR); printf(" process grid for a ");
    printf(GRID_SIZE_STR); printf(" local grid using %d x %d threads\n", NUM_PROC, NUM_THREADS);
#else
    printf("Solving for a "); printf(GRID_SIZE_STR); printf(" grid using %d threads\n", NUM_THREADS);
#endif

    init_grid();
    init_global_vars();
#if USE_MPI
    MPI_Iallreduce_Data_init(global_time_step_reduction);
#endif
    check_all_states(0);

    run_all_procs();  // Launch all threads of all MPI processes

    (_nr_pr == init_pid);  // Wait for the solver to complete (i.e. until all solver threads have terminated)
    check_all_states(NUM_CYCLES);
}
