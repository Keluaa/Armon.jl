
#ifndef NUM_PROC
#define NUM_PROC 4
#endif


#include "../mpi_spin.pml"


MPI_Barrier_Data barrier;
MPI_Ibarrier_Data ibarrier;


proctype mpi_proc(byte RANK)
{
    byte is_done;

    MPI_Barrier(barrier);

    // Two testing paths for diversity
    if
    :: (RANK % 2 == 0) -> {
        MPI_Assert_Ibarrier(ibarrier);
        MPI_Ibarrier(ibarrier);
    }
    :: (RANK % 2 == 1) -> {
        // Processes which did not yet start the barrier consider it as inactive
        MPI_Assert_Ibarrier(ibarrier);
        MPI_Test_Ibarrier(is_done, ibarrier);
        assert(is_done);
        MPI_Ibarrier(ibarrier);
    }
    fi

    // There is two ways to complete a barrier with MPI, both are equivalent, in that if the current
    // process knows that all other ranks have contributed, it now considers the barrier as completed
    // and a new one can be started.
    MPI_Wait_Ibarrier(ibarrier);
    MPI_Test_Ibarrier(is_done, ibarrier);
    assert(is_done);
    MPI_Assert_Ibarrier(ibarrier);

    // Barriers should be reusable
    MPI_Ibarrier(ibarrier);
    MPI_Wait_Ibarrier(ibarrier);
    MPI_Assert_Ibarrier(ibarrier);

    MPI_Barrier(barrier);
    MPI_Barrier(barrier);
}


init
{
    byte rank;

    MPI_Barrier_Data_init(barrier);
    MPI_Ibarrier_Data_init(ibarrier);

    atomic {
        for (rank : 0 .. NUM_PROC-1) {
            run mpi_proc(rank);
        }
    }
}
