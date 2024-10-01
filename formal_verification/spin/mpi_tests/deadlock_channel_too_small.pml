
#ifndef NUM_PROC
#define NUM_PROC 4
#endif

#define MPI_SPIN_STORE_MESSAGE 0


#include "../mpi_spin.pml"


chan ring = [NUM_PROC-2] of { byte, byte };


proctype mpi_proc(byte RANK)
{
    byte src_rank = (RANK > 0 -> RANK - 1 : NUM_PROC - 1);

    MPI_Request req;
    MPI_Request_init(req, ring, ring);

    // Simple synchronous MPI ring
    if
    :: (RANK < NUM_PROC-1) -> { MPI_Start(req, RANK, _); }
    :: else -> skip;
    fi
    if
    :: (RANK > 0) -> { MPI_Wait_Request(req, src_rank); }
    :: else -> skip;
    fi
}


init
{
    byte rank;
    atomic {
        for (rank : 0 .. NUM_PROC-1) {
            run mpi_proc(rank);
        }
    }
}
