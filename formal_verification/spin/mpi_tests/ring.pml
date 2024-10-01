
#ifndef NUM_PROC
#define NUM_PROC 4
#endif

#define MPI_SPIN_STORE_MESSAGE 0


#include "../mpi_spin.pml"


chan ring = [NUM_PROC*MPI_SPIN_P2P_CHAN_SIZE] of { byte, byte };


proctype mpi_proc(byte RANK)
{
    byte src_rank = (RANK > 0 -> RANK - 1 : NUM_PROC - 1);

    MPI_Request req;
    MPI_Request_init(req, ring, ring);

    // Important note: here there is no matching of the contents of the messages.
    // Some processes will receive message 'out-of-order'. For correctness some 'MPI_Barrier' should
    // be used here.

    // Simple synchronous MPI ring
    if
    :: (RANK < NUM_PROC-1) -> { MPI_Start(req, RANK, _); }
    :: else -> skip;
    fi
    if
    :: (RANK > 0) -> { MPI_Wait_Request(req, src_rank); }
    :: else -> skip;
    fi

    // Simple asynchronous MPI ring, repeated twice.
    // This works since the channel can contain at least NUM_PROC messages, which would be the maximum
    // number of active MPI request active at once.
    byte is_done;
    byte loops = 2;
    do
    :: (loops > 0) -> {
        MPI_Start(req, RANK, _);
        MPI_Test_Request(is_done, req, src_rank);
        if
        :: (!is_done) -> { MPI_Wait_Request(req, src_rank); }
        :: else -> skip;
        fi
        loops--;
    }
    :: (loops == 0) -> break;
    od
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
