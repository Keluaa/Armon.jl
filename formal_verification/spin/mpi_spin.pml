
#ifndef MPI_SPIN_P2P_CHAN_SIZE
// Alaways start at the lowest value possible, then slowly increase it, as it adds a lot of possible
// states to explore. A value of 0 means will impose blocking communications, which may be
// incompatible with some asynchronous MPI calls.
#define MPI_SPIN_P2P_CHAN_SIZE 1
#endif

#ifndef MPI_SPIN_CHECK_MESSAGE
#define MPI_SPIN_CHECK_MESSAGE 1
#endif

#ifndef MPI_SPIN_USE_RANDOM_RECEIVE
#define MPI_SPIN_USE_RANDOM_RECEIVE 0
#endif


// TODO: multi-threading tests to check exclusive usage of a request by a single thread


mtype:MPI_Obj_Type = { MPI_Request_Type, MPI_Barrier_Data_Type, MPI_Async_Barrier_Data_Type };


typedef MPI_Request {
    mtype:MPI_Obj_Type type;  // MPI_Request_Type
#if MPI_SPIN_CHECK_MESSAGE
    byte counter;
#endif
    chan send;
    chan recv;
};


inline MPI_Request_init(request, send_chan, recv_chan)
{
    d_step {
        request.type = MPI_Request_Type;
    #if MPI_SPIN_CHECK_MESSAGE
        request.counter = 0;
    #endif
        request.send = send_chan;
        request.recv = recv_chan;
    }
}


typedef MPI_Barrier_Data {
    mtype:MPI_Obj_Type type;  // MPI_Barrier_Data_Type
    byte contribs;
    bool completed;
};
#define MPI_Barrier_Data   MPI_Barrier_Data
#define MPI_Allreduce_Data MPI_Barrier_Data
#define MPI_Reduce_Data    MPI_Barrier_Data


inline MPI_Barrier_Data_init(barrier)
{
    d_step {
        barrier.type = MPI_Barrier_Data_Type;
    }
}
inline MPI_Allreduce_Data_init(barrier) { MPI_Barrier_Data_init(barrier); }
inline MPI_Reduce_Data_init(barrier)    { MPI_Barrier_Data_init(barrier); }


typedef MPI_Async_Barrier_Data {
    mtype:MPI_Obj_Type type;  // MPI_Async_Barrier_Data_Type
    byte contribs;
    bool completed;
    bool contributions[NUM_PROC];
};
#define MPI_Ibarrier_Data   MPI_Async_Barrier_Data
#define MPI_Iallreduce_Data MPI_Async_Barrier_Data
#define MPI_Ireduce_Data    MPI_Async_Barrier_Data


inline MPI_Ibarrier_Data_init(barrier)
{
    d_step {
        barrier.type = MPI_Async_Barrier_Data_Type;
    }
}
inline MPI_Iallreduce_Data_init(barrier) { MPI_Ibarrier_Data_init(barrier); }
inline MPI_Ireduce_Data_init(barrier)    { MPI_Ibarrier_Data_init(barrier); }


inline MPI_Start(request)
{
    d_step {
        assert(request.type == MPI_Request_Type);
        assert(nfull(request.send));
#if MPI_SPIN_CHECK_MESSAGE
        request.send!request.counter;
#else
        request.send!0;
#endif
    }
}


inline MPI_Request_Receive(request)
{
    atomic {  // cannot use 'd_step' here because of the random receive
        // 'request.recv' must be non-empty, otherwise it will block
#if MPI_SPIN_CHECK_MESSAGE
        byte recv_val;
#if MPI_SPIN_USE_RANDOM_RECEIVE
        request.recv??recv_val;
#else
        request.recv?recv_val;
#endif // MPI_SPIN_USE_RANDOM_RECEIVE
        // Verify that we received the correct message
        assert(recv_val == request.counter);
        request.counter++;
#else
        request.recv?_;
#endif // MPI_SPIN_CHECK_MESSAGE
    }
}


inline MPI_Test_Request(res, request)
{
    d_step {
        assert(request.type == MPI_Request_Type);
        // TODO: is this an invalid usage of `nempty` since it is combined with channel assertions?
        res = nempty(request.recv);
        if
        :: res  -> MPI_Request_Receive(request);
        :: else -> skip;
        fi
    }
}


inline MPI_Wait_Request(res, request)
{
    atomic {
        assert(request.type == MPI_Request_Type);
        (nempty(request.recv));
        MPI_Request_Receive(request);
    }
}


inline MPI_Barrier(barrier)
{
    // The barrier is done entirely in an atomic block to make sure other processes run only when
    // we are explicitly waiting for them.
    atomic {
        assert(barrier.type == MPI_Barrier_Data_Type);
        (!barrier.completed);  // wait until any previous barrier is fully completed
        barrier.contribs++;
        bool last_contrib = barrier.contribs == NUM_PROC;
        if
        :: (last_contrib) -> {
            barrier.completed = true;
            (barrier.contribs == 1);  // wait until all other processes have left MPI_Barrier
            // Reset the barrier
            barrier.contribs = 0;
            barrier.completed = false;
        }
        :: (barrier.completed) -> barrier.contribs--;
        fi
    }
}

// Since we are discarding values, all of these are simply barriers
inline MPI_Reduce(reduction)     { MPI_Barrier(reduction); }
inline MPI_Allreduce(reduction)  { MPI_Barrier(reduction); }


inline MPI_Ibarrier(barrier)
{
    // The barrier is done entirely in a d_step block, this implies that we suppose that only
    // one process can touch the barrier data at once. This allows to reduce the possible states.
    d_step {
        assert(barrier.type == MPI_Async_Barrier_Data_Type);
        assert(!barrier.completed);  // The user is required to make sure a barrier is completed before starting it again
        assert(!barrier.contributions[RANK]);

        barrier.contribs++;
        barrier.contributions[RANK] = true;  // requires RANK to be defined in the proctype

        if
        :: (barrier.contribs == NUM_PROC) -> {
            barrier.completed = true;
            // Now each process must acknowledge that the barrier has completed and remove their
            // contribution from it. Then the barrier will be able to start again.
            barrier.contribs--;
        }
        :: else -> skip;
        fi
    }
}

inline MPI_Ireduce(reduction)    { MPI_Ibarrier(reduction); }
inline MPI_Iallreduce(reduction) { MPI_Ibarrier(reduction); }


inline MPI_reset_barrier_contrib(barrier)
{
    // Make sure our contribution is only reset once
    // This doubly important with MPI+threads since two threads may try to reset it at the same time
    if
    :: (!(barrier.contributions[RANK])) -> {
        barrier.contributions[RANK] = false;
        barrier.contribs--;
        barrier.completed = (barrier.contribs != 0);  // reset once all processes have acknowledged it
    }
    :: else -> skip;
    fi
}


inline MPI_Test_Ibarrier(res, barrier)
{
    d_step {
        assert(barrier.type == MPI_Async_Barrier_Data_Type);
        if
        :: (barrier.completed) -> {
            // The barrier is active, and has completed.
            res = true;
            MPI_reset_barrier_contrib(barrier);
        }
        :: else -> {
            // Even if 'barrier.contribs > 0', what matters here is the two things that a MPI rank 
            // can know: if it contributed, and if all processes have reached the barrier.
            // The barrier is active only if all processes have reached it or if the current process
            // reached the barrier. Otherwise, behave as if the request is inactive.
            res = !(barrier.contributions[RANK]);
        }
        fi
    }
}

inline MPI_Test_Ireduce(reduction)    { MPI_Test_Ibarrier(reduction); }
inline MPI_Test_Iallreduce(reduction) { MPI_Test_Ibarrier(reduction); }


inline MPI_Wait_Ibarrier(barrier)
{
    atomic {
        assert(barrier.type == MPI_Async_Barrier_Data_Type);
        if
        :: (barrier.completed) -> MPI_reset_barrier_contrib(barrier);
        :: (barrier.contribs > 0) -> {
            // The barrier is active only if there is at least one contribution.
            (barrier.completed);  // Wait
            // Immediatly reset our contribution
            MPI_reset_barrier_contrib(barrier);
        }
        :: else -> skip;  // Otherwise, the barrier is inactive
        fi
    }
}

inline MPI_Wait_Ireduce(reduction)    { MPI_Wait_Ibarrier(reduction); }
inline MPI_Wait_Iallreduce(reduction) { MPI_Wait_Ibarrier(reduction); }


// Similar to MPI_Test, but will fail the assertion if the barrier is active, **without** changing its state
inline MPI_Assert_Ibarrier(barrier)
{
    d_step {
        assert(barrier.type == MPI_Async_Barrier_Data_Type);
        assert(!barrier.completed && !barrier.contributions[RANK]);
    }
}

inline MPI_Assert_Ireduce(reduction)    { MPI_Assert_Ibarrier(reduction); }
inline MPI_Assert_Iallreduce(reduction) { MPI_Assert_Ibarrier(reduction); }
