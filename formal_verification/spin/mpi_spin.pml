
/**
 * A basic MPI specification in SPIN.
 *
 * To reduce the complexity of the model, some optimizations were done:
 *  - Collectives do not communicate data. This means that a reduction is equivalent to a barrier.
 *    Therefore they can be implemented using bounded atomic values.
 *  - P2P communications only send a tag and some dummy value. Channels should therefore be declared
 *    as 'of { byte, byte }'. If 'MPI_SPIN_STORE_MESSAGE == 1', then the dummy value is a counter
 *    given by the user whose value is unique at each communication cycle, allowing to check that we
 *    received the correct message. All receptions are done randomly, no ordering of messages is
 *    guarenteed by this implementation: it is up to the user to do so, and assertions about received
 *    values are done using the field `MPI_Request.val`.
 *
 * Collective operations must use a global variable of type 'MPI_<collective>_Data' which stores
 * data used to implement it. It should be initialized once by 'MPI_<collective>_Data_init'.
 *
 * P2P operations must use a local variable of type 'MPI_Request' initialized with 'MPI_Request_init',
 * which take as parameters two channels, one for sending operations, another for reception. This
 * pattern allow to use channel assertions (e.g. two processes use two channels, messages go in a
 * single direction in each channel, 'xr/xs' assertions can therefore be used by each process).
 * P2P operations use a tag to match messages together. Any message with that tag in the channel can
 * be received: all receive operations are done with the '??' random receive operator. This will
 * receive messages out-of-order if no precautions are taken: while the MPI specification guarentees
 * *some* ordering of messages, it is preferrable to make sure that it is the application *itself*
 * which imposes the ordering, especially in non-trivial MPI+Thread usages.
 *
 * Some values need to be manually defined in order to use MPI calls:
 *  - The macro 'NUM_PROC' should be set to the number of MPI processes. It must be defined before
 *    including this file.
 *  - The local variable 'RANK' is supposed to be defined in each 'proctype' using MPI. It should be
 *    unique and between 0 to NUM_PROC-1.
 */

#ifndef MPI_SPIN_P2P_CHAN_SIZE
// Alaways start at the lowest value possible, then slowly increase it, as it adds a lot of possible
// states to explore. A value of 0 means will impose blocking communications, which may be
// incompatible with some asynchronous MPI calls.
#define MPI_SPIN_P2P_CHAN_SIZE 1
#endif

#ifndef MPI_SPIN_STORE_MESSAGE
#define MPI_SPIN_STORE_MESSAGE 1
#endif


// TODO: multi-threading tests to check exclusive usage of a request by a single thread


mtype:MPI_Obj_Type = { MPI_Request_Type, MPI_Barrier_Data_Type, MPI_Async_Barrier_Data_Type };


typedef MPI_Request {
    mtype:MPI_Obj_Type type;  // MPI_Request_Type
    chan send;
    chan recv;
#if MPI_SPIN_STORE_MESSAGE
    byte val;
#endif
};


inline MPI_Request_init(request, send_chan, recv_chan)
{
    d_step {
        request.type = MPI_Request_Type;
        request.send = send_chan;
        request.recv = recv_chan;
#if MPI_SPIN_STORE_MESSAGE
        request.val = 0;
#endif
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


typedef MPI_Contributions_Side {
    bool s[2];
}

typedef MPI_Async_Barrier_Data {
    mtype:MPI_Obj_Type type;  // MPI_Async_Barrier_Data_Type
    // All fields are doubled: depending on '.side[RANK]', we access either side. This allows one
    // rank to start a new barrier while the previous one is still going.
    byte contribs[2];
    bool completed[2];
    MPI_Contributions_Side contributions[NUM_PROC];
    bit side[NUM_PROC];
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


inline MPI_Start(request, tag, val)  // TODO: rename to Isend, 'val' can be '_' if 'MPI_SPIN_STORE_MESSAGE == 0'
{
    d_step {
        assert(request.type == MPI_Request_Type);
        assert(nfull(request.send));
#if MPI_SPIN_STORE_MESSAGE
        request.send!tag,val;
#else
        request.send!tag,0;
#endif
    }
}


inline MPI_Request_Receive(request, tag)
{
    atomic {  // cannot use 'd_step' here because of the random receive
        // 'request.recv' must be non-empty, otherwise it will block
#if MPI_SPIN_STORE_MESSAGE
        request.recv??eval(tag),request.val;
#else
        request.recv??eval(tag),_;
#endif
    }
}


inline MPI_Test_Request(res, request, tag)
{
    d_step {
        assert(request.type == MPI_Request_Type);
        res = request.recv??[eval(tag),_];
        if
        :: res  -> MPI_Request_Receive(request, tag);
        :: else -> skip;
        fi
    }
}


inline MPI_Wait_Request(request, tag)  // TODO: add an alias for Recv
{
    atomic {
        assert(request.type == MPI_Request_Type);
        (request.recv??[eval(tag),_]);
        MPI_Request_Receive(request, tag);
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
        if
        :: (barrier.contribs == NUM_PROC) -> {
            barrier.completed = true;
            (barrier.contribs == 1);  // wait until all other processes have left MPI_Barrier
            // Reset the barrier
            barrier.contribs = 0;
            barrier.completed = false;
        }
        :: else -> {
            (barrier.completed);
            barrier.contribs--;
        }
        fi
    }
}

// Since we are discarding values, all of these are simply barriers
inline MPI_Reduce(reduction)     { MPI_Barrier(reduction); }
inline MPI_Allreduce(reduction)  { MPI_Barrier(reduction); }


inline MPI_reset_barrier_contrib(barrier)
{
    // Make sure our contribution is only reset once
    // This doubly important with MPI+threads since two threads may try to reset it at the same time
    if
    :: (barrier.contributions[RANK].s[barrier.side[RANK]]) -> {
        barrier.contributions[RANK].s[barrier.side[RANK]] = false;
        barrier.contribs[barrier.side[RANK]]--;

        // reset once all processes have acknowledged it
        barrier.completed[barrier.side[RANK]] = (barrier.contribs[barrier.side[RANK]] != 0);

        // Switch side everytime we reset our contribution. It allows this rank to start another
        // barrier while a previous one is still going.
        barrier.side[RANK] = barrier.side[RANK] ^ 1;  // 0 to 1, 1 to 0
    }
    :: else -> skip;
    fi
}


inline MPI_Ibarrier(barrier)
{
    // The barrier is done entirely in a d_step block, this implies that we suppose that only
    // one process can touch the barrier data at once. This allows to reduce the possible states.
    d_step {
        assert(barrier.type == MPI_Async_Barrier_Data_Type);
        assert(!barrier.contributions[RANK].s[barrier.side[RANK]]);
        assert(!barrier.completed[barrier.side[RANK]]);

        barrier.contribs[barrier.side[RANK]]++;
        barrier.contributions[RANK].s[barrier.side[RANK]] = true;

        if
        :: (barrier.contribs[barrier.side[RANK]] == NUM_PROC) -> {
            barrier.completed[barrier.side[RANK]] = true;
            // Now each process must acknowledge that the barrier has completed and remove their
            // contribution from it. Then the barrier will be able to start again.
            MPI_reset_barrier_contrib(barrier);
        }
        :: else -> skip;
        fi
    }
}

inline MPI_Ireduce(reduction)    { MPI_Ibarrier(reduction); }
inline MPI_Iallreduce(reduction) { MPI_Ibarrier(reduction); }


inline MPI_Test_Ibarrier(res, barrier)
{
    d_step {
        assert(barrier.type == MPI_Async_Barrier_Data_Type);
        if
        :: (barrier.completed[barrier.side[RANK]]) -> {
            // The barrier is active, and has completed.
            res = true;
            MPI_reset_barrier_contrib(barrier);
        }
        :: else -> {
            // Even if 'barrier.contribs > 0', what matters here is the two things that a MPI rank 
            // can know: if it contributed, and if all processes have reached the barrier.
            // The barrier is active only if all processes have reached it or if the current process
            // reached the barrier. Otherwise, behave as if the request is inactive.
            res = !(barrier.contributions[RANK].s[barrier.side[RANK]]);
        }
        fi
    }
}

inline MPI_Test_Ireduce(res, reduction)    { MPI_Test_Ibarrier(res, reduction); }
inline MPI_Test_Iallreduce(res, reduction) { MPI_Test_Ibarrier(res, reduction); }


inline MPI_Wait_Ibarrier(barrier)
{
    atomic {
        assert(barrier.type == MPI_Async_Barrier_Data_Type);
        if
        :: (barrier.completed[barrier.side[RANK]]) -> MPI_reset_barrier_contrib(barrier);
        :: (!barrier.completed[barrier.side[RANK]] && barrier.contributions[RANK].s[barrier.side[RANK]]) -> {
            // The barrier is active only if the current process contributed to it.
            (barrier.completed[barrier.side[RANK]]);  // Wait
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
        assert(!barrier.contributions[RANK].s[barrier.side[RANK]]);
    }
}

inline MPI_Assert_Ireduce(reduction)    { MPI_Assert_Ibarrier(reduction); }
inline MPI_Assert_Iallreduce(reduction) { MPI_Assert_Ibarrier(reduction); }
