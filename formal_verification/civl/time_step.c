
#include "time_step.h"

void update_dt()
{
    enum TimeStepState prev_dt_state = global_dt_state;

    if (prev_dt_state == DT_AllContributed) {
#if USE_MPI
        global_dt_state = DT_DoingMPI;
        return;
#endif
    } else if (prev_dt_state == DT_WaitingForMPI) {
        // receive the MPI reduction result
        // TODO
//        assert(nempty(subdomain_neighbours_dt_reduction));
//        subdomain_neighbours_dt_reduction?_;
    } else {
        // TODO: another thread went here first? is this normal? this isn't in the Julia implementation!
        return;
    }

    dt_contributions = 0;
    global_dt_state = DT_Done;
}

enum TimeStepState wait_for_dt()
{
    bool cas_ok = false;
    enum TimeStepState old_val = DT_DoingMPI, new_val = DT_WaitingForMPI;
    CIVL_atomic_cas(cas_ok, global_dt_state, old_val, new_val);

    if (cas_ok) {
#if USE_MPI
        // TODO: wait until `subdomain_neighbours_dt_reduction` has a value + make sure only one thread waits on the request
#endif
        update_dt();
        return global_dt_state;  // TODO: is this correct? shouldn't we use a "return value" from `update_dt` instead?
    } else {
        return DT_WaitingForMPI;
    }
}

void contribute_to_dt(struct Block* block)
{
    byte current_contributions = 0;
    CIVL_atomic_incr_fetch(dt_contributions, current_contributions);

    if (current_contributions == TOTAL_BLOCKS) {
        bool cas_ok = false;
        enum TimeStepState old_val = DT_Ready, new_val = DT_AllContributed;
        CIVL_atomic_cas(cas_ok, global_dt_state, old_val, new_val);
        if (cas_ok) update_dt();
    }
}

void next_cycle()
{
    enum TimeStepState current_dt_state;

#if DO_TIME_STEP
    current_dt_state = global_dt_state;
    retry_next_dt:
    CIVL_assert(current_dt_state == DT_DoingMPI || current_dt_state == DT_Done);
    if (current_dt_state == DT_DoingMPI) {
        current_dt_state = wait_for_dt();
        goto retry_next_dt;
    } else if (current_dt_state == DT_Done) {
        global_dt_state = DT_Ready;
    }
#endif

    global_cycle++;
}

void next_time_step(struct Block* block)
{
#if DO_TIME_STEP
    enum TimeStepState dt_state = global_dt_state;
    retry_time_step:
    CIVL_assert(dt_state == DT_DoingMPI || dt_state == DT_Ready || dt_state == DT_Done);
    if (dt_state == DT_DoingMPI) {
        dt_state = wait_for_dt();
        CIVL_assert(dt_state != DT_DoingMPI);
        goto retry_time_step;
    } else if (dt_state == DT_Ready) {
        // local_time_step
        contribute_to_dt(block);
        block->must_wait = global_cycle == 0;  // The first cycle requires the time step before continuing
    } else if (dt_state == DT_Done) {
        block->must_wait = false;
    }
#else
    block->must_wait = false;
#endif
}
