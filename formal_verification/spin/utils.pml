
inline atomic_cas(success, val, expected_value, new_value)
{
    atomic {
        if
        :: (val == expected_value) -> { val = new_value; success = true; }
        :: else -> { success = false; }
        fi
    }
}


inline thread_barrier(global_counter, barrier_counter, threads_count)
{
    byte local_counter, old_barrier_counter;
    atomic {
        old_barrier_counter = barrier_counter;
        global_counter++;
        local_counter = global_counter;
    }
    // The remaining threads wait for the last thread to increment the barrier
    if
    :: (local_counter == threads_count) -> {
        global_counter = 0;
        barrier_counter++;
    }
    :: (old_barrier_counter < barrier_counter) -> skip;
    fi
}
