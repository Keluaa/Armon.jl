
#include "utils.h"

#ifndef _CIVL
void assert_fail(const char* test, const char* file, int line)
{
    fprintf(stderr, "ASSERT FAILED at %s line %d: %s\n", file, line, test);
}
#endif

void thread_barrier(AtomicVar int* join_counter, AtomicVar int* barrier_counter, int num_threads)
{
    int old_counter, join_count;
    CIVL_atomic_load(*barrier_counter, old_counter);
    CIVL_atomic_incr_fetch(*join_counter, join_count);  // join_count = ++join_counter

    if (join_count == num_threads) {
        // All threads reached the barrier
        *join_counter = 0;
        *barrier_counter += 1;
    } else {
        // Wait until the last thread resets the barrier
#ifdef _CIVL
        $when (old_counter != *barrier_counter);
#else
        while (old_counter == *barrier_counter) {}  // ugly busy wait
#endif
    }
}
