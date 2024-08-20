#ifndef ARMON_CIVL_UTILS_H
#define ARMON_CIVL_UTILS_H

#ifdef _CIVL

#include <civlc.cvh>

// Multiline macros messes with CIVL's parser, so everything is in a signle line :(

#define CIVL_atomic_cas(success, val, expected_value, new_value) $atomic { if ((val) == (expected_value)) { (val) = (new_value); (success) = true; } else { (success) = false; } }
#define CIVL_atomic_incr_fetch(val, res)     $atomic { (val)++;        (res) = (val); }
#define CIVL_atomic_add_fetch(val, arg, res) $atomic { (val) += (arg); (res) = (val); }
#define CIVL_atomic_store(val, arg)          (val) = (arg)
#define CIVL_atomic_load(val, res)           (res) = (val)

#define AtomicVar

#define CIVL_assert(expr) $assert(expr)
#define CIVL_assume(expr) $assume(expr)
#define CIVL_pure         $pure

typedef unsigned short byte;  // 'char' is not an integer in CIVL

#else

#include <stdatomic.h>

#define CIVL_atomic_cas(success, val, expected_value, new_value) \
    (success) = atomic_compare_exchange_strong(&(val), &(expected_value), (new_value))

#define CIVL_atomic_incr_fetch(val, res) \
    (res) = (atomic_fetch_add(&(val), 1) + 1)

#define CIVL_atomic_add_fetch(val, arg, res) \
    (res) = (atomic_fetch_add(&(val), (arg)) + (arg))

#define CIVL_atomic_store(val, arg) \
    atomic_store(&(val), (arg))

#define CIVL_atomic_load(val, res) \
    (res) = atomic_load(&(val))

#define AtomicVar _Atomic

#define CIVL_assert(expr) if (!(expr)) assert_fail(#expr, __FILE_NAME__, __LINE__)
#define CIVL_assume(expr) while (false)
#define CIVL_pure

void assert_fail(const char* test, const char* file, int line);

typedef unsigned char byte;

#endif //_CIVL

#include <stdbool.h>
#include <stdio.h>

void thread_barrier(AtomicVar int* join_counter, AtomicVar int* barrier_counter, int num_threads);

#endif //ARMON_CIVL_UTILS_H
