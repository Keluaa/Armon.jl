#ifndef ARMON_CIVL_VARS_H
#define ARMON_CIVL_VARS_H

#ifdef _CIVL

int get_num_threads();
int get_num_cycles();
int get_max_sweeps();
int get_grid_size_x();
int get_grid_size_y();
int get_proc_grid_x();
int get_proc_grid_y();
int get_total_blocks();
int get_total_remote_blocks();
int get_total_interfaces();

#define NUM_THREADS          get_num_threads()
#define NUM_CYCLES           get_num_cycles()
#define MAX_SWEEPS           get_max_sweeps()
#define GRID_SIZE_X          get_grid_size_x()
#define GRID_SIZE_Y          get_grid_size_y()
#define PROC_GRID_X          get_proc_grid_x()
#define PROC_GRID_Y          get_proc_grid_y()
#define TOTAL_BLOCKS         get_total_blocks()
#define TOTAL_REMOTE_BLOCKS  get_total_remote_blocks()
#define TOTAL_INTERFACES     get_total_interfaces()

#else

#ifndef NUM_THREADS
#define NUM_THREADS   4
#endif

#ifndef NUM_CYCLES
#define NUM_CYCLES    3
#endif

#ifndef MAX_SWEEPS
#define MAX_SWEEPS    2
#endif

#ifndef GRID_SIZE_X
#define GRID_SIZE_X   3
#endif
#ifndef GRID_SIZE_Y
#define GRID_SIZE_Y   3
#endif

#define TOTAL_BLOCKS         (GRID_SIZE_X * GRID_SIZE_Y)
#define TOTAL_REMOTE_BLOCKS  (2 * (GRID_SIZE_X + GRID_SIZE_Y))
#define TOTAL_INTERFACES     (2 * TOTAL_BLOCKS - GRID_SIZE_X - GRID_SIZE_Y)

#ifndef PROC_GRID_X
#define PROC_GRID_X 2
#endif

#ifndef PROC_GRID_Y
#define PROC_GRID_Y 2
#endif

#endif //_CIVL

#ifndef DO_TIME_STEP
#define DO_TIME_STEP      1
#endif
#ifndef DO_HALO_EXCHANGE
#define DO_HALO_EXCHANGE  1
#endif

#ifndef SIMPLE_XCHG
#define SIMPLE_XCHG 0
#elif SIMPLE_XCHG == 1 && !defined(_CIVL)
#error "SIMPLE_XCHG can only be used with CIVL"
#endif

#endif //ARMON_CIVL_VARS_H
