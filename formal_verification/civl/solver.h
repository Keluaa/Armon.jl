#ifndef ARMON_CIVL_SOLVER_H
#define ARMON_CIVL_SOLVER_H

#include "blocks.h"

bool solver_cycle_async(struct ThreadWorkload* thread_workload);
bool solver_thread(int tid, struct ThreadWorkload* thread_workload, struct BlockGrid* block_grid);

#endif //ARMON_CIVL_SOLVER_H
