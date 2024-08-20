#ifndef ARMON_CIVL_TIME_STEP_H
#define ARMON_CIVL_TIME_STEP_H

#include "blocks.h"

void update_dt();
enum TimeStepState wait_for_dt();
void contribute_to_dt(struct Block* block);
void next_cycle();
void next_time_step(struct Block* block);

#endif //ARMON_CIVL_TIME_STEP_H
