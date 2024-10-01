#ifndef ARMON_CIVL_DEBUG_H
#define ARMON_CIVL_DEBUG_H

#include "blocks.h"

const char* side_to_str(enum Side side);
const char* block_state_to_str(enum BlockState state);
const char* time_step_state_to_str(enum TimeStepState time_step_state);
const char* block_xchg_state_to_str(enum BlockExchangeState block_exchange_state);

void print_interface(struct BlockInterface* interface);
void print_block(struct Block* block);
void print_block_grid(struct BlockGrid* block_grid);
void print_global_state(void);
void pretty_print_solver_state(struct BlockGrid* block_grid);

bool check_interfaces(struct BlockGrid* block_grid);
void check_end_state(struct BlockGrid* block_grid);

#endif //ARMON_CIVL_DEBUG_H
