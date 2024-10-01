#ifndef ARMON_CIVL_BLOCK_INTERFACE_H
#define ARMON_CIVL_BLOCK_INTERFACE_H

#include "utils.h"

enum BlockExchangeState {
    XCHG_NotReady   = 0,
    XCHG_InProgress = 1,
    XCHG_Done       = 3,
};

struct BlockInterface {
#if SIMPLE_XCHG
    bool ready[2];
    enum BlockExchangeState bint_state;
#else
    AtomicVar byte int_state;
#endif
    bool is_done[2];
};

void block_interface_state(struct BlockInterface* bint, enum BlockExchangeState* bint_state, byte* bint_ready);
void interface_side_ready(struct BlockInterface* bint, byte ready_flag, enum BlockExchangeState* bint_state, byte* bint_ready);
bool interface_start_exchange(struct BlockInterface* bint, byte side_flag);
void interface_end_exchange(struct BlockInterface* bint);
bool interface_acknowledge_exchange(struct BlockInterface* bint, byte side_flag);
bool mark_ready_for_exchange(struct BlockInterface* bint, bool is_first_side, enum BlockExchangeState* new_state);
enum BlockExchangeState exchange_done(struct BlockInterface* bint);

#endif //ARMON_CIVL_BLOCK_INTERFACE_H
