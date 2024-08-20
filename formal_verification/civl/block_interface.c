
#include "block_interface.h"

CIVL_pure static enum BlockExchangeState exchange_state_from_val(byte val)
{
    CIVL_assume(0 <= val && val <= 3);
    switch (val) {
    case 0:  return XCHG_NotReady;
    case 1:  return XCHG_InProgress;
    case 3:  return XCHG_Done;
    default: return XCHG_NotReady;
    }
}

CIVL_pure static byte first_two_bits(byte val)
{
    CIVL_assume(0 <= val && val <= 15);
#ifdef _CIVL
    return val - ((val / 4) * 4);
#else
    return val & 3;
#endif
}

CIVL_pure static byte last_two_bits(byte val)
{
    // last two bits of a 4 bit value
    CIVL_assume(0 <= val && val <= 15);
#ifdef _CIVL
    return val / 4;
#else
    return (val & 3) >> 2;
#endif
}

CIVL_pure static byte two_bit_and(byte a, byte b)
{
    CIVL_assume(0 <= a <= 3 && 0 <= b <= 3);
    switch (a) {
    default:
    case 0: return 0;
    case 1: return (b == 1 || b == 3) ? 1 : 0;
    case 2: return (b == 2 || b == 3) ? 2 : 0;
    case 3: return b;
    }
}

CIVL_pure static byte two_bit_or(byte a, byte b)
{
    CIVL_assume(0 <= a <= 3 && 0 <= b <= 3);
    switch (a) {
    default:
    case 0: return b;
    case 1: return (b == 2 || b == 3) ? 3 : 1;
    case 2: return (b == 1 || b == 3) ? 3 : 2;
    case 3: return 3;
    }
}

#ifdef _CIVL
#define CIVL_atomic_or_fetch(val, arg, res) $atomic { (val) = two_bit_or((val), (arg)); (res) = (val); }
#else
#define CIVL_atomic_or_fetch(val, arg, res) \
    (res) = (atomic_fetch_or(&(val), (arg)) | (arg))
#endif

void block_interface_state(struct BlockInterface* bint, enum BlockExchangeState* bint_state, byte* bint_ready)
{
    byte bint_flags;
    CIVL_atomic_load(bint->int_state, bint_flags);
    *bint_state = exchange_state_from_val(last_two_bits(bint_flags));
    *bint_ready = first_two_bits(bint_flags);
}

void interface_side_ready(struct BlockInterface* bint, byte ready_flag, enum BlockExchangeState* bint_state, byte* bint_ready)
{
    byte flags = first_two_bits(ready_flag);
    byte bint_flags;
    CIVL_atomic_or_fetch(bint->int_state, flags, bint_flags);
    *bint_state = exchange_state_from_val(last_two_bits(bint_flags));
    *bint_ready = first_two_bits(bint_flags);
}

bool interface_start_exchange(struct BlockInterface* bint, byte side_flag)
{
    byte other_side_flag = side_flag == 1 ? 2 : 1;  // opposite side flag (eq to 'side_flag ^ 0b11')
    byte ready_state  = ((byte) XCHG_NotReady)   * 4 + 3;
    byte target_state = ((byte) XCHG_InProgress) * 4 + other_side_flag;
    bool success;
    CIVL_atomic_cas(success, bint->int_state, ready_state, target_state);
    return success;
}

void interface_end_exchange(struct BlockInterface* bint)
{
    byte done_flag = ((byte) XCHG_Done) * 4;
    byte _;
    CIVL_atomic_or_fetch(bint->int_state, done_flag, _);
}

bool interface_acknowledge_exchange(struct BlockInterface* bint, byte side_flag)
{
    byte current_state = (((byte) XCHG_Done)     * 4) + side_flag;
    byte target_state  = (((byte) XCHG_NotReady) * 4) + 0;
    bool success;
    CIVL_atomic_cas(success, bint->int_state, current_state, target_state);
    return success;
}

bool mark_ready_for_exchange(struct BlockInterface* bint, bool is_first_side, enum BlockExchangeState* new_state)
{
    enum BlockExchangeState bint_state;
    byte ready_flags;
    block_interface_state(bint, &bint_state, &ready_flags);
    byte side_flag = is_first_side ? 2 : 1;

    switch (bint_state) {
    case XCHG_InProgress:
        // The other block is still doing the exchange
        *new_state = bint_state;
        return false;
    case XCHG_Done:
        // One of the blocks did the exchange
        if (interface_acknowledge_exchange(bint, side_flag)) {
            // It was the other one, now the interface is reset and we can continue
            *new_state = XCHG_Done;
        } else {
            // It was this one, we are waiting for the other block to acknowledge it
            *new_state = XCHG_NotReady;
        }
        return false;
    default:
        break;
    }

    // If `bint_state` is `NotReady`, we can safely set our flag.
    CIVL_assert(bint_state == XCHG_NotReady);
    if (two_bit_and(ready_flags, side_flag) == 0) {
        interface_side_ready(bint, two_bit_or(ready_flags, side_flag), &bint_state, &ready_flags);
    }

    if (ready_flags == 3) {
        // Both sides are ready
        if (interface_start_exchange(bint, side_flag)) {
            // This block will do the exchange
            *new_state = XCHG_InProgress;
            return true;
        } else {
            // The other block will do the exchange
            *new_state = XCHG_InProgress;
            return false;
        }
    } else {
        // Wait for the other side to be ready
        *new_state = XCHG_NotReady;
        return false;
    }
}

enum BlockExchangeState exchange_done(struct BlockInterface* bint)
{
    interface_end_exchange(bint);
    return XCHG_Done;
}
