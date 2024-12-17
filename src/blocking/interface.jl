
"""
    BlockExchangeState

State of an interface between two blocks, controlling if a cell exchange can happen or did happen.
"""
@enumx BlockExchangeState::UInt8 begin
    "One of the blocks is not ready yet"
    NotReady=0b00
    "One of the blocks is performing the exchange"
    InProgress=0b01
    "The exchange is done"
    Done=0b11
    _err=0b10  # dummy error state
end


function exchange_state_from_val(v::UInt8)
    # Exception-less alternative to `BlockExchangeState.T(v)`
    if     v == 0b00 return BlockExchangeState.NotReady
    elseif v == 0b01 return BlockExchangeState.InProgress
    elseif v == 0b11 return BlockExchangeState.Done
    else             return BlockExchangeState._err
    end
end


function block_interface_state(bint::BlockInterface)
    bint_flags = @atomic bint.flags.x
    bint_state = exchange_state_from_val((bint_flags & 0b1100) >> 2)
    bint_ready = bint_flags & 0b0011
    return bint_state, bint_ready
end


function interface_side_ready!(bint::BlockInterface, ready_flag::UInt8)
    flags = (ready_flag & 0b0011)
    bint_flags = @atomic bint.flags.x |= flags
    bint_state = exchange_state_from_val((bint_flags & 0b1100) >> 2)
    bint_ready = bint_flags & 0b0011
    return bint_state, bint_ready
end


function interface_start_exchange!(bint::BlockInterface, side_flag::UInt8; for_MPI=false)
    if !for_MPI
        # Ready && 0b11  =>  InProgress & other_side_flag
        # Once the exchange is done, only the other side will be able to reset the interface.
        other_side_flag = side_flag ⊻ 0b11
        ready_state = (UInt8(BlockExchangeState.NotReady) << 2) | 0b0011
        target_state = (UInt8(BlockExchangeState.InProgress) << 2) | other_side_flag
        _, ok = @atomicreplace bint.flags.x ready_state => target_state
        return ok
    else
        @atomic bint.flags.x = (UInt8(BlockExchangeState.InProgress) << 2)
        return true
    end
end


function interface_end_exchange!(bint::BlockInterface; for_MPI=false)
    if !for_MPI
        done_flag = UInt8(BlockExchangeState.Done) << 2  # Requires that `BlockExchangeState.Done == 0b11`
        @atomic bint.flags.x |= done_flag
    else
        @atomic bint.flags.x = UInt8(BlockExchangeState.NotReady) << 2
    end
end


function interface_acknowledge_exchange!(bint::BlockInterface, side_flag::UInt8)
    # Will only work if `side_flag` is of the side which did not do the exchange
    current_state = (UInt8(BlockExchangeState.Done) << 2) | side_flag
    target_state = (UInt8(BlockExchangeState.NotReady) << 2) | 0b00
    _, ok = @atomicreplace bint.flags.x current_state => target_state
    return ok
end


function reset!(bint::BlockInterface)
    @atomic bint.flags.x = 0
    bint.is_left_side_done = false
    bint.is_right_side_done = false
end


"""
    mark_ready_for_exchange!(blk, side)

Mark `blk` in the interface along its `side` as ready for an exchange.

Return `true` if the `blk` should do the exchange, and the new [`BlockExchangeState`](@ref) of the interface.
"""
function mark_ready_for_exchange!(blk::LocalTaskBlock, side::Side.T)
    bint = blk.exchanges[side]
    bint_state, ready_flags = block_interface_state(bint)
    side_flag = first_side(side) ? 0b10 : 0b01

    if bint_state == BlockExchangeState.InProgress
        # The other block is still doing the exchange
        return false, bint_state
    elseif bint_state == BlockExchangeState.Done
        # One of the blocks did the exchange
        if interface_acknowledge_exchange!(bint, side_flag)
            # It was the other one, now the interface is reset and we can continue
            return false, BlockExchangeState.Done
        else
            # It was this one, we are waiting for the other block to acknowledge it
            return false, BlockExchangeState.NotReady
        end
    elseif bint_state != BlockExchangeState.NotReady
        error("concurrency error at $(Tuple(blk.pos)) along $side, wrong interface state: $bint_state")
    end

    # If `bint_state` is `NotReady`, we can safely set our flag.
    if (ready_flags & side_flag) == 0
        bint_state, ready_flags = interface_side_ready!(bint, ready_flags | side_flag)
    end

    if ready_flags == 0b11
        # Both sides are ready
        if interface_start_exchange!(bint, side_flag)
            # This block will do the exchange
            return true, BlockExchangeState.InProgress
        else
            # The other block will do the exchange
            return false, BlockExchangeState.InProgress
        end
    else
        # Wait for the other side to be ready
        return false, BlockExchangeState.NotReady
    end
end


"""
    exchange_done!(blk, side)

Mark the exchange of the interface of the `blk` along `side` as done.
The other block need to acknowledge this before the exchange state can be reset.
"""
function exchange_done!(blk::LocalTaskBlock, side::Side.T)
    # Mark the exchange as done for the other block to acknowledge it
    bint = blk.exchanges[side]
    interface_end_exchange!(bint)
    return BlockExchangeState.Done
end


"""
    is_side_done(blk::LocalTaskBlock, side::Side.T)

The value set by [`side_is_done!`](@ref).
"""
function is_side_done(blk::LocalTaskBlock, side::Side.T)
    bint = blk.exchanges[side]
    # left/right are reversed as sides are relative to the interface
    return first_side(side) ? bint.is_right_side_done : bint.is_left_side_done
end


"""
    side_is_done!(blk::LocalTaskBlock, side::Side.T, done::Bool=true)

Set the interface of `blk` along `side` as `done`. This is a non-atomic operation, meant only to
avoid repeating the exchange logic multiple times.
"""
function side_is_done!(blk::LocalTaskBlock, side::Side.T, done::Bool=true)
    bint = blk.exchanges[side]
    if first_side(side)
        bint.is_right_side_done = done
    else
        bint.is_left_side_done = done
    end
end
