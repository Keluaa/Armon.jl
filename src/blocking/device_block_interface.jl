
# Using UInt8 would be better for memory usage, but Atomic operations on small operands are not present on GPU
# Using `@enum` or `@enumx` isn't possible since we want to increment values with `+=`
module DeviceBlockInterfaceState
    const T = UInt32

    "No block is ready"
    const NotReady      :: T = 0
    "One block is ready"
    const OneReady      :: T = 1
    "Both blocks are ready, the first one to set this state will do the exchange"
    const BothReady     :: T = 2
    "One of the blocks is doing the exchange"
    const DoingExchange :: T = 3
    "The exchange is done, waiting for the other block to acknowledge this"
    const ExchangeDone  :: T = 4
end


module DeviceBlockInterfaceStatus
    const T = UInt32

    "The exchange isn't done"
    const NotReady        :: T = 0
    "The block is ready for the exchange, the interface at least [`DeviceBlockInterfaceState.OneReady`](@ref), we wait for the other block to do the exchange"
    const WaitForExchange :: T = 1
    "The block is ready and is doing the exchange, the interface is at [`DeviceBlockInterfaceState.DoingExchange`](@ref) too"
    const DoingExchange   :: T = 2
    "The block has completed the exchange, the interface is either [`DeviceBlockInterfaceState.ExchangeDone`](@ref) or [`DeviceBlockInterfaceState.NotReady`](@ref) if already reset"
    const Done            :: T = 3
end


struct GridInterfaces{
    InterfaceStates   <: AbstractArray{DeviceBlockInterfaceState.T},
    InterfaceStatuses <: AbstractArray{DeviceBlockInterfaceStatus.T}
}
    states   :: InterfaceStates    # state of interfaces, shared by two blocks
    statuses :: InterfaceStatuses  # status of a block regarding one of its interfaces (not shared/not atomic)
end

Adapt.@adapt_structure GridInterfaces


function GridInterfaces(grid_size::Dims{D}, device) where {D}
    block_sides_count = prod(grid_size) * 2D  # 2 sides per axis. For simplicity, includes also sides with no neighbours
    interface_count = sum(1:D) do d
        grid_size[d] == 1 && return 0
        int_grid = ntuple(D) do i; ifelse(i == d, grid_size[i] - 1, grid_size[i]) end
        return prod(int_grid)
    end

    array_type     = device_array_type(device)
    int_state      = array_type{DeviceBlockInterfaceState.T}(undef, interface_count)
    int_blk_status = array_type{DeviceBlockInterfaceStatus.T}(undef, block_sides_count)

    interfaces = GridInterfaces(int_state, int_blk_status)
    reset!(interfaces)

    return interfaces
end


function reset!(interfaces::GridInterfaces)
    fill!(interfaces.states, DeviceBlockInterfaceState.NotReady)
    fill!(interfaces.statuses, DeviceBlockInterfaceStatus.NotReady)
    return interfaces
end


function block_interface_index(grid_size::Dims{D}, block_pos::CartesianIndex{D}, side::Side.T) where {D}
    if side in first_sides()
        # The previous block along the side's axis stores the interface.
        block_pos = block_pos + CartesianIndex(offset_to(side))
        side = opposite_of(side)  # Since it is the opposite block of the interface, it is the opposite side
    end

    if block_pos ∉ CartesianIndices(grid_size) || (block_pos + CartesianIndex(offset_to(side))) ∉ CartesianIndices(grid_size)
        return -1  # An interface to (or from) a remote block (or a global boundary)
    end

    # Interfaces are stored contiguously by axis: all along X then all along Y.
    dim_idx = Int(axis_of(side))
    int_idx = 0
    for d in 1:D
        # `int_grid` is the grid of interfaces along axis `d`.
        # Each block stores the Right (or Top) interface, therefore the last blocks of the grid
        # along `d` will have no interface assigned, hence the interface grid is one block shorter
        # in that direction.
        int_grid = ntuple(D) do i; ifelse(i == d, grid_size[i] - 1, grid_size[i]) end
        if d == dim_idx
            # We know `block_pos` must be in `int_grid` because of the checks above
            int_idx += @inbounds LinearIndices(int_grid)[block_pos]
            break
        else
            int_idx += prod(int_grid)
        end
    end

    return int_idx
end


function base_block_interface_status_index(grid_size::Dims{D}, block_pos::CartesianIndex{D}) where {D}
    # Each status is unique to each block: there is `prod(grid_size) * 2 * D` statuses in total
    blk_idx = (LinearIndices(grid_size)[block_pos] - 1) * 2D  # 0-index
    # Then `blk_idx + Int(side)` would give the index of a side's status in `GridInterfaces.statuses`
    return blk_idx
end


function interface_exchange!(interfaces::GridInterfaces, interface_idx, block_status_idx)
    # Returns `is_xchg_completed, should_do_the_exchange`

    blk_status = interfaces.statuses[block_status_idx]
    if interface_idx < 0
        # Special case: negative indices indicate that there is no neighbour, but either a remote block
        # or a global boundary.
        # Those two can always be done since they do not depend on another local block.
        interfaces.statuses[block_status_idx] = DeviceBlockInterfaceStatus.Done  # "do it only once"
        return true, blk_status != DeviceBlockInterfaceStatus.Done

    elseif blk_status == DeviceBlockInterfaceStatus.NotReady
        # TODO: simple atomic loads are not supported, so we use an add
        int_state = (Atomix.@atomic interfaces.states[interface_idx] += 0)
        if int_state > DeviceBlockInterfaceState.BothReady
            # The previous exchange isn't completed: we must wait for the other block to acknowledge
            # it and reset the interface.
            return false, false
        end

        # Mark this block as ready
        new = (Atomix.@atomic interfaces.states[interface_idx] += DeviceBlockInterfaceState.T(1))
        # TODO: replace this `if` in case there is problems
        if new == DeviceBlockInterfaceState.OneReady
            # It is certain that the other block isn't ready, no need to attempt the CAS
            interfaces.statuses[block_status_idx] = DeviceBlockInterfaceStatus.WaitForExchange
            return false, false
        end

        # TODO: it may be safe to remove this CAS since the atomic add above might be acting as one on all GPUs (check)
        res = Atomix.@atomicreplace(
            interfaces.states[interface_idx],
            DeviceBlockInterfaceState.BothReady => DeviceBlockInterfaceState.DoingExchange
        )

        if res.success
            # This block will do the exchange, the other will wait
            interfaces.statuses[block_status_idx] = DeviceBlockInterfaceStatus.DoingExchange
            return true, true
        else
            interfaces.statuses[block_status_idx] = DeviceBlockInterfaceStatus.WaitForExchange
            return false, false
        end

    elseif blk_status == DeviceBlockInterfaceStatus.WaitForExchange
        # Is the exchange done? Reset the interface if so.
        # TODO: using a CAS for this is overkill, but might still be more performant than a load followed by a store => check
        res = Atomix.@atomicreplace(
            interfaces.states[interface_idx],
            DeviceBlockInterfaceState.ExchangeDone => DeviceBlockInterfaceState.NotReady
        )
        if res.success
            interfaces.statuses[block_status_idx] = DeviceBlockInterfaceStatus.Done
            return true, false
        else
            return false, false
        end

    elseif blk_status == DeviceBlockInterfaceStatus.DoingExchange
        # This is reached only when the exchange is completed by the current block.
        # TODO: no need for a acknowledge for this style of exchange protocol right?
        # TODO: simple atomic stores are not supported, so we use an atomic swap instead
        Atomix.@atomicswap interfaces.states[interface_idx] = DeviceBlockInterfaceState.ExchangeDone
        interfaces.statuses[block_status_idx] = DeviceBlockInterfaceStatus.Done
        return true, false

    else
        # DeviceBlockInterfaceStatus.Done
        return true, false
    end
end


function reset_exchange!(interfaces::GridInterfaces, block_status_idx)
    # Since the interface state is reset automatically, we only have to reset the block state
    # regarding that interface, which isn't shared and therefore no atomic operations are needed.
    interfaces.statuses[block_status_idx] = DeviceBlockInterfaceStatus.NotReady
    return
end
