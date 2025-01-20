
@inline function block_border_info(bsize::BlockSize, axis::Axis.T, side::Side.T)
    # TODO: use Int32? since blocks are small this could reduce register usage

    # Position of the first real cell of the border domain
    # Left:   x=G       y=0
    # Right:  x=Sx-G+1  y=0
    # Bottom: x=0       y=G
    # Top:    x=0       y=Sy-G+1
    pos = (
        ifelse(side in sides_along(Axis.X), ifelse(side in first_sides(), ghosts(bsize), real_block_size(bsize)[1] - ghosts(bsize) + 1), 0),
        ifelse(side in sides_along(Axis.Y), ifelse(side in first_sides(), ghosts(bsize), real_block_size(bsize)[2] - ghosts(bsize) + 1), 0),
    )
    base_index = lin_position(bsize, pos)

    # Stride between ghost cells on the `side`
    ghost_stride = stride_along(bsize, axis)
    ghost_stride = ifelse(side in first_sides(), -ghost_stride, ghost_stride)

    # Stride between real cells along the border side
    other_axis = axis == Axis.X ? Axis.Y : Axis.X
    stride = stride_along(bsize, other_axis)

    return base_index, ghost_stride, stride
end


@inline function border_exchange_indices(bsize::BlockSize, (base_index, ghost_stride, stride), i_side, i_ghost)
    # `i` is the position of the farthest real cell from the ghost border of the block
    i = base_index + i_side * stride

    # `ig` is the position of the farthest ghost cell from the ghost the border of the block
    ig = i + ghost_stride * (2 * ghosts(bsize) - 1)

    # Distance to the real/ghost cells border of the blocks
    j = (i_ghost - 1) * ghost_stride

    # "ghost cell index", "real cell index"
    return ig - j, i + j
end


function device_boundary_condition(
    KernelAbstractions.@context(),
    vars::NTuple{N, V}, factors::NTuple{N, T},  # TODO: using arrays instead of tuples could reduce register pressure
    bsize::BlockSize, axis::Axis.T, side::Side.T
) where {N, T, V <: AbstractArray{T}}
    border_info = block_border_info(bsize, axis, side)

    other_axis = axis == Axis.X ? Axis.Y : Axis.X
    side_length = real_block_size(bsize)[Int(other_axis)]
    xchg_domain = (side_length, ghosts(bsize), N)

    groupsize = prod(KernelAbstractions.@groupsize())
    tid = @index(Local, Linear)
    for i in tid:groupsize:prod(xchg_domain)
        # TODO: unroll the loop? there should be between 1 to 3 iterations max (for small blocks)

        # Since `N` is the last element of `xchg_domain`, it is the index which increases the slowest,
        # giving some opportunity for coalescing of memory accesses.
        xchg_idx = CartesianIndices(xchg_domain)[i]
        i_side, i_ghost, i_var = Tuple(xchg_idx)

        ghost_idx, real_idx = border_exchange_indices(bsize, border_info, i_side, i_ghost)

        var = @inbounds vars[i_var]
        factor = @inbounds factors[i_var]
        var[ghost_idx] = var[real_idx] * factor
    end
end


function device_block_exchange(
    KernelAbstractions.@context(),
    # TODO: using arrays instead of tuples could reduce register pressure
    vars₁::NTuple{N, V}, bsize₁::BlockSize,
    vars₂::NTuple{N, V}, bsize₂::BlockSize,
    axis::Axis.T, side₁::Side.T
) where {N, V}
    border_info₁ = block_border_info(bsize₁, axis, side₁)
    border_info₂ = block_border_info(bsize₂, axis, opposite_of(side₁))

    other_axis = axis == Axis.X ? Axis.Y : Axis.X
    side_length = real_block_size(bsize₁)[Int(other_axis)]
    xchg_domain = (side_length, ghosts(bsize₁), N)

    # There are `side_length * ghosts` cells to exchange, each with `N` variables.
    # Since the amount of cells can be much lower than the workgroup size (`56*4=224` vs. 1024), we
    # treat each variable as work to be divided among all threads, ensuring that the workload is
    # spread as evenly as possible.
    # Surprisingly, this compiles quite well to the GPU.
    # However, since the iteration space is a mess, this is certainly suboptimal, but with minimal
    # impact: there is never a lot of work to do here.
    groupsize = prod(KernelAbstractions.@groupsize())
    tid = @index(Local, Linear)
    for i in tid:groupsize:prod(xchg_domain)
        # TODO: unroll the loop? there should be between 1 to 3 iterations max (for small blocks)

        # Since `N` is the last element of `xchg_domain`, it is the index which increases the slowest,
        # giving some opportunity for coalescing of memory accesses.
        xchg_idx = CartesianIndices(xchg_domain)[i]
        i_side, i_ghost, i_var = Tuple(xchg_idx)

        ghost_idx₁, real_idx₁ = border_exchange_indices(bsize₁, border_info₁, i_side, i_ghost)
        ghost_idx₂, real_idx₂ = border_exchange_indices(bsize₂, border_info₂, i_side, i_ghost)

        var₁ = @inbounds vars₁[i_var]
        var₂ = @inbounds vars₂[i_var]

        # TODO: we might need to make loads and stores to the other block as atomic operations, as
        # they might be cached in an another SM. This would also mean that if the other block is on
        # another SM and doesn't do the exchange, then its up-to-date data might not be flushed and
        # we might read old data instead. This is very bad, but fixable by flushing the data before
        # marking the block as ready. By adding more steps to the exchange sequence, we might be
        # able to make it so that only the block which doesn't do the exchange flushes its data, but
        # this might not be more performant as more steps means more dancing around and kernel
        # launches, etc...
        # TODO: inbounds
        var₂[ghost_idx₂] = var₁[real_idx₁]
        var₁[ghost_idx₁] = var₂[real_idx₂]
    end
end


function device_block_packing(
    KernelAbstractions.@context(),
    # TODO: using arrays instead of tuples could reduce register pressure
    packed_array::P, vars::NTuple{N, V},
    bsize::BlockSize, side::Side.T,
    ::Val{packing}
) where {N, T, P <: AbstractArray{T}, V <: AbstractArray{T}, packing}
    border_info = block_border_info(bsize, axis, side)

    other_axis = axis == Axis.X ? Axis.Y : Axis.X
    side_length = real_block_size(bsize)[Int(other_axis)]
    xchg_domain = (side_length, ghosts(bsize), N)

    groupsize = prod(KernelAbstractions.@groupsize())
    tid = @index(Local, Linear)
    for i in tid:groupsize:prod(xchg_domain)
        # TODO: unroll the loop? there should be between 1 to 3 iterations max (for small blocks)

        # Since `N` is the last element of `xchg_domain`, it is the index which increases the slowest,
        # giving some opportunity for coalescing of memory accesses.
        xchg_idx = CartesianIndices(xchg_domain)[i]
        i_side, i_ghost, i_var = Tuple(xchg_idx)

        _, real_idx = border_exchange_indices(bsize, border_info, i_side, i_ghost)

        var = @inbounds vars[i_var]
        if packing
            packed_array[i] = var[real_idx]  # marshalling
        else
            var[real_idx] = packed_array[i]  # unmarshalling
        end
    end
end


function device_border_exchange(
    KernelAbstractions.@context(),
    grid::DeviceBlockGrid, block::DeviceLocalBlock, state::BasicSolverState,
    xchg_state
)
    side_1 = first_side(state.axis)
    side_2 = last_side(state.axis)

    # Let the thread 1 and 2 do the atomic logic for both interfaces, then all threads will cooperate
    # to perform the exchange.
    tid = @index(Local, Linear)
    if tid == 1
        interface_idx = block.interfaces_idx[Int(state.axis)][1]
        block_status_idx = block.base_status_idx + Int(side_1)
        xchg_state[1], xchg_state[2] = interface_exchange!(grid.interfaces, interface_idx, block_status_idx)
    elseif tid == 2
        interface_idx = block.interfaces_idx[Int(state.axis)][2]
        block_status_idx = block.base_status_idx + Int(side_2)
        xchg_state[3], xchg_state[4] = interface_exchange!(grid.interfaces, interface_idx, block_status_idx)
    else
        interface_idx = 0
        block_status_idx = 0
    end

    KernelAbstractions.@synchronize()

    for (side, do_xchg_idx) in zip((side_1, side_2), (2, 4))
        # Only do the border operations if this block will `do_xchg`

        !xchg_state[do_xchg_idx] && continue

        neighbour_pos = block.pos + CartesianIndex(offset_to(side))
        if in_grid(neighbour_pos, grid.grid_size)
            # Exchange between two local blocks
            neighbour_data, neighbour_size = @apply_device_block(
                neighbour_blk = grid[neighbour_pos],
            begin
                # Uniformize all block sizes to `DynamicBSize` to reduce the amount of combinasions
                # to compile for.
                neighbour_blk.data, DynamicBSize(block_size(neighbour_blk.size), ghosts(neighbour_blk.size))
            end)

            device_block_exchange(
                KernelAbstractions.@context(),
                comm_vars(block.data), block.size,
                comm_vars(neighbour_data), neighbour_size,
                state.axis, side
            )

            if tid == 1
                # Mark the exchange as done
                # To limit register pressure we recompute the interface indices
                interface_idx = block.interfaces_idx[Int(state.axis)][side in first_sides() ? 1 : 2]
                block_status_idx = block.base_status_idx + Int(side)
                interface_exchange!(grid.interfaces, interface_idx, block_status_idx)
            end
        else
            # Exchange with a remote block, or a global boundary
            remote_pos = block.pos + CartesianIndex(offset_to(side))
            remote_idx = grid.index_map[remote_pos + oneunit(CartesianIndex{2})]
            remote_blk = @inbounds grid.remote_blocks[remote_idx]

            if remote_blk.exists
                error("NYI")

                # Pack the border cells' data to the remote block's buffer for communication
                if remote_blk.on_device
                    device_block_packing(
                        KernelAbstractions.@context(),
                        remote_blk.buffer, comm_vars(block.data),
                        block.size, side, Val(true)
                    )
                    # TODO: how to trigger the communication?
                else
                    # TODO: the buffer is on the host, where to pack the data?
                    #  => if the buffer is mapped and the device allows it, we can write directly to it from the device
                    #  => otherwise we need a temporary array on the device
                    #    => one per stream: then the next blocks in the stream will not be able to pack their data
                    #    => one per block:  then we use more memory (and bandwidth?)
                    error("NYI")
                end

                # TODO: when do we unpack?
                #  => after packing, the block is stopped.
                #  => it is scheduled only after we know the communication is done
                #  => if we make it so that the GPU can notify the host that it can launch the
                #     communication, then there is some work to do about atomic ops and system scopes
                #     with CUDA (but is it available on other GPUs?)
                device_block_packing(
                    KernelAbstractions.@context(),
                    packed_array, comm_vars(block.data),
                    block.size, side, Val(false)
                )

                # TODO: when do we mark the exchange as completed?
                #  => what about another solver step, which does all of that?
            else
                # Global domain boundary condition
                u_factor, v_factor = boundary_condition(state.schemes.test_case, side)
                factors = (;
                    # TODO: this is not ideal, if `boundary_condition` could do the transformation
                    #   for us it would be better
                    ρ = eltype(block)(1),
                    u = eltype(block)(u_factor),
                    v = eltype(block)(v_factor),
                    E = eltype(block)(1),
                    p = eltype(block)(1),
                    c = eltype(block)(1),
                    g = eltype(block)(1),
                )
                device_boundary_condition(
                    KernelAbstractions.@context(),
                    comm_vars(block.data), Tuple(factors),
                    block.size, state.axis, side
                )
            end
        end
    end

    xchg_done = xchg_state[1] && xchg_state[3]
    if tid == 1 && xchg_done
        # Both sides are done, we can reset the exchange state of this block
        reset_exchange!(grid.interfaces, block.base_status_idx + Int(side_1))
        reset_exchange!(grid.interfaces, block.base_status_idx + Int(side_2))
    end

    return xchg_done
end

