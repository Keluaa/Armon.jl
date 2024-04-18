
@generic_kernel function boundary_conditions!(
    ρ::V, u::V, v::V, p::V, c::V, g::V, E::V,
    bsize::BlockSize, axis::Axis.T, side::Side.T,
    u_factor::T, v_factor::T
) where {T, V <: AbstractArray{T}}
    @kernel_init begin
        # `incr` is the stride of `axis` going towards the edge along `side`
        incr = stride_along(bsize, axis)
        incr = ifelse(side in first_sides(ndims(bsize)), -incr, incr)
    end

    i = @index_2D_lin()
    ig = i + incr  # index of the ghost cell

    # TODO: use `comm_vars()` and iterate it like `pack_to_array` ?
    for _ in 1:ghosts(bsize)
        ρ[ig] = ρ[i]
        u[ig] = u[i] * u_factor
        v[ig] = v[i] * v_factor
        p[ig] = p[i]
        c[ig] = c[i]
        g[ig] = g[i]
        E[ig] = E[i]

        i  -= incr
        ig += incr
    end
end


function boundary_conditions!(params::ArmonParameters{T}, state::SolverState, blk::LocalTaskBlock, side::Side.T) where {T}
    (u_factor::T, v_factor::T) = boundary_condition(state.test_case, side)
    domain = border_domain(blk.size, side)
    boundary_conditions!(params, block_device_data(blk), domain, blk.size, state.axis, side, u_factor, v_factor)
end


@kernel_function function vars_ghost_exchange(
    vars₁::NTuple{N, V}, i₁, ig₁, sg₁,
    vars₂::NTuple{N, V}, i₂, ig₂, sg₂,
    ghosts
) where {N, V}
    for gᵢ in 0:ghosts-1
        j₁ = gᵢ * sg₁
        j₂ = gᵢ * sg₂

        # Real cells of 1 to ghosts of 2
        # TODO: KernelAbstractions.Extras.@unroll ??
        for v in 1:N
            vars₂[v][ig₂ - j₂] = vars₁[v][i₁ + j₁]
        end

        # Real cells of 2 to ghosts of 1
        # TODO: KernelAbstractions.Extras.@unroll ??
        for v in 1:N
            vars₁[v][ig₁ - j₁] = vars₂[v][i₂ + j₂]
        end
    end
end


@generic_kernel function block_ghost_exchange(
    vars₁::NTuple{N, V},
    vars₂::NTuple{N, V},
    bsize::BlockSize, axis::Axis.T, side₁::Side.T
) where {N, V}
    @kernel_init begin
        side₂ = opposite_of(side₁)

        # Offsets going towards the ghost cells
        sg  = stride_along(bsize, axis)
        sg₁ = ifelse(side₁ in first_sides(ndims(bsize)), -sg, sg)
        sg₂ = -sg₁

        # Convertion from `i₁` to `i₂`, exploiting the fact that both blocks have the same size
        d₂ = stride_along(bsize, axis) * (real_size_along(bsize, axis) - 1)
        d₂ = ifelse(side₂ in first_sides(ndims(bsize)), -d₂, d₂)
    end

    i₁ = @index_2D_lin()
    i₂ = i₁ + d₂

    # `ig` is the position of the ghost cell at the border of the block
    ig₁ = i₁ + sg₁ * ghosts(bsize)
    ig₂ = i₂ + sg₂ * ghosts(bsize)

    # `i` is the position of the farthest real cell from the ghost border of the block
    i₁ -= sg₁ * (ghosts(bsize) - 1)
    i₂ -= sg₂ * (ghosts(bsize) - 1)

    vars_ghost_exchange(
        vars₁, i₁, ig₁, sg₁,
        vars₂, i₂, ig₂, sg₂,
        ghosts(bsize)
    )
end


function block_ghost_exchange(
    params::ArmonParameters, state::SolverState,
    blk₁::LocalTaskBlock{V, Size}, blk₂::LocalTaskBlock{V, Size}, side::Side.T
) where {V, Size <: StaticBSize}
    do_xchg, xchg_state = mark_ready_for_exchange!(blk₁, side)
    !do_xchg && return xchg_state

    # Exchange between two blocks with the same dimensions
    domain = border_domain(blk₁.size, side)
    block_ghost_exchange(params, domain,
        comm_vars(blk₁), comm_vars(blk₂),
        blk₁.size, state.axis, side
    )

    return exchange_done!(blk₁, side)
end


@generic_kernel function block_ghost_exchange(
    vars₁::NTuple{N, V}, bsize₁::BlockSize,
    vars₂::NTuple{N, V}, bsize₂::BlockSize,
    axis::Axis.T, side₁::Side.T
) where {N, V}
    @kernel_init begin
        side₂ = opposite_of(side₁)

        # Offsets going towards the ghost cells
        sg₁ = stride_along(bsize₁, axis)
        sg₂ = stride_along(bsize₂, axis)
        sg₁ = ifelse(side₁ in first_sides(ndims(bsize₁)), -sg₁, sg₁)
        sg₂ = ifelse(side₂ in first_sides(ndims(bsize₁)), -sg₂, sg₂)
    end

    i₁ = @index_2D_lin()

    # `bsize₁` and `bsize₂` are different, therefore such is the iteration domain. We translate the
    # `i₁` index to its reciprocal `i₂` on the other side using the nD index.
    I₁ = position(bsize₁, i₁)

    # TODO: cleanup
    I₂x = ifelse(side₂ in sides_along(Axis.X), ifelse(side₂ in first_sides(ndims(bsize₁)), 1, real_block_size(bsize₂)[1]), I₁[1])
    I₂y = ifelse(side₂ in sides_along(Axis.Y), ifelse(side₂ in first_sides(ndims(bsize₁)), 1, real_block_size(bsize₂)[2]), I₁[2])
    I₂ = (I₂x, I₂y)
    # TODO: Unreadable but efficient and dimension-agnostic?
    # I₂ = ifelse.(
    #     in.(side₂, (sides_along(Axis.X), sides_along(Axis.Y))),
    #     ifelse.(side₂ in first_sides(ndims(bsize₁)), 1, real_block_size(bsize₂)),
    #     I₁
    # )

    i₂ = lin_position(bsize₂, I₂)

    # `ig` is the position of the ghost cell at the border of the block
    ig₁ = i₁ + sg₁ * ghosts(bsize₁)
    ig₂ = i₂ + sg₂ * ghosts(bsize₁)

    # `i` is the position of the farthest real cell from the ghost border of the block
    i₁ -= sg₁ * (ghosts(bsize₁) - 1)
    i₂ -= sg₂ * (ghosts(bsize₁) - 1)

    vars_ghost_exchange(
        vars₁, i₁, ig₁, sg₁,
        vars₂, i₂, ig₂, sg₂,
        ghosts(bsize₁)
    )
end


function block_ghost_exchange(
    params::ArmonParameters, state::SolverState,
    blk₁::LocalTaskBlock{V}, blk₂::LocalTaskBlock{V}, side::Side.T
) where {V}
    do_xchg, xchg_state = mark_ready_for_exchange!(blk₁, side)
    !do_xchg && return xchg_state

    # Exchange between two blocks with (possibly) different dimensions, but the same length along `side`
    domain = border_domain(blk₁.size, side)
    block_ghost_exchange(params, domain,
        comm_vars(blk₁), blk₁.size,
        comm_vars(blk₂), blk₂.size,
        state.axis, side
    )

    return exchange_done!(blk₁, side)
end


@generic_kernel function pack_to_array!(
    bsize::BlockSize, side::Side.T, array::V, vars::NTuple{N, V}
) where {N, V}
    idx = @index_2D_lin()
    itr = @iter_idx()

    (i, i_g) = divrem(itr - 1, ghosts(bsize))
    i_arr = (i_g * real_face_size(bsize, side) + i) * N

    # TODO: KernelAbstractions.Extras.@unroll ??
    for v in 1:N
        array[i_arr+v] = vars[v][idx]
    end
end


@generic_kernel function unpack_from_array!(
    bsize::BlockSize, side::Side.T, array::V, vars::NTuple{N, V}
) where {N, V}
    idx = @index_2D_lin()
    itr = @iter_idx()

    (i, i_g) = divrem(itr - 1, ghosts(bsize))
    i_arr = (i_g * real_face_size(bsize, side) + i) * N

    # TODO: KernelAbstractions.Extras.@unroll ??
    for v in 1:N
        vars[v][idx] = array[i_arr+v]
    end
end


"""
    start_exchange(
        params::ArmonParameters,
        blk::LocalTaskBlock{D, H}, other_blk::RemoteTaskBlock{B}, side::Side.T
    ) where {D, H, B}

Start the exchange between one local block and a remote block from another sub-domain.
Returns `true` if the exchange is [`BlockExchangeState.Done`](@ref), `false` if
[`BlockExchangeState.InProgress`](@ref).
"""
function start_exchange(
    params::ArmonParameters,
    blk::LocalTaskBlock{D, H}, other_blk::RemoteTaskBlock{B}, side::Side.T
) where {D, H, B}
    buffer_are_on_device = D == B
    if !buffer_are_on_device
        # MPI buffers are not located where the up-to-date data is: we must to a copy first.
        device_to_host!(blk)
    end

    send_domain = border_domain(blk.size, side; single_strip=false)
    vars = comm_vars(blk; on_device=buffer_are_on_device)
    # TODO: run on host if `D != B`, or perform it on the device on a tmp array
    pack_to_array!(params, send_domain, blk.size, side, other_blk.send_buf.data, vars)

    wait(params)  # Wait for the copy to complete

    # TODO: use RMA with processes local to the node.
    MPI.Startall(other_blk.requests)

    return false
end


"""
    finish_exchange(
        params::ArmonParameters,
        blk::LocalTaskBlock{D, H}, other_blk::RemoteTaskBlock{B}, side::Side.T
    ) where {D, H, B}

Finish the exchange between one local block and a remote block from another sub-domain, if MPI
communications are done.
Returns `true` if the exchange is [`BlockExchangeState.Done`](@ref), `false` if
[`BlockExchangeState.InProgress`](@ref).
"""
function finish_exchange(
    params::ArmonParameters,
    blk::LocalTaskBlock{D, H}, other_blk::RemoteTaskBlock{B}, side::Side.T
) where {D, H, B}
    # Finish the exchange between one local block and a remote block from another sub-domain
    !MPI.Testall(other_blk.requests) && return false  # Still waiting

    recv_domain = ghost_domain(blk.size, side; single_strip=false)
    buffer_are_on_device = D == B
    vars = comm_vars(blk; on_device=buffer_are_on_device)
    # TODO: run on host if `D != B`
    unpack_from_array!(params, recv_domain, blk.size, side, other_blk.recv_buf.data, vars)

    if !buffer_are_on_device
        # MPI buffers are not where we want the data to be. Retreive the result of the exchange.
        host_to_device!(blk)
    end

    return true
end


function block_ghost_exchange(
    params::ArmonParameters, state::SolverState,
    blk::LocalTaskBlock{D, H}, other_blk::RemoteTaskBlock{B}, side::Side.T
) where {D, H, B}
    if other_blk.rank == -1
        # `other_blk` is fake, this is the border of the global domain
        boundary_conditions!(params, state, blk, side)
        return BlockExchangeState.Done
    end

    bint = blk.exchanges[Int(side)]
    bint_state = block_interface_state(bint)[1]

    # Exchange between one local block and a remote block from another sub-domain
    if bint_state == BlockExchangeState.NotReady
        exchange_ended = start_exchange(params, blk, other_blk, side)
        side_flag = side in first_sides() ? 0b10 : 0b01
        !exchange_ended && interface_start_exchange!(bint, side_flag; for_MPI=true)
    else
        exchange_ended = finish_exchange(params, blk, other_blk, side)
        exchange_ended && interface_end_exchange!(bint; for_MPI=true)
    end

    return exchange_ended ? BlockExchangeState.Done : BlockExchangeState.InProgress
end


"""
    block_ghost_exchange(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)

Handles communications between the `blk` neighbours, along the current `state.axis`.
If `blk` is on one of the edges of the grid, a remote exchange is performed with the neighbouring
[`RemoteTaskBlock`](@ref), or the global boundary conditions are applied.

Returns `true` if exchanges were not completed, and the block is waiting on another to be ready for
the exchange.
"""
function block_ghost_exchange(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    # TODO: skip the whole interface logic if both blocks are in the same thread

    # Exchange with the Left/Bottom neighbour
    side = first_side(state.axis)
    if !is_side_done(blk, side)
        other_blk = blk.neighbours[side]
        left_exchange_state = block_ghost_exchange(params, state, blk, other_blk, side)
        left_exchange_state == BlockExchangeState.Done && side_is_done!(blk, side, true)
    else
        left_exchange_state = BlockExchangeState.Done
    end

    # Exchange with the Right/Top neighbour
    other_side = opposite_of(side)
    if !is_side_done(blk, other_side)
        other_blk = blk.neighbours[other_side]
        right_exchange_state = block_ghost_exchange(params, state, blk, other_blk, other_side)
        right_exchange_state == BlockExchangeState.Done && side_is_done!(blk, other_side, true)
    else
        right_exchange_state = BlockExchangeState.Done
    end

    if left_exchange_state == right_exchange_state == BlockExchangeState.Done
        # Both sides are `Done`: this exchange step is finished, reset the interface.
        side_is_done!(blk, side, false)
        side_is_done!(blk, other_side, false)
        return false
    else
        return true  # Still waiting for neighbours
    end
end


function block_ghost_exchange(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    # We must repeatedly update all blocks' states until the exchanges are done, as they are designed
    # to work independantly and asynchronously in a state machine, which isn't the case here.
    waiting_for = ones(Bool, grid.grid_size)  # Cannot use `BitArray` because of multithreading
    while any(waiting_for)
        @iter_blocks for blk in grid
            if waiting_for[blk.pos] && !block_ghost_exchange(params, state, blk)
                waiting_for[blk.pos] = false
            end
        end
    end
end
