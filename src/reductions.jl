
@inline @fast function dtCFL_kernel_reduction(u::T, v::T, c::T, mask::T, dx::T, dy::T) where T
    # We need the absolute value of the divisor since the result of the max can be negative,
    # because of some IEEE 754 non-compliance since fast math is enabled when compiling this code
    # for GPU, e.g.: `@fastmath max(-0., 0.) == -0.`, while `max(-0., 0.) == 0.`
    # If the mask is 0, then: `dx / -0.0 == -Inf`, which will then make the result incorrect.
    return min(
        dx / abs(max(abs(u + c), abs(u - c)) * mask),
        dy / abs(max(abs(v + c), abs(v - c)) * mask)
    )
end


@inline @fast function dtCFL_kernel_reduction(u::T, v::T, c::T, dx::T, dy::T) where T
    # Mask-less version
    return min(
        dx / abs(max(abs(u + c), abs(u - c))),
        dy / abs(max(abs(v + c), abs(v - c)))
    )
end


@fast function dtCFL_kernel(params::ArmonParameters{T, CPU_HP}, state::SolverState, blk::LocalTaskBlock, Δx::NTuple{2, T}) where {T}
    # CPU reduction
    (; u, v, c) = block_device_data(blk)
    range = block_domain_range(blk.size, state.steps_ranges[Int(state.axis)].real_domain)

    if params.use_cache_blocking
        # Reduction exploiting multithreading from the caller
        res = typemax(T)
        for j in range.col, i in range.row .+ (j - 1)
            cell_dt = dtCFL_kernel_reduction(u[i], v[i], c[i], Δx...)
            res = min(res, cell_dt)
        end
        return res
    else
        # Reduction using explicit multithreading, since the caller isn't multithreaded
        threads_res = Vector{T}(undef, params.nthreads)
        threads_res .= typemax(T)

        @threaded for j in range.col
            tid = Threads.threadid()
            res = threads_res[tid]
            for i in range.row .+ (j - 1)
                cell_dt = dtCFL_kernel_reduction(u[i], v[i], c[i], Δx...)
                res = min(res, cell_dt)
            end
            threads_res[tid] = res
        end

        return minimum(threads_res)
    end
end


@generic_kernel function dtCFL_kernel(
    u::V, v::V, c::V, res::V, bsize::BlockSize, dx::T, dy::T
) where {T, V <: AbstractArray{T}}
    i = @index_2D_lin()
    mask = T(!is_ghost(bsize, i))
    res[i] = dtCFL_kernel_reduction(u[i], v[i], c[i], mask, dx, dy)
end


function dtCFL_kernel(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, Δx::NTuple{2})
    # GPU generic reduction
    range = block_domain_range(blk.size, state.steps_ranges[Int(state.axis)].full_domain)
    blk_data = block_device_data(blk)

    if params.use_two_step_reduction
        # Use a temporary array to store the partial reduction result. This is inefficient but can
        # be more performant for some GPU backends.
        # We avoid filling the whole array with `typemax(T)` by applying the kernel on the whole
        # array (`full_domain`) and by using `is_ghost(i)` as a mask.
        # TODO: There may be some synchronization issues with oneAPI.jl.
        dtCFL_kernel(params, blk_data, range, blk_data.work_1, blk.size, Δx...)
        wait(params)
        return reduce(min, blk_data.work_1)
    else
        # Direct reduction, which depends on a pre-computed `mask`
        lin_range = first(range):last(range)  # Reduce on a 1D range
        c_v    = @view blk_data.c[lin_range]
        u_v    = @view blk_data.u[lin_range]
        v_v    = @view blk_data.v[lin_range]
        mask_v = @view blk_data.mask[lin_range]

        # Since GPUArrays.jl doesn't define a clean way of reducing with arrays and scalars, we must
        # put `Δx` in a closure.
        dtCFL_kernel_reduction_func(u, v, c, mask) = dtCFL_kernel_reduction(u, v, c, mask, Δx...)

        return mapreduce(dtCFL_kernel_reduction_func, min, u_v, v_v, c_v, mask_v)
    end
end


function local_time_step(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    Δx = params.domain_size ./ params.global_grid
    return dtCFL_kernel(params, state, blk, Δx)
end


function local_time_step(params::ArmonParameters{T}, state::SolverState, grid::BlockGrid) where {T}
    mt_reduction = params.use_threading && params.use_cache_blocking
    threads_res = Vector{T}(undef, mt_reduction ? params.nthreads : 1)
    threads_res .= typemax(T)

    @iter_blocks for blk in grid
        blk_res = local_time_step(params, state, blk)

        tid = mt_reduction ? Threads.threadid() : 1
        threads_res[tid] = min(blk_res, threads_res[tid])
    end

    return minimum(threads_res)
end


"""
    next_time_step(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock; already_contributed=false)
    next_time_step(params::ArmonParameters, state::SolverState, grid::BlockGrid)

Compute the time step of the next cycle. This is done at the start of the current cycle.

Since the current cycle does not rely on an up-to-date time step, the time step reduction is done
fully asynchronously, including the global MPI reduction.
The accuracy cost of this optimisation is minimal, as the CFL condition prevents the time step from
being too large.
Additionally, we prevent the time step from increasing of more than +5% of the previous one.

For first cycle, if no initial time step is given, the time step of the next cycle is reused for the
initial cycle.

If `blk` is given, its contribution is only added to the `state.global_dt` (the [`GlobalTimeStep`](@ref)).
Passing the whole block `grid` will block until the new time step is computed.

Return if the time step for the next cycle cannot be computed, and we must wait before doing so.
"""
function next_time_step(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    !need_to_update_time_step(params, state) && return false

    # We may need to wait for other blocks/sub-domains beforehand
    !can_contribute_to_local_time_step(params, state.global_dt, state.cycle) && return true

    # Compute this block's contribution to the next cycle's time step
    local_dt = local_time_step(params, state, blk)
    contribute_to_local_time_step!(params, state.global_dt, local_dt)

    # Update the time step for the current cycle
    state.dt = state.global_dt.current_dt

    return false
end


function fetch_time_step(params::ArmonParameters, state::SolverState, ::LocalTaskBlock)
    # The block already contributed, and we were waiting for a new time step. Is it available?
    !need_to_update_time_step(params, state) && return true
    if (@atomic state.global_dt.state.x) == TimeStepState.AllDone
        state.dt = state.global_dt.current_dt
        return true
    else
        return false
    end
end


function next_time_step(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    !need_to_update_time_step(params, state) && return false

    # Compute the contribution of all blocks to the next cycle's time step
    @section "local_time_step" begin
        local_dt = local_time_step(params, state, grid)
    end

    @section "time_step_reduction" begin
        contribute_to_local_time_step!(params, state.global_dt, local_dt; all_blocks=true)
    end

    if state.dt == 0
        # Force a blocking wait for the first cycle's time step
        wait_for_time_step!(params, state.global_dt)
    end

    # Update the time step for this cycle
    state.dt = state.global_dt.current_dt

    return false
end


@inline @fast function conservation_vars_kernel_reduction(ρ::T, E::T, mask::T) where T
    return (
        ρ * mask,     # Mass
        ρ * E * mask  # Energy
    )
end


@inline @fast function conservation_vars_kernel_reduction(ρ::T, E::T) where T
    # Mask-less version
    return (
        ρ,     # Mass
        ρ * E  # Energy
    )
end


@fast function conservation_vars(params::ArmonParameters{T, CPU_HP}, blk::LocalTaskBlock) where {T}
    # CPU reduction
    (; ρ, E) = block_device_data(blk)
    range = block_domain_range(blk.size, blk.state.steps_ranges[Int(state.axis)].real_domain)

    if params.use_cache_blocking
        # Reduction exploiting multithreading from the caller
        res_mass = zero(T)
        res_energy = zero(T)
        for j in range.col, i in range.row .+ (j - 1)
            (res_mass, res_energy) = (res_mass, res_energy) .+ conservation_vars_kernel_reduction(ρ[i], E[i])
        end
    else
        # Reduction using explicit multithreading, since the caller isn't multithreaded
        threads_mass   = Vector{T}(undef, params.nthreads)
        threads_energy = Vector{T}(undef, params.nthreads)
        threads_mass   .= 0
        threads_energy .= 0

        @threaded for j in range.col
            tid = Threads.threadid()
            thread_mass = threads_mass[tid]
            thread_energy = threads_energy[tid]
            for i in range.row .+ (j - 1)
                cell_mass, cell_energy = conservation_vars_kernel_reduction(ρ[i], E[i])
                thread_mass += cell_mass
                thread_energy += cell_energy
            end
            threads_mass[tid] = thread_mass
            threads_energy[tid] = thread_energy
        end

        (res_mass, res_energy) = sum(threads_mass), sum(threads_energy)
    end

    ds = prod(params.domain_size ./ params.global_grid)
    res_mass   *= ds
    res_energy *= ds

    return res_mass, res_energy
end


@generic_kernel function conservation_vars(
    ρ::V, E::V, res_mass::V, res_energy::V, bsize::BlockSize
) where {V}
    i = @index_2D_lin()
    mask = eltype(V)(!is_ghost(bsize, i))
    (res_mass[i], res_energy[i]) = conservation_vars_kernel_reduction(ρ[i], E[i], mask)
end


function conservation_vars(params::ArmonParameters{T}, blk::LocalTaskBlock) where {T}
    # GPU generic reduction
    range = block_domain_range(blk.size, blk.state.steps_ranges[Int(blk.state.axis)].full_domain)
    blk_data = block_device_data(blk)

    if params.use_two_step_reduction
        # Use a temporary array to store the partial reduction results.
        # Same comments as for `dtCFL_kernel`.
        conservation_vars(params, blk_data, range, blk_data.work_1, blk_data.work_2, blk.size)
        wait(params)
        identity2(a, b) = (a, b)
        total_mass, total_energy = mapreduce(identity2, min, blk_data.work_1, blk_data.work_2;
            init=(zero(T), zero(T)))
    else
        lin_range = first(range):last(range)
        ρ_v    = @view blk_data.ρ[lin_range]
        E_v    = @view blk_data.E[lin_range]
        mask_v = @view blk_data.mask[lin_range]
        total_mass, total_energy = mapreduce(conservation_vars_kernel_reduction, .+, ρ_v, E_v, mask_v;
            init=(zero(T), zero(T)))
    end

    ds = prod(params.domain_size ./ params.global_grid)
    total_mass   *= ds
    total_energy *= ds

    return total_mass, total_energy
end


function conservation_vars(params::ArmonParameters{T}, grid::BlockGrid) where {T}
    mt_reduction = params.use_threading && params.use_cache_blocking
    threads_mass   = Vector{T}(undef, mt_reduction ? params.nthreads : 1)
    threads_energy = Vector{T}(undef, mt_reduction ? params.nthreads : 1)
    threads_mass   .= 0
    threads_energy .= 0

    # Sum the energy and mass over all blocks of each thread
    @iter_blocks for blk in grid
        tid = mt_reduction ? Threads.threadid() : 1
        (threads_mass[tid], threads_energy[tid]) =
            (threads_mass[tid], threads_energy[tid]) .+ conservation_vars(params, blk)
    end

    # Local reduction over all threads
    total_mass   = sum(threads_mass)
    total_energy = sum(threads_energy)

    # Global reduction over all MPI processes
    # TODO: this is wrong! we MUST ensure that the right thread associated with the communicator
    #  of `params.reduc_model` is doing the reduction
    #  => how do we run a task on a specific thread in Julia ?
    if params.use_MPI
        global_reduction = Communications.init_reduce_broadcast(params.reduc_model, +, Vector{T}, 2)

        send_buf = Communications.acquire_send_buffer!(global_reduction)
        send_buf .= (total_mass, total_energy)
        Communications.release_send_buffer!(global_reduction)

        recv_buf = Communications.acquire_recv_buffer!(global_reduction)
        (total_mass, total_energy) = recv_buf
        Communications.release_recv_buffer!(global_reduction)
    end

    return total_mass, total_energy
end
