
struct StateMachineKernelError <: Exception
    step :: SolverStep.T
end

function Base.showerror(io::IO, ex::StateMachineKernelError)
    println(io, "Invalid device state machine step: ", ex.step)
end


function is_gpu_thread_in_step_domain(state::BasicSolverState, bsize::BlockSize, step::SolverStep.T, I::NTuple{2})
    steps_ranges = state.steps_ranges[Int(state.axis)]

    if     step == SolverStep.TimeStep         return false  # TODO: full_domain, or real_domain? We can nicely optimize things here
    elseif step == SolverStep.NewSweep         return false
    elseif step == SolverStep.EOS              corners = steps_ranges.EOS
    elseif step == SolverStep.Exchange         return false  # TODO
    elseif step == SolverStep.Fluxes           corners = steps_ranges.fluxes
    elseif step == SolverStep.CellUpdate       corners = steps_ranges.cell_update
    elseif step == SolverStep.RemapAdvection   corners = steps_ranges.advection
    elseif step == SolverStep.RemapProjection  corners = steps_ranges.projection
    elseif step == SolverStep.EndCycle         return false
    else                                       return false
    end

    # A domain is represented by offsets to the bottom-left and top-right corners of the block size.
    # By some basic index shifting, we can guarentee that one cell will only be worked on by a single
    # kernel item/thread.
    # Kernel items outside of the domain have no work to do for that step.
    return all(corners[1] .≤ I .- ghosts(bsize) .- 1 .< real_block_size(bsize) .+ corners[2])
end


@kernel cpu=false function state_machine_kernel(data::BlockData, state::BasicSolverState, bsize::BlockSize, queue::DeviceStepQueue)
    # TODO: first attempt: kernels are launched with as many items as there is cells in the block
    #   => this means that we don't have to worry about tiling
    #   => this means that all items on the edges will do almost no work
    #   => this means that items are not a multiple of the GPU block size

    # TODO: UVM? Host memory accesses? since all memory accesses will be cached for a long time, they might be viable

    I = @index(Global, NTuple)
    idx = (;
        idx = @index(Global, Linear),
        lin_1D = 0,  # TODO: not sure, only used in some steps
        lin_2D = lin_position(bsize, I .- ghosts(bsize))
    )

    step_idx = @inbounds queue.status[2]
    KernelAbstractions.@synchronize()

@label step_loop
    # In order to spare some registers, we don't use a `for` loop but use `@goto` and an index stored
    # in `queue.status[2]`.
    step_idx > length(queue) && return
    step = queue.steps[step_idx]

    if !can_run_on_device(queue.device, step)
        throw(StateMachineKernelError(step))
    end

    # Each step may use a different domain, therefore some threads will have nothing to do for
    # some steps. They go immediately to the `__syncthreads` call.
    if !is_gpu_thread_in_step_domain(state, bsize, step, I)
        @goto nothing_to_do
    end

    if step == SolverStep.TimeStep
        # TODO: reduction? => possible for a single block, need to reuse the logic inside CUDA.mapreduce/AMDGPU.mapreduce
        # TODO: use `Armon.gpu_workgroup_reduction` for this
        # TODO: are we forced to use a single workgroup? couldn't we just place the results of
        # each workgroup in an array (which wouldn't be very large), and then perform the final
        # reduction on the CPU after a copy. Since the reduction isn't always required for the
        # rest of the steps, this is a very attractive option.
        # => however this comes at the cost of preventing any other type of reduction in this kernel,
        #    this can be viewed as an ad-hoc optimization...
        # => the only alternative is parsing the block by tiles of the workgroup size

    elseif step == SolverStep.NewSweep
        # TODO: start a new sweep here? or do it on the host?
        # TODO: since we have the domains, strides, and all data arrays, it is possible to do this
        # TODO: but we may be missing Δx or Δy

    elseif step == SolverStep.EOS
        if state.schemes.test_case isa Bizarrium
            @sub_kernel_call idx bizarrium_EOS!(data.ρ, data.u, data.v, data.E, data.p, data.c, data.g)
        else
            γ = eltype(state)(specific_heat_ratio(state.schemes.test_case))
            @sub_kernel_call idx perfect_gas_EOS!(γ, data.ρ, data.E, data.u, data.v, data.p, data.c, data.g)
        end

    elseif step == SolverStep.Exchange
        # TODO: boundary conditions use a very different domain, dertermine if `idx` is in the domain
        #   and compute the BC only if so, otherwise `@synchronize(is_idx_not_on_the_side)`
        # TODO: since there is one kernel item per cell, we could do the BC of both sides at once
        # TODO: u_factor, v_factor could also be computed from here
        # side = Side.Left
        # u_factor, v_factor = zero(eltype(state)), zero(eltype(state))
        # @sub_kernel_call idx boundary_conditions!(
        #     data.ρ, data.u, data.v, data.p, data.c, data.g, data.E,
        #     bsize, state.axis, side,
        #     u_factor, v_factor
        # )
        # TODO: local block exchange using atomic (but how to get the dimensions + data arrays of the neighbouring block?)

    elseif step == SolverStep.Fluxes
        s = stride_along(bsize, state.axis)
        uₐ = state.axis == Axis.X ? data.u : data.v
        if state.schemes.riemann_scheme isa RiemannGodunov
            @sub_kernel_call idx acoustic!(s, data.uˢ, data.pˢ, data.ρ, uₐ, data.p, data.c)
        elseif state.schemes.riemann_scheme isa RiemannGAD
            @sub_kernel_call idx acoustic_GAD!(
                s, state.dt, state.dx,
                data.uˢ, data.pˢ, data.ρ, uₐ, data.p, data.c,
                state.schemes.riemann_limiter
            )
        end

    elseif step == SolverStep.CellUpdate
        s = stride_along(bsize, state.axis)
        uₐ = state.axis == Axis.X ? data.u : data.v
        @sub_kernel_call idx cell_update!(s, state.dx, state.dt, data.uˢ, data.pˢ, data.ρ, uₐ, data.E)

    elseif step == SolverStep.RemapAdvection
        s = stride_along(bsize, state.axis)
        if state.schemes.projection_scheme isa EulerProjection
            @sub_kernel_call idx advection_first_order!(
                s, state.dt,
                data.uˢ, data.ρ, data.u, data.v, data.E,
                data.work_1, data.work_2, data.work_3, data.work_4
            )
        elseif state.schemes.projection_scheme isa Euler2ndProjection
            @sub_kernel_call idx advection_second_order!(
                s, state.dx, state.dt,
                data.uˢ, data.ρ, data.u, data.v, data.E,
                data.work_1, data.work_2, data.work_3, data.work_4
            )
        end

    elseif step == SolverStep.RemapProjection
        s = stride_along(bsize, state.axis)
        @sub_kernel_call idx euler_projection!(
            s, state.dx, state.dt, data.uˢ, data.ρ, data.u, data.v, data.E,
            data.work_1, data.work_2, data.work_3, data.work_4
        )
    end

@label nothing_to_do
    if idx.idx == 1  # TODO: ugly
        queue.status[2] += 1  # increment the position, as this step was completed
        # queue.status[3] = true  # we can always continue to the next step (for now)
    end
    KernelAbstractions.@synchronize()

    # step_idx = queue.status[2]  # TODO: using atomic add + load is mandatory here I think
    step_idx += 1
    @goto step_loop
end


function thread_position_to_block_index(thread_pos, tile_iter_idx, state, wrap_tile, step)
    tile_idx = thread_pos .+ tile_iter_idx .* wrap_tile

    # `corners` are offsets, and `I` would be an index in the real cells of the block
    corners = getfield(state.steps_ranges[Int(state.axis)], step)
    I = tile_idx .+ corners[1]

    in_bounds = all(I .< real_block_size(bsize) .+ corners[2])
    return I, in_bounds
end


macro tiled_2D_iter(step, step_call)
    return esc(quote
        for tile_iter_idx_y in 1:tile_count[2], tile_iter_idx_x in 1:tile_count[1]
            # Compute everything from `tile_iter_idx`, in order to minimize the amount of
            # memory dependancies across loop iterations.
            thread_pos = @index(Local, NTuple)
            I, in_bounds = thread_position_to_block_index(thread_pos, (tile_iter_idx_x, tile_iter_idx_y), state, wrap_tile, $step)
            if in_bounds
                idx = (;
                    idx = 0,
                    lin_1D = 0,
                    lin_2D = lin_position(bsize, I)
                )
                $step_call
            end
        end
    end)
end


@kernel cpu=false function tiled_block_iter(
    data::BlockData, state::BasicSolverState, queue::DeviceStepQueue, bsize::BlockSize,
    ::Val{wrap_tile}, ::Val{tile_count}
) where {wrap_tile, tile_count}
    # TODO: does using gotos instead of loops could improve register usage and performance?
    for step_idx in 1:length(queue)
        step = queue.steps[step_idx]

        if step == SolverStep.EOS
            @tiled_2D_iter :EOS begin
                if state.schemes.test_case isa Bizarrium
                    @sub_kernel_call idx bizarrium_EOS!(data.ρ, data.u, data.v, data.E, data.p, data.c, data.g)
                else
                    γ = eltype(state)(specific_heat_ratio(state.schemes.test_case))
                    @sub_kernel_call idx perfect_gas_EOS!(γ, data.ρ, data.E, data.u, data.v, data.p, data.c, data.g)
                end
            end

        elseif step == SolverStep.Fluxes
            @tiled_2D_iter :fluxes begin
                s = stride_along(bsize, state.axis)
                uₐ = state.axis == Axis.X ? data.u : data.v
                if state.schemes.riemann_scheme isa RiemannGodunov
                    @sub_kernel_call idx acoustic!(s, data.uˢ, data.pˢ, data.ρ, uₐ, data.p, data.c)
                elseif state.schemes.riemann_scheme isa RiemannGAD
                    @sub_kernel_call idx acoustic_GAD!(
                        s, state.dt, state.dx,
                        data.uˢ, data.pˢ, data.ρ, uₐ, data.p, data.c,
                        state.schemes.riemann_limiter
                    )
                end
            end

        elseif step == SolverStep.CellUpdate
            @tiled_2D_iter :cell_updsate begin
                s = stride_along(bsize, state.axis)
                uₐ = state.axis == Axis.X ? data.u : data.v
                @sub_kernel_call idx cell_update!(s, state.dx, state.dt, data.uˢ, data.pˢ, data.ρ, uₐ, data.E)
            end
    
        elseif step == SolverStep.RemapAdvection
            @tiled_2D_iter :advection begin
                s = stride_along(bsize, state.axis)
                if state.schemes.projection_scheme isa EulerProjection
                    @sub_kernel_call idx advection_first_order!(
                        s, state.dt,
                        data.uˢ, data.ρ, data.u, data.v, data.E,
                        data.work_1, data.work_2, data.work_3, data.work_4
                    )
                elseif state.schemes.projection_scheme isa Euler2ndProjection
                    @sub_kernel_call idx advection_second_order!(
                        s, state.dx, state.dt,
                        data.uˢ, data.ρ, data.u, data.v, data.E,
                        data.work_1, data.work_2, data.work_3, data.work_4
                    )
                end
            end
    
        elseif step == SolverStep.RemapProjection
            @tiled_2D_iter :projection begin
                s = stride_along(bsize, state.axis)
                @sub_kernel_call idx euler_projection!(
                    s, state.dx, state.dt, data.uˢ, data.ρ, data.u, data.v, data.E,
                    data.work_1, data.work_2, data.work_3, data.work_4
                )
            end
        end

        @synchronize()
    end

    if @index(Global, Cartesian) == CartesianIndex(1, 1)
        queue.status[2] = length(queue) + 1
    end
end


function process_queue!(queue::DeviceStepQueue, params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    basic_state = BasicSolverState(state)

    if !params.use_tiled_state_machine
        state_machine_func = state_machine_kernel(queue.device, params.workgroup_size)
        # TODO: the ndrange is quite important, the current choice is maybe sub-optimal since it includes all ghost cells
        state_machine_func(blk.device_data, basic_state, blk.size, queue; ndrange=block_size(blk))
    else
        # grid size == workgroup size  =>  exactly 1 workgroup per kernel
        tiled_block_iter_func = tiled_block_iter(queue.device, params.workgroup_size, params.workgroup_size)

        # split `block_size(blk)` into tiles
        # `(64, 64)` split into tiles of `workgroup_size`
        # => `(64, 64) .÷ (32, 32) = (2, 2)`
        # => each thread does:
        #  - `thread_pos .+ work_group_size .* (0, 0)`
        #  - `thread_pos .+ work_group_size .* (1, 0)`
        #  - `thread_pos .+ work_group_size .* (0, 0)`
        #  - `thread_pos .+ work_group_size .* (1, 1)`

        wrap_tile = params.workgroup_size
        tile_count = cld.(block_size(blk), wrap_tile)

        tiled_block_iter_func(
            blk.device_data, basic_state, queue, blk.size,
            Val(wrap_tile), Val(tile_count)
        )
    end

    # Place an event in the device stream in order to be able to know when the kernel has completed,
    # independantly of the status of the stream.
    put_kernel_event(queue.device, queue.event)
    return
end
