
struct StateMachineKernelError <: Exception
    task :: SolverStep.T
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


@kernel cpu=false function state_machine_kernel(data::BlockData, state::BasicSolverState, bsize::BlockSize, sub_tasks::AbstractArray{SolverStep.T})

    # TODO: first attempt: kernels are launched with as many items as there is cells in the block
    #   => this means that we don't have to worry about tiling
    #   => this means that all items on the edges will do almost no work
    #   => this means that items are not a multiple of the GPU block size

    # TODO: UVM? Host memory accesses? since all memory accesses will be cached for a long time, they might be viable

    I = @index(Global, NTuple)
    idx = (;
        idx = @index(Global, Linear),
        lin_1D = 0,  # TODO: not sure
        lin_2D = lin_position(bsize, I)  # TODO: check
    )

    for task in sub_tasks
        # Each step may use a different domain, therefore some threads will have nothing to do for
        # some steps. They go immediately to the `__syncthreads` call.
        if !is_gpu_thread_in_step_domain(state, bisze, task, I)
            @goto nothing_to_do
        end

        if task == SolverStep.TimeStep
            # TODO: reduction? => possible for a single block, need to reuse the logic inside CUDA.mapreduce/AMDGPU.mapreduce
            # TODO: use `Armon.gpu_workgroup_reduction` for this
            # TODO: are we forced to use a single workgroup? couldn't we just place the results of
            # each workgroup in an array (which wouldn't be very large), and then perform the final
            # reduction on the CPU after a copy. Since the reduction isn't always required for the
            # rest of the steps, this is a very attractive option.
            # => however this comes at the cost of preventing any other type of reduction in this kernel,
            #    this can be viewed as an ad-hoc optimization...
            # => the only alternative is parsing the block by tiles of the workgroup size
            throw(StateMachineKernelError(task))

        elseif task == SolverStep.NewSweep
            # TODO: start a new sweep here? or do it on the host?
            # TODO: since we have the domains, strides, and all data arrays, it is possible to do this
            # TODO: but we may be missing Δx or Δy
            throw(StateMachineKernelError(task))

        elseif task == SolverStep.EOS
            if state.schemes.test_case <: Bizarrium
                @sub_kernel_call idx bizarrium_EOS!(data.ρ, data.u, data.v, data.E, data.p, data.c, data.g)
            else
                γ = eltype(state)(specific_heat_ratio(state.schemes.test_case))
                @sub_kernel_call idx perfect_gas_EOS!(γ, data.ρ, data.E, data.u, data.v, data.p, data.c, data.g)
            end

        elseif task == SolverStep.Exchange
            throw(StateMachineKernelError(task))

            # TODO: boundary conditions use a very different domain, dertermine if `idx` is in the domain
            #   and compute the BC only if so, otherwise `@synchronize(is_idx_not_on_the_side)`
            # TODO: since there is one kernel item per cell, we could do the BC of both sides at once
            # TODO: u_factor, v_factor could also be computed from here
            side = Side.Left
            u_factor, v_factor = zero(eltype(state)), zero(eltype(state))
            @sub_kernel_call idx boundary_conditions!(
                data.ρ, data.u, data.v, data.p, data.c, data.g, data.E,
                bsize, state.axis, side,
                u_factor, v_factor
            )
            # TODO: local block exchange using atomic (but how to get the dimensions + data arrays of the neighbouring block?)

        elseif task == SolverStep.Fluxes
            s = stride_along(bsize, state.axis)
            uₐ = state.axis == Axis.X ? data.u : data.v
            if state.schemes.riemann_scheme <: RiemannGodunov
                @sub_kernel_call idx acoustic!(s, data.uˢ, data.pˢ, data.ρ, uₐ, data.p, data.c)
            elseif state.schemes.riemann_scheme <: RiemannGAD
                @sub_kernel_call idx acoustic_GAD!(
                    s, state.dt, state.dx,
                    data.uˢ, data.pˢ, data.ρ, uₐ, data.p, data.c,
                    state.schemes.riemann_limiter
                )
            else
                throw(StateMachineKernelError(task))
            end

        elseif task == SolverStep.CellUpdate
            s = stride_along(bsize, state.axis)
            uₐ = state.axis == Axis.X ? data.u : data.v
            @sub_kernel_call idx cell_update!(s, state.dx, state.dt, data.uˢ, data.pˢ, data.ρ, uₐ, data.E)

        elseif step == SolverStep.RemapAdvection
            s = stride_along(bsize, state.axis)
            if state.schemes.projection_scheme <: EulerProjection
                @sub_kernel_call idx advection_first_order!(
                    s, state.dt,
                    data.uˢ, data.ρ, data.u, data.v, data.E,
                    data.work_1, data.work_2, data.work_3, data.work_4
                )
            elseif state.schemes.projection_scheme <: Euler2ndProjection
                @sub_kernel_call idx advection_second_order!(
                    s, state.dx, state.dt,
                    data.uˢ, data.ρ, data.u, data.v, data.E,
                    data.work_1, data.work_2, data.work_3, data.work_4
                )
            else
                throw(StateMachineKernelError(task))
            end

        elseif task == SolverStep.RemapProjection
            s = stride_along(bsize, state.axis)
            @sub_kernel_call idx euler_projection!(
                s, state.dx, state.dt, data.uˢ, data.ρ, data.u, data.v, data.E,
                data.work_1, data.work_2, data.work_3, data.work_4
            )

        elseif task == SolverStep.EndCycle
            break

        else
            throw(StateMachineKernelError(task))
        end

@label nothing_to_do
        KernelAbstractions.@synchronize()
    end
end
