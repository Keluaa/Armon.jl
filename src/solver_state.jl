
@enumx TimeStepState::UInt32 begin
    "`current_dt` is up-to-date and blocks can contribute to `next_dt`"
    Ready
    "All blocks contributed to `next_dt`, the MPI reduction has started"
    DoingMPI
    "One thread is waiting for the MPI reduction to complete"
    WaitingForMPI
    "MPI is done: `current_dt` is up-to-date"
    Done
end


"""
    GlobalTimeStep

Holds all information about the current time and time step for the current solver cycle. This struct
is global and shared among all blocks.

When reaching `next_time_step`, blocks will contribute to the calculation of the next time step. The
last block doing so will start the MPI reduction. The first block reaching the start of the next
cycle will wait until this reduction is completed, updating the `GlobalTimeStep` when so.
"""
mutable struct GlobalTimeStep{T}
    state          :: Atomic{TimeStepState.T}
    cycle          :: Int
    time           :: T
    current_dt     :: T
    next_cycle_dt  :: T  # Result of the reduction for `next_dt`. `Inf` if not ready.
    next_dt        :: Atomic{T}  # Time step accumulator
    contributions  :: Atomic{Int}
    expected_count :: Int
    MPI_reduction  :: MPI.AbstractRequest
    MPI_buffer     :: MPI.RBuffer
end


function GlobalTimeStep(params::ArmonParameters{T}, grid::BlockGrid) where {T}
    init_dt = params.cst_dt ? params.Dt : zero(T)
    blocks_count = prod(grid.grid_size)
    return GlobalTimeStep{T}(
        Atomic(TimeStepState.Ready),
        0, zero(T),
        init_dt, typemax(T), Atomic(typemax(T)),
        Atomic(0), blocks_count,
        MPI.Request(), MPI.RBuffer(Ref{T}(), Ref{T}())
    )
end


time_step_state(global_dt::GlobalTimeStep) = @atomic global_dt.state.x


function contribute_to_dt!(params::ArmonParameters, global_dt::GlobalTimeStep{T}, dt::T; all_blocks=false) where {T}
    @atomic global_dt.next_dt.x min dt  # Atomic reduction

    contributed_blocks = all_blocks ? global_dt.expected_count : 1
    contributions = @atomic global_dt.contributions.x += contributed_blocks
    if contributions == global_dt.expected_count
        _, ok = @atomicreplace global_dt.state.x TimeStepState.Ready => TimeStepState.DoingMPI
        !ok && return TimeStepState.DoingMPI
        # All blocks have contributed, therefore `global_dt.current_dt` is outdated: we can update
        # it safely.
        return update_dt!(params, global_dt)
    else
        return TimeStepState.Ready
    end
end


function wait_for_dt!(params::ArmonParameters, global_dt::GlobalTimeStep)
    _, ok = @atomicreplace global_dt.state.x TimeStepState.DoingMPI => TimeStepState.WaitingForMPI
    !ok && return TimeStepState.WaitingForMPI

    # Since this thread started working on a block without the time step for the new cycle, we
    # consider that all blocks of that thread are in the same state, therefore loosing no time by
    # using a blocking wait here. Only a single thread will wait.
    params.use_MPI && wait(global_dt.MPI_reduction)
    return update_dt!(params, global_dt, true)
end


function update_dt!(params::ArmonParameters, global_dt::GlobalTimeStep{T}, reduction_done=false) where {T}
    if params.use_MPI
        if !reduction_done
            local_dt = @atomicswap global_dt.next_dt.x = typemax(T)
            global_dt.MPI_buffer.senddata[] = local_dt
            global_dt.MPI_buffer.recvdata[] = zero(T)
            IAllreduce!(global_dt.MPI_buffer, MPI.Op(min, T), params.cart_comm, global_dt.MPI_reduction)
            return TimeStepState.DoingMPI
        else
            new_dt = global_dt.MPI_buffer.recvdata[]
        end
    else
        new_dt = @atomicswap global_dt.next_dt.x = typemax(T)
    end

    previous_dt = global_dt.current_dt

    if (!isfinite(new_dt) || new_dt ≤ 0)
        solver_error(:time, "Invalid time step for cycle $(global_dt.cycle): $new_dt")
    elseif previous_dt == 0
        new_dt = params.cfl * new_dt
    else
        # CFL condition and maximum increase per cycle of the time step
        new_dt = convert(T, min(params.cfl * new_dt, 1.05 * previous_dt))
    end

    global_dt.next_cycle_dt = new_dt

    if global_dt.current_dt == 0
        # The current time step needs to be initialized
        global_dt.current_dt = global_dt.next_cycle_dt
    end

    @atomic global_dt.contributions.x = 0
    @atomic global_dt.state.x = TimeStepState.Done
    return TimeStepState.Done
end


function next_cycle!(params::ArmonParameters, global_dt::GlobalTimeStep{T}) where {T}
    global_dt.cycle += 1
    global_dt.time += global_dt.current_dt

    if time_step_state(global_dt) == TimeStepState.DoingMPI
        wait_for_dt!(params, global_dt)
    end

    dt_state = time_step_state(global_dt)
    if dt_state != TimeStepState.Done
        error("expected time step to be done, got: $dt_state")
    end

    @atomic global_dt.state.x = TimeStepState.Ready
    global_dt.current_dt = global_dt.next_cycle_dt
    global_dt.next_cycle_dt = typemax(T)
end


"""
    SolverState

Object containing all non-constant parameters needed to run the solver, as well as type-parameters
needed to avoid runtime dispatch.

This object is local to a block (or set of blocks): multiple blocks could be at different steps of
the solver at once.
"""
mutable struct SolverState{T, Splitting, Riemann, RiemannLimiter, Projection, TestCase}
    dx::T  # Space step along the current axis
    dt::T  # Scaled time step for the current cycle
    axis::Axis
    axis_splitting_idx::Int
    splitting::Splitting
    riemann_scheme::Riemann
    riemann_limiter::RiemannLimiter
    projection_scheme::Projection
    test_case::TestCase
    global_dt::GlobalTimeStep{T}
    steps_ranges::StepsRanges

    function SolverState{T}(
        splitting::S, riemann::R, limiter::RL, projection::P, test_case::TC, global_dt, steps_ranges
    ) where {
        T, S <: SplittingMethod, R <: RiemannScheme, RL <: Limiter, P <: ProjectionScheme, TC <: TestCase
    }
        return new{T, S, R, RL, P, TC}(
            zero(T), zero(T), X_axis, 1,
            splitting, riemann, limiter, projection, test_case,
            global_dt, steps_ranges
        )
    end
end


function SolverState(params::ArmonParameters{T}, global_dt::GlobalTimeStep{T}) where {T}
    return SolverState{T}(
        params.axis_splitting,
        params.riemann_scheme, params.riemann_limiter,
        params.projection_scheme,
        params.test,
        global_dt,
        first(params.steps_ranges)
    )
end


function next_axis_sweep!(params::ArmonParameters, state::SolverState)
    if state.axis_splitting_idx == 0
        iter_val = iterate(split_axes(state))
    else
        iter_val = iterate(split_axes(state), state.axis_splitting_idx)
    end

    if isnothing(iter_val)
        state.axis_splitting_idx = 0
        return true
    else
        ((axis, dt_factor), state.axis_splitting_idx) = iter_val
        update_solver_state!(params, state, axis, dt_factor)
        return false
    end
end


function update_solver_state!(params::ArmonParameters, state::SolverState, axis::Axis, dt_factor)
    i_ax = Int(axis)
    state.dx = params.domain_size[i_ax] / params.global_grid[i_ax]
    state.dt = state.global_dt.current_dt * dt_factor
    state.axis = axis
    state.steps_ranges = params.steps_ranges[i_ax]
end
