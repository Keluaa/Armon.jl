
@enumx TimeStepState::UInt32 begin
    "`current_dt` is up-to-date and blocks can contribute to `next_dt`"
    Ready
    "All blocks contributed to `next_dt`, one thread will start the MPI reduction"
    AllContributed
    "The MPI reduction has started"
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

    function GlobalTimeStep{T}() where {T}
        return new{T}(
            Atomic(TimeStepState.Ready),
            0, zero(T),
            zero(T), typemax(T), Atomic(typemax(T)),
            Atomic(0), 0,
            MPI.Request(), MPI.RBuffer(Ref{T}(), Ref{T}())
        )
    end
end


time_step_state(global_dt::GlobalTimeStep) = @atomic global_dt.state.x
time_step_state!(global_dt::GlobalTimeStep, state::TimeStepState.T) = @atomic global_dt.state.x = state
function replace_time_step_state!(global_dt::GlobalTimeStep, transition::Pair{TimeStepState.T, TimeStepState.T})
    _, ok = @atomicreplace global_dt.state.x transition
    return ok
end


function reset!(global_dt::GlobalTimeStep{T}, params::ArmonParameters{T}, block_count) where {T}
    time_step_state!(global_dt, TimeStepState.Ready)
    global_dt.cycle = 0
    global_dt.time = zero(T)
    global_dt.current_dt = params.cst_dt ? params.Dt : zero(T)
    global_dt.next_cycle_dt = typemax(T)
    @atomic global_dt.next_dt.x = typemax(T)
    @atomic global_dt.contributions.x = 0
    global_dt.expected_count = block_count
end


function contribute_to_dt!(params::ArmonParameters, global_dt::GlobalTimeStep{T}, dt::T; all_blocks=false) where {T}
    @atomic global_dt.next_dt.x min dt  # Atomic reduction

    contributed_blocks = all_blocks ? global_dt.expected_count : 1
    contributions = @atomic global_dt.contributions.x += contributed_blocks
    if contributions == global_dt.expected_count
        if !replace_time_step_state!(global_dt, TimeStepState.Ready => TimeStepState.AllContributed)
            return TimeStepState.AllContributed
        end

        # All blocks have contributed, therefore `global_dt.current_dt` is outdated: we can update
        # it safely.
        return update_dt!(params, global_dt)
    else
        return TimeStepState.Ready
    end
end


function wait_for_dt!(params::ArmonParameters, global_dt::GlobalTimeStep)
    if !replace_time_step_state!(global_dt, TimeStepState.DoingMPI => TimeStepState.WaitingForMPI)
        return TimeStepState.WaitingForMPI
    end

    # Since this thread started working on a block without the time step for the new cycle, we
    # consider that all blocks of that thread are in the same state, therefore loosing no time by
    # using a blocking wait here. Only a single thread will wait.
    params.use_MPI && wait(global_dt.MPI_reduction)
    return update_dt!(params, global_dt)
end


function update_dt!(params::ArmonParameters, global_dt::GlobalTimeStep{T}) where {T}
    state = time_step_state(global_dt)
    if state == TimeStepState.AllContributed
        local_dt = @atomicswap global_dt.next_dt.x = typemax(T)
        if params.use_MPI
            global_dt.MPI_buffer.senddata[] = local_dt
            global_dt.MPI_buffer.recvdata[] = typemax(T)
            IAllreduce!(global_dt.MPI_buffer, MPI.MIN, params.cart_comm, global_dt.MPI_reduction)
            time_step_state!(global_dt, TimeStepState.DoingMPI)
            return TimeStepState.DoingMPI
        else
            new_dt = local_dt
        end
    elseif state == TimeStepState.WaitingForMPI
        new_dt = global_dt.MPI_buffer.recvdata[]
    else
        error("unexpected time step state: $state")
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
    time_step_state!(global_dt, TimeStepState.Done)
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

    time_step_state!(global_dt, TimeStepState.Ready)
    global_dt.current_dt = global_dt.next_cycle_dt
    global_dt.next_cycle_dt = typemax(T)
end


"""
    SolverStep

Enumeration of each state a [`LocalTaskBlock`](@ref) can be in.
[`block_state_machine`](@ref) advances this state.
"""
@enumx SolverStep begin
    NewCycle
    TimeStep
    InitTimeStep
    NewSweep
    EOS
    Exchange
    Fluxes
    CellUpdate
    Remap
    EndCycle
    ErrorState
end


"""
    SolverState

Object containing all non-constant parameters needed to run the solver, as well as type-parameters
needed to avoid runtime dispatch.

This object is local to a block (or set of blocks): multiple blocks could be at different steps of
the solver at once.
"""
mutable struct SolverState{T, Splitting, Riemann, RiemannLimiter, Projection, TestCase}
    step               :: SolverStep.T  # Solver step the associated block is at. Unused if `params.async_cycle == false`
    dx                 :: T    # Space step along the current axis
    dt                 :: T    # Scaled time step for the current cycle
    axis               :: Axis
    axis_splitting_idx :: Int
    cycle              :: Int  # Local cycle of the block
    splitting          :: Splitting
    riemann_scheme     :: Riemann
    riemann_limiter    :: RiemannLimiter
    projection_scheme  :: Projection
    test_case          :: TestCase
    global_dt          :: GlobalTimeStep{T}
    steps_ranges       :: StepsRanges

    function SolverState{T}(
        splitting::S, riemann::R, limiter::RL, projection::P, test_case::TC, global_dt, steps_ranges
    ) where {
        T, S <: SplittingMethod, R <: RiemannScheme, RL <: Limiter, P <: ProjectionScheme, TC <: TestCase
    }
        return new{T, S, R, RL, P, TC}(
            SolverStep.NewCycle, zero(T), zero(T), X_axis, 1, 0,
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


function start_cycle(state::SolverState)
    # If `cycle > global_dt.cycle` then we must wait for the other blocks to finish the previous cycle.
    return state.cycle == state.global_dt.cycle
end


function end_cycle!(state::SolverState)
    state.cycle += 1
end


solver_step(state::SolverState) = state.step
finished_cycle(state::SolverState) = state.cycle == state.global_dt.cycle && state.step == SolverStep.NewCycle


function reset!(state::SolverState{T}) where {T}
    state.step = SolverStep.NewCycle
    state.dx = zero(T)
    state.dt = zero(T)
    state.axis = X_axis
    state.axis_splitting_idx = 1
    state.cycle = 0
end
