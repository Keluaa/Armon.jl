
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

    if params.cst_dt
        global_dt.current_dt = global_dt.next_cycle_dt = params.Dt
        return
    end

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
@enumx SolverStep::UInt8 begin
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


# Bit flags for each of the variables of a block
const STEPS_VARS_FLAGS = (;
    x      = UInt16(1) << 0,
    ρ      = UInt16(1) << 1,
    u      = UInt16(1) << 2,
    E      = UInt16(1) << 3,
    p      = UInt16(1) << 4,
    c      = UInt16(1) << 5,
    g      = UInt16(1) << 6,
    uˢ     = UInt16(1) << 7,
    pˢ     = UInt16(1) << 8,
    work_1 = UInt16(1) << 9,
    work_2 = UInt16(1) << 10,
    work_3 = UInt16(1) << 11,
    mask   = UInt16(1) << 12,
)


# OR of all dimensional vars' flags
steps_dimensional_vars_flags() = |(getindex.(Ref(STEPS_VARS_FLAGS), dim_vars())...)


# Flags for the arrays used by kernels: `count_ones(SOLVER_STEPS_VARS[step][2])` represents the
# number of arrays the kernels of `step` can bring into the cache. If `SOLVER_STEPS_VARS[step][1] == true`
# then the kernel uses only one of the array of all dimensional variables.
# TODO: deduce them from kernel + steps definitions?
const SOLVER_STEPS_VARS = let SVF = STEPS_VARS_FLAGS
    Dict{SolverStep.T, Tuple{Bool, UInt16}}(
        SolverStep.NewCycle     => (false, 0),
        SolverStep.TimeStep     => (false, SVF.u | SVF.c),
        SolverStep.InitTimeStep => (false, SVF.u | SVF.c),
        SolverStep.NewSweep     => (false, 0),
        SolverStep.EOS          => (false, SVF.ρ | SVF.E | SVF.u | SVF.p | SVF.c | SVF.g),
        SolverStep.Exchange     => (false, SVF.ρ | SVF.E | SVF.u | SVF.p | SVF.c | SVF.g),
        SolverStep.Fluxes       => (true,  SVF.ρ | SVF.p | SVF.c | SVF.uˢ| SVF.pˢ),
        SolverStep.CellUpdate   => (true,  SVF.ρ | SVF.E | SVF.uˢ| SVF.pˢ),
        SolverStep.Remap        => (false, SVF.ρ | SVF.E | SVF.u | SVF.uˢ| SVF.work_1 | SVF.work_2 | SVF.work_3),
        SolverStep.EndCycle     => (false, 0),
        SolverStep.ErrorState   => (false, 0),
    )
end


"""
    ThreadLogEvent

Info about a thread, emitted after a call to `solver_cycle_async`, which itself calls [`block_state_machine`](@ref).
"""
struct ThreadLogEvent
    cycle             :: Int16
    blk_count         :: Int16
    stop_count        :: Int16
    mpi_waits         :: Int16
    step_count        :: Int32
    no_progress_count :: Int32
    wait_time         :: Float64
    cycle_time        :: Float64
end


"""
    BlockLogEvent

Info about a block, emitted after a call to [`block_state_machine`](@ref) which successfully advanced
the internal state of the block, only if `params.log_blocks == true`.
"""
struct BlockLogEvent
    cycle           :: Int16   # Cycle at which the event occured
    tid             :: UInt16  # Thread which processed the block
    axis            :: Axis.T  # Final axis of the block
    new_state       :: SolverStep.T  # Final state of the block
    steps_count     :: UInt8   # Number of solver steps done
    steps_vars      :: UInt16  # Flag with a 1 when a variable was used
    steps_var_count :: Int16   # Number of times all variables were used
    tid_blk_idx     :: Int32   # Number of blocks processed by the thread before this event
    stalls          :: Int64   # Number of times `block_state_machine` was called without completing a step
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
    axis               :: Axis.T
    axis_splitting_idx :: Int
    cycle              :: Int  # Local cycle of the block
    splitting          :: Splitting
    riemann_scheme     :: Riemann
    riemann_limiter    :: RiemannLimiter
    projection_scheme  :: Projection
    test_case          :: TestCase
    global_dt          :: GlobalTimeStep{T}
    steps_ranges       :: StepsRanges
    blk_logs           :: Vector{BlockLogEvent}
    total_stalls       :: Int

    function SolverState{T}(
        splitting::S, riemann::R, limiter::RL, projection::P, test_case::TC, global_dt, steps_ranges, log_size
    ) where {
        T, S <: SplittingMethod, R <: RiemannScheme, RL <: Limiter, P <: ProjectionScheme, TC <: TestCase
    }
        blk_logs = Vector{BlockLogEvent}()
        log_size > 0 && sizehint!(blk_logs, log_size)
        return new{T, S, R, RL, P, TC}(
            SolverStep.NewCycle, zero(T), zero(T), Axis.X, 1, 0,
            splitting, riemann, limiter, projection, test_case,
            global_dt, steps_ranges, blk_logs, 0
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
        first(params.steps_ranges),
        params.estimated_blk_log_size
    )
end


function next_axis_sweep!(params::ArmonParameters, state::SolverState)
    if state.axis_splitting_idx == 0
        iter_val = iterate(split_axes(state, ndims(params)))
    else
        iter_val = iterate(split_axes(state, ndims(params)), state.axis_splitting_idx)
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


function update_solver_state!(params::ArmonParameters, state::SolverState, axis::Axis.T, dt_factor)
    state.dx = params.domain_size[axis] / params.global_grid[axis]
    state.dt = state.global_dt.current_dt * dt_factor
    state.axis = axis
    state.steps_ranges = params.steps_ranges[axis]
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
    state.axis = Axis.X
    state.axis_splitting_idx = 1
    state.cycle = 0
    empty!(state.blk_logs)
    state.total_stalls = 0
end


"""
    BLOCK_LOG_THREAD_LOCAL_STORAGE::Dict{UInt16, Int32}

Incremented by 1 every time a `BlockLogEvent` is created in a thread, i.e. each time a block has
solver kernels applied to it through [`block_state_machine`](@ref).

Since only differences between values are interesting, no need to reset it.
"""
const BLOCK_LOG_THREAD_LOCAL_STORAGE = Dict{UInt16, Int32}()


function BlockLogEvent(blk_state::SolverState, new_state::SolverStep.T, steps_count, steps_vars, steps_var_count)
    tid = convert(UInt16, Threads.threadid())
    steps_count = convert(UInt8, steps_count)
    tid_block_event_counter = BLOCK_LOG_THREAD_LOCAL_STORAGE[tid] += 1
    stalls = blk_state.total_stalls
    blk_state.total_stalls = 0
    return BlockLogEvent(
        blk_state.cycle, tid, blk_state.axis, new_state,
        steps_count, steps_vars, steps_var_count,
        tid_block_event_counter, stalls
    )
end


push_log!(state::SolverState, blk_log::BlockLogEvent) = push!(state.blk_logs, blk_log)
