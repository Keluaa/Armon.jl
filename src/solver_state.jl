
@enumx TimeStepState::UInt32 begin
    "Blocks can to contribute to the next cycle's time step"
    LocalReady
    "One thread is starting the global reduction for the next cycle's time step"
    GlobalStart
    "The global reduction is in progress"
    GlobalInProgress
    "The global reduction is complete"
    GlobalDone
    "The next cycle's time step is available"
    AllDone
end


"""
    GlobalTimeStep

Holds all information about the current time and time step for the current solver cycle. This struct
is global and shared among all blocks.

When reaching `next_time_step`, blocks will contribute to the calculation of the time step for the
next cycle. Once it is done, the global reduction among all MPI processes will start. Starting and
completing the MPI reduction is done only by the main thread if `params.thread_split_comm == true`.
"""
mutable struct GlobalTimeStep{T}
    state          :: Atomic{TimeStepState.T}
    state_lock     :: Atomic{Int}
    cycle          :: Int
    time           :: T
    current_dt     :: T
    next_cycle_dt  :: T  # Result of the reduction for `next_dt`. `Inf` if not ready.
    next_dt        :: Atomic{T}  # Time step accumulator
    contributions  :: Atomic{Int}
    expected_count :: Int
    reduction_data :: AbstractCommunication{Vector{T}}

    function GlobalTimeStep{T}(params::ArmonParameters) where {T}
        # TODO: wrap the model with a Communications.ThreadCollective, then remove most atomic/thread related logic
        model = params.use_MPI ? params.reduc_model : Communications.NoCommunicationModel()
        reduc_data = Communications.init_reduce_broadcast(model, MPI.MIN, Vector{T}, 1)
        return new{T}(
            Atomic(TimeStepState.LocalReady), Atomic(0),
            0, zero(T), params.cst_dt ? params.Dt : zero(T),
            typemax(T), Atomic(typemax(T)),
            Atomic(0), 0,
            reduc_data
        )
    end
end


function reset!(global_dt::GlobalTimeStep{T}, params::ArmonParameters{T}, block_count) where {T}
    @atomic global_dt.state.x = TimeStepState.LocalReady
    @atomic global_dt.state_lock.x = 0
    global_dt.cycle = 0
    global_dt.time = zero(T)
    global_dt.current_dt = params.cst_dt ? params.Dt : zero(T)
    global_dt.next_cycle_dt = typemax(T)
    @atomic global_dt.next_dt.x = typemax(T)
    @atomic global_dt.contributions.x = 0
    global_dt.expected_count = block_count
end


function can_touch_global_mpi_time_step(params::ArmonParameters)
    # `thread_split_comm` imposes that all communications started by a thread are tested and completed
    # by the same thread. For simplicity, the main thread is in charge of handling the global MPI
    # reduction for the time step.
    return !params.use_MPI || !params.thread_split_comm || Threads.threadid() == 1
end


function advance_time_step_state!(
    params::ArmonParameters{T}, global_dt::GlobalTimeStep{T};
    wait_until_global_done=false, force=false
) where {T}
    # We use a state machine to control the time step state transitions and actions, as it has a bit
    # of complex logic to be thread-safe, MPI-safe, and support additional constraints enforced by
    # `params.thread_split_comm`.

    # The global lock ensures that only a single thread updates the state. Most MPI operations on the
    # same request are thread-unsafe, so the lock has two purposes.
    if force
        Communications.wait_acquire_atomic_lock!(global_dt.state_lock)
    else
        # Always use non-blocking locks by default
        locked = Communications.try_acquire_atomic_lock!(global_dt.state_lock)
        !locked && return (@atomic global_dt.state.x)
    end

    state = @atomic global_dt.state.x

    @label next_state
    new_state = state

    if state == TimeStepState.LocalReady
        contributions = @atomic global_dt.contributions.x
        if contributions == global_dt.expected_count
            # All blocks have contributed: the local (block-wise) reduction is done
            if params.use_MPI
                new_state = TimeStepState.GlobalStart
            else
                global_dt.next_cycle_dt = @atomic global_dt.next_dt.x
                new_state = TimeStepState.GlobalDone
            end
        end

    elseif state == TimeStepState.GlobalStart
        if can_touch_global_mpi_time_step(params)
            send_buf = Communications.acquire_send_buffer!(global_dt.reduction_data)
            send_buf[1] = @atomic global_dt.next_dt.x
            Communications.release_send_buffer!(global_dt.reduction_data)
            new_state = TimeStepState.GlobalInProgress
        end

    elseif state == TimeStepState.GlobalInProgress
        if can_touch_global_mpi_time_step(params)
            global_done = Communications.recv_completed(global_dt.reduction_data)
            if wait_until_global_done
                Communications.wait_recv_completed(global_dt.reduction_data)
                global_done = true
            end

            if global_done
                recv_buf = Communications.acquire_recv_buffer!(global_dt.reduction_data)
                global_dt.next_cycle_dt = recv_buf[1]
                Communications.release_recv_buffer!(global_dt.reduction_data)
                new_state = TimeStepState.GlobalDone
            end
        end

    elseif state == TimeStepState.GlobalDone
        # Apply the CFL condition
        prev_Δt = global_dt.current_dt
        next_Δt = global_dt.next_cycle_dt

        if (!isfinite(next_Δt) || next_Δt ≤ 0)
            solver_error(:time, "Invalid next time step for cycle $(global_dt.cycle): $next_Δt")
        elseif prev_Δt == 0
            # First cycle time step initialization
            next_Δt = params.cfl * next_Δt
        else
            # CFL condition and maximum increase per cycle of the time step
            next_Δt = convert(T, min(params.cfl * next_Δt, 1.05 * prev_Δt))
        end

        global_dt.next_cycle_dt = next_Δt

        if global_dt.current_dt == 0
            # The current time step needs to be initialized (first cycle)
            global_dt.current_dt = global_dt.next_cycle_dt
        end

        # Reset the local contributions
        @atomic global_dt.next_dt.x = typemax(T)
        @atomic global_dt.contributions.x = 0

        new_state = TimeStepState.AllDone

    elseif state == TimeStepState.AllDone
        # Nothing more to do. It is now up to `next_cycle!` to reset the state.

    else
        error("unknown state: $state")
    end

    if new_state != state
        @atomic global_dt.state.x = new_state
        state = new_state
        @goto next_state
    end

    Communications.release_atomic_lock!(global_dt.state_lock)
    return state
end


function loop_until_time_step_available(params::ArmonParameters, global_dt::GlobalTimeStep)
    !can_touch_global_mpi_time_step(params) && return 0

    # This exist only to avoid a very specific and rare deadlock, where the main thread doesn't
    # contribute to the local time step (as it has no assigned blocks), and we are in the first
    # cycle, where we must wait for the time step to be available before finishing the cycle.
    # Not my proudest lines of code...

    Δt_state = advance_time_step_state!(params, global_dt; force=true, wait_until_global_done=true)
    loop_count = 1
    while Δt_state ≠ TimeStepState.AllDone
        GC.safepoint()
        µs_to_wait = 2^clamp(loop_count, 1, 13)
        # Avoid to use Julia's `sleep` as this is supposed to be called in a multithreaded loop
        Libc.systemsleep(µs_to_wait * 1e-6)
        loop_count += 1
        Δt_state = advance_time_step_state!(params, global_dt; force=true, wait_until_global_done=true)
    end

    return loop_count
end


function contribute_to_local_time_step!(params::ArmonParameters, global_dt::GlobalTimeStep{T}, dt::T; all_blocks=false) where {T}
    @atomic global_dt.next_dt.x min dt  # Atomic reduction

    # Either the whole grid or a single block contributed
    contributed_blocks = all_blocks ? global_dt.expected_count : 1
    contributions = @atomic global_dt.contributions.x += contributed_blocks

    if contributions == global_dt.expected_count
        # Try to advance the state only when all local blocks have contributed
        return advance_time_step_state!(params, global_dt)
    else
        return TimeStepState.LocalReady
    end
end


function can_contribute_to_local_time_step(params::ArmonParameters, global_dt::GlobalTimeStep, current_cycle)
    global_dt.cycle != current_cycle && return false

    dt_state = @atomic global_dt.state.x
    if dt_state != TimeStepState.LocalReady
        # The time step state is maybe not up-to-date with the MPI state, etc...
        dt_state = wait_for_time_step!(params, global_dt)
    end

    return dt_state == TimeStepState.LocalReady
end


wait_for_time_step!(params, global_dt) = advance_time_step_state!(params, global_dt; wait_until_global_done=true)


function next_cycle!(params::ArmonParameters, global_dt::GlobalTimeStep{T}) where {T}
    if params.cst_dt
        global_dt.current_dt = global_dt.next_cycle_dt = params.Dt
    else
        time_step_state = wait_for_time_step!(params, global_dt)
        if time_step_state != TimeStepState.AllDone
            error("expected time step to be done, got: $time_step_state")
        end
    end

    global_dt.time += global_dt.current_dt
    global_dt.cycle += 1

    # Reset the time step reduction state
    @atomic global_dt.state.x = TimeStepState.LocalReady
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
    RemapAdvection
    RemapProjection
    EndCycle
    ErrorState
end


# Bit flags for each of the variables of a block
const STEPS_VARS_FLAGS = (;
    x      = 0b0000_0000_0000_0001,
    y      = 0b0000_0000_0000_0010,
    ρ      = 0b0000_0000_0000_0100,
    u      = 0b0000_0000_0000_1000,
    v      = 0b0000_0000_0001_0000,
    E      = 0b0000_0000_0010_0000,
    p      = 0b0000_0000_0100_0000,
    c      = 0b0000_0000_1000_0000,
    g      = 0b0000_0001_0000_0000,
    uˢ     = 0b0000_0010_0000_0000,
    pˢ     = 0b0000_0100_0000_0000,
    work_1 = 0b0000_1000_0000_0000,
    work_2 = 0b0001_0000_0000_0000,
    work_3 = 0b0010_0000_0000_0000,
    work_4 = 0b0100_0000_0000_0000,
    mask   = 0b1000_0000_0000_0000,
)

# Flags for the arrays used by kernels: `count_ones(SOLVER_STEPS_VARS[step][2])` represents the
# number of arrays the kernels of `step` can bring into the cache. If `SOLVER_STEPS_VARS[step][1] == true`
# then the kernel uses one of `STEPS_VARS_FLAGS.u` or `STEPS_VARS_FLAGS.v` depending on the current
# axis.
# TODO: deduce them from kernel + steps definitions?
const SOLVER_STEPS_VARS = Dict{SolverStep.T, Tuple{Bool, UInt16}}(
    SolverStep.NewCycle        => (false, 0),
    SolverStep.TimeStep        => (false, STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.c),
    SolverStep.InitTimeStep    => (false, STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.c),
    SolverStep.NewSweep        => (false, 0),
    SolverStep.EOS             => (false, STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.g),
    SolverStep.Exchange        => (false, STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.g),
    SolverStep.Fluxes          => (true,  STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.pˢ),
    SolverStep.CellUpdate      => (true,  STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.pˢ),
    SolverStep.RemapAdvection  => (false, STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.work_1 | STEPS_VARS_FLAGS.work_2 | STEPS_VARS_FLAGS.work_3 | STEPS_VARS_FLAGS.work_4),
    SolverStep.RemapProjection => (false, STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.work_1 | STEPS_VARS_FLAGS.work_2 | STEPS_VARS_FLAGS.work_3 | STEPS_VARS_FLAGS.work_4),
    SolverStep.EndCycle        => (false, 0),
    SolverStep.ErrorState      => (false, 0),
)


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
    SolverSchemes

The different numerical schemes to use in the solver, and their parameters.
"""
struct SolverSchemes{Splitting, Riemann, RiemannLimiter, Projection, TestCase}
    splitting         :: Splitting
    riemann_scheme    :: Riemann
    riemann_limiter   :: RiemannLimiter
    projection_scheme :: Projection
    test_case         :: TestCase

    function SolverSchemes(
        splitting::S, riemann::R, limiter::RL, projection::P, test_case::TC
    ) where {
        S <: SplittingMethod, R <: RiemannScheme, RL <: Limiter, P <: ProjectionScheme, TC <: TestCase
    }
        return new{S, R, RL, P, TC}(splitting, riemann, limiter, projection, test_case)
    end
end


function SolverSchemes(params::ArmonParameters)
    return SolverSchemes(
        params.axis_splitting,
        params.riemann_scheme, params.riemann_limiter,
        params.projection_scheme,
        params.test
    )
end


"""
    SolverState

Object containing all non-constant parameters needed to run the solver, as well as type-parameters
needed to avoid runtime dispatch.

This object is local to a block (or set of blocks): multiple blocks could be at different steps of
the solver at once.
"""
mutable struct SolverState{T, Schemes <: SolverSchemes, StepsRangesArray <: AbstractArray{StepsRanges}, Queue <: AbstractStepQueue}
    step               :: SolverStep.T  # Solver step the associated block is at. Unused if `params.async_cycle == false`
    dx                 :: T    # Space step along the current axis
    dt                 :: T    # Scaled time step for the current cycle
    axis               :: Axis.T
    axis_splitting_idx :: Int
    cycle              :: Int  # Local cycle of the block
    schemes            :: Schemes
    global_dt          :: GlobalTimeStep{T}
    steps_ranges       :: Vector{StepsRanges}
    device_ranges      :: StepsRangesArray
    queue              :: Queue  # Queue to schedule the solver steps to the device
    blk_logs           :: Vector{BlockLogEvent}
    total_stalls       :: Int

    function SolverState{T}(schemes::Schemes, global_dt, steps_ranges, device_steps_ranges, queue::Queue, log_size) where {T, Schemes, Queue}
        blk_logs = Vector{BlockLogEvent}()
        log_size > 0 && sizehint!(blk_logs, log_size)
        return new{T, Schemes, typeof(device_steps_ranges), Queue}(
            SolverStep.NewCycle, zero(T), zero(T), Axis.X, 1, 0,
            schemes, global_dt, steps_ranges, device_steps_ranges, queue,
            blk_logs, 0
        )
    end
end


function SolverState(params::ArmonParameters{T}, global_dt::GlobalTimeStep{T}) where {T}
    schemes = SolverSchemes(params)
    if params.use_step_queue
        queue = StepQueue(params.device, params.step_queue_capacity)
    else
        queue = NoQueue(params.device)
    end
    return SolverState{T}(
        schemes, global_dt, params.steps_ranges, params.device_steps_ranges, queue,
        params.estimated_blk_log_size
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


function update_solver_state!(params::ArmonParameters, state::SolverState, axis::Axis.T, dt_factor)
    i_ax = Int(axis)
    state.dx = params.domain_size[i_ax] / params.global_grid[i_ax]
    state.dt = state.global_dt.current_dt * dt_factor
    state.axis = axis
end


function need_to_update_time_step(params::ArmonParameters, state::SolverState)
    params.cst_dt && return false
    state.dt == 0 && return true  # first cycle
    if params.dt_on_even_cycles
        return iseven(state.global_dt.cycle)
    else
        return true
    end
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
    BasicSolverState

Immutable version of a lightweight [`SolverState`](@ref), for use in GPU kernels.
"""
struct BasicSolverState{T, Schemes <: SolverSchemes, StepsRangesArray <: AbstractArray{StepsRanges}}
    dx                :: T
    dt                :: T
    axis              :: Axis.T
    schemes           :: Schemes
    steps_ranges      :: StepsRangesArray
end

BasicSolverState(solver_state::SolverState) =
    BasicSolverState(solver_state.dx, solver_state.dt, solver_state.axis, solver_state.schemes, solver_state.device_ranges)

Base.eltype(::ObjOrType{BasicSolverState{T}}) where {T} = T

Adapt.adapt_structure(to, bss::BasicSolverState) =
    BasicSolverState(bss.dx, bss.dt, bss.axis, bss.schemes, Adapt.adapt(to, bss.steps_ranges))


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
