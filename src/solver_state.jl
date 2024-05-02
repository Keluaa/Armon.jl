
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


struct BlockLogEvent
    cycle           :: Int32   # Cycle at which the event occured
    tid             :: UInt8   # Thread which processed the block
    axis            :: Axis.T  # Final axis of the block
    new_state       :: SolverStep.T  # Final state of the block
    steps_count     :: UInt8   # Number of solver steps done
    steps_vars      :: UInt16  # Flag with a 1 when a variable was used
    steps_var_count :: Int16   # Number of times all variables were used
    tid_blk_idx     :: Int32   # Number of blocks processed by the thread before this event
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
            global_dt, steps_ranges, blk_logs
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
    state.axis = Axis.X
    state.axis_splitting_idx = 1
    state.cycle = 0
    empty!(state.blk_logs)
end


"""
    BLOCK_LOG_THREAD_LOCAL_STORAGE::Dict{Int, Int}

Incremented by 1 every time a `BlockLogEvent` is created in a thread, i.e. each time a block has
solver kernels applied to it through [`block_state_machine`](@ref).

Since only differences between values are interesting, no need to reset it.
"""
const BLOCK_LOG_THREAD_LOCAL_STORAGE = Dict{Int, Int}(tid => 0 for tid in Threads.nthreads())


function BlockLogEvent(blk_state::SolverState, new_state::SolverStep.T, steps_count, steps_vars, steps_var_count)
    tid = convert(UInt8, Threads.threadid())
    steps_count = convert(UInt8, steps_count)
    tid_block_event_counter = BLOCK_LOG_THREAD_LOCAL_STORAGE[tid] += 1
    return BlockLogEvent(
        blk_state.cycle, tid, blk_state.axis, new_state,
        steps_count, steps_vars, steps_var_count,
        tid_block_event_counter
    )
end


push_log!(state::SolverState, blk_log::BlockLogEvent) = push!(state.blk_logs, blk_log)


struct BlockGridLog
    blk_logs           :: Array{Vector{BlockLogEvent}}
    mean_blk_cells     :: Float64  # Mean number of cells in all blocks
    mean_vars_per_cell :: Float64  # Mean number of variables in all cells
    var_data_type_size :: Int      # Byte size of variables' data type
end


mutable struct BlockGridLogStats
    tot_blk                       :: Int
    tot_events                    :: Int
    tot_steps                     :: Int

    # Number of blocks with an inconsistent processing thread (>0 invalidates most measurements)
    inconsistent_threads          :: Int

    # Number of events per block
    min_events_per_blk            :: Int
    max_events_per_blk            :: Int
    mean_events_per_blk           :: Float64

    # Number of solver steps per event
    min_steps_per_event           :: Int
    max_steps_per_event           :: Int
    mean_steps_per_event          :: Float64

    # Number of blocks processed by the thread between events of the same block
    min_blk_before_per_event      :: Int
    max_blk_before_per_event      :: Int
    mean_blk_before_per_event     :: Float64

    # Number of unique variables used during the event
    min_indep_vars_per_event      :: Int
    max_indep_vars_per_event      :: Int
    mean_indep_vars_per_event     :: Float64

    # Total number of variables used during the event
    min_vars_per_event            :: Int
    max_vars_per_event            :: Int
    mean_vars_per_event           :: Float64

    # Number of times an event stopped at each solver step
    steps_stats                   :: Dict{SolverStep.T, Int}

    mean_blk_cells                :: Float64
    mean_vars_per_cell            :: Float64
    var_data_type_size            :: Int

    BlockGridLogStats() = new(
        0, 0, 0, 0,
        typemax(Int), typemin(Int), zero(Float64),
        typemax(Int), typemin(Int), zero(Float64),
        typemax(Int), typemin(Int), zero(Float64),
        typemax(Int), typemin(Int), zero(Float64),
        typemax(Int), typemin(Int), zero(Float64),
        Dict{SolverStep.T, Int}(step => 0 for step in instances(SolverStep.T)),
        zero(Float64), zero(Float64), 0
    )
end


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

# Those are the arrays used by kernels.
# Therefore `count_ones(SOLVER_STEPS_VARS[step])` represents the number of arrays the kernels of `step`
# can bring into the cache.
# TODO: deduce them from kernel + steps definitions?
const SOLVER_STEPS_VARS = Dict{SolverStep.T, UInt16}(
    SolverStep.NewCycle     => 0,
    SolverStep.TimeStep     => STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.c,
    SolverStep.InitTimeStep => STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.c,
    SolverStep.NewSweep     => 0,
    SolverStep.EOS          => STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.g,
    SolverStep.Exchange     => STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.g,
    SolverStep.Fluxes       => STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.pˢ,
    SolverStep.CellUpdate   => STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.pˢ,
    SolverStep.Remap        => STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.work_1 | STEPS_VARS_FLAGS.work_2 | STEPS_VARS_FLAGS.work_3 | STEPS_VARS_FLAGS.work_4,
    SolverStep.EndCycle     => 0,
    SolverStep.ErrorState   => 0,
)


function accumulate_grid_stats!(grid_stats::BlockGridLogStats, ::CartesianIndex, blk_events::Vector{BlockLogEvent})
    grid_stats.tot_blk += 1

    event_count = length(blk_events)
    grid_stats.tot_events += event_count
    grid_stats.min_events_per_blk = min(grid_stats.min_events_per_blk, event_count)
    grid_stats.max_events_per_blk = max(grid_stats.max_events_per_blk, event_count)

    min_steps = typemax(Int)
    max_steps = typemin(Int)
    tot_steps = 0
    expected_tid = first(blk_events).tid
    inconsistent_tids = 0
    prev_tid_blk_idx = first(blk_events).tid_blk_idx
    min_blk_before = typemax(Int)
    max_blk_before = typemin(Int)
    tot_blk_before = 0
    min_indep_vars_per_event = typemax(Int)
    max_indep_vars_per_event = typemin(Int)
    tot_indep_vars_per_event = 0
    min_vars_per_event = typemax(Int)
    max_vars_per_event = typemin(Int)
    tot_vars_per_event = 0
    for event in blk_events
        tot_steps += event.steps_count
        min_steps = min(min_steps, event.steps_count)
        max_steps = max(max_steps, event.steps_count)

        indep_var_count = count_ones(event.steps_vars)
        tot_indep_vars_per_event += indep_var_count
        min_indep_vars_per_event = min(min_indep_vars_per_event, indep_var_count)
        max_indep_vars_per_event = max(max_indep_vars_per_event, indep_var_count)

        tot_vars_per_event += event.steps_var_count
        min_vars_per_event = min(min_vars_per_event, event.steps_var_count)
        max_vars_per_event = max(max_vars_per_event, event.steps_var_count)

        if event.tid == expected_tid
            blk_before = event.tid_blk_idx - prev_tid_blk_idx
            if blk_before != 0  # 0 if it is the first event of the block, but we only care about differences
                min_blk_before = min(min_blk_before, blk_before)
                max_blk_before = max(max_blk_before, blk_before)
                tot_blk_before += blk_before
                prev_tid_blk_idx = event.tid_blk_idx
            end
        else
            inconsistent_tids += 1
        end

        grid_stats.steps_stats[event.new_state] += 1
    end

    grid_stats.inconsistent_threads += inconsistent_tids
    grid_stats.min_steps_per_event = min(grid_stats.min_steps_per_event, min_steps)
    grid_stats.max_steps_per_event = max(grid_stats.max_steps_per_event, max_steps)
    grid_stats.tot_steps += tot_steps
    grid_stats.min_blk_before_per_event = min(grid_stats.min_blk_before_per_event, min_blk_before)
    grid_stats.max_blk_before_per_event = max(grid_stats.max_blk_before_per_event, max_blk_before)
    grid_stats.mean_blk_before_per_event += tot_blk_before
    grid_stats.min_indep_vars_per_event = min(grid_stats.min_indep_vars_per_event, min_indep_vars_per_event)
    grid_stats.max_indep_vars_per_event = max(grid_stats.max_indep_vars_per_event, max_indep_vars_per_event)
    grid_stats.mean_indep_vars_per_event += tot_indep_vars_per_event
    grid_stats.min_vars_per_event = min(grid_stats.min_vars_per_event, min_vars_per_event)
    grid_stats.max_vars_per_event = max(grid_stats.max_vars_per_event, max_vars_per_event)
    grid_stats.mean_vars_per_event += tot_vars_per_event

    return grid_stats
end


function analyse_log_stats(f, grid_log::BlockGridLog)
    for pos in eachindex(Base.IndexCartesian(), grid_log.blk_logs)
        !isassigned(grid_log.blk_logs, pos) && continue
        f(pos, grid_log.blk_logs[pos])
    end
end


function analyse_log_stats(grid_log::BlockGridLog)
    grid_stats = BlockGridLogStats()
    analyse_log_stats((args...) -> accumulate_grid_stats!(grid_stats, args...), grid_log)

    grid_stats.mean_events_per_blk = grid_stats.tot_events / grid_stats.tot_blk
    grid_stats.mean_steps_per_event = grid_stats.tot_steps / grid_stats.tot_events
    grid_stats.mean_blk_before_per_event /= grid_stats.tot_events - grid_stats.tot_blk  # `N - 1` since we only care about differences
    grid_stats.mean_indep_vars_per_event /= grid_stats.tot_events
    grid_stats.mean_vars_per_event /= grid_stats.tot_events

    grid_stats.mean_blk_cells = grid_log.mean_blk_cells
    grid_stats.mean_vars_per_cell = grid_log.mean_vars_per_cell
    grid_stats.var_data_type_size = grid_log.var_data_type_size

    return grid_stats
end


function as_SI_magnitude(x)
    prefixes = ["", "k", "M", "G", "T", "P"]
    mag = floor(Int, log(1000, abs(x)))
    mag = clamp(mag, 0, length(prefixes) - 1)
    prefix = prefixes[mag + 1]
    return x / 1000^mag, prefix
end


function Base.show(io::IO, ::MIME"text/plain", grid_stats::BlockGridLogStats)
    println(io, "BlockGrid solve stats:")
    println(io, " - total blocks\t\t\t", grid_stats.tot_blk)
    println(io, " - total events\t\t\t", grid_stats.tot_events)
    println(io, " - total steps \t\t\t", grid_stats.tot_steps)
    printstyled(io, " - inconsistent threads\t\t", grid_stats.inconsistent_threads, "\n";
        color=grid_stats.inconsistent_threads == 0 ? :normal : :red)
    println(io, " - events per block\t\t", @sprintf("%5.2f", grid_stats.mean_events_per_blk), ", ",
        grid_stats.min_events_per_blk, "..", grid_stats.max_events_per_blk, "\t(mean, min..max)")
    println(io, " - steps per event\t\t", @sprintf("%5.2f", grid_stats.mean_steps_per_event), ", ",
        grid_stats.min_steps_per_event, "..", grid_stats.max_steps_per_event, "\t(mean, min..max)")
    println(io, " - vars per event\t\t", @sprintf("%5.2f", grid_stats.mean_vars_per_event), ", ",
        grid_stats.min_vars_per_event, "..", grid_stats.max_vars_per_event, "\t(mean, min..max)")
    println(io, " - indep vars per event\t\t", @sprintf("%5.2f", grid_stats.mean_indep_vars_per_event), ", ",
        grid_stats.min_indep_vars_per_event, "..", grid_stats.max_indep_vars_per_event, "\t(mean, min..max)")
    println(io, " - blocks before event\t\t", @sprintf("%5.2f", grid_stats.mean_blk_before_per_event), ", ",
        grid_stats.min_blk_before_per_event, "..", grid_stats.max_blk_before_per_event, "\t(mean, min..max)")
    println(io, " - mean block cells\t\t", @sprintf("%5.2f", grid_stats.mean_blk_cells))
    println(io, " - mean vars per cell\t\t", @sprintf("%5.2f", grid_stats.mean_vars_per_cell), ", ",
        grid_stats.var_data_type_size, " bytes per var")

    mean_blk_size = grid_stats.mean_blk_cells * grid_stats.mean_vars_per_cell * grid_stats.var_data_type_size
    mean_evt_size = grid_stats.mean_blk_cells * grid_stats.mean_vars_per_event * grid_stats.var_data_type_size
    mean_evt_size_indep = grid_stats.mean_blk_cells * grid_stats.mean_indep_vars_per_event * grid_stats.var_data_type_size
    mean_bytes_before_event = grid_stats.mean_blk_before_per_event * mean_evt_size
    mean_bytes_before_event_indep = grid_stats.mean_blk_before_per_event * mean_evt_size_indep

    println(io, " - mean block size\t\t", @sprintf("%6.2f %sB", as_SI_magnitude(mean_blk_size)...))
    println(io, " - mean event size\t\t", @sprintf("%6.2f %sB", as_SI_magnitude(mean_evt_size)...), ", ",
        @sprintf("%6.2f %sB", as_SI_magnitude(mean_evt_size_indep)...), "\t(total, unique)")
    println(io, " - mean bytes before event\t", @sprintf("%6.2f %sB", as_SI_magnitude(mean_bytes_before_event)...), ", ",
        @sprintf("%6.2f %sB", as_SI_magnitude(mean_bytes_before_event_indep)...), "\t(total, unique)")

    print(io, " - steps stop count")
    for (step, cnt) in filter(≠(0) ∘ last, grid_stats.steps_stats |> collect |> sort!)
        println(io)
        print(io, "    - $step = $cnt")
    end
end
