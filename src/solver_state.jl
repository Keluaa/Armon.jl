
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


struct BlockGridLog
    blk_logs           :: Array{Vector{BlockLogEvent}}
    blk_sizes          :: Array{NTuple{2, Int}}
    ghosts             :: Int
    mean_blk_cells     :: Float64  # Mean number of cells in all blocks
    mean_vars_per_cell :: Float64  # Mean number of variables in all cells
    var_data_type_size :: Int      # Byte size of variables' data type
end


mutable struct LogStat{T}
    min  :: T
    max  :: T
    tot  :: T
    mean :: Float64

    LogStat{T}() where {T} = new{T}(typemax(T), typemin(T), zero(T), zero(Float64))
end


function Base.:*(ls::LogStat{T}, x::V) where {T, V}
    r_ls = LogStat{typeof(ls.min * x)}()
    r_ls.min = ls.min * x; r_ls.max = ls.max * x; r_ls.tot = ls.tot * x; r_ls.mean = ls.mean * x
    return r_ls
end

function Base.:*(ls₁::LogStat{T}, ls₂::LogStat{V}) where {T, V}
    r_ls = LogStat{typeof(ls₁.min * ls₂.min)}()
    r_ls.min = ls₁.min * ls₂.min; r_ls.max = ls₁.max * ls₂.max
    r_ls.tot = ls₁.tot * ls₂.tot; r_ls.mean = ls₁.mean * ls₂.mean
    return r_ls
end


mutable struct BlockGridLogThreadStats
    tot_blk        :: Int
    tot_events     :: Int
    tot_bytes      :: Int  # Size (in bytes) of all real cell main variables in all of the thread's blocks
    tot_used_bytes :: Int  # Abount of bytes which were used at least once by the thread's blocks
    tot_stalls     :: Int  # Number of calls to `block_state_machine` which didn't progress a block's state

    BlockGridLogThreadStats() = new(0, 0, 0, 0, 0)
end


mutable struct BlockGridLogStats
    tot_blk                 :: Int

    # Number of blocks with an inconsistent processing thread (>0 invalidates most measurements)
    inconsistent_threads    :: Int

    events_per_blk          :: LogStat{Int}  # Number of events per block
    steps_per_event         :: LogStat{Int}  # Number of solver steps per event
    blk_before_per_event    :: LogStat{Int}  # Number of blocks processed by the thread between events of the same block
    indep_vars_per_event    :: LogStat{Int}  # Number of unique variables used during the event
    indep_bytes_per_event   :: LogStat{Int}  # Bytes of unique variables used during the event
    vars_per_event          :: LogStat{Int}  # Total number of variables used during the event
    stalls_per_event        :: LogStat{Int}  # Number of calls to `block_state_machine` which didn't progress the block's state

    event_size              :: LogStat{Float64}  # Total bytes used during a block event
    event_indep_size        :: LogStat{Float64}  # Unique bytes used during a block event
    bytes_prev_blk          :: LogStat{Float64}  # Total bytes processed by the thread between the same block
    indep_bytes_prev_blk    :: LogStat{Float64}  # Unique bytes processed by the thread between the same block

    # Number of times an event stopped at each solver step
    steps_stats             :: Dict{SolverStep.T, Int}

    # Thread stats
    threads_stats           :: Dict{UInt16, BlockGridLogThreadStats}
    active_threads          :: Int
    blk_per_thread          :: LogStat{Int}
    events_per_thread       :: LogStat{Int}
    stalls_per_thread       :: LogStat{Int}      # Stalls for all of the thread's blocks
    bytes_prev_thread       :: LogStat{Float64}  # Total bytes processed by the thread between each solver iteration
    indep_bytes_prev_thread :: LogStat{Float64}  # Unique bytes processed by the thread between each solver interation

    # Stats from the original `BlockGrid`
    blk_ghosts              :: Int
    mean_blk_cells          :: Float64
    mean_vars_per_cell      :: Float64
    var_data_type_size      :: Int
    mean_blk_size           :: Float64

    BlockGridLogStats() = new(
        0, 0,
        LogStat{Int}(), LogStat{Int}(), LogStat{Int}(), LogStat{Int}(),
        LogStat{Int}(), LogStat{Int}(), LogStat{Int}(),
        LogStat{Float64}(), LogStat{Float64}(),
        LogStat{Float64}(), LogStat{Float64}(),
        Dict{SolverStep.T, Int}(step => 0 for step in instances(SolverStep.T)),
        Dict{UInt16, BlockGridLogThreadStats}(), 0,
        LogStat{Int}(), LogStat{Int}(), LogStat{Int}(),
        LogStat{Float64}(), LogStat{Float64}(),
        0, zero(Float64), zero(Float64), 0, zero(Float64)
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

# Flags for the arrays used by kernels: `count_ones(SOLVER_STEPS_VARS[step][2])` represents the
# number of arrays the kernels of `step` can bring into the cache. If `SOLVER_STEPS_VARS[step][1] == true`
# then the kernel uses one of `STEPS_VARS_FLAGS.u` or `STEPS_VARS_FLAGS.v` depending on the current
# axis.
# TODO: deduce them from kernel + steps definitions?
const SOLVER_STEPS_VARS = Dict{SolverStep.T, Tuple{Bool, UInt16}}(
    SolverStep.NewCycle     => (false, 0),
    SolverStep.TimeStep     => (false, STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.c),
    SolverStep.InitTimeStep => (false, STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.c),
    SolverStep.NewSweep     => (false, 0),
    SolverStep.EOS          => (false, STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.g),
    SolverStep.Exchange     => (false, STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.g),
    SolverStep.Fluxes       => (true,  STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.p | STEPS_VARS_FLAGS.c | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.pˢ),
    SolverStep.CellUpdate   => (true,  STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.pˢ),
    SolverStep.Remap        => (false, STEPS_VARS_FLAGS.ρ | STEPS_VARS_FLAGS.E | STEPS_VARS_FLAGS.u | STEPS_VARS_FLAGS.v | STEPS_VARS_FLAGS.uˢ| STEPS_VARS_FLAGS.work_1 | STEPS_VARS_FLAGS.work_2 | STEPS_VARS_FLAGS.work_3 | STEPS_VARS_FLAGS.work_4),
    SolverStep.EndCycle     => (false, 0),
    SolverStep.ErrorState   => (false, 0),
)


function accumulate_grid_stats!(stat::LogStat, value)
    stat.min = min(stat.min, value)
    stat.max = max(stat.max, value)
    stat.tot += value
end


function accumulate_grid_stats!(gs::BlockGridLogStats, ::CartesianIndex, blk_events::Vector{BlockLogEvent}, blk_size)
    gs.tot_blk += 1

    event_count = length(blk_events)
    accumulate_grid_stats!(gs.events_per_blk, event_count)

    block_cells = prod(blk_size .- gs.blk_ghosts)  # We only consider real cells
    tot_blk_bytes = 0
    used_blk_vars = UInt16(0)
    tot_stalls = 0

    expected_tid = first(blk_events).tid
    prev_tid_blk_idx = first(blk_events).tid_blk_idx
    for event in blk_events
        accumulate_grid_stats!(gs.steps_per_event, event.steps_count)

        used_blk_vars |= event.steps_vars
        indep_var_count = count_ones(event.steps_vars)
        accumulate_grid_stats!(gs.indep_vars_per_event, indep_var_count)

        event_bytes = indep_var_count * block_cells * gs.var_data_type_size
        tot_blk_bytes += event_bytes
        accumulate_grid_stats!(gs.indep_bytes_per_event, event_bytes)

        accumulate_grid_stats!(gs.vars_per_event, event.steps_var_count)

        tot_stalls += event.stalls
        accumulate_grid_stats!(gs.stalls_per_event, event.stalls)

        if event.tid == expected_tid
            blk_before = event.tid_blk_idx - prev_tid_blk_idx
            prev_tid_blk_idx = event.tid_blk_idx
            if blk_before != 0  # 0 if it is the first event of the block, but we only care about differences
                accumulate_grid_stats!(gs.blk_before_per_event, blk_before)
            end
        else
            gs.inconsistent_threads += 1
        end

        gs.steps_stats[event.new_state] += 1
    end

    tot_used_blk_bytes = count_ones(used_blk_vars) * block_cells * gs.var_data_type_size

    ts = get!(BlockGridLogThreadStats, gs.threads_stats, expected_tid)
    ts.tot_blk += 1
    ts.tot_events += event_count
    ts.tot_bytes += tot_blk_bytes
    ts.tot_used_bytes += tot_used_blk_bytes
    ts.tot_stalls += tot_stalls

    return gs
end


function analyse_log_stats(f, grid_log::BlockGridLog)
    for pos in eachindex(Base.IndexCartesian(), grid_log.blk_logs)
        !isassigned(grid_log.blk_logs, pos) && continue
        f(pos, grid_log.blk_logs[pos], grid_log.blk_sizes[pos])
    end
end


function analyse_log_stats(grid_log::BlockGridLog)
    gs = BlockGridLogStats()
    gs.blk_ghosts = grid_log.ghosts
    gs.mean_blk_cells = grid_log.mean_blk_cells
    gs.mean_vars_per_cell = grid_log.mean_vars_per_cell
    gs.var_data_type_size = grid_log.var_data_type_size
    gs.mean_blk_size = gs.mean_blk_cells * gs.mean_vars_per_cell * gs.var_data_type_size

    analyse_log_stats((args...) -> accumulate_grid_stats!(gs, args...), grid_log)

    tot_events = gs.events_per_blk.tot
    gs.events_per_blk.mean = tot_events / gs.tot_blk
    gs.steps_per_event.mean = gs.steps_per_event.tot / tot_events
    gs.blk_before_per_event.mean = gs.blk_before_per_event.tot / (tot_events - gs.tot_blk)  # `N - 1` since we only care about differences
    gs.indep_vars_per_event.mean = gs.indep_vars_per_event.tot / tot_events
    gs.indep_bytes_per_event.mean = gs.indep_bytes_per_event.tot / tot_events
    gs.vars_per_event.mean = gs.vars_per_event.tot / tot_events
    gs.stalls_per_event.mean = gs.stalls_per_event.tot / tot_events

    gs.active_threads = length(gs.threads_stats)
    for (_, ts) in gs.threads_stats
        accumulate_grid_stats!(gs.blk_per_thread, ts.tot_blk)
        accumulate_grid_stats!(gs.events_per_thread, ts.tot_events)
        accumulate_grid_stats!(gs.stalls_per_thread, ts.tot_stalls)
    end
    gs.blk_per_thread.mean    = gs.blk_per_thread.tot    / gs.active_threads
    gs.events_per_thread.mean = gs.events_per_thread.tot / gs.active_threads
    gs.stalls_per_thread.mean = gs.stalls_per_thread.tot / gs.active_threads

    gs.event_size = gs.vars_per_event * gs.mean_blk_cells * gs.var_data_type_size
    gs.event_indep_size = gs.indep_vars_per_event * gs.mean_blk_cells * gs.var_data_type_size

    gs.bytes_prev_blk = gs.blk_before_per_event * gs.event_size
    gs.indep_bytes_prev_blk = gs.blk_before_per_event * gs.event_indep_size

    gs.bytes_prev_thread = gs.bytes_prev_blk * gs.blk_per_thread
    gs.indep_bytes_prev_thread = gs.indep_bytes_prev_blk * gs.blk_per_thread

    return gs
end


function as_SI_magnitude(x)
    prefixes = ["", "k", "M", "G", "T", "P"]
    mag = floor(Int, log(1000, abs(x)))
    mag = clamp(mag, 0, length(prefixes) - 1)
    prefix = prefixes[mag + 1]
    return x / 1000^mag, prefix
end


function Base.show(io::IO, ls::LogStat)
    print(io, @sprintf("%6.2f", ls.mean), ", ", ls.min, "..", ls.max, "\t(mean, min..max)")
end


function Base.show(io::IO, ::MIME"text/plain", gs::BlockGridLogStats)
    println(io, "BlockGrid solve stats:")
    println(io, " - total blocks\t\t\t", gs.tot_blk)
    println(io, " - total events\t\t\t", gs.events_per_blk.tot)
    println(io, " - total steps \t\t\t", gs.steps_per_event.tot)
    printstyled(io, " - inconsistent threads\t\t", gs.inconsistent_threads, "\n";
        color=gs.inconsistent_threads == 0 ? :normal : :red)
    println(io, " - events per block\t\t", gs.events_per_blk)
    println(io, " - steps per event\t\t", gs.steps_per_event)
    println(io, " - vars per event\t\t", gs.vars_per_event)
    println(io, " - indep vars per event\t\t", gs.indep_vars_per_event)
    println(io, " - blocks before event\t\t", gs.blk_before_per_event)
    println(io, " - mean block cells\t\t", @sprintf("%5.2f", gs.mean_blk_cells))
    println(io, " - mean vars per cell\t\t", @sprintf("%5.2f", gs.mean_vars_per_cell), ", ",
        gs.var_data_type_size, " bytes per var")

    println(io, " - blocks per thread\t\t", gs.blk_per_thread, ", ", gs.active_threads, " active threads")
    println(io, " - events per thread\t\t", gs.events_per_thread)

    println(io, " - mean block size\t\t", @sprintf("%6.2f %sB", as_SI_magnitude(gs.mean_blk_size)...))
    println(io, " - mean event size\t\t", @sprintf("%6.2f %sB", as_SI_magnitude(gs.event_size.mean)...), ", ",
        @sprintf("%6.2f %sB", as_SI_magnitude(gs.event_indep_size.mean)...), "\t(total, unique)")
    println(io, " - mean bytes before block\t", @sprintf("%6.2f %sB", as_SI_magnitude(gs.bytes_prev_blk.mean)...), ", ",
        @sprintf("%6.2f %sB", as_SI_magnitude(gs.indep_bytes_prev_blk.mean)...), "\t(total, unique)")

    println(io, " - mean bytes per thread iter\t",
        @sprintf("%6.2f %sB", as_SI_magnitude(gs.bytes_prev_thread.mean)...), ", ",
        @sprintf("%6.2f %sB", as_SI_magnitude(gs.indep_bytes_prev_thread.mean)...), "\t(total, unique)")

    println(io, " - stalls per event\t\t", gs.stalls_per_event)
    println(io, " - stalls per thread\t\t", gs.stalls_per_thread)

    print(io, " - steps stop count")
    for (step, cnt) in filter(≠(0) ∘ last, gs.steps_stats |> collect |> sort!)
        println(io)
        print(io, "    - $step = $cnt")
    end
end
