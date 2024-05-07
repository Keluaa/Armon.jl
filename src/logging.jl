
struct BlockGridLog
    blk_logs           :: Array{Vector{BlockLogEvent}}
    blk_sizes          :: Array{NTuple{2, Int}}
    ghosts             :: Int
    mean_blk_cells     :: Float64  # Mean number of cells in all blocks
    mean_vars_per_cell :: Float64  # Mean number of variables in all cells
    var_data_type_size :: Int      # Byte size of variables' data type
end


"""
    collect_logs(grid::BlockGrid)

Collect all [`BlockLogEvent`](@ref) of the `grid` into a single object (a `BlockGridLog`).
See [`analyse_log_stats`](@ref) to print metrics about all blocks.
"""
function collect_logs(grid::BlockGrid{T}) where {T}
    logs = Array{Vector{BlockLogEvent}}(undef, grid.grid_size)
    blk_sizes = Array{NTuple{2, Int}}(undef, grid.grid_size)
    tot_blk_cells = 0
    for blk in all_blocks(grid)
        logs[blk.pos] = solver_state(blk).blk_logs
        blk_sizes[blk.pos] = real_block_size(blk)
        tot_blk_cells += prod(real_block_size(blk))
    end
    mean_blk_cells = tot_blk_cells / prod(grid.grid_size)
    mean_vars_per_cell = length(main_vars())
    data_type_size = sizeof(T)
    return BlockGridLog(logs, blk_sizes, ghosts(grid), mean_blk_cells, mean_vars_per_cell, data_type_size)
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


"""
    BlockGridLogStats

Metrics about the blocks of a full solver excution.
"""
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


"""
    analyse_log_stats(f, grid_log::BlockGridLog)

Call `f` for each block of the `grid_log`, passing the block's position (a `CartesianIndex`), a
`Vector` of [`BlockLogEvent`](@ref)s, and the size of the block.
"""
function analyse_log_stats(f, grid_log::BlockGridLog)
    for pos in eachindex(Base.IndexCartesian(), grid_log.blk_logs)
        !isassigned(grid_log.blk_logs, pos) && continue
        f(pos, grid_log.blk_logs[pos], grid_log.blk_sizes[pos])
    end
end


"""
    analyse_log_stats(grid_log::BlockGridLog)

Crunch all data of `grid_log` into tangible metrics contained in a [`BlockGridLogStats`](@ref).
"""
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
