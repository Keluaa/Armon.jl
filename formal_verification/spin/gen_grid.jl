
using TOML


const USAGE = """
Script to generate a static specification of a N-D block grid for the SPIN specification.

Usage:
    julia ./gen_grid.jl <grid_file.toml>

See './grids/default_grid_3x3_4threads.toml' for an explaination of each option.
"""

const DEFAULT_ARGS = Dict(
    "grid"            => [3, 3],
    "periodic"        => [false, false],
    "num_threads"     => 4,
    "use_mpi"         => false,
    "processes"       => [1, 1],
    "real_block_xchg" => false,
    "num_cycles"      => 3,
    "max_sweeps"      => 2,
    "do_time_step"    => true,
    "do_halo_xchg"    => true,
)

const REMOTE_BLOCK   = 254
const NULL_BLOCK     = 255
const NULL_INTERFACE = 255
const NULL_RANK      = -1


function parse_options_file(filename)
    options = try
        open(TOML.tryparse, filename, "r")
    catch e
        println("Could not parse '$filename'")
        rethrow(e)
    end

    grid_dim = length(options["grid"])
    if grid_dim != 2
        DEFAULT_ARGS["periodic"]  = map(Returns(false), 1:grid_dim)
        DEFAULT_ARGS["processes"] = map(Returns(false), 1:grid_dim)
    end

    unknown_options = setdiff(keys(options), keys(DEFAULT_ARGS))
    !isempty(unknown_options) && error("unknown options: $(join(unknown_options, ", "))")

    options = merge(DEFAULT_ARGS, options)

    if !(length(options["grid"]) == length(options["periodic"]) == length(options["processes"]))
        error("incoherent dimension for 'grid' ($(length(options["grid"]))), \
               'periodic' ($(length(options["periodic"]))) \
               and 'processes' ($(length(options["processes"])))")
    end

    options["grid"] = Tuple(options["grid"])
    options["processes"] = Tuple(options["processes"])

    return options
end


function parse_arguments()
    length(ARGS) != 1 && error("Expected 1 argument. Usage:\n" * USAGE)
    ARGS[1] in ("-h", "--help") && (println(USAGE); exit())
    return parse_options_file(ARGS[1])
end


mutable struct P2PChannel
    idx   :: Int
    usage :: Int
end

P2PChannel(idx) = P2PChannel(idx, 0)


mutable struct RemoteBlock{Dim}
    idx       :: Int
    rank      :: Int
    pos       :: NTuple{Dim}
    tag       :: Int
    send_chan :: P2PChannel
    recv_chan :: P2PChannel
end

RemoteBlock(idx, rank, pos, tag) = RemoteBlock{length(pos)}(idx, rank, pos, tag, P2PChannel(-1), P2PChannel(-1))


mutable struct Interface
    idx :: Int
end


mutable struct Block{Dim, Neigh}
    idx        :: Int
    pos        :: NTuple{Dim, Int}
    neighbours :: NTuple{Neigh, Union{Block{Dim, Neigh}, RemoteBlock{Dim}}}
    interfaces :: NTuple{Neigh, Interface}
    tid        :: Int

    Block{D}(idx, pos) where {D} = new{D, 2D}(idx, pos)
    Block(idx, pos) = Block{length(pos)}(idx, pos)
end


mutable struct GridTopo{Dim, Neigh}
    pos               :: NTuple{Dim, Int}
    rank              :: Int
    neighbours        :: NTuple{Neigh, Int}
    threads           :: Int
    size              :: NTuple{Dim, Int}
    periodicity       :: NTuple{Dim, Bool}
    blocks            :: Vector{Block{Dim, Neigh}}
    interfaces        :: Vector{Interface}
    remote_blocks     :: Vector{RemoteBlock{Dim}}
    workloads_offset  :: Int
    threads_workloads :: Vector{Vector{Int}}

    function GridTopo{Dim, Neigh}(pos, rank, threads, grid_size, periodicity) where {Dim, Neigh}
        neighbours = ntuple(Returns(NULL_RANK), Neigh)
        threads_workloads = [Int[] for _ in 1:threads]
        return new{Dim, Neigh}(pos, rank, neighbours, threads, grid_size, periodicity, [], [], [], 0, threads_workloads)
    end
end


struct ProcessGridTopo{Dim, Neigh}
    use_mpi         :: Bool
    proc_size       :: NTuple{Dim, Int}
    periodicity     :: NTuple{Dim, Bool}
    channels        :: Vector{P2PChannel}
    processes       :: Vector{GridTopo{Dim, Neigh}}

    ProcessGridTopo(dim, proc_size, periodicity, use_mpi) = new{dim, 2*dim}(use_mpi, proc_size, periodicity, [], [])
end

Base.ndims(::ProcessGridTopo{Dim}) where {Dim} = Dim


struct SolverOptions
    real_block_xchg :: Bool
    num_cycles      :: Int
    max_sweeps      :: Int
    do_time_step    :: Bool
    do_halo_xchg    :: Bool
end


is_in_grid(grid_size, pos)     = all(1 .≤ Tuple(pos) .≤ Tuple(grid_size))
is_in_grid(grid_size, pos, ax) =    (1 .≤ Tuple(pos) .≤ Tuple(grid_size))[ax]

offset_along(ax, dim) = ntuple(i -> i == ax ? 1 : 0, dim)


function build_grid_topo!(grid::GridTopo{dim, neigh}) where {dim, neigh}
    null_interface = Interface(NULL_INTERFACE)
    remote_interface = Interface(REMOTE_BLOCK)
    null_neighbour = Block(NULL_BLOCK, ntuple(Returns(0), dim))

    resize!(grid.blocks, prod(grid.size))
    for pos in CartesianIndices(grid.size)
        idx = LinearIndices(grid.size)[pos]
        grid.blocks[idx] = Block(idx, pos)
    end

    for pos in CartesianIndices(grid.size)
        idx = LinearIndices(grid.size)[pos]
        blk_neighbours = Vector{Union{Block{dim, neigh}, RemoteBlock{dim}}}(undef, neigh)
        blk_interfaces = Vector{Interface}(undef, neigh)
        for ax in 1:dim, is_backwards in (true, false)
            side_idx       = (ax - 1) * 2 + (is_backwards ? 1 : 2)
            other_side_idx = (ax - 1) * 2 + (is_backwards ? 2 : 1)
            offset = offset_along(ax, dim) .* (is_backwards ? -1 : 1)
            neighbour_pos = pos + CartesianIndex(offset)

            if grid.periodicity[ax] && !is_in_grid(grid.size, neighbour_pos, ax)
                neighbour_pos = neighbour_pos + CartesianIndex(grid.size .* offset .* -1)
            end

            if is_in_grid(grid.size, neighbour_pos, ax)
                neighbour_idx = LinearIndices(grid.size)[neighbour_pos]
                blk_neighbours[side_idx] = grid.blocks[neighbour_idx]
                if isdefined(grid.blocks[neighbour_idx], :interfaces)
                    blk_interfaces[side_idx] = grid.blocks[neighbour_idx].interfaces[other_side_idx]
                else
                    interface = Interface(length(grid.interfaces) + 1)
                    push!(grid.interfaces, interface)
                    blk_interfaces[side_idx] = interface
                end
            elseif grid.neighbours[side_idx] != NULL_RANK
                # the index along the interface between the processes uniquely identifies the remote block
                tag = pos[mod1(ax + 1, dim)]
                remote_block = RemoteBlock(length(grid.remote_blocks) + 1, grid.neighbours[side_idx], neighbour_pos, tag)
                push!(grid.remote_blocks, remote_block)
                blk_neighbours[side_idx] = remote_block
                blk_interfaces[side_idx] = remote_interface
            else
                blk_neighbours[side_idx] = null_neighbour
                blk_interfaces[side_idx] = null_interface
            end
        end
        grid.blocks[idx].neighbours = Tuple(blk_neighbours)
        grid.blocks[idx].interfaces = Tuple(blk_interfaces)
    end
end


function assign_workloads!(grid)
    block_count = length(grid.blocks)
    blocks_per_thread = fld(block_count, grid.threads)
    remaining_blocks = block_count - grid.threads * blocks_per_thread

    # Assign to the n-th thread the `(n:n+1) .* blocks_per_thread` blocks.
    # The first `remaining_blocks` threads have one more block to even out the extra workload.
    for tid in 1:grid.threads
        prev_tids_blocks = blocks_per_thread * (tid - 1)
        tid_blocks = blocks_per_thread
        if tid > remaining_blocks
            prev_tids_blocks += remaining_blocks
        else
            prev_tids_blocks += tid - 1
            tid_blocks += 1
        end

        workload = [blk_idx for blk_idx in (1:tid_blocks) .+ prev_tids_blocks]
        for idx in workload
            grid.blocks[idx].tid = tid
        end

        grid.threads_workloads[tid] = workload
    end
end


function assign_channels!(proc_grid::ProcessGridTopo)
    # Two MPI processes with P2P communications between each other will share one channel.
    # Each channel will be big enough to handle messages from both sides.
    # Having one channel per remote block to remote block communication creates too many channels,
    # causing many issues, mainly having to explore depths far too large (>10^8).
    tot_channels = 0
    channels_topo = Dict(rank => Dict{Int, P2PChannel}() for rank in 1:length(proc_grid.processes))  # rank to {other rank to channel idx}
    for grid in proc_grid.processes
        for neigh_grid_rank in Iterators.filter(!=(NULL_RANK), grid.neighbours)
            neigh_grid = proc_grid.processes[neigh_grid_rank]
            neigh_grid.rank != neigh_grid_rank && error("wrong rank: expected=$neigh_grid_rank, got=$(neigh_grid.rank)")

            chan = get!(channels_topo[grid.rank], neigh_grid_rank) do
                tot_channels += 1
                chan = P2PChannel(tot_channels)
                push!(proc_grid.channels, chan)
                channels_topo[neigh_grid_rank][grid.rank] = chan
                return chan
            end

            for blk in grid.remote_blocks
                blk.rank != neigh_grid_rank && continue
                blk.send_chan.idx != -1 && continue  # the channels are already assigned

                neigh_blk_idx = findfirst(ob -> ob.rank == grid.rank && ob.tag == blk.tag, neigh_grid.remote_blocks)
                isnothing(neigh_blk_idx) && error("could not match remote block of grid $(grid.pos) (rank $(grid.rank)) at $(blk.pos) with grid at $(neigh_grid.pos) (rank $neigh_grid_rank)")
                neigh_blk = neigh_grid.remote_blocks[neigh_blk_idx]

                blk.send_chan = blk.recv_chan = neigh_blk.recv_chan = neigh_blk.send_chan = chan
                chan.usage += 2  # 1 per block
            end
        end

        # Check if all remote blocks were assigned
        for blk in grid.remote_blocks
            blk.rank == NULL_RANK && continue
            blk.send_chan.idx == -1 && error("unassigned remote block in $(grid.pos) at $(blk.pos) with rank $(blk.rank)")
        end
    end
end


function build_process_grid_topo(options)
    use_mpi = options["use_mpi"]
    threads = options["num_threads"]
    grid_size = options["grid"]
    proc_grid_size = options["processes"]
    proc_grid_periodicity = options["periodic"] .&& proc_grid_size .>  1  # periodicity between processes
    grid_periodicity      = options["periodic"] .&& proc_grid_size .== 1  # periodicity in each grid
    dim = length(grid_size)
    !use_mpi && prod(proc_grid_size) > 1 && error("'use_mpi' must be 'true' when using more than 1 MPI process")

    tot_blocks        = 0
    tot_interfaces    = 0
    tot_remote_blocks = 0
    proc_grid = ProcessGridTopo(dim, proc_grid_size, Tuple(proc_grid_periodicity), use_mpi)
    for (rank, proc_pos) in enumerate(CartesianIndices(proc_grid_size))
        # Get the ranks of neighbouring processes
        # 1D order: left, right
        # 2D order: left, right, bottom, top
        neighbour_procs = ntuple(2*dim) do i
            ax = (i - 1) ÷ 2 + 1
            offset = offset_along(ax, dim)
            mod1(i, 2) == 1 && (offset = offset .* -1)
            neighbour_pos = proc_pos + CartesianIndex(offset)

            if proc_grid_periodicity[ax] && !is_in_grid(proc_grid_size, neighbour_pos, ax)
                neighbour_pos = neighbour_pos + CartesianIndex(proc_grid_size .* offset .* -1)
            end

            if is_in_grid(proc_grid_size, neighbour_pos)
                return LinearIndices(proc_grid_size)[neighbour_pos]
            else
                return NULL_RANK
            end
        end

        grid = GridTopo{dim, 2*dim}(proc_pos, rank, threads, grid_size, Tuple(grid_periodicity))
        grid.neighbours = neighbour_procs
        build_grid_topo!(grid)
        assign_workloads!(grid)

        # Shift all indices to make them unique among all processes
        foreach(blk -> blk.idx += tot_blocks,        grid.blocks)
        foreach(int -> int.idx += tot_interfaces,    grid.interfaces)
        foreach(rmt -> rmt.idx += tot_remote_blocks, grid.remote_blocks)
        grid.workloads_offset = tot_blocks
        tot_blocks        += length(grid.blocks)
        tot_interfaces    += length(grid.interfaces)
        tot_remote_blocks += length(grid.remote_blocks)

        push!(proc_grid.processes, grid)
    end

    assign_channels!(proc_grid)

    return proc_grid
end


function get_solver_options(options)
    real_block_xchg = options["real_block_xchg"]
    num_cycles      = options["num_cycles"]
    max_sweeps      = options["max_sweeps"]
    do_time_step    = options["do_time_step"]
    do_halo_xchg    = options["do_halo_xchg"]
    return SolverOptions(real_block_xchg, num_cycles, max_sweeps, do_time_step, do_halo_xchg)
end


write_c_array(io::IO, array) = (print(io, '{'); join(io, array, ", "); print(io, '}'))
function write_idx_array(io::IO, array)
    print(io, '{')
    join(io, lpad.(ifelse.(in.(array, Ref((NULL_BLOCK, NULL_INTERFACE, REMOTE_BLOCK))), array, array .- 1), 3), ", ")
    print(io, '}')
end


function write_block(io::IO, block::Block, rank)
    # "{<idx>, <rank>, <pos>, <neighbours>, <interfaces>, <tid>},"
    print(io, '{', lpad(block.idx - 1, 3), ", ", rank - 1, ", ")
    write_c_array(io, block.pos .- 1)
    print(io, ", ")
    write_idx_array(io, getfield.(block.neighbours, :idx))
    print(io, ", ")
    write_idx_array(io, getfield.(block.interfaces, :idx))
    print(io, ", ", block.tid - 1, '}')
end


function write_process(io::IO, grid::GridTopo)
    # Get all remote blocks (and their channels) this thread will interact with
    remote_blocks = RemoteBlock[]
    send_channels = String[]
    recv_channels = String[]
    for blk in grid.blocks, neigh_blk in blk.neighbours
        !(neigh_blk isa RemoteBlock) && continue
        push!(remote_blocks, neigh_blk)
        push!(send_channels, "channel_" * string(neigh_blk.send_chan.idx-1))
        push!(recv_channels, "channel_" * string(neigh_blk.recv_chan.idx-1))
    end

    if !isempty(remote_blocks)
        channels_inits = map(enumerate(remote_blocks)) do (i, blk)
            "MPI_Request_init(remote_blocks[$(blk.idx-1)].req, $(send_channels[i]), $(recv_channels[i]));"
        end

        requests_inits = map(remote_blocks) do blk
            is_low_side = grid.rank < blk.rank
            if is_low_side
                "remote_blocks[$(blk.idx-1)].send_tag = $(blk.tag*2); remote_blocks[$(blk.idx-1)].recv_tag = $(blk.tag*2-1);"
            else
                "remote_blocks[$(blk.idx-1)].send_tag = $(blk.tag*2-1); remote_blocks[$(blk.idx-1)].recv_tag = $(blk.tag*2);"
            end
        end

        channels_inits = join(channels_inits, "\n    ") * "\n\n    " * join(requests_inits, "\n    ")
    else
        channels_inits = "// no channels for rank $(grid.rank-1)"
    end

    return channels_inits
end


function write_all_processes(io::IO, proc_grid::ProcessGridTopo)
    max_channel_usage = maximum(chan -> chan.usage, proc_grid.channels; init=0)
    if max_channel_usage ≥ 255
        error("Maximum channel usage ($max_channel_usage) exceeds the maximum value of a byte (255)")
    end

    num_pad   = floor(Int, log10(max(length(proc_grid.channels) - 1, 1))) + 1
    usage_pad = floor(Int, log10(max(max_channel_usage, 1))) + 1
    for chan in proc_grid.channels
        println(io, "chan channel_$(rpad(chan.idx-1, num_pad)) = [$(rpad(chan.usage, usage_pad))*MPI_SPIN_P2P_CHAN_SIZE] of { byte, byte };")
    end
    println(io)

    # Defines all proctypes for SPIN, and a function to initialize them
    ranks_inits = String[]
    for grid in proc_grid.processes
        rank_init = write_process(io, grid)
        rank_init = ":: (RANK == $(grid.rank-1)) -> {\n    $rank_init\n}"
        rank_init = replace(rank_init, '\n' => "\n    ")
        push!(ranks_inits, rank_init)
    end
    println(io, """
    inline init_rank()
    {
        d_step {
            if
            $(replace(join(ranks_inits, "\n    "), '\n' => "\n    "))
            :: else -> assert(false);
            fi
        }
    }
    """)

    # The proctype common to all processes, and a function to run all processes and threads
    print(io, """
    proctype solver_rank_thread(byte RANK; byte TID)
    {
        // This proctype is the base type of all threads of all MPI ranks
        // 'RANK' and 'TID' can be considered as alaways-defined variables in all non-init functions
        init_rank();
        solver_thread();
    }

    inline run_all_procs()
    {
        byte rank, tid;
        for (rank : 0 .. NUM_PROC-1) {
            for (tid : 0 .. NUM_THREADS-1) {
                run solver_rank_thread(rank, tid);
            }
        }
    }
    """)
end


function write_proc_grid(io_h::IO, io_c::IO, proc_grid::ProcessGridTopo)
    max_work            = maximum(g -> maximum(length, g.threads_workloads), proc_grid.processes)
    total_blocks        = sum(g -> length(g.blocks),        proc_grid.processes)
    total_interfaces    = sum(g -> length(g.interfaces),    proc_grid.processes)
    total_remote_blocks = sum(g -> length(g.remote_blocks), proc_grid.processes)

    interface_idx_limit = max(NULL_INTERFACE, REMOTE_BLOCK)
    total_blocks        ≥ NULL_BLOCK          && error("too many blocks! limit=$NULL_BLOCK, got=$total_blocks")
    total_interfaces    ≥ interface_idx_limit && error("too many interfaces! limit=$interface_idx_limit, got=$total_interfaces")
    total_remote_blocks ≥ NULL_BLOCK          && error("too many remote blocks! limit=$NULL_BLOCK, got=$total_remote_blocks")

    # The header is included "as is" in the Promela source, therefore it cannot have C definitions (a bit stupid I know)
    print(io_h, """
    #define PROC_GRID_STR       "$(proc_grid.proc_size)"
    #define GRID_SIZE_STR       "$(first(proc_grid.processes).size)"
    #define DIMS                $(ndims(proc_grid))
    #define NUM_NEIGHBOURS      $(2*ndims(proc_grid))
    #define TOTAL_BLOCKS        $(max(total_blocks, 1))
    #define TOTAL_INTERFACES    $(max(total_interfaces, 1))
    #define TOTAL_REMOTE_BLOCKS $(max(total_remote_blocks, 1))

    #define REMOTE_BLOCK   $(REMOTE_BLOCK)
    #define NULL_BLOCK     $(NULL_BLOCK)
    #define NULL_INTERFACE $(NULL_INTERFACE)
    #define NULL_RANK      $(NULL_RANK)

    #define USE_MPI       $(Int(proc_grid.use_mpi))
    #define NUM_PROC      $(prod(proc_grid.proc_size))
    #define NUM_THREADS   $(first(proc_grid.processes).threads)
    #define TOTAL_THREADS $(sum(g -> g.threads, proc_grid.processes))
    #define MAX_WORKLOAD  $max_work
    """)

    # Note: "byte" in Promela is converted to "uchar", itself an alias to "unsigned char"
    print(io_c, """
    typedef struct BlockTopo {
        uchar idx;                         // index in 'grid_topology'
        uchar rank;                        // assiociated MPI rank
        uchar pos[DIMS];                   // (x,y) position in the grid
        uchar neighbours[NUM_NEIGHBOURS];  // indexes of the neighbouring blocks ($(NULL_BLOCK) if none)
        uchar interfaces[NUM_NEIGHBOURS];  // indexes of the interfaces to the neighbouring blocks ($(NULL_INTERFACE) if none, $(REMOTE_BLOCK) if remote block)
        uchar tid;                         // thread associated to the block (local to the MPI rank)
    } BlockTopo;

    const BlockTopo grid_topology[TOTAL_BLOCKS] = {
    """)

    for grid in proc_grid.processes, blk in grid.blocks
        print(io_c, "    ")
        write_block(io_c, blk, grid.rank)
        println(io_c, ',')
    end

    print(io_c, """
    };

    typedef struct ThreadWorkload {
        uchar rank;
        uchar tid;
        uchar num_blocks;
        uchar blocks[MAX_WORKLOAD];  // indexes of blocks assigned to the thread
    } ThreadWorkload;

    const ThreadWorkload threads_workload[TOTAL_THREADS] = {
    """)

    for grid in proc_grid.processes, (tid, workload) in enumerate(grid.threads_workloads)
        norm_workload = map(1:max_work) do i
            if i ≤ length(workload)
                workload[i] + grid.workloads_offset
            else
                NULL_BLOCK
            end
        end
        print(io_c, "    { $(grid.rank-1), $(tid-1), $(length(workload)), ")
        write_idx_array(io_c, norm_workload)
        println(io_c, "},")
    end

    println(io_c, "};")
end


function write_solver_options(io_h::IO, solver::SolverOptions)
    println(io_h, """
    // solver options
    #define REAL_BLOCK_XCHG  $(Int(solver.real_block_xchg))
    #define NUM_CYCLES       $(solver.num_cycles)
    #define MAX_SWEEPS       $(solver.max_sweeps)
    #define DO_TIME_STEP     $(Int(solver.do_time_step))
    #define DO_HALO_EXCHANGE $(Int(solver.do_halo_xchg))
    """)
end


function write_proc_grid(filename, proc_grid::ProcessGridTopo, solver_options::SolverOptions)
    open(filename * ".h", "w") do header_file
        println(header_file, "#ifndef _PROC_GRID_H")
        println(header_file, "#define _PROC_GRID_H\n")

        open(filename * ".c", "w") do c_file
            println(c_file, "#include \"$(last(splitpath(filename)) * ".h")\"\n")
            write_proc_grid(header_file, c_file, proc_grid)
        end

        println(header_file)
        write_solver_options(header_file, solver_options)

        println(header_file, "#endif // _PROC_GRID_H")
    end

    open(filename * ".pml", "w") do promela_file
        write_all_processes(promela_file, proc_grid)
    end
end


print_solver_stats(proc_grid::ProcessGridTopo, solver::SolverOptions) = print_solver_stats(stdout, proc_grid, solver)
function print_solver_stats(io::IO, proc_grid::ProcessGridTopo, solver::SolverOptions)
    tot_blocks  = sum(g -> length(g.blocks), proc_grid.processes)
    tot_remotes = sum(g -> length(g.remote_blocks), proc_grid.processes)
    tot_threads = sum(g -> g.threads, proc_grid.processes)

    println(io, "Grid of ", proc_grid.proc_size, " processes ($(length(proc_grid.processes)) total):")
    println(io, " - $tot_blocks blocks")
    println(io, " - $tot_remotes remote blocks")
    println(io, " - $tot_threads threads (proctypes)")
    println(io, " - $(length(proc_grid.channels)) channels for exchanges")

    println(io, "Solver:")
    println(io, " - $(solver.num_cycles) cycles of $(solver.max_sweeps) sweeps ($(solver.num_cycles*solver.max_sweeps) total)")
    println(io, " - time step reduction: $(solver.do_time_step)")
    print(io,   " - halo exchange: ")
    if solver.do_halo_xchg && !solver.real_block_xchg
        print(io, "approximative")
    else
        print(io, solver.real_block_xchg)
    end
    println(io)
end


if !isinteractive()
    options = parse_arguments()
    proc_grid = build_process_grid_topo(options)
    solver_options = get_solver_options(options)
    write_proc_grid("grid_definition", proc_grid, solver_options)
    print_solver_stats(proc_grid, solver_options)
end
