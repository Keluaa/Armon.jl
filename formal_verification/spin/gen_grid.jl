
using TOML


const USAGE = """
Script to generate a static specification of a 2D block grid for the SPIN specification.

Usage:
    julia ./gen_grid.jl <grid_file.toml>

See './grids/default_grid_3x3_4threads.toml' for an explaination of each option.
"""

const DEFAULT_ARGS = Dict(
    "grid"        => [3, 3],
    "periodic"    => [false, false],
    "num_threads" => 4,
    "use_mpi"     => false,
    "processes"   => [1, 1],
)

const REMOTE_BLOCK   = 254
const NULL_BLOCK     = 255
const NULL_INTERFACE = 255


function parse_arguments()
    length(ARGS) != 1 && error("Expected 1 argument. Usage:\n" * USAGE)
    ARGS[1] in ("-h", "--help") && (println(USAGE); exit())

    options = try
        open(TOML.tryparse, ARGS[1], "r")
    catch e
        println("Could not parse '$(ARGS[1])'")
        rethrow(e)
    end

    grid_dim = length(options["grid"])
    if grid_dim != 2
        DEFAULT_ARGS["periodic"]  = map(Returns(false), 1:grid_dim)
        DEFAULT_ARGS["processes"] = map(Returns(false), 1:grid_dim)
    end

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


mutable struct Interface
    idx :: Int
end


mutable struct Block{Dim, Neigh}
    idx        :: Int
    pos        :: NTuple{Dim, Int}
    neighbours :: NTuple{Neigh, Block}
    interfaces :: NTuple{Neigh, Interface}
    tid        :: Int

    Block{D}(idx, pos) where {D} = new{D, 2D}(idx, pos)
    Block(idx, pos) = Block{length(pos)}(idx, pos)
end


mutable struct RemoteBlock
    idx       :: Int
    rank      :: Int
    other_idx :: Int
end


is_in_grid(grid_size, pos)     =  1 .≤ Tuple(pos) .≤ Tuple(grid_size)
is_in_grid(grid_size, pos, ax) = (1 .≤ Tuple(pos) .≤ Tuple(grid_size))[ax]

offset_along(ax, dim) = ntuple(i -> i == ax ? 1 : 0, dim)


function build_grid_topo(grid_size, periodicity)
    dim = length(grid_size)
    null_interface = Interface(0)
    null_neighbour = Block(0, ntuple(Returns(0), dim))

    interfaces = Vector{Interface}()
    blocks = Vector{Block{dim, 2*dim}}(undef, prod(grid_size))
    for pos in CartesianIndices(grid_size)
        idx = LinearIndices(grid_size)[pos]
        blocks[idx] = Block(idx, pos)
    end

    for pos in CartesianIndices(grid_size)
        idx = LinearIndices(grid_size)[pos]
        blk_neighbours = Vector{Block{dim, 2*dim}}(undef, 2*dim)
        blk_interfaces = Vector{Interface}(undef, 2*dim)
        for ax in 1:dim, is_backwards in (true, false)
            side_idx = (ax - 1) * 2 + (is_backwards ? 1 : 2)
            other_side_idx = side_idx + (is_backwards ? 1 : -1)
            offset = offset_along(ax, dim) .* (is_backwards ? -1 : 1)
            neighbour_pos = pos + CartesianIndex(offset)

            if periodicity[ax] && !is_in_grid(grid_size, neighbour_pos, ax)
                neighbour_pos = neighbour_pos + CartesianIndex(grid_size .* offset .* -1)
            end

            if is_in_grid(grid_size, neighbour_pos, ax)
                neighbour_idx = LinearIndices(grid_size)[neighbour_pos]
                blk_neighbours[side_idx] = blocks[neighbour_idx]
                if isdefined(blocks[neighbour_idx], :interfaces)
                    blk_interfaces[side_idx] = blocks[neighbour_idx].interfaces[other_side_idx]
                else
                    interface = Interface(length(interfaces) + 1)
                    push!(interfaces, interface)
                    blk_interfaces[side_idx] = interface
                end
            else
                blk_neighbours[side_idx] = null_neighbour
                blk_interfaces[side_idx] = null_interface
            end
        end
        blocks[idx].neighbours = Tuple(blk_neighbours)
        blocks[idx].interfaces = Tuple(blk_interfaces)
    end

    return blocks, interfaces
end


function assign_workloads(blocks, num_threads)
    block_count = length(blocks)
    blocks_per_thread = fld(block_count, num_threads)
    remaining_blocks = block_count - num_threads * blocks_per_thread

    # Assign to the n-th thread the `(n:n+1) .* blocks_per_thread` blocks.
    # The first `remaining_blocks` threads have one more block to even out the extra workload.
    threads_workload = map(1:num_threads) do tid
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
            blocks[idx].tid = tid
        end

        return workload
    end

    return threads_workload
end


function build_process_grid_topo(options)
    if !options["use_mpi"]
        blocks, interfaces = build_grid_topo(options["grid"], options["periodic"])
        threads_workload = assign_workloads(blocks, options["num_threads"])
        return blocks, interfaces, RemoteBlock[], threads_workload
    end

    proc_grid = options["processes"]
    proc_grid_periodicity = options["periodic"] .&& proc_grid .>  1  # periodicity between processes
    grid_periodicity      = options["periodic"] .&& proc_grid .== 1  # periodicity in each grid

    all_blocks        = Vector{Block{dim, 2*dim}}()
    all_interfaces    = Vector{Interface}()
    all_remote_blocks = Vector{RemoteBlock}()
    all_workloads     = Int[][]

    # TODO

    return all_blocks, all_interfaces, all_remote_blocks, all_workloads
end


write_c_array(io::IO, array) = (print(io, '{'); join(io, array, ", "); print(io, '}'))
function write_idx_array(io::IO, array, neg_value)
    print(io, '{')
    join(io, lpad.(ifelse.(array .== -1, neg_value, array), 3), ", ")
    print(io, '}')
end


function write_grid_topo(io_h::IO, io_c::IO, proc_grid, grid_size, blocks, remote_blocks, interfaces)
    # The header is included "as is" in the Promela source, therefore it cannot have C definitions (a bit stupid I know)
    println(io_h, """
    #define PROC_GRID_STR       "$(proc_grid)"
    #define GRID_SIZE_STR       "$(grid_size)"
    #define DIMS                $(length(grid_size))
    #define NUM_NEIGHBOURS      $(2*length(grid_size))
    #define TOTAL_BLOCKS        $(max(length(blocks), 1))
    #define TOTAL_INTERFACES    $(max(length(interfaces), 1))
    #define TOTAL_REMOTE_BLOCKS $(max(length(remote_blocks), 1))

    #define REMOTE_BLOCK   $(REMOTE_BLOCK)
    #define NULL_BLOCK     $(NULL_BLOCK)
    #define NULL_INTERFACE $(NULL_INTERFACE)
    """)

    # Note: "byte" in Promela is converted to "uchar", itself an alias to "unsigned char"
    print(io_c, """
    typedef struct BlockTopo {
        uchar pos[DIMS];                   // (x,y) position in the grid
        uchar neighbours[NUM_NEIGHBOURS];  // indexes of the neighbouring blocks ($(NULL_BLOCK) if none)
        uchar interfaces[NUM_NEIGHBOURS];  // indexes of the interfaces to the neighbouring blocks ($(NULL_INTERFACE) if none)
        uchar tid;                         // thread associated to the block
    } BlockTopo;

    const BlockTopo grid_topology[TOTAL_BLOCKS] = {
    """)

    for block in blocks
        # "{<pos>, <neighbours>, <interfaces>, <tid>},"
        print(io_c, "    {")
        write_c_array(io_c, block.pos .- 1)
        print(io_c, ", ")
        write_idx_array(io_c, getfield.(block.neighbours, :idx) .- 1, NULL_BLOCK)
        print(io_c, ", ")
        write_idx_array(io_c, getfield.(block.interfaces, :idx) .- 1, NULL_INTERFACE)
        println(io_c, ", ", block.tid - 1, "},")
    end

    println(io_c, "};\n")
end


function write_workload(io_h::IO, io_c::IO, proc_grid, threads_workload::Vector{Vector{Int}})
    max_work = maximum(length, threads_workload)

    println(io_h, """
    #define NUM_PROC     $(prod(proc_grid))
    #define NUM_THREADS  $(length(threads_workload))
    #define MAX_WORKLOAD $max_work
    """)

    print(io_c, """
    typedef struct ThreadWorkload {
        uchar tid;
        uchar num_blocks;
        uchar blocks[MAX_WORKLOAD];  // indexes of blocks assigned to the thread
    } ThreadWorkload;

    const ThreadWorkload threads_workload[NUM_THREADS] = {
    """)

    for (tid, workload) in enumerate(threads_workload)
        norm_workload = map(1:max_work) do i; get(workload, i, 0) end
        print(io_c, "    { $(tid-1), $(length(workload)), ")
        write_idx_array(io_c, norm_workload .- 1, NULL_BLOCK)
        println(io_c, "},")
    end

    println(io_c, "};")
end


function write_grid_to_file(filename, proc_grid, grid_size, blocks, interfaces, remote_blocks, threads_workload)
    open(filename * ".h", "w") do header_file
        open(filename * ".c", "w") do c_file
            write_grid_topo(header_file, c_file, proc_grid, grid_size, blocks, remote_blocks, interfaces)
            write_workload(header_file, c_file, proc_grid, threads_workload)
        end
    end
end


if !isinteractive()
    options = parse_arguments()
    write_grid_to_file("grid_definition", options["processes"], options["grid"], build_process_grid_topo(options)...)
end
