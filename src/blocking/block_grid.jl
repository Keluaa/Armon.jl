
"""
    BlockGrid{T, DeviceA, HostA, BufferA, Ghost, BlockSize, Device, SolverState}

Stores [`TaskBlock`](@ref)s on the `Device` and host memory, in a grid.

[`LocalTaskBlock`](@ref) are stored separately depending on if they have a [`StaticBSize`](@ref) of
`BlockSize` (in `blocks`) or if they have a [`DynamicBSize`](@ref) (in `edge_blocks`).

Blocks have `Ghost` cells padding their real cells. This is included in their [`block_size`](@ref).
A block cannot have a number of real cells along an axis smaller than the number of ghost cells,
unless there is only a single block in the grid along that axis.

"Edge blocks" are blocks located on the right and/or top edge of the grid. They exist in order to
handle domains with dimensions which are not multiples of the block size.

`DeviceA` and `HostA` are `AbstractArray` types for the device and host respectively.

`BufferA` is the type of storage used for MPI buffers. MPI buffers are homogenous: they are either
all on the host or all on the device.
"""
struct BlockGrid{
    T,
    DeviceArray <: AbstractArray{T},
    HostArray   <: AbstractArray{T},
    BufferArray <: AbstractArray{T},
    Ghost,
    BS          <: StaticBSize{<:Any, Ghost},
    SState      <: SolverState,
    Device
}
    grid_size          :: NTuple{2, Int}  # Size of the grid, including all local blocks
    static_sized_grid  :: NTuple{2, Int}  # Size of the grid of statically sized local blocks
    cell_size          :: NTuple{2, Int}  # Number of real cells in each direction
    edge_size          :: NTuple{2, Int}  # Number of real cells in edge blocks in each direction (only along non-edge directions)
    periodic           :: NTuple{2, Bool}  # If the grid is periodic along an axis
    device             :: Device
    global_dt          :: GlobalTimeStep{T}
    blocks             :: Vector{LocalTaskBlock{DeviceArray, HostArray, BS, SState}}
    edge_blocks        :: Vector{LocalTaskBlock{DeviceArray, HostArray, DynamicBSize{Ghost}, SState}}
    remote_blocks      :: Vector{RemoteTaskBlock{BufferArray}}
    threads_workload   :: Vector{Vector{CartesianIndex{2}}}  # `tid => block index` map for all threads, distributing each block to each thread
    threads_logs       :: Vector{Vector{ThreadLogEvent}}
end


function BlockGrid(params::ArmonParameters{T}) where {T}
    grid_size, static_sized_grid, remainder_block_size = grid_dimensions(params)
    cell_size = params.N
    static_sized_block_count = prod(static_sized_grid)
    dyn_sized_block_count = prod(grid_size) - static_sized_block_count

    device_array = device_array_type(params.device){T, 1}
    host_array = host_array_type(params.device){T, 1}

    # A bit janky, but required as `device_array` or `host_array` might still be incomplete
    device_array = Core.Compiler.return_type(device_array, Tuple{UndefInitializer, Int})
    host_array = Core.Compiler.return_type(host_array, Tuple{UndefInitializer, Int})

    global_dt = GlobalTimeStep{T}()
    state_type = typeof(SolverState(params, global_dt))

    ghost = params.nghost

    # Container for blocks with a static size
    static_size = StaticBSize(params.block_size, ghost)
    blocks = Vector{LocalTaskBlock{device_array, host_array, typeof(static_size), state_type}}(undef, static_sized_block_count)

    # Container for blocks on the edges, with a non-uniform size
    edge_blocks = Vector{LocalTaskBlock{device_array, host_array, DynamicBSize{ghost}, state_type}}(undef, dyn_sized_block_count)

    # Container for remote blocks, neighbours of blocks on the edges. Corners are excluded.
    buffer_array = params.gpu_aware ? device_array : host_array
    grid_perimeter = sum(grid_size) * length(grid_size)  # (nx+ny) * 2
    remote_blocks = Vector{RemoteTaskBlock{buffer_array}}(undef, grid_perimeter)

    threads_workload = thread_workload_distribution(params)

    log_size = params.log_blocks ? min(params.maxcycle, 1000) : 0
    threads_logs = map(1:Threads.nthreads()) do _
        logs = Vector{ThreadLogEvent}()
        sizehint!(logs, log_size)
        return logs
    end

    # Main grid container
    edge_size = remainder_block_size .- 2*ghost
    periodic_grid = params.periodic .&& params.proc_dims .== 1  # `true` only if the grid itself should be periodic
    grid = BlockGrid{
        T, device_array, host_array, buffer_array,
        ghost, typeof(static_size),
        state_type, typeof(params.device)
    }(
        grid_size, static_sized_grid, cell_size, edge_size, periodic_grid, params.device, global_dt,
        blocks, edge_blocks, remote_blocks, threads_workload, threads_logs
    )

    # Allocate all local and remote blocks
    # Non-static (edge) blocks are placed on the right and top sides.
    # Remote blocks are placed on the edge of the grid.
    # Multithreading is necessary here in order to guarentee that no array is shared between two
    # NUMA node (when we move the pages afterward), which can happen when allocations are done
    # sequentially.
    inner_grid = CartesianIndices(static_sized_grid)
    device_kwargs = alloc_device_kwargs(params)
    host_kwargs = alloc_host_kwargs(params)
    @threaded :outside_kernel for _ in 1:length(threads_workload)
        tid = Threads.threadid()
        for pos in threads_workload[tid]
            # Static block for the inner grid, edge block otherwise
            local_block = if pos in inner_grid
                idx = block_idx(grid, pos)
                blk_size = static_size
                blk_state = SolverState(params, global_dt)
                blocks[idx] = eltype(blocks)(blk_size, pos, blk_state, device_kwargs, host_kwargs)
            else
                idx = edge_block_idx(grid, pos)
                blk_size = DynamicBSize(
                    ifelse.(Tuple(pos) .== grid_size .&& remainder_block_size .!= 0,
                        remainder_block_size,
                        params.block_size
                    ),
                    ghost
                )
                blk_state = SolverState(params, global_dt)
                edge_blocks[idx] = eltype(edge_blocks)(blk_size, pos, blk_state, device_kwargs, host_kwargs)
            end

            # Create the neighbouring remote blocks if we are at the border of the grid
            is_border_block = !in_grid(2, pos, grid_size .- 1)
            is_border_block && for side in instances(Side.T)
                # Position of the neighbouring block in the grid
                remote_blk_pos = pos + CartesianIndex(offset_to(side))
                in_grid(remote_blk_pos, grid_size) && continue

                # Neighbour MPI rank. If this side is periodic and the proc grid side is 1 along
                # this axis then it is the rank of this sub-domain. In this case we still create the
                # remote block but the local block will not connect to it.
                neighbour = neighbour_at(params, side)

                remote_block = if has_neighbour(params, side) && neighbour != params.rank
                    # The buffer must be the same size as the side of our block which is a neighbour to...
                    buffer_size = real_face_size(local_block.size, side)
                    # ...for each variable to communicate of each ghost cell
                    buffer_size *= length(comm_vars()) * params.nghost

                    global_pos = CartesianIndex(params.cart_coords .+ offset_to(side))  # pos in the cart_comm

                    RemoteTaskBlock{buffer_array}(buffer_size, remote_blk_pos, neighbour, global_pos, params.cart_comm, side)
                else
                    # "Fake" remote block for non-existant neighbour at the edge of the global domain
                    RemoteTaskBlock{buffer_array}(remote_blk_pos)
                end

                remote_blocks[remote_block_idx(grid, remote_blk_pos)] = remote_block
            end
        end
    end

    # Initialize all block neighbours references and exchanges
    for idx in CartesianIndex(1, 1):CartesianIndex(grid_size)
        left_idx   = idx + CartesianIndex(offset_to(Side.Left))
        right_idx  = idx + CartesianIndex(offset_to(Side.Right))
        bottom_idx = idx + CartesianIndex(offset_to(Side.Bottom))
        top_idx    = idx + CartesianIndex(offset_to(Side.Top))

        this_block   = block_at(grid, idx)
        left_block   = block_at(grid, left_idx)
        right_block  = block_at(grid, right_idx)
        bottom_block = block_at(grid, bottom_idx)
        top_block    = block_at(grid, top_idx)
        this_block.neighbours = Neighbours{TaskBlock}((left_block, right_block, bottom_block, top_block))

        # TODO: there might be some optimisations possible when a block is connected to itself, which
        # allow skipping most of the block interface logic

        # Blocks sharing a side must share the same `BlockInterface`
        this_block.exchanges = Neighbours{BlockInterface}((
            isdefined(left_block,   :exchanges) ? left_block.exchanges[Int(Side.Right)] : BlockInterface(),
            isdefined(right_block,  :exchanges) ? right_block.exchanges[Int(Side.Left)] : BlockInterface(),
            isdefined(bottom_block, :exchanges) ? bottom_block.exchanges[Int(Side.Top)] : BlockInterface(),
            isdefined(top_block,    :exchanges) ? top_block.exchanges[Int(Side.Bottom)] : BlockInterface(),
        ))

        for blk in this_block.neighbours
            blk isa RemoteTaskBlock || continue
            blk.neighbour = this_block
        end
    end

    return grid
end


"""
    EdgeBlockRegions(grid::BlockGrid; real_sizes=false)
    EdgeBlockRegions(
        grid_size::NTuple{D, Int}, static_sized_grid::NTuple{D, Int},
        block_size::NTuple{D, Int}, remainder_block_size::NTuple{D, Int}; ghosts=0
    )

Iterator over edge blocks, their positions and their size in each edge region of the `grid`.

`ghosts` only affects the size of edge blocks: `B .- 2*ghosts`. If `real_sizes == true`, then `ghosts`
is `ghosts(grid)`, therefore the real block size of edge blocks is returned.

There is a maximum of `2^D-1` edge regions in a grid
([see this explanation](https://en.wikipedia.org/wiki/Binomial_theorem#Geometric_explanation)).
When `grid_size[i] != static_sized_grid[i]`, there cannot be an edge region, hence there are exactly
`2^sum(grid_size .!= static_sized_grid)-1` regions.

```jldoctest
julia> collect(Armon.EdgeBlockRegions((5, 5), (4, 4), (32, 32), (16, 16)))
3-element Vector{Tuple{Int64, CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}, Tuple{Int64, Int64}}}:
 (4, CartesianIndices((5:5, 1:4)), (16, 32))
 (4, CartesianIndices((1:4, 5:5)), (32, 16))
 (1, CartesianIndices((5:5, 5:5)), (16, 16))

julia> collect(Armon.EdgeBlockRegions((5, 5, 5), (4, 5, 4), (32, 32, 32), (16, 16, 16)))
3-element Vector{Tuple{Int64, CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, Tuple{Int64, Int64, Int64}}}:
 (20, CartesianIndices((5:5, 1:5, 1:4)), (16, 32, 32))
 (20, CartesianIndices((1:4, 1:5, 5:5)), (32, 32, 16))
 (5, CartesianIndices((5:5, 1:5, 5:5)), (16, 32, 16))
```
"""
struct EdgeBlockRegions{D, AE}
    axis_edges :: AE
    iter       :: Iterators.ProductIterator{AE}
    grid_size  :: NTuple{D, Int}
    static_sized_grid    :: NTuple{D, Int}
    remainder_block_size :: NTuple{D, Int}
    block_size :: NTuple{D, Int}
    ghosts     :: Int

    function EdgeBlockRegions(
        grid_size::NTuple{D, Int}, static_sized_grid::NTuple{D, Int},
        block_size::NTuple{D, Int}, remainder_block_size::NTuple{D, Int};
        ghosts::Int=0
    ) where {D}
        # `false` => when axis `d` has no edge region, `true` => when axis `d` has an edge region
        # By iterating over `(false, true)` for axes with a non-empty edge, and `(false,)` without
        # an edge, we cover all possible regions.
        axis_edges = ifelse.(grid_size .!= static_sized_grid, Ref((false, true)), Ref((false,)))
        iter = Iterators.product(axis_edges...)
        return new{D, typeof(axis_edges)}(
            axis_edges, iter,
            grid_size, static_sized_grid, remainder_block_size, block_size, ghosts
        )
    end
end

EdgeBlockRegions(grid::BlockGrid; real_sizes=false) =
    EdgeBlockRegions(
        grid.grid_size, grid.static_sized_grid, static_block_size(grid), grid.edge_size .+ 2*ghosts(grid);
        ghosts=real_sizes ? ghosts(grid) : 0
    )

Base.eltype(::Type{<:EdgeBlockRegions{D}}) where {D} = Tuple{Int, CartesianIndices{D, NTuple{D, UnitRange{Int}}}, NTuple{D, Int}}
Base.length(it::EdgeBlockRegions) = 2^sum(it.grid_size .!= it.static_sized_grid) - 1
Base.size(it::EdgeBlockRegions) = (length(it),)

function Base.iterate(ebr::EdgeBlockRegions, state=0)
    if state == 0
        s0 = iterate(ebr.iter)  # Discard the first element, which is all `false`, i.e. the static sized region
        s0 === nothing && return nothing
        _, state = s0
    end

    s = iterate(ebr.iter, state)
    s === nothing && return nothing
    side_states, s = s

    # Use `grid_size - 1` over `static_sized_grid` to handle the edge case with multiple edge blocks
    # but no static blocks.
    region_block_count = prod(ifelse.(side_states, 1, max.(ebr.static_sized_grid, ebr.grid_size .- 1)))
    first_pos = CartesianIndex(ifelse.(side_states, ebr.grid_size, 1))
    last_pos  = CartesianIndex(ifelse.(side_states, ebr.grid_size, max.(ebr.static_sized_grid, ebr.grid_size .- 1)))
    region_block_size = ifelse.(side_states, ebr.remainder_block_size, ebr.block_size) .- 2*ebr.ghosts

    return (region_block_count, first_pos:last_pos, region_block_size), s
end


"""
    RemoteBlockRegions(grid::BlockGrid)
    RemoteBlockRegions(grid_size::NTuple{D, Int})

Iterator of all remote block positions in each region (the faces of the `grid`).
There is always `2*D` regions.

```jldoctest
julia> collect(Armon.RemoteBlockRegions((5, 3)))
4-element Vector{Tuple{Armon.Side.T, CartesianIndices{2, Tuple{UnitRange{Int64}, UnitRange{Int64}}}}}:
 (Armon.Side.Left, CartesianIndices((0:0, 1:3)))
 (Armon.Side.Right, CartesianIndices((6:6, 1:3)))
 (Armon.Side.Bottom, CartesianIndices((1:5, 0:0)))
 (Armon.Side.Top, CartesianIndices((1:5, 4:4)))
```
"""
struct RemoteBlockRegions{D}
    grid_size :: NTuple{D, Int}

    RemoteBlockRegions(grid_size::NTuple{D, Int}) where {D} = new{D}(grid_size)
end

RemoteBlockRegions(grid::BlockGrid) = RemoteBlockRegions(grid.grid_size)

Base.eltype(::Type{<:RemoteBlockRegions{D}}) where {D} = Tuple{Side.T, CartesianIndices{D, NTuple{D, UnitRange{Int}}}}
Base.length(::RemoteBlockRegions{D}) where {D} = 2*D
Base.size(it::RemoteBlockRegions) = (length(it),)

function Base.iterate(rbr::RemoteBlockRegions{D}, state=0) where {D}
    if state == 0
        # TODO: dimension agnostic
        s = iterate(instances(Side.T))
    else
        s = iterate(instances(Side.T), state)
    end

    s === nothing && return nothing
    side, s = s

    o = offset_to(side)
    first_idx = CartesianIndex(ifelse.(o .== 0,             1, ifelse.(o .> 0, rbr.grid_size .+ 1, 0)))
    last_idx  = CartesianIndex(ifelse.(o .== 0, rbr.grid_size, ifelse.(o .> 0, rbr.grid_size .+ 1, 0)))

    return (side, first_idx:last_idx), s
end


"""
    grid_dimensions(params::ArmonParameters)
    grid_dimensions(block_size::NTuple{D, Int}, domain_size::NTuple{D, Int}, ghost::Int) where {D}

Returns the dimensions of the grid in the form `(grid_size, static_sized_grid, remainder_block_size)`
from the `block_size` (the size of blocks in the `static_sized_grid`), the `domain_size` (number of
real cells) and the number of `ghost` cells, common to all blocks.

`grid_size` is the `static_sized_grid` including the edge blocks.
Edge blocks along the axis `d` have a size of `remainder_block_size[d]` along `d` only if they are
at the edge of the grid along `d`, or `block_size` otherwise.

`block_size` includes the `ghost` cells in its dimensions, which must all be greater than `2*ghost`.

If `prod(block_size) == 0`, then `block_size` is ignored and the grid is made of only a single block
of size `domain_size`.

!!! note

    Blocks must not be smaller than `ghost`, therefore edge blocks might be made bigger than
    `block_size`.

!!! note
    
    In case `domain_size` is smaller than `block_size .- 2*ghost` along any axis, the grid will
    contain only edge blocks.
"""
grid_dimensions(params::ArmonParameters) = grid_dimensions(params.block_size, params.N, params.nghost)

function grid_dimensions(block_size::NTuple{D, Int}, domain_size::NTuple{D, Int}, ghost::Int) where {D}
    if prod(block_size) == 0
        # No cache blocking: only one dynamic block for the whole grid
        return ntuple(Returns(1), D), ntuple(Returns(0), D), domain_size .+ 2*ghost
    end

    if any(block_size .< 3*ghost)
        # 2 ghost regions in all axes + the real cell region should have enough ghost cells to
        # prevent one block from depending on cells two (or more) blocks away, therefore `≥3*ghost`.
        solver_error(:config, "block size $block_size is too small for $ghost ghosts cells: \
                               it should be `≥3*ghosts` for all axes, got: $block_size .< $(3*ghost)")
    end

    # `block_size` includes the number of ghost cells, while `domain_size` is in real cells.
    real_block_size = block_size .- 2*ghost  # number of real cells in a block in each dim

    grid_size = domain_size .÷ real_block_size
    remainder = domain_size .% real_block_size

    # Edge-case: the remainder is smaller than the number of ghost cells. This creates problems during
    # the halo-exchange since ghosts cells at the edge of the static sized grid depend on all real cells
    # of the neighbouring edge block BUT ALSO the ghosts of that same edge block, themselves depending
    # on a remote block (or on the global border conditions).
    # To solve this, those remaining cells are merged with the neighbouring static sized blocks, by
    # creating more edge blocks with a size bigger than `block_size`.
    remainder_too_small = 0 .< remainder .< ghost .&& grid_size .> 0
    remainder = remainder .+ real_block_size .* remainder_too_small

    # This deals with the edge case where the `domain_size` is smaller than `block_size` in one axis,
    # resulting in all blocks to be edge blocks, but with a grid potentially bigger than `1×1`.
    grid_size = (domain_size .- remainder) .÷ real_block_size

    static_sized_grid = grid_size  # Grid of blocks with a static size of `block_size`
    prod(static_sized_grid) == 0 && (static_sized_grid = ntuple(Returns(0), D))

    grid_size = grid_size .+ (remainder .> 0)  # Include edge blocks in the `grid_size`

    if prod(static_sized_grid) == 0
        # Edge case: for grids made only of edge blocks, axes with `remainder[ax] == 0` should instead
        # be `real_block_size`.
        remainder = ifelse.(remainder .== 0, real_block_size, remainder)
    end

    if any(grid_size .< 1)
        solver_error(:config, "could not partition $domain_size domain using $block_size blocks with $ghost ghost cells")
    end

    # (x, y) size of edge blocks, including ghost cells. A block for a X axis edge will have a size
    # of (x, block_size.y), and for the Y axis edge the size would be (block_size.x, y). In the edge
    # corner the block would be of size (x, y). Axes without an edge region have a remainder of 0.
    remainder_block_size = ifelse.(remainder .> 0, remainder .+ 2 * ghost, 0)

    # Sanity checks
    # `static blocks + edge blocks = whole grid`
    if prod(static_sized_grid) > 0
        valid_grid = all(static_sized_grid .+ (remainder_block_size .> 0) .== grid_size)
    else
        valid_grid = all(remainder_block_size .> 0)
    end
    # `real cells in static blocks + real cells in edge blocks = real domain`
    static_block_cells = prod(static_sized_grid) * prod(real_block_size)
    r = EdgeBlockRegions(grid_size, static_sized_grid, block_size, remainder_block_size; ghosts=ghost)
    edge_block_cells = sum(first.(r) .* prod.(last.(r)))
    valid_domain = static_block_cells + edge_block_cells == prod(domain_size)

    if !(valid_grid && valid_domain)
        err_msg = "failed to create a valid grid"
        !valid_grid   && (err_msg *= " (invalid grid)")
        !valid_domain && (err_msg *= " (invalid domain)")
        solver_error(:config, "$err_msg: \
                               block_size=$block_size, domain_size=$domain_size, ghost=$ghost, \
                               grid_size=$grid_size, static_sized_grid=$static_sized_grid, \
                               remainder_block_size=$remainder_block_size")
    end

    return grid_size, static_sized_grid, remainder_block_size
end


"""
    block_idx(grid::BlockGrid, idx::CartesianIndex)

Linear index in `grid.blocks` of the block at `idx` in the statically-sized grid.
"""
block_idx(grid::BlockGrid, idx::CartesianIndex) =
    LinearIndices(CartesianIndices(grid.static_sized_grid))[idx]


"""
    edge_block_idx(grid::BlockGrid, idx::CartesianIndex)

Linear index in `grid.edge_blocks` of the block at `idx`, along the (dynamically-sized) edges of the
`grid`.
"""
function edge_block_idx(grid::BlockGrid, idx::CartesianIndex, ::Val{:slow})
    offset = 0
    for (block_count, blocks_pos, _) in EdgeBlockRegions(grid)
        if idx in blocks_pos
            pos_idx = idx - first(blocks_pos) + oneunit(idx)  # `block_pos[pos_idx] == idx`
            return LinearIndices(blocks_pos)[pos_idx] + offset
        end
        offset += block_count
    end
    error("Block index $idx is not at the edge of the grid $grid")
end


@generated function edge_block_idx(grid::BlockGrid, idx::CartesianIndex{N}) where {N}
    if N > 8
        # Don't exceed 255 if statements. Above that, unrolling the loop seems too heavy to be more efficient.
        return quote edge_block_idx(grid, idx, Val{:slow}()) end
    end

    # Here we "simply" unroll the iteration over `EdgeBlockRegions` of the slow version of `edge_block_idx`
    # There is a lot of additional optimizations the compiler is able to do in this version, but the
    # most important part is the fact that this version requires 0 allocations.
    # For N=2, we go from >1µs to <10ns in execution time: a x100 improvement. Since this function is
    # can be used very frequently in `solver_async_cycle`, it is a worthwhile optimization.

    edge_regions = repeat([(false, true)], N)
    iter = Iterators.product(edge_regions...)
    iter = Iterators.drop(iter, 1)

    expr = quote
        grid_size = grid.grid_size
        static_sized_grid = grid.static_sized_grid
        offset = 0
    end

    for side_states in iter
        cond = Expr(:&&)  # crazy fact: the && expression gives coherent results even with 0 or 1 arguments
        for (i, side_state) in enumerate(side_states)
            !side_state && continue
            push!(cond.args, :(grid_size[$i] != static_sized_grid[$i]))
        end

        side_states_t = Tuple(side_states)
        test_edge_region = quote
            if $cond
                # Equivalent to a subset of `Base.iterate(::EdgeBlockRegions)` (see above)
                first_pos = CartesianIndex(ifelse.($side_states_t, grid_size, 1))
                last_pos  = CartesianIndex(ifelse.($side_states_t, grid_size, max.(static_sized_grid, grid_size .- 1)))
                edge_block_range = first_pos:last_pos
                if idx in edge_block_range
                    pos_idx = idx - first(edge_block_range) + oneunit(idx)  # `block_pos[pos_idx] == idx`
                    return LinearIndices(edge_block_range)[pos_idx] + offset
                else
                    offset += length(edge_block_range)
                end
            end
        end

        push!(expr.args, test_edge_region)
    end

    push!(expr.args, quote
        error("Block index $idx is not at the edge of the grid $grid")
    end)

    return expr
end


"""
    remote_block_idx(grid::BlockGrid, idx::CartesianIndex)

Linear index in `grid.remote_blocks` of the remote block at `idx` in the `grid`.
"""
function remote_block_idx(grid::BlockGrid, idx::CartesianIndex)
    offset = 0
    for (_, blocks_pos) in RemoteBlockRegions(grid)
        if idx in blocks_pos
            pos_idx = idx - first(blocks_pos) + oneunit(idx)  # `block_pos[pos_idx] == idx`
            return LinearIndices(blocks_pos)[pos_idx] + offset
        end
        offset += length(blocks_pos)
    end
    error("Block index $idx is not on the remote edges of the grid")
end


mirror_idx(grid::BlockGrid, idx::CartesianIndex) = CartesianIndex(mirror_idx(grid, Tuple(idx)))
function mirror_idx(grid::BlockGrid, idx::NTuple{2, Int})
    return ifelse.(grid.periodic,
        ifelse.(idx .== 0, grid.grid_size,
            ifelse.(idx .== grid.grid_size .+ 1, 1, idx)
        ), idx)
end


"""
    block_at(grid::BlockGrid, idx::CartesianIndex)

The [`TaskBlock`](@ref) at position `idx` in the `grid`.
"""
function block_at(grid::BlockGrid, idx::CartesianIndex)
    try
        if in_grid(idx, grid.static_sized_grid)
            return grid.blocks[block_idx(grid, idx)]
        elseif in_grid(idx, grid.grid_size)
            return grid.edge_blocks[edge_block_idx(grid, idx)]
        else
            periodic_axes = grid.periodic .&& (Tuple(idx) .== 0 .|| Tuple(idx) .== grid.grid_size .+ 1)
            if any(periodic_axes)
                # At a periodic edge of a grid, we return the local block of the other side instead
                # of a remote block, since the local blocks shouldn't link to them.
                return block_at(grid, mirror_idx(grid, idx))
            else
                return grid.remote_blocks[remote_block_idx(grid, idx)]
            end
        end
    catch e
        # TODO: very rarely, I get some UndefRefError and I don't know why
        !(e isa UndefRefError) && rethrow(e)
        kind_str = if in_grid(idx, grid.static_sized_grid)
            "static"
        elseif in_grid(idx, grid.grid_size)
            "edge"
        else
            "remote"
        end
        println(stderr, "Error while accessing ", Tuple(idx), " as a ", kind_str, " block")
        rethrow(e)
    end
end


"""
    block_size_at(grid::BlockGrid, idx)
    block_size_at(idx, grid_size, static_sized_grid, block_size, remainder_block_size, ghosts)

Theoretical size of block at `idx` in `grid`.
If `ghosts == 0`, ghosts cells will be included in the size.
"""
function block_size_at(grid::BlockGrid, idx)
    return block_size_at(idx,
        grid.grid_size, grid.static_sized_grid, static_block_size(grid),
        grid.edge_size .+ 2*ghosts(grid), ghosts(grid)
    )
end


function block_size_at(idx, grid_size, static_sized_grid, block_size, remainder_block_size, ghosts)
    if in_grid(idx, static_sized_grid)
        return block_size .- 2*ghosts
    else
        edge_block_size = ifelse.(
            in_grid.(Ref(idx), max.(static_sized_grid, grid_size .- 1), instances(Axis.T)),  # TODO: replace with `axes_of(dim)`
            block_size, remainder_block_size
        )
        return edge_block_size .- 2*ghosts
    end
end


device_array_type(::ObjOrType{BlockGrid{<:Any, D}}) where {D} = D
host_array_type(::ObjOrType{BlockGrid{<:Any, <:Any, H}}) where {H} = H
buffer_array_type(::ObjOrType{BlockGrid{<:Any, <:Any, <:Any, B}}) where {B} = B
ghosts(::ObjOrType{BlockGrid{<:Any, <:Any, <:Any, <:Any, Ghost}}) where {Ghost} = Ghost
static_block_size(::ObjOrType{BlockGrid{<:Any, <:Any, <:Any, <:Any, G, BS}}) where {G, BS} = block_size(BS)
real_block_size(::ObjOrType{BlockGrid{<:Any, <:Any, <:Any, <:Any, G, BS}}) where {G, BS} = real_block_size(BS)


"""
    all_blocks(grid::BlockGrid)

Simple iterator over all blocks of the grid, excluding [`RemoteTaskBlock`](@ref)s.
"""
all_blocks(grid::BlockGrid) = Iterators.flatten((grid.blocks, grid.edge_blocks))


"""
    device_is_host(::BlockGrid{T, D, H})
    device_is_host(::Type{<:BlockGrid{T, D, H}})

`true` if the device is the host, i.e. device blocks and host blocks are the same (and `D == H`).
"""
device_is_host(::ObjOrType{BlockGrid{<:Any, D, H}}) where {D, H} = D === H


"""
    buffers_on_device(::BlockGrid)
    buffers_on_device(::Type{<:BlockGrid})

`true` if the communication buffers are stored on the device, allowing direct transfers without
passing through the host (GPU-aware communication).
"""
buffers_on_device(::ObjOrType{BlockGrid{<:Any, D, H, B}}) where {D, H, B} = D === B


function reset!(grid::BlockGrid, params::ArmonParameters)
    reset!(grid.global_dt, params, prod(grid.grid_size))
    foreach(empty!, grid.threads_logs)
    for blk in all_blocks(grid)
        reset!(blk)
    end
end


"""
    first_state(grid::BlockGrid)

A [`SolverState`](@ref) which can be used as a global state when outside of a solver cycle.
It belongs to the first device block.
"""
first_state(grid::BlockGrid) = first(all_blocks(grid)).state


function ThreadLogEvent(grid::BlockGrid, tid, step_count, no_progress_count, stop_count, mpi_waits, wait_time, cycle_time)
    tw = grid.threads_workload[tid]
    cycle = isempty(tw) ? first_state(grid).cycle : block_at(grid, first(tw)).state.cycle
    blk_count = min(length(tw), typemax(Int16))
    return ThreadLogEvent(
        cycle, blk_count,
        stop_count, mpi_waits,
        step_count, no_progress_count,
        wait_time, cycle_time
    )
end

push_log!(grid::BlockGrid, tid, thread_log::ThreadLogEvent) = push!(grid.threads_logs[tid], thread_log)


"""
    memory_required(params::ArmonParameters)

`(device_memory, host_memory)` required for `params`.

MPI buffers size are included in the appropriate field depending on `params.gpu_aware`.
`params.use_MPI` and `params.neighbours` is taken into account.

If `device_is_host`, then, `device_memory` only includes memory required by data arrays and MPI buffers.
"""
function memory_required(params::ArmonParameters{T}) where {T}
    device_array = device_array_type(params.device){T, 1}
    host_array = host_array_type(params.device){T, 1}

    # A bit janky, but required as `device_array` or `host_array` might still be incomplete
    device_array = Core.Compiler.return_type(device_array, Tuple{UndefInitializer, Int})
    host_array = Core.Compiler.return_type(host_array, Tuple{UndefInitializer, Int})
    device_is_host = host_array == device_array
    buffer_array = params.gpu_aware ? device_array : host_array

    solver_state_type = typeof(SolverState(params, GlobalTimeStep{T}()))

    arrays_byte_count, MPI_buffer_byte_count, host_overhead = memory_required(
        params.N, params.block_size, params.nghost,
        device_array, host_array, buffer_array,
        solver_state_type
    )

    device_memory = arrays_byte_count
    host_memory = host_overhead + (device_is_host ? device_memory : arrays_byte_count)

    if params.use_MPI
        live_neighbours_factor = neighbour_count(params) / length(params.neighbours)
        if params.gpu_aware
            device_memory += Int(MPI_buffer_byte_count * live_neighbours_factor)
        else
            host_memory += Int(MPI_buffer_byte_count * live_neighbours_factor)
        end
    end

    return device_memory, host_memory
end


"""
    memory_required(N::NTuple{2, Int}, block_size::NTuple{2, Int}, ghost::Int, data_type)
    memory_required(N::NTuple{2, Int}, block_size::NTuple{2, Int}, ghost::Int,
        device_array_type, host_array_type, buffer_array_type[, solver_state_type])

Compute the number of bytes needed to allocate all blocks. If only `data_type` is given, then
`device_array_type`, `host_array_type` and `buffer_array_type` default to `Vector{T}`.
`solver_state_type` defaults to `SolverState{T, #= default schemes and test =#}`.

In order of returned values:
 1. Amount of bytes needed for all arrays on the device. This amount is also required on the host
    when the host and device are not the same.
 2. Amount of bytes needed for all MPI buffers, if the sub-domain has neighbours on all of its sides.
    If `params.gpu_aware`, then this memory is allocated on the device.
 3. Amount of bytes needed on the host memory for all block objects, excluding array data and buffers.
    This memory is always allocated on the host.

```
res = memory_required((1000, 1000), (64, 64), 4, CuArray{Float64}, Vector{Float64}, Vector{Float64})
device_memory = res[1]
host_memory = res[3] + (device_is_host ? device_memory : res[1])
if params.gpu_aware
    device_memory += res[2]
else
    host_memory += res[2]
end
```
"""
function memory_required(
    N::NTuple{2, Int}, block_size::NTuple{2, Int}, ghost::Int,
    ::Type{DeviceArray}, ::Type{HostArray}, ::Type{BufferArray}, ::Type{SState}
) where {
    T,
    DeviceArray <: AbstractArray{T},
    HostArray   <: AbstractArray{T},
    BufferArray <: AbstractArray{T},
    SState      <: SolverState
}
    grid_size, static_sized_grid, remainder_block_size = grid_dimensions(block_size, N, ghost)

    # Static blocks
    cell_count = prod(static_sized_grid) * prod(block_size)

    # Edge blocks
    for (block_count, _, edge_block_size) in EdgeBlockRegions(grid_size, static_sized_grid, block_size, remainder_block_size)
        cell_count += block_count * prod(edge_block_size)
    end

    arrays_byte_count = cell_count * length(block_vars()) * sizeof(T)

    # Size of `TaskBlock`s objects. They are only stored on the host memory.
    static_block_size = StaticBSize(block_size, ghost)
    static_block_count = prod(static_sized_grid)
    sizeof_static_block = sizeof(LocalTaskBlock{DeviceArray, HostArray, typeof(static_block_size), SState})
    blocks_overhead = static_block_count * sizeof_static_block

    edge_block_size = DynamicBSize(remainder_block_size, ghost)  # Dummy size
    edge_block_count = prod(grid_size) - static_block_count
    sizeof_edge_block = sizeof(LocalTaskBlock{DeviceArray, HostArray, typeof(edge_block_size), SState})
    blocks_overhead += edge_block_count * sizeof_edge_block

    remote_block_count = sum((length∘last).(RemoteBlockRegions(grid_size)))
    blocks_overhead += remote_block_count * sizeof(RemoteTaskBlock{BufferArray})

    # MPI Buffers size
    domain_perimeter = sum(N) * 2  # TODO: dimension agnostic
    MPI_buffer_byte_count = domain_perimeter * #= send+recv =# 2 * length(comm_vars()) * ghost * sizeof(T)

    return arrays_byte_count, MPI_buffer_byte_count, blocks_overhead
end

memory_required(N::Tuple, block_size::Tuple, ghost::Int, device_array, host_array, buffer_array) =
    memory_required(N, block_size, ghost, device_array, host_array, buffer_array,
        # Default `SolverState` for a good enough estimation
        SolverState{T, GodunovSplitting, RiemannGodunov, MinmodLimiter, EulerProjection, Sod})

memory_required(N::Tuple, block_size::Tuple, ghost::Int, ::Type{T}) where {T} =
    memory_required(N, block_size, ghost, Vector{T}, Vector{T}, Vector{T})


"""
    device_to_host!(grid::BlockGrid)

Copies device data of all blocks to the host data. A no-op if the device is the host.
"""
function device_to_host!(grid::BlockGrid{<:Any, D, H}) where {D, H}
    for blk in all_blocks(grid)
        device_to_host!(blk)
    end
end

device_to_host!(::BlockGrid{<:Any, D, D}) where {D} = nothing


"""
    device_to_host!(grid::BlockGrid)

Copies host data of all blocks to the device data. A no-op if the device is the host.
"""
function host_to_device!(grid::BlockGrid{<:Any, D, H}) where {D, H}
    for blk in all_blocks(grid)
        host_to_device!(blk)
    end
end

host_to_device!(::BlockGrid{<:Any, D, D}) where {D} = nothing


"""
    move_pages(grid::BlockGrid)

Move the pages of all blocks of the `grid`, including remote blocks, to the NUMA node of the thread
which is in charge of working on that block.
"""
function move_pages(grid::BlockGrid)
    numa_map = tid_to_numa_node_map()
    for (tid, blks_pos) in enumerate(grid.threads_workload), blk_pos in blks_pos
        target_numa = numa_map[tid]
        blk = block_at(grid, blk_pos)
        move_pages(blk, target_numa)

        # Make sure the MPI buffers are as close as the data they will be interacting with
        for neighbour in blk.neighbours
            !(neighbour isa RemoteTaskBlock) && continue
            move_pages(neighbour, target_numa)
        end
    end
end


"""
    lock_pages(grid::BlockGrid)

Locks the pages of all blocks of the `grid`, including remote blocks.
"""
function lock_pages(grid::BlockGrid)
    foreach(lock_pages, grid.blocks)
    foreach(lock_pages, grid.edge_blocks)
    foreach(lock_pages, grid.remote_blocks)
end


"""
    block_pos_containing_cell(grid::BlockGrid, pos::Union{CartesianIndex, NTuple})

Returns two `CartesianIndex`es: the first is the position of the block containing the cell at `pos`,
the second is the position of the cell in that block.
"""
function block_pos_containing_cell(grid::BlockGrid, pos)
    if !in_grid(pos, grid.cell_size)
        error("real cell position $pos is outside of the grid: $(grid.cell_size)")
    elseif prod(static_block_size(grid)) == 0
        # No blocking: there is only a single block in the grid
        return CartesianIndex{length(pos)}(1), CartesianIndex(pos)
    end

    block_pos = (Tuple(pos) .- 1) .÷ real_block_size(grid) .+ 1
    if !in_grid(block_pos, grid.static_sized_grid)
        # In the edge blocks
        block_pos = clamp.(block_pos, 1, grid.grid_size)
    end

    cell_pos = Tuple(pos) .- (block_pos .- 1) .* real_block_size(grid)
    return CartesianIndex(block_pos), CartesianIndex(cell_pos)
end


"""
    block_origin(grid::BlockGrid, pos, include_ghosts=false)

A `Tuple` of the position of the cell at the bottom left corner of the [`LocalTaskBlock`](@ref) at
`pos` in the `grid`.

If `include_ghosts == true`, then the cell position includes all ghost cells of the `grid`.

`pos` can be any of `Integer`, `NTuple{N, Integer}` or `CartesianIndex`.
"""
function block_origin(grid::BlockGrid, pos, include_ghosts=false)
    bs = include_ghosts ? static_block_size(grid) : real_block_size(grid)
    if in_grid(pos, grid.static_sized_grid)
        return bs .* (Tuple(pos) .- 1) .+ 1
    elseif any(grid.static_sized_grid .== 0)
        return ifelse.(grid.grid_size .> 1, bs .* (Tuple(pos) .- 1) .+ 1, 1)
    else
        return ifelse.(
            in_grid.(Ref(Tuple(pos)), Ref(grid.static_sized_grid), instances(Axis.T)),
            bs .* (Tuple(pos) .- 1) .+ 1,
            bs .* grid.static_sized_grid .+ 1
        )
    end
end


function print_grid_dimensions(
    io::IO, grid_size::Tuple, static_grid::Tuple, static_block_size::Tuple,
    cell_size::Tuple, ghost; pad=20
)
    grid_str = join(grid_size, '×')
    block_count = prod(grid_size)

    static_grid_str = join(static_grid, '×')
    static_block_str = join(static_block_size, '×')
    static_block_count = prod(static_grid)
    static_cell_count = prod(static_block_size)

    real_size = join(static_block_size .- 2*ghost, '×')
    real_count = prod(static_block_size .- 2*ghost)
    ghost_count = static_cell_count - real_count
    real_ghost_ratio = @sprintf("%.03g%%", ghost_count / static_cell_count * 100)

    edge_block_count = prod(grid_size) - static_block_count
    edge_pos = String[]
    for axis in instances(Axis.T)
        static_grid[Int(axis)] ≥ grid_size[Int(axis)] && continue
        push!(edge_pos, lowercasefirst(string(last_side(axis))))
    end
    if !isempty(edge_pos)
        edge_pos_str = "at the " * join(edge_pos, ", ", " and ") * " edge"
        length(edge_pos) > 1 && (edge_pos_str *= "s")
    else
        edge_pos_str = "none"
    end

    static_block_cells = static_block_count * prod(static_block_size .- 2*ghost)
    edge_block_cells = prod(cell_size) - static_block_cells
    edge_cell_ratio = @sprintf("%.03g%%", edge_block_cells / prod(cell_size) * 100)
    edge_block_ratio = @sprintf("%.03g%%", edge_block_count / block_count * 100)

    remote_block_count = 2*sum(grid_size)
    remote_buffers_size = 2*sum(cell_size) * ghost

    static_block_ratio = @sprintf("%.03g%%", static_block_count / block_count * 100)

    print_parameter(io, pad, "block size", "$static_block_str cells ($static_cell_count total)")

    print_parameter(io, pad, "grid", "$grid_str blocks ($block_count total)")
    if static_block_count == 0 && edge_block_count == 1
        # No blocking: there is only a single block in the grid
        total_cells = prod(cell_size .+ 2*ghost)
        ghost_count = total_cells - prod(cell_size)
        real_ghost_ratio = @sprintf("%.03g%%", ghost_count / total_cells * 100)
        print_parameter(io, pad, "static grid", "0 static blocks")
        print_parameter(io, pad, "edge grid", "1 edge block, containing all cells, \
            with $ghost ghost cells ($ghost_count total, $real_ghost_ratio of the block)")
    else
        print_parameter(io, pad, "static grid", "$static_grid_str static blocks \
            ($static_block_count total, $static_block_ratio of all blocks)")
        print_parameter(io, pad, "static block",
            "$real_size real cells ($real_count total), \
            with $ghost ghost cells ($ghost_count total, $real_ghost_ratio of the block)")
        print_parameter(io, pad, "edge grid", "$edge_block_count edge blocks ($edge_block_ratio of all blocks)")
        print_parameter(io, pad, "edge blocks", "$edge_pos_str, containing $edge_cell_ratio of all real cells")
    end
    print_parameter(io, pad, "remote grid", "$remote_block_count remote blocks, \
        containing up to $remote_buffers_size cells"; nl=false)
end


function Base.show(io::IO, ::MIME"text/plain", grid::BlockGrid{T, D, H, B, Ghost, BS, Device};
    pad=16
) where {T, D, H, B, Ghost, BS, Device}
    println(io, "BlockGrid:")
    print_grid_dimensions(io, grid.grid_size, grid.static_sized_grid, block_size(BS), grid.cell_size, Ghost; pad)
    println()

    remote_dev_str = buffers_on_device(grid) ? "device" : "host"
    print_parameter(io, pad, "remote buffers", "stored on the $remote_dev_str")
    print_parameter(io, pad, "device", grid.device)
    print_parameter(io, pad, "device array", D)
    print_parameter(io, pad, "host array", D == H ? "same as device" : H; nl=false)
end


function Base.show(io::IO, grid::BlockGrid{T, D, H, B, G, BS}) where {T, D, H, B, G, BS}
    grid_str = join(grid.grid_size, '×')
    bs_str = join(block_size(BS), '×')
    print(io, "BlockGrid{$T, $D, $H, $B}($grid_str, bs: $bs_str, ghost: $G)")
end
