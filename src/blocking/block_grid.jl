
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
    cell_size          :: NTuple{2, Int}  # Number of cells in each direction
    edge_size          :: NTuple{2, Int}  # Number of cells in edge blocks in each direction (only along non-edge directions)
    device             :: Device
    global_dt          :: GlobalTimeStep{T}
    blocks             :: Vector{LocalTaskBlock{DeviceArray, HostArray, BS, SState}}
    edge_blocks        :: Vector{LocalTaskBlock{DeviceArray, HostArray, DynamicBSize{Ghost}, SState}}
    remote_blocks      :: Vector{RemoteTaskBlock{BufferArray}}
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

    # Main grid container
    grid = BlockGrid{
        T, device_array, host_array, buffer_array,
        ghost, typeof(static_size),
        state_type, typeof(params.device)
    }(
        grid_size, static_sized_grid, cell_size, remainder_block_size, params.device, global_dt,
        blocks, edge_blocks, remote_blocks
    )

    # Allocate all local blocks
    # Non-static blocks are placed on the right and top sides.
    inner_grid = CartesianIndices(static_sized_grid)
    device_kwargs = alloc_device_kwargs(params)
    host_kwargs = alloc_host_kwargs(params)
    for pos in CartesianIndices(grid_size)
        if pos in inner_grid
            blks = blocks
            idx = block_idx(grid, pos)
            blk_size = static_size
        else
            blks = edge_blocks
            idx = edge_block_idx(grid, pos)
            blk_size = DynamicBSize(ifelse.(Tuple(pos) .== grid_size, remainder_block_size, params.block_size), ghost)
        end

        blk_state = SolverState(params, global_dt)
        blks[idx] = eltype(blks)(blk_size, pos, blk_state, device_kwargs, host_kwargs)
    end

    # Allocate all remote blocks
    # X (Left, Right) then Y (Bottom, Top)
    remote_i = 1
    for axis in (Axis.X, Axis.Y), side in sides_along(axis)
        # Range of all blocks along `side`
        blocks_x = axis == Axis.X ? (side == Side.Left   ? (1:1) : (grid_size[1]:grid_size[1])) : (1:grid_size[1])
        blocks_y = axis == Axis.Y ? (side == Side.Bottom ? (1:1) : (grid_size[2]:grid_size[2])) : (1:grid_size[2])

        for our_pos in CartesianIndices((blocks_x, blocks_y))
            pos = our_pos + CartesianIndex(offset_to(side))

            if has_neighbour(params, side)
                # The buffer must be the same size as the side of our block which is a neighbour to...
                buffer_size = size_along(block_at(grid, our_pos).size, side)
                # ...for each variable to communicate of each ghost cell
                buffer_size *= length(comm_vars()) * params.nghost

                neighbour = neighbour_at(params, side)  # rank
                global_pos = CartesianIndex(params.cart_coords .+ offset_to(side))  # pos in the cart_comm

                block = RemoteTaskBlock{buffer_array}(buffer_size, pos, neighbour, global_pos, params.cart_comm)
            else
                # "Fake" remote block for non-existant neighbour at the edge of the global domain
                block = RemoteTaskBlock{buffer_array}(pos)
            end

            remote_blocks[remote_i] = block
            remote_i += 1
        end
    end

    # Initialize all block neighbours references
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

        for blk in this_block.neighbours
            blk isa RemoteTaskBlock || continue
            blk.neighbour = this_block
        end
    end

    return grid
end


"""
    grid_dimensions(params::ArmonParameters)
    grid_dimensions(block_size::NTuple{D, Int}, domain_size::NTuple{D, Int}, ghost::Int) where {D}

Returns the dimensions of the grid in the form `(grid_size, static_sized_grid, remainder_block_size)`
from the `block_size` (the size of blocks in the `static_sized_grid`), the `domain_size` (number of
real cells) and the number of `ghost` cells, common to all blocks.

`grid_size` is the `static_sized_grid` including the edge blocks.
Edge blocks along the axis `d` have a size of `remainder_block_size[d]` along `d`, and `block_size`
for the other axes.

`block_size` includes the `ghost` cells in its dimensions, which must all be greater than `2*ghost`.

If `prod(block_size) == 0`, then `block_size` is ignored and the grid is made of only a single block
of size `domain_size`.
"""
grid_dimensions(params::ArmonParameters) = grid_dimensions(params.block_size, params.N, params.nghost)

function grid_dimensions(block_size::NTuple{D, Int}, domain_size::NTuple{D, Int}, ghost::Int) where {D}
    if prod(block_size) == 0
        # No cache blocking: only one dynamic block for the whole grid
        return ntuple(Returns(1), D), ntuple(Returns(0), D), domain_size .+ 2*ghost
    end

    # `block_size` includes the number of ghost cells, while `domain_size` is in real cells.
    real_block_size = block_size .- 2*ghost  # number of real cells in a block in each dim
    if prod(real_block_size) < 1
        solver_error(:config, "block size $block_size is too small with $ghost ghost cells: $real_block_size")
    end

    grid_size = domain_size .÷ real_block_size
    remainder = domain_size .% real_block_size

    prod(grid_size) == 0 && (grid_size = ntuple(Returns(0), D))

    # Edge-case: the remainder is smaller than the number of ghost cells. This creates problems during
    # the halo-exchange since ghosts cells at the edge of the static sized grid depend on all real cells
    # of the neighbouring edge block BUT ALSO the ghosts of that same edge block, themselves depending
    # on a remote block (or on the global border conditions).
    # To solve this, those remaining cells are merged with the neighbouring static sized blocks.
    remainder_too_small = 0 .< remainder .< ghost .&& grid_size .> 0
    if any(remainder_too_small)
        remainder = remainder .+ real_block_size .* remainder_too_small
        grid_size = grid_size .- 1 .* remainder_too_small
    end

    static_sized_grid = grid_size  # Grid of blocks with a static size
    prod(remainder) > 0 && (grid_size = grid_size .+ (remainder .> 0))

    if prod(grid_size) < 1
        solver_error(:config, "could not partition $domain_size domain into $block_size blocks with $ghost ghost cells")
    end

    # (x, y) size of dynamic blocks, including ghost cells. A block for a X axis edge will have a
    # size of (x, block_size.y), and for the Y axis edge the size would be (block_size.x, y). In the
    # edge corner the block would be of size (x, y).
    if prod(remainder) == 0
        remainder_block_size = ntuple(Returns(0), D)
    elseif prod(static_sized_grid) > 0
        remainder_block_size = remainder .+ 2 * ghost
    else
        # Edge-case for when the block size is bigger than the domain
        remainder_block_size = domain_size .+ 2 * ghost
    end

    return grid_size, static_sized_grid, remainder_block_size
end


"Linear index of a block in the statically-sized grid"
block_idx(grid::BlockGrid, idx::CartesianIndex{2}) =
    LinearIndices(CartesianIndices(grid.static_sized_grid))[idx]


"Linear index of a block at the (dynamically-sized) edges of the grid"
function edge_block_idx(grid::BlockGrid, idx::CartesianIndex{2})
    # A 2D grid has 2 edge regions:
    #  - Last X column (including last block of last Y row if they overlap)
    #  - Last Y row
    # They are stored contigously.
    i = 0
    offset = 0

    if grid.static_sized_grid[1] < grid.grid_size[1]
        if grid.static_sized_grid[1] < idx[1]
            offset = idx[2]
        else
            offset = grid.static_sized_grid[2] + 1
        end
    end

    if grid.static_sized_grid[2] < grid.grid_size[2]
        if grid.static_sized_grid[2] < idx[2]
            i = idx[1] + offset - 1  # -1 to account for the overlap
        else
            i = offset
        end
    else
        i = offset
    end

    i == 0 && error("Index $idx is not at the edge of the grid $grid")
    return i
end


"Linear index of a remote block at the edges of the grid"
function remote_block_idx(grid::BlockGrid, idx::CartesianIndex{2})
    offset = 0
    for axis in instances(Axis.T), side in sides_along(axis)
        ai = Int(axis)  # Axis index

        side_length = grid.grid_size[ai]
        edge_length = prod(grid.grid_size) ÷ side_length

        side_idx = side in first_sides() ? 1 : side_length
        side_idx += offset_to(side)[ai]

        edge_idx = Tuple(idx)[mod1(ai + 1, 2)]  # TODO: not dimension-agnostic :(

        if Tuple(idx)[ai] == side_idx && 1 ≤ edge_idx ≤ edge_length
            return offset + edge_idx
        else
            offset += edge_length
        end
    end

    throw(ArgumentError("Block index $idx is not on the remote edges of the grid"))
end


function block_at(grid::BlockGrid, idx::CartesianIndex{2})
    if in_grid(idx, grid.static_sized_grid)
        return grid.blocks[block_idx(grid, idx)]
    elseif in_grid(idx, grid.grid_size)
        return grid.edge_blocks[edge_block_idx(grid, idx)]
    else
        return grid.remote_blocks[remote_block_idx(grid, idx)]
    end
end


device_array_type(::ObjOrType{BlockGrid{<:Any, D}}) where {D} = D
host_array_type(::ObjOrType{BlockGrid{<:Any, <:Any, H}}) where {H} = H
buffer_array_type(::ObjOrType{BlockGrid{<:Any, <:Any, <:Any, B}}) where {B} = B
ghosts(::ObjOrType{BlockGrid{<:Any, <:Any, <:Any, <:Any, Ghost}}) where {Ghost} = Ghost
static_block_size(::ObjOrType{BlockGrid{<:Any, <:Any, <:Any, <:Any, G, BS}}) where {G, BS} = BS


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

    # Edge blocks on the right, along the Y axis (excluding the top-right corner)
    if static_sized_grid[1] < grid_size[1]
        cell_count += static_sized_grid[2] * (remainder_block_size[1] * block_size[2])
    end

    # Edge blocks on the top, along the X axis (excluding the top-right corner)
    if static_sized_grid[2] < grid_size[2]
        cell_count += static_sized_grid[1] * (remainder_block_size[2] * block_size[1])
    end

    # Edge block on the top-right corner
    if static_sized_grid[1] < grid_size[1] && static_sized_grid[2] < grid_size[2]
        cell_count += prod(remainder_block_size)
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

    remote_block_count = sum(grid_size) * 2
    blocks_overhead += remote_block_count * sizeof(RemoteTaskBlock{BufferArray})

    # MPI Buffers size
    domain_perimeter = sum(N) * 2
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
    real_ghost_ratio = @sprintf("%.02g%%", ghost_count / static_cell_count * 100)

    edge_block_count = prod(grid_size) - static_block_count
    edge_pos = []
    static_grid[1] < grid_size[1] && push!(edge_pos, "right")
    static_grid[2] < grid_size[2] && push!(edge_pos, "top")
    if !isempty(edge_pos)
        edge_pos_str = "at the " * join(edge_pos, ", ", " and ") * " edge"
        length(edge_pos) > 1 && (edge_pos_str *= "s")
    else
        edge_pos_str = "none"
    end

    static_block_cells = static_block_count * prod(static_block_size .- 2*ghost)
    edge_block_cells = prod(cell_size) - static_block_cells
    edge_cell_ratio = @sprintf("%.02g%%", edge_block_cells / prod(cell_size) * 100)
    edge_block_ratio = @sprintf("%.02g%%", edge_block_count / block_count * 100)

    remote_block_count = 2*sum(grid_size)
    remote_buffers_size = 2*sum(cell_size) * ghost

    static_block_ratio = @sprintf("%.02g%%", static_block_count / block_count * 100)

    print_parameter(io, pad, "block size", "$static_block_str cells ($static_cell_count total)")

    print_parameter(io, pad, "grid", "$grid_str blocks ($block_count total)")
    if static_block_count == 0 && edge_block_count == 1
        # No blocking: there is only a single block in the grid
        total_cells = prod(cell_size .+ 2*ghost)
        ghost_count = total_cells - prod(cell_size)
        real_ghost_ratio = @sprintf("%.02g%%", ghost_count / total_cells * 100)
        print_parameter(io, pad, "static grid", "0 static blocks")
        print_parameter(io, pad, "edge grid", "1 edge block, containing all cells, \
            with $ghost ghost cells ($ghost_count total, $real_ghost_ratio of the block)")
    else
        print_parameter(io, pad, "static grid", "$static_grid_str static blocks \
            ($static_block_count total, $static_block_ratio)")
        print_parameter(io, pad, "static block",
            "$real_size real cells ($real_count total), \
            with $ghost ghost cells ($ghost_count total, $real_ghost_ratio of the block)")
        print_parameter(io, pad, "edge grid", "$edge_block_count edge blocks ($edge_block_ratio)")
        print_parameter(io, pad, "edge blocks", "$edge_pos_str, containing $edge_cell_ratio of all real cells")
    end
    print_parameter(io, pad, "remote grid", "$remote_block_count remote blocks, \
        containing $remote_buffers_size cells (total)"; nl=false)
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
