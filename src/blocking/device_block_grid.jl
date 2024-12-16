
"""
    DeviceLocalBlock

Mirror of a [`LocalTaskBlock`](@ref) on the device (GPU).

The structure can be entirely stored and manipulated from the device.
"""
struct DeviceLocalBlock{D <: AbstractArray, Size <: BlockSize}
    size            :: Size
    pos             :: CartesianIndex{2}          # position in the grid
    interfaces_idx  :: NTuple{2, NTuple{2, Int}}  # indexes of the block's interfaces (per axis, then per side)
    base_status_idx :: Int                        # index of the status of the first side of the block
    data            :: BlockData{D}
end

Adapt.@adapt_structure DeviceLocalBlock


function DeviceLocalBlock(block::LocalTaskBlock{D, <:Any, Size}, grid_size, device) where {D, Size}
    base_status_idx = base_block_interface_status_index(grid_size, block.pos)
    interfaces_idx = ntuple(2) do i
        axis = instances(Axis.T)[i]
        return (
            block_interface_index(grid_size, block.pos, first_side(axis)),
            block_interface_index(grid_size, block.pos, last_side(axis)),
        )
    end

    # Data arrays to pointers, as `isbitstype(Array) == false`, which is unsupported on GPU.
    # Since arrays are owned by the host `block`, they cannot be GC'ed.
    kernel_data = Adapt.adapt(device_converter(device), block.device_data)

    return DeviceLocalBlock{D, Size}(
        block.size, block.pos, interfaces_idx, base_status_idx, kernel_data
    )
end


"""
    DeviceRemoteBlock

Mirror of a [`RemoteTaskBlock`](@ref) on the device (GPU).

The communication logic isn't present on the device, as it is intrinsically the host's job: MPI or
NCCL (and others) communication are all initiated from a host call.
"""
struct DeviceRemoteBlock{B <: AbstractArray}
    pos       :: CartesianIndex{2}
    exists    :: Bool  # `false` if there is no remote block, and it a border of a global domain
    on_device :: Bool
    buffer    :: B     # empty if `(on_device || exists) == false`
end

Adapt.@adapt_structure DeviceRemoteBlock


function DeviceRemoteBlock(block::RemoteTaskBlock, device)
    if block.rank != -1 && block.on_device
        # TODO: which buffer to take? how about send/recv buffers? what about communication styles with global buffers?
        error("NYI")
        buffer = 0
    else
        # Empty dummy buffer
        buffer = device_array_type(device){eltype(block), 1}()
    end
    kernel_buffer = Adapt.adapt(device_converter(device), buffer)
    return DeviceRemoteBlock(block.pos, block.rank != -1, block.on_device, kernel_buffer)
end


"""
    DeviceBlockGrid

Mirror of a [`BlockGrid`](@ref) such that it can be stored and used from the device (e.g. a GPU).

All data arrays are replaced with device pointers: they are immutable and therefore GPU-compatible.

`@atomic` fields of block interfaces, which require a mutable container on the host, are instead all
stored in the same array, which itself is indexed using atomic operations.

It is necessarily backed by a host [`BlockGrid`](@ref), as device array pointers do not own their
data from the perspective of the GC: it is the host grid which would prevent them from being
discarded.
"""
struct DeviceBlockGrid{
    T,
    DeviceArray <: AbstractArray{T},
    Device,
    Ghost,
    BS <: StaticBSize{<:Any, Ghost},
    IndexMap <: AbstractArray{UInt32, 2},
    Interfaces <: GridInterfaces,
    SB_Container <: AbstractArray{DeviceLocalBlock{DeviceArray, BS}},
    EB_Container <: AbstractArray{DeviceLocalBlock{DeviceArray, DynamicBSize{Ghost}}},
    RB_Container <: AbstractArray{DeviceRemoteBlock{DeviceArray}}
} <: AbstractBlockGrid{T, Ghost, BS, Device}
    # Same fields as for `BlockGrid`
    grid_size         :: NTuple{2, Int}
    static_sized_grid :: NTuple{2, Int}
    cell_size         :: NTuple{2, Int}
    edge_size         :: NTuple{2, Int}

    device            :: Device

    # 2D array from block position to its linear index in its respective container
    # Remote blocks are included, therefore it should be indexed with an offset of 1 for the 2D index
    index_map         :: IndexMap  # TODO: use the 2-3 LSB of an index to identify the block kind, leaving >1e8 blocks per kind which is plenty

    # Interfaces between the blocks of the grid
    interfaces        :: Interfaces

    # Device arrays of device blocks
    blocks            :: SB_Container
    edge_blocks       :: EB_Container
    remote_blocks     :: RB_Container
end

Adapt.@adapt_structure DeviceBlockGrid


device(grid::DeviceBlockGrid) = grid.device
grid_sizes(grid::DeviceBlockGrid) =
    (; grid=grid.grid_size, static_grid=grid.static_sized_grid, real_cells=grid.cell_size, edge=grid.edge_size)


function DeviceBlockGrid(::Type{T}, device::Device, static_size::StaticBSize, dyn_size_t::Type{<:DynamicBSize}, grid_sizes) where {T, Device}
    device_array = device_array_type(device)

    SB_Container = device_array{DeviceLocalBlock{device_array, typeof(static_size)}, 1}
    EB_Container = device_array{DeviceLocalBlock{device_array, dyn_size_t}, 1}
    RB_Container = device_array{DeviceRemoteBlock{device_array}, 1}

    IndexMap = device_array{UInt32, 2}
    index_map = IndexMap(undef, grid_sizes.grid .+ 1)

    interfaces = GridInterfaces(grid_sizes.grid, device)

    dev_blocks = SB_Container(undef, size(host_grid.blocks))
    dev_edge_blocks = EB_Container(undef, size(host_grid.edge_blocks))
    dev_remote_blocks = RB_Container(undef, size(host_grid.remote_blocks))

    return DeviceBlockGrid{
        T, D, Device, ghosts(static_size), typeof(static_size), IndexMap, typeof(interfaces),
        SB_Container, EB_Container, RB_Container
    }(
        grid_sizes.grid, grid_sizes.static_grid, grid_sizes.real_cells, grid_sizes.edge,
        device, index_map, interfaces,
        dev_blocks, dev_edge_blocks, dev_remote_blocks
    )
end


function init_device_block_grid!(device_grid::DeviceBlockGrid, host_grid::AbstractBlockGrid)
    # Working on device memory without a kernel is annoying: we do everything on the host then copy
    # to the device. The alternative of initializing directly on the device is not possible as most
    # host data structures are mutable or `!isbitstype`, and therefore cannot be sent on the device.
    host_index_map = zeros(UInt32, size(device_grid.index_map))
    host_blocks = Vector{eltype(device_grid.blocks)}(undef, size(device_grid.blocks))
    host_edge_blocks = Vector{eltype(device_grid.edge_blocks)}(undef, size(device_grid.edge_blocks))
    host_remote_blocks = Vector{eltype(device_grid.remote_blocks)}(undef, size(device_grid.remote_blocks))

    # Local blocks
    for (idx, blk) in enumerate(host_grid.blocks)
        host_index_map[blk.pos + one(CartesianIndex{2})] = idx
        host_blocks[idx] = DeviceLocalBlock(blk, sizes.grid, dev)
    end

    for (idx, blk) in enumerate(host_grid.edge_blocks)
        host_index_map[blk.pos + one(CartesianIndex{2})] = idx
        host_edge_blocks[idx] = DeviceLocalBlock(blk, sizes.grid, dev)
    end

    # Remote blocks
    for (idx, blk) in enumerate(host_grid.remote_blocks)
        host_index_map[blk.pos + one(CartesianIndex{2})] = idx
        host_edge_blocks[idx] = DeviceRemoteBlock(blk, dev)
    end

    copyto!(device_grid.index_map, host_index_map)
    copyto!(device_grid.blocks, host_blocks)
    copyto!(device_grid.edge_blocks, host_edge_blocks)
    copyto!(device_grid.remote_blocks, host_remote_blocks)

    return device_grid
end


function put_block_grid_on_device(grid::DeviceBlockGrid)
    # `kernel_grid` is a device-code compatible version of `grid`, itself a copy of
    # a `BlockGrid` on the device's memory
    kernel_grid = Adapt.adapt(device_converter(grid.device), grid)

    # `sizeof(kernel_grid)` is large, we want to minimize kernel launch overhead by storing it once
    # in global memory. This way we only have to pass a pointer with `pointer(kernel_grid_arr)` when
    # launching the kernel, which can then retrieve the `kernel_grid` with `unsafe_load`.
    kernel_grid_arr = device_array_type(grid.device){typeof(kernel_grid), 0}(undef)
    kernel_grid_arr[1:1] .= Ref(kernel_grid)  # equiv. to `kernel_grid_arr[1] = kernel_grid` but without the `allowscalar` checks

    return kernel_grid_arr
end


function block_kind(grid::DeviceBlockGrid, pos::CartesianIndex)
    if in_grid(pos, grid.static_sized_grid)
        return :static
    elseif in_grid(pos, grid.grid_size)
        return :edge
    else
        return :remote
    end
end


"""
    @apply_device_block blk = grid[pos] begin
        # body
    end

    @apply_device_block blk = pointer(grid, pos) begin
        # body
    end

Indexes `grid` at `pos` to retrive `blk`, then applies `body`.
`body` is duplicated for each kind of block (static or edge): therefore everything is
type-safe and no runtime dispatch is performed in `body`.

Indexing a remote block will result in an runtime error.

Using `pointer(grid, pos)` in the first macro argument will turn `blk` into a pointer to the block
in the grid.
"""
macro apply_device_block(grid_get_expr, expr)
    if @capture(grid_get_expr, blk_ = grid_[pos_])
        get_block = :(@inbounds $grid.blocks[block_lin_idx])
        get_edge_block = :(@inbounds $grid.edge_blocks[block_lin_idx])

    elseif @capture(grid_get_expr, blk_ = pointer(grid_, pos_))
        get_block = :(pointer($grid.blocks, block_lin_idx))
        get_edge_block = :(pointer($grid.edge_blocks, block_lin_idx))

    else
        error("Expected a `blk = grid[pos]` expression, got: ", grid_get_expr)
    end

    return esc(quote
        kind = block_kind($grid, $pos)
        block_lin_idx = $grid.index_map[$pos + one(CartesianIndex{2})]
        if kind === :static
            $blk = $get_block
            $expr
        elseif kind === :edge
            $blk = $get_edge_block
            $expr
        else kind === :remote
            error("unexpected index to a remote block")
        end
    end)
end

