
"""
    TaskBlock{V}

Abstract block used for cache blocking.
"""
abstract type TaskBlock{V <: AbstractArray} end

array_type(::TaskBlock{V}) where {V} = V
Base.eltype(::TaskBlock{V}) where {V} = eltype(V)


"""
    BlockData{V}

Holds the variables of all cells of a [`LocalTaskBlock`](@ref).
"""
struct BlockData{V}
    x      :: V
    y      :: V
    ρ      :: V
    u      :: V
    v      :: V
    E      :: V
    p      :: V
    c      :: V
    g      :: V
    uˢ     :: V
    pˢ     :: V
    work_1 :: V
    work_2 :: V
    work_3 :: V
    work_4 :: V
    mask   :: V  # TODO: remove ??

    function BlockData{V}(size; kwargs...) where {V}
        vars = ntuple(length(block_vars())) do i
            var = block_vars()[i]
            label = string(var)
            return V(undef, size; alloc_array_kwargs(; label, kwargs...)...)
        end
        return new{V}(vars...)
    end
end


block_vars() = (:x, :y, :ρ, :u, :v, :E, :p, :c, :g, :uˢ, :pˢ, :work_1, :work_2, :work_3, :work_4, :mask)
main_vars()  = (:x, :y, :ρ, :u, :v, :E, :p, :c, :g, :uˢ, :pˢ)  # Variables synchronized between host and device
saved_vars() = (:x, :y, :ρ, :u, :v,     :p                  )  # Variables saved to/read from I/O
comm_vars()  = (        :ρ, :u, :v, :E, :p, :c, :g          )  # Variables exchanged between ghost cells

get_vars(data::BlockData, vars) = getfield.(Ref(data), vars)
block_vars(data::BlockData) = get_vars(data, block_vars())
main_vars(data::BlockData)  = get_vars(data, main_vars())
saved_vars(data::BlockData) = get_vars(data, saved_vars())
comm_vars(data::BlockData)  = get_vars(data, comm_vars())


"""
    LocalTaskBlock{D, H, Size, SolverState} <: TaskBlock{V}

Container of `Size` and variables of type `D` on the device and `H` on the host.
Part of a [`BlockGrid`](@ref).

The block stores its own solver state, allowing it to run all solver steps independantly of all
other blocks, apart from steps requiring synchronization.
"""
mutable struct LocalTaskBlock{D, H, Size <: BlockSize, SState <: SolverState} <: TaskBlock{D}
    state        :: SState             # Solver state for the block
    exchange_age :: Atomic{Int}        # Incremented every time an exchange is completed
    exchanges    :: Neighbours{Atomic{BlockExchangeState.T}}  # State of ghost cells exchanges for each side
    size         :: Size               # Size (in cells) of the block
    pos          :: CartesianIndex{2}  # Position in the local block grid
    neighbours   :: Neighbours{TaskBlock}
    device_data  :: BlockData{D}
    host_data    :: BlockData{H}  # Host data uses the same arrays as device data if `D == H`
    # TODO: device? storing the associated GPU stream, or CPU cores (or maybe not, to allow relocations?)

    function LocalTaskBlock{D, H, Size, SState}(size::Size, pos, blk_state::SState, device_kwargs, host_kwargs) where {D, H, Size, SState}
        # `neighbours` is set afterwards, when all blocks are created.
        block = new{D, H, Size, SState}(
            blk_state, Atomic(0),
            Neighbours(Atomic{BlockExchangeState.T}, BlockExchangeState.NotReady),
            size, pos #= undef =#
        )

        cell_count = prod(block_size(size))
        if D != H
            block.device_data = BlockData{D}(cell_count; device_kwargs...)
            block.host_data   = BlockData{H}(cell_count; host_kwargs...)
        else
            block.host_data = block.device_data = BlockData{D}(cell_count; device_kwargs...)
        end

        return block
    end
end


block_size(blk::LocalTaskBlock) = block_size(blk.size)
real_block_size(blk::LocalTaskBlock) = real_block_size(blk.size)
ghosts(blk::LocalTaskBlock) = ghosts(blk.size)
solver_state(blk::LocalTaskBlock) = blk.state

block_device_data(blk::LocalTaskBlock) = blk.device_data
block_host_data(blk::LocalTaskBlock) = blk.host_data
block_data(blk::LocalTaskBlock; on_device=true) = on_device ? block_device_data(blk) : block_host_data(blk)

get_vars(blk::LocalTaskBlock, vars; on_device=true) = get_vars(block_data(blk; on_device), vars)
block_vars(blk::LocalTaskBlock; on_device=true) = block_vars(block_data(blk; on_device))
main_vars(blk::LocalTaskBlock; on_device=true)  = main_vars(block_data(blk; on_device))
saved_vars(blk::LocalTaskBlock; on_device=true) = saved_vars(block_data(blk; on_device))
comm_vars(blk::LocalTaskBlock; on_device=true)  = comm_vars(block_data(blk; on_device))


exchange_age(blk::LocalTaskBlock) = @atomic blk.exchange_age.x
incr_exchange_age!(blk::LocalTaskBlock) = @atomic blk.exchange_age.x += 1

exchange_state(blk::LocalTaskBlock, side::Side.T) = @atomic blk.exchanges[Int(side)].x
exchange_state!(blk::LocalTaskBlock, side::Side.T, state::BlockExchangeState.T) = @atomic blk.exchanges[Int(side)].x = state
function replace_exchange_state!(blk::LocalTaskBlock, side::Side.T, transition::Pair{BlockExchangeState.T, BlockExchangeState.T})
    _, ok = @atomicreplace blk.exchanges[Int(side)].x transition
    return ok
end


function reset!(blk::LocalTaskBlock)
    reset!(blk.state)
    @atomic blk.exchange_age.x = 0
    for exchange_state in blk.exchanges
        @atomic exchange_state.x = BlockExchangeState.NotReady
    end
end


"""
    device_to_host!(blk::LocalTaskBlock)

Copies the device data of `blk` to the host data. A no-op if the device is the host.
"""
function device_to_host!(blk::LocalTaskBlock{D, H}) where {D, H}
    for (dst_var, src_var) in zip(main_vars(blk.host_data), main_vars(blk.device_data))
        # TODO: ensure this is done in the correct thread/device
        # KernelAbstractions.copyto!(blk.device, dst_var, src_var)
        copyto!(dst_var, src_var)
    end
end

device_to_host!(::LocalTaskBlock{H, H}) where {H} = nothing


"""
    device_to_host!(blk::LocalTaskBlock)

Copies the host data of `blk` to its device data. A no-op if the device is the host.
"""
function host_to_device!(blk::LocalTaskBlock{D, H}) where {D, H}
    for (dst_var, src_var) in zip(main_vars(blk.device_data), main_vars(blk.host_data))
        # TODO: ensure this is done in the correct thread/device
        # KernelAbstractions.copyto!(blk.device, dst_var, src_var)
        copyto!(dst_var, src_var)
    end
end

host_to_device!(::LocalTaskBlock{H, H}) where {H} = nothing


function Base.show(io::IO, blk::LocalTaskBlock)
    pos_str = join(Tuple(blk.pos), ',')
    size_str = join(block_size(blk), '×')
    print(io, "LocalTaskBlock(at ($pos_str) of size $size_str, state: $(blk.state.step))")
end

#
# RemoteTaskBlock
#

"""
    RemoteTaskBlock{B} <: TaskBlock{B}

Block located at the border of a [`BlockGrid`](@ref), containing MPI buffer of type `B` for
communication with other [`BlockGrid`](@ref)s.
"""
mutable struct RemoteTaskBlock{B} <: TaskBlock{B}
    pos        :: CartesianIndex{2}  # Position in the local block grid
    neighbour  :: LocalTaskBlock     # Remote blocks are on the edges of the sub-domain: there can only be one real neighbour
    rank       :: Int                # `-1` if the remote block has no MPI rank to communicate with
    global_pos :: CartesianIndex{2}  # Rank position in the Cartesian process grid
    send_buf   :: MPI.Buffer{B}
    recv_buf   :: MPI.Buffer{B}
    requests   :: MPI.UnsafeMultiRequest

    function RemoteTaskBlock{B}(size, pos, rank, global_pos, comm) where {B}
        # `neighbour` is set afterwards, when all blocks are created.
        block = new{B}(pos)
        block.rank = rank
        block.global_pos = global_pos
        block.send_buf = MPI.Buffer(B(undef, size))
        block.recv_buf = MPI.Buffer(B(undef, size))
        block.requests = MPI.UnsafeMultiRequest(2)  # We always keep a reference to the buffers, therefore it is safe
        MPI.Send_init(block.send_buf, comm, block.requests[1]; dest=rank)
        MPI.Recv_init(block.recv_buf, comm, block.requests[2]; source=rank)
        return block
    end

    function RemoteTaskBlock{B}(pos) where {B}
        # Constructor for an non-existant task block, found at the edges of the global domain, where
        # there is no neighbouring MPI rank.
        block = new{B}(pos)
        block.rank = -1
        block.global_pos = CartesianIndex(0, 0)
        block.send_buf = MPI.Buffer(B(undef, 0))
        block.recv_buf = MPI.Buffer(B(undef, 0))
        block.requests = MPI.UnsafeMultiRequest(0)
        return block
    end
end


function Base.show(io::IO, blk::RemoteTaskBlock)
    pos_str = join(Tuple(blk.pos), ',')
    if blk.rank == -1
        print(io, "RemoteTaskBlock(at ($pos_str), to: nowhere)")
    else
        global_str = join(Tuple(blk.global_pos), ',')
        print(io, "RemoteTaskBlock(at ($pos_str), to: process $(blk.rank) at ($global_str))")
    end
end
