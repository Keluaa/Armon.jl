
"""
    TaskBlock{V, Dim}

Abstract block used for cache blocking.
"""
abstract type TaskBlock{V <: AbstractArray, Dim} end

array_type(::TaskBlock{V}) where {V} = V
Base.eltype(::TaskBlock{V}) where {V} = eltype(V)
Base.ndims(::ObjOrType{TaskBlock{V, Dim}}) where {V, Dim} = Dim


"""
    BlockData{V, D}

Holds the variables of all cells of a [`LocalTaskBlock`](@ref).

Variables are separated into scalar variables and dimensional variables.
There is one array per scalar variable, and `D` arrays per dimensional variable.
"""
struct BlockData{V, D}
    scalar_vars :: @NamedTuple begin
        ρ      :: V
        E      :: V
        p      :: V
        c      :: V
        g      :: V
        uˢ     :: V
        pˢ     :: V
        work_1 :: V
        work_2 :: V
        mask   :: V  # TODO: remove ??    
    end

    dim_vars :: @NamedTuple begin
        x      :: NTuple{D, V}
        u      :: NTuple{D, V}
        work_3 :: NTuple{D, V}
    end

    function BlockData{V, D}(size; kwargs...) where {V, D}
        scalars = map(scalar_vars()) do var
            label = string(var)
            return var => V(undef, size; alloc_array_kwargs(; label, kwargs...)...)
        end
        dims = map(dim_vars()) do var
            return var => ntuple(D) do d
                label = "$(var)_$d"
                return V(undef, size; alloc_array_kwargs(; label, kwargs...)...)
            end
        end
        return new{V, D}(NamedTuple(scalars), NamedTuple(dims))
    end
end


scalar_vars() = fieldnames(fieldtype(BlockData{Any, 1}, :scalar_vars))  # Variables with one scalar value
dim_vars()    = fieldnames(fieldtype(BlockData{Any, 1}, :dim_vars))     # Variables with one component (array) per dimension
block_vars()  = (scalar_vars()..., dim_vars()...)       # All variables
main_vars()   = (:x, :ρ, :u, :E, :p, :c, :g, :uˢ, :pˢ)  # Variables synchronized between host and device
saved_vars()  = (:x, :ρ, :u,     :p                  )  # Variables saved to/read from I/O
comm_vars()   = (    :ρ, :u, :E, :p, :c, :g          )  # Variables exchanged between ghost cells, and on which boundary condition are applied


num_arrays_per_cell(dim) = length(scalar_vars()) + length(dim_vars()) * dim
num_arrays_for_vars(vars, dim) = length(vars) + count(in(dim_vars()), vars) * (dim-1)
num_arrays_per_comm(dim) = num_arrays_for_vars(comm_vars(), dim)


"""
    var_arrays(data::BlockData{V, D})
    var_arrays(data::BlockData{V, D}, vars)

A `Tuple` of arrays of type `V` from the block `data`.
By default all arrays in `data` are returned.

`vars` allows to keep only the arrays of the given variable names: it can be a `Tuple` of
`Symbol`s, or a function returning a `Tuple` of `Symbol`s (such as `main_vars`, `saved_vars`
or `comm_vars`), in which case optimal non-allocating code may be generated, as the size of
the tuple can be inferred.
"""
function var_arrays(data::BlockData, var::Symbol)
    var in scalar_vars() && return (getfield(data.scalar_vars, var),)
    var in dim_vars()    && return getfield(data.dim_vars, var)
    return ()
end

var_arrays(data::BlockData) = var_arrays(data, block_vars())
function var_arrays(data::BlockData{V, D}, all_vars::Vars) where {V, D, Vars}
    # Allow passing functions like `comm_vars` or `main_vars` to generate optimal code
    vars_names = Vars <: Function ? all_vars() : all_vars
    splat_concat(t1, t2) = (t1..., t2...)
    return foldl(splat_concat, var_arrays.(Ref(data), vars_names))
end


"""
    var_arrays_names(data::BlockData)
    var_arrays_names(data::BlockData, vars)
    var_arrays_names(::Type{BlockData})
    var_arrays_names(::Type{BlockData}, vars)

Same as [`var_arrays`](@ref), but returns a `Tuple` of unique `Symbol`s for each array.
"""
function var_arrays_names(::Val{Dim}, var::Symbol) where {Dim}
    var in scalar_vars() && return (var,)
    var in dim_vars()    && return ntuple(d -> Symbol(var, :_, d), Dim)
    return ()
end

var_arrays_names(bd::ObjOrType{BlockData}) = var_arrays_names(bd, block_vars())
function var_arrays_names(::ObjOrType{<:BlockData{V, D}}, all_vars::Vars) where {V, D, Vars}
    vars_names = Vars <: Function ? all_vars() : all_vars
    splat_concat(t1, t2) = (t1..., t2...)
    return foldl(splat_concat, var_arrays_names.(Val(D), vars_names))
end


main_arrays(data::BlockData)  = var_arrays(data, main_vars)
saved_arrays(data::BlockData) = var_arrays(data, saved_vars)
comm_arrays(data::BlockData)  = var_arrays(data, comm_vars)

Base.ndims(::ObjOrType{<:BlockData{<:Any, D}}) where {D} = D


"""
    LocalTaskBlock{D, H, Dim, Size, SolverState} <: TaskBlock{V, Dim}

Container of `Size` and variables of type `D` on the device and `H` on the host.
Part of a [`BlockGrid`](@ref).

The block stores its own solver state, allowing it to run all solver steps independantly of all
other blocks, apart from steps requiring synchronization with neighbours or the whole grid.
"""
mutable struct LocalTaskBlock{D, H, Dim, Size <: BlockSize{Dim}, SState <: SolverState} <: TaskBlock{D, Dim}
    state        :: SState               # Solver state for the block
    size         :: Size                 # Size (in cells) of the block
    pos          :: CartesianIndex{Dim}  # Position in the local block grid
    neighbours   :: Neighbours{TaskBlock, Dim}
    exchanges    :: Neighbours{BlockInterface, Dim}  # State of ghost cells exchanges for each side
    device_data  :: BlockData{D, Dim}
    host_data    :: BlockData{H, Dim}  # Host data uses the same arrays as device data if `D == H`
    # TODO: device? storing the associated GPU stream, or CPU cores (or maybe not, to allow relocations?)

    function LocalTaskBlock{D, H, Dim, Size, SState}(
        size::Size, pos, blk_state::SState, device_kwargs, host_kwargs
    ) where {D, H, Dim, Size <: BlockSize{Dim}, SState}
        # `exchanges` and `neighbours` are set afterwards, after all blocks are created.
        block = new{D, H, Dim, Size, SState}(blk_state, size, pos, #= undef =#)
        cell_count = prod(block_size(size))
        if D != H
            block.device_data = BlockData{D, Dim}(cell_count; device_kwargs...)
            block.host_data   = BlockData{H, Dim}(cell_count; host_kwargs...)
        else
            block.host_data = block.device_data = BlockData{D, Dim}(cell_count; device_kwargs...)
        end
        return block
    end
end


block_size(blk::LocalTaskBlock) = block_size(blk.size)
real_block_size(blk::LocalTaskBlock) = real_block_size(blk.size)
ghosts(blk::LocalTaskBlock) = ghosts(blk.size)
solver_state(blk::LocalTaskBlock) = blk.state
Base.ndims(blk::LocalTaskBlock) = ndims(blk.size)

block_device_data(blk::LocalTaskBlock) = blk.device_data
block_host_data(blk::LocalTaskBlock) = blk.host_data
block_data(blk::LocalTaskBlock; on_device=true) = on_device ? block_device_data(blk) : block_host_data(blk)

var_arrays(blk::LocalTaskBlock, vars; on_device=true) = var_arrays(block_data(blk; on_device), vars)
main_arrays(blk::LocalTaskBlock; on_device=true)  = main_arrays(block_data(blk; on_device))
saved_arrays(blk::LocalTaskBlock; on_device=true) = saved_arrays(block_data(blk; on_device))
comm_arrays(blk::LocalTaskBlock; on_device=true)  = comm_arrays(block_data(blk; on_device))


function reset!(blk::LocalTaskBlock)
    reset!(blk.state)
    foreach(reset!, blk.exchanges)
end


"""
    device_to_host!(blk::LocalTaskBlock)

Copies the device data of `blk` to the host data. A no-op if the device is the host.
"""
function device_to_host!(blk::LocalTaskBlock{D, H}) where {D, H}
    for (dst_var, src_var) in zip(main_arrays(blk.host_data), main_arrays(blk.device_data))
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
    for (dst_var, src_var) in zip(main_arrays(blk.device_data), main_arrays(blk.host_data))
        # TODO: ensure this is done in the correct thread/device
        # KernelAbstractions.copyto!(blk.device, dst_var, src_var)
        copyto!(dst_var, src_var)
    end
end

host_to_device!(::LocalTaskBlock{H, H}) where {H} = nothing


function move_pages(blk::LocalTaskBlock, target_node)
    for array in var_arrays(blk)
        move_pages(array, target_node)
    end
end


lock_pages(blk::LocalTaskBlock) = foreach(lock_pages, var_arrays(blk))


function Base.show(io::IO, blk::LocalTaskBlock)
    pos_str = join(Tuple(blk.pos), ',')
    size_str = join(block_size(blk), '×')
    print(io, "LocalTaskBlock(at ($pos_str) of size $size_str, state: $(blk.state.step))")
end

#
# RemoteTaskBlock
#

"""
    RemoteTaskBlock{B, Dim} <: TaskBlock{B, Dim}

Block located at the border of a [`BlockGrid`](@ref), containing MPI buffers of type `B` for
communication with another [`BlockGrid`](@ref) in a distant process.
"""
mutable struct RemoteTaskBlock{B, Dim} <: TaskBlock{B, Dim}
    pos        :: CartesianIndex{Dim}  # Position in the local block grid
    neighbour  :: LocalTaskBlock       # Remote blocks are on the edges of the sub-domain: there can only be one real neighbour
    rank       :: Int                  # `-1` if the remote block has no MPI rank to communicate with
    global_pos :: CartesianIndex{Dim}  # Rank position in the Cartesian process grid
    send_buf   :: MPI.Buffer{B}
    recv_buf   :: MPI.Buffer{B}
    requests   :: MPI.UnsafeMultiRequest

    function RemoteTaskBlock{B, Dim}(size, pos, rank, global_pos, comm, side_idx) where {B, Dim}
        # `neighbour` is set afterwards, when all blocks are created.
        block = new{B, Dim}(pos)
        block.rank = rank
        block.global_pos = global_pos
        block.send_buf = MPI.Buffer(B(undef, size))
        block.recv_buf = MPI.Buffer(B(undef, size))
        block.requests = MPI.UnsafeMultiRequest(2)  # We always keep a reference to the buffers, therefore it is safe

        # Because two ranks may have several comms at once, we must use tags. They must match at both sides.
        # Since both ranks share the same (flat) side and block size, a unique index could be the block
        # position along the communication side.
        # TODO: this is not enough to support arbitrary block distributions (non-flat sides), there
        #   could be tag collisions in this case
        #   hashes are not a viable alternative, as `MPI.tab_ub()` is only 2^15 at min (2^23 for OpenMPI)
        MPI.Send_init(block.send_buf, comm, block.requests[1]; dest=rank, tag=side_idx)
        MPI.Recv_init(block.recv_buf, comm, block.requests[2]; source=rank, tag=side_idx)

        return block
    end

    function RemoteTaskBlock{B, Dim}(pos) where {B, Dim}
        # Constructor for an non-existant remote block, found at the edges of the global domain, where
        # there is no neighbouring MPI rank.
        block = new{B, Dim}(pos)
        block.rank = -1
        block.global_pos = zero(CartesianIndex{Dim})
        block.send_buf = MPI.Buffer(B(undef, 0))
        block.recv_buf = MPI.Buffer(B(undef, 0))
        block.requests = MPI.UnsafeMultiRequest(0)
        return block
    end
end


function move_pages(blk::RemoteTaskBlock, target_node)
    blk.rank == -1 && return
    # MPI communications might not have happened, therefore pages might not be placed on a NUMA node
    # yet, so it is safer to touch them first.
    touch_pages(array_pages(blk.send_buf.data))
    touch_pages(array_pages(blk.recv_buf.data))
    move_pages(blk.send_buf.data, target_node)
    move_pages(blk.recv_buf.data, target_node)
end


function lock_pages(blk::RemoteTaskBlock)
    blk.rank == -1 && return
    lock_pages(blk.send_buf.data)
    lock_pages(blk.recv_buf.data)
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
