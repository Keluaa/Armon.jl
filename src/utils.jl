
const ObjOrType = Union{T, Type{<:T}} where T


mutable struct Atomic{T}
    @atomic x::T
end


"""
    Axis

Enumeration of the axes of the domain
"""
@enumx Axis::UInt8 X=1 Y=2

next_axis(axis::Axis.T) = Axis.T(mod1(Int(axis) + 1, length(instances(Axis.T))))


"""
    Side

Enumeration of the sides of the domain
"""
@enumx Side Left=1 Right=2 Bottom=3 Top=4

@inline sides_along(dir::Axis.T) = dir == Axis.X ? (Side.Left, Side.Right) : (Side.Bottom, Side.Top)

@inline first_side(dir::Axis.T) = dir == Axis.X ? Side.Left : Side.Bottom
@inline first_sides() = (Side.Left, Side.Bottom)

@inline last_side(dir::Axis.T) = dir == Axis.X ? Side.Right : Side.Top
@inline last_sides() = (Side.Right, Side.Top)

@inline function axis_of(side::Side.T)
    # TODO: `Axis.T(((Integer(side) - 1) >> 1) + 1)` ??
    if     side == Side.Left   return Axis.X
    elseif side == Side.Right  return Axis.X
    elseif side == Side.Bottom return Axis.Y
    elseif side == Side.Top    return Axis.Y
    else                       return Axis.X  # Should not happen, here only for type-stability
    end
end

@inline function opposite_of(side::Side.T)
    if     side == Side.Left   return Side.Right
    elseif side == Side.Right  return Side.Left
    elseif side == Side.Bottom return Side.Top
    elseif side == Side.Top    return Side.Bottom
    else                       return Side.Left  # Should not happen, here only for type-stability
    end
end

@inline function offset_to(side::Side.T)
    if     side == Side.Left   return (-1,  0)
    elseif side == Side.Right  return ( 1,  0)
    elseif side == Side.Bottom return ( 0, -1)
    elseif side == Side.Top    return ( 0,  1)
    else                       return ( 0,  0)  # Should not happen, here only for type-stability
    end
end


@inline function side_from_offset(offset::Tuple)
    if     offset[1] < 0 return Side.Left
    elseif offset[1] > 0 return Side.Right
    elseif offset[2] < 0 return Side.Bottom
    elseif offset[2] > 0 return Side.Top
    else                 return Side.Left  # Should not happen, here only for type-stability (responsability of the caller)
    end
end


"""
    CPU_HP

Device tag for the high-performance CPU backend using multithreading (with Polyester.jl) and
vectorisation.
"""
struct CPU_HP end


"""
    SolverException

Thrown when the solver encounters an invalid state.

The `category` field can be used to distinguish between error types without inspecting the error
message:
 - `:config`: a problem in the solver configuration, usually thrown when constructing `ArmonParameters`
 - `:cpp`: a C++ exception thrown by the C++ Kokkos backend
 - `:time`: an invalid time step

`ErrorException`s thrown by the solver represent internal errors.
"""
struct SolverException <: Exception
    category::Symbol
    msg::String
end


solver_error(category::Symbol, msg::String) = throw(SolverException(category, msg))

function solver_error(category::Symbol, msgs::Vararg{Any, N}) where {N}
    throw(SolverException(category, Base.string(msgs...)))
end


function Base.showerror(io::IO, ex::SolverException)
    print(io, "SolverException ($(ex.category)): ", ex.msg)
end


disp_blk(blk, var; on_device=true) = reshape(getfield(block_data(blk; on_device), var), block_size(blk))'
disp_real_blk(blk, var; on_device=true) = view(disp_blk(blk, var; on_device)', (.+).(Base.oneto.(real_block_size(blk.size)), ghosts(blk.size))...)'
disp_grid_state(grid) = permutedims(reshape(solver_step.(solver_state.(all_blocks(grid))), grid.grid_size))
disp_mirror_y(A) = view(A, size(A, 1):-1:1, :)  # places the bottom-left cell at the bottom-left of the display


function IAllreduce!(rbuf::MPI.RBuffer, op::Union{MPI.Op, MPI.MPI_Op}, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request())
    @assert MPI.isnull(req)
    # int MPI_Allreduce(const void* sendbuf, void* recvbuf, int count,
    #                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
    #                   MPI_Request* req)
    MPI.API.MPI_Iallreduce(rbuf.senddata, rbuf.recvdata, rbuf.count, rbuf.datatype, op, comm, req)
    MPI.setbuffer!(req, rbuf)
    return req
end

IAllreduce!(rbuf::MPI.RBuffer, op, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request()) =
    IAllreduce!(rbuf, MPI.Op(op, eltype(rbuf)), comm, req)
IAllreduce!(sendbuf, recvbuf, op, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request()) =
    IAllreduce!(MPI.RBuffer(sendbuf, recvbuf), op, comm, req)

# inplace
IAllreduce!(rbuf, op, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request()) =
    IAllreduce!(MPI.IN_PLACE, rbuf, op, comm, req)
