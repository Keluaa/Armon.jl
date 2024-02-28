
const ObjOrType = Union{T, Type{<:T}} where T


"""
    Axis

Enumeration of the axes of the domain
"""
@enum Axis X_axis=1 Y_axis=2

next_axis(axis::Axis) = Axis(mod1(Int(axis) + 1, length(instances(Axis))))


"""
    Side

Enumeration of the sides of the domain
"""
@enum Side Left=1 Right=2 Bottom=3 Top=4

@inline sides_along(dir::Axis) = dir == X_axis ? (Left, Right) : (Bottom, Top)

@inline first_side(dir::Axis) = dir == X_axis ? Left : Bottom
@inline first_sides() = (Left, Bottom)

@inline last_side(dir::Axis) = dir == X_axis ? Right : Top
@inline last_sides() = (Right, Top)

@inline function axis_of(side::Side)
    # TODO: `Axis(((Integer(side) - 1) >> 1) + 1)` ??
    if     side == Left   return X_axis
    elseif side == Right  return X_axis
    elseif side == Bottom return Y_axis
    elseif side == Top    return Y_axis
    else                  return X_axis  # Should not happen, here only for type-stability
    end
end

@inline function opposite_of(side::Side)
    if     side == Left   return Right
    elseif side == Right  return Left
    elseif side == Bottom return Top
    elseif side == Top    return Bottom
    else                  return Left  # Should not happen, here only for type-stability
    end
end

@inline function offset_to(side::Side)
    if     side == Left   return (-1,  0)
    elseif side == Right  return ( 1,  0)
    elseif side == Bottom return ( 0, -1)
    elseif side == Top    return ( 0,  1)
    else                  return ( 0,  0)  # Should not happen, here only for type-stability
    end
end


@inline function side_from_offset(offset::Tuple)
    if     offset[1] < 0 return Left
    elseif offset[1] > 0 return Right
    elseif offset[2] < 0 return Bottom
    elseif offset[2] > 0 return Top
    else                 return Left  # Should not happen, here only for type-stability (responsability of the caller)
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


disp_blk(blk, var) = reshape(getfield(blk, var), block_size(blk))'
disp_real_blk(blk, var) = view(disp_blk(blk, var)', (.+).(Base.oneto.(real_block_size(blk.size)), ghosts(blk.size))...)'
disp_mirror_y(A) = view(A, size(A, 1):-1:1, :)  # places the bottom-left cell at the bottom-left of the display
