
const ObjOrType = Union{T, Type{<:T}} where T


mutable struct Atomic{T}
    @atomic x::T
end


"""
    Axis

Enumeration of the axes of a domain. Accepts values from `1` (X axis) to `126`.
Only the first 3 axes are explicitly named.
"""
module Axis
    primitive type T <: Base.Enums.Enum{UInt8} 8 end 

    function T(a::Integer)
        !(1 ≤ a ≤ 126) && throw(ArgumentError("Axis index must be ≥ 1 and ≤ 126"))
        return Core.bitcast(T, Base.convert(UInt8, a))
    end

    (::Type{UInt8})(a::T) = UInt8(Core.bitcast(UInt8, a))::UInt8
    Base.cconvert(::Type{UInt8}, a::T) = UInt8(a)::UInt8

    Base.typemin(::Type{T}) = T(UInt8(1))
    Base.typemax(::Type{T}) = T(UInt8(126))  # Nearest multiple of two before `typemax(UInt8)÷2`

    axis_name(a::UInt8) = 1 ≤ a ≤ 3 ? (:X, :Y, :Z)[a] : Symbol(:Axis_, a)
    axis_name(a::T) = axis_name(UInt8(a))
    Base.Symbol(a::T) = axis_name(a)
    Base.Enums._symbol(a::T) = axis_name(a)

    Base.show(io::IO, ::MIME"text/plain", a::T) = print(io, a, "::Axis.T = ", Integer(a))

    function Base.show(io::IO, ::MIME"text/plain", st::Type{T})
        println(io, "Enum Axis:")
        for x in 1:7
            println(io, " Axis.", Symbol(T(x)), " = ", x)
        end
        print(io, " ... (up to Axis.$(typemax(T)))")
    end

    let type_hash = hash(T)
        Base.Enums._enum_hash(a::T, h::UInt) = Base.Enums.hash(type_hash, Base.Enums.hash(Integer(a), h))
    end

    Base.instances(::Type{T}) = error("instances of Axis should not be iterated, use `axes_of(dim)`")
    Base.Enums.namemap(::Type{T}) = error("names of Axis should not be iterated, use `Symbol.(axes_of(dim))`")

    const X = T(1)
    const Y = T(2)
    const Z = T(3)
end


"""
    axes_of(::Val{dim})
    axes_of(dim::Integer)

Axes for all dimensions up to `dim`.
    
```jldoctest
julia> Armon.axes_of(3)
(Armon.Axis.X, Armon.Axis.Y, Armon.Axis.Z)
```
"""
axes_of(::Val{dim}) where {dim} = ntuple(Axis.T, dim)
axes_of(dim::Integer) = axes_of(Val(dim))


"""
    next_axis(axis::Axis.T, ::Val{dim})
    next_axis(axis::Axis.T, dim::Integer)

The axis after `axis`, mod `dim`.

```jldoctest
julia> Armon.next_axis(Armon.Axis.X, 3)
Y::Axis.T = 2

julia> Armon.next_axis(Armon.Axis.Z, 3)
X::Axis.T = 1
```
"""
next_axis(axis::Axis.T, ::Val{dim}) where {dim} = Axis.T(mod1(UInt(axis) + 1, dim))
next_axis(axis::Axis.T, dim::Integer) = next_axis(axis, Val(dim))


offset_to(axis::Axis.T) = ntuple(i -> i == Integer(axis) ? 1 : 0, dim)


"""
    Side

Enumeration of the sides of a domain. Accepts values from `1` (Left side) to `254`.
Only the sides of the first 3 axes are explicitly named.
"""
module Side
    primitive type T <: Base.Enums.Enum{UInt8} 8 end 

    function T(s::Integer)
        !(1 ≤ s ≤ 254) && throw(ArgumentError("Side index must be ≥ 1 and ≤ 254"))
        return Core.bitcast(T, Base.convert(UInt8, s))
    end

    (::Type{UInt8})(s::T) = UInt8(Core.bitcast(UInt8, s))::UInt8
    Base.cconvert(::Type{UInt8}, s::T) = UInt8(s)::UInt8

    Base.typemin(::Type{T}) = T(UInt8(1))
    Base.typemax(::Type{T}) = T(UInt8(254))  # Nearest multiple of two before `typemax(UInt8)`

    side_name(s::UInt8) = 1 ≤ s ≤ 6 ? (:Left, :Right, :Bottom, :Top, :Back, :Front)[s] : Symbol(:Side_, s)
    side_name(s::T) = side_name(UInt8(s))
    Base.Symbol(s::T) = side_name(s)
    Base.Enums._symbol(s::T) = side_name(s)

    Base.show(io::IO, ::MIME"text/plain", s::T) = print(io, s, "::Side.T = ", Integer(s))

    function Base.show(io::IO, ::MIME"text/plain", st::Type{T})
        println(io, "Enum Side:")
        for x in 1:7
            println(io, " Side.", Symbol(T(x)), " = ", x)
        end
        print(io, " ... (up to Side.$(typemax(T)))")
    end

    let type_hash = hash(T)
        Base.Enums._enum_hash(s::T, h::UInt) = Base.Enums.hash(type_hash, Base.Enums.hash(Integer(s), h))
    end

    Base.instances(::Type{T}) = error("instances of Side should not be iterated, use `sides_of(dim)`")
    Base.Enums.namemap(::Type{T}) = error("names of Side should not be iterated, use `Symbol.(sides_of(dim))`")

    const Left   = T(1)
    const Right  = T(2)
    const Bottom = T(3)
    const Top    = T(4)
    const Back   = T(5)
    const Front  = T(6)
end


"""
    sides_of(::Val{dim})
    sides_of(dim::Integer)

Sides for all dimensions up to `dim`.

```jldoctest
julia> Armon.sides_of(2)
(Armon.Side.Left, Armon.Side.Right, Armon.Side.Bottom, Armon.Side.Top)
```
"""
sides_of(::Val{dim}) where {dim} = ntuple(Side.T, dim*2)
sides_of(dim::Integer) = sides_of(Val(dim))


"""
    sides_along(ax::Axis.T)

Sides before and after `ax`.

```jldoctest
julia> Armon.sides_along(Armon.Axis.X)
(Armon.Side.Left, Armon.Side.Right)

julia> Armon.sides_along(Armon.Axis.Z)
(Armon.Side.Back, Armon.Side.Front)
```
"""
sides_along(ax::Axis.T) = (Side.T(2*UInt8(ax)-1), Side.T(2*UInt8(ax)))


"""
    first_side(ax::Axis.T)

The first side of `sides_along(ax)`. Always at the lowest coordinates of `ax`.
"""
first_side(ax::Axis.T) = Side.T(2*UInt8(ax)-1)


"""
    first_side(s::Side.T)

`true` if `s` is the first side of an axis.
"""
first_side(s::Side.T) = isodd(Integer(s))


"""
    first_sides(::Val{dim})
    first_sides(dim::Integer)

Tuple of the first sides of all axes up to `dim`.

```jldoctest
julia> Armon.first_sides(3)
(Armon.Side.Left, Armon.Side.Bottom, Armon.Side.Back)
```
"""
first_sides(::Val{dim}) where {dim} = ntuple(d -> Side.T(2*d-1), dim)
first_sides(dim::Integer) = first_sides(Val(dim))


"""
    last_side(ax::Axis.T)

The last side of `sides_along(ax)`. Always at the highest coordinates of `ax`.
"""
last_side(ax::Axis.T) = Side.T(2*UInt8(ax))


"""
    last_side(ax::Axis.T)

`true` if `s` is the last side of an axis.
"""
last_side(s::Side.T) = iseven(Integer(s))


"""
    last_sides(::Val{dim})
    last_sides(dim::Integer)

Tuple of the last sides of all axes up to `dim`.

```jldoctest
julia> Armon.last_sides(3)
(Armon.Side.Right, Armon.Side.Top, Armon.Side.Front)
```
"""
last_sides(::Val{dim}) where {dim} = ntuple(d -> Side.T(2*d), dim)
last_sides(dim::Integer) = last_sides(Val(dim))


"""
    axis_of(s::Side.T)

[`Axis`](@ref) this [`Side`](@ref) belongs to.

```jldoctest
julia> Armon.axis_of(Armon.Side.Left)
X::Axis.T = 1
```
"""
axis_of(s::Side.T) = Axis.T((Integer(s) - 1) >> 1 + 1)


"""
    opposite_of(s::Side.T)

The [`Side`](@ref) opposite of `s`.

```jldoctest
julia> Armon.opposite_of(Armon.Side.Left)
Right::Side.T = 2
```
"""
opposite_of(s::Side.T) = first_side(s) ? Side.T(Integer(s)+1) : Side.T(Integer(s)-1)


"""
    offset_of(s::Side.T, dim::Integer)

An `NTuple{dim, Int}` offset going towards `s`.

```jldoctest
julia> Armon.offset_to(Armon.Side.Left, 2)
(-1, 0)

julia> Armon.offset_to(Armon.Side.Top, 3)
(0, 1, 0)
```
"""
function offset_to(s::Side.T, dim::Integer)
    i_ax = Integer(axis_of(s))
    if first_side(s)
        return ntuple(i -> i == i_ax ? -1 : 0, dim)
    else
        return ntuple(i -> i == i_ax ?  1 : 0, dim)
    end
end


"""
    side_from_offset(offset::NTuple{D, <:Integer})

Opposite of [`offset_to`](@ref): deduces the [`Side`](@ref) from the offset.
If `offset` is all zeros, `Side.Left` is returned.

```jldoctest
julia> Armon.side_from_offset((-1, 0))
Left::Side.T = 1

julia> Armon.side_from_offset((0, 1, 0))
Top::Side.T = 4
```
"""
function side_from_offset(offset::NTuple{D, <:Integer}) where {D}
    for (i_ax, o) in enumerate(offset)
        o < 0 && return first_side(Axis.T(i_ax))
        o > 0 && return last_side(Axis.T(i_ax))
    end
    return Side.Left  # Should not happen, here only for type-stability (responsability of the caller)
end


"""
    Neighbours{T, D}

Immutable `Vector`-like object storing one `T` for each of [`sides_of(D)`](@ref).
Can be indexed using [`Side`](@ref) (returns a `T`), [`Axis`](@ref) or `Int` (returns a `NTuple{2, T}`,
one `T` for each side of the axis).

```jldoctest
julia> n = Armon.Neighbours((ax, s) -> s, 2)
2-element Armon.Neighbours{Armon.Side.T, 2}:
 (Armon.Side.Left, Armon.Side.Right)
 (Armon.Side.Bottom, Armon.Side.Top)

julia> n[Armon.Side.Left]
Left::Side.T = 1

julia> n[Armon.Axis.X]
(Armon.Side.Left, Armon.Side.Right)
```
"""
struct Neighbours{T, D} <: AbstractVector{NTuple{2, T}}
    val :: NTuple{D, NTuple{2, T}}

    Neighbours{T, D}(val) where {T, D} = new{T, D}(val)
    Neighbours(val::NTuple{D, NTuple{2, T}}) where {T, D} = Neighbours{T, D}(val)
end

Neighbours(f::Base.Callable, axes::Tuple{Vararg{Axis.T}}) = Neighbours(f.(axes))
Neighbours(f::Base.Callable, dim::Integer) = Neighbours(ax -> f.(ax, sides_along(ax)), axes_of(dim))

Base.@propagate_inbounds Base.getindex(neigh::Neighbours, axis::Axis.T) = neigh.val[Integer(axis)]
Base.@propagate_inbounds Base.getindex(neigh::Neighbours, side::Side.T) = neigh[axis_of(side)][1+last_side(side)]
Base.@propagate_inbounds Base.getindex(neigh::Neighbours, i::Int) = neigh.val[i]

Base.IndexStyle(::Type{<:Neighbours}) = Base.IndexLinear()
Base.size(::Neighbours{T, D}) where {T, D} = (D,)


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
