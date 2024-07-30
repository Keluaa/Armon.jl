
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

    @noinline dim_error(d::Integer) = throw(ArgumentError("Axis index must be ≥ 1 and ≤ 127, got: $d"))
    @inline check_dim(d::Integer) = ((1 ≤ d ≤ 127) || dim_error(d); UInt8(d))

    function T(a::Integer)
        check_dim(a)
        return Core.bitcast(T, Base.convert(UInt8, a))
    end

    (::Type{UInt8})(a::T) = UInt8(Core.bitcast(UInt8, a))::UInt8
    Base.cconvert(::Type{UInt8}, a::T) = UInt8(a)::UInt8
    unsafe_cast(a::UInt8) = Core.bitcast(T, a)

    Base.typemin(::Type{T}) = T(UInt8(1))
    Base.typemax(::Type{T}) = T(UInt8(127))  # `typemax(UInt8)÷2` for compat with `Side`

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


# Make Axis.T a valid index to any array
Base.checkindex(::Type{Bool}, inds::AbstractUnitRange, axis::Axis.T) = Base.checkindex(Bool, inds, Int(axis))
Base.to_index(axis::Axis.T) = Int(axis)


"""
    axes_of(::Val{dim})
    axes_of(dim::Integer)

Axes for all dimensions up to `dim`.

```jldoctest
julia> Armon.axes_of(3)
(Armon.Axis.X, Armon.Axis.Y, Armon.Axis.Z)
```
"""
axes_of(dim::Integer) = ntuple(Axis.T, dim)
axes_of(::Val{dim}) where {dim} = axes_of(dim)


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
next_axis(axis::Axis.T, dim::Integer) = Axis.unsafe_cast(mod1(UInt8(axis) + 0x1, Axis.check_dim(dim)))
next_axis(axis::Axis.T, ::Val{dim}) where {dim} = next_axis(axis, dim)


offset_to(axis::Axis.T, ::Val{dim}, mag=1) where {dim} = offset_to(axis, dim, mag)
function offset_to(axis::Axis.T, dim::Integer, mag::T=1) where {T}
    i_ax = Int64(axis)
    return ntuple(i -> i == i_ax ? mag : zero(T), dim)
end


"""
    Side

Enumeration of the sides of a domain. Accepts values from `1` (Left side) to `254`.
Only the sides of the first 3 axes are explicitly named.
"""
module Side
    primitive type T <: Base.Enums.Enum{UInt8} 8 end

    @noinline side_error(s::Integer) = throw(ArgumentError("Side index must be ≥ 1 and ≤ 254, got: $s"))
    @inline side_check(s::Integer) = (1 ≤ s ≤ 254) || side_error(s)

    function T(s::Integer)
        side_check(s)
        return Core.bitcast(T, Base.convert(UInt8, s))
    end

    (::Type{UInt8})(s::T) = UInt8(Core.bitcast(UInt8, s))::UInt8
    Base.cconvert(::Type{UInt8}, s::T) = UInt8(s)::UInt8
    unsafe_cast(s::UInt8) = Core.bitcast(T, s)

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

    # The names respect the right-hand rule, with the Z axis facing you, and the Y axis going up.
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
sides_of(dim::Integer) = ntuple(Side.T, dim*2)
sides_of(::Val{dim}) where {dim} = sides_of(dim)


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
sides_along(ax::Axis.T) = (Side.unsafe_cast(0x2*UInt8(ax)-0x1), Side.unsafe_cast(0x2*UInt8(ax)))


"""
    first_side(ax::Axis.T)

The first side of `sides_along(ax)`. Always at the lowest coordinates of `ax`.
"""
first_side(ax::Axis.T) = Side.unsafe_cast(0x2*UInt8(ax)-0x1)


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
(Armon.Side.Left, Armon.Side.Bottom, Armon.Side.Front)
```
"""
first_sides(dim::Integer) = ntuple(d -> Side.T(2*d-1), dim)
first_sides(::Val{dim}) where {dim} = first_sides(dim)


"""
    last_side(ax::Axis.T)

The last side of `sides_along(ax)`. Always at the highest coordinates of `ax`.
"""
last_side(ax::Axis.T) = Side.unsafe_cast(0x2*UInt8(ax))


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
last_sides(dim::Integer) = ntuple(d -> Side.T(2*d), dim)
last_sides(::Val{dim}) where {dim} = last_sides(dim)


"""
    axis_of(s::Side.T)

[`Axis`](@ref) this [`Side`](@ref) belongs to.

```jldoctest
julia> Armon.axis_of(Armon.Side.Left)
X::Axis.T = 1
```
"""
axis_of(s::Side.T) = Axis.unsafe_cast((UInt8(s) - 0x1) >> 1 + 0x1)


"""
    opposite_of(s::Side.T)

The [`Side`](@ref) opposite of `s`.

```jldoctest
julia> Armon.opposite_of(Armon.Side.Left)
Right::Side.T = 2
```
"""
opposite_of(s::Side.T) = Side.unsafe_cast(((UInt8(s) - 0x1) ⊻ 0x1) + 0x1)


"""
    offset_of(s::Side.T, dim::Integer, magnitude=1)

An `NTuple{dim, Int}` offset going towards `s`.

```jldoctest
julia> Armon.offset_to(Armon.Side.Left, 2)
(-1, 0)

julia> Armon.offset_to(Armon.Side.Top, 3)
(0, 1, 0)

julia> Armon.offset_to(Armon.Side.Top, 3, 5)
(0, 5, 0)
```
"""
offset_to(s::Side.T, ::Val{dim}, mag=1) where {dim} = offset_to(s, dim, mag)
function offset_to(s::Side.T, dim::Integer, mag::T=1) where {T}
    i_ax = Int64(axis_of(s))
    if first_side(s)
        return ntuple(i -> i == i_ax ? -mag : zero(T), dim)
    else
        return ntuple(i -> i == i_ax ?  mag : zero(T), dim)
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
Can be indexed using [`Side`](@ref) (returns a `T`), [`Axis`](@ref) or `Int`
(returns a `NTuple{2, T}`, one `T` for each side of the axis).

Custom accessors for each axis and side of the first 3 dimensions are also defined.

```jldoctest
julia> n = Armon.Neighbours((ax, s) -> s, 2)
2-element Armon.Neighbours{Armon.Side.T, 2}:
 (Armon.Side.Left, Armon.Side.Right)
 (Armon.Side.Bottom, Armon.Side.Top)

julia> n[Armon.Side.Left]
Left::Side.T = 1

julia> n[Armon.Axis.X]
(Armon.Side.Left, Armon.Side.Right)

julia> n.Left
Left::Side.T = 1

julia> n.X
(Armon.Side.Left, Armon.Side.Right)
```
"""
struct Neighbours{T, D} <: AbstractVector{NTuple{2, T}}
    val :: NTuple{D, NTuple{2, T}}

    Neighbours{T, D}(val) where {T, D} = new{T, D}(val)
    Neighbours(val::NTuple{D, NTuple{2, T}}) where {T, D} = Neighbours{T, D}(val)
    Neighbours(val::NTuple{0}) = Neighbours{Nothing, 0}(val)
end

function Neighbours(f::Base.Callable, dim::Union{Integer, Val})
    axis_generator = ax -> f.(ax, sides_along(ax))
    axes = axes_of(dim)
    return Neighbours(axis_generator.(axes))
end

function Neighbours(f::Base.Callable, dim::Union{Integer, Val}, ::Type{T}) where {T}
    axis_generator = ax -> (f.(ax, sides_along(ax))::NTuple{2, T})
    axes = axes_of(dim)
    return Neighbours{T, length(axes)}(axis_generator.(axes))
end

Base.@propagate_inbounds Base.getindex(neigh::Neighbours, side::Side.T) = neigh[axis_of(side)][1+last_side(side)]
Base.@propagate_inbounds Base.getindex(neigh::Neighbours, i::Integer)   = neigh.val[i]

Base.IndexStyle(::Type{<:Neighbours}) = Base.IndexLinear()
Base.size(::Neighbours{T, D}) where {T, D} = (D,)

function Base.getproperty(neigh::Neighbours, sym::Symbol)
    if     sym === :X      return neigh[Axis.X]
    elseif sym === :Y      return neigh[Axis.Y]
    elseif sym === :Z      return neigh[Axis.Z]
    elseif sym === :Left   return neigh[Side.Left]
    elseif sym === :Right  return neigh[Side.Right]
    elseif sym === :Bottom return neigh[Side.Bottom]
    elseif sym === :Top    return neigh[Side.Top]
    elseif sym === :Back   return neigh[Side.Back]
    elseif sym === :Front  return neigh[Side.Front]
    else                   return getfield(neigh, sym)
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


@noinline solver_error(category::Symbol, msg::String) = throw(SolverException(category, msg))

@noinline function solver_error(category::Symbol, msgs::Vararg{Any, N}) where {N}
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
