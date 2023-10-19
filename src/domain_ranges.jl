
import Base: oneto, length, size, axes, isempty, first, last, in, show

#
# Basic indexing macros
#

"""
    @indexing_vars(params)

Brings the parameters needed for the `@i` macro into the current scope.
"""
macro indexing_vars(params)
    return esc(quote
        (; index_start, row_length, col_length, idx_row, idx_col) = $(params)
    end)
end

"""
    @i(i, j)

Converts the two-dimensional indexes `i` and `j` to a mono-dimensional index.

```julia
    idx = @i(i, j)
```
"""
macro i(i, j)
    return esc(quote
        index_start + $(j) * idx_row + $(i) * idx_col
    end)
end

#
# Range utilities
#

shift(r::AbstractUnitRange,   n::Int = 1) = r .+ n
shift(r::OrdinalRange,        n::Int = 1) = r .+ (step(r) * n)
expand(r::AbstractUnitRange,  n::Int = 1) = first(r):(last(r)+n)
expand(r::OrdinalRange,       n::Int = 1) = first(r):step(r):(last(r)+step(r)*n)
prepend(r::AbstractUnitRange, n::Int = 1) = (first(r)-n):last(r)
prepend(r::OrdinalRange,      n::Int = 1) = (first(r)-step(r)*n):step(r):last(r)
inflate(r::AbstractUnitRange, n::Int = 1) = (first(r)-n):(last(r)+n)
inflate(r::OrdinalRange,      n::Int = 1) = (first(r)-step(r)*n):step(r):(last(r)+step(r)*n)

function connect_ranges(r1::AbstractUnitRange, r2::AbstractUnitRange)
    last(r1) > first(r2) && return connect_ranges(r2, r1)
    !isempty(intersect(r1, r2)) && error("$r1 and $r2 intersect")
    last(r1) + 1 != first(r2) && error("$r1 doesn't touch $r2")
    return first(r1):last(r2)
end

function connect_ranges(r1::OrdinalRange, r2::OrdinalRange)
    last(r1) > first(r2) && return connect_ranges(r2, r1)
    !isempty(intersect(r1, r2)) && error("$r1 and $r2 intersect")
    step(r1) != step(r2) && error("$r1 and $r2 have different steps")
    last(r1) + step(r1) != first(r2) && error("$r1 doesn't touch $r2")
    return first(r1):step(r1):last(r2)
end

#
# DomainRange: Two dimensional range to index a 2D array stored with contiguous rows
#

struct DomainRange
    col::StepRange{Int, Int}
    row::StepRange{Int, Int}
end

DomainRange() = DomainRange(1:1:0, 1:1:0)

length(dr::DomainRange) = length(dr.col) * length(dr.row)
size(dr::DomainRange) = (length(dr.col), length(dr.row))
axes(dr::DomainRange) = (oneto(length(dr.col)), oneto(length(dr.row)))

isempty(dr::DomainRange) = length(dr) == 0

first(dr::DomainRange) = first(dr.col) + first(dr.row) - 1
last(dr::DomainRange)  = last(dr.col)  + last(dr.row)  - 1

function in(x::Integer, dr::DomainRange)
    first(dr) <= x <= last(dr) || return false
    ix = x - first(dr.col) + 1
    id = fld(ix, step(dr.col))
    ix -= id * step(dr.col)
    return ix in dr.row
end

shift(dr::DomainRange,   n::Int = 1) = DomainRange(dr.col, shift(dr.row, n))
prepend(dr::DomainRange, n::Int = 1) = DomainRange(dr.col, prepend(dr.row, n))
expand(dr::DomainRange,  n::Int = 1) = DomainRange(dr.col, expand(dr.row, n))
inflate(dr::DomainRange, n::Int = 1) = DomainRange(dr.col, inflate(dr.row, n))

@inline function apply_along_direction(dr::DomainRange, dir::Axis, f, args...)
    if dir == X_axis
        return DomainRange(dr.col, f(dr.row, args...))
    else
        return DomainRange(f(dr.col, args...), dr.row)
    end
end

shift_dir(dr::DomainRange, dir::Axis, n::Int = 1)   = apply_along_direction(dr, dir, shift, n)
prepend_dir(dr::DomainRange, dir::Axis, n::Int = 1) = apply_along_direction(dr, dir, prepend, n)
expand_dir(dr::DomainRange, dir::Axis, n::Int = 1)  = apply_along_direction(dr, dir, expand, n)
inflate_dir(dr::DomainRange, dir::Axis, n::Int = 1) = apply_along_direction(dr, dir, inflate, n)

direction_length(dr::DomainRange, dir::Axis) = dir == X_axis ? length(dr.row) : length(dr.col)

linear_range(dr::DomainRange) = first(dr):last(dr)

show(io::IO, dr::DomainRange) = print(io, "DomainRange{$(dr.col), $(dr.row)}")

"""
    StepsRanges

Holds indexing information for all steps of the solver.
"""
mutable struct StepsRanges
    direction::Axis
    real_domain::DomainRange
    full_domain::DomainRange

    EOS::DomainRange
    fluxes::DomainRange
    cell_update::DomainRange
    advection::DomainRange
    projection::DomainRange

    outer_lb_EOS::DomainRange
    outer_rt_EOS::DomainRange
    outer_lb_fluxes::DomainRange
    outer_rt_fluxes::DomainRange
    inner_EOS::DomainRange
    inner_fluxes::DomainRange
end


function StepsRanges()
    StepsRanges(
        X_axis,
        DomainRange(), DomainRange(),
        DomainRange(), DomainRange(), DomainRange(), DomainRange(), DomainRange(),
        DomainRange(), DomainRange(), DomainRange(), DomainRange(), DomainRange(), DomainRange(),
    )
end
