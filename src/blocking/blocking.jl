
"""
    BlockSize

Dimensions of a [`LocalTaskBlock`](@ref).
"""
abstract type BlockSize end


"""
    StaticBSize{S, Ghost} <: BlockSize

A [`BlockSize`](@ref) of size `S` and `Ghost` cells. `S` is embedded in the type, this reduces the
amount of memory in the parameters of every kernel, as well as allowing the compiler to make some
minor optizations in indexing expressions.

Ghost cells are included in `S`: there are `S .- 2*Ghost` real cells.
"""
struct StaticBSize{S, Ghost} <: BlockSize
    StaticBSize{S, G}() where {S, G} = new{S, G}()
end

StaticBSize(s::Tuple{Vararg{Integer}}, g) = StaticBSize{s, g}()

"`NTuple` of the dimensions of a [`BlockSize`](@ref)"
block_size(::ObjOrType{StaticBSize{S}}) where {S} = S

"Number of ghost cells of a [`BlockSize`](@ref)"
ghosts(::ObjOrType{StaticBSize{S, G}}) where {S, G} = G

"`NTuple` of the dimensions of a [`BlockSize`](@ref), excluding ghost cells"
real_block_size(::ObjOrType{StaticBSize{S, G}}) where {S, G} = S .- 2*G

Base.ndims(::ObjOrType{StaticBSize}) = 2


"""
    DynamicBSize{Ghost} <: BlockSize

Similar to [`StaticBSize`](@ref), but for blocks with a less-than-ideal size: block size is therefore
not stored in the type. This results in less compilation when testing for different domain sizes with
a constant [`StaticBSize`](@ref).

The number of `Ghost` cells is still embedded in the type, as it can simplify some indexing expressions
and for coherency.
"""
struct DynamicBSize{Ghost} <: BlockSize
    s::NTuple{2, UInt16}  # TODO: should it really be UInt16? What about 1D?

    DynamicBSize{G}(s) where {G} = new{G}(s)
end

DynamicBSize(s::Tuple{Vararg{Integer}}, g) = DynamicBSize{g}(Base.convert(Tuple{Vararg{UInt16}}, s))

block_size(bs::DynamicBSize) = bs.s
ghosts(::ObjOrType{DynamicBSize{G}}) where {G} = G
real_block_size(bs::DynamicBSize{G}) where {G} = bs.s .- 2*G
Base.ndims(::ObjOrType{DynamicBSize}) = 2


"""
    block_domain_range(bsize::BlockSize, corners)
    block_domain_range(bsize::BlockSize, bottom_left::Tuple, top_right::Tuple)

A [`DomainRange`](@ref) built from offsets from the corners of `bsize`.

`block_domain_range(bsize, (0, 0), (0, 0))` is the domain of all real cells in the block.
`block_domain_range(bsize, (-g, -g), (g, g))` would be the domain of all cells (real cells + `g`
ghost cells) in the block.
"""
block_domain_range(bsize::BlockSize, corners::NTuple{2, Dims{2}}) = block_domain_range(bsize, corners...)
function block_domain_range(bsize::BlockSize, bottom_left::Dims{2}, top_right::Dims{2})
    ghost = ghosts(bsize)
    row = block_size(bsize)[1]

    block_start = row * (ghost - 1) + ghost
    block_idx(I) = block_start + I[2] * row + I[1]

    first_I = bottom_left .+ (1, 1)
    last_I  = top_right   .+ real_block_size(bsize)

    col_range = block_idx(first_I):row:block_idx(last_I)
    row_range = Base.OneTo(length(first_I[1]:last_I[1]))
    return DomainRange(col_range, row_range)
end


"""
    position(bsize::BlockSize, i)

N-dim position of the `i`-th cell in the block.

If `1 ≤ position(bsize, i)[d] ≤ block_size(bsize)[d]` then the cell is not a ghost cell along the `d`
dimension. See [`is_ghost`](@ref).
"""
position(bsize::BlockSize, i) = position(bsize, block_size(bsize), i)

position(bsize::BlockSize, ::NTuple{1}, i) = i - ghosts(bsize)

function position(bsize::BlockSize, size::NTuple{2}, i::I) where {I}
    row_length = I(size[1])
    iy = (i - one(I)) ÷ row_length
    ix = (i - one(I)) % row_length
    # g=5
    # 6 -> (1, 0)
    return (ix, iy) .- ghosts(bsize) .+ 1
end

function position(bsize::BlockSize, size::NTuple{3}, i::I) where {I}
    row_length   = I(size[1])
    plane_length = I(size[2]) * row_length

    iz      = (i - one(I)) ÷ plane_length
    i_plane = (i - one(I)) % plane_length

    iy = i_plane ÷ row_length
    ix = i_plane % row_length

    return (ix, iy, iz) .- ghosts(bsize) .+ 1
end


"""
    lin_position(bsize::BlockSize, I)

From the `NTuple` (e.g. returned from [`position`](@ref)), return the linear index in the block.
`lin_position(bsize, position(bsize, i)) == i`.
"""
lin_position(bsize::BlockSize, I::NTuple{1}) = I + ghosts(bsize)

function lin_position(bsize::BlockSize, I::NTuple{2})
    return (I[2] + ghosts(bsize) - 1) * block_size(bsize)[1] + (I[1] + ghosts(bsize))
end

function lin_position(bsize::BlockSize, I::NTuple)
    return sum((I .+ (ghosts(bsize) - 1)) .* Base.size_to_strides(1, block_size(bsize)...)) + 1
end


"""
    border_domain(bsize::BlockSize, side::Side)

[`DomainRange`](@ref) of the real cells along `side`. It includes only one "strip" of cells, that is
`length(border_domain(bsize, side)) == size_along(bsize, side)`.
"""
function border_domain(bsize::BlockSize, side::Side)
    rsize = real_block_size(bsize)

    if side == Left
        bl_corner = (0, 0)
        tr_corner = (1 - rsize[1], 0)
    elseif side == Right
        bl_corner = (rsize[1] - 1, 0)
        tr_corner = (0, 0)
    elseif side == Bottom
        bl_corner = (0, 0)
        tr_corner = (0, 1 - rsize[2])
    elseif side == Top
        bl_corner = (0, rsize[2] - 1)
        tr_corner = (0, 0)
    end

    return block_domain_range(bsize, bl_corner, tr_corner)
end


"""
    ghost_domain(params::ArmonParameters, side::Side)

[`DomainRange`](@ref) of all ghosts cells of `side`, excluding the corners of the block.
"""
function ghost_domain(bsize::BlockSize, side::Side)
    g = ghosts(bsize)
    domain = border_domain(bsize, side)
    domain = shift_dir(domain, axis_of(side), side in first_sides() ? -g : g)
    domain = expand_dir(domain, axis_of(side), g - 1)
    return domain
end


stride_along(bsize::BlockSize, axis::Axis)    = Base.size_to_strides(1, block_size(bsize)...)[Int(axis)]
size_along(bsize::BlockSize, axis::Axis)      = block_size(bsize)[Int(axis)]
size_along(bsize::BlockSize, side::Side)      = size_along(bsize, axis_of(side))
real_size_along(bsize::BlockSize, axis::Axis) = real_block_size(bsize)[Int(axis)]
real_size_along(bsize::BlockSize, side::Side) = real_size_along(bsize, axis_of(side))


"""
    is_ghost(bsize::BlockSize, i, o=0)

`true` if the `i`-th cell of the block is a ghost cell, `false` otherwise.

`o` would be a "ring" index: `o == 1` excludes the first ring of ghost cells, etc.
"""
is_ghost(bsize::BlockSize, i, o=0) = !in_grid(1 - o, position(bsize, i), real_block_size(bsize) .+ o)


"""
    in_grid(idx, size)
    in_grid(start, idx, size)

`true` if each axis of `idx` is between `start` and `size`. `start` defaults to `1`.

Argument types can any mix of `Integer`, `Tuple` or `CartesianIndex`.

```julia
julia> in_grid(1, (1, 2), 2)  # same as `in_grid((1, 1), (1, 2), (2, 2))`
true

julia> in_grid((3, 1), (3, 2))
true

julia> in_grid((1, 3), (3, 2))
false
```
"""
in_grid(idx, grid) = in_grid(1, idx, grid)
in_grid(start, idx, grid) = all(Tuple(start) .≤ Tuple(idx) .≤ Tuple(grid))

"""
    in_grid(idx, grid, axis::Axis)
    in_grid(start, idx, grid, axis::Axis)

Same as `in_grid(start, idx, grid)`, but only checks along `axis`.
"""
in_grid(idx, grid, axis::Axis) = in_grid(1, idx, grid, axis)
in_grid(start, idx, grid, axis::Axis) = (Tuple(start) .≤ Tuple(idx) .≤ Tuple(grid))[Int(axis)]


const Neighbours = @NamedTuple{left::T, right::T, bottom::T, top::T} where {T}

function (::Type{Neighbours})(f::Base.Callable, default)
    # Default constructor, e.g: `Neighbours(Int, 0)` => Neighbours{Int}(left=0, right=0, ...)
    values = ntuple(_ -> f(default), fieldcount(Neighbours))
    return Neighbours{eltype(values)}(values)
end


@enumx BlockExchangeState::UInt8 begin
    "There are steps to do before border cells are ready"
    NotReady=0
    "Border cells are ready, waiting for the other side to be `Ready` as well"
    Ready=1
    "One of the blocks is performing the exchange"
    InProgress=2
    "The exchange is done"
    Done=3
end


include("blocks.jl")
include("block_grid.jl")


"""
    @iter_blocks for blk in all_blocks(grid)
        # body...
    end

Applies the body of the for-loop in to all blocks of the `grid`.

The body is duplicated for inner and edge blocks, ensuring type-inference.

If `params.use_multithreading`, then an attempt will be made at equilibrating workload among threads.

```julia
# Iterate on all blocks of the grid
@iter_blocks for blk in all_blocks(grid)
    some_function(blk)
end
```
"""
macro iter_blocks(expr)
    !(expr isa Expr && expr.head === :for) && error("expected for-loop")
    block_var = expr.args[1].args[1]
    block_range = expr.args[1].args[2]
    body = expr.args[2]

    if @capture(block_range, f_(grid_var_))
        if f === :all_blocks || f === all_blocks
            inner_blocks_range = :($grid_var.blocks)
            edge_blocks_range  = :($grid_var.edge_blocks)
        else
            error("expected `all_blocks(grid)`, got: $block_range")
        end
    else
        error("wrong style of block iteration: $block_range")
    end

    return esc(quote
        Armon.@section "Inner blocks" begin
            Armon.@threaded for $block_var in $inner_blocks_range
                $body
            end
        end

        Armon.@section "Edge blocks" begin
            Armon.@threaded for $block_var in $edge_blocks_range
                $body
            end
        end
    end)
end


"""
    BlockRowIterator(grid::BlockGrid; kwargs...)
    BlockRowIterator(grid::BlockGrid, blk::LocalTaskBlock; kwargs...)
    BlockRowIterator(grid::BlockGrid, sub_grid; global_ghosts=false, all_ghosts=false)

Iterate the rows of all blocks of the `grid`, row by row (and not block by block).
This allows to iterate the cells of the `grid` as if it was a single block.

Giving `blk` will return an iterator on the rows of the block.

`sub_grid` defaults to the whole grid: it is a `Tuple` of iterables, one for each axis.

If `global_ghosts == true`, then the ghost cells of at the border of the global domain are also returned.
If `all_ghosts == true`, then the ghost cells of at the border of all blocks are also returned.

```julia
julia> for (blk, row_range) in BlockRowIterator(grid; all_ghosts=true)
           println(blk.pos, " - ", row_range)
       end
(1, 1) - 1:64
(2, 1) - 1:64
(1, 1) - 65:128
(2, 1) - 65:128
(1, 1) - 129:192
(2, 1) - 129:192
...
```
"""
struct BlockRowIterator
    grid::BlockGrid
    row_iter::Iterators.ProductIterator
    global_ghosts::Bool
    all_ghosts::Bool
end


BlockRowIterator(grid::BlockGrid; kwargs...) =
    BlockRowIterator(grid, Base.oneto.(grid.grid_size); kwargs...)

function BlockRowIterator(grid::BlockGrid, blk::LocalTaskBlock; kwargs...)
    # Only iterate the grid at the position of `blk`
    grid_iter = UnitRange.(Tuple(blk.pos), Tuple(blk.pos))
    return BlockRowIterator(grid, grid_iter; kwargs...)
end

function BlockRowIterator(grid::BlockGrid, sub_grid; global_ghosts=false, all_ghosts=false)
    # `sub_grid` is a `Tuple` of iterables
    bs = block_size(static_block_size(grid))
    bs = max.(bs, grid.edge_size)  # Include dimensions of edge blocks in case some are bigger than static blocks
    row_count = (1, bs[2:end]...)

    # Black magic: by intertwining `row_count` and `sub_grid` we create an iterator over all rows
    # and blocks of the grid, in the same order as if we where iterating cell by cell.
    # Edge blocks are handled at iteration time.
    blk_rows_iter = Base.oneto.(row_count)
    row_iter = Iterators.product(Iterators.flatten(zip(blk_rows_iter, sub_grid))...)

    global_ghosts && error("global_ghosts NYI")  # TODO

    return BlockRowIterator(grid, row_iter, global_ghosts, all_ghosts)
end


Base.IteratorSize(::BlockRowIterator) = Base.SizeUnknown()
Base.eltype(::BlockRowIterator) = Tuple{LocalTaskBlock, UnitRange}


function Base.iterate(iter::BlockRowIterator, row_iter_state=0)
    if row_iter_state == 0
        row_iter_state = iterate(iter.row_iter)
    else
        @label next_row
        row_iter_state = iterate(iter.row_iter, row_iter_state)
    end

    row_iter_state === nothing && return nothing
    row_iter_val, row_iter_state = row_iter_state

    # De-intertwine `row_iter_val`
    row_idx = row_iter_val[1:2:end]
    blk_idx = row_iter_val[2:2:end]

    blk = block_at(iter.grid, CartesianIndex(blk_idx))
    blk_size = block_size(blk)
    block_row_count = (1, Base.tail(blk_size)...)

    if blk.size isa DynamicBSize
        # `row_idx` might be outside of `blk`
        !all(row_idx .≤ block_row_count) && @goto next_row
    end

    on_global_edge = false
    if iter.all_ghosts
        # Iterate through all rows
    elseif iter.global_ghosts && !in_grid(2, blk_idx, iter.grid.grid_size .- 1)
        # `blk` is on the global grid edge, is the row also on the edge?
        # TODO
        on_global_edge = true
    else
        # Keep rows with real cells
        if !all(1 .≤ (row_idx[2:end] .- ghosts(blk.size)) .≤ (block_row_count[2:end] .- 2*ghosts(blk.size)))
            @goto next_row
        end
    end

    row_lin_idx = sum(Base.size_to_strides(1, block_row_count...) .* (row_idx .- 1)) + 1
    row_length = blk_size[1]
    row_range = (((row_lin_idx-1) * row_length)+1):(row_lin_idx * row_length)

    if iter.all_ghosts
        # Keep all cells of the row
    elseif iter.global_ghosts && on_global_edge
        # TODO: add ghosts to the range
    else
        # Keep only the real cells of the row
        g = ghosts(blk)
        row_range = (first(row_range) + g):(last(row_range) - g)
    end

    next_element = (blk, row_idx, row_range)
    return next_element, row_iter_state
end
