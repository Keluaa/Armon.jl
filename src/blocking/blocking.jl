
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
    border_domain(bsize::BlockSize, side::Side.T; single_strip=true)

[`DomainRange`](@ref) of the real cells along `side`.

If `single_strip == true`,  it includes only one "strip" of cells, that is
`length(border_domain(bsize, side)) == size_along(bsize, side)`.
Otherwise, there are `ghosts(bsize)` strips of cells: all real cells which would be exchanged with
another block along `side`.
"""
function border_domain(bsize::BlockSize, side::Side.T; single_strip=true)
    rsize = real_block_size(bsize)

    if side == Side.Left
        bl_corner = (0, 0)
        tr_corner = (1 - rsize[1], 0)
    elseif side == Side.Right
        bl_corner = (rsize[1] - 1, 0)
        tr_corner = (0, 0)
    elseif side == Side.Bottom
        bl_corner = (0, 0)
        tr_corner = (0, 1 - rsize[2])
    elseif side == Side.Top
        bl_corner = (0, rsize[2] - 1)
        tr_corner = (0, 0)
    end

    domain = block_domain_range(bsize, bl_corner, tr_corner)
    single_strip && return domain
    if side in first_sides()
        return expand_dir(domain, axis_of(side), ghosts(bsize) - 1)
    else
        return prepend_dir(domain, axis_of(side), ghosts(bsize) - 1)
    end
end


"""
    ghost_domain(bsize::BlockSize, side::Side.T; single_strip=true)

[`DomainRange`](@ref) of all ghosts cells of `side`, excluding the corners of the block.

If `single_strip == true`, then the domain is only 1 cell thick, positionned at the furthest ghost
cell from the real cells.
Otherwise, there are `ghosts(bsize)` strips of cells: all ghost cells which would be exchanged with
another block along `side`.
"""
function ghost_domain(bsize::BlockSize, side::Side.T; single_strip=true)
    domain = border_domain(bsize, side)
    domain = shift_dir(domain, axis_of(side), side in first_sides() ? -ghosts(bsize) : ghosts(bsize))
    single_strip && return domain
    if side in first_sides()
        return expand_dir(domain, axis_of(side), ghosts(bsize) - 1)
    else
        return prepend_dir(domain, axis_of(side), ghosts(bsize) - 1)
    end
end


stride_along(bsize::BlockSize, axis::Axis.T)    = Base.size_to_strides(1, block_size(bsize)...)[Int(axis)]
size_along(bsize::BlockSize, axis::Axis.T)      = block_size(bsize)[Int(axis)]
size_along(bsize::BlockSize, side::Side.T)      = size_along(bsize, axis_of(side))
real_size_along(bsize::BlockSize, axis::Axis.T) = real_block_size(bsize)[Int(axis)]
real_size_along(bsize::BlockSize, side::Side.T) = real_size_along(bsize, axis_of(side))

face_size(bsize::BlockSize, side::Axis.T)      = prod(block_size(bsize)) ÷ size_along(bsize, side)
face_size(bsize::BlockSize, side::Side.T)      = face_size(bsize, axis_of(side))
real_face_size(bsize::BlockSize, side::Axis.T) = prod(real_block_size(bsize)) ÷ real_size_along(bsize, side)
real_face_size(bsize::BlockSize, side::Side.T) = real_face_size(bsize, axis_of(side))


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
    in_grid(idx, grid, axis::Axis.T)
    in_grid(start, idx, grid, axis::Axis.T)

Same as `in_grid(start, idx, grid)`, but only checks along `axis`.
"""
in_grid(idx, grid, axis::Axis.T) = in_grid(1, idx, grid, axis)
in_grid(start, idx, grid, axis::Axis.T) = (Tuple(start) .≤ Tuple(idx) .≤ Tuple(grid))[Int(axis)]


const Neighbours = @NamedTuple{left::T, right::T, bottom::T, top::T} where {T}

function (::Type{Neighbours})(f::Base.Callable, default)
    # Default constructor, e.g: `Neighbours(Int, 0)` => Neighbours{Int}(left=0, right=0, ...)
    values = ntuple(_ -> f(default), fieldcount(Neighbours))
    return Neighbours{eltype(values)}(values)
end


"""
    BlockInterface

Represents the interface between two neighbouring [`TaskBlock`](@ref)s.
It synchronizes the state of the blocks to make sure the halo exchange happens when both blocks are
ready, and that the operation is done by only one of the blocks.
"""
mutable struct BlockInterface
    """
    0b00XX: ready flag, 0bXX00: exchange state

    Ready flag:
     - `0b00`: no side ready
     - `0b10`: left side ready
     - `0b01`: right side ready
     - `0b11`: both sides ready
    """
    flags :: Atomic{UInt8}
    # Non-atomics to be used only by their respective blocks, to avoid repeating the same exchange
    # multiple times in the same step.
    is_left_side_done :: Bool
    is_right_side_done :: Bool

    BlockInterface() = new(Atomic{UInt8}(0), false, false)
end


include("blocks.jl")
include("block_grid.jl")
include("interface.jl")
include("block_tree.jl")


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
    BlockRowIterator(grid::BlockGrid, blk₁_pos, blk₂_pos; kwargs...)
    BlockRowIterator(grid::BlockGrid, sub_domain::NTuple{2, CartesianIndex}; kwargs...)
    BlockRowIterator(grid::BlockGrid, row_iter::CartesianIndices; global_ghosts=false, all_ghosts=false)

Iterate the rows of all blocks of the `grid`, row by row (and not block by block).
This allows to iterate the cells of the `grid` as if it was a single block.

Elements are tuples of `(block, global_row_idx, row_range)`. `row_range` is the range of cells in
`block` for the current row.

Giving `blk` will return an iterator on the rows of the block.

Giving `blk₁_pos` and `blk₂_pos` will return an iterator over all rows between those blocks.

Giving `sub_domain` will return an iterator including only the cells contained in `sub_domain`.
`sub_domain` is a cuboid defined by the position of the cells in the whole domain of `grid`, using
its lower and upper corners.

`row_iter` is a iterator over global row indices.

If `global_ghosts == true`, then the ghost cells of at the border of the global domain are also returned.
If `all_ghosts == true`, then the ghost cells of at the border of all blocks are also returned.

```jldoctest
julia> params = ArmonParameters(; N=(24, 8), nghost=4, block_size=(20, 12), use_MPI=false);

julia> grid = BlockGrid(params);

julia> for (blk, row_idx, row_range) in Armon.BlockRowIterator(grid)
           println(Tuple(blk.pos), " - ", row_idx, " - ", row_range)
       end
(1, 1) - (1, 1) - 85:96
(2, 1) - (2, 1) - 85:96
(1, 1) - (1, 2) - 105:116
(2, 1) - (2, 2) - 105:116
(1, 1) - (1, 3) - 125:136
(2, 1) - (2, 3) - 125:136
(1, 1) - (1, 4) - 145:156
(2, 1) - (2, 4) - 145:156
(1, 2) - (1, 5) - 85:96
(2, 2) - (2, 5) - 85:96
(1, 2) - (1, 6) - 105:116
(2, 2) - (2, 6) - 105:116
(1, 2) - (1, 7) - 125:136
(2, 2) - (2, 7) - 125:136
(1, 2) - (1, 8) - 145:156
(2, 2) - (2, 8) - 145:156
```
"""
struct BlockRowIterator  # TODO: use have Dim as type param + use NTuple 
    grid          :: BlockGrid
    row_iter      :: CartesianIndices
    rows_per_blk  :: CartesianIndex
    global_ghosts :: Bool
    all_ghosts    :: Bool
end


function row_range_from_corners(grid, bl_blk_pos, tr_blk_pos, include_ghosts, last_offset=nothing)
    bs = include_ghosts ? static_block_size(grid) : real_block_size(grid)
    static_row_count = (1, bs[2:end]...)

    first_row_pos = block_origin(grid, bl_blk_pos, include_ghosts)
    first_row_pos = ((first_row_pos[1] - 1) ÷ ifelse(bs[1] == 0, 1, bs[1]) + 1, first_row_pos[2:end]...)

    last_row_pos = block_origin(grid, tr_blk_pos, include_ghosts)
    last_row_pos = ((last_row_pos[1] - 1) ÷ ifelse(bs[1] == 0, 1, bs[1]) + 1, last_row_pos[2:end]...)

    if last_offset === nothing
        # The default is to iterate over all rows, so we want the number of rows in `tr_blk_pos`
        edge_size = (1, (grid.edge_size[2:end] .+ (include_ghosts ? 2*ghosts(grid) : 0))...)
        last_offset = ifelse.(
            in_grid.(Ref(tr_blk_pos), Ref(grid.static_sized_grid), instances(Axis.T)),
            static_row_count,
            edge_size
        )
    end
    last_row_pos = last_row_pos .+ last_offset .- 1

    return first_row_pos, last_row_pos
end


BlockRowIterator(grid::BlockGrid; kwargs...) =
    BlockRowIterator(grid, 1, grid.grid_size; kwargs...)

BlockRowIterator(grid::BlockGrid, blk::LocalTaskBlock; kwargs...) =
    BlockRowIterator(grid, blk.pos, blk.pos; kwargs...)

function BlockRowIterator(grid::BlockGrid, bl_blk_pos, tr_blk_pos; global_ghosts=false, all_ghosts=false)
    first_row_pos, last_row_pos = row_range_from_corners(
        grid, Tuple(bl_blk_pos), Tuple(tr_blk_pos), global_ghosts || all_ghosts
    )
    row_iter = CartesianIndex(first_row_pos):CartesianIndex(last_row_pos)
    return BlockRowIterator(grid, row_iter; global_ghosts, all_ghosts)
end

function BlockRowIterator(grid::BlockGrid, sub_domain::NTuple{2, CartesianIndex}; global_ghosts=false, all_ghosts=false)
    if global_ghosts || all_ghosts
        sub_domain_with_ghosts = sub_domain
        # Build `sub_domain` solely to find which blocks we need
        sub_domain_bl = CartesianIndex(clamp.(Tuple(sub_domain[1]), 1, grid.cell_size))
        sub_domain_tr = CartesianIndex(clamp.(Tuple(sub_domain[2]), 1, grid.cell_size))
        sub_domain = (sub_domain_bl, sub_domain_tr)
    end

    # Only iterate the grid for the rows included in `sub_domain` (in real cells)
    bl_blk_pos, bl_cell_pos = block_pos_containing_cell(grid, sub_domain[1])
    tr_blk_pos, tr_cell_pos = block_pos_containing_cell(grid, sub_domain[2])

    if global_ghosts || all_ghosts
        # Shift the original `sub_domain` to global indexing, relative to the origin of the block
        g = ghosts(grid)
        bl_cell_pos = Tuple(sub_domain_with_ghosts[1]) .- (Tuple(bl_blk_pos) .- 1) .* real_block_size(grid) .+ g
        tr_cell_pos = Tuple(sub_domain_with_ghosts[2]) .- (Tuple(tr_blk_pos) .- 1) .* real_block_size(grid) .+ g
    else
        bl_cell_pos = Tuple(bl_cell_pos)
        tr_cell_pos = Tuple(tr_cell_pos)
    end

    bl_row_pos = (1, bl_cell_pos[2:end]...)
    tr_row_pos = (1, tr_cell_pos[2:end]...)

    first_row_pos, last_row_pos = row_range_from_corners(
        grid, Tuple(bl_blk_pos), Tuple(tr_blk_pos), global_ghosts || all_ghosts, tr_row_pos
    )
    first_row_pos = first_row_pos .+ bl_row_pos .- 1

    row_iter = CartesianIndex(first_row_pos):CartesianIndex(last_row_pos)
    return BlockRowIterator(grid, row_iter; global_ghosts, all_ghosts)
end

function BlockRowIterator(grid::BlockGrid, row_iter::CartesianIndices; global_ghosts=false, all_ghosts=false)
    if global_ghosts || all_ghosts
        row_count = (1, static_block_size(grid)[2:end]...)
    else
        row_count = (1, real_block_size(grid)[2:end]...)
    end
    row_count = ifelse.(row_count .≤ 0, grid.edge_size, row_count)
    return BlockRowIterator(grid, row_iter, CartesianIndex(row_count), global_ghosts, all_ghosts)
end


Base.IteratorSize(::Type{BlockRowIterator}) = Base.SizeUnknown()
Base.eltype(::Type{BlockRowIterator}) = Tuple{LocalTaskBlock, NTuple{2, Int}, UnitRange}


function Base.iterate(iter::BlockRowIterator, row_iter_state=0)
    if row_iter_state == 0
        row_iter_state = iterate(iter.row_iter)
    else
        @label next_row
        row_iter_state = iterate(iter.row_iter, row_iter_state)
    end

    row_iter_state === nothing && return nothing
    global_row_idx, row_iter_state = row_iter_state
    global_row_idx = Tuple(global_row_idx)

    # Block position from the global row index
    blk_pos = (global_row_idx .- 1) .÷ Tuple(iter.rows_per_blk) .+ 1
    if !in_grid(blk_pos, iter.grid.static_sized_grid)
        # In the edge blocks
        blk_pos = clamp.(blk_pos, 1, iter.grid.grid_size)
    end

    blk = block_at(iter.grid, CartesianIndex(blk_pos))
    blk_size = block_size(blk)
    block_row_count = (1, Base.tail(blk_size)...)
    g = ghosts(blk)

    # Local row index
    row_idx = global_row_idx .- (blk_pos .- 1) .* Tuple(iter.rows_per_blk)
    if !(iter.all_ghosts || iter.global_ghosts)
        row_idx = (row_idx[1], (row_idx[2:end] .+ g)...)
    end

    on_global_edge = false
    left_ghosts = false
    right_ghosts = false
    if iter.global_ghosts
        if in_grid(2, blk_pos, iter.grid.grid_size .- 1)
            # Keep rows with real cells
            if !all(g .< row_idx[2:end] .≤ (block_row_count[2:end] .- g))
                @goto next_row
            end
        else
            # `blk` is on the global grid edge. Keep only rows containing real cells or ghost cells on the edge
            on_global_edge = true

            # Include all rows containing real cells
            include_row = all(g .< row_idx[2:end] .≤ (block_row_count[2:end] .- g))

            # Include all rows with ghosts along the global edge
            include_row |= all(blk_pos[2:end] .== 1 .&& row_idx[2:end] .≤ block_row_count[2:end] .- g)
            include_row |= all(blk_pos[2:end] .== iter.grid.grid_size[2:end] .&& g .< row_idx[2:end])

            !include_row && @goto next_row
        end
    end

    row_lin_idx = sum(Base.size_to_strides(1, block_row_count...) .* (row_idx .- 1)) + 1
    row_length = blk_size[1]
    row_range = (((row_lin_idx-1) * row_length)+1):(row_lin_idx * row_length)

    if iter.all_ghosts
        # Keep all cells of the row
    elseif iter.global_ghosts && on_global_edge
        # Keep only ghost cells on the global edge
        left_ghosts = blk_pos[1] == 1
        right_ghosts = blk_pos[1] == iter.grid.grid_size[1]
        row_range = (first(row_range) + g * !left_ghosts):(last(row_range) - g * !right_ghosts)
    else
        # Keep only the real cells of the row
        row_range = (first(row_range) + g):(last(row_range) - g)
    end

    next_element = (blk, global_row_idx, row_range)
    return next_element, row_iter_state
end
