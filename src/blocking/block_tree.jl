
"""
    BlockTree(grid::BlockGrid, branches_per_level::Vector{Int})

A tree-like structure containing `sub_blocks` (other `BlockTree`s) as well as blocks
([`LocalTaskBlock`](@ref)s). It contains `branches_per_level[1]` branches, themselves containing
`branches_per_level[2]` branches each, etc...
The last level contains `branches_per_level[end]` [`LocalTaskBlock`](@ref)s, and no branches.

`depth(bt) == 0` is the root node.

Technically, there is no topology constraint to which block may be part of the tree: for now blocks
are naively assigned during the construction of the tree.

`branches_per_level` might specify a larger amount of blocks than there is in `grid`. In this case
the extra blocks are not added to the tree, and empty branches are pruned.
"""
struct BlockTree{D, H, BS <: BlockSize, SState <: SolverState, Ghost} <: TaskBlock{D}
    sub_blocks  :: Vector{BlockTree{D, H, BS, SState, Ghost}}
    blocks      :: Vector{LocalTaskBlock{D, H, BS, SState}}
    edge_blocks :: Vector{LocalTaskBlock{D, H, DynamicBSize{Ghost}, SState}}
    depth       :: Int   # 0: root meta-block

    function BlockTree{D, H, BS, SState, Ghost}(depth::Int) where {D, H, BS, SState, Ghost}
        return new{D, H, BS, SState, Ghost}([], [], [], depth)
    end
end


function BlockTree(grid::BlockGrid, branches_per_level::Vector{Int})
    if prod(branches_per_level) < prod(grid.grid_size)
        solver_error(:config, "BlockTree does not have enough sub blocks per level \
                               ($branches_per_level, total: $(prod(branches_per_level))) \
                               to cover the whole grid $(grid.grid_size), total: $(prod(grid.grid_size))")
    end

    bt = BlockTree(grid, branches_per_level, (), 1)

    # Basic sanity check
    if tree_block_count(bt) ≠ prod(grid.grid_size)
        # solver_error(:config, "failed to create a BlockTree, block count mismatch: $(tree_block_count(bt)) ≠ $(prod(grid.grid_size))")
        @warn "failed to create a BlockTree, block count mismatch: $(tree_block_count(bt)) ≠ $(prod(grid.grid_size))"
    end

    return bt
end


function BlockTree(grid::BlockGrid{<:Any, D, H, <:Any, Ghost, BS, SState, <:Any},
    branches_per_level::Vector{Int}, branch_pos::Tuple{Vararg{Int}}, depth::Int
) where {D, H, Ghost, BS, SState}
    bt = BlockTree{D, H, BS, SState, Ghost}(depth - 1)
    sub_blocks_count = branches_per_level[depth]

    if depth < length(branches_per_level)
        resize!(bt.sub_blocks, sub_blocks_count)
        map!(bt.sub_blocks, 1:sub_blocks_count) do i
            BlockTree(grid, branches_per_level, (branch_pos..., i), depth+1)
        end

        # Prune empty branches
        filter!(bt.sub_blocks) do b
            !isempty(b.blocks) || !isempty(b.edge_blocks) || !isempty(b.sub_blocks)
        end
    else
        # Only the last level contains normal blocks
        strides = Base.size_to_strides(reverse(branches_per_level)...)
        global_pos = sum(strides .* reverse(branch_pos .- 1))

        first_blk_i = global_pos + 1
        last_blk_i  = min(first_blk_i + sub_blocks_count - 1, prod(grid.grid_size))
        blk_i_range = first_blk_i:last_blk_i

        sizehint!(bt.blocks, length(blk_i_range))
        for blk_i in blk_i_range
            blk_pos = CartesianIndices(grid.grid_size)[blk_i]
            if in_grid(blk_pos, grid.static_sized_grid)
                blk = grid.blocks[block_idx(grid, blk_pos)]
                push!(bt.blocks, blk)
            else
                blk = grid.edge_blocks[edge_block_idx(grid, blk_pos)]
                push!(bt.edge_blocks, blk)
            end
        end
        resize!(bt.blocks, length(bt.blocks))
    end

    return bt
end


"""
    all_blocks(bt::BlockTree)

Iterator over every [`LocalTaskBlock`](@ref)s contained in `bt` (recursive).
"""
all_blocks(bt::BlockTree) = Iterators.flatten((bt.blocks, bt.edge_blocks, all_blocks.(bt.sub_blocks)...))


"""
    is_leaf(bt::BlockTree)

`true` if `bt` has no children, and contains only [`LocalTaskBlock`](@ref)s.
"""
is_leaf(bt::BlockTree) = isempty(bt.sub_blocks)


"""
    depth(bt::BlockTree)

Depth of `bt` in the parent tree.
"""
depth(bt::BlockTree) = bt.depth


"""
    tree_block_count(bt::BlockTree)

Number of [`LocalTaskBlock`](@ref)s contained in `bt` (recursive).
"""
function tree_block_count(bt::BlockTree)
    # TODO: why does `mapreduce(tree_block_count, +, bt.sub_blocks; init=blk_count)` allocates so much,
    # and why `sum(tree_block_count.(bt.sub_blocks); init=0)` creates systematic compilation time?
    blk_count = length(bt.blocks) + length(bt.edge_blocks)
    for sb in bt.sub_blocks
        blk_count += tree_block_count(sb)
    end
    return blk_count
end


static_block_size(::ObjOrType{BlockTree{D, H, BS}}) where {D, H, BS} = block_size(BS)
real_block_size(::ObjOrType{BlockTree{D, H, BS}}) where {D, H, BS} = real_block_size(BS)
ghosts(::ObjOrType{BlockTree{D, H, BS, SS, Ghost}}) where {D, H, BS, SS, Ghost} = Ghost


"""
    visit_all_tree_block(f, bt::BlockTree; topdown=true, min_depth=0, max_depth=typemax(Int))

Recursively calls `f(pos::Tuple, mb_child::BlockTree)` for `bt` and its children.

If `topdown == true`, then `f` is called on parents before their children.

`f` is called only for [`BlockTree`](@ref)s with a depth `≥ min_depth` and `≤ max_depth`, defaulting
to the whole `BlockTree`.
Children with a depth `> max_depth` will not be recursed into.
"""
function visit_all_tree_block(f, bt::BlockTree, pos=(); topdown=true, min_depth=0, max_depth=typemax(Int))
    topdown && bt.depth ≥ min_depth && f(pos, bt)
    bt.depth < max_depth && for (i, sb) in enumerate(bt.sub_blocks)
        visit_all_tree_block(f, sb, (pos..., i); topdown, min_depth, max_depth)
    end
    !topdown && bt.depth ≥ min_depth && f(pos, bt)
    return
end


"""
    iter_tree_block(f, bt::BlockTree)

Calls `f(blk::LocalTaskBlock)` for all [`LocalTaskBlock`](@ref)s in `bt`.
"""
function iter_tree_block(f, bt::BlockTree)
    foreach(f, bt.blocks)
    foreach(f, bt.edge_blocks)
    foreach(sub_blk -> iter_tree_block(f, sub_blk), bt.sub_blocks)
end


reset!(bt::BlockTree) = iter_tree_block(reset!, bt)

device_to_host!(::BlockTree{H, H}) where {H} = nothing
host_to_device!(::BlockTree{H, H}) where {H} = nothing

device_to_host!(bt::BlockTree{D, H}) where {D, H} = iter_tree_block(device_to_host!, bt)
host_to_device!(bt::BlockTree{D, H}) where {D, H} = iter_tree_block(host_to_device!, bt)


"""
    apply_until_true(f, bt::BlockTree, sub_blocks_range=eachindex(bt.sub_blocks), repeats=2)

Applies recursively `f` on all [`LocalTaskBlock`](@ref)s of `bt`. As long as `f` returns `false`,
it will be repeatedly called up to `repeats^max_depth` times for the block (`max_depth` being the
depth of the `LocalTaskBlock` relative to `bt`).
"""
function apply_until_true(f, bt::BlockTree, sub_blocks_range=eachindex(bt.sub_blocks), repeats=2)
    reps = 1
    val = false
    if is_leaf(bt)
        while !val && reps ≤ repeats
            val = mapreduce(f, &, bt.blocks; init=true)
            val = mapreduce(f, &, bt.edge_blocks; init=val)
            reps += 1
        end
    else
        states = falses(length(sub_blocks_range))
        while !val && reps ≤ repeats
            val = true
            for (i, sub_blk_i) in enumerate(sub_blocks_range)
                states[i] && continue  # Avoid excessive exponential recursion by skipping sub blocks which are already done
                val &= states[i] = apply_until_true(f, bt.sub_blocks[sub_blk_i])
            end
            reps += 1
            val = val || all(states)
        end
    end
    return val
end


function Base.show(io::IO, @nospecialize(bt::BlockTree))
    Base.show(io, typeof(bt))
    if is_leaf(bt)
        print(io, "(blocks: ", length(bt.blocks), ", edge blocks: ", length(bt.edge_blocks), ")")
    else
        print(io, "(sub blocks: ", length(bt.sub_blocks), ")")
    end
end


function Base.show(io::IO, ::MIME"text/plain", @nospecialize(bt::BlockTree))
    print(io, "BlockTree:")
    visit_all_tree_block(bt) do pos, sub_mb
        print(io, '\n', "  "^sub_mb.depth, pos, " - ")
        if is_leaf(sub_mb)
            total_blocks = length(sub_mb.blocks) + length(sub_mb.edge_blocks)
            print(io, total_blocks, " block", total_blocks != 1 ? "s" : "", ", at [")

            max_in_line = 8
            blks_pos = map(blk -> Tuple(blk.pos), Iterators.take(all_blocks(sub_mb), max_in_line))

            if total_blocks > max_in_line
                foreach(blk_pos -> print(io, blk_pos, ", "), blks_pos)
                print(io, "... ")
                blks_pos = map(blk -> Tuple(blk.pos), Iterators.drop(all_blocks(sub_mb), max(max_in_line, total_blocks - max_in_line)))
            end

            foreach(blk_pos -> print(io, blk_pos, ", "), blks_pos[1:end-1])
            print(io, blks_pos[end], ']')
        else
            sub_count = length(sub_mb.sub_blocks)
            print(io, sub_count, sub_count == 1 ? " child" : " children")
        end
    end
end


"""
    block_levels_for_machine(params::ArmonParameters; include_ghosts=false, kwargs...)
    block_levels_for_machine(
        data_type::Type, block_count::Int, block_size::Tuple, vars_per_block::Int;
        levels=[:L3_cache, :L2_cache, :L1_cache]
    )

Attempts to compute efficient sizes for branches of a [`BlockTree`](@ref) for the current machine
using its topology and the cache sizes.
"""
function block_levels_for_machine(params::ArmonParameters; include_ghosts=false, kwargs...)
    block_size = params.block_size
    include_ghosts && (block_size = block_size .- 2*params.nghost)
    grid_size, _, _ = grid_dimensions(params)
    var_count = length(main_vars())  # Not all variables, as they are never all used in a single kernel
    return block_levels_for_machine(data_type(params), prod(grid_size), block_size, var_count; kwargs...)
end


function block_levels_for_machine(
    data_type::Type, block_count::Int, block_size::Tuple, vars_per_block::Int;
    levels=[:L3_cache, :L2_cache, :L1_cache]
)
    data_size = sizeof(data_type)
    useful_block_size = prod(block_size) * vars_per_block * data_size  # In bytes

    all_levels = [:socket, :numa, :L3_cache, :L2_cache, :L1_cache]
    wrong_levels = setdiff(levels, all_levels)
    !isempty(wrong_levels) && error("unknown topology levels: $wrong_levels")

    # TODO: more generic approach to handle any topology?
    topo_info = Hwloc.getinfo(; list_all=true)
    levels_per_prev = Int[]
    blocks_per_level = Int[]
    for level in all_levels
        prev_count = prod(levels_per_prev; init=1)
        entity_count = if level === :socket
            topo_info[:Package]
        elseif level === :numa
            topo_info[:NUMANode]
        elseif level === :L3_cache
            topo_info[:L3Cache]
        elseif level === :L2_cache
            topo_info[:L2Cache]
        elseif level === :L1_cache
            topo_info[:L1Cache]
        else
            error("unknown topology level: $level")
        end
        entity_count == 0 && error("Hwloc found 0 $level")
        push!(levels_per_prev, entity_count ÷ prev_count)
        level in levels && push!(blocks_per_level, last(levels_per_prev))
    end

    flc_size = if last(levels) === :L3_cache
        Hwloc.cachesize(:L3)
    elseif last(levels) === :L2_cache
        Hwloc.cachesize(:L2)
    elseif last(levels) === :L1_cache
        Hwloc.cachesize(:L1)
    else
        error("deepest level must be a cache, got: $(last(levels))")
    end
    blocks_per_level[end] *= cld(flc_size, useful_block_size)

    # Grow the tree to guarrentee there will be enough blocks to cover the entire grid
    main_branches = cld(block_count, prod(blocks_per_level))
    pushfirst!(blocks_per_level, main_branches)

    return blocks_per_level
end
