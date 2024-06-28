
function sort_blocks_by_perimeter_first!(threads_workload, blk_grid, grid_size)
    function is_block_at_thread_perimeter(blk_pos, tid)
        for side in instances(Side.T)  # TODO: replace by `sides_of(length(blk_pos))`
            neighbour_pos = blk_pos + CartesianIndex(offset_to(side))
            !in_grid(neighbour_pos, grid_size) && return true
            blk_grid[neighbour_pos] != tid && return true
        end
        return false
    end

    for (tid, thread_workload) in enumerate(threads_workload)
        sort!(thread_workload; by=blk_pos->!is_block_at_thread_perimeter(blk_pos, tid))  # `false` first
    end

    return threads_workload
end


"""
    simple_workload_distribution(threads, grid_size)

Basic distribution of `B = prod(grid_size)` blocks into groups of `B ÷ threads`, with the remaining
blocks distributed evenly.
"""
function simple_workload_distribution(threads, grid_size)
    block_count = prod(grid_size)
    blocks_per_thread = fld(block_count, threads)
    remaining_blocks = block_count - threads * blocks_per_thread

    # Assign to the n-th thread the `(n:n+1) .* blocks_per_thread` blocks.
    # The first `remaining_blocks` threads have one more block to even out the extra workload.
    threads_workload = map(1:threads) do tid
        prev_tids_blocks = blocks_per_thread * (tid - 1)
        tid_blocks = blocks_per_thread
        if tid > remaining_blocks
            prev_tids_blocks += remaining_blocks
        else
            prev_tids_blocks += tid - 1
            tid_blocks += 1
        end

        return [
            CartesianIndices(grid_size)[blk_idx]
            for blk_idx in (1:tid_blocks) .+ prev_tids_blocks
        ]
    end

    return threads_workload
end


"""
    block_grid_from_workload(grid_size, threads_workload)

Convenience function to convert a `threads_workload` (result of [`thread_workload_distribution`](@ref))
into an `Array` of `grid_size`, with each element assigned to the `tid` given by the distribution.

This makes it easy to visualize the efficiency of the distribution.
"""
function block_grid_from_workload(grid_size, threads_workload)
    blk_grid = zeros(Int, grid_size)
    for (tid, thread_workload) in enumerate(threads_workload)
        blk_grid[thread_workload] .= tid
    end
    return blk_grid
end


function grid_to_scotch_graph(grid_size; weighted=false,
    static_sized_grid=nothing, block_size=nothing, remainder_block_size=nothing, ghosts=0
)
    # `block_size` includes `ghosts`. Therefore for the block weights to include only real cells: `ghosts > 0`

    dim = length(grid_size)
    n_vertices = Scotch.SCOTCH_Num(prod(grid_size))

    # Parse through the regions of the n-cube. This algorithm is just so useful, I love it very much.
    # See this representation: https://en.wikipedia.org/wiki/Binomial_theorem#Geometric_explanation
    n_edges = 0
    grid_regions = ntuple(Returns((false, true)), dim)
    for region in Iterators.product(grid_regions...)
        region_size = ifelse.(region, 1, grid_size .- 1)
        region_order = sum(region)
        # `dim - region_order` is the number of edges which we haven't counted yet for each vertices
        # of this region.
        n_edges += prod(region_size) * (dim - region_order)
    end

    # Compressed storage format of the graph for the grid. Same as for METIS.
    xadj      = Vector{Scotch.SCOTCH_Num}(undef, n_vertices + 1)
    adjncy    = Vector{Scotch.SCOTCH_Num}(undef, 2 * n_edges)
    v_weights = weighted ? Vector{Scotch.SCOTCH_Num}(undef, n_vertices) : nothing

    xadj[1] = 1
    adjncy_idx = 1
    for (vn, v) in enumerate(CartesianIndices(grid_size))
        ne = 0
        for side in instances(Side.T)  # TODO: replace with `sides_of(dim)`
            v_neighbour = v + CartesianIndex(offset_to(side))
            !in_grid(v_neighbour, grid_size) && continue
            i_neighbour = LinearIndices(grid_size)[v_neighbour]
            adjncy[adjncy_idx] = i_neighbour
            adjncy_idx += 1
            ne += 1
        end
        xadj[vn + 1] = xadj[vn] + ne

        if weighted
            # The weight is the number of (real) cells in the block
            v_weights[vn] = prod(block_size_at(v, grid_size, static_sized_grid, block_size, remainder_block_size, ghosts))
        end
    end

    return Scotch.graph_build(xadj, adjncy; v_weights)
end


function partition_cost(parts, grid_size, partition)
    part_tmp = Vector{eltype(partition)}(undef, parts)
    costs = Vector{eltype(partition)}(undef, prod(grid_size))
    for (i, p) in enumerate(partition)
        blk_pos = CartesianIndices(grid_size)[i]
        out_connections = 0
        part_tmp .= 0
        part_tmp_i = 0
        for side in instances(Side.T)
            neighbour_pos = blk_pos + CartesianIndex(offset_to(side))
            !in_grid(neighbour_pos, grid_size) && continue
            neighbour_i = LinearIndices(grid_size)[neighbour_pos]
            neighbour_p = partition[neighbour_i]
            neighbour_p == p && continue
            out_connections += 1
            if !(neighbour_p in (@view part_tmp[1:part_tmp_i]))
                part_tmp_i += 1
                part_tmp[part_tmp_i] = neighbour_p
            end
        end
        # cost = number of neighbours of another part + number of different neighbouring parts
        costs[i] = -(out_connections + part_tmp_i)  # negative as higher cost means harder to move
    end
    return costs
end


"""
    scotch_grid_partition(
        threads, grid_size;
        strategy=:default, workload_tolerance=0, repart=false, retries=10, weighted=false,
        static_sized_grid=nothing, block_size=nothing, remainder_block_size=nothing, ghosts=0
    )

Split `grid_size` to the `threads`.

`strategy` is passed to [`Scotch.strat_flags`](https://keluaa.github.io/Scotch.jl/dev/#Scotch.strat_flags).

`workload_tolerance` is the acceptable workload uneveness among the partitions.
Giving a few blocks of margin (e.g. at least `1/prod(grid_size)`) is preferrable.

`weighted == true` will distribute blocks while taking into account the number of real cells they have.
In this case all parameters of the grid must be present: `static_sized_grid`, `block_size`,
`remainder_block_size` and `ghosts` must be given (obtained with e.g. [`grid_dimensions`](@ref)).

The partitioning is random, hence results may vary. To counterbalance this, giving `retries > 0` will
repeat the partitioning `retries` times and keep the best one.
"""
function scotch_grid_partition(
    threads, grid_size;
    strategy=:default, workload_tolerance=0, repart=false, retries=10, weighted=false,
    static_sized_grid=nothing, block_size=nothing, remainder_block_size=nothing, ghosts=0
)
    graph = grid_to_scotch_graph(grid_size;
        weighted, static_sized_grid, block_size, remainder_block_size, ghosts
    )

    # TODO: with `Scotch.graph_map` it is possible to partition the grid according to a topology of threads, is it better?
    # -> by adding weights to the topology it might be possible to reproduce the processor topology to account for e.g. NUMA or caches
    # -> cores sharing the same L3      => dense graph + weight 1
    # -> cores on the same NUMA         => dense graph + weight 2?
    # -> cores on the same socket       => dense graph + weight 4?
    # -> cores on the different sockets => dense graph + weight 8?

    # TODO: for larger grids, using graph coarsening might be necessary (+ it may help the solver to reach better solutions)

    # TODO: results are random, impose the RNG seed or do something else (repeatedly call the solver N times and keep the best?)
    strat = Scotch.strat_build(:graph_map; strategy, parts=threads, imbalance_ratio=Float64(workload_tolerance))
    partition = Scotch.graph_part(graph, threads, strat)

    if repart
        cost_factor = 1.0
        costs = partition_cost(threads, grid_size, partition)
        partition = Scotch.graph_repart(graph, threads, partition, cost_factor, costs, strat)
    end

    threads_workload = map(1:threads) do tid
        indices = findall(==(tid), partition)
        return CartesianIndices(grid_size)[indices]
    end

    if retries > 0
        # Keep the distribution with the most even workload first, as it is the best way to avoid
        # stalls, and then among those keep the one with the smallest partitions perimeters, to
        # minimize the amount of exchanges between threads.
        block_weights  = graph.v_weights
        best_workload  = threads_workload
        best_eveness   = weighted ? workload_eveness(best_workload, block_weights, grid_size) : workload_eveness(best_workload)
        best_perimeter = total_workload_perimeter(best_workload)
        for _ in 1:retries
            # TODO: reuse the same graph for the retries, as introduces unnecessary overhead
            new_threads_workload = scotch_grid_partition(threads, grid_size; strategy, workload_tolerance, repart, retries=0)
            new_threads_workload == best_workload && continue

            new_eveness = weighted ? workload_eveness(new_threads_workload, block_weights, grid_size) : workload_eveness(new_threads_workload)
            best_eveness < new_eveness && continue

            new_perimeter = total_workload_perimeter(new_threads_workload)
            if (new_eveness < best_eveness) || (new_perimeter < best_perimeter)
                threads_workload = new_threads_workload
                best_workload  = new_threads_workload
                best_eveness   = new_eveness
                best_perimeter = new_perimeter
            end
        end
        return best_workload
    else
        return threads_workload
    end
end


function thread_workload_distribution(params::ArmonParameters; threads=nothing, kwargs...)
    thread_count = @something threads (params.use_threading ? Threads.nthreads() : 1)
    grid_size, static_sized_grid, remainder_block_size = grid_dimensions(params)
    simple = params.workload_distribution === :simple
    scotch = params.workload_distribution in (:scotch, :sorted_scotch, :weighted_sorted_scotch)
    perimeter_first = params.workload_distribution in (:sorted_scotch, :weighted_sorted_scotch)
    merged_kw = merge(params.distrib_params, Dict(kwargs...))
    if params.workload_distribution === :weighted_sorted_scotch
        return thread_workload_distribution(thread_count, grid_size;
            simple, scotch, perimeter_first,
            weighted=true, static_sized_grid, block_size=params.block_size, remainder_block_size, ghosts=params.nghost,
            merged_kw...
        )
    else
        return thread_workload_distribution(thread_count, grid_size; simple, scotch, perimeter_first, merged_kw...)
    end
end


"""
    thread_workload_distribution(params::ArmonParameters; threads=nothing)
    thread_workload_distribution(
        threads::Int, grid_size::Tuple;
        scotch=true, simple=false, perimeter_first=false, kwargs...
    )

Distribute each block in `grid_size` among the `threads`, as evenly as possible.

With `simple == true`, blocks are distributed with [`simple_workload_distribution`](@ref).

With `scotch == true`, the [`Scotch`](https://gitlab.inria.fr/scotch/scotch) graph partitioning solver
is used to better split the grid.
`kwargs` are passed to [`scotch_grid_partition`](@ref).

If `perimeter_first == true`, the resulting distribution will have blocks sorted in way that will
place neighbours of other threads' blocks first in the list.
By doing so, communications between threads may be overlapped more frequently.
"""
function thread_workload_distribution(
    threads::Int, grid_size::Tuple;
    scotch=true, simple=false, perimeter_first=false, kwargs...
)
    if simple
        threads_workload = simple_workload_distribution(threads, grid_size)
    elseif scotch
        threads_workload = scotch_grid_partition(threads, grid_size; kwargs...)
        if perimeter_first
            blk_grid = block_grid_from_workload(grid_size, threads_workload)
            sort_blocks_by_perimeter_first!(threads_workload, blk_grid, grid_size)
        end
    else
        error("unknown workload distribution, expected `simple == true` or `scotch == true`")
    end
    return threads_workload
end


function total_workload_perimeter(threads_workload)
    return map(threads_workload) do thread_workload
        perimeter = 0
        for blk_pos in thread_workload, side in instances(Side.T)
            neighbour_pos = blk_pos + CartesianIndex(offset_to(side))
            perimeter += !(neighbour_pos in thread_workload)
        end
        return perimeter
    end
end


function workload_eveness(threads_workload)
    workloads = length.(threads_workload)
    total_items = sum(workloads)
    threads = length(threads_workload)

    ideal_workload = total_items ÷ threads
    remainder = total_items % threads

    uneveness = 0
    for workload in workloads
        if workload == ideal_workload
            # ok
        elseif workload == ideal_workload + (remainder > 0)
            remainder -= 1
        else
            uneveness += abs(ideal_workload - workload)
        end
    end

    return uneveness + remainder  # 0 if it is optimal
end


function workload_eveness(threads_workload, block_weights, grid_size)
    mean_weight = sum(block_weights) / length(threads_workload)
    uneveness = 0.0
    for thread_workload in threads_workload
        blk_indexes = LinearIndices(grid_size)[thread_workload]
        uneveness += sum((block_weights[blk_indexes] .- mean_weight).^2)
    end
    return uneveness
end
