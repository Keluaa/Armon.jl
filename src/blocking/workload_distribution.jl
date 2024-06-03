
# From https://stackoverflow.com/a/40700741
function largest_2_factors(n)
    # TODO: 3D+
    # n = a*b  =>  a ≤ √n ≤ b  =>  choose (a, b) such that a - b is minimized
    a = round(Int, sqrt(n))
    while n % a > 0
        a -= 1
    end
    b = n ÷ a
    return a, b
end


function surface_perimeter_ratio(g)
    # TODO: 3D+
    surface = prod(g)
    perimeter = sum(g) .* 2
    return surface / perimeter
end


function smoothest_square_subdivision(grid_size, k; optimize=true)
    # F = largest_2_factors(k)  # TODO: remove?
    N = prod(grid_size)
    # r = N % k  # extra elements
    # G = (N - r) ÷ k  # group size

    r = 0
    G = (N + (k - N % k)) ÷ k  # group size

    # Find the group dimensions with the largest group size and biggest surface/perimeter ratio.
    # Starting from an initial group size, we only decrease it to keep a minimum of `k` groups.
    best_G = G
    best_Gd = init_Gd = largest_2_factors(G)
    best_Gd_score = surface_perimeter_ratio(best_Gd)
    for Gi in G-1:-1:max(1, G-k, optimize ? 0 : G)
        Gd = largest_2_factors(Gi)
        Gd_score = surface_perimeter_ratio(Gd)
        if Gd_score > best_Gd_score  # No ties to break: we would always keep the bigger `best_G`, the previous one
            best_G = Gi
            best_Gd = Gd
            best_Gd_score = Gd_score
        end
    end
    r += (prod(init_Gd) - prod(best_Gd)) * k

    best_Gd = Tuple(sort(collect(best_Gd); rev=true))  # Place the biggest axes first
    F = cld.(grid_size, best_Gd)
    prod(F) < k && (F = F .+ 1)

    return F, best_G, best_Gd, r
end


function expanding_chocolate_chips_cake(grid_size::Dims{D}, k::Int) where {D}
    ideal_size = ceil(Int, k^(1/D))
    ideal_size = min(ideal_size, minimum(grid_size))

    expansion_factor = minimum(grid_size .÷ ideal_size)  # at least 1
    group_grid = ntuple(Returns(ideal_size), D)
    group_size = ntuple(Returns(expansion_factor), D)

    changed = true
    while changed && any(group_grid .* group_size .< grid_size)
        changed = false
        # Continue expanding in each direction until we fill the `grid_size`
        for axis in instances(Axis.T)  # TODO: replace by `axis_of(D)`
            new_group_size = group_size .+ offset_to(axis)
            if all(group_grid .* new_group_size .< grid_size .+ new_group_size)
                group_size = new_group_size
                changed = true
            end
        end
    end

    return group_grid, group_size
end


function square_block_distribution(threads, grid_size::Dims{D}) where {D}
    # group_grid, _, group_size, _ = smoothest_square_subdivision(grid_size, threads)
    # @show F, G, Gd, r
    # TODO: fix this, but how?
    # any(group_size .> grid_size) && error("invalid group size with $threads groups: $group_size vs grid size: $grid_size")

    group_grid, group_size = expanding_chocolate_chips_cake(grid_size, threads)

    threads_workload = Vector{Vector{CartesianIndex{D}}}(undef, threads)
    groups_pos = CartesianIndices(group_grid)
    blk_grid = zeros(Int, grid_size)
    for tid in 1:threads
        group_pos = groups_pos[tid]
        thread_work = CartesianIndex{D}[]
        sizehint!(thread_work, prod(group_size))
        threads_workload[tid] = thread_work

        # Trivial distribution of the first `G` elements
        group_origin = CartesianIndex(group_size .* (Tuple(group_pos) .- 1))
        for blk_group_pos in CartesianIndices(group_size)
            blk_pos = group_origin + blk_group_pos
            !in_grid(blk_pos, grid_size) && continue  # TODO: debug??
            push!(thread_work, blk_pos)
            blk_grid[blk_pos] = tid
        end

        isempty(thread_work) && error("$threads $grid_size led ")
    end

    # Parse through the grid, assigning the unassigned blocks to nearby groups.
    # Since we do `group_origin = group_size .* group_pos`, then they must be at the edges.
    # R_size = grid_size .% group_size
    # Magic code similar to `EdgeBlockRegions` to iterate over edges of a n-grid
    # edge_iter = Iterators.product(ifelse.(R_size .> 0, Ref((false, true)), Ref((false,)))...)
    # edge_iter = Iterators.drop(edge_iter, 1)  # drop the all `false` element
    # for edge_axes in edge_iter
        # `R_size` controls the 'width' of the edge
        # first_pos = CartesianIndex(ifelse.(edge_axes, grid_size .- R_size .+ 1, 1))
        # last_pos  = CartesianIndex(ifelse.(edge_axes, grid_size, grid_size .- R_size))
        # for blk_pos in first_pos:last_pos
        for blk_pos in CartesianIndices(grid_size)
            blk_grid[blk_pos] != 0 && continue  # already assigned
            # `first_pos` is the corner closest to the grid's origin, therefore it always has one
            # assigned neighbouring block. Hence we can be sure that at any stage in the algorithm,
            # we would have one assigned block in the neighbours of `blk_pos`.
            # Here we choose the neighbouring group with the least amount of assigned blocks.
            min_workload = (typemax(Int), 0)
            for side in instances(Side.T)  # TODO: replace by `sides_of(D)`
                neighbour_pos = blk_pos + CartesianIndex(offset_to(side))
                !in_grid(neighbour_pos, grid_size) && continue
                tid = blk_grid[neighbour_pos]
                tid == 0 && continue
                min_workload = min(min_workload, (length(threads_workload[tid]), tid))
            end
            _, tid = min_workload
            tid == 0 && error("at `$(blk_pos)`, no assigned neighbours")
            push!(threads_workload[tid], blk_pos)
            blk_grid[blk_pos] = tid
        end
    # end

    return blk_grid, threads_workload
end


function diffuse_workload!(
    threads_workload::Vector{Vector{CartesianIndex{D}}}, blk_grid::Array{Int, D}, grid_size::Dims{D};
    workload_tolerance=1.01, max_iter=1000
) where {D}
    # All threads should handle a maximum of `min_workload*workload_tolerance` blocks.
    mean_workload = round(prod(grid_size) / length(threads_workload))  # Ideal expected workload
    workload_tolerance = max(1 + 1 / mean_workload, workload_tolerance)  # 5% or 1 whole additional block

    function move_to_group(blk::CartesianIndex, prev_tid, new_tid)
        deleteat!(threads_workload[prev_tid], findfirst(==(blk), threads_workload[prev_tid]))
        push!(threads_workload[new_tid], blk)
        blk_grid[blk] = new_tid
    end

    function move_cost(blk::CartesianIndex, prev_tid, new_tid)
        prev_tid == new_tid && return 0
    
        # The cost is the added (`D-1` dimensional) perimeter length to `new_tid`.
        perimeter_increment = 0
        for side in instances(Side.T)  # TODO: replace by `sides_of(D)`
            neighbour_pos = blk + CartesianIndex(offset_to(side))
            # !in_grid(neighbour_pos, grid_size) && continue
            # perimeter_increment += blk_grid[neighbour_pos] != new_tid ? 1 : -1

            # TODO: this seems to improve results
            if in_grid(neighbour_pos, grid_size)
                perimeter_increment += blk_grid[neighbour_pos] != new_tid ? 1 : -1
            else
                perimeter_increment += 1
            end
        end

        workload_increment = 0  # TODO: get the weight of what we move, e.g. 1 for a static block, something else for an edge block
        return perimeter_increment + workload_increment
    end

    # "Diffuse" the workload by moving blocks from one thread to another, until the condition is met for all threads.
    # Candidates are sorted by ascending target tid workload then ascending move cost.
    move_candidates = Vector{Tuple{Int, Int, Int, Int, CartesianIndex{D}}}()
    performed_changes = true
    iter = 0
    while performed_changes && iter < max_iter
        performed_changes = false
        iter += 1
        min_workload = minimum(length.(threads_workload))
        ideal_workload = ceil(Int, min_workload * workload_tolerance)

        # TODO: debug
        # println("After $(iter-1) iterations (min: $min_workload, max: $(maximum(length.(threads_workload)))):")
        # Main.disp_block_distrib(grid_size, threads_workload)
        # println()

        # Build a list of possible blocks to move to neighbouring threads
        empty!(move_candidates)
        for (tid, thread_workload) in enumerate(threads_workload)
            length(thread_workload) ≤ ideal_workload && continue

            for blk_pos in thread_workload
                for side in instances(Side.T)  # TODO: replace by `sides_of(D)`
                    neighbour_pos = blk_pos + CartesianIndex(offset_to(side))
                    !in_grid(neighbour_pos, grid_size) && continue

                    neighbour_tid = blk_grid[neighbour_pos]
                    neighbour_tid == tid && continue
                    neighbour_workload = length(threads_workload[neighbour_tid])
                    length(thread_workload) ≤ neighbour_workload + 1 && continue  # `+1` after the move

                    # TODO: privilégier `length(thread_workload) > neighbour_workload + 1` avant `length(thread_workload) > neighbour_workload`
                    cost = move_cost(blk_pos, tid, neighbour_tid)
                    candidate = (neighbour_workload, cost, tid, neighbour_tid, blk_pos)
                    insert!(move_candidates, searchsortedfirst(move_candidates, candidate), candidate)
                    break
                end
            end
        end

        while !isempty(move_candidates)
            # Move a block to another neighbouring block with a lower amount of workload                
            _, _, prev_tid, new_tid, blk_pos = popfirst!(move_candidates)
            length(threads_workload[prev_tid]) ≤ ideal_workload && continue
            move_to_group(blk_pos, prev_tid, new_tid)
            performed_changes = true
            break  # TODO: this seems to improve results

            # Update the move cost of the other candidates
            changed_move_candidates = false
            for (i, (neighbour_workload, cost, current_tid, neighbour_tid, pos)) in enumerate(move_candidates)
                if !(neighbour_tid in (prev_tid, new_tid) || sum(abs, Tuple(blk_pos - pos)) == 1)
                    continue  # Only update if same tid or is a neighbour of the block we just moved
                end
                neighbour_workload = length(threads_workload[neighbour_tid])
                cost = move_cost(pos, current_tid, neighbour_tid)
                move_candidates[i] = (neighbour_workload, cost, current_tid, neighbour_tid, pos)
                changed_move_candidates = true
            end
            changed_move_candidates && sort!(move_candidates)
        end
    end

    # TODO: remove?
    if iter == max_iter && max_iter ≥ 20
        @warn "workload diffusion did not converge after $iter iterations"
    else
        @info "workload diffusion convered after $iter iterations"
    end

    return blk_grid, threads_workload
end


function diffuse_workload_matrix!(
    threads_workload::Vector{Vector{CartesianIndex{D}}}, blk_grid::Array{Int, D}, grid_size::Dims{D};
    workload_tolerance=1.01, max_iter=1000
) where {D}
    # All threads should handle a maximum of `min_workload*workload_tolerance` blocks.
    mean_workload = round(prod(grid_size) / length(threads_workload))  # Ideal expected workload
    workload_tolerance = max(1 + 1 / mean_workload, workload_tolerance)  # 5% or 1 whole additional block

    function move_to_group(blk::CartesianIndex, prev_tid, new_tid)
        deleteat!(threads_workload[prev_tid], findfirst(==(blk), threads_workload[prev_tid]))
        push!(threads_workload[new_tid], blk)
        blk_grid[blk] = new_tid
    end

    function move_cost(blk::CartesianIndex, prev_tid, new_tid)
        prev_tid == new_tid && return 0
    
        # The cost is the added (`D-1` dimensional) perimeter length to `new_tid`.
        perimeter_increment = 0
        for side in instances(Side.T)  # TODO: replace by `sides_of(D)`
            neighbour_pos = blk + CartesianIndex(offset_to(side))
            # !in_grid(neighbour_pos, grid_size) && continue
            # perimeter_increment += blk_grid[neighbour_pos] != new_tid ? 1 : -1

            # TODO: this seems to improve results
            if in_grid(neighbour_pos, grid_size)
                perimeter_increment += blk_grid[neighbour_pos] != new_tid ? 1 : -1
            else
                perimeter_increment += 1
            end
        end

        workload_increment = 0  # TODO: get the weight of what we move, e.g. 1 for a static block, something else for an edge block
        return perimeter_increment + workload_increment
    end


    function update_gradient!(workload_gradient)
        min_workload = minimum(length.(threads_workload))
        min_ideal_workload = floor(Int, min_workload * workload_tolerance)
        max_ideal_workload =  ceil(Int, min_workload * workload_tolerance)

        # Initial fill
        for (tid, thread_workload) in enumerate(threads_workload)
            gradient = if min_ideal_workload < length(thread_workload)
                min_ideal_workload - length(thread_workload)  # negative (needs more work)
            elseif length(thread_workload) ≤ max_ideal_workload
                0  # should not change
            else
                max_ideal_workload - length(thread_workload)  # positive (needs less work)
            end

            workload_gradient[tid] = gradient
        end

        # Spread the gradient to neighbours with a basic 
        for (tid, threads_workload) in enumerate(threads_workload)
            workload_gradient[tid] != 0.0 && continue
            neighbours = Int[]
            for blk_pos in threads_workload, side in instances(Side.T)  # TODO: replace by `sides_of(D)`
                neighbour_pos = blk_pos + CartesianIndex(offset_to(side))
                !in_grid(neighbour_pos, grid_size) && continue
                neighbour_tid = blk_grid[neighbour_pos]
                neighbour_tid == tid && continue
                neighbour_tid in neighbours
            end
        end
    end

    workload_gradient = zeros(Float64, length(threads_workload))

    # "Diffuse" the workload by moving blocks from one thread to another, until the condition is met for all threads.
    # Candidates are sorted by ascending target tid workload then ascending move cost.
    move_candidates = Vector{Tuple{Int, Int, Int, Int, CartesianIndex{D}}}()
    performed_changes = true
    iter = 0
    while performed_changes && iter < max_iter
        performed_changes = false
        iter += 1
        min_workload = minimum(length.(threads_workload))
        ideal_workload = ceil(Int, min_workload * workload_tolerance)

        # TODO: debug
        # println("After $(iter-1) iterations (min: $min_workload, max: $(maximum(length.(threads_workload)))):")
        # Main.disp_block_distrib(grid_size, threads_workload)

        # Build a list of possible blocks to move to neighbouring threads
        empty!(move_candidates)
        for (tid, thread_workload) in enumerate(threads_workload)
            length(thread_workload) ≤ ideal_workload && continue

            for blk_pos in thread_workload
                for side in instances(Side.T)  # TODO: replace by `sides_of(D)`
                    neighbour_pos = blk_pos + CartesianIndex(offset_to(side))
                    !in_grid(neighbour_pos, grid_size) && continue

                    neighbour_tid = blk_grid[neighbour_pos]
                    neighbour_tid == tid && continue
                    neighbour_workload = length(threads_workload[neighbour_tid])
                    length(thread_workload) ≤ neighbour_workload && continue

                    cost = move_cost(blk_pos, tid, neighbour_tid)
                    candidate = (neighbour_workload, cost, tid, neighbour_tid, blk_pos)
                    insert!(move_candidates, searchsortedfirst(move_candidates, candidate), candidate)
                    break
                end
            end
        end

        while !isempty(move_candidates)
            # Move a block to another neighbouring block with a lower amount of workload                
            _, _, prev_tid, new_tid, blk_pos = popfirst!(move_candidates)
            length(threads_workload[prev_tid]) ≤ ideal_workload && continue
            move_to_group(blk_pos, prev_tid, new_tid)
            performed_changes = true
            break  # TODO: this seems to improve results, if we keep this we can remove `move_candidates` entierely

            # Update the move cost of the other candidates
            changed_move_candidates = false
            for (i, (neighbour_workload, cost, current_tid, neighbour_tid, pos)) in enumerate(move_candidates)
                if !(neighbour_tid in (prev_tid, new_tid) || sum(abs, Tuple(blk_pos - pos)) == 1)
                    continue  # Only update if same tid or is a neighbour of the block we just moved
                end
                neighbour_workload = length(threads_workload[neighbour_tid])
                cost = move_cost(pos, current_tid, neighbour_tid)
                move_candidates[i] = (neighbour_workload, cost, current_tid, neighbour_tid, pos)
                changed_move_candidates = true
            end
            changed_move_candidates && sort!(move_candidates)
        end
    end

    # TODO: remove?
    if iter == max_iter && max_iter ≥ 20
        @warn "workload diffusion did not converge after $iter iterations"
    end

    return blk_grid, threads_workload
end


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
        strategy=:default, workload_tolerance=0.03, repart=false, retries=0, weighted=false,
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
    strategy=:default, workload_tolerance=0.03, repart=false, retries=0, weighted=false,
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
    strat = Scotch.strat_build(:graph_map; strategy, parts=threads, imbalance_ratio=workload_tolerance)
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
            new_threads_workload = scotch_grid_partition(threads, grid_size; strategy, workload_tolerance, repart)
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
    perimeter_first = params.workload_distribution in (:sorted_square, :sorted_scotch, :weighted_sorted_scotch)
    if params.workload_distribution === :weighted_sorted_scotch
        return thread_workload_distribution(thread_count, grid_size;
            simple, scotch, perimeter_first,
            weighted=true, static_sized_grid, block_size=params.block_size, remainder_block_size, ghosts=params.nghost,
            kwargs...
        )
    else
        return thread_workload_distribution(thread_count, grid_size; simple, scotch, perimeter_first, kwargs...)
    end
end


"""
    thread_workload_distribution(params::ArmonParameters; threads=nothing)
    thread_workload_distribution(
        threads::Int, grid_size::Tuple;
        scotch=true, simple=false, diffuse=true, perimeter_first=false, kwargs...
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
    scotch=true, simple=false, diffuse=true, perimeter_first=false, kwargs...
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
        blk_grid, threads_workload = square_block_distribution(threads, grid_size)
        diffuse && diffuse_workload!(threads_workload, blk_grid, grid_size; kwargs...)
        perimeter_first && sort_blocks_by_perimeter_first!(threads_workload, blk_grid, grid_size)
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
