
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


function scotch_graph_part_fixed(graph, parts, initial_partition, args...; repart=false)
    # There is a bug/weird behaviour with `SCOTCH_graphPartFixed` (and others) where fixed vertices
    # are treated as if they didn't exist, yet SCOTCH still checks that their values are between `0`
    # and `parts`. Therefore only the subgraph of non-fixed vertices is partitioned into `parts`.
    # To work around this, we suppose that the parts in fixed vertices are complete, and do the
    # partitioning for `parts - "fixed parts"` instead.
    # Because of the constraint that fixed vertices should be between `0` and `parts`, we fix them
    # all to `0` and correct their value afterward.
    initial_partition_parts = Set{Scotch.SCOTCH_Num}()
    partition::Vector{Scotch.SCOTCH_Num} = map(initial_partition) do part
        part == -1 && return part
        push!(initial_partition_parts, part)
        return 0
    end

    actual_parts = parts - length(initial_partition_parts)
    if actual_parts == 0
        # No partititoning to do: avoid calling SCOTCH since I don't know if it would work
    elseif repart
        partition = Scotch.graph_repart(graph, actual_parts, partition, args...; partition, fixed=true)
    else
        partition = Scotch.graph_part(graph, actual_parts, args...; partition, fixed=true)
    end

    # Map from `sub_graph part index` to `whole graph part index`
    # Since part numbers of the fixed vertices can be any value from `0` to `parts`, we need to
    # transform the newly mapped vertices to a part index which isn't used by the fixed vertices.
    part_mapping = Vector{Scotch.SCOTCH_Num}(undef, actual_parts)
    next_part = 0
    for i in 1:actual_parts
        while next_part in initial_partition_parts
            next_part += 1
        end
        part_mapping[i] = next_part
        next_part += 1
    end

    for (i, (p, init_p)) in enumerate(zip(partition, initial_partition))
        if init_p == -1
            partition[i] = part_mapping[p + 1]
        else
            partition[i] = initial_partition[i]
        end
    end

    return partition
end


"""
    scotch_grid_partition(
        threads, grid_size;
        strategy=:default, workload_tolerance=0, repart=false, retries=10, weighted=false,
        static_sized_grid=nothing, block_size=nothing, remainder_block_size=nothing, ghosts=0,
        initial_partition=nothing
    )

Split `grid_size` to the `threads`.

`strategy` is passed to [`Scotch.strat_flags`](https://keluaa.github.io/Scotch.jl/dev/#Scotch.strat_flags).

`workload_tolerance` is the acceptable workload uneveness among the partitions.
Giving a few blocks of margin (e.g. at least `1/prod(grid_size)`) is preferrable.

`weighted == true` will distribute blocks while taking into account the number of real cells they have.
In this case all parameters of the grid must be present: `static_sized_grid`, `block_size`,
`remainder_block_size` and `ghosts` must be given (obtained with e.g. [`grid_dimensions`](@ref)).

`initial_partition` is the initial partition containing pre-assigned vertices. `-1` is
for movable vertices, while other values are fixed vertices (0-indexed).

The partitioning is random, hence results may vary. To counterbalance this, giving `retries > 0` will
repeat the partitioning `retries` times and keep the best one.
"""
function scotch_grid_partition(
    graph::Scotch.Graph, strat::Scotch.Strat, threads, grid_size;
    repart=false, retries=10, weighted=false, initial_partition=nothing
)
    if isnothing(initial_partition)
        partition = Scotch.graph_part(graph, threads, strat)
    else
        partition = scotch_graph_part_fixed(graph, threads, initial_partition, strat)
    end

    if repart
        cost_factor = 1.0
        costs = partition_cost(threads, grid_size, partition)
        if isnothing(initial_partition)
            partition = Scotch.graph_repart(graph, threads, partition, cost_factor, costs, strat)
        else
            partition = scotch_graph_part_fixed(graph, threads, initial_partition, cost_factor, costs, strat; repart=true)
        end
    end

    threads_workload = map(1:threads) do tid
        indices = findall(==(tid - 1), partition)
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
            new_threads_workload = scotch_grid_partition(graph, strat, threads, grid_size;repart, retries=0, weighted, initial_partition)
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


function side_workload_constraints(partition, side, grid_size::NTuple{N, Int}) where {N}
    # Build the constraints for the sub-domain along `side`

    # Build an iterator for all elements along `side`
    # The iterator should behave the same as the one on the receive side.
    # TODO: better dimension agnosticity?
    axis = Int(axis_of(side))
    ax_pos = side in first_sides() ? 1 : grid_size[axis]
    side_iter = CartesianIndices(Tuple(
        ifelse.(1:N .== axis, Ref(ax_pos:ax_pos), Base.OneTo.(grid_size))
    ))
    side_root = first(side_iter)  # transferred block positions are relative to this corner

    side_constraints = Vector{Tuple{CartesianIndex{N}, Scotch.SCOTCH_Num}}()
    side_size = prod(grid_size) ÷ grid_size[axis]
    sizehint!(side_constraints, side_size)

    for (tid, workload) in enumerate(partition)
        isdisjoint(workload, side_iter) && continue
        # This thread has blocks on the side of the grid: add all blocks to the constraints.
        # While we could only transfer those on the edge, passing whole groups gives make things
        # easier for SCOTCH, which will not have to handle incomplete groups.
        # From what I could see, having incomplete groups in constraints gives low quality output,
        # with some disjoint groups and uneven workloads, no matter the amount of retries.
        for blk_pos in workload
            rel_blk_pos = blk_pos - side_root
            if rel_blk_pos[axis] != 0
                # Since we are sending relative block positions to the other side, we need to mirror
                # the axis of the side (X for left and right, etc...): this way the other side can
                # use all positions directly.
                rel_blk_pos = CartesianIndex(Tuple(ifelse.(1:N .== axis, Tuple(-rel_blk_pos), Tuple(rel_blk_pos))))
            end
            push!(side_constraints, (rel_blk_pos, tid - 1))  # tid-1 as SCOTCH uses 0-indices
        end
    end

    return side_constraints
end


function side_workload_constraints!(
    initial_partition, side, grid_size::NTuple{N, Int},
    side_constraints::Vector{Tuple{CartesianIndex{N}, Scotch.SCOTCH_Num}}
) where {N}
    # Update the `initial_partition` with the constraints from the sub-domain along `side`.
    # Return the new partition and the fixed groups added to it.

    # Build an iterator for all elements along `side`
    # TODO: better dimension agnosticity?
    axis = Int(axis_of(side))
    ax_pos = side in first_sides() ? 1 : grid_size[axis]
    side_iter = CartesianIndices(Tuple(
        ifelse.(1:length(grid_size) .== axis, Ref(ax_pos:ax_pos), Base.OneTo.(grid_size))
    ))
    side_root = first(side_iter)  # transferred block positions are relative to this corner

    fixed_groups = Set{Scotch.SCOTCH_Num}()
    for (rel_blk_pos, tid) in side_constraints
        blk_pos = rel_blk_pos + side_root

        # The neighbouring sub-domain has no information about the dimensions of our block grid.
        # This means it can send block positions which aren't part of our grid (this shouldn't
        # concern blocks along the exact side).
        !(blk_pos in CartesianIndices(grid_size)) && continue

        initial_partition[LinearIndices(grid_size)[blk_pos]] = tid
        push!(fixed_groups, tid)
    end

    return initial_partition, fixed_groups
end


function scotch_grid_matched_partition(grid_graph, strat, threads, grid_size, comm::MPI.Comm; kwargs...)
    # The partitioning is global among the processes of the communicator.
    # Neighbour processes will need the same threads associated with the blocks along sides they share.

    # The sub-domain at the center of the grid is the one which initiates the partitioning.
    # We could also start from anywhere else, but the center is the most optimal one, as it is
    # the barycenter of the grid: it minimizes the number of waits required to complete the
    # whole partitioning.
    cart_size, _, coords = MPI.Cart_get(comm)
    center_pos = fld.(cart_size .- 1, 2)
    is_center = all(center_pos .== coords)
    dim = length(grid_size)

    # Get the neighbours which we will communicate with
    # TODO: dimension agnostic
    neighbours = [
        Side.Left   => MPI.Cart_shift(comm, 0, -1)[2],
        Side.Right  => MPI.Cart_shift(comm, 0,  1)[2],
        Side.Bottom => MPI.Cart_shift(comm, 1, -1)[2],
        Side.Top    => MPI.Cart_shift(comm, 1,  1)[2]
    ]

    center_offset = coords .- center_pos
    center_dist = sum(abs.(center_offset))
    recv_neighbours = filter(neighbours) do (side, rank)
        rank_dist = sum(abs.(center_offset .+ offset_to(side)))
        return rank != MPI.PROC_NULL && rank_dist ≤ center_dist
    end
    send_neighbours = filter(neighbours) do (side, rank)
        rank_dist = sum(abs.(center_offset .+ offset_to(side)))
        return rank != MPI.PROC_NULL && rank_dist > center_dist
    end

    if is_center
        # The center domain has no constraints
        partition = scotch_grid_partition(grid_graph, strat, threads, grid_size; kwargs...)
    else
        # Wait for the neighbours closer (Manhattan distance) to the center to communicate their
        # edge partitions

        # Get the number of constrained blocks each neighbour wants to send
        reqs = MPI.MultiRequest(length(recv_neighbours))
        side_constraints_sizes = Vector{Ref{Int}}(undef, 2*dim)
        for (req, (side, rank)) in zip(reqs, recv_neighbours)
            side_constraints_sizes[Int(side)] = Ref(-1)
            MPI.Irecv!(side_constraints_sizes[Int(side)], comm, req; source=rank, tag=0)
        end
        MPI.Waitall(reqs)

        # Receive all constrained blocks
        constraints = Dict{Side.T, Vector{Tuple{CartesianIndex{dim}, Scotch.SCOTCH_Num}}}()
        for (req, (side, rank)) in zip(reqs, recv_neighbours)
            side_size = side_constraints_sizes[Int(side)][]
            constraints[side] = valtype(constraints)(undef, side_size)
            MPI.Irecv!(constraints[side], comm, req; source=rank, tag=1)
        end
        MPI.Waitall(reqs)

        # Build the initial partition with the constraints of each neighbour
        initial_partition = Vector{Scotch.SCOTCH_Num}(undef, prod(grid_size))
        initial_partition .= -1  # all vertices have no constraints by default
        fixed_groups = Set{Scotch.SCOTCH_Num}()
        for (side, side_cstr) in constraints
            _, side_groups = side_workload_constraints!(initial_partition, side, grid_size, side_cstr)
            union!(fixed_groups, side_groups)
        end

        partition = scotch_grid_partition(grid_graph, strat, threads, grid_size; initial_partition, kwargs...)
    end

    # Send our edge partitions to the neighbours farther from the center
    reqs = MPI.MultiRequest(2*length(send_neighbours))
    for (i, (side, rank)) in enumerate(send_neighbours)
        constraint = side_workload_constraints(partition, side, grid_size)

        # Send the length and the data at once, but with a different tag to not anger the MPI gods
        MPI.Isend(length(constraint), comm, reqs[2*i-1]; dest=rank, tag=0)
        MPI.Isend(constraint, comm, reqs[2*i]; dest=rank, tag=1)
    end

    MPI.Waitall(reqs)  # Don't leave any requests incomplete
    MPI.Barrier(comm)  # Just to be safe, maybe unneeded

    return partition
end


function scotch_grid_partition(
    threads, grid_size;
    strategy=:default, workload_tolerance=0, weighted=false,
    static_sized_grid=nothing, block_size=nothing, remainder_block_size=nothing, ghosts=0,
    match_neighbour_domains::Union{Bool, MPI.Comm}=false,
    kwargs...
)
    # TODO: for larger grids, using graph coarsening might be necessary (+ it may help the solver to reach better solutions)
    graph = grid_to_scotch_graph(grid_size;
        weighted, static_sized_grid, block_size, remainder_block_size, ghosts
    )
    # TODO: with `Scotch.graph_map` it is possible to partition the grid according to a topology of threads, is it better?
    # -> by adding weights to the topology it might be possible to reproduce the processor topology to account for e.g. NUMA or caches
    # -> cores sharing the same L3      => dense graph + weight 1
    # -> cores on the same NUMA         => dense graph + weight 2?
    # -> cores on the same socket       => dense graph + weight 4?
    # -> cores on the different sockets => dense graph + weight 8?

    strat = Scotch.strat_build(:graph_map; strategy, parts=threads, imbalance_ratio=Float64(workload_tolerance))

    if match_neighbour_domains == false
        return scotch_grid_partition(graph, strat, threads, grid_size; weighted, kwargs...)
    else
        return scotch_grid_matched_partition(graph, strat, threads, grid_size, match_neighbour_domains; weighted, kwargs...)
    end
end


function thread_workload_distribution(params::ArmonParameters; threads=nothing, kwargs...)
    thread_count = @something threads params.nthreads
    grid_size, static_sized_grid, remainder_block_size = grid_dimensions(params)
    simple = params.workload_distribution === :simple
    scotch = params.workload_distribution in (:scotch, :sorted_scotch, :weighted_sorted_scotch)
    perimeter_first = params.workload_distribution in (:sorted_scotch, :weighted_sorted_scotch)
    match_neighbour_domains = params.thread_split_comm ? params.cart_comm : false
    merged_kw = merge(params.distrib_params, Dict(kwargs...))
    if params.workload_distribution === :weighted_sorted_scotch
        return thread_workload_distribution(thread_count, grid_size;
            simple, scotch, perimeter_first,
            weighted=true, static_sized_grid, block_size=params.block_size, remainder_block_size, ghosts=params.nghost,
            match_neighbour_domains, merged_kw...
        )
    else
        return thread_workload_distribution(thread_count, grid_size;
            simple, scotch, perimeter_first, match_neighbour_domains, merged_kw...
        )
    end
end


"""
    thread_workload_distribution(params::ArmonParameters; threads=nothing)
    thread_workload_distribution(
        threads::Int, grid_size::Tuple;
        scotch=true, simple=false, perimeter_first=false, check=true, match_neighbour_domains=false,
        kwargs...
    )

Distribute each block in `grid_size` among the `threads`, as evenly as possible.

With `simple == true`, blocks are distributed with [`simple_workload_distribution`](@ref).

With `scotch == true`, the [`Scotch`](https://gitlab.inria.fr/scotch/scotch) graph partitioning solver
is used to better split the grid.
`kwargs` are passed to [`scotch_grid_partition`](@ref).

If `perimeter_first == true`, the resulting distribution will have blocks sorted in way that will
place neighbours of other threads' blocks first in the list.
By doing so, communications between threads may be overlapped more frequently.

If `check == true`, the resulting distribution is checked to ensure that all blocks are assigned to
a thread exactly once.

If `match_neighbour_domains == comm`, then we have an additional constraint: two neighbouring
sub-domains must share the same distribution of blocks and threads on the side they share.
The partitioning problem becomes global, and is resolved by starting from the center domain (without
any constraints), which communicates the distribution on each side to its neighbours, becoming a
partitioning constraint. They each solve their distribution, and the process continues until it
reaches all sub-domains.
`comm` is a cartesian `MPI.Comm`unicator, restraining the algorithm to ranks part of `comm` only.
It is incompatible with the `simple == true`.
"""
function thread_workload_distribution(
    threads::Int, grid_size::Tuple;
    scotch=true, simple=false, perimeter_first=false, check=true, match_neighbour_domains=false,
    kwargs...
)
    if simple
        if match_neighbour_domains != false
            error("simple workload distribution is incompatible with `match_neighbour_domains`")
        end
        threads_workload = simple_workload_distribution(threads, grid_size)
    elseif scotch
        threads_workload = scotch_grid_partition(threads, grid_size; match_neighbour_domains, kwargs...)
    else
        error("unknown workload distribution, expected `simple == true` or `scotch == true`")
    end

    if check || perimeter_first
        blk_grid = block_grid_from_workload(grid_size, threads_workload)
    end
    perimeter_first && sort_blocks_by_perimeter_first!(threads_workload, blk_grid, grid_size)
    check && check_workload(threads_workload, blk_grid)
    return threads_workload
end


function check_workload(threads_workload, blk_grid)
    unassigned_blocks = count(==(0), blk_grid)
    if unassigned_blocks > 0
        plurial = unassigned_blocks > 0 ? "blocks are" : "block is"
        error("invalid block distribution: $unassigned_blocks $plurial not assigned to a thread")
    end

    for (tid, workload) in enumerate(threads_workload), blk_pos in workload
        blk_grid[blk_pos] == tid && continue
        error("invalid block distribution: block $(Tuple(blk_pos)) is assigned to threads $(blk_grid[blk_pos]) and $tid")
    end

    return true
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
