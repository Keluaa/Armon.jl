
import NUMA
using ThreadPinning

@testset "NUMA" begin
    if NUMA.numa_available()
        pinthreads(:cores)
        tid_map = Armon.tid_to_numa_node_map()
    else
        tid_map = Int[]
    end


    function check_for_page_sharing(all_the_pages; disp_first=0)
        # Check if any arrays are sharing a memory page.
        # This would cause issues if e.g. the array are in different blocks, themselves in different
        # NUMA groups: then it would be impossible to efficently place all pages.
        # The check is exhaustive so it runs in 𝓞(N²).
        sharing_count = 0
        for (i, (numa_id, page_range)) in enumerate(all_the_pages),
                (other_numa_id, other_page_range) in all_the_pages[i+1:end]
            if numa_id != other_numa_id && !isdisjoint(page_range, other_page_range)
                sharing_count += 1
                if disp_first > 0
                    disp_first -= 1
                    println("Those two page ranges are both on the NUMA node $numa_id and $other_numa_id: \n\
                              - $page_range\n\
                              - $other_page_range")
                end
            end
        end
        return sharing_count
    end


    function get_all_grid_pages(grid)
        # All `(expected_numa_id, page_range)` pairs for all arrays of every block of the `grid`
        all_pages = Vector{Tuple{Int, StepRange{Ptr{Float64}, Int}}}()
        for (tid, blks_pos) in enumerate(grid.threads_workload), blk_pos in blks_pos
            blk = Armon.block_at(grid, blk_pos)
            numa_id = tid_map[tid]
            for var in Armon.block_vars(blk)
                push!(all_pages, (numa_id, Armon.array_pages(var)))
            end
        end
        return all_pages
    end


    function numa_blk_grid(grid)
        # Returns a NxM grid of which NUMA node is set to the first page of the first variable of all blocks
        numa_blks = Array{Int}(undef, grid.grid_size)
        for blk_pos in eachindex(Base.IndexCartesian(), numa_blks)
            blk = Armon.block_at(grid, blk_pos)
            var = Armon.get_vars(blk, (:x,)) |> only
            idx_mid = (lastindex(var) + firstindex(var)) ÷ 2
            numa_node = NUMA.which_numa_node(var, idx_mid)
            numa_blks[blk_pos] = numa_node
        end
        return numa_blks
    end


    function check_block_numa(blk, test_domain, expected_tid)
        vars = Armon.block_vars()
        blk_range = Armon.block_domain_range(blk.size, test_domain)
        target_numa = tid_map[expected_tid]
        err_count = 0
        first_err_var = nothing
        for (var_name, var) in zip(vars, Armon.get_vars(blk, vars))
            # Only check the first and last index of the array
            first_numa = NUMA.which_numa_node(var, first(blk_range))
            last_numa  = NUMA.which_numa_node(var, last(blk_range))
            if !(first_numa == last_numa == target_numa)
                println("at block $(Tuple(blk.pos)) for var $var_name \
                         expected $target_numa for tid $expected_tid, got:\n\
                          - first_numa = $first_numa (idx=$(first(blk_range)))\n\
                          - last_numa  = $last_numa (idx=$(last(blk_range)))\n")
                err_count == 0 && (first_err_var = var_name)
                err_count += 1
            end
        end
        return err_count, first_err_var
    end


    function check_grid_numa(grid, test_domain)
        err_count = 0
        first_err_blk_pos = nothing
        first_err_var = nothing
        for (tid, blks_pos) in enumerate(grid.threads_workload), blk_pos in blks_pos
            page_err_count, fail_var_name = if Armon.in_grid(blk_pos, grid.static_sized_grid)
                blk = grid.blocks[Armon.block_idx(grid, blk_pos)]
                check_block_numa(blk, test_domain, tid)
            else
                blk = grid.edge_blocks[Armon.edge_block_idx(grid, blk_pos)]
                check_block_numa(blk, test_domain, tid)
            end

            if page_err_count > 0
                if err_count == 0
                    first_err_blk_pos = blk_pos
                    first_err_var = fail_var_name
                end
                err_count += page_err_count
            end
        end
        return err_count, first_err_blk_pos, first_err_var
    end


    if NUMA.numa_available()
        @testset "move_pages" begin
            GC.gc(true)
            # It is a bit difficult to find an array size for which the kernel will always give us
            # pages for which `move_pages` will fail (easiest is EFAULT -> page not present).
            # Those values are quite arbirairy. You may ignore a test failure.
            byte_count = round(Int, min(Sys.maxrss(), Sys.free_memory() * 0.20))
            byte_count = min(byte_count, 50_000_000)
            v = Vector{Int64}(undef, byte_count ÷ sizeof(Int64))
            v_pages = Armon.array_pages(v)
            @test_throws "could not move" Armon.move_pages(v_pages, 1)
            Armon.touch_pages(v_pages)
            @test Armon.move_pages(v_pages, 1) === nothing
        end

        n = round(Int, sqrt(10000 * Threads.nthreads()))  # 100x100 per core
        N = (n, n)
        params = get_reference_params(:Sod_circ, Float64;
            N, block_size=(64, 64), maxcycle=5,
            use_threading=true, use_cache_blocking=true, async_cycle=true,
            workload_distribution=:scotch, numa_aware=true, lock_memory=true
        )

        grid = Armon.BlockGrid(params)
        test_domain = first(params.steps_ranges).full_domain

        # Page positions should be properly enforced after initialization
        Armon.init_test(params, grid)
        @test check_grid_numa(grid, test_domain) == (0, nothing, nothing)
        @test check_for_page_sharing(get_all_grid_pages(grid)) == 0

        # Pages shouldn't move after some solver iterations
        Armon.time_loop(params, grid)
        @test check_grid_numa(grid, test_domain) == (0, nothing, nothing)
    else
        @warn "NUMA is unavailable"
        @test false skip=true
    end
end
