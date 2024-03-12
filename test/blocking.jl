@testset "Blocking" begin

@testset "BlockGrid" begin
    @testset "$block_size" for block_size in (
            ( 32,  32),  # Normal mix of static and edge blocks
            ( 16,  48),  # Uneven mix of static and edge blocks
            ( 57,  57),  # Edge blocks too small, merged with closest static blocks
            ( 64,  57),  # Same but uneven
            (106, 106),  # Edge blocks too small, merged with closest static blocks (only edge blocks left)
            (108, 108),  # Perfect match: only a single static block
            (  0,   0),  # No blocking, only edge blocks
        )
        nghost = 4
        ref_params = get_reference_params(:Sod, Float64; N=(100, 100), block_size, nghost)
        grid = Armon.BlockGrid(ref_params)

        @testset "Grid" begin
            if prod(block_size) > 0
                @test prod(grid.grid_size) == length(grid.blocks) + length(grid.edge_blocks)
                grid_size, static_grid, remainder_size = Armon.grid_dimensions(ref_params)
                @test all(remainder_size .> nghost .|| remainder_size .== 0)
            else
                @test isempty(grid.blocks)
                @test !isempty(grid.edge_blocks)
                @test prod(grid.grid_size) == length(grid.edge_blocks)
            end

            @test Armon.device_array_type(grid) == Armon.host_array_type(grid) == Armon.buffer_array_type(grid) == Vector{Float64}
            @test Armon.device_is_host(grid)
            @test Armon.buffers_on_device(grid)
            @test Armon.ghosts(grid) == nghost
            @test Armon.block_size(Armon.static_block_size(grid)) == block_size
        end

        @testset "Domain" begin
            real_cell_count = 0
            all_cell_count = 0
            blk_bytes = 0
            for blk in Armon.all_blocks(grid)
                real_cell_count += prod(Armon.real_block_size(blk.size))
                all_cell_count  += prod(Armon.block_size(blk.size))
                blk_bytes += sizeof(blk)
            end

            @test real_cell_count == prod(ref_params.N)
            if prod(grid.grid_size) == 1
                @test all_cell_count  == prod(ref_params.N .+ 2*ref_params.nghost)
            end

            total_mem = all_cell_count * sizeof(Float64) * length(Armon.block_vars())
            remote_blocks_overhead = sum(sizeof.(grid.remote_blocks))
            device_memory, host_memory = Armon.memory_required(ref_params)
            @test total_mem == device_memory
            @test host_memory == total_mem + blk_bytes + remote_blocks_overhead
        end

        @testset "Neighbours" begin
            for blk in Iterators.flatten((grid.blocks, grid.edge_blocks, grid.remote_blocks))
                if blk isa Armon.RemoteTaskBlock
                    if blk.neighbour isa Armon.RemoteTaskBlock
                        @test blk.neighbour.neighbour == blk
                    else
                        opposite_side = Armon.side_from_offset(Tuple(blk.pos .- blk.neighbour.pos))
                        @test blk.neighbour.neighbours[Int(opposite_side)] == blk
                    end
                else
                    for (side, neighbour) in zip(instances(Armon.Side.T), blk.neighbours)
                        if neighbour isa Armon.RemoteTaskBlock
                            @test neighbour.neighbour == blk
                        else 
                            opposite_side = Armon.opposite_of(side)
                            @test neighbour.neighbours[Int(opposite_side)] == blk
                        end
                    end
                end
            end
        end

        @testset "Block position" begin
            for pos in CartesianIndex(1, 1):CartesianIndex(grid.grid_size)
                blk = Armon.block_at(grid, pos)
                @test blk.pos == pos
                @test blk isa Armon.LocalTaskBlock
            end

            # Remote blocks
            for pos in Iterators.flatten((
                    CartesianIndex(                  0, 1):CartesianIndex(                  0, grid.grid_size[2]),  # left
                    CartesianIndex(grid.grid_size[1]+1, 1):CartesianIndex(grid.grid_size[1]+1, grid.grid_size[2]),  # right
                    CartesianIndex(1,                   0):CartesianIndex(grid.grid_size[1],                   0),  # bottom
                    CartesianIndex(1, grid.grid_size[2]+1):CartesianIndex(grid.grid_size[1], grid.grid_size[2]+1))) # top
                blk = Armon.block_at(grid, pos)
                @test blk.pos == pos
                @test blk isa Armon.RemoteTaskBlock
            end
        end
    end
end


@testset "Block size" begin
    @testset "$bsize" for bsize in (
            Armon.StaticBSize((64, 64), 5),
            Armon.StaticBSize((37, 39), 5),
            Armon.DynamicBSize((64, 37), 5),
            Armon.DynamicBSize((37, 39), 1),
            Armon.DynamicBSize((37, 39), 0))

        @test Armon.ghosts(bsize) in (5, 1, 0)
        @test Armon.real_block_size(bsize) == Armon.block_size(bsize) .- 2*Armon.ghosts(bsize)
        @test ndims(bsize) == 2

        @testset "Full domain" begin
            g = Armon.ghosts(bsize)
            all_cells = Armon.block_domain_range(bsize, (-g, -g), (g, g))
            @test size(all_cells) == Armon.block_size(bsize)

            pos_ok = 0
            lin_pos_ok = 0
            ghost_ok = 0
            for (ij, j) in enumerate(all_cells.col), (ii, i) in enumerate(all_cells.row .+ (j - 1))
                I = (ii - g, ij - g)
                pos_ok     += Armon.position(bsize, i) == I
                lin_pos_ok += Armon.lin_position(bsize, I) == i
                ghost_ok   += Armon.is_ghost(bsize, i) == (any(I .≤ 0) || any(I .> Armon.real_block_size(bsize)))
            end
            @test pos_ok     == length(all_cells)
            @test lin_pos_ok == length(all_cells)
            @test ghost_ok   == length(all_cells)
        end

        @testset "Real domain" begin
            real_cells = Armon.block_domain_range(bsize, (0, 0), (0, 0))            
            @test size(real_cells) == Armon.real_block_size(bsize)

            pos_ok = 0
            lin_pos_ok = 0
            ghost_ok = 0
            for (ij, j) in enumerate(real_cells.col), (ii, i) in enumerate(real_cells.row .+ (j - 1))
                I = (ii, ij)
                pos_ok     += Armon.position(bsize, i) == I
                lin_pos_ok += Armon.lin_position(bsize, I) == i
                ghost_ok   += !Armon.is_ghost(bsize, i)
            end
            @test pos_ok     == length(real_cells)
            @test lin_pos_ok == length(real_cells)
            @test ghost_ok   == length(real_cells)
        end

        @testset "$side domain" for side in instances(Armon.Side.T)
            border = Armon.border_domain(bsize, side)
            expected_size = Armon.real_size_along(bsize, Armon.next_axis(Armon.axis_of(side)))
            @test length(border) == expected_size
            if Armon.axis_of(side) == Armon.Axis.X
                @test size(border) == (1, expected_size)
            else
                @test size(border) == (expected_size, 1)
            end

            ghost_border = Armon.ghost_domain(bsize, side)
            @test length(ghost_border) == expected_size * Armon.ghosts(bsize)
            if Armon.axis_of(side) == Armon.Axis.X
                @test size(ghost_border) == (Armon.ghosts(bsize), expected_size)
            else
                @test size(ghost_border) == (expected_size, Armon.ghosts(bsize))
            end
        end

        @test Armon.stride_along(bsize, Armon.Axis.X) == abs(Armon.lin_position(bsize, (1, 1)) - Armon.lin_position(bsize, (2, 1)))
        @test Armon.stride_along(bsize, Armon.Axis.Y) == abs(Armon.lin_position(bsize, (1, 1)) - Armon.lin_position(bsize, (1, 2)))
        @test Armon.size_along(bsize, Armon.Axis.X) == Armon.size_along(bsize, Armon.Side.Left) == Armon.size_along(bsize, Armon.Side.Right) == Armon.block_size(bsize)[1]
        @test Armon.size_along(bsize, Armon.Axis.Y) == Armon.size_along(bsize, Armon.Side.Bottom) == Armon.size_along(bsize, Armon.Side.Top) == Armon.block_size(bsize)[2]
    end
end


@testset "Row iterator" begin
    @testset "$(join(N, '×')) - $(join(B, '×')) - $g" for (g, N, B) in (
        (5, (100, 100), (32, 32)),
        (5, ( 47, 100), (17, 37)),
        (4, ( 96,  96), (32, 32)),  # No edge blocks
        (4, ( 16,  16), (32, 32)),  # Only edge blocks
        (3, ( 16,  33), ( 0,  0)),  # No blocking, only edge blocks
        (4, (100, 100), (57, 57)),  # Edge blocks bigger than static blocks
        (0, (100, 100), (32, 32)),  # 0 ghosts
    )
        params = ArmonParameters(;
            test=:DebugIndexes, nghost=5, N, block_size=B,
            use_MPI=false, data_type=Float64
        )
        params.nghost = g  # We must do this after the constructor to avoid checks with the schemes
        Armon.compute_steps_ranges(params)
        grid = Armon.BlockGrid(params)
        Armon.init_test(params, grid)

        i = 1
        fail_pos = nothing
        for (blk, row_idx, row_range) in Armon.BlockRowIterator(grid)
            expected_length = Armon.real_size_along(blk.size, Armon.Axis.X)
            blk_data = Armon.block_host_data(blk)
            if length(row_range) == expected_length && all((i:i+expected_length-1) .== blk_data.ρ[row_range])
                i += expected_length
            else
                fail_pos = (blk.pos, row_idx, row_range, blk_data.ρ[row_range])
                break
            end
        end
        @test fail_pos === nothing
        @test i - 1 == prod(N)
    end
end

end
