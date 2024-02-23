@testset "Blocking" begin

@testset "BlockGrid" begin
    @testset "$block_size" for block_size in ((32, 32), (16, 48), (110, 110))
        ref_params = get_reference_params(:Sod, Float64; N=(100, 100), block_size)
        grid = Armon.BlockGrid(ref_params)

        @testset "Neighbours" begin
            for blk in Iterators.flatten((grid.device_blocks, grid.device_edge_blocks, grid.remote_blocks))
                if blk isa Armon.RemoteTaskBlock
                    if blk.neighbour isa Armon.RemoteTaskBlock
                        @test blk.neighbour.neighbour == blk
                    else
                        opposite_side = Armon.side_from_offset(Tuple(blk.pos .- blk.neighbour.pos))
                        @test blk.neighbour.neighbours[Int(opposite_side)] == blk
                    end
                else
                    for (side, neighbour) in zip(instances(Armon.Side), blk.neighbours)
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
                blk = Armon.device_block(grid, pos)
                host_blk = Armon.host_block(grid, pos)
                @test blk.pos == pos
                @test blk isa Armon.LocalTaskBlock
                @test blk.mirror == host_blk == blk
            end

            # Remote blocks
            for pos in Iterators.flatten((
                    CartesianIndex(                  0, 1):CartesianIndex(                  0, grid.grid_size[2]),  # left
                    CartesianIndex(grid.grid_size[1]+1, 1):CartesianIndex(grid.grid_size[1]+1, grid.grid_size[2]),  # right
                    CartesianIndex(1,                   0):CartesianIndex(grid.grid_size[1],                   0),  # bottom
                    CartesianIndex(1, grid.grid_size[2]+1):CartesianIndex(grid.grid_size[1], grid.grid_size[2]+1))) # top
                blk = Armon.device_block(grid, pos)
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

        @testset "$side domain" for side in instances(Armon.Side)
            border = Armon.border_domain(bsize, side)
            expected_size = Armon.real_size_along(bsize, Armon.next_axis(Armon.axis_of(side)))
            @test length(border) == expected_size
            if Armon.axis_of(side) == Armon.X_axis
                @test size(border) == (1, expected_size)
            else
                @test size(border) == (expected_size, 1)
            end

            ghost_border = Armon.ghost_domain(bsize, side)
            @test length(ghost_border) == expected_size * Armon.ghosts(bsize)
            if Armon.axis_of(side) == Armon.X_axis
                @test size(ghost_border) == (Armon.ghosts(bsize), expected_size)
            else
                @test size(ghost_border) == (expected_size, Armon.ghosts(bsize))
            end
        end

        @test Armon.stride_along(bsize, Armon.X_axis) == abs(Armon.lin_position(bsize, (1, 1)) - Armon.lin_position(bsize, (2, 1)))
        @test Armon.stride_along(bsize, Armon.Y_axis) == abs(Armon.lin_position(bsize, (1, 1)) - Armon.lin_position(bsize, (1, 2)))
        @test Armon.size_along(bsize, Armon.X_axis) == Armon.size_along(bsize, Armon.Left) == Armon.size_along(bsize, Armon.Right) == Armon.block_size(bsize)[1]
        @test Armon.size_along(bsize, Armon.Y_axis) == Armon.size_along(bsize, Armon.Bottom) == Armon.size_along(bsize, Armon.Top) == Armon.block_size(bsize)[2]
    end
end


@testset "Row iterator" begin
    @testset "$(join(N, '×')) - $(join(B, '×')) - $g" for (g, N, B) in (
        (5, (100, 100), (32, 32)),
        (5, ( 47, 100), (17, 37)),
        (4, ( 96,  96), (32, 32)),  # No edge blocks
        (4, ( 16,  16), (32, 32)),  # Only edge blocks
        (0, (100, 100), (32, 32)),  # 0 ghosts
    )
        params = ArmonParameters(;
            nghost=5, N, block_size=B,
            debug_indexes=true, use_MPI=false, data_type=Float64
        )
        params.nghost = g  # We must do this after the constructor to avoid checks with the schemes
        Armon.update_steps_ranges(params)
        grid = Armon.BlockGrid(params)
        Armon.init_test(params, grid)

        i = 1
        fail_pos = nothing
        for (blk, row_idx, row_range) in Armon.BlockRowIterator(grid)
            expected_length = Armon.real_size_along(blk.size, Armon.X_axis)
            if length(row_range) == expected_length && all((i:i+expected_length-1) .== blk.ρ[row_range])
                i += expected_length
            else
                fail_pos = (blk.pos, row_idx, row_range, blk.ρ[row_range])
                break
            end
        end
        @test fail_pos === nothing
    end
end

end
