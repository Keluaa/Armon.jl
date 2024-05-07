
@testset "logging" begin
    ref_params = get_reference_params(:Sod_circ, Float64; maxcycle=10,
        use_cache_blocking=true, async_cycle=true, log_blocks=true
    )

    @test ref_params.log_blocks
    @test ref_params.estimated_blk_log_size == 10 * 2

    data = BlockGrid(ref_params)
    stats = armon(ref_params)

    @test !isnothing(stats.grid_log)
    grid_stats = Armon.analyse_log_stats(stats.grid_log)

    @test grid_stats.tot_blk == prod(data.grid_size)
    @test grid_stats.inconsistent_threads == 0
    @test length(grid_stats.threads_stats) == grid_stats.active_threads == Threads.nthreads()

    @test occursin(string(grid_stats.stalls_per_thread), sprint(Base.show, MIME"text/plain"(), grid_stats))

    @test_throws Armon.SolverException get_reference_params(:Sod_circ, Float64; maxcycle=10, log_blocks=true)
end
