
using Printf
import Armon: Axis, X_axis, Y_axis


function cmp_cpu_with_reference(test::Symbol, type::Type; options...)
    ref_params = get_reference_params(test, type; options...)
    dt, cycles, data = run_armon_reference(ref_params)
    T = data_type(ref_params)
    ref_data = BlockGrid(ref_params)

    differences_count, max_diff = compare_with_reference_data(
        ref_params, dt, cycles,
        data, ref_data,
        save_diff=WRITE_FAILED
    )

    if differences_count > 0 && WRITE_FAILED
        file_name = "test_$(Armon.test_name(ref_params.test))_$(T)"
        open(file_name, "w") do file
            write_reference_data(ref_params, file, data, dt, cycles; more_vars=(:work_1,))
        end
    end

    if !(test in (:Bizarrium, :Sedov))
        @test differences_count == 0
        @test max_diff == 0
    end
end


function axis_invariance(test::Symbol, type::Type, axis::Axis; options...)
    ref_params = get_reference_params(test, type; options...)
    _, _, data = run_armon_reference(ref_params)

    atol = abs_tol(type, ref_params.test)
    rtol = rel_tol(type, ref_params.test)

    vars = setdiff(Armon.main_vars(), (:x, :y))
    vars_errors = Dict(vars .=> 0)

    for blk in Armon.host_blocks(data)
        bsize = blk.size
        blk_size = Armon.block_size(bsize)
        real_blk_size = Armon.real_block_size(bsize)
        real_blk_range = (.+).(Base.oneto.(real_blk_size), Armon.ghosts(bsize))

        r        = ntuple(i -> i == Int(axis) ? real_blk_size[i][1:end-1] : Colon(), ndims(bsize))
        r_offset = ntuple(i -> i == Int(axis) ? real_blk_size[i][2:end]   : Colon(), ndims(bsize))

        for var in vars
            v_data = getfield(blk, var)
            v_data = reshape(v_data, blk_size)  # 1D to n-D array
            v_data = view(v_data, real_blk_range...)  # the real (non-ghost) data
    
            # Substract each axis with its neighbour => if the axis invariance is ok it should be 0
            errors_count = count((!isapprox).(v_data[r...], v_data[r_offset...]; atol, rtol))
            vars_errors[var] += errors_count
        end
    end

    @testset "$var" for var in vars
        @test vars_errors[var] == 0
    end
end


function uninit_vars_propagation(test, type; options...)
    ref_params = get_reference_params(test, type; options...)

    data = BlockGrid(ref_params)
    Armon.init_test(ref_params, data)

    big_val = type == Float32 ? 1e30 : 1e100
    vars = setdiff(Armon.main_vars(), (:x, :y, :mask))
    for blk in Armon.host_blocks(data)
        for i in 1:prod(Armon.block_size(blk.size))
            !Armon.is_ghost(blk.size, i) && continue  # only set `big_val` to ghost cells
            for var in vars
                getfield(blk, var)[i] = big_val
            end
        end
    end

    dt, cycles, _, _ = Armon.time_loop(ref_params, data)

    ref_data = BlockGrid(ref_params)
    differences_count, max_diff = compare_with_reference_data(
        ref_params, dt, cycles, data, ref_data;
        save_diff=WRITE_FAILED
    )

    if differences_count > 0 && WRITE_FAILED
        file_name = "test_$(Armon.test_name(ref_params.test))_$(data_type(ref_params))_uninit_vars"
        open(file_name, "w") do file
            write_reference_data(ref_params, file, data, dt, cycles; more_vars=(:work_1,))
        end
    end

    @test differences_count == 0
    @test max_diff == 0
end


@testset "Convergence" begin
    @testset "Reference" begin
        @testset "$test with $type" for type in (Float32, Float64),
                                        test in (:Sod, :Sod_y, :Sod_circ, :Bizarrium, :Sedov)
            cmp_cpu_with_reference(test, type; use_threading=false, use_simd=false)
        end
    end

    @testset "Axis invariance" begin
        @testset "$test" for (test, axis) in ([:Sod, Y_axis], [:Sod_y, X_axis], [:Bizarrium, Y_axis])
            axis_invariance(test, Float64, axis; use_threading=false, use_simd=false)
        end
    end

    @testset "Uninitialized values propagation" begin
        uninit_vars_propagation(:Sod_circ, Float64; use_threading=false, use_simd=false)
    end
end
