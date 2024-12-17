
using Kokkos
if !Kokkos.is_initialized()
    Kokkos.set_omp_vars()
    Kokkos.set_backends(split(TEST_KOKKOS_BACKEND, ',') .|> strip)
    Kokkos.initialize()
end


function run_armon_cpp_reference(ref_params::ArmonParameters)
    data = BlockGrid(ref_params)
    Armon.init_test(ref_params, data)
    _, dt, cycles, _, _ = Armon.time_loop(ref_params, data)
    Armon.device_to_host!(data)
    return dt, cycles, data
end


function cmp_cpp_with_reference_for(type, test; kwargs...)
    debug_kernels = true
    armon_cpp_lib_src = TEST_KOKKOS_PATH
    use_md_iter = true
    use_simd = !use_md_iter
    ref_params = get_reference_params(test, type; use_kokkos=true, debug_kernels, armon_cpp_lib_src, use_md_iter, use_simd, kwargs...)
    dt, cycles, data = run_armon_cpp_reference(ref_params)
    ref_data = BlockGrid(ref_params)

    differences_count, max_diff = compare_with_reference_data(
        ref_params, dt, cycles,
        data, ref_data;
        save_diff=WRITE_FAILED
    )

    if differences_count > 0 && WRITE_FAILED
        file_name = "test_kokkos_$(Armon.test_name(ref_params.test))_$(Armon.data_type(ref_params))"
        open(file_name, "w") do file
            write_reference_data(ref_params, file, data, dt, cycles; more_vars=(:work_1,))
        end
    end

    @test differences_count == 0
    @test max_diff == 0
end


@testset "Kokkos" begin
    @testset "Reference" begin
        @testset "$test with $type" for type in (Float64,),
                                        test in (:Sod, :Sod_y, :Sod_circ, :Bizarrium, :Sedov)
            cmp_cpp_with_reference_for(type, test)
        end
    end
end
