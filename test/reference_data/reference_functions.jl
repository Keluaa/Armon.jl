
using Printf
using Test
import Armon: @i, @indexing_vars, ArmonData, TestCase, init_test, time_loop 
import Armon: write_data_to_file, read_data_from_file


function get_reference_params(test::Symbol, type::Type; overriden_options...)
    ref_options = Dict(
        :ieee_bits => sizeof(type)*8,
        :test => test, :scheme => :GAD, :projection => :euler_2nd, :riemann_limiter => :minmod,
        :nghost => 5, :nx => 100, :ny => 100,
        :cfl => 0,
        :maxcycle => 1000, :maxtime => 0,  # Always run until reaching the default maximum time of the test
        :silent => 5, :write_output => false, :measure_time => false,
        :use_MPI => false, :async_comms => false
        # TODO: disable SIMD + threading
    )
    merge!(ref_options, overriden_options)
    ArmonParameters(; ref_options...)
end


function run_armon_reference(ref_params::ArmonParameters{T}) where T
    data = ArmonDualData(ref_params)
    init_test(ref_params, data)
    dt, cycles, _ = time_loop(ref_params, data)
    return dt, cycles, data
end


function get_reference_data_file_name(test::TestCase, type::Type)
    test_name = typeof(test).name.name
    return joinpath(@__DIR__, "ref_$(test_name)_$(sizeof(type)*8)bits.csv")
end


function write_reference_data(ref_params::ArmonParameters{T}, ref_file::IO, ref_data::ArmonDualData, 
        dt::T, cycles::Int) where T
    (; nx, ny) = ref_params

    @printf(ref_file, "%#.15g, %d\n", dt, cycles)

    col_range = 1:ny
    row_range = 1:nx
    write_data_to_file(ref_params, ref_data, col_range, row_range, ref_file)
end


function read_reference_data(ref_params::ArmonParameters{T}, ref_file::IO, 
        ref_data::ArmonData{V}) where {T, V <: AbstractArray{T}}
    (; nx, ny) = ref_params
    @indexing_vars(ref_params)

    ref_dt = parse(T, readuntil(ref_file, ','))
    ref_cycles = parse(Int, readuntil(ref_file, '\n'))

    col_range = 1:ny
    row_range = 1:nx
    read_data_from_file(ref_params, ref_data, col_range, row_range, ref_file)

    return ref_dt, ref_cycles
end


# TODO: GPU/CPU still fails, but it is most likely that is is a problem with comparison and tolerance
abs_tol(::Type{Float64}) = 1e-13
abs_tol(::Type{Float32}) = 1e-6
abs_tol(::Flt) where {Flt <: AbstractFloat} = abs_tol(Flt)
rel_tol(::Type{Flt}) where {Flt <: AbstractFloat} = 4*eps(Flt)
rel_tol(::Flt) where {Flt <: AbstractFloat} = rel_tol(Flt)


function count_differences(ref_params::ArmonParameters{T}, 
        data::ArmonData{V}, ref_data::ArmonData{V};
        atol=abs_tol(T), rtol=rel_tol(T)) where {T, V <: AbstractArray{T}}
    (; nx, ny) = ref_params
    @indexing_vars(ref_params)

    differences_count = 0
    for j in 1:ny
        row_range = @i(1,j):@i(nx,j)
        for field in saved_variables()
            ref_row = @view getfield(ref_data, field)[row_range]
            cur_row = @view getfield(data, field)[row_range]
            diff_count = sum(.~ isapprox.(ref_row, cur_row; atol, rtol))
            differences_count += diff_count
            (diff_count > 0) && @debug begin
                max_diff = maximum(abs.((ref_row .- cur_row) .* (.~ isapprox.(ref_row, cur_row; atol, rtol))))
                "Row $j has $diff_count differences in '$field' with the reference. Max diff=$max_diff"
            end
        end
    end

    return differences_count
end


function compare_with_reference_data(ref_params::ArmonParameters{T}, dt::T, cycles::Int, 
        data::ArmonData{V}, ref_data::ArmonData{V}) where {T, V <: AbstractArray{T}}
    ref_file_name = get_reference_data_file_name(ref_params.test, T)

    atol = abs_tol(T)
    rtol = rel_tol(T)

    open(ref_file_name, "r") do ref_file
        ref_dt, ref_cycles = read_reference_data(ref_params, ref_file, ref_data)
        @test ref_dt ??? dt atol=atol rtol=rtol
        @test ref_cycles == cycles
    end

    return count_differences(ref_params, data, ref_data; atol, rtol)
end
