
using Printf
using Test
import Armon: TestCase


function get_reference_params(test::Symbol, type::Type; overriden_options...)
    ref_options = Dict(
        :data_type => type,
        :test => test, :scheme => :GAD, :projection => :euler_2nd, :riemann_limiter => :minmod,
        :nghost => 4, :N => (100, 100),
        :cfl => 0,
        :maxcycle => 1000, :maxtime => 0,  # Run until reaching the default maximum time of the test by default
        :silent => 5, :write_output => false, :measure_time => false,
        :use_MPI => false
    )
    merge!(ref_options, overriden_options)
    ArmonParameters(; ref_options...)
end


function run_armon_reference(ref_params::ArmonParameters)
    data = BlockGrid(ref_params)
    Armon.init_test(ref_params, data)
    dt, cycles, _, _ = Armon.time_loop(ref_params, data)
    Armon.device_to_host!(data)
    return dt, cycles, data
end


function get_reference_data_file_name(test::TestCase, type::Type)
    test_name = nameof(typeof(test))
    return joinpath(@__DIR__, "ref_$(test_name)_$(sizeof(type)*8)bits.csv")
end


function write_reference_data(
    ref_params::ArmonParameters, ref_file::IO, ref_data::BlockGrid,
    dt, cycles::Int; options...
)
    @printf(ref_file, "%#.15g, %d\n", dt, cycles)
    Armon.write_blocks_to_file(ref_params, ref_data, ref_file; options...)
end


function read_reference_data(ref_params::ArmonParameters, ref_file::IO, ref_data::BlockGrid; options...)
    ref_dt = parse(Armon.data_type(ref_params), readuntil(ref_file, ','))
    ref_cycles = parse(Int, readuntil(ref_file, '\n'))
    Armon.read_data_from_file(ref_params, ref_data, ref_file; options...)
    return ref_dt, ref_cycles
end


# TODO: GPU/CPU still fails, but it is most likely that is is a problem with comparison and tolerance
abs_tol(::Type{Float64}, ::Type{<:TestCase}) = 1e-13
abs_tol(::Type{Float32}, ::Type{<:TestCase}) = 1e-5
rel_tol(::Type{Float64}, ::Type{<:TestCase}) = 4*eps(Float64)
rel_tol(::Type{Float32}, ::Type{<:TestCase}) = 20*eps(Float32)

abs_tol(t::Type, tc::TestCase) = abs_tol(t, typeof(tc))
rel_tol(t::Type, tc::TestCase) = rel_tol(t, typeof(tc))

abs_tol(::Flt) where {Flt <: AbstractFloat} = abs_tol(Flt, TestCase)
rel_tol(::Flt) where {Flt <: AbstractFloat} = rel_tol(Flt, TestCase)

no_zero(x::Flt) where {Flt <: AbstractFloat} = ifelse(iszero(x), nextfloat(zero(Flt)), x)


function count_differences(
    ref_params::ArmonParameters{T}, grid::BlockGrid,
    blk::Armon.LocalTaskBlock{V, Size}, ref_blk::Armon.LocalTaskBlock{V, Size};
    atol=abs_tol(T, ref_params.test), rtol=rel_tol(T, ref_params.test), save_diff=false
) where {T, V <: AbstractArray{T}, Size}
    diff_var = blk.work_1
    save_diff && (diff_var .= 0)

    differences_count = 0
    max_diff = zero(T)
    for (_, row_idx, row_range) in Armon.BlockRowIterator(grid, blk), field in Armon.saved_vars()
        ref_row = @view getfield(ref_blk, field)[row_range]
        cur_row = @view getfield(blk, field)[row_range]

        diff_count = sum((!isapprox).(ref_row, cur_row; atol, rtol))
        differences_count += diff_count
        diff_count == 0 && continue

        if save_diff
            diff_var[row_range] .= (!isapprox).(ref_row, cur_row; atol, rtol)
        end

        row_max_diff = maximum(
            abs.((ref_row .- cur_row) ./ no_zero.(ref_row)) .*
            ((!isapprox).(ref_row, cur_row; atol, rtol))
        )
        max_diff = max(max_diff, row_max_diff)

        @debug begin
            blk_str = join(blk.pos, '×')
            row_str = join(row_idx, '×')
            "In block $blk_str row $row_str has $diff_count differences in '$field' with the reference. Max diff=$row_max_diff"
        end
    end

    return differences_count, max_diff
end


function count_differences(
    ref_params::ArmonParameters{T}, grid::BlockGrid, ref_grid::BlockGrid;
    device_blocks=false, kwargs...
) where {T}
    differences_count = 0
    max_diff = zero(T)
    blk_iter = device_blocks ? Armon.device_blocks : Armon.host_blocks
    for (blk, ref_blk) in zip(blk_iter(grid), blk_iter(ref_grid))
        blk_diff, blk_max_diff = count_differences(ref_params, grid, blk, ref_blk; kwargs...)
        differences_count += blk_diff
        max_diff += blk_max_diff
    end
    return differences_count, max_diff
end


function compare_with_reference_data(
    ref_params::ArmonParameters{T}, dt::T, cycles::Int, 
    grid::BlockGrid, ref_grid::BlockGrid; options...
) where {T}
    ref_file_name = get_reference_data_file_name(ref_params.test, T)

    atol = abs_tol(T, ref_params.test)
    rtol = rel_tol(T, ref_params.test)

    open(ref_file_name, "r") do ref_file
        ref_dt, ref_cycles = read_reference_data(ref_params, ref_file, ref_grid)
        @test ref_dt ≈ dt atol=atol rtol=rtol
        @test ref_cycles == cycles
    end

    return count_differences(ref_params, grid, ref_grid; atol, rtol, options...)
end
