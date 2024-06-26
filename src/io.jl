
# TODO: use hdf5 for much more efficient read/write

function write_blocks_to_file(
    params::ArmonParameters, grid::BlockGrid, file::IO, row_iter_params...;
    global_ghosts=false, all_ghosts=false, for_3D=true, vars=saved_vars()
)
    p = params.output_precision
    format = Printf.Format(join(repeat(["%#$(p+7).$(p)e"], length(vars)), ", ") * "\n")

    # Write cells in the correct ascending (X, Y, Z) order, combining the cells of all blocks
    prev_row_idx = nothing
    for (blk, row_idx, row_range) in BlockRowIterator(grid, row_iter_params...; global_ghosts, all_ghosts)
        row_idx = row_idx[2:end]
        if prev_row_idx != row_idx && !isnothing(prev_row_idx)
            for_3D && println(file)  # Separate rows to use pm3d plotting with gnuplot
        end

        blk_vars = get_vars(blk, vars; on_device=false)
        # TODO: center the positions of the cells
        for idx in row_range
            Printf.format(file, format, getindex.(blk_vars, idx)...)
        end

        prev_row_idx = row_idx
    end
end


function read_data_from_file(
    ::ArmonParameters{T}, grid::BlockGrid, file::IO, row_iter_params...;
    global_ghosts=false, all_ghosts=false, vars=saved_vars()
) where {T}
    for (blk, _, row_range) in BlockRowIterator(grid, row_iter_params...; global_ghosts, all_ghosts)
        blk_vars = get_vars(blk, vars; on_device=false)
        for idx in row_range
            for var in blk_vars[1:end-1]
                var[idx] = parse(T, readuntil(file, ','))
            end
            blk_vars[end][idx] = parse(T, readuntil(file, '\n'))
        end
    end
end


function build_file_path(params::ArmonParameters, file_name::String)
    file_path = joinpath(params.output_dir, file_name)

    if params.is_root && !isdir(params.output_dir)
        mkpath(params.output_dir)
    end

    if params.use_MPI
        coords_str = join(params.cart_coords, '×')
        params.use_MPI && (file_path *= "_$coords_str")
    end

    return file_path
end


function write_sub_domain_file(
    params::ArmonParameters, data::BlockGrid, file_name::String;
    no_msg=false, options...
)
    output_file_path = build_file_path(params, file_name)
    open(output_file_path, "w") do file
        write_blocks_to_file(params, data, file; global_ghosts=params.write_ghosts, options...)
    end

    if !no_msg && params.is_root && params.silent < 2
        println("\nWrote to files $(output_file_path)_*x*")
    end
end


function read_sub_domain_file!(
    params::ArmonParameters, data::BlockGrid, file_name::String; options...
)
    output_file_path = build_file_path(params, file_name)
    open(output_file_path, "r") do file
        read_data_from_file(params, data, file; global_ghosts=params.write_ghosts, options...)
    end
end


function write_time_step_file(params::ArmonParameters, state::SolverState, file_name::String)
    file_path = build_file_path(params, file_name)

    p = params.output_precision
    format = Printf.Format("%#$(p+7).$(p)e\n")

    open(file_path, "w") do file
        Printf.format(file, format, state.global_dt.current_dt)
    end
end


function read_time_step_file(params::ArmonParameters{T}, file_name::String) where {T}
    file_path = build_file_path(params, file_name)

    open(file_path, "r") do file
        return parse(T, readchomp(file))
    end
end

#
# Comparison functions
#

function compare_block(
    params::ArmonParameters, ref_blk::LocalTaskBlock, our_blk::LocalTaskBlock, label::String;
    vars=saved_vars()
)
    different = false

    real_static_bsize = params.block_size .- 2*params.nghost
    blk_global_pos = params.N_origin .- 1 .+ (Tuple(our_blk.pos) .- 1) .* real_static_bsize

    ref_vars = get_vars(ref_blk, vars; on_device=false)
    our_vars = get_vars(our_blk, vars; on_device=false)
    for (var, ref_var, our_var) in zip(vars, ref_vars, our_vars)
        diff_mask = (!isapprox).(ref_var, our_var; rtol=params.comparison_tolerance)
        !params.write_ghosts && (diff_mask .*= (!is_ghost).(Ref(our_blk.size), 1:prod(block_size(our_blk))))

        diff_count = sum(diff_mask)
        diff_count == 0 && continue

        !different && println("At $label, in block $(our_blk.pos):")
        different = true
        print("  $diff_count differences found in $var")

        if diff_count ≤ 200
            println(" (ref ≢ current)")
            for (idx, mask) in enumerate(diff_mask)
                !mask && continue
                I = position(our_blk.size, idx)
                gI = I .+ blk_global_pos .- 1

                val_diff = ref_var[idx] - our_var[idx]
                diff_ulp = val_diff / eps(ref_var[idx])
                abs(diff_ulp) > 1e10 && (diff_ulp = Inf)

                pos_str  = join((@sprintf("%3d", i) for i in I ), ',')
                gpos_str = join((@sprintf("%3d", i) for i in gI), ',')
                @printf("   - %5d (%s | %s): %12.5g ≢ %12.5g (%12.5g, ulp: %8g)\n",
                    idx, pos_str, gpos_str, ref_var[idx], our_var[idx], val_diff, diff_ulp)
            end
        else
            println()
        end
    end

    return different
end


function compare_data(
    params::ArmonParameters, ref_data::BlockGrid, our_data::BlockGrid, label::String;
    vars=saved_vars()
)
    different = false
    for (ref_blk, our_blk) in zip(all_blocks(ref_data), all_blocks(our_data))
        different |= compare_block(params, ref_blk, our_blk, label; vars)
    end
    return different
end


function compare_with_file(
    params::ArmonParameters, grid::BlockGrid, file_name::String, label::String
)
    ref_data = BlockGrid(params)
    read_sub_domain_file!(params, ref_data, file_name)
    different = compare_data(params, ref_data, grid, label)

    if params.use_MPI
        different = MPI.Allreduce(different, |, params.cart_comm)
    end

    return different
end


function step_checkpoint(params::ArmonParameters, state::SolverState, grid::BlockGrid, step_label::String)
    !params.compare && return false

    wait(params)
    device_to_host!(grid)
    wait(params)

    step_file_name = params.output_file * @sprintf("_%03d_%s", state.global_dt.cycle, step_label)
    if state.global_dt.cycle == 0 && step_label == "time_step"
        axis = Axis.X
    else
        axis = state.axis
    end
    step_file_name *= "_" * string(axis)[1]

    if params.is_ref
        if step_label == "time_step"
            write_time_step_file(params, state, step_file_name)
        else
            write_sub_domain_file(params, grid, step_file_name; no_msg=true)
        end

        return false
    else
        if step_label == "time_step"
            ref_dt = read_time_step_file(params, step_file_name)
            different = !isapprox(ref_dt, state.dt; rtol=params.comparison_tolerance)
            if different
                @printf("Time step difference: ref Δt = %.18f, Δt = %.18f, diff = %.18f\n",
                        ref_dt, state.dt, ref_dt - state.dt)
            end
        else
            different = compare_with_file(params, grid, step_file_name, step_label)
        end

        if different
            write_sub_domain_file(params, grid, step_file_name * "_diff"; no_msg=true)
            println("Difference file written to $(step_file_name)_diff")
        end

        return different
    end
end
