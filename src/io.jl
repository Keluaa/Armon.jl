
# TODO: use hdf5 for much more efficient read/write

function write_blocks_to_file(
    params::ArmonParameters, grid::BlockGrid, file::IO;
    global_ghosts=false, all_ghosts=false, for_3D=true, more_vars=()
)
    vars_to_write = tuple(saved_vars()..., more_vars...)

    p = params.output_precision
    format = Printf.Format(join(repeat(["%#$(p+7).$(p)e"], length(vars_to_write)), ", ") * "\n")

    # Write cells in the correct ascending (X, Y, Z) order, combining the cells of all blocks
    prev_row_idx = nothing
    for (blk, row_idx, row_range) in BlockRowIterator(grid; global_ghosts, all_ghosts, device_blocks=false)
        if prev_row_idx != row_idx && !isnothing(prev_row_idx)
            for_3D && println(file)  # Separate rows to use pm3d plotting with gnuplot
        end

        vars = getfield.(Ref(blk), vars_to_write)
        # TODO: center the positions of the cells
        for idx in row_range
            Printf.format(file, format, getindex.(vars, idx)...)
        end

        prev_row_idx = row_idx
    end
end


function read_data_from_file(
    ::ArmonParameters{T}, grid::BlockGrid, file::IO;
    global_ghosts=false, all_ghosts=false, more_vars=()
) where {T}
    vars_to_read = tuple(saved_vars()..., more_vars...)

    for (blk, row_range) in BlockRowIterator(grid; global_ghosts, all_ghosts, device_blocks=false)
        vars = getfield.(Ref(blk), vars_to_read)
        for idx in row_range
            for var in vars[1:end-1]
                var[idx] = parse(T, readuntil(file, ','))
            end
            vars[end][idx] = parse(T, readuntil(file, '\n'))
        end
    end
end


function build_file_path(params::ArmonParameters, file_name::String)
    (; output_dir, use_MPI, is_root, cart_coords) = params

    file_path = joinpath(output_dir, file_name)

    if is_root && !isdir(output_dir)
        mkpath(output_dir)
    end

    if use_MPI
        (cx, cy) = cart_coords
        params.use_MPI && (file_path *= "_$(cx)x$(cy)")
    end

    return file_path
end


function write_sub_domain_file(
    params::ArmonParameters, data::BlockGrid, file_name::String;
    no_msg=false, options...
)
    output_file_path = build_file_path(params, file_name)
    open(output_file_path, "w") do file
        global_ghosts = params.write_ghosts
        write_blocks_to_file(params, data, file; global_ghosts, options...)
    end

    if !no_msg && params.is_root && params.silent < 2
        println("\nWrote to files $(output_file_path)_*x*")
    end
end


function write_time_step_file(params::ArmonParameters, file_name::String)
    file_path = build_file_path(params, file_name)

    p = params.output_precision
    format = Printf.Format("%#$(p+7).$(p)e\n")

    open(file_path, "w") do file
        Printf.format(file, format, params.curr_cycle_dt)
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

    for var in vars
        ref_var = getfield(ref_blk, var)
        our_var = getfield(our_blk, var)

        diff_mask = (!isapprox).(ref_var, our_var; rtol=params.comparison_tolerance)
        !params.write_ghosts && (diff_mask .*= is_ghost.(Ref(our_blk.size), 1:prod(block_size(our_blk))))

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

                val_diff = ref_var[idx] - our_var[idx]
                diff_ulp = val_diff / eps(ref_var[idx])
                abs(diff_ulp) > 1e10 && (diff_ulp = Inf)

                pos_str = join(map(i -> @sprintf("%3d"), I), ',')
                @printf("   - %5d (%s): %10.5g ≢ %10.5g (%11.5g, ulp: %8g)\n",
                    idx, pos_str, ref_var[idx], our_var[idx], val_diff, diff_ulp)
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
    for (ref_blk, our_blk) in zip(host_blocks(ref_data), host_blocks(our_data))
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


function step_checkpoint(params::ArmonParameters, grid::BlockGrid, step_label::String)
    !params.compare && return false

    wait(params)
    device_to_host!(grid)
    wait(params)

    step_file_name = params.output_file * @sprintf("_%03d_%s", params.cycle, step_label)
    step_file_name *= "_" * string(params.current_axis)[1]

    if params.is_ref
        if step_label == "time_step"
            write_time_step_file(params, step_file_name)
        else
            write_sub_domain_file(params, grid, step_file_name; no_msg=true)
        end

        return false
    else
        if step_label == "time_step"
            ref_dt = read_time_step_file(params, step_file_name)
            different = isapprox(ref_dt, params.curr_cycle_dt; rtol=params.comparison_tolerance)
            @printf("Time step difference: ref Δt = %.18f, Δt = %.18f, diff = %.18f\n",
                    ref_dt, params.curr_cycle_dt, ref_dt - params.curr_cycle_dt)
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
