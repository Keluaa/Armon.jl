
function compare_block(
    params::ArmonParameters, ref_blk::LocalTaskBlock, our_blk::LocalTaskBlock, label::String;
    vars=saved_vars
)
    different = false

    real_static_bsize = params.block_size .- 2*params.nghost
    blk_global_pos = params.N_origin .- 1 .+ (Tuple(our_blk.pos) .- 1) .* real_static_bsize

    var_names = var_arrays_names(ref_blk, vars)
    ref_vars = var_arrays(ref_blk, vars; on_device=false)
    our_vars = var_arrays(our_blk, vars; on_device=false)
    for (var, ref_var, our_var) in zip(var_names, ref_vars, our_vars)
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
    vars=saved_vars
)
    different = false
    for (ref_blk, our_blk) in zip(all_blocks(ref_data), all_blocks(our_data))
        different |= compare_block(params, ref_blk, our_blk, label; vars)
    end
    return different
end


function compare_with_file(params::ArmonParameters, grid::BlockGrid, file_name::String, label::String)
    ref_data = BlockGrid(params)
    read_sub_domain_file!(params, ref_data, file_name)
    different = compare_data(params, ref_data, grid, label)

    if params.use_MPI
        different = MPI.Allreduce(different, |, params.cart_comm)
    end

    return different
end


function write_time_step_file(params::ArmonParameters, state::SolverState, file_name::String)
    file_path = build_file_path(CSVSolverIO, file_name, params, state.cycle)

    p = 17  # enough for exact decimal representation for Float64
    format = Printf.Format("%#$(p+7).$(p)e\n")

    open(file_path, "w") do file
        Printf.format(file, format, state.global_dt.current_dt)
    end

    return
end


function read_time_step_file(params::ArmonParameters{T}, file_name::String) where {T}
    file_path = build_file_path(CSVSolverIO, file_name, params, state.cycle)
    open(file_path, "w") do file
        return parse(T, readchomp(file))
    end
end


function step_checkpoint(params::ArmonParameters, state::SolverState, grid::BlockGrid, step_label::String)
    !params.compare && return false

    wait(params)
    device_to_host!(grid)
    wait(params)

    if state.global_dt.cycle == 0 && step_label == "time_step"
        axis = Axis.X
    else
        axis = state.axis
    end
    step_file_name = params.output_file * "_" * string(axis)[1] * "_" * step_label

    if params.is_ref
        if step_label == "time_step"
            write_time_step_file(params, state, step_file_name)
        else
            write_sub_domain_file(params, grid, step_file_name)
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
            write_sub_domain_file(params, grid, step_file_name * "_diff")
            println("Difference file written to $(step_file_name)_diff")
        end

        return different
    end
end
