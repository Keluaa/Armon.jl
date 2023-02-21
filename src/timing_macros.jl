
in_warmup_cycle = false
set_warmup(b::Bool) = global in_warmup_cycle = b
is_warming_up()::Bool = in_warmup_cycle

#
# Timings storage
#

axis_time_contrib = Dict{Axis, Dict{String, Float64}}()
total_time_contrib = Dict{String, Float64}()
const time_contrib_lock = ReentrantLock()


function add_common_time(label, time)
    lock(time_contrib_lock) do
        if !haskey(total_time_contrib, label)
            global total_time_contrib[label] = time
        else
            global total_time_contrib[label] += time
        end
    end
end


function add_axis_time(axis, label, time)
    lock(time_contrib_lock) do
        if !haskey(axis_time_contrib, axis)
            global axis_time_contrib[axis] = Dict{String, Float64}()
        end

        if !haskey(axis_time_contrib[axis], label)
            global axis_time_contrib[axis][label] = time
        else
            global axis_time_contrib[axis][label] += time
        end
    end
end


#
# Hardware counters storage
#

axis_hw_counters = nothing
const hw_counters_lock = ReentrantLock()


function init_hw_counter(axis, label)
    inner_dict_type = Dict{LinuxPerf.EventType, LinuxPerf.Counter}
    lock(hw_counters_lock) do
        if isnothing(axis_hw_counters)
            global axis_hw_counters = Dict{Axis, Dict{String, inner_dict_type}}()
        end
        if !haskey(axis_hw_counters, axis)
            global axis_hw_counters[axis] = Dict{String, inner_dict_type}()
        end
        if !haskey(axis_hw_counters[axis], label)
            global axis_hw_counters[axis][label] = inner_dict_type()
        end
    end
    return
end

#
# Macros 
#

function build_time_expr(label, common_time_only, expr; use_wait=true, exclude_from_total=false)
    return esc(quote
        if params.measure_time && !is_warming_up()
            if params.measure_hw_counters
                var"_$(label)_hw" = init_bench(params.hw_counters_options)
                enable_bench(var"_$(label)_hw")
            end

            @static if $(use_wait)
                var"_$(label)_res"   = $(expr)
                var"_$(label)_time"  = @elapsed wait(var"_$(label)_res")
                var"_$(label)_time" *= 1e9
            else
                var"_$(label)_start" = time_ns()
                var"_$(label)_res"   = $(expr)
                var"_$(label)_end"   = time_ns()
                var"_$(label)_time"  = var"_$(label)_end" - var"_$(label)_start"
            end

            if params.measure_hw_counters
                disable_bench(var"_$(label)_hw")
                var"_$(label)_hw_stats" = close_bench(var"_$(label)_hw") |> reduce_threads_stats!
                sum_stats!(axis_hw_counters[$(common_time_only) ? X_axis : params.current_axis][$(label)], var"_$(label)_hw_stats")
            end

            @static if $(!common_time_only)
                add_axis_time(params.current_axis, $(label), var"_$(label)_time")
            end
            @static if $(!exclude_from_total)
                add_common_time($(label), var"_$(label)_time")
            end
            
            var"_$(label)_res"
        elseif params.measure_hw_counters && is_warming_up()
            init_hw_counter(params.current_axis, $(label))
            $(expr)
        else
            $(expr)
        end
    end)
end


function extract_function_name(expr)
    if expr.head == :call
        function_name = expr.args[1]
    elseif isa(expr.args[2], Expr) && expr.args[2].head == :call
        function_name = expr.args[2].args[1]
    else
        error("Could not find the function name of the provided expression")
    end
    return string(function_name)
end


macro time_event(label, expr)   return build_time_expr(label, false, expr) end
macro time_event_c(label, expr) return build_time_expr(label, true,  expr) end
macro time_expr(label, expr)    return build_time_expr(label, false, expr; use_wait=false) end
macro time_expr_c(label, expr)  return build_time_expr(label, true,  expr; use_wait=false) end
macro time_event_a(label, expr) return build_time_expr(label, false,  expr; exclude_from_total=true) end
macro time_expr_a(label, expr)  return build_time_expr(label, false,  expr; exclude_from_total=true, use_wait=false) end
macro time_event(expr)          return build_time_expr(extract_function_name(expr), false, expr) end
macro time_event_c(expr)        return build_time_expr(extract_function_name(expr), true,  expr) end
macro time_expr(expr)           return build_time_expr(extract_function_name(expr), false, expr; use_wait=false) end
macro time_expr_c(expr)         return build_time_expr(extract_function_name(expr), true,  expr; use_wait=false) end
macro time_event_a(expr)        return build_time_expr(extract_function_name(expr), false, expr; exclude_from_total=true) end
macro time_expr_a(expr)         return build_time_expr(extract_function_name(expr), false, expr; exclude_from_total=true, use_wait=false) end


@static if VERSION >= v"1.8"
    compile_time_before_ns() = Base.cumulative_compile_time_ns()[1]
    compile_time_after_ns() = Base.cumulative_compile_time_ns()[1]
else
    compile_time_before_ns() = Base.cumulative_compile_time_ns_before()
    compile_time_after_ns() = Base.cumulative_compile_time_ns_after()
end


# Equivalent to `@time` but with a better output
macro pretty_time(expr)
    function_name = extract_function_name(expr)
    return esc(quote
        if params.is_root && params.silent <= 3
            # Same structure as `@time` (see `@macroexpand @time`), using some undocumented functions.
            gc_info_before = Base.gc_num()
            time_before = Base.time_ns()
            compile_time_before = compile_time_before_ns()

            $(expr)

            compile_time_after = compile_time_after_ns()
            time_after = Base.time_ns()
            gc_info_after = Base.gc_num()
            
            elapsed_time = time_after - time_before
            compile_time = compile_time_after - compile_time_before
            gc_diff = Base.GC_Diff(gc_info_after, gc_info_before)

            allocations_size = gc_diff.allocd / 1e3
            allocations_count = Base.gc_alloc_count(gc_diff)
            gc_time = gc_diff.total_time
    
            println("\nTime info for $($(function_name)):")
            if allocations_count > 0
                @printf(" - %d allocations for %g kB\n", 
                    allocations_count, convert(Float64, allocations_size))
            end
            if gc_time > 0
                @printf(" - GC:      %10.5f ms (%5.2f%%)\n", 
                    gc_time / 1e6, gc_time / elapsed_time * 100)
            end
            if compile_time > 0
                @printf(" - Compile: %10.5f ms (%5.2f%%)\n", 
                    compile_time / 1e6, compile_time / elapsed_time * 100)
            end
            @printf(" - Total:   %10.5f ms\n", elapsed_time / 1e6)
        else
            $(expr)
        end
    end)
end
