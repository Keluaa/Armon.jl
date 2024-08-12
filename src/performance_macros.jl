
const use_std_lib_threads = @load_preference("use_std_lib_threads", false)
const use_fast_math = @load_preference("use_fast_math", true)


macro kernel_function(func)
    func = :(@inline $func)
    func = use_fast_math ? :(@fastmath $func) : func
    return esc(func)
end


macro threads(expr)
    if use_std_lib_threads
        return esc(quote
            Threads.@threads :static $(expr)
        end)
    else
        return esc(quote
            Armon.@batch $(expr)
        end)
    end
end


function make_threaded_loop(expr::Expr; choice=:dynamic)
    with = :(@inbounds Armon.@threads $(expr))
    without = :(@inbounds $(expr))

    if choice == :dynamic
        return quote
            if params.use_threading
                $(with)
            else
                $(without)
            end
        end
    elseif choice == :with
        return with
    elseif choice == :without
        return without
    else
        error("Unknown 'choice' value: $choice")
    end
end


"""
    @threaded(expr)

Allows to enable/disable multithreading of the loop depending on `params.use_threading`.

```julia
    @threaded for i = 1:n
        y[i] = log10(x[i]) + x[i]
    end
```
"""
macro threaded(expr)
    return esc(make_threaded_loop(expr))
end
