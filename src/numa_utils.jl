
"""
    array_pages(A::Ptr, A_length)
    array_pages(A::DenseArray)

Iterator over the memory pages (aligned to `PAGE_SIZE`) used by the array `A`.
"""
function array_pages(A::Ptr{T}, A_length) where {T}
    page_size = NUMA.numa_pagesize()
    A_end = A + (A_length - 1) * sizeof(T)
    A -= UInt(A) % page_size  # Align to the first page
    return A:page_size:A_end
end

array_pages(A::DenseArray) = array_pages(pointer(A), length(A))


function move_pages_status_to_str(status)
    # See https://linux.die.net/man/2/move_pages
    if     0 ≤ status              return "ok"
    elseif status == -Libc.EACCES  return "EACCES"
    elseif status == -Libc.EBUSY   return "EBUSY"
    elseif status == -Libc.EFAULT  return "EFAULT"
    elseif status == -Libc.EIO     return "EIO"
    elseif status == -Libc.EINVAL  return "EINVAL"
    elseif status == -Libc.ENOENT  return "ENOENT"
    elseif status == -Libc.ENOMEM  return "ENOMEM"
    else                           return "unknown page status code: $status"
    end
end


"""
    move_pages(pages::Vector{Ptr}, target_node; retries=4)
    move_pages(pages::OrdinalRange{Ptr, Int}, target_node)
    move_pages(A::DenseArray, target_node)

Move all pages (of either the array `A` or from the `pages` iterable) to the `target_node` (1-index).

The move status of all pages is checked afterwards.
No error is thrown only if all pages are on `target_node` at the end of the move.

In case some pages cannot be moved because of EBUSY status, the move is retried `retries`-times
until it succeeds, with a short `sleep` in-between attempts.
"""
function move_pages(pages::Vector{Ptr{T}}, target_node; retries=4) where {T}
    target_nodes = similar(pages, Cint)
    fill!(target_nodes, Cint(target_node - 1))
    page_statuses = similar(pages, Cint)
    ret = NUMA.LibNuma.move_pages(0, length(pages), pages, target_nodes, page_statuses, NUMA.LibNuma.MPOL_MF_MOVE)
    if ret != 0
        Libc.errno() == Libc.ENOENT && return  # all pages already on the target node
        error("move_pages failed: $(Libc.errno())")
    end

    move_fail = count(page_statuses .!= Cint(target_node - 1))
    if move_fail > 0
        busy_pages_idx = findall(==(Cint(-Libc.EBUSY)), page_statuses)
        if length(busy_pages_idx) == move_fail && retries > 0
            # For EBUSY failures, waiting for a bit fixes most issues
            sleep(0.050)  # 50ms
            move_pages(pages[busy_pages_idx], target_node; retries=retries-1)
            move_fail -= length(busy_pages_idx)
        end

        if move_fail > 0
            probl_pages_idx = findall(!=(Cint(target_node - 1)), page_statuses)
            probl_pages_idx = probl_pages_idx[1:min(50, end)]  # limit the error output to 50 pages
            pages_ptr = pages[probl_pages_idx]
            pages_reason = move_pages_status_to_str.(page_statuses[probl_pages_idx])
            probl_strs = Ref(" - ") .* string.(pages_ptr) .* Ref(" : ") .* pages_reason
            move_fail > length(probl_pages_idx) && push!(probl_strs, "   (too many pages to display)")
            error("could not move $move_fail of the $(length(pages)) pages:\n$(join(probl_strs, '\n'))")
        end
    end

    return
end

move_pages(pages::OrdinalRange{Ptr{T}, Int}, target_node) where {T} = move_pages(collect(pages), target_node)
move_pages(A::DenseArray, target_node) = move_pages(array_pages(A), target_node)


"""
    lock_pages(ptr, len)
    lock_pages(pages::OrdinalRange{Ptr, Int})
    lock_pages(A::DenseArray)

Lock the `pages` ranging from `ptr` to `ptr+len` on the RAM, preventing the kernel from moving it
around or putting it in swap memory.
"""
function lock_pages(ptr::Ptr, len)
    ret = ccall(:mlock, Cint, (Ptr{Cvoid}, Csize_t), Ptr{Cvoid}(ptr), len)
    ret != 0 && error("mlock failed: $(Libc.errno())")
    return
end

lock_pages(pages::OrdinalRange{Ptr{T}, Int}) where {T} = lock_pages(first(pages), UInt(last(pages) - first(pages)))
lock_pages(A::DenseArray) = lock_pages(array_pages(A))


"""
    unlock_pages(ptr, len)
    unlock_pages(pages::OrdinalRange{Ptr, Int})
    unlock_pages(A::DenseArray)

Unlock the `pages`, opposite of [`lock_pages`](@ref).
"""
function unlock_pages(ptr, len)
    ret = ccall(:munlock, Cint, (Ptr{Cvoid}, Csize_t), Ptr{Cvoid}(ptr), len)
    ret != 0 && error("munlock failed: $(Libc.errno())")
    return
end

unlock_pages(pages::OrdinalRange{Ptr{T}, Int}) where {T} = unlock_pages(first(pages), UInt(last(pages) - first(pages)))
unlock_pages(A::DenseArray) = unlock_pages(array_pages(A))


function volatile_touch(ptr::Ptr)
    # We simply want to touch the memory, without creating any side effects.
    # To do this one would simply make a load followed by a store at the same location, however this
    # would be optimized away by the compiler. The best way to prevent this is to use volatile
    # operations, which cannot be optimized by the compiler.
    # There is no `volatile` keyword in Julia, therefore we must rely on inline LLVM IR.
    Core.Intrinsics.llvmcall("""
        %ptr = inttoptr i64 %0 to i64*
        %val = load volatile i64, i64* %ptr
        store volatile i64 %val, i64* %ptr
        ret void
    """, Cvoid, Tuple{Ptr{Int64}}, Ptr{Int64}(ptr))
end


"""
    touch_pages(pages)

Parses over each of `pages` and "touch" them by loading a value then storing it again.

The (only) use case is to manually trigger the first-touch policy and/or force the kernel to
physically allocate all pages.
"""
touch_pages(pages) = foreach(volatile_touch, pages)


function tid_to_numa_node_map()
    # Julia thread ID to its NUMA node ID
    numa_map = Vector{Int}(undef, Threads.nthreads())
    Threads.@threads :static for _ in 1:Threads.nthreads()
        tid = Threads.threadid()
        numa_map[tid] = NUMA.current_numa_node()
    end
    return numa_map
end
