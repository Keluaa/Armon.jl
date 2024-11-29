
"""
    create_device(::Val{:device_name})

Create a device object from its name.

Default devices:
 - `:CPU`: the CPU backend of `KernelAbstractions.jl`
 - `:CPU_HP`: low overhead multithreading with `Polyester.jl`

Extensions:
 - `:Kokkos`: the default `Kokkos.jl` device
 - `:CUDA`: the `CUDA.jl` backend of `KernelAbstractions.jl`
 - `:ROCM`: the `AMDGPU.jl` backend of `KernelAbstractions.jl`
 - `:oneAPI`: the `oneAPI.jl` backend of `KernelAbstractions.jl`
"""
create_device(::Val{:CPU}) = CPU()
create_device(::Val{:CPU_HP}) = CPU_HP()


struct CPUThreadInfo <: ThreadInfo
    tid :: Int
end

thread_info(params::ArmonParameters, tid = Threads.threadid()) = params.threads_info[tid]


"""
    setup_task_for_device(params::ArmonParameters, tid)

Since each multi-threaded section may create a new task, the task local info of GPU backends is
reset.

GPU backends should therefore set the right device and stream to use in this function, ensuring the
same stream is consistently used by the same thread throughout all solver iterations.
"""
function setup_task_for_device(::ArmonParameters{<:Any, <:Union{CPU, CPU_HP}}, tid)
    # nothing to do on CPU
end


"""
    Base.wait(params::ArmonParameters, tid)

Wait for the completion of all kernels launched in the stream assigned to the thread `tid`.

This operation is blocking. Most GPU backends use a mix of spin-loops and GPU API calls to perform
the synchronization as fast as possible.

To synchronize the whole device, use [`Base.wait(::ArmonParameters)`](@ref).
"""
function Base.wait(::ArmonParameters{<:Any, <:Union{CPU, CPU_HP}}, tid)
    # CPU backends are synchronous
end


"""
    Base.wait(params::ArmonParameters)

Wait for the completion of all kernels launched on the `params.device`.

This operation is blocking.
"""
function Base.wait(::ArmonParameters{<:Any, <:Union{CPU, CPU_HP}})
    # CPU backends are synchronous
end

function Base.wait(params::ArmonParameters{<:Any, <:GPU})
    KernelAbstractions.synchronize(params.device)
end


Base.copyto!(::ArmonParameters{<:Any, <:Union{CPU, CPU_HP}}, dst, src) = copyto!(dst, src)


"""
    memory_info(params)

The total and free memory the current process can store on the `params.device`.
"""
function memory_info(params::ArmonParameters)
    mem_info = device_memory_info(params.device)
    # TODO: MPI support
    return mem_info
end


"""
    device_memory_info(device)

The total and free memory on the device, in bytes.
"""
function device_memory_info(::Union{CPU_HP, CPU})
    return (
        total = UInt64(Sys.total_physical_memory()),
        free  = UInt64(Sys.free_physical_memory())
    )
end
