module ArmonAMDGPU

using Armon
using KernelAbstractions
import AMDGPU
import AMDGPU: ROCBackend


Armon.create_device(::Val{:ROCM}) = ROCBackend()
Armon.device_array_type(::ROCBackend) = AMDGPU.ROCArray


mutable struct ROCThreadInfo <: Armon.ThreadInfo
    tid    :: Int
    device :: AMDGPU.HIPDevice
    stream :: AMDGPU.HIPStream
end


function Armon.init_backend(params::ArmonParameters, ::ROCBackend; options...)
    device = AMDGPU.device()
    for tid in 1:params.nthreads
        stream = AMDGPU.HIPStream()
        params.threads_info[tid] = ROCThreadInfo(tid, device, stream)
    end

    params.backend_options = Armon.EmptyParams()
    return options
end


function Armon.setup_task_for_device(params::ArmonParameters{<:Any, <:ROCBackend}, tid)
    # AMDGPU.jl automatically creates a new stream when using its API in a new task.
    # We want to reuse the streams across tasks in order to track them easily, especially when
    # viewing GPU activity.
    thread_info::ROCThreadInfo = Armon.thread_info(params, tid)
    AMDGPU.device!(thread_info.device)
    AMDGPU.stream!(thread_info.stream)
    return
end


function Base.wait(params::ArmonParameters{<:Any, <:ROCBackend}, tid)
    thread_info::ROCThreadInfo = Armon.thread_info(params, tid)
    AMDGPU.synchronize(thread_info.stream)
    return
end


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, <:ROCBackend})
    Armon.print_parameter(io, pad, "GPU", true, nl=false)
    println(io, ": ROCm (block size: ", join(p.block_size, '×'), ")")
end


function Armon.device_memory_info(::ROCBackend)
    @static if pkgversion(AMDGPU) ≥ v"0.5"
        free, total = AMDGPU.Runtime.Mem.info()
    else
        free_p = Ref{UInt64}()
        total_p = Ref{UInt64}()
        ccall((:hipMemGetInfo, AMDGPU.libhip),
            AMDGPU.HIP.hipError_t, (Ptr{Csize_t}, Ptr{Csize_t}),
            free_p, total_p)
        free = free_p[]
        total = total_p[]
    end

    return (
        total = UInt64(total),
        free  = UInt64(free)
    )
end


# TODO: profiling with rocprofile with roctracer library

end
