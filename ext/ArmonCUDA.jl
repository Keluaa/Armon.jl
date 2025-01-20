module ArmonCUDA

using Armon
using KernelAbstractions
import CUDA
import CUDA: CUDABackend


Armon.create_device(::Val{:CUDA}) = CUDABackend()
Armon.device_array_type(::CUDABackend) = CUDA.CuArray
Armon.device_converter(::CUDABackend) = CUDA.KernelAdaptor()


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, <:CUDABackend})
    Armon.print_parameter(io, pad, "GPU", true, nl=false)
    device = first(p.threads_info).device
    println(io, ": CUDA, workgroup size: ", join(p.workgroup_size, '×'), ", ",
        p.nthreads, " streams, device n°", CUDA.deviceid(device))
end


function Armon.device_memory_info(::CUDABackend)
    free, total = CUDA.Mem.info()
    return (
        total = UInt64(total),
        free  = UInt64(free)
    )
end


mutable struct CuThreadInfo <: Armon.ThreadInfo
    tid    :: Int
    device :: CUDA.CuDevice
    stream :: CUDA.CuStream
end


function Armon.init_backend(params::ArmonParameters, ::CUDABackend; options...)
    armon_nvtx = Base.get_extension(Armon, :ArmonNVTX)

    device = CUDA.device()  # TODO: would this allow us to use multiple devices from the same process?
    for tid in 1:params.nthreads
        stream = CUDA.create_stream()
        if !isnothing(armon_nvtx)
            armon_nvtx.name_stream(stream, "Armon stream $tid")
        end
        params.threads_info[tid] = CuThreadInfo(tid, device, stream)
    end

    params.backend_options = Armon.EmptyParams()
    return options
end


function Armon.setup_task_for_device(params::ArmonParameters{<:Any, <:CUDABackend}, tid)
    # CUDA.jl automatically creates a new stream when using its API in a new task.
    # We want to reuse the streams across tasks in order to track them easily, especially when
    # viewing GPU activity in Nsight Systems.
    thread_info::CuThreadInfo = Armon.thread_info(params, tid)
    CUDA.device!(thread_info.device)
    CUDA.stream!(thread_info.stream)
    return
end


Armon.create_kernel_event(::CUDABackend) = CUDA.CuEvent(CUDA.EVENT_DISABLE_TIMING)
Armon.put_kernel_event(::CUDABackend, event) = CUDA.record(event)
Armon.query_kernel_event(::CUDABackend, event) = CUDA.isdone(event)


function Base.wait(params::ArmonParameters{<:Any, <:CUDABackend}, tid)
    thread_info::CuThreadInfo = Armon.thread_info(params, tid)
    CUDA.synchronize(thread_info.stream)
    return
end


function Armon.lock_pages(::CUDABackend, ptr::Ptr, len)
    len == 0 && return
    # TODO: passing the `CUDA.MEMHOSTREGISTER_DEVICEMAP` flag would allow the GPU to access the host
    # memory from the GPU seamlessly.
    flags = 0
    CUDA.register(CUDA.HostMemory, ptr, len, flags)
    return
end


function Armon.unlock_pages(::CUDABackend, ptr::Ptr, len)
    len == 0 && return
    # TODO: passing the `CUDA.MEMHOSTREGISTER_DEVICEMAP` flag would allow the GPU to access the host
    # memory from the GPU seamlessly.
    CUDA.unregister(CUDA.HostMemory(CUDA.context(), ptr, len))
    return
end


function Base.copyto!(::ArmonParameters{<:Any, <:CUDABackend}, dst::AbstractArray, src::AbstractArray)
    # Asynchronous copy using the current active CuStream. This requires that host arrays to be
    # pinned.
    # The implementation is similar to `CUDAKernels.copyto!`, but without the expensive `CUDA.pin`
    # on host arrays.
    GC.@preserve dst src begin
        dst_ptr = pointer(dst_var)
        src_ptr = pointer(src_var)
        unsafe_copyto!(dst_ptr, src_ptr, length(dst); async=true)
    end
end


function cuda_kernel_start(_, _)
    # Equivalent to CUDA.@profile
    CUDA.Profile.start()
    return nothing
end


function cuda_kernel_end(_, _, _)
    CUDA.Profile.stop()
end


function __init__()
    Armon.register_kernel_callback(Armon.KernelCallback((
        :CUDA_kernels,
        cuda_kernel_start,
        cuda_kernel_end
    )))
end


# TODO: this is temporary type-piracy in order to support `@atomicswap`
# See this issue: https://github.com/JuliaConcurrent/Atomix.jl/issues/49
@inline function Armon.Atomix.modify!(ref::Armon.Atomix.IndexableRef{<:CUDA.CuDeviceArray}, op::OP, x, order) where {OP}
    x = convert(eltype(ref), x)
    ptr = Armon.Atomix.pointer(ref)
    begin
        old = if op === (+)
            CUDA.atomic_add!(ptr, x)
        elseif op === (-)
            CUDA.atomic_sub!(ptr, x)
        elseif op === (&)
            CUDA.atomic_and!(ptr, x)
        elseif op === (|)
            CUDA.atomic_or!(ptr, x)
        elseif op === xor
            CUDA.atomic_xor!(ptr, x)
        elseif op === min
            CUDA.atomic_min!(ptr, x)
        elseif op === max
            CUDA.atomic_max!(ptr, x)
        elseif op === Armon.Atomix.right
            CUDA.atomic_xchg!(ptr, x)
        else
            error("not implemented")
        end
    end
    return old => op(old, x)
end

end
