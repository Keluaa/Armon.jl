module ArmonCUDA

using Armon
using KernelAbstractions
import CUDA
import CUDA: CUDABackend


Armon.create_device(::Val{:CUDA}) = CUDABackend()
Armon.device_array_type(::CUDABackend) = CUDA.CuArray


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, <:CUDABackend})
    Armon.print_parameter(io, pad, "GPU", true, nl=false)
    println(io, ": CUDA (block size: ", join(p.block_size, '×'), ")")
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
    device :: CuDevice
    stream :: CuStream
end


function Armon.init_backend(params::ArmonParameters, ::CUDABackend; options...)
    device = CUDA.device()  # TODO: would this allow us to use multiple devices from the same process?
    for tid in 1:params.nthreads
        stream = CUDA.create_stream()  # TODO: name the stream with NVTX.jl : `NVTX.nvtxNameCuStreamA`
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


function Base.wait(params::ArmonParameters{<:Any, <:CUDABackend}, tid)
    thread_info::CuThreadInfo = Armon.thread_info(params, tid)
    CUDA.synchronize(thread_info.stream)
    return
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

end
