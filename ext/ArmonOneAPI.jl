module ArmonOneAPI

using Armon
using KernelAbstractions
import oneAPI
import oneAPI: oneAPIBackend


Armon.create_device(::Val{:oneAPI}) = oneAPIBackend()
Armon.device_array_type(::oneAPIBackend) = oneAPI.oneArray


mutable struct oneAPIThreadInfo <: Armon.ThreadInfo
    tid     :: Int
    driver  :: oneAPI.ZeDriver
    device  :: oneAPI.ZeDevice
    context :: oneAPI.ZeContext
    queue   :: oneAPI.ZeCommandQueue
end


function Armon.init_backend(params::ArmonParameters, ::oneAPIBackend; options...)
    device = oneAPI.device()
    driver = oneAPI.driver()
    context = oneAPI.context()
    for tid in 1:params.nthreads
        queue = oneAPI.ZeCommandQueue(context, device)
        params.threads_info[tid] = oneAPIThreadInfo(tid, driver, device, context, queue)
    end

    params.backend_options = Armon.EmptyParams()
    return options
end


function Armon.setup_task_for_device(params::ArmonParameters{<:Any, <:oneAPIBackend}, tid)
    # oneAPI.jl automatically creates a new stream when using its API in a new task.
    # We want to reuse the streams across tasks in order to track them easily, especially when
    # viewing GPU activity.
    thread_info::oneAPIThreadInfo = Armon.thread_info(params, tid)
    oneAPI.driver!(thread_info.driver)
    oneAPI.device!(thread_info.device)
    oneAPI.context!(thread_info.context)
    # TODO: oneAPI.jl does not provide a `global_queue!`, most likely because it is a global queue,
    # not a local queue... does this mean using multiple queues per device isn't supported?
    task_local_storage((:ZeCommandQueue, thread_info.context, thread_info.device), thread_info.queue)
    return
end


function Base.wait(params::ArmonParameters{<:Any, <:oneAPIBackend}, tid)
    thread_info::oneAPIThreadInfo = Armon.thread_info(params, tid)
    oneAPI.oneL0.synchronize(thread_info.queue)
    return
end


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, <:oneAPIBackend})
    Armon.print_parameter(io, pad, "GPU", true, nl=false)
    println(io, ": oneAPI (block size: ", join(p.block_size, '×'), ")")
end


function Armon.device_memory_info(::oneAPIBackend)
    # TODO: I think it might be impossible to know how much memory is free, but total memory?
    return (
        total = UInt64(0),
        free  = UInt64(0)
    )
end

end
