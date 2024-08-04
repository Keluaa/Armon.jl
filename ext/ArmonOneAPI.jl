module ArmonOneAPI

using Armon
using KernelAbstractions
import oneAPI
import oneAPI: oneAPIBackend


Armon.create_device(::Val{:oneAPI}) = oneAPIBackend()
Armon.device_array_type(::oneAPIBackend) = oneAPI.oneArray


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, Dim, <:oneAPIBackend}) where {Dim}
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
