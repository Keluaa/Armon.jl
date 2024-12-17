module ArmonAMDGPU

using Armon
using KernelAbstractions
import AMDGPU
import AMDGPU: ROCBackend


Armon.create_device(::Val{:ROCM}) = ROCBackend()
Armon.device_array_type(::ROCBackend) = AMDGPU.ROCArray


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, Dim, <:ROCBackend}) where {Dim}
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
