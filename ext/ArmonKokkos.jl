module ArmonKokkos

using Armon
import Armon: ArmonData, Side, solver_error
import Armon: NoLimiter, MinmodLimiter, SuperbeeLimiter
import Armon: Sod, Sod_y, Sod_circ, Bizarrium, Sedov

isdefined(Base, :get_extension) ? (import Kokkos) : (import ..Kokkos)
import Kokkos: CMakeKokkosProject, option!
using MPI


"Mirror of `Range`, in 'src/kernels/indexing.h'"
mutable struct Range{Idx <: Integer}
    start::Idx
    var"end"::Idx  # exclusive
end


"Mirror of `InnerRange1D`, in 'src/kernels/indexing.h'"
mutable struct InnerRange1D{Idx <: Integer}
    start::Idx
    step::Idx
end


"Mirror of `InnerRange2D`, in 'src/kernels/indexing.h'"
mutable struct InnerRange2D{Idx <: Integer, UIdx <: Integer}
    main_range_start::Idx
    main_range_step::Idx
    row_range_start::Idx
    row_range_length::UIdx
end


struct ArmonKokkosParams{Idx <: Integer, UIdx <: Integer} <: Armon.BackendParams
    project::CMakeKokkosProject
    lib::Kokkos.CLibrary
    limiter_index::Cint
    test_case_index::Cint
    range::Range{Idx}
    inner_range_1D::InnerRange1D{Idx}
    inner_range_2D::InnerRange2D{Idx, UIdx}
end


function Armon.create_device(::Val{:Kokkos})
    !Kokkos.is_initialized() && solver_error(:config, "Kokkos should be initialized before creating ArmonParameters")
    return Base.invokelatest(Kokkos.DEFAULT_DEVICE_SPACE)
end


function limiter_type_to_int(params::ArmonParameters)
    if     params.riemann_limiter isa NoLimiter       return 0
    elseif params.riemann_limiter isa MinmodLimiter   return 1
    elseif params.riemann_limiter isa SuperbeeLimiter return 2
    else
        solver_error(:config, "This limiter is not recognized by armon_cpp")
    end
end


function test_case_to_int(params::ArmonParameters)
    if     params.test isa Sod       return 0
    elseif params.test isa Sod_y     return 1
    elseif params.test isa Sod_circ  return 2
    elseif params.test isa Bizarrium return 3
    elseif params.test isa Sedov     return 4
    else
        solver_error(:config, "This test case is not recognized by armon_cpp")
    end
end


function raise_cpp_exception(kernel::Cstring, msg::Cstring)
    kernel_str = Base.unsafe_string(kernel)
    msg_str = Base.unsafe_string(msg)
    return solver_error(:cpp, "C++ exception in kernel $kernel_str: $msg_str")
end


const Idx_type = Ref{DataType}(Cint)
get_idx_type() = Idx_type[]


function Armon.init_backend(params::ArmonParameters, ::Kokkos.ExecutionSpace;
    armon_cpp_lib_src = nothing, cmake_options = [], kokkos_options = nothing,
    debug_kernels = false, use_md_iter = false,
    options...
)
    !Kokkos.is_initialized() && solver_error(:config, "Kokkos has not yet been initialized")

    # Assume we are in the ArmonBenchmark project
    armon_cpp_lib_src = @something armon_cpp_lib_src joinpath(@__DIR__, "..", "..", "kokkos")
    !isdir(armon_cpp_lib_src) && solver_error(:config, "Invalid path to the Armon C++ library source: $armon_cpp_lib_src")

    armon_cpp_lib_build = joinpath(Kokkos.KOKKOS_BUILD_DIR, "armon_kokkos")
    armon_cpp = CMakeKokkosProject(armon_cpp_lib_src, "src/kernels/libarmon_kernels";
        target="armon_kernels", build_dir=armon_cpp_lib_build, cmake_options, kokkos_options)
    option!(armon_cpp, "ENABLE_DEBUG_BOUNDS_CHECK", debug_kernels)
    option!(armon_cpp, "USE_SINGLE_PRECISION", Armon.data_type(params) == Float32; prefix="")
    option!(armon_cpp, "TRY_ALL_CALLS", debug_kernels; prefix="")
    option!(armon_cpp, "CHECK_VIEW_ORDER", debug_kernels; prefix="")
    option!(armon_cpp, "USE_SIMD_KERNELS", params.use_simd; prefix="")
    option!(armon_cpp, "USE_MD_ITER", use_md_iter; prefix="")
    is_NVTX_loaded = !isnothing(Base.get_extension(Armon, :ArmonNVTX))
    option!(armon_cpp, "USE_NVTX", is_NVTX_loaded; prefix="")

    if params.use_MPI
        # Prevent concurrent compilation
        params.is_root && Kokkos.compile(armon_cpp)
        MPI.Barrier(params.global_comm)
    else
        Kokkos.compile(armon_cpp)
    end
    kokkos_lib = Kokkos.load_lib(armon_cpp)

    # Initialize the library and its wrapper

    cpp_exception_handler = cglobal(Kokkos.get_symbol(kokkos_lib, :raise_exception_handler), Ptr{Ptr{Cvoid}})
    unsafe_store!(cpp_exception_handler, @cfunction(raise_cpp_exception, Cvoid, (Cstring, Cstring)))

    cpp_flt_size = ccall(Kokkos.get_symbol(kokkos_lib, :flt_size), Cint, ())
    T = Armon.data_type(params)
    if cpp_flt_size != sizeof(T)
        solver_error(:config, "eltype size mismatch: expected $(sizeof(T)) bytes (for $T), got $cpp_flt_size bytes")
    end

    cpp_idx_size = ccall(Kokkos.get_symbol(kokkos_lib, :idx_size), Cint, ())
    cpp_is_uidx_signed = ccall(Kokkos.get_symbol(kokkos_lib, :is_uidx_signed), Cuchar, ()) |> Bool

    if cpp_idx_size == 8
        Idx = Int64
    elseif cpp_idx_size == 4
        Idx = Int32
    else
        solver_error(:config, "unknown index type size: $cpp_idx_size bytes")
    end

    UIdx = cpp_is_uidx_signed ? Idx : unsigned(Idx)
    Idx_type[] = Idx

    params.backend_options = ArmonKokkosParams{Idx, UIdx}(
        armon_cpp, kokkos_lib,
        limiter_type_to_int(params), test_case_to_int(params),
        Range{Idx}(zero(Idx), zero(Idx)),
        InnerRange1D{Idx}(zero(Idx), zero(Idx)),
        InnerRange2D{Idx, UIdx}(zero(Idx), zero(Idx), zero(Idx), zero(UIdx))
    )

    return options
end


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, <:Kokkos.ExecutionSpace})
    Armon.print_parameter(io, pad, "use_kokkos", true)
    Armon.print_parameter(io, pad, "device", nameof(Kokkos.main_space_type(p.device)))
    Armon.print_parameter(io, pad, "memory", nameof(Kokkos.main_space_type(Kokkos.memory_space(p.device))))
end


function Armon.device_memory_info(exec::Kokkos.ExecutionSpace)
    if exec isa Kokkos.Cuda || exec isa Kokkos.HIP
        free, total = Kokkos.BackendFunctions.memory_info()
        return (
            total = UInt64(total),
            free  = UInt64(free)
        )
    elseif exec isa Kokkos.Serial || exec isa Kokkos.OpenMP
        return Armon.device_memory_info(Armon.CPU_HP())
    else
        error("`device_memory_info` for $(Kokkos.main_space_type(exec)) NYI")
    end
end


function Base.wait(::ArmonParameters{<:Any, <:Kokkos.ExecutionSpace})
    Kokkos.fence()
end

#
# Array allocation
#

function Armon.host_array_type(::Kokkos.ExecutionSpace)
    return Kokkos.View{T, D,
        Kokkos.array_layout(Kokkos.DEFAULT_HOST_SPACE),
        Kokkos.DEFAULT_HOST_MEM_SPACE
    } where {T, D}
end


function Armon.device_array_type(device::Kokkos.ExecutionSpace)
    return Kokkos.View{T, D,
        Kokkos.array_layout(device),
        Kokkos.memory_space(device)
    } where {T, D}
end

#
# Custom reduction kernels
#

@generated function Armon.dtCFL_kernel(
    params::ArmonParameters{T, <:Kokkos.ExecutionSpace}, data::ArmonData{V}, range, dx, dy
) where {T, V <: AbstractArray{T}}
    quote
        cpp_range = params.backend_options.range
        cpp_range.start = 0
        cpp_range.end = length(range)

        inner_range_1D = params.backend_options.inner_range_1D
        inner_range_1D.start = first(range) - 1  # to 0-index
        inner_range_1D.step = 1

        return ccall(Kokkos.get_symbol(params.backend_options.lib, :dt_CFL),
            T, (Ptr{Cvoid}, Ptr{Cvoid},
                T, T, Ref{V}, Ref{V}, Ref{V}, Ref{V}),
            pointer_from_objref(cpp_range), pointer_from_objref(inner_range_1D),
            dx, dy, data.umat, data.vmat, data.cmat, data.domain_mask,
        )
    end
end


@generated function Armon.conservation_vars_kernel(
    params::ArmonParameters{T, <:Kokkos.ExecutionSpace}, data::ArmonData{V}, range
) where {T, V <: Kokkos.View{T}}
    quote
        cpp_range = params.backend_options.range
        cpp_range.start = 0
        cpp_range.end = length(range)

        inner_range_1D = params.backend_options.inner_range_1D
        inner_range_1D.start = first(range) - 1  # to 0-index
        inner_range_1D.step = 1

        total_mass_ref = Ref{T}(0)
        total_energy_ref = Ref{T}(0)
        ccall(Kokkos.get_symbol(params.backend_options.lib, :conservation_vars),
            Cvoid, (Ptr{Cvoid}, Ptr{Cvoid},
                    T, Ref{V}, Ref{V}, Ref{V},
                    Ref{T}, Ref{T}),
            pointer_from_objref(cpp_range), pointer_from_objref(inner_range_1D),
            params.dx, data.rho, data.Emat, data.domain_mask,
            total_mass_ref, total_energy_ref
        )
        return total_mass_ref[], total_energy_ref[]
    end
end

#
# Copies
#

function Armon.copy_to_send_buffer!(data::ArmonDualData{D, H, <:Kokkos.ExecutionSpace},
    array::D, buffer::H
) where {D, H}
    array_data = Kokkos.subview(array, 1:length(buffer))
    Kokkos.deep_copy(Armon.device_type(data), buffer, array_data)
end


function Armon.copy_from_recv_buffer!(data::ArmonDualData{D, H, <:Kokkos.ExecutionSpace},
    array::D, buffer::H
) where {D, H}
    array_data = Kokkos.subview(array, 1:length(buffer))
    Kokkos.deep_copy(Armon.device_type(data), array_data, buffer)
end


function Armon.get_send_comm_array(data::ArmonDualData{H, H, <:Kokkos.ExecutionSpace},
    side::Side
) where H
    # Here only for method disambiguation
    return Base.invoke(Armon.get_send_comm_array, Tuple{ArmonDualData{H, H, <:Any}, Side}, data, side)
end


function Armon.get_send_comm_array(data::ArmonDualData{D, H, <:Kokkos.ExecutionSpace},
    side::Side
) where {D, H}
    # Kokkos functions (like deep_copy, etc...) called on the comm array require an actual
    # Kokkos.View, not a view.
    array_view = Base.invoke(Armon.get_send_comm_array, Tuple{ArmonDualData{D, H, <:Any}, Side}, data, side)
    return Kokkos.subview(parent(array_view), parentindices(array_view))
end

end
