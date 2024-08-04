module ArmonKokkos

using Armon
import Armon: solver_error
import Kokkos
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
    if     params.riemann_limiter isa Armon.NoLimiter       return 0
    elseif params.riemann_limiter isa Armon.MinmodLimiter   return 1
    elseif params.riemann_limiter isa Armon.SuperbeeLimiter return 2
    else
        solver_error(:config, "This limiter is not recognized by armon_cpp")
    end
end


function test_case_to_int(params::ArmonParameters)
    if     params.test isa Armon.Sod       return 0
    elseif params.test isa Armon.Sod_y     return 1
    elseif params.test isa Armon.Sod_circ  return 2
    elseif params.test isa Armon.Bizarrium return 3
    elseif params.test isa Armon.Sedov     return 4
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
    debug_kernels = false, use_md_iter = 0,
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
    option!(armon_cpp, "USE_SIMD_KERNELS", use_md_iter == 0 && params.use_simd; prefix="")
    option!(armon_cpp, "USE_2D_ITER", use_md_iter == 1; prefix="")
    option!(armon_cpp, "USE_MD_ITER", use_md_iter >= 2; prefix="")
    option!(armon_cpp, "BALANCE_MD_ITER", use_md_iter == 3; prefix="")
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


function Armon.print_device_info(io::IO, pad::Int, p::ArmonParameters{<:Any, Dim, <:Kokkos.ExecutionSpace}) where {Dim}
    Armon.print_parameter(io, pad, "use_kokkos", true)

    device_str = Kokkos.main_space_type(p.device) |> nameof |> string
    if p.device isa Kokkos.Cuda || p.device isa Kokkos.HIP
        device_id = Kokkos.BackendFunctions.device_id(p.device)
        device_str *= " (GPU index: $device_id)"
    end
    Armon.print_parameter(io, pad, "device", device_str)

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


function Base.wait(::ArmonParameters{<:Any, Dim, <:Kokkos.ExecutionSpace}) where {Dim}
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

# TODO: fix

@generated function Armon.dtCFL_kernel(
    params::ArmonParameters{T, Dim, <:Kokkos.ExecutionSpace}, data::BlockGrid{V}, range, dx, dy
) where {T, Dim, V <: AbstractArray{T}}
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
    params::ArmonParameters{T, Dim, <:Kokkos.ExecutionSpace}, data::BlockGrid{V}, range
) where {T, Dim, V <: Kokkos.View{T}}
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

function Base.copyto!(dst_blk::LocalTaskBlock{A, Size}, src_blk::LocalTaskBlock{B, Size}) where {A <: Kokkos.View, B <: Kokkos.View, Size}
    if state(dst_blk) != Done || state(src_blk) != Done
        error("Both destination and source blocks must be `Done` to copy")
    elseif block_size(dst_blk) != block_size(src_blk)
        error("Destination and source blocks have different sizes: $(block_size(dst_blk)) != $(block_size(src_blk))")
    end

    for (dst_var, src_var) in zip(main_arrays(dst_blk), main_arrays(src_blk))
        Kokkos.deep_copy(Armon.device_type(dst_blk), dst_var, src_var)
    end
end

end
