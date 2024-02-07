
"""
    ArmonData{V}

Generic array holder for all variables and temporary variables used throughout the solver.
"""
struct ArmonData{V <: AbstractArray}
    # TODO: rename some variables (rho -> ρ, umat -> u, work_array_1 -> work_1, ustar -> uˢ, etc...)
    x            :: V
    y            :: V
    rho          :: V
    umat         :: V
    vmat         :: V
    Emat         :: V
    pmat         :: V
    cmat         :: V
    gmat         :: V
    ustar        :: V
    pstar        :: V
    work_array_1 :: V
    work_array_2 :: V
    work_array_3 :: V
    work_array_4 :: V
    domain_mask  :: V
end


ArmonData(params::ArmonParameters{T}) where T = ArmonData(T, params.nbcell, params.comm_array_size)
ArmonData(type::Type, size::Int64, comm_size::Int64) = ArmonData(Vector{type}, size, comm_size)

function ArmonData(array::Type{V}, size::Int64, comm_size::Int64; kwargs...) where {V <: AbstractArray}
    a1 = array(undef, size; alloc_array_kwargs(; label="x", kwargs...)...)

    # For very small domains (~20×20), the communication array might be bigger than a normal data array
    comm_size = max(size, comm_size)

    complete_array_type = typeof(a1)
    return ArmonData{complete_array_type}(
        a1,
        array(undef, size; alloc_array_kwargs(; label="y", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="rho", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="umat", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="vmat", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="Emat", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="pmat", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="cmat", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="gmat", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="ustar", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="pstar", kwargs...)...),
        array(undef, comm_size; alloc_array_kwargs(; label="work_array_1", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="work_array_2", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="work_array_3", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="work_array_4", kwargs...)...),
        array(undef, size; alloc_array_kwargs(; label="domain_mask", kwargs...)...)
    )
end


array_type(::ArmonData{V}) where V = V


main_variables() = (:x, :y, :rho, :umat, :vmat, :Emat, :pmat, :cmat, :gmat, :ustar, :pstar, :domain_mask)
main_variables(data::ArmonData; more_vars=()) = map(f -> getfield(data, f), Iterators.flatten((main_variables(), more_vars)))

saved_variables() = (:x, :y, :rho, :umat, :vmat, :pmat)
saved_variables(data::ArmonData; more_vars=()) = map(f -> getfield(data, f), Iterators.flatten((saved_variables(), more_vars)))


"""
    memory_required(params::ArmonParameters)

Compute the number of bytes needed on the device to allocate all data arrays.

While the result is precise, it does not account for additional memory required by MPI buffers and
the solver, as well as resources shared among MPI processes (i.e. CPU memory).
"""
memory_required(params::ArmonParameters{T}) where T = memory_required(params.nbcell, T)

function memory_required(N, float_type)
    field_count = fieldcount(ArmonData{Vector{float_type}})
    floats = field_count * N
    return floats * sizeof(float_type)
end


"""
    ArmonDualData{DeviceArray, HostArray, BufferArray, Device}

Holds two version of `ArmonData`, one for the `Device` and one for the host, as well as the buffers
necessary for the halo exchange.

If the host and device are the same, the `device_data` and `host_data` fields point to the same data.

`device` might be a `KernelAbstractions.Device`, `AMDGPU.ROCDevice` or `Kokkos.ExecutionSpace`.
"""
struct ArmonDualData{DeviceArray <: AbstractArray, HostArray <: AbstractArray, BufferArray <: AbstractArray, Device}
    device       :: Device
    device_data  :: ArmonData{DeviceArray}
    host_data    :: ArmonData{HostArray}
    comm_buffers :: Dict{Side, NamedTuple{(:send, :recv), NTuple{2, MPI.Buffer{BufferArray}}}}
    requests     :: Dict{Side, NamedTuple{(:send, :recv), NTuple{2, MPI.AbstractRequest}}}
end


function ArmonDualData(params::ArmonParameters{T}) where T
    device_array = device_array_type(params.device){T, 1}
    host_array = host_array_type(params.device){T, 1}

    device_data = ArmonData(device_array, params.nbcell, params.comm_array_size; alloc_device_kwargs(params)...)
    if host_array == device_array
        host_data = device_data
    else
        host_data = ArmonData(host_array, params.nbcell, params.comm_array_size; alloc_host_kwargs(params)...)
    end

    # By changing the types, we change which methods of `copy_to_send_buffer` and `copy_from_recv_buffer`
    # gets called => GPU-awareness with very low code footprint
    if params.gpu_aware
        buf_array = device_array
        buf_type = array_type(device_data)
    else
        buf_array = host_array
        buf_type = array_type(host_data)
    end

    # In case we don't use MPI: since there is no neighbours, no array is allocated.
    comm_buffers = Dict{Side, NamedTuple{(:send, :recv), NTuple{2, MPI.Buffer{buf_type}}}}()
    requests = Dict{Side, NamedTuple{(:send, :recv), NTuple{2, MPI.AbstractRequest}}}()
    for side in instances(Side)
        has_neighbour(params, side) || continue
        neighbour = neighbour_at(params, side)
        comm_buffers[side] = (
            send = MPI.Buffer(buf_array(undef, params.comm_array_size)),
            recv = MPI.Buffer(buf_array(undef, params.comm_array_size))
        )
        requests[side] = (
            send = MPI.Send_init(comm_buffers[side].send, params.cart_comm; dest=neighbour),
            recv = MPI.Recv_init(comm_buffers[side].recv, params.cart_comm; source=neighbour)
        )
    end

    return ArmonDualData{array_type(device_data), array_type(host_data), buf_type, typeof(params.device)}(
        params.device, device_data, host_data, comm_buffers, requests
    )
end


device_type(data::ArmonDualData) = data.device
device(data::ArmonDualData) = data.device_data
host(data::ArmonDualData) = data.host_data

"""
    iter_send_requests(data::ArmonDualData)

Iterator over all active MPI send requests
"""
iter_send_requests(data::ArmonDualData) = 
    Iterators.map(p -> first(p) => first(last(p)), 
        Iterators.filter(!MPI.isnull ∘ first ∘ last, data.requests))

"""
    iter_recv_requests(data::ArmonDualData)

Iterator over all active MPI receive requests
"""
iter_recv_requests(data::ArmonDualData) = 
    Iterators.map(p -> first(p) => last(last(p)), 
        Iterators.filter(!MPI.isnull ∘ last ∘ last, data.requests))

"""
    buffers_on_device(data::ArmonDualData)

`true` if the MPI buffers are on the device, i.e. no copy from device to host is needed.
This is the case when working on the CPU only or when using GPU-aware MPI.
"""
buffers_on_device(::ArmonDualData{D, H, B}) where {D, H, B} = D == B

send_buffer(data::ArmonDualData, side::Side) = data.comm_buffers[side].send
recv_buffer(data::ArmonDualData, side::Side) = data.comm_buffers[side].recv

# MPI buffers are on the device (or everything is on the host) => no copy needed
get_send_comm_array(data::ArmonDualData{D, H, D}, side::Side) where {D, H} = send_buffer(data, side).data
get_recv_comm_array(data::ArmonDualData{D, H, D}, side::Side) where {D, H} = recv_buffer(data, side).data

function get_send_comm_array(data::ArmonDualData{D, H, B}, side::Side) where {D, H, B}
    D == B && return  # The buffer is on the device
    # Generic case: device to host copies may be needed
    if side == Left
        comm_array = device(data).work_array_1
    elseif side == Right
        comm_array = device(data).work_array_2
    elseif side == Bottom
        comm_array = device(data).work_array_3
    else
        comm_array = device(data).work_array_4
    end
    array_range = Base.OneTo(length(send_buffer(data, side).data))
    return view(comm_array, array_range)
end

get_recv_comm_array(data::ArmonDualData{D, H, B}, side::Side) where {D, H, B} = get_send_comm_array(data, side)


"""
    device_to_host!(data::ArmonDualData)

Copies all `main_variables()` from the device to the host data. A no-op if the host is the device.
"""
device_to_host!(::ArmonDualData{D, D}) where D = nothing


"""
    host_to_device!(data::ArmonDualData)

Copies all `main_variables()` from the host to the device data. A no-op if the host is the device.
"""
host_to_device!(::ArmonDualData{D, D}) where D = nothing


function device_to_host!(data::ArmonDualData{D, H}) where {D, H}
    for var in main_variables()
        copyto!(getfield(host(data), var), getfield(device(data), var))
    end
end


function host_to_device!(data::ArmonDualData{D, H}) where {D, H}
    for var in main_variables()
        copyto!(getfield(device(data), var), getfield(host(data), var))
    end
end


# The device data is already in the send/recv buffer (see `get_send_comm_array` and `get_recv_comm_array`)
# For CPU, `D == H`
copy_to_send_buffer!(::ArmonDualData{D, H, D}, ::D, ::Side)   where {D, H} = nothing
copy_to_send_buffer!(::ArmonDualData{D, H, D}, ::D, ::D)      where {D <: AbstractArray, H} = nothing
copy_from_recv_buffer!(::ArmonDualData{D, H, D}, ::D, ::Side) where {D, H} = nothing
copy_from_recv_buffer!(::ArmonDualData{D, H, D}, ::D, ::D)    where {D <: AbstractArray, H} = nothing


function copy_to_send_buffer!(data::ArmonDualData{D, H, B}, array::D, side::Side) where {D, H, B}
    buffer_data = send_buffer(data, side).data
    copy_to_send_buffer!(data, array, buffer_data)
end


function copy_to_send_buffer!(data::ArmonDualData{D, H, B}, array::D, buffer::B) where {D, H, B <: AbstractArray}
    D == B && return  # Buffers are already on the device
    array_data = view(array, 1:length(buffer))
    KernelAbstractions.copyto!(device_type(data), buffer, array_data)
end


function copy_from_recv_buffer!(data::ArmonDualData{D, H, B}, array::D, side::Side) where {D, H, B}
    buffer_data = recv_buffer(data, side).data
    copy_from_recv_buffer!(data, array, buffer_data)
end


function copy_from_recv_buffer!(data::ArmonDualData{D, H, B}, array::D, buffer::B) where {D, H, B <: AbstractArray}
    D == B && return  # Buffers are already on the device
    array_data = view(array, 1:length(buffer))
    KernelAbstractions.copyto!(device_type(data), array_data, buffer)
end


ArmonDataOrDual = Union{ArmonData, ArmonDualData}
