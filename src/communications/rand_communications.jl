
"""
    RandCommunicationModel(; s=5e4, m=500)

Communication model without communications: the process communicates with itself.

Unlike [`NoCommunicationModel`](@ref), communications choose a random completion time after they are
started. This allows to test the behaviour of communication overlap with compute.

The communication time follows an [exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution):
`x = m + s Z` with `Z` an exponential distribution of scale `1`.
Hence, the mean communication time would be `m + s` (by default `50` µs), and the minimum is `m` (by
default `500` ns). All values are in nanoseconds.
"""
struct RandCommunicationModel <: AbstractCommunicationModel
    s :: Float64
    m :: Float64
end

RandCommunicationModel(; s=50e3, m=500) = RandCommunicationModel(s, m)

buffer_type(::ObjOrType{RandCommunicationModel}, ::Type{A}) where {A} = A
is_thread_safe(::ObjOrType{RandCommunicationModel}) = true

function Base.show(io::IO, model::RandCommunicationModel)
    print(io, "RandCommunicationModel(s=", model.s, ", m=", model.m, ")")
end


function random_completion_time(model::RandCommunicationModel)
    nano_time = model.m + model.s * randexp()
    # Safely convert from Float64 to UInt64
    nano_time > typemax(UInt64) && (nano_time %= typemax(UInt64))
    return time_ns() + trunc(UInt64, nano_time)  # we don't care about overflows
end


mutable struct RandCommunication{A} <: AbstractCommunication{A}
    model         :: RandCommunicationModel
    buffer        :: A
    send_end_time :: UInt64
    recv_end_time :: UInt64
end

unsafe_send_buffer(c::RandCommunication) = (c.buffer,)
unsafe_recv_buffer(c::RandCommunication) = (c.buffer,)
exchange_position(::RandCommunication) = (; rank=1, side=1, side_pos=1)

function init_exchange(
    model::RandCommunicationModel,
    rank, side, side_pos, array_type, buffer_size, total_side_buffer_size
)
    buffer = array_type(undef, buffer_size)
    current_time = time_ns()
    recv_time = random_completion_time(model)
    return RandCommunication{typeof(buffer)}(model, buffer, current_time, recv_time)
end

function init_reduce_broadcast(model::RandCommunicationModel, reduction_op, array_type, count)
    buffer = array_type(undef, count)
    current_time = time_ns()
    recv_time = random_completion_time(model)
    return RandCommunication{typeof(buffer)}(model, buffer, current_time, recv_time)
end


function try_acquire_send_buffer!(c::RandCommunication)
    !send_completed(c) && return nothing
    return c.buffer
end

function acquire_send_buffer!(c::RandCommunication)
    wait_send_completed(c)
    return c.buffer
end

function release_send_buffer!(c::RandCommunication)
    c.send_end_time = random_completion_time(c.model)
    return nothing
end

send_completed(c::RandCommunication) = time_ns() ≥ c.send_end_time
function wait_send_completed(c::RandCommunication)
    while !send_completed(c)
        ccall(:jl_cpu_pause, Cvoid, ())
    end
    return true
end


function try_acquire_recv_buffer!(c::RandCommunication)
    !recv_completed(c) && return nothing
    return c.buffer
end

function acquire_recv_buffer!(c::RandCommunication)
    wait_recv_completed(c)
    return c.buffer
end

function release_recv_buffer!(c::RandCommunication)
    c.recv_end_time = random_completion_time(c.model)
    return nothing
end

recv_completed(c::RandCommunication) = time_ns() ≥ c.recv_end_time
function wait_recv_completed(c::RandCommunication)
    while !recv_completed(c)
        ccall(:jl_cpu_pause, Cvoid, ())
    end
    return true
end


finalize_comm!(::RandCommunication) = nothing
