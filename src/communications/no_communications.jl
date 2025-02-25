
"""
    NoCommunicationModel

Communication model without communications: the process communicates with itself.
"""
struct NoCommunicationModel <: AbstractCommunicationModel end

buffer_type(::ObjOrType{NoCommunicationModel}, ::Type{A}) where {A} = A
is_thread_safe(::ObjOrType{NoCommunicationModel}) = true

function Base.show(io::IO, ::NoCommunicationModel)
    print(io, "NoCommunicationModel()")
end


struct DummyCommunication{A} <: AbstractCommunication{A}
    model::NoCommunicationModel
    # The send and receive buffers are the same, therefore the communication is equivalent to sending
    # an email to yourself.
    buffer::A
end

unsafe_send_buffer(c::DummyCommunication) = (c.buffer,)
unsafe_recv_buffer(c::DummyCommunication) = (c.buffer,)
exchange_position(::DummyCommunication) = (; rank=1, side=1, side_pos=1)

function init_exchange(
    model::NoCommunicationModel,
    rank, side, side_pos, array_type, buffer_size, total_side_buffer_size
)
    buffer = array_type(undef, buffer_size)
    return DummyCommunication{typeof(buffer)}(model, buffer)
end

function init_reduce_broadcast(model::NoCommunicationModel, reduction_op, array_type, count)
    buffer = array_type(undef, count)
    return DummyCommunication{typeof(buffer)}(model, buffer)
end


try_acquire_send_buffer!(c::DummyCommunication) = c.buffer
acquire_send_buffer!(c::DummyCommunication) = c.buffer
release_send_buffer!(::DummyCommunication) = nothing

send_completed(::DummyCommunication) = true
wait_send_completed(::DummyCommunication) = true

try_acquire_recv_buffer!(c::DummyCommunication) = c.buffer
acquire_recv_buffer!(c::DummyCommunication) = c.buffer
release_recv_buffer!(::DummyCommunication) = nothing

recv_completed(::DummyCommunication) = true
wait_recv_completed(::DummyCommunication) = true

finalize_comm!(::DummyCommunication) = nothing
