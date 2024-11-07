
"""
    EmptyCommunicationModel

Communication model where exactly nothing is communicated. For use as a placeholder.
"""
struct EmptyCommunicationModel <: AbstractCommunicationModel end

buffer_type(::ObjOrType{EmptyCommunicationModel}, ::Type{A}) where {A} = A
is_thread_safe(::ObjOrType{EmptyCommunicationModel}) = true

function Base.show(io::IO, ::EmptyCommunicationModel)
    print(io, "EmptyCommunicationModel()")
end


struct EmptyCommunication{A} <: AbstractCommunication{A} end

unsafe_send_buffer(::EmptyCommunication) = ()
unsafe_recv_buffer(::EmptyCommunication) = ()
exchange_position(::EmptyCommunication) = (; rank=1, side=1, side_pos=1)

init_exchange(::EmptyCommunicationModel, _, _, _, array_type, _, _) = EmptyCommunication{array_type}()
init_reduce_broadcast(::EmptyCommunicationModel, _, array_type, _) = EmptyCommunication{array_type}()

try_acquire_send_buffer!(::EmptyCommunication{A}) where {A} = eltype(A)[]
acquire_send_buffer!(::EmptyCommunication{A}) where {A} = eltype(A)[]
release_send_buffer!(::EmptyCommunication) = nothing

send_completed(::EmptyCommunication) = true
wait_send_completed(::EmptyCommunication) = true

try_acquire_recv_buffer!(::EmptyCommunication{A}) where {A} = eltype(A)[]
acquire_recv_buffer!(::EmptyCommunication{A}) where {A} = eltype(A)[]
release_recv_buffer!(::EmptyCommunication) = nothing

recv_completed(::EmptyCommunication) = true
wait_recv_completed(::EmptyCommunication) = true

finalize_comm!(::EmptyCommunication) = nothing
