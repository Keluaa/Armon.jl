
"""
    MPISyncCommunicationModel

Synchronous communication model.
Sends and receives are done at the same time using `MPI.Sendrecv`, preventing deadlocks.
"""
struct MPISyncCommunicationModel <: AbstractCommunicationModel
    comm::MPI.Comm
end

buffer_type(::Type{MPISyncCommunicationModel}, ::Type{A}) where {A} = A
is_async(::ObjOrType{MPISyncCommunicationModel}) = false

function Base.show(io::IO, ::MPISyncCommunicationModel)
    print(io, "MPISyncCommunicationModel()")
end


struct MPISyncP2P{A} <: AbstractCommunication{A}
    model       :: MPISyncCommunicationModel
    send_buffer :: MPI.Buffer{A}
    recv_buffer :: MPI.Buffer{A}
    rank        :: Int
    side        :: Int
    tag         :: Int
end

MPISyncP2P(m, send::A, recv::A, rank, side, tag) where {A} =
    new{A}(m, MPI.Buffer(send), MPI.Buffer(recv), rank, side, tag)

unsafe_send_buffer(c::MPISyncP2P) = (c.send_buffer.data,)
unsafe_recv_buffer(c::MPISyncP2P) = (c.recv_buffer.data,)

exchange_position(c::MPISyncP2P) = (; rank=c.rank, side=c.side, side_pos=c.tag)

function init_exchange(
    model::MPISyncCommunicationModel,
    rank, side, side_pos, array_type, buffer_size, total_side_buffer_size
)
    send_buf = array_type(undef, buffer_size)
    recv_buf = array_type(undef, buffer_size)
    return MPISyncP2P(model, send_buf, recv_buf, rank, side, side_pos)
end

try_acquire_send_buffer!(c::MPISyncP2P) = c.send_buffer.data
acquire_send_buffer!(c::MPISyncP2P) = c.send_buffer.data

function release_send_buffer!(c::MPISyncP2P)
    MPI.Sendrecv!(c.send_buffer, c.recv_buffer, c.model.comm;
        dest=c.rank, sendtag=c.tag,
        source=c.rank, recvtag=c.tag
    )
end

send_completed(::MPISyncP2P) = true
wait_send_completed(::MPISyncP2P) = true

try_acquire_recv_buffer!(c::MPISyncP2P) = c.recv_buffer.data
acquire_recv_buffer!(c::MPISyncP2P) = c.recv_buffer.data
release_recv_buffer!(::MPISyncP2P) = nothing

recv_completed(::MPISyncP2P) = true
wait_recv_completed(::MPISyncP2P) = true


struct MPISyncCollective{A} <: AbstractCommunication{A}
    model       :: MPISyncCommunicationModel
    send_buffer :: MPI.Buffer{A}
    recv_buffer :: MPI.Buffer{A}
    op          :: MPI.Op
end

MPISyncCollective(model, send::A, recv::A, op) where {A} = new{A}(model, MPI.Buffer(send), MPI.Buffer(recv), MPI.Op(op))

unsafe_send_buffer(c::MPISyncCollective) = (c.send_buffer.data,)
unsafe_recv_buffer(c::MPISyncCollective) = (c.recv_buffer.data,)


function init_reduce_broadcast(model::MPISyncCommunicationModel, reduction_op, array_type, count)
    send_buf = array_type(undef, count)
    recv_buf = array_type(undef, count)
    return MPISyncCollective(model, send_buf, recv_buf, reduction_op)
end

try_acquire_send_buffer!(c::MPISyncCollective) = c.send_buffer.data
acquire_send_buffer!(c::MPISyncCollective) = c.send_buffer.data
release_send_buffer!(c::MPISyncCollective) = MPI.Allreduce!(c.send_buf, c.recv_buf, c.op, c.model.comm)

send_completed(::MPISyncCollective) = true
wait_send_completed(::MPISyncCollective) = true

try_acquire_recv_buffer!(c::MPISyncCollective) = c.recv_buffer.data
acquire_recv_buffer!(c::MPISyncCollective) = c.recv_buffer.data
release_recv_buffer!(::MPISyncCollective) = nothing

recv_completed(::MPISyncCollective) = true
wait_recv_completed(::MPISyncCollective) = true
