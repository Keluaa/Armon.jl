
"""
    MPIAsyncCommunicationModel

Asynchronous communication model. All communications are done asynchronously, without thread-safety
mechanisms.

Options:
 - `persistant_reduction::Bool = true`: use persistant requests for reduction operations (`MPI_Allreduce_init`)
   It defaults to `false` as persistant collective operations are not correctly supported by most
   MPI implementations.

!!! warn

    While buffers are obviously not thread-safe, be aware that `MPI.Test` and `MPI.Wait` are also not
    thread-safe when used on the same request. Starting a request from a thread and testing/waiting on
    it in another will also invoke borderline unpredictable behaviour from the MPI implementation.
    See [`MPIAsyncSafeCommunicationModel`](@ref) for a thread-safe alternative.
"""
mutable struct MPIAsyncCommunicationModel <: AbstractCommunicationModel
    comm::MPI.Comm
    persistant_reduction::Bool
end

MPIAsyncCommunicationModel(comm::MPI.Comm; persistant_reduction::Bool=false) =
    MPIAsyncCommunicationModel(comm, persistant_reduction)

buffer_type(::ObjOrType{MPIAsyncCommunicationModel}, ::Type{A}) where {A} = A

function Base.show(io::IO, model::MPIAsyncCommunicationModel)
    print(io, "MPIAsyncCommunicationModel(; persistant_reduction=", model.persistant_reduction, ")")
end


mutable struct MPIAsyncP2P{A} <: AbstractCommunication{A}
    model        :: MPIAsyncCommunicationModel
    send_buffer  :: MPI.Buffer{A}
    recv_buffer  :: MPI.Buffer{A}
    send_request :: MPI.Request
    recv_request :: MPI.Request
    rank         :: Int
    side         :: Int
    side_pos     :: Int
end

function MPIAsyncP2P(model, send::A, recv::A, rank, side, side_pos) where {A}
    c = MPIAsyncP2P{A}(
        model,
        MPI.Buffer(send), MPI.Buffer(recv),
        MPI.Request(), MPI.Request(),
        rank, side, side_pos
    )

    finalizer(c) do c_obj
        # Cancel the active receive request, as it is always active otherwise
        if !MPI.Finalized() && !MPI.Test(c_obj.recv_request)
            MPI.Cancel!(c_obj.recv_request)
        end
    end

    return c
end

unsafe_send_buffer(c::MPIAsyncP2P) = (c.send_buffer.data,)
unsafe_recv_buffer(c::MPIAsyncP2P) = (c.recv_buffer.data,)

exchange_position(c::MPIAsyncP2P) = (; rank=c.rank, side=c.side, side_pos=c.side_pos)

function init_exchange(
    model::MPIAsyncCommunicationModel,
    rank, side, side_pos, array_type, buffer_size, total_side_buffer_size
)
    send_buf = array_type(undef, buffer_size)
    recv_buf = array_type(undef, buffer_size)
    p2p_data = MPIAsyncP2P(model, send_buf, recv_buf, rank, side, side_pos)

    MPI.Send_init(p2p_data.send_buffer, model.comm, p2p_data.send_request; dest=rank, tag=side_pos)
    MPI.Recv_init(p2p_data.recv_buffer, model.comm, p2p_data.recv_request; source=rank, tag=side_pos)

    MPI.Start(p2p_data.recv_request)

    return p2p_data
end


try_acquire_send_buffer!(c::MPIAsyncP2P) = MPI.Test(c.send_request) ? c.send_buffer.data : nothing

function acquire_send_buffer!(c::MPIAsyncP2P)
    buf = try_acquire_send_buffer!(c)
    !isnothing(buf) && return buf
    MPI.Wait(c.send_request)
    return c.send_buffer.data
end

release_send_buffer!(c::MPIAsyncP2P) = MPI.Start(c.send_request)

send_completed(c::MPIAsyncP2P) = MPI.Test(c.send_request)
wait_send_completed(c::MPIAsyncP2P) = (MPI.Wait(c.send_request); true)


try_acquire_recv_buffer!(c::MPIAsyncP2P) = MPI.Test(c.recv_request) ? c.recv_buffer.data : nothing

function acquire_recv_buffer!(c::MPIAsyncP2P)
    buf = try_acquire_recv_buffer!(c)
    !isnothing(buf) && return buf
    MPI.Wait(c.recv_request)
    return c.recv_buffer.data
end

release_recv_buffer!(c::MPIAsyncP2P) = MPI.Start(c.recv_request)

recv_completed(c::MPIAsyncP2P) = MPI.Test(c.recv_request)
wait_recv_completed(c::MPIAsyncP2P) = (MPI.Wait(c.recv_request); true)


struct MPIAsyncCollective{A} <: AbstractCommunication{A}
    model         :: MPIAsyncCommunicationModel
    send_buffer   :: MPI.Buffer{A}
    recv_buffer   :: MPI.Buffer{A}
    request       :: MPI.Request
    op            :: MPI.Op
end

function MPIAsyncCollective(model, send::A, recv::A, op) where {A}
    mpi_op = op isa MPI.Op ? op : MPI.Op(op, eltype(send))
    return MPIAsyncCollective{A}(
        model,
        MPI.Buffer(send), MPI.Buffer(recv),
        MPI.Request(), mpi_op
    )
end

unsafe_send_buffer(c::MPIAsyncCollective) = (c.send_buffer.data,)
unsafe_recv_buffer(c::MPIAsyncCollective) = (c.recv_buffer.data,)


function init_reduce_broadcast(model::MPIAsyncCommunicationModel, reduction_op, array_type, count)
    send_buf = array_type(undef, count)
    recv_buf = array_type(undef, count)
    comm = MPIAsyncCollective(model, send_buf, recv_buf, reduction_op)
    if model.persistant_reduction
        MPI_Allreduce_init(
            comm.send_buf, comm.recv_buf,
            count, MPI.Datatype(eltype(array_type)), comm.op,
            model.comm, MPI.Info(), comm.req
        )
    end
    return comm
end

send_completed(c::MPIAsyncCollective) = MPI.Test(c.request)
wait_send_completed(c::MPIAsyncCollective) = (MPI.Wait(c.request); true)

try_acquire_send_buffer!(c::MPIAsyncCollective) = send_completed(c) ? c.send_buffer.data : nothing
acquire_send_buffer!(c::MPIAsyncCollective) = (wait_send_completed(c); c.send_buffer.data)

function release_send_buffer!(c::MPIAsyncCollective)
    if c.model.persistant_reduction
        MPI.Start(c.request)
    else
        MPI_IAllreduce!(c.send_buffer.data, c.recv_buffer.data, c.op, c.model.comm, c.request)
    end
end

recv_completed(c::MPIAsyncCollective) = send_completed(c)
wait_recv_completed(c::MPIAsyncCollective) = wait_send_completed(c)

try_acquire_recv_buffer!(c::MPIAsyncCollective) = recv_completed(c) ? c.recv_buffer.data : nothing
acquire_recv_buffer!(c::MPIAsyncCollective) = (wait_recv_completed(c); c.recv_buffer.data)
release_recv_buffer!(::MPIAsyncCollective) = nothing
