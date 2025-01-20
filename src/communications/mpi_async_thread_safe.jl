
"""
    MPIAsyncSafeCommunicationModel

Asynchronous communication model. All communications are done asynchronously, using thread-safety
mechanisms which prevent several threads from operating on a communication (MPI request and/or buffer)
at the same time.

See [`MPIAsyncCommunicationModel`](@ref) for a thread-unsafe alternative, which doesn't use any
global locks.

Options:
 - `persistant_reduction::Bool = false`: use persistant requests for reduction operations (`MPI_Allreduce_init`)
   It defaults to `false` as persistant collective operations are not correctly supported by most
   MPI implementations.
"""
mutable struct MPIAsyncSafeCommunicationModel <: AbstractCommunicationModel
    comm::MPI.Comm
    persistant_reduction::Bool
end

MPIAsyncSafeCommunicationModel(comm::MPI.Comm; persistant_reduction::Bool=false) =
    MPIAsyncSafeCommunicationModel(comm, persistant_reduction)

buffer_type(::ObjOrType{MPIAsyncSafeCommunicationModel}, ::Type{A}) where {A} = A
is_thread_safe(::ObjOrType{MPIAsyncSafeCommunicationModel}) = true

function Base.show(io::IO, model::MPIAsyncSafeCommunicationModel)
    print(io, "MPIAsyncSafeCommunicationModel(; persistant_reduction=", model.persistant_reduction, ")")
end


mutable struct MPIAsyncSafeP2P{A} <: AbstractCommunication{A}
    model        :: MPIAsyncSafeCommunicationModel
    send_buffer  :: MPI.Buffer{A}
    recv_buffer  :: MPI.Buffer{A}
    requests     :: NTuple{2, @NamedTuple{
        send     :: MPI.Request,
        recv     :: MPI.Request
    }}
    # The MPI communication ordering is guareenteed using those indexes to `requests`
    send_state   :: Int
    recv_state   :: Int  # Note: to make sends match with recvs, `recv_state` starts as 2 and `send_state` starts at 1
    # The buffers and requests atomicity is guareenteed with those atomics, by storing the thread ID
    # using them. `0` means that they are currently not owned.
    send_lock    :: Atomic{Int}
    recv_lock    :: Atomic{Int}
    rank         :: Int
    side         :: Int
    side_pos     :: Int
end

function MPIAsyncSafeP2P(model, send::A, recv::A, rank, side, side_pos) where {A}
    c = MPIAsyncSafeP2P{A}(model, MPI.Buffer(send), MPI.Buffer(recv), (
        (; send=MPI.Request(), recv=MPI.Request()),
        (; send=MPI.Request(), recv=MPI.Request()),
    ), 1, 2, Atomic{Int}(0), Atomic{Int}(0), rank, side, side_pos)

    finalizer(c) do c_obj
        # Cancel the active receive request, as they are always active otherwise
        if !MPI.Finalized() && !MPI.Test(c_obj.requests[c_obj.recv_state].recv)
            MPI.Cancel!(c_obj.requests[c_obj.recv_state].recv)
        end
    end

    return c
end

unsafe_send_buffer(c::MPIAsyncSafeP2P) = (c.send_buffer.data,)
unsafe_recv_buffer(c::MPIAsyncSafeP2P) = (c.recv_buffer.data,)

exchange_position(c::MPIAsyncSafeP2P) = (; rank=c.rank, side=c.side, side_pos=c.side_pos)


function init_exchange(
    model::MPIAsyncSafeCommunicationModel,
    rank, side, side_pos, array_type, buffer_size, total_side_buffer_size
)
    send_buf = array_type(undef, buffer_size)
    recv_buf = array_type(undef, buffer_size)
    p2p_data = MPIAsyncSafeP2P(model, send_buf, recv_buf, rank, side, side_pos)

    tag_1 = (side_pos << 1)
    tag_2 = (side_pos << 1) | 1

    # Same buffers, different tags. By alternating which request is started/waited upon, we can enforce
    # an ordering of the permanent communications without using barriers between communication cycles.
    # Without this, if a remote process were to send two successive communications for the same request,
    # the order in which they are receive could be undetermined.
    MPI.Send_init(p2p_data.send_buffer, model.comm, p2p_data.requests[1].send; dest=rank,   tag=tag_1)
    MPI.Recv_init(p2p_data.recv_buffer, model.comm, p2p_data.requests[1].recv; source=rank, tag=tag_1)

    MPI.Send_init(p2p_data.send_buffer, model.comm, p2p_data.requests[2].send; dest=rank,   tag=tag_2)
    MPI.Recv_init(p2p_data.recv_buffer, model.comm, p2p_data.requests[2].recv; source=rank, tag=tag_2)

    # One receive request is always active when this rank isn't touching the receive buffer.
    MPI.Start(p2p_data.requests[p2p_data.recv_state].recv)

    return p2p_data
end


function try_acquire_send_buffer!(c::MPIAsyncSafeP2P)
    !try_acquire_atomic_lock!(c.send_lock) && return nothing
    return MPI.Test(c.requests[c.send_state].send) ? c.send_buffer.data : nothing
end

function acquire_send_buffer!(c::MPIAsyncSafeP2P)
    wait_acquire_atomic_lock!(c.send_lock)  # will throw after 10 seconds of waiting
    MPI.Wait(c.requests[c.send_state].send)
    return c.send_buffer.data
end

function release_send_buffer!(c::MPIAsyncSafeP2P)
    c.send_state = mod1(c.send_state + 1, 2)
    MPI.Start(c.requests[c.send_state].send)
    release_atomic_lock!(c.send_lock)
end

function send_completed(c::MPIAsyncSafeP2P)
    return atomic_lock!(c.send_lock) do
        MPI.Test(c.requests[c.send_state].send)
    end
end

function wait_send_completed(c::MPIAsyncSafeP2P)
    return atomic_lock!(c.send_lock) do
        MPI.Wait(c.requests[c.send_state].send)
        return true
    end
end


function try_acquire_recv_buffer!(c::MPIAsyncSafeP2P)
    !try_acquire_atomic_lock!(c.recv_lock) && return nothing
    return MPI.Test(c.requests[c.recv_state].recv) ? c.recv_buffer.data : nothing
end

function acquire_recv_buffer!(c::MPIAsyncSafeP2P)
    wait_acquire_atomic_lock!(c.recv_lock)  # will throw after 10 seconds of waiting
    MPI.Wait(c.requests[c.recv_state].recv)
    return c.recv_buffer.data
end

function release_recv_buffer!(c::MPIAsyncSafeP2P)
    c.recv_state = mod1(c.recv_state + 1, 2)
    MPI.Start(c.requests[c.recv_state].recv)
    release_atomic_lock!(c.recv_lock)
end

function recv_completed(c::MPIAsyncSafeP2P)
    return atomic_lock!(c.recv_lock) do
        MPI.Test(c.requests[c.recv_state].recv)
    end
end

function wait_recv_completed(c::MPIAsyncSafeP2P)
    return atomic_lock!(c.recv_lock) do
        MPI.Wait(c.requests[c.recv_state].recv)
        return true
    end
end


struct MPIAsyncSafeCollective{A} <: AbstractCommunication{A}
    model         :: MPIAsyncSafeCommunicationModel
    send_buffer   :: MPI.Buffer{A}
    recv_buffer   :: MPI.Buffer{A}
    request       :: MPI.Request
    request_lock  :: Atomic{Int}
    send_lock     :: Atomic{Int}
    recv_lock     :: Atomic{Int}
    op            :: MPI.Op
end

function MPIAsyncSafeCollective(model, send::A, recv::A, op) where {A}
    mpi_op = op isa MPI.Op ? op : MPI.Op(op, eltype(send))
    return MPIAsyncSafeCollective{A}(
        model,
        MPI.Buffer(send), MPI.Buffer(recv),
        MPI.Request(), Atomic{Int}(0),
        Atomic{Int}(0), Atomic{Int}(0),
        mpi_op
    )
end

unsafe_send_buffer(c::MPIAsyncSafeCollective) = (c.send_buffer.data,)
unsafe_recv_buffer(c::MPIAsyncSafeCollective) = (c.recv_buffer.data,)


function init_reduce_broadcast(model::MPIAsyncSafeCommunicationModel, reduction_op, array_type, count)
    send_buf = array_type(undef, count)
    recv_buf = array_type(undef, count)
    comm = MPIAsyncSafeCollective(model, send_buf, recv_buf, reduction_op)
    if model.persistant_reduction
        MPI_Allreduce_init(
            send_buf, recv_buf,
            count, MPI.Datatype(eltype(array_type)), comm.op,
            model.comm, MPI.Info(), comm.request
        )
    end
    return comm
end


function try_acquire_send_buffer!(c::MPIAsyncSafeCollective)
    !try_acquire_atomic_lock!(c.send_lock) && return nothing
    ok = atomic_lock!(c.request_lock) do
        MPI.Test(c.request)
    end
    if ok
        return c.send_buffer.data
    else
        release_atomic_lock!(c.send_lock)
        return nothing
    end
end

function acquire_send_buffer!(c::MPIAsyncSafeCollective)
    wait_acquire_atomic_lock!(c.send_lock)
    wait_for_atomic_lock!(c.request_lock) do
        MPI.Wait(c.request)
    end
    return c.send_buffer.data
end

function release_send_buffer!(c::MPIAsyncSafeCollective)
    if c.model.persistant_reduction
        MPI.Start(c.request)
    else
        MPI_IAllreduce!(c.send_buffer.data, c.recv_buffer.data, c.op, c.model.comm, c.request)
    end
    release_atomic_lock!(c.send_lock)
end

function send_completed(c::MPIAsyncSafeCollective)
    return atomic_lock!(c.request_lock) do
        MPI.Test(c.request)
    end
end

function wait_send_completed(c::MPIAsyncSafeCollective)
    return atomic_lock!(c.request_lock) do
        MPI.Wait(c.request)
        return true
    end
end


function try_acquire_recv_buffer!(c::MPIAsyncSafeCollective)
    !try_acquire_atomic_lock!(c.recv_lock) && return nothing
    ok = atomic_lock!(c.request_lock) do
        MPI.Test(c.request)
    end
    if ok
        return c.recv_buffer.data
    else
        release_atomic_lock!(c.recv_lock)
        return nothing
    end
end

function acquire_recv_buffer!(c::MPIAsyncSafeCollective)
    wait_acquire_atomic_lock!(c.recv_lock)
    wait_for_atomic_lock!(c.request_lock) do
        MPI.Wait(c.request)
    end
    return c.recv_buffer.data
end

release_recv_buffer!(c::MPIAsyncSafeCollective) = release_atomic_lock!(c.recv_lock)

recv_completed(c::MPIAsyncSafeCollective) = send_completed(c)
wait_recv_completed(c::MPIAsyncSafeCollective) = wait_send_completed(c)
