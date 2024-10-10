
"""
    RMACommunicationModel

Asynchronous communication model using MPI's Remote Memory Access interface.
   
Point-to-point operations are implemented using active RMA.

Options:
 - `homogenous::Bool = true`: if all ranks use the same data type and layout
 - `origin_is_source::Bool = true`: use `MPI_Put`, `MPI_Get` otherwise
"""
struct RMACommunicationModel
    comm             :: MPI.Comm
    homogenous       :: Bool  # If all ranks use the same data type and layout
    origin_is_source :: Bool  # If the origin rank is the source of the data (`MPI.Put!`, or `MPI.Get!` otherwise)
end

function RMACommunicationModel(comm::MPI.Comm; homogenous::Bool=true, origin_is_source::Bool=true)
    error("NYI")
    RMACommunicationModel(comm, homogenous, origin_is_source)
end

buffer_type(::ObjOrType{RMACommunicationModel}, ::Type{A}) where {A} = A
is_thread_safe(::ObjOrType{RMACommunicationModel}) = false
supports_collectives(::ObjOrType{RMACommunicationModel}) = false  # TODO: implement it with passive aggregation on a root rank, which then broadcasts its result

function Base.show(io::IO, model::RMACommunicationModel)
    print(io, "RMACommunicationModel(; homogenous=", model.homogenous, ", origin_is_source=", model.origin_is_source, ")")
end


struct OneWayCommunication{A} <: AbstractCommunication{A}
    model         :: RMACommunicationModel
    target_rank   :: Int  # either the source or destination based on `model.origin_is_source`
    side          :: Int
    side_pos      :: Int
    xchg_buffer   :: MPI.Buffer{A}
    window_buffer :: MPI.Buffer{A}
    window        :: MPI.Win
    comm_group    :: MPI.Group
end

unsafe_send_buffer(owc::OneWayCommunication) = (owc.xchg_buffer.data,)
unsafe_recv_buffer(owc::OneWayCommunication) = (owc.window_buffer.data,)

exchange_position(c::OneWayCommunication) = (; rank=c.target_rank, side=c.side, side_pos=c.side_pos)


function init_exchange(
    model::RMACommunicationModel,
    rank, side, side_pos, array_type, buffer_size, total_side_buffer_size
)
    xchg_buffer = array_type(undef, buffer_size)
    win_buffer  = array_type(undef, buffer_size)

    # TODO: since remote blocks are created asynchronously, it IS needed to have synchronization AT
    # window creation. This way tags in window operations will be implicit. However tags must be
    # explicitly used at window creation.

    # TODO: instead it might be smarter to have a single window per rank side, but idk how MPI implements
    #   concurrent window accesses at non-overlapping places

    # See the MPI spec v4.1, section 12.2.1, for a detailed explaination on the kwargs
    window = MPI.Win_create(win_buffer, model.comm;
        no_locks=true,  # no passive RMA
        same_size=model.homogenous,  # same window for both parties
        same_disp_unit=model.homogenous  # same data type for both parties
    )

    global_group = MPI.Comm_group(model.comm)
    comm_group = MPI.Group_incl(global_group, Int32[rank])

    return OneWayCommunication{array_type}(
        model, rank, side, side_pos,
        MPI.Buffer(xchg_buffer), MPI.Buffer(win_buffer),
        window, comm_group
    )
end


function start_send_epoch(owc::OneWayCommunication; wait=false)
    if owc.model.origin_is_source
        MPI_Win_start(owc.group, owc.window)
        return true
    else
        if wait
            MPI_Win_wait(owc.window)
        elseif !MPI_Win_test(owc.window)
            return false  # previous epoch is not yet completed
        end
        MPI_Win_post(owc.group, owc.window)
        return true
    end
end


function complete_send_epoch(owc::OneWayCommunication)
    if owc.model.origin_is_source
        # It is the source which initiates the communication
        MPI.Put!(owc.xchg_buffer, owc.window; rank=owc.target_rank)
        MPI_Win_complete(owc.window)
    else
        # It is the destination which initiates the communication
        # Instead of closing the window immediately, we keep it open, and it is subsequent calls to
        # `MPI_Win_test()`/`MPI_Win_wait()` which will close it.
    end
end


function start_recv_epoch(owc::OneWayCommunication; wait=false)
    if owc.model.origin_is_source
        if wait
            MPI_Win_wait(owc.window)
        elseif !MPI_Win_test(owc.window)
            return false  # previous epoch is not yet completed
        end
        MPI_Win_post(owc.group, owc.window)
        return true
    else
        MPI_Win_start(owc.group, owc.window)
        return true
    end
end


function complete_recv_epoch(owc::OneWayCommunication)
    if owc.model.origin_is_source
        # It is the source which initiates the communication
        # Instead of closing the window immediately, we keep it open, and it is subsequent calls to
        # `MPI_Win_test()`/`MPI_Win_wait()` which will close it.
    else
        # It is the destination which initiates the communication
        MPI.Get!(owc.xchg_buffer, owc.window; rank=owc.target_rank)
        MPI_Win_complete(owc.window)
    end
end


# TODO: I messed up which buffer goes where
try_acquire_send_buffer!(owc::OneWayCommunication) = start_send_epoch(owc) ? owc.xchg_buffer.data : nothing
function acquire_send_buffer!(owc::OneWayCommunication)
    start_send_epoch(owc; wait=true)
    return owc.xchg_buffer.data
end

release_send_buffer!(owc::OneWayCommunication) = complete_send_epoch(owc)


try_acquire_recv_buffer!(owc::OneWayCommunication) = start_recv_epoch(owc) ? owc.xchg_buffer.data : nothing
function acquire_recv_buffer!(owc::OneWayCommunication)
    start_recv_epoch(owc; wait=true)
    return owc.xchg_buffer.data
end

release_recv_buffer!(owc::OneWayCommunication) = complete_recv_epoch(owc)


# TODO
function send_completed end
function wait_send_completed end

function recv_completed end
function wait_recv_completed end

