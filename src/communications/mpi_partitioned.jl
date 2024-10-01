
mutable struct MPIPartitionedP2PGlobalInfo{A}
    send_buffer      :: MPI.Buffer{A}  # Global send buffer
    recv_buffer      :: MPI.Buffer{A}  # Global recv buffer
    send_request     :: PartitionedRequest
    recv_request     :: PartitionedRequest
    send_lock        :: Atomic{Int}   # Restricts the access to `send_request` to a single thread
    send_started     :: Atomic{Bool}  # If `MPI_Start` was called on `send_request`, and it did not complete yet
    recv_count       :: Atomic{Int}   # Number of partitions which acknowledged the arrival of their partition.
    send_cycle       :: Atomic{Int}   # Global cycle of the send request. Incremented every time `MPI_Start` is called.
    recv_cycle       :: Atomic{Int}   # Global cycle of the recv request. Incremented every time `MPI_Start` is called.
    total_partitions :: Int
    rank             :: Int
    side             :: Int
    tag              :: Int
end


"""
    MPIPartitionedCommunicationModel

Thread-safe communication using MPI partitioned communications, available in MPI 4.1.

Options:
 - `partition_size::Int`: number of elements in each partition (mandatory)

!!! note

    Point-to-point exchanges are matched by the combinaison of the `rank` and `side`, therefore there
    should be only one exchange per side.

!!! note

    Point-to-point communication supports buffer sizes (the total data to send) which are not a multiple
    of the `partition_size` (the data to send per partition). However it is expected that the extra
    data (`buffer_size % partition_size`) would be stored in the last partition, therefore the local
    `buffer_size` to give should be different from the others.

!!! warn

    Point-to-point communication call `MPI_Start` when a thread need it, therefore any thread might
    start and complete (with `MPI_Test` or `MPI_Wait`) the partitioned request. This is the tricky
    part, made thread-safe with global locks: only one thread can touch the request at once. Even
    though it is allowed by the MPI standard, it might be incompatible with your implementation.
"""
struct MPIPartitionedCommunicationModel
    comm           :: MPI.Comm
    partition_size :: Int
    # Dict of all initialized partitioned communications with the other rank+side
    partitions     :: Dict{Tuple{Int, Int}, MPIPartitionedP2PGlobalInfo}
    lock           :: ReentrantLock  # lock for dict accesses
end

MPIPartitionedCommunicationModel(comm::MPI.Comm; partition_size) =
    MPIPartitionedCommunicationModel(comm, partition_size, Dict(), ReentrantLock())

buffer_type(::ObjOrType{MPIPartitionedCommunicationModel}, ::Type{A}) where {T, N, A <: AbstractArray{T, N}} =
    SubArray{T, N, A, Tuple{UnitRange{Int}}, true}
is_thread_safe(::ObjOrType{MPIPartitionedCommunicationModel}) = true
supports_collectives(::ObjOrType{MPIPartitionedCommunicationModel}) = false
uses_global_buffers(::ObjOrType{AbstractCommunicationModel}) = true

function Base.show(io::IO, model::MPIPartitionedCommunicationModel)
    print(io, "MPIPartitionedCommunicationModel(; partition_size=", model.partition_size, ")")
end


struct MPIPartitionedP2P{S} <: AbstractCommunication{S}
    model              :: MPIPartitionedCommunicationModel
    global_info        :: MPIPartitionedP2PGlobalInfo
    side_pos           :: Int
    partitions         :: UnitRange{Int}  # Range of partitions associated with this object
    local_send_buffer  :: S  # local view to the global send buffer data
    local_recv_buffer  :: S  # local view to the global recv buffer data
    # About global and local cycles (`model.send_cycle`/`model.recv_cycle` vs. `c.send_cycle`/`c.recv_cycle`):
    #  - one global cycle is completed when MPI_Test returns `true`, i.e. the communication has been
    #    completed and we can start contributing to the partitions once again
    #  - one local cycle is completed when `MPI_Pready` is called on the partition for a send operation,
    #    or when `MPI_Parrived` returns `true` for a receive operation.
    #  - when a cycle is completed, it is incremented. Therefore, in order for a partition to be
    #    accessible, the local cycle must be equal to the global cycle.
    # This way, any thread can start, test and wait on the partitioned request, and it is possible to
    # limit the contension on global locks.
    send_cycle         :: Int  # Local cycle of the send side of the partitions. Incremented every time `MPI_Pready` is called.
    recv_cycle         :: Int  # Local cycle of the recv side of the partitions. Incremented every time `MPI_Parrived` returns `true` for all `partitions`
end


# Return the global buffers instead of the local ones, as it may be more useful, instead of having
# to parse through all partitions
unsafe_send_buffer(c::MPIPartitionedP2P) = (c.global_info.send_buffer.data,)
unsafe_recv_buffer(c::MPIPartitionedP2P) = (c.global_info.recv_buffer.data,)

exchange_position(c::MPIPartitionedP2P) = (; rank=c.global_info.rank, side=c.global_info.rank, side_pos=c.side_pos)

function build_partitioned_communication(model::MPIPartitionedCommunicationModel, rank, side, array_type, required_buffer_size)
    (; partition_size) = model

    total_partitions = ceil(Int, required_buffer_size / partition_size)
    real_buffer_size = total_partitions * partition_size

    send_buffer = MPI.Buffer(array_type(undef, real_buffer_size))
    recv_buffer = MPI.Buffer(array_type(undef, real_buffer_size))

    # The tag for each model/communicator is the number of partitioned communications already
    # initialized with the other rank.
    tag = count((r, _) -> r == rank, keys(model.partitions))

    send_req = MPI_Psend_init(send_buffer, model.comm; partitions=total_partitions, count=partition_size, dest=rank, tag)
    recv_req = MPI_Precv_init(recv_buffer, model.comm; partitions=total_partitions, count=partition_size, source=rank, tag)

    # The receive request is always active, until we need to access the receive buffer
    MPI.Start(recv_req)

    global_c_info = MPIPartitionedP2PGlobalInfo{array_type}(
        send_buffer, recv_buffer,
        send_req, recv_req,
        Atomic{Int}(0),
        Atomic{Bool}(false),
        Atomic{Int}(0),
        Atomic{Int}(0), Atomic{Int}(0),
        total_partitions, rank, side, tag
    )

    finalizer(global_c_info) do c_obj
        # Cancel any active receive request
        if !MPI.Test(c_obj.recv_request)
            MPI.Cancel!(c_obj.recv_request)
        end
    end

    return global_c_info
end


function init_exchange(
    model::MPIPartitionedCommunicationModel,
    rank, side, side_pos, array_type, buffer_size, total_side_buffer_size
)
    # Get (or create) the partition communication for this combinaison of rank+side
    partition_info = lock(model.lock) do
        partition_key = (rank, side)
        return get!(model.partitions, partition_key) do
            return build_partitioned_communication(model, rank, side, array_type, total_side_buffer_size)
        end
    end

    partition_pos = side_pos
    if model.partition_size != buffer_size
        # Split this buffer into multiple partitions
        num_partitions = cld(buffer_size, model.partition_size)
    else
        num_partitions = 1
    end

    partition_range = (1:num_partitions) .+ (partition_pos - 1)
    if num_partitions > 1 && last(partition_range) != partition_info.total_partitions
        # In order for this to be correct, it must be the last buffer, otherwise our communication
        # buffer would be discontinuous if the buffer size isn't a multiple of the partition size,
        # and the next partitions would also need to be shifted.
        error("partitions with a size ≠ `partition_size` must be at the end of the global buffer")
    end

    global_buffer_offset = (first(partition_range) - 1) * model.partition_size
    local_send_buffer = view(partition_info.send_buffer, (1:buffer_size) .+ global_buffer_offset)
    local_recv_buffer = view(partition_info.recv_buffer, (1:buffer_size) .+ global_buffer_offset)

    expected_buffer_type = buffer_type(model, array_type)
    @assert expected_buffer_type == typeof(local_send_buffer)

    return MPIPartitionedP2P{expected_buffer_type}(
        model, partition_info, side_pos, partition_range,
        local_send_buffer, local_recv_buffer,
        0, 0
    )
end


function try_acquire_send_buffer!(c::MPIPartitionedP2P)
    # We avoid the global lock if we can. Since the current cycle cannot complete without the this
    # partition from being ready, this is thread-safe.
    (@atomic c.global_info.send_cycle) == c.send_cycle && return c.local_send_buffer
    # Otherwise we must explicitly check the request
    return atomic_lock!(c.global_info.send_lock) do  # non-blocking lock
        # Another thread could have done the job while we were acquiring the lock, so we must check again
        (@atomic c.global_info.send_cycle) == c.send_cycle && return c.local_send_buffer
        completed = MPI.Test(c.global_info.send_request)
        if completed
            @atomic c.global_info.send_cycle += 1
            @atomic c.global_info.send_started = false
        end
        return completed ? c.local_send_buffer : nothing
    end
end


function acquire_send_buffer!(c::MPIPartitionedP2P)
    # Same as `try_acquire_send_buffer!`
    (@atomic c.global_info.send_cycle) == c.send_cycle && return c.local_send_buffer
    wait_for_atomic_lock!(c.global_info.send_lock) do
        (@atomic c.global_info.send_cycle) == c.send_cycle && return c.local_send_buffer
        MPI.Wait(c.global_info.send_request)
        @atomic c.global_info.send_cycle += 1
        @atomic c.global_info.send_started = false
    end
    return c.local_send_buffer
end


function release_send_buffer!(c::MPIPartitionedP2P)
    # The request needs to be started before we can call `MPI_Pready`
    if !(@atomic c.global_info.send_started)
        # We impose a wait here for convenience.
        # Note that as per the MPI spec, we cannot simplify things by immediately starting the
        # request after it completed, as a requirement for proper call to `MPI_Finalize`, ALL
        # requests must be inactive. This could also cause implementation (MPI-side) problems if
        # several calls to the solver would be made in the same Julia session, as some requests
        # would remain active forever.
        wait_for_atomic_lock!(c.global_info.send_lock) do 
            (@atomic c.global_info.send_started) && return  # thread-safety, etc...
            MPI.Start(c.global_info.send_request)
            @atomic c.global_info.send_started = true
        end
    end

    MPI_Pready!(c.global_info.send_request, c.partitions)
    c.send_cycle += 1  # Now this partition must wait for the request to complete before contributing again
    return nothing
end


wait_send_completed(c::MPIPartitionedP2P) = send_completed(c, true)
function send_completed(c::MPIPartitionedP2P, wait=false)
    if !(@atomic c.global_info.send_started)
        # `MPI.Test` will return `true` in all cases as the request is inactive, but:
        #   - if `c.global_info.send_cycle == c.send_cycle` => all partitions contributed and the
        #     MPI send has already completed
        #   - if `c.global_info.send_cycle != c.send_cycle` => this partition must have started the
        #     request before calling `MPI_Pready`: this isn't a possible scenario
        # Therefore it makes sense to always return `true`.
        return true
    end

    return atomic_lock!(c.global_info.send_lock) do  # non-blocking lock
        # Another thread could have done the job while we were acquiring the lock, so we must check again
        previously_not_completed = (@atomic c.global_info.send_cycle) != c.send_cycle
        if !previously_not_completed
            # If another thread completed the request and started the next one, we mustn't test it as
            # the result wouldn't be relevent for this partition.
            return true
        end

        if wait
            MPI.Wait(c.global_info.send_request)
            completed = true
        else
            completed = MPI.Test(c.global_info.send_request)
        end

        if completed
            # Begin a new communication cycle
            @atomic c.global_info.send_cycle += 1
            @atomic c.global_info.send_started = false
        end

        return completed
    end
end


function try_acquire_recv_buffer!(c::MPIPartitionedP2P)
    # Unlike with the send request, the receive request is started after initialization.
    # Since it is always active, there is no need to `MPI_Test` it.
    return recv_completed(c) ? c.local_recv_buffer : nothing
end


function acquire_recv_buffer!(c::MPIPartitionedP2P)
    # We could use `MPI_Wait` on the whole partitioned request, but this would enforce the use of a
    # global lock in `release_recv_buffer!` to prevent concurrent access to the request, plus one
    # global lock here to make sure only a single thread waits on the request.
    # This is too heavy, and stupid given that `MPI_Parrived` is thread-safe (unlike most MPI
    # functions using requests): instead we call `MPI_Parrived` repeatedly until it returns true.
    wait_recv_completed(c)
    return c.local_recv_buffer
end


function release_recv_buffer!(c::MPIPartitionedP2P)
    # Mark the associated partitions as done.
    current_value = (@atomic c.global_info.recv_count += length(c.partitions))

    if current_value == c.global_info.total_partitions
        # Once all partitions are done, we are safe to start the receive operation again.
        # The atomic `+=` acts as a CAS, so we are thread-safe: no need for a global lock.
        # Since all partitions arrived, the `MPI_Wait` is here only to mark the completed active
        # request as inactive. `MPI_Test` would have done the same, but `MPI_Wait` is safer. No
        # actual wait should take place here.
        MPI.Wait(c.global_info.recv_request)
        MPI.Start(c.global_info.recv_request)
        @atomic c.global_info.recv_count = 0
        @atomic c.global_info.recv_cycle += 1
    end

    c.recv_cycle += 1
    return nothing
end


function recv_completed(c::MPIPartitionedP2P)
    # If the previous exchange is not yet complete then `MPI_Parrived` is irrelevant
    (@atomic c.global_info.recv_cycle) != c.recv_cycle && return false
    return MPI_Parrived(c.global_info.recv_request, c.partitions)
end


function wait_recv_completed(c::MPIPartitionedP2P)
    # Wait 10s max, poll every 1ms. Very suboptimal, but waiting by design is suboptimal.
    # There is a fast path in `timedwait` if the condition is already `true`.
    res = timedwait(10.0; pollint=0.001) do
        return recv_completed(c)
    end
    res === :time_out && timeout_error && wait_lock_timeout(10.0)  # unhelpful error
    return true
end
