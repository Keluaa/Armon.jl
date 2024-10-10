
function try_acquire_atomic_lock!(lock::Atomic{Int})
    tid = Threads.threadid()
    (@atomic :monotonic lock.x) == tid && return true  # already acquired
    _, locked = @atomicreplace lock.x 0 => tid
    return locked
end

@noinline wait_lock_timeout(timeout) = error("could not acquire the lock after $timeout seconds")

function wait_acquire_atomic_lock!(lock::Atomic{Int}; timeout=10.0, pollint=0.001, timeout_error=true)
    try_acquire_atomic_lock!(lock) && return true
    res = timedwait(timeout; pollint) do
        return try_acquire_atomic_lock!(lock)
    end
    res === :time_out && timeout_error && wait_lock_timeout(timeout)
    return res === :ok
end

release_atomic_lock!(lock::Atomic{Int}) = @atomic lock.x = 0


function atomic_lock!(f, lock::Atomic{Int})
    locked = try_acquire_atomic_lock!(lock)
    !locked && return false
    ret = f()
    release_atomic_lock!(lock)
    return ret
end


function wait_for_atomic_lock!(f, lock::Atomic{Int}; kwargs...)
    wait_acquire_atomic_lock!(lock; kwargs..., timeout_error=true)
    ret = f()
    release_atomic_lock!(lock)
    return ret
end

#
# API for MPI_Iallreduce
#

function MPI_IAllreduce!(rbuf::MPI.RBuffer, op::Union{MPI.Op, MPI.MPI_Op}, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request())
    @assert MPI.isnull(req)
    # int MPI_Allreduce(const void* sendbuf, void* recvbuf, int count,
    #                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
    #                   MPI_Request* req)
    MPI.API.MPI_Iallreduce(rbuf.senddata, rbuf.recvdata, rbuf.count, rbuf.datatype, op, comm, req)
    MPI.setbuffer!(req, rbuf)
    return req
end

MPI_IAllreduce!(rbuf::MPI.RBuffer, op, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request()) =
    MPI_IAllreduce!(rbuf, MPI.Op(op, eltype(rbuf)), comm, req)
MPI_IAllreduce!(sendbuf, recvbuf, op, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request()) =
    MPI_IAllreduce!(MPI.RBuffer(sendbuf, recvbuf), op, comm, req)

# inplace
MPI_IAllreduce!(rbuf, op, comm::MPI.Comm, req::MPI.AbstractRequest=MPI.Request()) =
    MPI_IAllreduce!(MPI.IN_PLACE, rbuf, op, comm, req)


function MPI_Allreduce_init(send_buf, recv_buf, count, datatype, op, comm, info, req)
    MPI.API.@mpichk ccall(
        (:MPIX_Allreduce_init, MPI.API.libmpi),
        Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Cint, MPI.API.MPI_Datatype, MPI.API.MPI_Op, MPI.Comm, MPI.API.MPI_Info, Ptr{MPI.API.MPI_Request}),
        send_buf, recv_buf, count, datatype, op, comm, info, req
    )
    return req
end

#
# API for MPI 4 partitioned communications
#

function MPI_partitioned_p2p_supported()
    try
        cglobal((:MPI_Psend_init, MPI.API.libmpi))
        return true
    catch _
        return false
    end
end


mutable struct PartitionedRequest <: MPI.AbstractRequest
    val::MPI.API.MPI_Request
end

function PartitionedRequest()
    req = PartitionedRequest(MPI.API.MPI_REQUEST_NULL[])
    return finalizer(MPI.free, req)
end

MPI.setbuffer!(req::PartitionedRequest, val) = nothing

Base.cconvert(::Type{MPI.API.MPI_Request}, request::PartitionedRequest) = request
Base.unsafe_convert(::Type{MPI.API.MPI_Request}, request::PartitionedRequest) = request.val
Base.unsafe_convert(::Type{Ptr{MPI.API.MPI_Request}}, request::PartitionedRequest) = convert(Ptr{MPI.API.MPI_Request}, pointer_from_objref(request))


function MPI_Parrived(req::PartitionedRequest, partition::Integer)
    flag = Ref{Cint}(0)
    # int MPI_Parrived(MPI_Request request, int partition, int *flag)
    MPI.API.@mpichk ccall(
        (:MPI_Parrived, MPI.API.libmpi),
        Cint, (MPI.API.MPI_Request, Cint, Ptr{Cint}),
        req, partition - 1, flag
    )

    @static if MPI.MPI_LIBRARY == "OpenMPI" && MPI.MPI_LIBRARY_VERSION < v"5-"
        # The initial implementation of partitioned communications in OpenMPI 4 didn't include a
        # progress call in `MPI_Parrived`, leading to systematic deadlocks in receive loops.
        # See https://github.com/open-mpi/ompi/pull/10077
        # The fix here is the same fix as in `MPI_Parrived` in the above PR.
        if flag[] == 0
            ccall((:opal_progress, MPI.API.libmpi), Cint, ())
        end
    end

    return flag[] != 0
end

# MPI_Parrived can be called multiple times on a partition, unlike with MPI_Test, so this is correct
MPI_Parrived(req::PartitionedRequest, partitions) = all(i -> MPI_Parrived(req, i), partitions)


function MPI_Pready!(req::PartitionedRequest, partition::Integer)
    # int MPI_Pready(int partitions, MPI_Request request)
    MPI.API.@mpichk ccall(
        (:MPI_Pready, MPI.API.libmpi),
        Cint, (Cint, MPI.API.MPI_Request),
        partition - 1, req
    )
    return req
end


function MPI_Pready!(req::PartitionedRequest, partition_range::UnitRange{<:Integer})
    isempty(partition_range) && return req  # handles cases where `first(range) > last(range)`
    # int MPI_Pready_range(int partition_low, int partition_high, MPI_Request request)
    MPI.API.@mpichk ccall(
        (:MPI_Pready_range, MPI.API.libmpi),
        Cint, (Cint, Cint, MPI.API.MPI_Request),
        first(partition_range) - 1, last(partition_range) - 1, req
    )
    return req
end


function MPI_Pready!(req::PartitionedRequest, partition_list::Vector{<:Integer})
    # int MPI_Pready_list(int length, int partition_list[], MPI_Request request)
    MPI.API.@mpichk ccall(
        (:MPI_Pready_list, MPI.API.libmpi),
        Cint, (Cint, Cint, MPI.API.MPI_Request),
        length(partition_list), Cint.(partition_list .- 1), req
    )
    return req
end


function MPI_Precv_init(
    buf::MPI.Buffer, comm::MPI.Comm, req::MPI.AbstractRequest=PartitionedRequest();
    partitions::Integer, count::Integer, source::Integer=MPI.ANY_SOURCE, tag::Integer=MPI.ANY_TAG, infokws...
)
    MPI_Precv_init(buf.data, partitions, count, buf.datatype, source, tag, comm, MPI.Info(infokws...), req)
    MPI.setbuffer!(req, buf)
    return req
end

function MPI_Precv_init(buf, partitions, count, datatype, source, tag, comm, info, req)
    # int MPI_Precv_init(void* buf, int partitions, MPI_Count count,
    #                    MPI_Datatype datatype, int source, int tag, MPI_Comm comm,
    #                    MPI_Info info, MPI_Request *request)
    MPI.API.@mpichk ccall(
        (:MPI_Precv_init, MPI.API.libmpi),
        Cint,
        (Ptr{Cvoid}, Cint, Cint, MPI.API.MPI_Datatype, Cint, Cint, MPI.API.MPI_Comm, MPI.API.MPI_Info, Ptr{MPI.API.MPI_Request}),
        buf, partitions, count, datatype, source, tag, comm, info, req
    )
end


function MPI_Psend_init(
    buf::MPI.Buffer, comm::MPI.Comm, req::MPI.AbstractRequest=PartitionedRequest();
    partitions::Integer, count::Integer, dest::Integer, tag::Integer=0, infokws...
)
    MPI_Psend_init(buf.data, partitions, count, buf.datatype, dest, tag, comm, MPI.Info(infokws...), req)
    MPI.setbuffer!(req, buf)
    return req
end

function MPI_Psend_init(buf, partitions, count, datatype, dest, tag, comm, info, req)
    # int MPI_Psend_init(const void* buf, int partitions, MPI_Count count,
    #                    MPI_Datatype datatype, int dest, int tag, MPI_Comm comm,
    #                    MPI_Info info, MPI_Request *req)
    MPI.API.@mpichk ccall(
        (:MPI_Psend_init, MPI.API.libmpi),
        Cint,
        (Ptr{Cvoid}, Cint, Cint, MPI.API.MPI_Datatype, Cint, Cint, MPI.API.MPI_Comm, MPI.API.MPI_Info, Ptr{MPI.API.MPI_Request}),
        buf, partitions, count, datatype, dest, tag, comm, info, req
    )
end

#
# API for MPI active RMA
#

window_assert_value(assert::Integer) = Cint(assert)

function window_assert_value(assert::Symbol, constraints::NTuple{N, Symbol}) where {N}
    # See MPI spec section 12.5.5 for an explaination of each flag
    if assert ∉ constraints
        error("expected $(join(constraints, ", ", " or ")), got: $assert")
    elseif assert === :none        return Cint(0)
    elseif assert === :no_check    return MPI.API.MPI_MODE_NOCHECK[]
    elseif assert === :no_store    return MPI.API.MPI_MODE_NOSTORE[]
    elseif assert === :no_put      return MPI.API.MPI_MODE_NOPUT[]
    elseif assert === :no_precede  return MPI.API.MPI_MODE_NOPRECEDE[]
    elseif assert === :no_succeed  return MPI.API.MPI_MODE_NOSUCCEED[]
    else                           return Cint(0)
    end
end

function window_assert_value(asserts::Tuple{Vararg{Symbol}}, constraints::NTuple{N, Symbol}) where {N}
    return |(window_assert_value.(asserts, Ref(constraints))...)
end

function window_assert_value(asserts::Vector{Symbol}, constraints::NTuple{N, Symbol}) where {N}
    return |(window_assert_value.(asserts, Ref(constraints))...)
end


function MPI_Win_start(group::MPI.Group, window::MPI.Win; asserts=:none)
    assert_val = window_assert_value(asserts, (:none, :no_check))
    MPI.API.MPI_Win_start(group, assert_val, window)
    return nothing
end


function MPI_Win_post(group::MPI.Group, window::MPI.Win; asserts=:none)
    assert_val = window_assert_value(asserts, (:none, :no_check, :no_store, :no_put))
    MPI.API.MPI_Win_post(group, assert_val, window)
    return nothing
end


function MPI_Win_complete(window::MPI.Win)
    MPI.API.MPI_Win_complete(window)
    return nothing
end


function MPI_Win_test(window::MPI.Win)
    flag = Ref{Cint}(0)
    MPI.API.MPI_Win_test(window, flag)
    return flag[] != 0
end


function MPI_Win_wait(window::MPI.Win)
    MPI.API.MPI_Win_wait(window)
    return nothing
end


function MPI_Win_fence(window::MPI.Win; asserts=:none)
    assert_val = window_assert_value(asserts, (:none, :no_store, :no_put, :no_precede, :no_succeed))
    MPI.API.MPI_Win_fence(assert_val, window)
    return nothing
end

