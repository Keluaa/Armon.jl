
"""
    ThreadsCommunicationModel{InterProcessModel}

Communication model for which an inter-process communication is triggered once a certain amount of
threads have contributed.

The communication begins once all threads have contributed to the send **and** receive buffers, that
is all of them have completed all read/write operations on the buffers.
This ensures that a communication epoch has elapsed on all threads before starting a new one.
"""
mutable struct ThreadsCommunicationModel{InterProcessModel <: AbstractCommunicationModel} <: AbstractCommunicationModel
    inter_process_model :: InterProcessModel
    contributors        :: Int
    threads_mask        :: BitVector  # 1: non-contributing thread, 0: contributing thread
end

function ThreadsCommunicationModel(ipm_model, threads)
    threads_mask = trues(Threads.nthreads())  # By default a thread will not contribute
    if threads isa Integer
        contributors = Int(threads)
        threads_mask[1:contributors] .= 0  # Only the first `threads` will be contributing
    else
        contributors = length(threads)
        threads_mask[threads] .= 0
    end
    return ThreadsCommunicationModel(ipm_model, contributors, threads_mask)
end


buffer_type(::ObjOrType{ThreadsCommunicationModel{IPM}}, ::Type{A}) where {IPM, A} = buffer_type(IPM, A)
is_thread_safe(::ObjOrType{ThreadsCommunicationModel}) = true
supports_point_to_point(::ObjOrType{ThreadsCommunicationModel}) = false

function Base.show(io::IO, model::ThreadsCommunicationModel{IPM}) where {IPM}
    print(io, "ThreadsCommunicationModel{$IPM}(; contributors=$(model.contributors))")
end


struct ThreadCollective{A, InterProcessModel, IPM_Collective} <: AbstractCommunication{A}
    model              :: ThreadsCommunicationModel{InterProcessModel}
    IPM_collective     :: IPM_Collective
    send_buffer        :: A
    recv_buffer        :: A
    contributions      :: BitVector
    # Once both `send_contributions` and `recv_contributions` are all ones, the `IPM_Collective` is started
    send_contributions :: BitVector
    recv_contributions :: BitVector
    all_contributed    :: Atomic{UInt8}
    send_lock          :: Atomic{Bool}
    recv_lock          :: Atomic{Bool}
end

function ThreadCollective(model::ThreadsCommunicationModel{IPM}, IPM_collective, send::A, recv::A) where {IPM, A}
    send_contributions = copy(model.threads_mask)
    recv_contributions = copy(model.threads_mask)
    return ThreadCollective{A, typeof(F), IPM, typeof(IPM_collective)}(
        model, IPM_collective,
        send, recv,
        send_contributions, recv_contributions, Atomic{UInt8}(0b00),
        Atomic{Bool}(false), Atomic{Bool}(false)
    )
end

unsafe_send_buffer(c::ThreadCollective) = (c.send_buffer,)
unsafe_recv_buffer(c::ThreadCollective) = (c.recv_buffer,)


function init_reduce_broadcast(model::ThreadsCommunicationModel, reduction_op, array_type, count)
    IPM_collective = init_reduce_broadcast(model.inter_process_model, reduction_op, array_type, count)
    send_buf = only(unsafe_send_buffer(IPM_collective))
    recv_buf = only(unsafe_recv_buffer(IPM_collective))
    return ThreadCollective(model, IPM_collective, send_buf, recv_buf)
end


function try_start_new_epoch(c::ThreadCollective, side::Symbol)
    sides_state = if side === :send
        @atomic c.all_contributed.x |= 0b01
    else#if side === :recv
        @atomic c.all_contributed.x |= 0b10
    end
    sides_state != 0b11 && return  # The other side is not done yet

    # Avoid edge-cases relating to the collective finishing immediately by reseting this variable
    # first.
    @atomic c.all_contributed.x = 0b00

    # All threads have contributed: trigger the inter-process collective
    acquire_send_buffer!(c.IPM_collective)  # May call `MPI.Wait` (or similar)
    release_send_buffer!(c.IPM_collective)
end


function wait_start_new_epoch(c::ThreadCollective, side::Symbol; timeout=10.0, pollint=0.001, timeout_error=true)
    # Wait until the current epoch is finished.
    side_lock = side === :send ? c.send_lock : c.recv_lock

    res = timedwait(timeout; pollint) do
        if try_acquire_atomic_lock!(side_lock)
            if all(c.send_contributions) && (side === :send ? send_completed(c.IPM_collective) : recv_completed(c.IPM_collective))
                release_atomic_lock!(side_lock)
                return false
            else
                return true  # note that we don't release the lock
            end
        end
    end

    if res === :time_out
        release_atomic_lock!(side_lock)
        timeout_error && wait_lock_timeout(timeout)
        return false
    end

    side_contribs = side === :send ? c.send_contributions : c.recv_contributions
    tid = Threads.threadid()
    if side_contribs[tid]
        # Then contributions have not been reset, it is up to this thread to do it
        copy!(c.send_contributions, c.model.threads_mask)
    end

    return true  # the `side_lock` is still acquired at this point
end


function try_acquire_send_buffer!(c::ThreadCollective)
    !try_acquire_atomic_lock!(c.send_lock) && return nothing
    tid = Threads.threadid()
    if c.send_contributions[tid]
        if all(c.send_contributions) && send_completed(c.IPM_collective)
            # This is the beginning of a new epoch
            copy!(c.send_contributions, c.model.threads_mask)  # reset the contributions
        else
            # The thread already contributed for this communication epoch, or the underlying
            # collective has not finished sending the previous send buffer.
            release_atomic_lock!(c.send_lock)
            return nothing
        end
    end
    return c.send_buffer
end


function acquire_send_buffer!(c::ThreadCollective)
    wait_send_completed(c)
    wait_acquire_atomic_lock!(c.send_lock)
    tid = Threads.threadid()
    if c.send_contributions[tid]
        if all(c.send_contributions) && send_completed(c.IPM_collective)
            # This is the beginning of a new epoch
            copy!(c.send_contributions, c.model.threads_mask)  # reset the contributions
        else
            # This thread is trying to contribute to a new epoch, while the current one isn't completed.
            release_atomic_lock!(c.send_lock)  # Other threads will need it to advance the current epoch's state
            wait_start_new_epoch(c, :send)  # When this returns, the send lock is acquired again
        end
    end
    return c.send_buffer
end


function release_send_buffer!(c::ThreadCollective)
    tid = Threads.threadid()
    c.send_contributions[tid] = true
    if all(c.send_contributions)
        try_start_new_epoch(c, :send)
    end
    release_atomic_lock!(c.send_lock)
end


function try_acquire_recv_buffer!(c::ThreadCollective)
    !try_acquire_atomic_lock!(c.recv_lock) && return nothing
    tid = Threads.threadid()
    if c.recv_contributions[tid]
        if all(c.recv_contributions) && recv_completed(c.IPM_collective)
            # This is the beginning of a new epoch
            copy!(c.recv_contributions, c.model.threads_mask)  # reset the contributions
        else
            # The thread already contributed for this communication epoch, or the underlying
            # collective has not finished receiving the previous recv buffer.
            release_atomic_lock!(c.recv_lock)
            return nothing
        end
    end
    return c.recv_buffer
end

function acquire_recv_buffer!(c::ThreadCollective)
    wait_recv_completed(c)
    wait_acquire_atomic_lock!(c.recv_lock)
    tid = Threads.threadid()
    if c.recv_contributions[tid]
        if all(c.recv_contributions) && recv_completed(c.IPM_collective)
            # This is the beginning of a new epoch
            copy!(c.recv_contributions, c.model.threads_mask)  # reset the contributions
        else
            # This thread is trying to contribute to a new epoch, while the current one isn't completed.
            release_atomic_lock!(c.recv_lock)  # Other threads will need it to advance the current epoch's state
            wait_start_new_epoch(c, :recv)  # When this returns, the recv lock is acquired again
        end
    end
    return c.recv_buffer
end

function release_recv_buffer!(c::ThreadCollective)
    tid = Threads.threadid()
    c.recv_contributions[tid] = true
    if all(c.recv_contributions)
        try_start_new_epoch(c, :recv)
    end
    release_atomic_lock!(c.recv_lock)
end


function send_completed(c::ThreadCollective)
    # TODO
end

function wait_send_completed(c::ThreadCollective)
    # TODO
end


function recv_completed(c::ThreadCollective)
    # TODO
end

function wait_recv_completed(c::ThreadCollective)
    # TODO
end
