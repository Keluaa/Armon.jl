module Communications

using Random
using MPI
using ..Armon: Atomic, ObjOrType

export AbstractCommunicationModel, AbstractCommunication


"""
    AbstractCommunicationModel

Represents a communication model.
"""
abstract type AbstractCommunicationModel end


"""
    buffer_type(model::AbstractCommunicationModel, array_type)

The buffer type used by the `model` for the given `array_type`.
"""
function buffer_type end


"""
    is_thread_safe(model::AbstractCommunicationModel)
    is_thread_safe(model::Type{AbstractCommunicationModel})

`true` if all operations (excluding `unsafe` ones) on [`AbstractCommunication`](@ref) built with
`model` can be called concurrently by different threads.

!!! warning

    [`try_acquire_send_buffer!`](@ref)/[`try_acquire_recv_buffer!`](@ref) should be followed by a
    call to [`release_send_buffer!`](@ref)/[`release_recv_buffer!`](@ref) from the **same** thread,
    respectively (it also applies to [`acquire_send_buffer!`](@ref)/[`acquire_recv_buffer!`](@ref)).

!!! note

    Collective operations can be launched/waited upon from any thread, but only one thread should
    start the communication. It the user's responsibility to ensure this.
"""
is_thread_safe(::ObjOrType{AbstractCommunicationModel}) = false


"""
    is_async(::AbstractCommunicationModel)
    is_async(::Type{AbstractCommunicationModel})

`true` if the `AbstractCommunicationModel` is asynchronous.
"""
is_async(::ObjOrType{AbstractCommunicationModel}) = true


"""
    supports_point_to_point(::AbstractCommunicationModel)
    supports_point_to_point(::Type{AbstractCommunicationModel})

`true` if the `AbstractCommunicationModel` implements [`init_exchange`](@ref).
"""
supports_point_to_point(::ObjOrType{AbstractCommunicationModel}) = true


"""
    supports_collectives(::AbstractCommunicationModel)
    supports_collectives(::Type{AbstractCommunicationModel})

`true` if the `AbstractCommunicationModel` implements [`init_reduce_broadcast`](@ref).
"""
supports_collectives(::ObjOrType{AbstractCommunicationModel}) = true


"""
    uses_global_buffers(::AbstractCommunicationModel)
    uses_global_buffers(::Type{AbstractCommunicationModel})

`true` if the `AbstractCommunicationModel` uses global communication buffers, therefore calling
[`unsafe_send_buffer`](@ref), [`unsafe_recv_buffer`](@ref) or [`unsafe_buffers`](@ref) on different
`AbstractCommunication` instances could return the same values.
Those three methods are the only ones affected by global buffer usage, others should return views on
local buffers.
"""
uses_global_buffers(::ObjOrType{AbstractCommunicationModel}) = false


"""
    AbstractCommunication{A}

Base type representing any communication between processes, using array buffers of type `A`.
"""
abstract type AbstractCommunication{A} end


is_thread_safe(c::AbstractCommunication) = is_thread_safe(c.model)
is_async(c::AbstractCommunication) = is_async(c.model)
uses_global_buffers(c::AbstractCommunication) = uses_global_buffers(c.model)


"""
    init_exchange(
        model::AbstractCommunicationModel,
        rank::Int, side::Int, side_pos::Int,
        array_type::Type, buffer_size::Int, total_side_buffer_size::Int
    )

Initialize an exchange between a local block and a remote block in another process.

`side_pos` must be a unique number among all communications between the current process and `rank`
along the `side` of the communication.

The exchange would use buffers of `array_type` and `buffer_size`.
The sum of all `buffer_size`s of the exchanges along `side` is `total_side_buffer_size`.

Returns a [`AbstractCommunication`](@ref), whose exact concrete type depends on `model`.

!!! info

    The order of operations in an exchange is important, as it is expected that sends are always
    followed by a receive:
     - first [`try_acquire_send_buffer!`](@ref) (or [`acquire_send_buffer!`](@ref)), then [`release_send_buffer!`](@ref)
     - second [`try_acquire_recv_buffer!`](@ref) (or [`acquire_recv_buffer!`](@ref)), then [`release_recv_buffer!`](@ref)
"""
function init_exchange end


"""
    init_reduce_broadcast(model::AbstractCommunicationModel, reduction_op, array_type, count)

Initialize a reduction operation whose result is broadcasted to all processes.

The `reduction_op` is applied on `count` values of `type`.

Returns a [`AbstractCommunication`](@ref), whose exact concrete type depends on `model`.
"""
function init_reduce_broadcast end


"""
    exchange_position(c::AbstractCommunication)

Return `(; rank, side, side_pos)` as given to the constructor of `c` in [`init_exchange`](@ref).
"""
function exchange_position end


"""
    unsafe_send_buffer(comm::AbstractCommunication{B})

A `Tuple{Vararg{B}}` containing all send buffers (usually the only one) of `comm`.

This is threads/MPI unsafe: no mechanisms ensure that the buffers are not in use.
The send buffers might be the same as the receive buffers.
"""
function unsafe_send_buffer end


"""
    unsafe_recv_buffer(comm::AbstractCommunication{B})

A `Tuple{Vararg{B}}` containing all receive buffers (usually the only one) of `comm`.

This is threads/MPI unsafe: no mechanisms ensure that the buffers are not in use.
The receive buffers might be the same as the send buffers.
"""
function unsafe_recv_buffer end


"""
    unsafe_buffers(comm::AbstractCommunication{B})

An iterator over all send and receive buffers of `comm`.
Some buffers can be iterated over more than once, e.g. where the send and receive buffers are the same.

See [unsafe_send_buffer](@ref) and [unsafe_recv_buffer](@ref).
"""
unsafe_buffers(c::AbstractCommunication) = Iterators.flatten((unsafe_send_buffer(c), unsafe_recv_buffer(c)))


"""
    try_acquire_send_buffer!(comm::AbstractCommunication)

Return the `send_buffer` if it can be written to by the current process, implying that no communication
using it is in progress. Otherwise, `nothing` is returned.

`send_buffer` may change from one call to another.

It is required to call [`release_send_buffer!`](@ref) as soon as the buffer is ready to be sent.

```julia
try_acquire_send_buffer!(comm)
# fill the buffer...
release_send_buffer!(comm)  # perform the send operation
```
"""
function try_acquire_send_buffer! end


"""
    acquire_send_buffer!(comm::AbstractCommunication)

Same as [`try_acquire_send_buffer!`](@ref), but blocks until the buffer is available.
Can be expensive for thread-safe and/or asynchronous communications.
"""
function acquire_send_buffer! end


"""
    release_send_buffer!(comm::AbstractCommunication)

Mark the `send_buffer` returned by [`try_acquire_send_buffer!`](@ref) as unwritable.
The underlying mechanism of `comm` is then free to send the buffer to the remote process whenever it
wants to.
"""
function release_send_buffer! end


"""
    try_acquire_recv_buffer!(comm::AbstractCommunication)

Return the `recv_buffer` if it can be written to by the current process, implying that no communication
using it is in progress. Otherwise, `nothing` is returned.

`recv_buffer` may change from one call to another.

It is required to call [`release_recv_buffer!`](@ref) as soon as the buffer is ready to be receive
new data.

```julia
try_acquire_recv_buffer!(comm)
# retreive data from the buffer...
release_recv_buffer!(comm)  # we are now ready to receive more data
```
"""
function try_acquire_recv_buffer! end


"""
    acquire_recv_buffer!(comm::AbstractCommunication)

Same as [`try_acquire_recv_buffer!`](@ref), but blocks until the buffer is available.
Can be expensive for thread-safe and/or asynchronous communications.
"""
function acquire_recv_buffer! end


"""
    release_recv_buffer!(comm::AbstractCommunication)

Mark the `recv_buffer` returned by [`try_acquire_recv_buffer!`](@ref) as unreadable.
The underlying mechanism of `comm` is then free to write to the buffer.
"""
function release_recv_buffer! end


"""
    send_completed(comm::AbstractCommunication)

`true` if the send side of `comm` is completed.

For thread-safe [`AbstractCommunicationModel`](@ref)s, if another thread acquired a global lock
needed to check if `comm` is completed, `false` is returned, regardless of the actual underlying
state of `comm`.
"""
function send_completed end


"""
    wait_send_completed(comm::AbstractCommunication)

Block until the send side of `comm` is completed.

For thread-safe [`AbstractCommunicationModel`](@ref)s, if another thread is waiting for `comm`,
`false` is returned (`true` otherwise).
"""
function wait_send_completed end


"""
    recv_completed(comm::AbstractCommunication)

`true` if the receive side of `comm` is completed.

For thread-safe [`AbstractCommunicationModel`](@ref)s, if another thread acquired a global lock
needed to check if `comm` is completed, `false` is returned, regardless of the actual underlying
state of `comm`.

Some models might use permanently active receive requests, which are started at initialization.
In this case, even if no communication has been started, `false` will be returned, unlike with
[`send_completed`](@ref).
"""
function recv_completed end


"""
    wait_recv_completed(comm::AbstractCommunication)

Block until the receive side of `comm` is completed.

For thread-safe [`AbstractCommunicationModel`](@ref)s, if another thread is waiting for `comm`,
`false` is returned (`true` otherwise).
"""
function wait_recv_completed end


"""
    finalize_comm!(comm::AbstractCommunication)

Finalize `comm` and free global ressources associated with it.
`comm` is unusable afterward, much like a call to `finalize`.

This call is thread-safe and supposes that no other thread is using any buffers, or starting another
communication.

For some models (such as [`MPIPartitionedCommunicationModel`](@ref)), this is a necessary step to
ensure no requests are kept active, possibly affecting future MPI calls by mismatching requests.
"""
function finalize_comm! end


include("mpi_extra.jl")
include("empty_communication.jl")
include("no_communications.jl")
include("rand_communications.jl")
include("mpi_sync_communications.jl")
include("mpi_async_communications.jl")
include("mpi_async_thread_safe.jl")
include("mpi_rma.jl")
include("mpi_partitioned.jl")
include("threads_collectives.jl")


"""
    communication_model(name::Symbol, comm::MPI.Comm; kwargs...)

Builds a new [`AbstractCommunicationModel`](@ref) from its name and an existing `MPI.Comm`unicator.

`kwargs` are passed to the model's constructor as-is.

Possible values for `name`:
 - `:async_safe`   [`MPIAsyncSafeCommunicationModel`](@ref)
 - `:async`        [`MPIAsyncCommunicationModel`](@ref)
 - `:sync`         [`MPISyncCommunicationModel`](@ref)
 - `:rma`          [`RMACommunicationModel`](@ref)
 - `:partitioned`  [`MPIPartitionedCommunicationModel`](@ref)
 - `:no_comms`     [`NoCommunicationModel`](@ref)
 - `:rand_comms`   [`RandCommunicationModel`](@ref)
"""
function communication_model(name::Symbol, comm::MPI.Comm; kwargs...)
    if name === :async_safe
        return MPIAsyncSafeCommunicationModel(comm; kwargs...)
    elseif name === :async
        return MPIAsyncCommunicationModel(comm; kwargs...)
    elseif name === :sync
        return MPISyncCommunicationModel(comm; kwargs...)
    elseif name === :rma
        return RMACommunicationModel(comm; kwargs...)
    elseif name === :partitioned
        if !MPI_partitioned_p2p_supported()
            error("MPI partitioned communications are not supported by the current implementation")
        end
        return MPIPartitionedCommunicationModel(comm; kwargs...)
    elseif name === :no_comms
        return NoCommunicationModel()
    elseif name === :rand_comms
        return RandCommunicationModel(; kwargs...)
    else
        error("Unknown communication model type: $name")
    end
end

end
