
"""
    NoQueue{Device} <: AbstractStepQueue{Device}

An always-empty queue: scheduled steps are executed immediately.
"""
struct NoQueue{Device} <: AbstractStepQueue{Device}
    device :: Device
end

Base.length(::NoQueue) = 0
Base.isempty(::NoQueue) = true
is_full(::NoQueue) = false
is_done(::NoQueue) = true
is_running(q::NoQueue) = false  # all steps are synchronous
Base.empty!(q::NoQueue) = q

function planify_step!(queue::NoQueue, params, state, blk, step)
    # Without a queue, steps are executed immediately
    return true, perform_step(queue.device, step, params, state, blk)
end

process_queue!(::NoQueue, params::ArmonParameters, state::SolverState, blk::LocalTaskBlock) = nothing
update_queue_status!(::NoQueue) = true


"""
    StepQueue(device, expected_max_num_steps)

A queue for [`SolverStep`](@ref)s: calling [`planify_step!`](@ref) would enqueue steps, and
[`process_queue!`](@ref) would execute steps in the same order.

The `device` is where steps are executed. It can be a `CPU` or `GPU`.

`expected_max_num_steps` is maximum capacity of the queue.
"""
mutable struct StepQueue{Device, DeviceQueue} <: AbstractStepQueue{Device}
    device       :: Device
    steps        :: Vector{SolverStep.T}
    length       :: Int
    pos          :: Int
    can_continue :: Bool
    device_queue :: DeviceQueue  # A `DeviceStepQueue` if `device` is a `GPU`
end

function StepQueue(device, expected_max_num_steps)
    if device isa Union{CPU, CPU_HP}
        # On CPU, `process_queue!` will process each step sequentially and synchronously.
        device_queue = NoQueue(device)
    else
        device_queue = DeviceStepQueue(device, expected_max_num_steps)
    end
    queue = StepQueue(device, Vector{SolverStep.T}(undef, expected_max_num_steps), 0, 1, true, device_queue)
    !(device isa Union{CPU, CPU_HP}) && copyto!(device_queue, queue)
    return queue
end

Base.length(q::StepQueue) = q.length
Base.isempty(q::StepQueue) = q.length == 0
is_full(q::StepQueue) = q.length == length(q.steps)
is_done(q::StepQueue) = q.pos > q.length
is_running(q::StepQueue) = is_running(q.device_queue)


function Base.empty!(q::StepQueue)
    q.length = 0
    q.pos = 1
    q.can_continue = true
end


@inbounds function Base.push!(queue::StepQueue, step::SolverStep.T)
    !can_run_on_device(queue.device, step) && return false
    is_full(queue) && return false
    queue.length += 1
    queue.steps[queue.length] = step
    return true
end


"""
    planify_step!(queue, params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, step::SolverStep.T)

Add `step` to the `queue`.

Return `(step_queued, step_executed)`:
- `step_queued` is `true` if the step was added to the `queue`
- `step_executed` is `true` if the step was added to the `queue` and executed immediately (always the
case for [`NoQueue`](@ref))
"""
function planify_step!(queue::StepQueue, params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, step::SolverStep.T)
    if can_run_on_device(queue.device, step)
        return push!(queue, step), false
    else
        # Execute the step without planning
        !is_done(queue) && return false, false  # all previous steps must be done beforehand
        return planify_step!(NoQueue(queue.device), params, state, blk, step)
    end
end


function process_queue!(queue::StepQueue, params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    is_done(queue) && return  # no step to process

    if !(queue.device_queue isa NoQueue)
        if queue.pos == 1
            # Initial step: send all steps to the device
            # TODO: we will repeatedly send the steps if the first one blocks for any reason
            copyto!(queue.device_queue, queue)
        end
        process_queue!(queue.device_queue, params, state, blk)
        return
    end

    pos = queue.pos
    can_continue = true
    for outer pos in queue.pos:queue.length
        can_continue = perform_step(queue.device, @inbounds(queue.steps[pos]), params, state, blk)
        !can_continue && break
    end

    # If `can_continue == true`, then all steps in the queue have been completed, so make `pos > length`
    # in order to mark the queue as done.
    # If `can_continue == false`, then the last step couldn't be completed: the next call to `process_queue!`
    # should retry that same step.
    queue.pos = pos + can_continue
    queue.can_continue = can_continue
    return
end


function update_queue_status!(q::StepQueue)
    q.device_queue isa NoQueue && return true
    is_running(q.device_queue) && return false
    copyto!(q, q.device_queue)  # Retreive the status of the device to the host
    return true
end


"""
    DeviceStepQueue(device, expected_max_num_steps)

Same as [`StepQueue`](@ref), but usable from a GPU kernel.

Steps cannot be scheduled using [`planify_step!`](@ref): this must be done with a [`StepQueue`](@ref)
first, then `copyto!` can be used to send the steps to the device.
"""
struct DeviceStepQueue{Device, StepArray, StatusArray, Event} <: AbstractStepQueue{Device}
    device :: Device
    steps  :: StepArray
    # The status array is required in order to return multiple values from the device
    #  1: number of steps (may be lower than `length(steps)`)
    #  2: index in `steps`
    #  3: `≠ 1` if we cannot continue (i.e. must wait before applying the next steps)
    status :: StatusArray
    event  :: Event
end

function DeviceStepQueue(device, expected_max_num_steps)
    steps  = device_array_type(device){SolverStep.T}(undef, expected_max_num_steps)
    status = device_array_type(device){Int}(undef, 3)
    event = create_kernel_event(device)
    return DeviceStepQueue(device, steps, status, event)
end

Adapt.adapt_structure(to, q::DeviceStepQueue) =
    DeviceStepQueue(q.device, Adapt.adapt(to, q.steps), Adapt.adapt(to, q.status), nothing)

# Those 4 methods can only be called from the device
Base.length(q::DeviceStepQueue) = @inbounds(q.status[1])
Base.isempty(q::DeviceStepQueue) = length(q) == 0
is_full(q::DeviceStepQueue) = length(q) == length(q.steps)
is_done(q::DeviceStepQueue) = @inbounds(q.status[2]) > length(q)

function is_running(q::DeviceStepQueue)
    isnothing(q.event) && return true
    # TODO: this is somewhat expensive (>1µs), make sure to not abuse it
    return !query_kernel_event(q.device, q.event)
end


function Base.copyto!(dst::DeviceStepQueue, src::StepQueue)
    # TODO: async?
    copyto!(dst.steps, src.steps)
    copyto!(dst.status, [src.length, src.pos, Int(src.can_continue)])
    return dst
end

function Base.copyto!(dst::StepQueue, src::DeviceStepQueue)
    # Don't retreive the steps, as it is up to the host to send them, the other way is wasteful.
    dst.length, dst.pos, dst.can_continue = Array(src.status)
    return dst
end


function process_queue!(queue::DeviceStepQueue, params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    basic_state = BasicSolverState(state)
    state_machine_func = state_machine_kernel(queue.device, params.workgroup_size)

    # TODO: the ndrange is quite important, the current choice is maybe sub-optimal since it includes all ghost cells
    state_machine_func(blk.device_data, basic_state, blk.size, queue; ndrange=block_size(blk))

    # Place an event in the device stream in order to be able to know when the kernel has completed,
    # independantly of the status of the stream.
    put_kernel_event(queue.device, queue.event)
    return
end


function perform_step(device, step::SolverStep.T, params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    # TODO: move this elsewhere
    must_wait = false
    if step == SolverStep.TimeStep
        must_wait = next_time_step(params, state, blk)

    elseif step == SolverStep.InitTimeStep
        must_wait = fetch_time_step(params, state, blk)

    elseif step == SolverStep.NewSweep
        must_wait = next_axis_sweep!(params, state)

    elseif step == SolverStep.EOS
        update_EOS!(params, state, blk)

    elseif step == SolverStep.Exchange
        must_wait = block_ghost_exchange(params, state, blk)

    elseif step == SolverStep.Fluxes
        numerical_fluxes!(params, state, blk)

    elseif step == SolverStep.CellUpdate
        cell_update!(params, state, blk)

    elseif step == SolverStep.RemapAdvection
        advection_fluxes!(params, state, blk)

    elseif step == SolverStep.RemapProjection
        euler_projection!(params, state, blk)

    else
        error("cannot perform step: ", step)
    end
    return !must_wait
end
