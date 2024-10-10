
using MPI
using ThreadPinning

if !MPI.Initialized()
    MPI.Init(; threadlevel=:multiple)
end

MPI.Barrier(MPI.COMM_WORLD)
global_rank = MPI.Comm_rank(MPI.COMM_WORLD)
global_size = MPI.Comm_size(MPI.COMM_WORLD)
is_root = global_rank == 0

# Local thread pinning
# This ensures that two Julia threads will run concurrently, without the help of the kernel.
# While it is important for performance, the issue here is that we need it for deadlock detection,
# as if two Julia threads are on the same core and one is busy waiting, no detection can take place.
node_local_comm = MPI.Comm_split_type(MPI.COMM_WORLD, MPI.COMM_TYPE_SHARED, global_rank)
local_rank = MPI.Comm_rank(node_local_comm)
cores = (0:Threads.nthreads()-1) .+ local_rank * Threads.nthreads()
pinthreads(cores; warn=false)

Comms = Armon.Communications


function check_for_core_overlaps()
    # If two threads of different ranks are pinned to the same core, then deadlocks are near unavoidable.
    rank_cores = collect(cores)
    is_local_root = MPI.Comm_rank(node_local_comm) == 0
    all_local_cores = MPI.Gather(rank_cores, node_local_comm)
    all_local_ranks = MPI.Gather(global_rank, node_local_comm)  # local to global rank
    (!is_local_root || allunique(all_local_cores)) && return

    node_name = gethostname()
    core_count = maximum(all_local_cores)
    if core_count > ThreadPinning.ncores()
        # We are pinning to hyperthreads as there is too many threads per rank: same deadlock issue.
        println("[$global_rank] ranks $(join(all_local_ranks, ", ", " and ")) at node '$node_name' \
                 use too many threads and are therefore pinned to hyperthreads.")
    end

    assigned_cores = zeros(Int, core_count)
    for (i, core) in enumerate(all_local_cores)
        local_rank = mod1(i, length(rank_cores))
        if assigned_cores[core] == 0
            assigned_cores[core] = local_rank
        else
            println("[$global_rank] core n°$core of node '$node_name' is already assigned to rank \
                     $(all_local_ranks[local_rank])")
        end
    end

    MPI.Abort(MPI.COMM_WORLD, 3)
end


function print_backtrace_and_abort(error, bt=backtrace())
    # Print the error message as single string, to avoid interleaving of messages by different MPI ranks
    rank_str = "[$global_rank] caught an error: "
    err_str = "ERROR: " * sprint(showerror, error; context=stdout)
    bt_str = sprint(Base.show_backtrace, bt; context=stdout)
    error_msg = rank_str * "\n" * err_str * "\n" * bt_str * "\n\n"

    # `print` uses global locks and therefore can yield to the scheduler. Since we only have to
    # print a single simple string, we can use `jl_printf` instead (it calls C's `vasprintf`).
    ccall(:jl_printf, Cint, (Ptr{Cvoid}, Cstring), stdout, error_msg)

    # Before aborting, let other MPI ranks print their backtrace if they detected the deadlock as
    # well. Once again, no `Base.sleep` as it yields to Julia's scheduler.
    Libc.systemsleep(0.25)
    MPI.Abort(MPI.COMM_WORLD, 2)
end


function detect_deadlock(f; label=nothing, timeout=2)
    # `f` can use multiple Julia threads, but must leave 1 (or 2?) available. `timeout` is in seconds.
    # TODO: for some very obscure reason, I get deadlocks with `timedwait` (and IO like `println`)
    #   when I only have 2 threads. 3 threads seems to be the minimum for this to work.
    if Threads.nthreads() < 3
        is_root && @warn "cannot detect deadlocks: 3 Julia threads or more are required" maxlog=1
        f()
        return
    end

    # The interesting backtrace is here, not in the deadlock guard task.
    # It would also be interesting to print the backtrace of the worker task, but it is not possible
    # as it may be in Julia or C code. When using GDB, `jl_backtracet` can print to `stdout` the
    # backtrace of any task, but it calls `printf` for each line, which is not ideal for MPI, and is
    # unsafe if any Julia code is running, which we cannot know.
    bt = backtrace()

    # We must disable the GC as if the guard task starts a GC pass while the other is waiting
    # forever in some C code, we will never get out of the deadlock.
    GC.enable(false)

    ok = Armon.Atomic{Bool}(false)
    guard_started = Threads.Event()
    @sync begin
        worker_task = Threads.@spawn begin
            wait(guard_started)
            f()
            @atomic ok.x = true
        end

        Threads.@spawn begin
            notify(guard_started)

            # In order to prevent this task from ever yielding to the Julia's scheduler, we cannot
            # use any task-programing construct here, as it would create an opportunity for the
            # scheduler to place only worker tasks on all threads. If this happens, there is no
            # longer any guard which can trigger the abort when too much time has elapsed.
            # This also includes all IO operations, such as `print`, as they involve a global stream
            # lock, which can yield to the scheduler.
            ns_timeout = timeout * 1e9
            start = time_ns()
            timed_out = true
            while true
                if (@atomic ok.x) || istaskfailed(worker_task)
                    timed_out = false
                    break
                elseif (time_ns() - start) > ns_timeout
                    break
                end
            end

            if timed_out
                label_str = isnothing(label) ? "" : (" at '" * label * "'")
                print_backtrace_and_abort(ErrorException("deadlock$label_str"), bt)
            end
        end
    end

    GC.enable(true)

    return
end


function MPI_try(f)
    try
        return f()
    catch e
        print_backtrace_and_abort(e, catch_backtrace())
    end
end


function test_point_to_point(model, test_name)
    # Basic 2-way ring exchange. Should work with any number of ranks ≥ 1
    prev_rank = mod(global_rank - 1, global_size)
    next_rank = mod(global_rank + 1, global_size)
    side_pos = 1
    array_type = Vector{Int}
    buffer_size = 10
    total_buffer_size = buffer_size
    tmp_buf_prev = Vector{Int}(undef, buffer_size)
    tmp_buf_next = Vector{Int}(undef, buffer_size)
    tmp_buf_prev .= -1
    tmp_buf_next .= -1

    if prev_rank == next_rank
        # Two exchanges with the same rank require using a different side, otherwise communications
        # could interact badly with each other for some models.
        side_1 = 1
        side_2 = 2
    else
        side_1 = side_2 = 1
    end

    xchgs = Any[nothing, nothing]

    detect_deadlock(; label="exchange init of $test_name") do
        xchgs[1] = Comms.init_exchange(model, prev_rank, side_1, side_pos, array_type, buffer_size, total_buffer_size)
        xchgs[2] = Comms.init_exchange(model, next_rank, side_2, side_pos, array_type, buffer_size, total_buffer_size)
    end

    xchg_prev, xchg_next = xchgs
    @test !isnothing(xchg_prev) && !isnothing(xchg_next)

    if !Comms.is_async(model) && isodd(global_rank)
        # Swap 'prev' and 'next' so that even ranks exchange with the previous rank first, and odd
        # ranks exchange with the next rank first, preventing any deadlock for synchronous communications.
        prev_rank,    next_rank    = next_rank,    prev_rank
        tmp_buf_prev, tmp_buf_next = tmp_buf_next, tmp_buf_prev
        xchg_prev,    xchg_next    = xchg_next,    xchg_prev
    end

    @test all(length.(Comms.unsafe_send_buffer(xchg_prev)) .== buffer_size)
    @test all(length.(Comms.unsafe_recv_buffer(xchg_prev)) .== buffer_size)
    @test all(length.(Comms.unsafe_send_buffer(xchg_next)) .== buffer_size)
    @test all(length.(Comms.unsafe_recv_buffer(xchg_next)) .== buffer_size)

    @test Comms.send_completed(xchg_prev)
    @test Comms.send_completed(xchg_next)
    if !Comms.is_async(model)
        # `recv_completed` is expected to work the same as `send_completed` only for synchronous communications
        @test Comms.recv_completed(xchg_prev)
        @test Comms.recv_completed(xchg_next)
    end

    # Note: it is expected that receives are always preceded by a send.
    detect_deadlock(; label="send to rank prev $prev_rank for $test_name") do
        buf = Comms.acquire_send_buffer!(xchg_prev)
        buf .= global_rank
        Comms.release_send_buffer!(xchg_prev)
    end

    detect_deadlock(; label="send to rank next $next_rank for $test_name") do
        while (buf = Comms.try_acquire_send_buffer!(xchg_next); isnothing(buf)) end
        buf .= global_rank
        Comms.release_send_buffer!(xchg_next)
    end

    detect_deadlock(; label="receive from prev rank $prev_rank for $test_name") do
        buf = Comms.acquire_recv_buffer!(xchg_prev)
        tmp_buf_prev .= buf
        Comms.release_recv_buffer!(xchg_prev)
    end

    detect_deadlock(; label="receive from next rank $next_rank for $test_name") do
        buf = Comms.acquire_recv_buffer!(xchg_next)
        tmp_buf_next .= buf
        Comms.release_recv_buffer!(xchg_next)
    end

    # TODO: try_acquire_send_buffer!
    # TODO: try_acquire_recv_buffer!
    # TODO: wait_recv_completed (but NOT for the test before the last recv_completed, as receive request might be permanent)

    expected_buf_prev = zeros(Int, size(tmp_buf_prev))
    expected_buf_next = zeros(Int, size(tmp_buf_next))
    if model isa Comms.NoCommunicationModel
        # Here we receive the data we sent
        expected_buf_prev .= global_rank
        expected_buf_next .= global_rank
    else
        expected_buf_prev .= prev_rank
        expected_buf_next .= next_rank
    end

    @test tmp_buf_prev == expected_buf_prev
    @test tmp_buf_next == expected_buf_next

    detect_deadlock(; label="wait send completed for prev rank $prev_rank") do
        Comms.wait_send_completed(xchg_prev)
    end
    @test Comms.send_completed(xchg_prev)
    
    detect_deadlock(; label="wait send completed for next rank $next_rank") do
        Comms.wait_send_completed(xchg_next)
    end
    @test Comms.send_completed(xchg_next)

    if !Comms.is_async(model)
        # Trivial for synchronous models, so no deadlock detection needed
        @test Comms.recv_completed(xchg_prev)
        @test Comms.recv_completed(xchg_next)
    end

    # Important: since we keep the same tags for each exchange, the next MPI exchange might collide
    # with this one if it uses permanently active receive requests (or similar).
    # By finalizing the object, any active request will be cancelled.
    finalize(xchg_prev)
    finalize(xchg_next)
end


function test_point_to_point_multithreaded(model, test_name)
    # TODO
end


function test_collective(model, test_name)
    # TODO
end


function test_collective_multithreaded(comm_model, comm_model_name)
    # TODO
end


function test_model(comm_model, comm_model_name)
    if Comms.supports_point_to_point(comm_model)
        MPI_try() do
            test_point_to_point(comm_model, comm_model_name)
        end

        if Comms.is_thread_safe(comm_model)
            MPI_try() do
                test_point_to_point_multithreaded(comm_model, comm_model_name)
            end
        end
    end

    if Comms.supports_collectives(comm_model)
        MPI_try() do
            test_collective(comm_model, comm_model_name)
        end

        if Comms.is_thread_safe(comm_model)
            MPI_try() do
                test_collective_multithreaded(comm_model, comm_model_name)
            end
        end
    end
end


@testset "Communications" begin
    check_for_core_overlaps()

    @testset "$comm_model_name model" for (comm_model_name, comm_model_kwargs) in (
        (:no_comms, (;)),
        (:sync, (;)),
        (:async, (;)),
        (:async_safe, (;)),
        (:rma, (;)),
        (:partitioned, (; partition_size=10)),
    )
        comm_model = try
            Comms.communication_model(comm_model_name, MPI.COMM_WORLD; comm_model_kwargs...)
        catch e
            # The model cannot be created (e.g. partitioned comms unsupported)
            @test true skip=true
            continue
        end

        @testset let MPI_rank = global_rank
            for r in 1:1000
                MPI.Barrier(MPI.COMM_WORLD)
                test_model(comm_model, comm_model_name)
            end
        end
    end
end
