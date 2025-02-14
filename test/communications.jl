
using MPI
using ThreadPinning

include("deadlock_watcher.jl")

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
    !is_local_root && return

    if allunique(all_local_cores)
        if is_root && length(all_local_cores) == ncores()
            # All cores of the node with the root rank are used. This is problematic as the deadlock
            # watcher process could have trouble having enough CPU time to do its job correctly.
            DeadlockWatcher.setup_test_suite()
            DeadlockWatcher.lower_process_priority()  # make the root rank have lower priority over the deadlock process
        end
        return
    end

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


function MPI_try(f)
    try
        return f()
    catch e
        print_backtrace_and_abort(e, catch_backtrace())
    end
end


THREAD_BARRIER_VAR = Comms.Atomic{Int}(0)
BARRIER_CYCLE = Comms.Atomic{Int}(0)

function thread_barrier(tid, nthreads)
    cycle = @atomic BARRIER_CYCLE.x
    @atomic THREAD_BARRIER_VAR.x += 1

    guard = 0
    while (@atomic THREAD_BARRIER_VAR.x) != nthreads && (@atomic BARRIER_CYCLE.x) == cycle
        ccall(:jl_cpu_pause, Cvoid, ())
        GC.safepoint()
        guard += 1
        if guard % 1_000_000 == 0
            @error "stuck in barrier n°$cycle: $(@atomic THREAD_BARRIER_VAR.x)"
            exit(1)
        end
    end

    # reset
    if tid == 1
        @atomic THREAD_BARRIER_VAR.x = 0
        @atomic BARRIER_CYCLE.x += 1
    else
        guard = 0
        while (@atomic BARRIER_CYCLE.x) == cycle
            ccall(:jl_cpu_pause, Cvoid, ())
            GC.safepoint()
            guard += 1
            if guard % 1_000_000 == 0
                @error "stuck in post barrier n°$cycle: $(@atomic THREAD_BARRIER_VAR.x)"
                exit(1)
            end
        end
    end
end


test_point_to_point(model, test_name) = test_point_to_point(model, test_name, 1, 1)

function test_point_to_point(model, test_name, tid, nthreads)
    # Basic 2-way ring exchange. Should work with any number of ranks ≥ 1
    prev_rank = mod(global_rank - 1, global_size)
    next_rank = mod(global_rank + 1, global_size)
    side_pos = tid
    array_type = Vector{Int}
    buffer_size = 10
    total_buffer_size = buffer_size * nthreads
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

    @detect_deadlocks "exchange init of $test_name" 1 begin
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

    expected_size = Comms.uses_global_buffers(model) ? total_buffer_size : buffer_size
    @test all(length.(Comms.unsafe_send_buffer(xchg_prev)) .== expected_size)
    @test all(length.(Comms.unsafe_recv_buffer(xchg_prev)) .== expected_size)
    @test all(length.(Comms.unsafe_send_buffer(xchg_next)) .== expected_size)
    @test all(length.(Comms.unsafe_recv_buffer(xchg_next)) .== expected_size)

    if nthreads == 1
        # `send_completed` may return `false` if another thread acquired a global lock, so those tests
        # are single-thread only
        @test Comms.send_completed(xchg_prev)
        @test Comms.send_completed(xchg_next)
    end

    if !Comms.is_async(model)
        # `recv_completed` is expected to work the same as `send_completed` only for synchronous communications
        @test Comms.recv_completed(xchg_prev)
        @test Comms.recv_completed(xchg_next)
    end

    for rep in 1:10
        # Note: it is expected that receives are always preceded by a send.
        @detect_deadlocks "send to rank prev $prev_rank for $test_name" 1 begin
            buf = Comms.acquire_send_buffer!(xchg_prev)
            buf .= global_rank + (tid - 1) * 1000 + rep * 100000
            Comms.release_send_buffer!(xchg_prev)
        end

        @detect_deadlocks "send to rank next $next_rank for $test_name" 1 begin
            while (buf = Comms.try_acquire_send_buffer!(xchg_next); isnothing(buf))
                ccall(:jl_cpu_pause, Cvoid, ())
                GC.safepoint()
            end
            buf .= global_rank + (tid - 1) * 1000 + rep * 100000
            Comms.release_send_buffer!(xchg_next)
        end

        @detect_deadlocks "receive from prev rank $prev_rank for $test_name" 1 begin
            # deadlock
            # mpiexec -np 2 --oversubscribe -- julia --project=.. -t 7 --color=yes ./runtests.jl comms
            buf = Comms.acquire_recv_buffer!(xchg_prev)
            tmp_buf_prev .= buf
            Comms.release_recv_buffer!(xchg_prev)
        end

        @detect_deadlocks "receive from next rank $next_rank for $test_name" 1 begin
            while (buf = Comms.try_acquire_recv_buffer!(xchg_next); isnothing(buf))
                ccall(:jl_cpu_pause, Cvoid, ())
                GC.safepoint()
            end
            tmp_buf_next .= buf
            Comms.release_recv_buffer!(xchg_next)
        end

        expected_buf_prev = zeros(Int, size(tmp_buf_prev))
        expected_buf_next = zeros(Int, size(tmp_buf_next))
        if model isa Comms.NoCommunicationModel
            # Here we receive the data we sent
            expected_buf_prev .= global_rank + (tid - 1) * 1000 + rep * 100000
            expected_buf_next .= global_rank + (tid - 1) * 1000 + rep * 100000
        else
            expected_buf_prev .= prev_rank + (tid - 1) * 1000 + rep * 100000
            expected_buf_next .= next_rank + (tid - 1) * 1000 + rep * 100000
        end

        @test tmp_buf_prev == expected_buf_prev
        @test tmp_buf_next == expected_buf_next

        @detect_deadlocks "wait send completed for prev rank $prev_rank" 1 begin
            guard = 0
            while !Comms.send_completed(xchg_prev)
                ccall(:jl_cpu_pause, Cvoid, ())
                GC.safepoint()
                guard += 1
                if guard % 1_000_000 == 0
                    # deadlock
                    # mpiexec -np 2 --oversubscribe -- julia --project=.. -t 7 --color=yes ./runtests.jl comms
                    Comms.log("stuck waiting send (1)"; extra=Comms.part_state_str(xchg_prev))
                end
            end
        end

        if nthreads > 1
            # `send_completed` may return `false` if another thread acquired a global lock
            @detect_deadlocks "test send completed for prev rank $prev_rank" 1 begin
                guard = 0
                while !Comms.send_completed(xchg_prev)
                    ccall(:jl_cpu_pause, Cvoid, ())
                    GC.safepoint()
                    guard += 1
                    if guard % 1_000_000 == 0
                        Comms.log("stuck waiting send (2)"; extra=Comms.part_state_str(xchg_prev))
                    end
                end
            end
        else
            @test Comms.send_completed(xchg_prev)
        end

        @detect_deadlocks "wait send completed for next rank $next_rank" 1 begin
            Comms.wait_send_completed(xchg_next)
        end

        if nthreads > 1
            # `send_completed` may return `false` if another thread acquired a global lock
            @detect_deadlocks "test send completed for next rank $next_rank" 1 begin
                while !Comms.send_completed(xchg_next)
                    ccall(:jl_cpu_pause, Cvoid, ())
                    GC.safepoint()
                end
            end
        else
            @test Comms.send_completed(xchg_next)
        end

        if !Comms.is_async(model)
            # Trivial for synchronous models, so no deadlock detection needed
            @test Comms.recv_completed(xchg_prev)
            @test Comms.recv_completed(xchg_next)
        end

        thread_barrier(tid, nthreads)
    end

    return xchg_prev, xchg_next
end


function test_point_to_point_multithreaded(model, test_name)
    @detect_deadlocks "parallel $test_name with $(Threads.nthreads()) threads" 20 :no_gc begin
        Threads.@threads :static for tid in 1:Threads.nthreads()
            MPI_try() do
                xchgs = test_point_to_point(model, test_name, tid, Threads.nthreads())

                thread_barrier(tid, Threads.nthreads())

                # Important: since we keep the same tags for each exchange, the next MPI exchange
                # might collide with this one if it uses permanently active receive requests (or similar).
                Comms.finalize_comm!.(xchgs)
                # Comms.finalize_comm!(xchgs)
            end
        end
    end
end


function test_collective(model, test_name)
    # Basic global reduction. Should work with any number of ranks.
    reduction_op = max
    array_type = Vector{Int}
    buffer_size = global_size
    tmp_buf = Vector{Int}(undef, buffer_size)
    tmp_buf .= -1

    reduc = nothing
    @detect_deadlocks "exchange init of $test_name" 1 begin
        reduc = Comms.init_reduce_broadcast(model, reduction_op, array_type, buffer_size)
    end
    @test !isnothing(reduc)

    @test all(length.(Comms.unsafe_send_buffer(reduc)) .== buffer_size)
    @test all(length.(Comms.unsafe_recv_buffer(reduc)) .== buffer_size)

    @test Comms.send_completed(reduc)
    @test Comms.send_completed(reduc)
    if !Comms.is_async(model)
        # `recv_completed` is expected to work the same as `send_completed` only for synchronous communications
        @test Comms.recv_completed(reduc)
        @test Comms.recv_completed(reduc)
    end

    # Note: it is expected that receives are always preceded by a send.
    @detect_deadlocks "collective send for $test_name" 1 begin
        buf = Comms.acquire_send_buffer!(reduc)
        buf .= 0
        buf[global_rank + 1] = global_rank
        Comms.release_send_buffer!(reduc)
    end

    @detect_deadlocks "collective recv for $test_name" 1 begin
        buf = Comms.acquire_recv_buffer!(reduc)
        tmp_buf .= buf
        Comms.release_recv_buffer!(reduc)
    end

    expected_buf = zeros(Int, size(tmp_buf))
    if model isa Comms.NoCommunicationModel
        # Here we receive the data we sent
        expected_buf[global_rank + 1] = global_rank
    else
        expected_buf .= 0:global_size-1
    end

    @test tmp_buf == expected_buf

    @detect_deadlocks "collective (try) send for $test_name" 1 begin
        while (buf = Comms.try_acquire_send_buffer!(reduc); isnothing(buf)) end
        buf .= 0
        buf[global_rank + 1] = global_rank
        Comms.release_send_buffer!(reduc)
    end

    @detect_deadlocks "collective (try) recv for $test_name" 1 begin
        while (buf = Comms.try_acquire_recv_buffer!(reduc); isnothing(buf)) end
        tmp_buf .= buf
        Comms.release_recv_buffer!(reduc)
    end

    @test tmp_buf == expected_buf

    # TODO: wait_recv_completed (but NOT for the test before the last recv_completed, as receive request might be permanent)

    @detect_deadlocks "collective wait send completed" 1 begin
        Comms.wait_send_completed(reduc)
    end
    @test Comms.send_completed(reduc)

    if !Comms.is_async(model)
        # Trivial for synchronous models, so no deadlock detection needed
        @test Comms.recv_completed(reduc)
    end

    return reduc
end


function test_collective_multithreaded(comm_model, comm_model_name)
    # TODO
end


function test_model(comm_model, comm_model_name)
    if Comms.supports_point_to_point(comm_model)
        MPI_try() do
            xchgs = test_point_to_point(comm_model, comm_model_name)
            # Important: since we keep the same tags for each exchange, the next MPI exchange might
            # collide with this one if it uses permanently active receive requests (or similar).
            Comms.finalize_comm!.(xchgs)
        end

        if Comms.is_thread_safe(comm_model)
            MPI_try() do
                test_point_to_point_multithreaded(comm_model, comm_model_name)
            end
        end
    end

    if Comms.supports_collectives(comm_model)
        MPI_try() do
            xchg = test_collective(comm_model, comm_model_name)
            Comms.finalize_comm!(xchg)
        end

        if Comms.is_thread_safe(comm_model)
            MPI_try() do
                test_collective_multithreaded(comm_model, comm_model_name)
            end
        end
    end
end


@testset "Communications" verbose=true begin
    check_for_core_overlaps()

    @testset "$comm_model_name model" verbose=true for (comm_model_name, comm_model_kwargs) in (
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
            is_root && @warn "Model $comm_model_name cannot be tested" maxlog=1
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
