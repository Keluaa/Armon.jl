module DeadlockWatcher

using MPI
import InterProcessCommunication: IPC

export @detect_deadlocks

# TODO: sending SIGUSR1 seems to create a profile report with a backtrace, would this be useful?
# TODO: in Julia 1.11, we can interrupt threads: then in the handler, we could interrupt all other
#   threads and retrieve their backtrace?


struct TestSection
    id      :: Int
    timeout :: Int  # expected number of time (in µs) the section can take at maximum
    line    :: Int
    file    :: String
    label   :: String

    # members not shared with the deadlock watcher process
    no_gc     :: Bool
    backtrace :: Union{Nothing, Vector{Union{Ptr{Nothing}, Base.InterpreterIP}}}
end

TestSection(id, timeout, lnn::LineNumberNode, label; no_gc=false, bt=nothing) =
    TestSection(id, timeout, lnn.line, string(something(lnn.file, "<unknown file>")), label, no_gc, bt)


const ENABLE_DEADLOCK_WATCHER = parse(Bool, get(ENV, "DEADLOCK_WATCHER", "true"))

# TODO: add a PID after the keys in order to make them unique and allow multiple test runners at once
#  => e.g. '/<pid>/name'
const TESTS_STATUS_KEY = "/tests_status"
const SEM_TESTS_STATUS = Ref{IPC.Semaphore}()

const WATCHER_STATUS_KEY = "/watcher_status"
const SEM_WATCHER_STATUS = Ref{IPC.Semaphore}()

const SHM_TESTS_STATUS_KEY = "/tests_section"
const SHM_TESTS_STATUS_LEN = 4096
const SHM_TESTS_STATUS = Ref{IPC.SharedMemory}()

const WATCHER_PROCESS = Ref{Base.Process}()
const TESTS_PROCESS_PID = Ref{IPC.ProcessId}()

const DEADLOCK_SIGNAL = IPC.SIGUSR2

const SECTION_ID = Ref(0)
const SECTION_STACK = Vector{TestSection}()

const ORIGINAL_SIGACTION = Ref{IPC.SigAction}()


function load_state_from_shm(shm::IPC.SharedMemory)
    base_ptr = pointer(shm)

    id      = unsafe_load(Ptr{Int}(base_ptr));                    base_ptr += sizeof(Int)
    timeout = unsafe_load(Ptr{Int}(base_ptr));                    base_ptr += sizeof(Int)
    line    = unsafe_load(Ptr{Int}(base_ptr));                    base_ptr += sizeof(Int)

    file_str_len  = unsafe_load(Ptr{Int}(base_ptr));              base_ptr += sizeof(Int)
    label_str_len = unsafe_load(Ptr{Int}(base_ptr));              base_ptr += sizeof(Int)

    file  = unsafe_string(Ptr{UInt8}(base_ptr), file_str_len);    base_ptr += file_str_len
    label = unsafe_string(Ptr{UInt8}(base_ptr), label_str_len);   base_ptr += label_str_len

    return TestSection(id, timeout, LineNumberNode(line, file), label)
end


function store_state_to_shm!(shm::IPC.SharedMemory, section::TestSection)
    base_ptr = pointer(shm)

    unsafe_store!(Ptr{Int}(base_ptr), section.id);       base_ptr += sizeof(Int)
    unsafe_store!(Ptr{Int}(base_ptr), section.timeout);  base_ptr += sizeof(Int)
    unsafe_store!(Ptr{Int}(base_ptr), section.line);     base_ptr += sizeof(Int)

    file_len  = max(0, min(ncodeunits(section.file),  sizeof(shm) - 5*sizeof(Int)))
    label_len = max(0, min(ncodeunits(section.label), sizeof(shm) - 5*sizeof(Int) - file_len))
    unsafe_store!(Ptr{Int}(base_ptr), file_len);         base_ptr += sizeof(Int)
    unsafe_store!(Ptr{Int}(base_ptr), label_len);        base_ptr += sizeof(Int)

    unsafe_copyto!(Ptr{UInt8}(base_ptr), pointer(section.file), file_len);   base_ptr += file_len
    unsafe_copyto!(Ptr{UInt8}(base_ptr), pointer(section.label), label_len); base_ptr += label_len

    return
end


function remove_ipc_files()
    # When the main process is terminated brutally, the semaphores and shared memory files are not
    # properly cleaned up.
    # This will happen when `kill` or `MPI.Abort` are used.
    isassigned(SEM_TESTS_STATUS)   && finalize(SEM_TESTS_STATUS[])
    isassigned(SEM_WATCHER_STATUS) && finalize(SEM_WATCHER_STATUS[])
    isassigned(SHM_TESTS_STATUS)   && finalize(SHM_TESTS_STATUS[])

    rm(IPC.Semaphore, TESTS_STATUS_KEY)
    rm(IPC.Semaphore, WATCHER_STATUS_KEY)
    IPC.shmrm(SHM_TESTS_STATUS_KEY)
end


function is_test_process_alive()
    # When a parent process dies, its children are not killed automatically, and their parent process
    # changes.
    # See https://stackoverflow.com/a/2035683
    return IPC.getppid() == TESTS_PROCESS_PID[]
end


function wait_for_new_deadlock_section()
    # Wait for the test suite to enter a deadlock-protected region.
    # Pause every second to check if the test process is still alive.
    while true
        new_section = try
            timedwait(SEM_TESTS_STATUS[], 1)
            true
        catch e
            !(e isa IPC.TimeoutError) && rethrow(e)
            false
        end

        if new_section
            return
        elseif !is_test_process_alive()
            exit()
        end
    end
end


function monitor_for_deadlocks()
    while true
        wait_for_new_deadlock_section()

        test_section = load_state_from_shm(SHM_TESTS_STATUS[])
        test_section.id == -1 && break  # The main process wants to stop

        # Notify the test process that we have read the section state
        IPC.post(SEM_WATCHER_STATUS[])

        try
            timedwait(SEM_TESTS_STATUS[], test_section.timeout / 1e6)
        catch e
            # There is most likely a deadlock in the test process
            !(e isa IPC.TimeoutError) && rethrow(e)

            # Because test sections can be nested, we need to load it again to get the top-most section.
            test_section = load_state_from_shm(SHM_TESTS_STATUS[])

            handle_deadlock(test_section)

            break
        end
    end
end


function handle_deadlock(section::TestSection)
    if !is_test_process_alive()
        print("Deadlock watcher: test process exited unexpectedly\n")
        exit()
    end

    print("""
    !!!
    Detected a deadlock in process $(TESTS_PROCESS_PID[].value) after $(section.timeout / 1e6) sec
    Section n°$(section.id) at $(section.file):$(section.line), '$(section.label)'
    Invoking deadlock handler...
    !!!
    """)

    # Invoke the `deadlock_signal_handler` in the deadlocked process
    IPC.sigqueue(TESTS_PROCESS_PID[], DEADLOCK_SIGNAL)

    # Wait 5 sec for the process to handle the signal.
    ok = timedwait(() -> !is_test_process_alive(), 5; pollint=0.01)
    if ok === :timed_out
        # If the process is still alive, we must kill it
        print("""
        !!!
        Process $(TESTS_PROCESS_PID[].value) is not responding.
        Terminating...
        !!!
        """)
        try
            IPC.sigqueue(TESTS_PROCESS_PID[], IPC.SIGKILL)
        catch e
            !(e isa SystemError) && rethrow(e)
        end
    end

    try
        # If we kill the process owning the semaphores and shared memory, we must remove them manually.
        # This is also the case when the process encountered an error in its handler.
        remove_ipc_files()
    catch e
        !(e isa IPC.SystemError) && rethrow(e)
        # The files were properly disposed of by the main process
    end
end


function deadlock_signal_handler(signum::Cint)
    signum != DEADLOCK_SIGNAL && return
    # Note: if the GC was running when the signal handler started, then it is very likely that we
    # will deadlock if a GC safepoint is triggered (which is almost guarenteed).
    # Adding SIGSEGV to the signal's handler mask does not help, as a hot safepoint will instead
    # terminate the process when triggered.
    # Maybe there exist a GC-safe way of using signal handlers, but I don't know about it.
    # For now, disabling the GC at the begining of sections is the best option.

    if MPI.Initialized()
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
    else
        rank = 0
    end
    tid = Threads.threadid()

    # Print the error message as single string, to avoid interleaving of messages by different MPI
    # ranks.
    # `print` uses global locks and therefore can yield to the scheduler. Since we only have to
    # print a single simple string, we can use `jl_printf` instead (it calls C's `vasprintf`).
    # We minimize the number of possible problems by using classic C formats instead of Julia's
    # string convertions.
    format_complete   = "[rank %d - thread %d] Caught an deadlock at section n°%d in %s:%d, '%s'\n%s\n\n"
    format_no_section = "[rank %d - thread %d] Caught an deadlock at <no active section>\n\n"

    section = get(SECTION_STACK, length(SECTION_STACK), nothing)
    if isnothing(section)
        ccall(:jl_printf, Cint, (Ptr{Cvoid}, Cstring, Cint, Cint), stdout, format_no_section, rank, tid)
    else
        if isnothing(section.backtrace)
            bt_str = "<no backtrace>"
        else
            bt_str = sprint(Base.show_backtrace, section.backtrace; context=stdout)
        end

        ccall(:jl_printf, Cint,
            (Ptr{Cvoid}, Cstring, Cint, Cint, Cint, Cstring, Cint, Cstring, Cstring),
            stdout, format_complete, rank, tid, section.id, section.file, section.line, section.label, bt_str
        )
    end

    if MPI.Initialized()
        remove_ipc_files()
        MPI.Abort(MPI.COMM_WORLD, 1)
    end

    exit(1)
end

const DEADLOCK_HANDLER_FPTR = @cfunction(deadlock_signal_handler, Cvoid, (Cint,))


function begin_deadlock_protected_section(id, timeout, lnn::LineNumberNode, label; no_gc=false, section_backtrace=nothing)
    if Threads.threadid() != 1 || (MPI.Initialized() && MPI.Comm_rank(MPI.COMM_WORLD) != 0)
        no_gc && GC.enable(false)
        return  # Only the main thread of the root rank communicates with the watcher process
    end

    if !isassigned(SEM_TESTS_STATUS)
        # Start the deadlock watcher process
        setup_test_suite()
    end

    timeout = round(Int, timeout * 1e6)  # seconds to µs
    test_section = TestSection(id, timeout, lnn, label; no_gc, bt=section_backtrace)
    push!(SECTION_STACK, test_section)

    if length(SECTION_STACK) == 1
        # Start the deadlock protected section
        new_sigaction = IPC.SigAction(DEADLOCK_HANDLER_FPTR, IPC.SigSet(), 0)
        old_sigaction = IPC.SigAction()
        old_sigaction = IPC.sigaction!(DEADLOCK_SIGNAL, new_sigaction, old_sigaction)

        ORIGINAL_SIGACTION[] = old_sigaction

        store_state_to_shm!(SHM_TESTS_STATUS[], test_section)
        IPC.post(SEM_TESTS_STATUS[])

        # Wait for the watcher process to read the section
        timedwait(SEM_WATCHER_STATUS[], 1)
    end

    no_gc && GC.enable(false)

    return
end


function end_deadlock_protected_section(; no_gc=false)
    no_gc && GC.enable(true)

    if Threads.threadid() != 1 || isempty(SECTION_STACK)
        return
    end

    pop!(SECTION_STACK)
    if isempty(SECTION_STACK)
        IPC.post(SEM_TESTS_STATUS[])
        IPC.sigaction(DEADLOCK_SIGNAL, ORIGINAL_SIGACTION[])
    end

    return
end


"""
    @detect_deadlocks section_name [timeout [options...]] body

Terminates the Julia process (including all other MPI ranks, if any) if `body` does not complete
within `timeout` seconds (support floats, with a µs resolution, defaults to 5 seconds).

This is done using an external process, therefore it cannot be affected by Julia's task scheduler,
global IO locks, or the GC.
It uses Inter-Process Communications for communication with that watcher process, therefore it is
very reactive (~10-40 µs to create a new top-level section).

In MPI applications, only the root rank will start and communicate with the watcher process.
Sections in other ranks are ignored.

In multithreaded applications, only the main thread will communicate with the watcher process.
A good practice would be to wrap the whole parallel section (e.g. `Threads.@threads`) with
`@detect_deadlocks`.

MPI+threads applications are supported.

`section_name` and `timeout` can be variables. `section_name` supports string interpolation.

Options can be:
 - `:no_gc` to disable the Garbage Collector during the whole section
 - `:no_bt` to not include a backtrace to the begining of the section in the error message

Nesting of sections is supported. Options are not inherited (i.e. you need to use `:no_gc` in all
nested sections).
Only information about the top-most section is displayed when a deadlock is detected.

```julia
for i in 1:10
    @detect_deadlocks "my section \$i" 1.5 :no_gc begin
        # body...
    end
end
```
"""
macro detect_deadlocks(section_name, args...)
    options..., body = args

    if !ENABLE_DEADLOCK_WATCHER
        return esc(body)
    end

    if isempty(options)
        timeout = 5  # 5 seconds by default
    else
        timeout, options... = options
    end

    options = Set(map(options) do opt
        opt isa QuoteNode && return opt.value
        return value
    end)

    no_gc = :no_gc in options
    no_bt = :no_bt in options
    setdiff!(options, (:no_gc, :no_bt))
    if !isempty(options)
        err_str = "unknown options to `@detect_deadlocks`: $(join(options, ", "))"
        return esc(quote error($err_str) end)
    end

    bt = no_bt ? nothing : Expr(:call, backtrace)
    id = (SECTION_ID[] += 1)
    lnn = QuoteNode(__source__)

    return esc(quote
        $begin_deadlock_protected_section($id, $timeout, $lnn, $section_name; no_gc=$no_gc, section_backtrace=$bt)
        $(Expr(:tryfinally,
            body,
            quote  # finally
                $end_deadlock_protected_section(; no_gc=$no_gc) 
            end
        ))
    end)
end


function launch_watcher_process()
    watcher_cmd = `$(Base.julia_cmd()) -t 1 --project=$(Base.active_project()) $(@__FILE__)`
    watcher_cmd = pipeline(watcher_cmd; stdout=stdout, stderr=stderr)
    return run(watcher_cmd; wait=false)
end


function setup_watcher_process()
    Threads.nthreads() > 1 && error("deadlock_watcher.jl must only use a single thread")

    # The PID of the test suite's process
    TESTS_PROCESS_PID[] = IPC.getppid()

    SEM_TESTS_STATUS[] = IPC.Semaphore(TESTS_STATUS_KEY)
    SEM_WATCHER_STATUS[] = IPC.Semaphore(WATCHER_STATUS_KEY)
    SHM_TESTS_STATUS[] = IPC.SharedMemory(SHM_TESTS_STATUS_KEY)

    @info "Deadlock watcher process is running..."

    # Readyness handshake
    IPC.post(SEM_WATCHER_STATUS[])
    ok = timedwait(() -> SEM_WATCHER_STATUS[][] == 0, 1; pollint=0.01)
    ok === :timed_out && error("test process did not acknowledge the watcher process")

    monitor_for_deadlocks()
end


"""
    lower_process_priority()

If there is not enough cores for both the application process(es) and the watcher process, then
it is recommended to lower the priority of the application, so that is it more likely that the
watcher could get some CPU time to check if any deadlock has occured.
    
This is especially important when there is some hot loops in the application which would never let
the linux scheduler yield to a process with the same priority.

This must be done *after* [`@detect_deadlocks`](@ref) (or [`setup_test_suite`](@ref)) was called,
otherwise the watcher will still have the same priority.
"""
function lower_process_priority()
    PRIO_PROCESS = 0
    niceness = ccall(:getpriority, Cint, (Cint, IPC.ProcessId), PRIO_PROCESS, IPC.ProcessId(0))
    niceness += 1
    ok = ccall(:setpriority, Cint, (Cint, IPC.ProcessId, Cint), PRIO_PROCESS, IPC.ProcessId(0), niceness)
    ok != 0 && systemerror("setpriority")
end


function setup_test_suite()
    !ENABLE_DEADLOCK_WATCHER && return

    if isassigned(SEM_TESTS_STATUS) || isassigned(SHM_TESTS_STATUS)
        error("deadlock watcher already setup")
    end

    # Remove the IPC files if they exist already
    try
        remove_ipc_files()
    catch e
        !(e isa SystemError) && rethrow(e)
    end

    SEM_TESTS_STATUS[] = IPC.Semaphore(TESTS_STATUS_KEY, 0)
    SEM_WATCHER_STATUS[] = IPC.Semaphore(WATCHER_STATUS_KEY, 0)
    SHM_TESTS_STATUS[] = IPC.SharedMemory(SHM_TESTS_STATUS_KEY, SHM_TESTS_STATUS_LEN)

    WATCHER_PROCESS[] = launch_watcher_process()

    # Wait for the watcher process to be ready
    timedwait(SEM_WATCHER_STATUS[], 10)

    return
end


function stop_watcher_process()
    !ENABLE_DEADLOCK_WATCHER && return
    # ID of -1 means 'stop'
    stop_section = TestSection(-1, 0, LineNumberNode(0), "")
    store_state_to_shm!(SHM_TESTS_STATUS[], stop_section)
    IPC.post(SEM_TESTS_STATUS[])  # wake up the process
end

end


if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    DeadlockWatcher.setup_watcher_process()
else
    using .DeadlockWatcher
end
