
!isinteractive() && @info "Loading..."

using Armon
import Armon: SolverException
using Test

WRITE_FAILED = parse(Bool, get(ENV, "WRITE_FAILED", "false"))
NO_MPI = parse(Bool, get(ENV, "NO_MPI", "false"))
MPI_PROCS = parse(Int, get(ENV, "MPI_PROCS", "0"))

if !NO_MPI
    using MPI
    MPI.Init(; threadlevel=:multiple)
    is_root = MPI.Comm_rank(MPI.COMM_WORLD) == 0
    world_size = MPI.Comm_size(MPI.COMM_WORLD)
    if is_root && world_size > 1
        println("Testing using $world_size processes")
    end
else
    is_root = true
end


include(joinpath(@__DIR__, "reference_data/reference_functions.jl"))


if isinteractive()
    menu = """
    Tests available:
     - all            All tests below
     - short          Equivalent to 'quality, stability, convergence, conservation'
     - quality        Code quality
     - stability      Type stability
     - domains        Domain 2D indexing
     - blocking       Blocking tests
     - convergence    Convergence to the reference solutions
     - conservation   Check that the energy and mass for each are kept constant throughout a lot of cycles.
     - GPU            Equivalence of the GPU backends (CUDA & ROCm) with the CPU
     - kokkos         Equivalence of the Kokkos backend with the Julia CPU backend
     - MPI            Equivalence with the single domain case and asynchronous communications.
                      If 'TEST_KOKKOS_MPI=true' or 'kokkos' is in the test list, MPI tests with
                      Kokkos will also be performed. Same for GPU backends.

    Separate multiple test sets with a comma.

    Choice: """
    printstyled(stdout, menu; color=:light_green)
    raw_main_options = readline()
elseif haskey(ENV, "TEST_FILTER")
    raw_main_options = ENV["TEST_FILTER"]
else
    raw_main_options = isempty(ARGS) ? "short" : join(ARGS, ", ")
end

main_options = split(raw_main_options, ',') .|> strip .|> lowercase
filter!(!isempty, main_options)
main_options = main_options .|> Symbol |> union

if :all in main_options
    expanded_options = [:quality, :stability, :domains, :blocking, :convergence, :conservation,
                        :kokkos, :gpu, :mpi]
elseif :short in main_options
    expanded_options = [:quality, :stability, :convergence, :conservation]
else
    expanded_options = []
end

deleteat!(main_options, findall(opt -> opt == :all || opt == :short, main_options))
append!(main_options, expanded_options)
union!(main_options)

ENV["TEST_KOKKOS_MPI"] = parse(Bool, get(ENV, "TEST_KOKKOS_MPI", "false")) || (:kokkos in main_options)
TEST_KOKKOS_BACKEND = ENV["TEST_KOKKOS_BACKEND"] = get(ENV, "TEST_KOKKOS_BACKEND", "OPENMP")
TEST_KOKKOS_PATH = get(ENV, "TEST_KOKKOS_PATH", nothing)


function do_tests(tests_to_do)
    is_root && println("Testing: ", join(tests_to_do, ", "))

    !is_root && (Test.TESTSET_PRINT_ENABLE[] = false)

    run_file(file_name) = include(joinpath(@__DIR__, file_name))

    ts = @testset "Armon.jl" begin
        for test in tests_to_do
            if !is_root
                if test === :mpi            run_file("mpi.jl")
                else
                    # the test is for a single process only
                end
            elseif test === :quality        run_file("code_quality.jl")
            elseif test === :stability      run_file("type_stability.jl")
            elseif test === :domains        run_file("domains.jl")
            elseif test === :blocking       run_file("blocking.jl")
            elseif test === :convergence    run_file("convergence.jl")
            elseif test === :conservation   run_file("conservation.jl")
            elseif test === :gpu            run_file("gpu.jl")
            elseif test === :kokkos         run_file("kokkos.jl")
            elseif test === :mpi
                if MPI_PROCS > 0 && world_size < MPI_PROCS
                    @info "Launching $MPI_PROCS MPI sub-processes"
                    project_dir = abspath(@__DIR__, "..")
                    mpi_cmd = `$(mpiexec()) -n $MPI_PROCS $(Base.julia_cmd()) --project=$(project_dir) --color=yes $(@__FILE__)`
                    @testset "MPI" begin
                        @test success(mpi_cmd)
                    end
                else
                    run_file("mpi.jl")
                end
            else
                error("Unknown test set: $test")
            end

            # TODO: susceptibility test comparing a result with different rounding modes
        end
    end

    if is_root && isinteractive()
        Test.print_test_results(ts)
    end

    kokkos_id = Base.PkgId(Base.UUID("3296cea9-b0de-4b57-aba0-ce554b517c3b"), "Kokkos")
    if !isinteractive() && haskey(Base.loaded_modules, kokkos_id)
        Base.loaded_modules[kokkos_id].finalize()
    end
end


!isempty(main_options) && do_tests(main_options)
