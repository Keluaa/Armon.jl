module Armon

using Printf
using Polyester
using ThreadPinning
using KernelAbstractions
using KernelsToolkit
using MPI
using MacroTools
using NUMA
using TimerOutputs
using Preferences
using EnumX
using Scotch

export ArmonParameters, BlockGrid, SolverStats, armon

# Forward declarations
abstract type TestCase end
abstract type Limiter end
abstract type RiemannScheme end
abstract type ProjectionScheme end
abstract type SplittingMethod end

@kernels_metadata

include("utils.jl")
include("numa_utils.jl")
include("parameters.jl")
include("tests.jl")
include("solver_state.jl")
include("profiling.jl")
include("performance_macros.jl")
include("blocking/blocking.jl")
include("kernels.jl")
include("reductions.jl")
include("limiters.jl")
include("riemann_schemes.jl")
include("projection_schemes.jl")
include("axis_splitting.jl")
include("halo_exchange.jl")
include("io.jl")
include("io_csv.jl")
include("compare.jl")
include("logging.jl")
include("solver.jl")

function __init__()
    @register_all_kernels
end

end
