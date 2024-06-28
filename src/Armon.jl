module Armon

using Printf
using Polyester
using ThreadPinning
using KernelAbstractions
using MPI
using MacroTools
using NUMA
using TimerOutputs
using Preferences
using EnumX
using Scotch

export ArmonParameters, BlockGrid, SolverStats, armon, data_type, memory_required
export device_to_host!, host_to_device!

# Forward declarations
abstract type Limiter end
abstract type RiemannScheme end
abstract type ProjectionScheme end
abstract type SplittingMethod end

include("utils.jl")
include("numa_utils.jl")
include("domain_ranges.jl")
include("tests.jl")
include("parameters.jl")
include("solver_state.jl")
include("profiling.jl")
include("generic_kernel.jl")
include("blocking/blocking.jl")
include("kernels.jl")
include("reductions.jl")
include("limiters.jl")
include("riemann_schemes.jl")
include("projection_schemes.jl")
include("axis_splitting.jl")
include("halo_exchange.jl")
include("io.jl")
include("logging.jl")
include("solver.jl")

end
