module Armon

using EnumX
using Hwloc
using KernelAbstractions
using MPI
using MacroTools
using Polyester
using Preferences
using Printf
using ThreadPinning
using TimerOutputs

export ArmonParameters, BlockGrid, SolverStats, armon, data_type, memory_required
export device_to_host!, host_to_device!

# Forward declarations
abstract type Limiter end
abstract type RiemannScheme end
abstract type ProjectionScheme end
abstract type SplittingMethod end

include("utils.jl")
include("domain_ranges.jl")
include("tests.jl")
include("parameters.jl")
include("solver_state.jl")
include("blocking/blocking.jl")
include("profiling.jl")
include("generic_kernel.jl")
include("kernels.jl")
include("reductions.jl")
include("limiters.jl")
include("riemann_schemes.jl")
include("projection_schemes.jl")
include("axis_splitting.jl")
include("halo_exchange.jl")
include("io.jl")
include("solver.jl")

end
