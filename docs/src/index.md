```@meta
CurrentModule = Armon
```

# Armon

Documentation for [Armon](https://github.com/Keluaa/Armon.jl).

Armon is a 2D CFD solver for compressible non-viscous fluids, using the finite volume method.

It was made to explore Julia's capabilities in HPC and for performance portability: it should
perform very well on any CPU and GPU.
Domain decomposition using MPI is supported.

The twin project [Armon-Kokkos](https://github.com/Keluaa/Armon-Kokkos) is a mirror of the core of
this solver (with much less options) written in C++ using the Kokkos library.
It is possible to reuse kernels from that solver in this one, using the
[Kokkos.jl](https://github.com/Keluaa/Kokkos.jl) package.

## Parameters and entry point

```@docs
ArmonParameters
armon
SolverStats
StepsRanges
data_type
```

## Grid and blocks

```@docs
BlockGrid
grid_dimensions
TaskBlock
LocalTaskBlock
BlockData
RemoteTaskBlock
var_arrays
var_arrays_names
device_to_host!
host_to_device!
buffers_on_device
device_is_host
all_blocks
block_idx
edge_block_idx
remote_block_idx
EdgeBlockRegions
RemoteBlockRegions
@iter_blocks
block_pos_containing_cell
block_origin
block_at
block_size_at
move_pages(::BlockGrid)
lock_pages(::BlockGrid)
```

### Block size and iteration

```@docs
BlockSize
StaticBSize
DynamicBSize
block_size
real_block_size
ghosts
border_domain
ghost_domain
block_domain_range
position
lin_position
in_grid
is_ghost
BlockRowIterator
DomainRange
```

## Block states

```@docs
SolverState
first_state
SolverStep
block_state_machine
```

### Time step reduction

```@docs
next_time_step
GlobalTimeStep
TimeStepState.WaitingForMPI
TimeStepState.Done
TimeStepState.Ready
TimeStepState.DoingMPI
TimeStepState.AllContributed
```

### Block exchanges

```@docs
block_ghost_exchange
start_exchange
finish_exchange
BlockInterface
mark_ready_for_exchange!
exchange_done!
side_is_done!
is_side_done
BlockExchangeState
BlockExchangeState.NotReady
BlockExchangeState.InProgress
BlockExchangeState.Done
```

### Block distribution

```@docs
thread_workload_distribution
simple_workload_distribution
scotch_grid_partition
block_grid_from_workload
```

## Device and backends

```@docs
CPU_HP
create_device
init_backend
device_memory_info
memory_info
memory_required
```

## Kernels

```@docs
@generic_kernel
@kernel_init
@kernel_options
@index_1D_lin
@index_2D_lin
@iter_idx
@simd_loop
@simd_threaded_iter
@simd_threaded_loop
@threaded
@threads
```

## Logging

```@docs
BlockLogEvent
ThreadLogEvent
collect_logs
analyse_log_stats
BlockGridLogStats
BLOCK_LOG_THREAD_LOCAL_STORAGE
```
