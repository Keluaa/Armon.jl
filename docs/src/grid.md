```@meta
CurrentModule = Armon
```

# Grid

```@docs
BlockGrid
grid_dimensions
```

## Grid indexing

```@docs
block_idx
edge_block_idx
remote_block_idx
block_pos_containing_cell
block_origin
block_at
block_size_at
```

## Iteration

```@docs
all_blocks
@iter_blocks
EdgeBlockRegions
RemoteBlockRegions
BlockRowIterator
```

## Host/Device transfers

```@docs
device_to_host!
host_to_device!
buffers_on_device
device_is_host
```

## NUMA

```@docs
move_pages(::BlockGrid)
lock_pages(::BlockGrid)
```

## Block distribution

```@docs
thread_workload_distribution
simple_workload_distribution
scotch_grid_partition
block_grid_from_workload
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
