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
