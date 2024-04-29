```@meta
CurrentModule = Armon
```

# Blocking

## Grid and blocks

```@docs
BlockGrid
grid_dimensions
TaskBlock
LocalTaskBlock
BlockData
RemoteTaskBlock
device_to_host!
host_to_device!
buffers_on_device
device_is_host
all_blocks(::BlockGrid)
block_idx
edge_block_idx
remote_block_idx
EdgeBlockRegions
RemoteBlockRegions
@iter_blocks
block_pos_containing_cell
block_origin
block_at
```

## Block size and iteration

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

## Block Tree

```@docs
BlockTree
all_blocks(::BlockTree)
is_leaf
depth
tree_block_count
visit_all_tree_block
iter_tree_block
apply_until_true
block_levels_for_machine
```

## Block states

```@docs
SolverState
first_state
SolverStep
block_state_machine
```

## Time step reduction

```@docs
next_time_step
GlobalTimeStep
TimeStepState.WaitingForMPI
TimeStepState.Done
TimeStepState.Ready
TimeStepState.DoingMPI
TimeStepState.AllContributed
```

## Block exchanges

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
