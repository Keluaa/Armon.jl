```@meta
CurrentModule = Armon
```

# Block states

```@docs
SolverState
first_state
SolverStep
block_state_machine
solver_cycle_async
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
