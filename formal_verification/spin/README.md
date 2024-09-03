
# Formal Verification of Armon.jl using SPIN

The [SPIN model checker](https://github.com/nimble-code/Spin) is used to perform an approximate
verification of the structure of parallelism in the solver.

From a description of the block grid topology across all MPI processes, multi-threading behaviour
and MPI calls are reproduced by the Promela model.
Only the structure and logic surrounding the compute kernels of the actual solver are replicated:
no floating-point arithmetic is done, and the verification is independent of the kernels' implementation.

SPIN is used to check if no deadlock or assertion violation can happen during any execution of the
solver, taking into account any possible order of operations (atomic ops, MPI calls, etc...).
However, since the state space is too large, this verification is partial: only a subset of all
reachable states are explored.

To maximize the coverage of reachable states in a reasonable time period, the
[`Swarm`](https://github.com/nimble-code/Swarm) tool is used.

## Verifying the solver

The model requires the topology of the block grid to verify. It can be generated using the
`gen_grid.jl` Julia script by giving it a file containing the parameters of the grid (size of the
grid, number of threads, number of processes, etc...):

```shell
julia ./gen_grid.jl ./grids/default_grid_3x3_4threads.toml
```

See `./grids/default_grid_3x3_4threads.toml` to find information about each parameter.

### Swarm verification

The verification can then be launched with the command `make swarm`.
By default it uses 128 processes with 1GiB of memory for one hour.
Different parameters can be given:

```shell
make clean_swarm  # remove any preexisting swarm file
# use 64 processes with 512MB for 30 min
make SWARM_PROC=64 SWARM_TIME=0.5 SWARM_MEM=512M swarm
```

Use `make help` for an explaination of the different `makefile` targets and variables.

Among the many verification runs done by SPIN, if any finds a problem, a trail file is generated in
`./swarm_check/trails/'`. Use `make <same options> swarm_trail` to generate the files needed to
inspect them. An empty folder after a swarm run means no issues where found.

Note that while a swarm run is deterministic (i.e. repeated calls to `make swarm` will yield the
same results), changing the Swarm parameters will change how the state space is explored and may
give different results.

### Non-progress verification

`make check_non_progress` will verify the absence of non-progress cycles (deadlocks).
This verification is incompatible with Swarm therefore it can only be done by a single process,
making it very slow.
