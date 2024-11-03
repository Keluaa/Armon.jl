```@meta
CurrentModule = Armon
```

# I/O

## Interface

```@docs
AbstractSolverIO
supports_mpi
supports_threads
supports_temporal_data
domain_writer
```

## CSV

Writing to CSV files is supported, as a basic output format.
Pass `:csv` to [`domain_writer`](@ref) or [`ArmonParameters`](@ref) to use it.

When using MPI, one file per rank is written.

Options:
- `precision`: number of digits to write. Default is enough for an exact representation.

## HDF5

HDF5 output is supported, using the VTKHDF file structure to allow direct visualization in Paraview.
Pass `:hdf5` to [`domain_writer`](@ref) or [`ArmonParameters`](@ref) to use it.
The HDF5 support is an extension, loaded only if the [`HDF5.jl`](https://github.com/JuliaIO/HDF5.jl) package is loaded.

Depending on the configuration of your HDF5 libraries:
- MPI support is required, as it will enable multiple ranks to read and write to the same file.
- Thread-safety is optional: it allows blocks to be saved to the same file concurrently, without needing a synchronization of all threads.
