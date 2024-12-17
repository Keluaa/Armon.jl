```@meta
CurrentModule = Armon
```

# I/O

## General interface

```@docs
write_sub_domain_file
read_sub_domain_file!
```

## Abstract I/O interface

This interface allows the solver to support multiple output formats.

```@docs
AbstractSolverIO
supports_mpi
supports_threads
supports_temporal_data
file_extension
domain_writer
domain_reader
Base.close(::AbstractSolverIO)
write_domain_to_file
write_block_to_file
read_domain_from_file
read_block_from_file
```

## CSV

Writing to CSV files is supported as a basic output format.
Pass `:csv` to [`domain_writer`](@ref) or [`ArmonParameters`](@ref) to use it.

When using MPI, one file per rank is written.

Options:

- `precision`: number of digits (base 10) to write. Default is enough for an exact representation.
- `all_ghosts`: read/write all ghost cells to the file. Defaults to `false`. A file written with ghost cells
  can only be read back with ghost cells, and vice-versa.
- `global_ghosts`: same as `all_ghosts`, but only for the ghost cells at the border of the global domain.
- `for_3D`: write one blank line between cell rows, for compatiblity with pm3d of Gnuplot. Defaults to `true`.
- `vars`: tuple of `Symbol`s of variables to write. Defaults to `saved_vars()`.

## HDF5

HDF5 output is supported, using the VTKHDF file structure to allow direct visualization in Paraview.
Pass `:hdf5` to [`domain_writer`](@ref) or [`ArmonParameters`](@ref) to use it.
The HDF5 support is an extension loaded only if the [`HDF5.jl`](https://github.com/JuliaIO/HDF5.jl) package is loaded.

Depending on the configuration of your HDF5 libraries:

- MPI support is required, as it will enable multiple ranks to read and write to the same file.
- Thread-safety is optional: it allows blocks to be saved to the same file concurrently, without needing a synchronization of all threads.

Options:

- `vars`: tuple of `Symbol`s of variables to write. Defaults to `saved_vars()`.
