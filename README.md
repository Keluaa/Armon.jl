# Armon.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://keluaa.github.io/Armon.jl/dev)
[![Build Status](https://github.com/Keluaa/Armon.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Keluaa/Armon.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Keluaa/Armon.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Keluaa/Armon.jl)

Armon.jl is an experimental CFD solver for compressible non-viscous fluids, using a finite volume method.
It supports domains in any dimension, although it is mainly used and tested for 1, 2 and 3D simulations.

It was developed to explore Julia's capabilities in HPC and for performance portability: it should
perform very well on any CPU (tested on Nvidia Grace CPU, AMD Milan 7763, Intel Skylake Xeon Platinum 8168...)
and GPU (tested on Nvidia A100, AMD MI250, Intel Max 1100...) without any additional configuration.
Domain decomposition is done using MPI, and cache-blocking is performed in a single process for
higher performance on CPU.

The twin project [Armon-Kokkos](https://github.com/Keluaa/Armon-Kokkos) is a mirror of the core of
this solver (in 2D, with much less options, without MPI or cache-blocking) written in C++ using the Kokkos library.
It is possible to reuse kernels from that solver in this one, using the
[Kokkos.jl](https://github.com/Keluaa/Kokkos.jl) package.

Armon.jl was created in order to create a portable CFD benchmark.
The [`benchmark`](https://github.com/Keluaa/Armon.jl/tree/benchmark) branch is the state of the
repository for which most of the results for an upcoming paper were obtained.

## Installation

Armon.jl is an unregistered Julia package: to install it, clone this repository
(`using Pkg; Pkg.develop(url="https://github.com/Keluaa/Armon.jl")`).

## Usage

Parameters are detailed in the help of `ArmonParameters`.

```julia
julia> using Armon

julia> params = ArmonParameters(use_MPI=false, options...);

julia> armon(params);  # Run the solver
```
