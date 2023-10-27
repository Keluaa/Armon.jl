# Armon

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://keluaa.github.io/Armon.jl/dev)
[![Build Status](https://github.com/Keluaa/Armon.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Keluaa/Armon.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Keluaa/Armon.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Keluaa/Armon.jl)

Armon is a 2D CFD solver for compressible non-viscous fluids, using the finite volume method.

It was made to explore Julia's capabilities in HPC and for performance portability: it should
perform very well on any CPU and GPU.
Domain decomposition using MPI is supported.

The twin project [Armon-Kokkos](https://github.com/Keluaa/Armon-Kokkos) is a mirror of the core of
this solver (with much less options) written in C++ using the Kokkos library.
It is possible to reuse kernels from that solver in this one, using the
[Kokkos.jl](https://github.com/Keluaa/Kokkos.jl) package.
