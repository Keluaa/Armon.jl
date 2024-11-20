```@meta
CurrentModule = Armon
```

# Utilities

```@docs
Axis
Side
SolverException
@section
```

## NUMA utilities

```@docs
array_pages
touch_pages
move_pages(::Vector{Ptr{T}}, ::Any) where T
lock_pages(::Ptr, ::Any)
unlock_pages
```
