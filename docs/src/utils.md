```@meta
CurrentModule = Armon
```

## Axis

```@docs
Axis
axes_of
next_axis
```

## Side

```@docs
Side
sides_of
sides_along
first_side
first_sides
last_side
last_sides
axis_of
opposite_of
offset_to
side_from_offset
Neighbours
```

## NUMA utilities

```@docs
array_pages
touch_pages
move_pages(::Vector{Ptr{T}}, ::Any) where T
lock_pages(::Ptr, ::Any)
unlock_pages
```

## Other

```@docs
SolverException
@section
```
