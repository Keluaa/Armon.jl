
struct SequentialSplitting     <: SplittingMethod end
struct GodunovSplitting        <: SplittingMethod end
struct StrangSplitting         <: SplittingMethod end
struct SinglePassSplitting{Ax} <: SplittingMethod end

splitting_from_name(::Val{:Sequential})    = SequentialSplitting()
splitting_from_name(::Val{:Godunov})       = GodunovSplitting()
splitting_from_name(::Val{:SequentialSym}) = GodunovSplitting()
splitting_from_name(::Val{:Strang})        = StrangSplitting()
splitting_from_name(::Val{:X_only})        = SinglePassSplitting{Axis.X}()
splitting_from_name(::Val{:Y_only})        = SinglePassSplitting{Axis.Y}()
splitting_from_name(::Val{:Z_only})        = SinglePassSplitting{Axis.Z}()

splitting_from_name(::Val{s}) where {s} = solver_error(:config, "Unknown splitting method: '$s'")
splitting_from_name(s::Symbol) = splitting_from_name(Val(s))

Base.show(io::IO, ::SequentialSplitting) = print(io, "Sequential (X, Y ; X, Y)")
Base.show(io::IO, ::GodunovSplitting)    = print(io, "Godunov (X, Y ; Y, X)")
Base.show(io::IO, ::StrangSplitting)     = print(io, "Strang (½X, Y, ½X ; ½Y, X, ½Y)")
Base.show(io::IO, ::SinglePassSplitting{Ax}) where {Ax} = print(io, "Single pass ($Ax ; $Ax)")

# TODO: passing `dim` as a type in `Val(dim)` would allow to generate more optimal code and prevent any allocations
split_axes(state::SolverState{T}, dim) where {T} = split_axes(state.splitting, T, dim, state.global_dt.cycle)

function split_axes(::SequentialSplitting, ::Type{T}, dim, _) where {T}
    return ntuple(d -> (Axis.T(d), one(T)), dim)  # `(X, 1); (Y, 1); (Z, 1)` etc...
end

function split_axes(::GodunovSplitting, ::Type{T}, dim, cycle) where {T}
    split_params = ntuple(d -> (Axis.T(d), one(T)), dim)  # `(X, 1); (Y, 1); (Z, 1)` etc...
    # Circular shift on a tuple
    first_axes = split_params[(dim - cycle % dim + 1):end]
    last_axes  = split_params[1:(dim - cycle % dim)]
    return ((first_axes..., last_axes...))
end

function split_axes(::StrangSplitting, ::Type{T}, dim, cycle) where {T}
    # Circular permutation of the axes
    first_axes  = axes_of(dim)[(dim - cycle % dim + 1):end]
    last_axes   = axes_of(dim)[1:(dim - cycle % dim)]
    strang_axes = (first_axes..., last_axes...)

    # The last axis is swept only once and has a time factor of 1, while the others are
    # swept again in the reverse order.
    split_params = ntuple(d -> (strang_axes[d], ifelse(d == dim, one(T), one(T)/2)), dim)
    
    # `½X, Y, ½X` in 2D, or `½X ½Y Z ½Y ½X` in 3D (for cycle 0)
    return (Base.front(split_params)..., split_params[end], reverse(Base.front(split_params))...,)
end

function split_axes(::SinglePassSplitting{Ax}, ::Type{T}, _, _) where {Ax, T}
    return ((Ax, one(T)),)
end
