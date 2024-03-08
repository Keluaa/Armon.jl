
struct SequentialSplitting     <: SplittingMethod end
struct GodunovSplitting        <: SplittingMethod end
struct StrangSplitting         <: SplittingMethod end
struct SinglePassSplitting{Ax} <: SplittingMethod end

splitting_from_name(::Val{:Sequential})    = SequentialSplitting()
splitting_from_name(::Val{:Godunov})       = GodunovSplitting()
splitting_from_name(::Val{:SequentialSym}) = GodunovSplitting()
splitting_from_name(::Val{:Strang})        = StrangSplitting()
splitting_from_name(::Val{:X_only})        = SinglePassSplitting{X_axis}()
splitting_from_name(::Val{:Y_only})        = SinglePassSplitting{Y_axis}()

splitting_from_name(::Val{s}) where {s} = solver_error(:config, "Unknown splitting method: '$s'")
splitting_from_name(s::Symbol) = splitting_from_name(Val(s))

Base.show(io::IO, ::SequentialSplitting) = print(io, "Sequential (X, Y ; X, Y)")
Base.show(io::IO, ::GodunovSplitting)    = print(io, "Godunov (X, Y ; Y, X)")
Base.show(io::IO, ::StrangSplitting)     = print(io, "Strong (½X, Y, ½X ; ½Y, X, ½Y)")
Base.show(io::IO, ::SinglePassSplitting{Ax}) where {Ax} = print(io, "Single pass ($Ax ; $Ax)")

split_axes(state::SolverState{T}) where {T} = split_axes(state.splitting, T, state.global_dt.cycle)

function split_axes(::SequentialSplitting, ::Type{T}, _) where {T}
    return ((X_axis, T(1.0)), (Y_axis, T(1.0)))
end

function split_axes(::GodunovSplitting, ::Type{T}, cycle) where {T}
    if iseven(cycle)
        return ((X_axis, T(1.0)), (Y_axis, T(1.0)))
    else
        return ((Y_axis, T(1.0)), (X_axis, T(1.0)))
    end
end

function split_axes(::StrangSplitting, ::Type{T}, cycle) where {T}
    if iseven(cycle)
        return ((X_axis, T(0.5)), (Y_axis, T(1.0)), (X_axis, T(0.5)))
    else
        return ((Y_axis, T(0.5)), (X_axis, T(1.0)), (Y_axis, T(0.5)))
    end
end

function split_axes(::SinglePassSplitting{Ax}, ::Type{T}, _) where {Ax, T}
    return ((Ax, T(1.0)),)
end
