
abstract type TestCase end

abstract type TwoStateTestCase <: TestCase end
struct Sod       <: TwoStateTestCase end
struct Sod_y     <: TwoStateTestCase end
struct Sod_circ  <: TwoStateTestCase end
struct Bizarrium <: TwoStateTestCase end
struct Sedov{T}  <: TwoStateTestCase
    r::T
end

create_test(::NTuple, ::Type{Test}) where {Test <: TestCase} = Test()

function create_test(Δx::NTuple, ::Type{Sedov})
    T = eltype(Δx)
    r_Sedov::T = hypot(Δx...) / sqrt(2)
    return Sedov{T}(r_Sedov)
end

test_from_name(::Val{:Sod})       = Sod
test_from_name(::Val{:Sod_y})     = Sod_y
test_from_name(::Val{:Sod_circ})  = Sod_circ
test_from_name(::Val{:Bizarrium}) = Bizarrium
test_from_name(::Val{:Sedov})     = Sedov

test_from_name(::Val{s}) where s = solver_error(:config, "Unknown test case: '$s'")
test_from_name(s::Symbol) = test_from_name(Val(s))

test_name(::Test) where {Test <: TestCase} = nameof(Test)

default_domain_size(::Type{<:TestCase}) = (1, 1)
default_domain_size(::Type{Sedov}) = (2, 2)

default_domain_origin(::Type{<:TestCase}) = (0, 0)
default_domain_origin(::Type{Sedov}) = (-1, -1)

default_CFL(::Union{Sod, Sod_y, Sod_circ}) = 0.95
default_CFL(::Bizarrium) = 0.6
default_CFL(::Sedov) = 0.7

default_max_time(::Union{Sod, Sod_y, Sod_circ}) = 0.20
default_max_time(::Bizarrium) = 80e-6
default_max_time(::Sedov) = 1.0

specific_heat_ratio(::TestCase) = 7/5

is_conservative(::TestCase) = true
is_conservative(::Bizarrium) = false

has_source_term(::TestCase) = false

Base.show(io::IO, ::Sod)       = print(io, "Sod shock tube")
Base.show(io::IO, ::Sod_y)     = print(io, "Sod shock tube (along the Y axis)")
Base.show(io::IO, ::Sod_circ)  = print(io, "Sod shock tube (cylindrical symmetry around the Z axis)")
Base.show(io::IO, ::Bizarrium) = print(io, "Bizarrium")
Base.show(io::IO, ::Sedov)     = print(io, "Sedov")

test_region_high(x::Tuple{Vararg{T}}, ::Sod)       where {T} = x[1] ≤ 0.5
test_region_high(x::Tuple{Vararg{T}}, ::Sod_y)     where {T} = x[2] ≤ 0.5
test_region_high(x::Tuple{Vararg{T}}, ::Sod_circ)  where {T} = sum((x .- T(0.5)).^2) ≤ T(0.09)  # radius of 0.3 
test_region_high(x::Tuple{Vararg{T}}, ::Bizarrium) where {T} = x[1] ≤ 0.5
test_region_high(x::Tuple{Vararg{T}}, s::Sedov{T}) where {T} = sum(x.^2) ≤ s.r^2


struct InitTestParamsTwoState{T}
    high_ρ::T
    low_ρ::T
    high_E::T
    low_E::T
    high_u::T
    low_u::T
    high_v::T
    low_v::T

    function InitTestParamsTwoState(;
        high_ρ::T, low_ρ::T, high_E::T, low_E::T, high_u::T, low_u::T, high_v::T, low_v::T
    ) where {T}
        new{T}(high_ρ, low_ρ, high_E, low_E, high_u, low_u, high_v, low_v)
    end
end


function init_test_params(::Union{Sod, Sod_y, Sod_circ}, ::Type{T}) where {T}
    return InitTestParamsTwoState(
        high_ρ = T(1.),
         low_ρ = T(0.125),
        high_E = T(2.5),
         low_E = T(2.0),
        high_u = zero(T),
         low_u = zero(T),
        high_v = zero(T),
         low_v = zero(T)
    )
end

function init_test_params(::Bizarrium, ::Type{T}) where {T}
    return InitTestParamsTwoState(
        high_ρ = T(1.42857142857e+4),
         low_ρ = T(10000.),
        high_E = T(4.48657821135e+6),
         low_E = T(0.5 * 250^2),
        high_u = zero(T),
         low_u = T(250.),
        high_v = zero(T),
         low_v = zero(T)
    )
end

function init_test_params(p::Sedov, ::Type{T}) where {T}
    return InitTestParamsTwoState(
        high_ρ = T(1.),
         low_ρ = T(1.),
        high_E = T((1/1.033)^5 / (π * p.r^2)),  # E so that the blast wave reaches r=1 at t=1 (E is spread in a circle of radius `p.r`)
         low_E = T(2.5e-14),
        high_u = zero(T),
         low_u = zero(T),
        high_v = zero(T),
         low_v = zero(T)
    )
end


@enum BC FreeFlow Dirichlet


struct Boundaries
    left::BC
    right::BC
    bottom::BC
    top::BC

    Boundaries(; left, right, bottom, top) = new(left, right, bottom, top)
end


function Base.getindex(bounds::Boundaries, side::Side.T)
    return if side == Side.Left
        bounds.left
    elseif side == Side.Right
        bounds.right
    elseif side == Side.Bottom
        bounds.bottom
    else
        bounds.top
    end
end


function boundary_condition(test, side::Side.T)::NTuple{2, Int}
    condition = boundary_condition(test)[side]
    if condition == FreeFlow
        return (1, 1)
    else  # if condition == Dirichlet
        if side in (Side.Left, Side.Right)
            return (-1, 1)  # mirror along X
        else
            return (1, -1)  # mirror along Y
        end
    end
end


function boundary_condition(::Sod)
    return Boundaries(
        left   = Dirichlet,
        right  = Dirichlet,
        bottom = FreeFlow,
        top    = FreeFlow
    )
end


function boundary_condition(::Sod_y)
    return Boundaries(
        left   = FreeFlow,
        right  = FreeFlow,
        bottom = Dirichlet,
        top    = Dirichlet
    )
end


function boundary_condition(::Sod_circ)
    return Boundaries(
        left   = Dirichlet,
        right  = Dirichlet,
        bottom = Dirichlet,
        top    = Dirichlet
    )
end


function boundary_condition(::Bizarrium)
    return Boundaries(
        left   = Dirichlet,
        right  = FreeFlow,
        bottom = Dirichlet,
        top    = Dirichlet
    )
end


function boundary_condition(::Sedov)
    return Boundaries(
        left   = FreeFlow,
        right  = FreeFlow,
        bottom = FreeFlow,
        top    = FreeFlow
    )
end

#
# DebugIndexes
#

struct DebugIndexes <: TestCase end

test_from_name(::Val{:DebugIndexes}) = DebugIndexes

default_CFL(::DebugIndexes) = 0
default_max_time(::DebugIndexes) = 0

Base.show(io::IO, ::DebugIndexes) = print(io, "DebugIndexes")

function boundary_condition(::DebugIndexes)
    return Boundaries(
        left   = Dirichlet,
        right  = Dirichlet,
        bottom = Dirichlet,
        top    = Dirichlet
    )
end
