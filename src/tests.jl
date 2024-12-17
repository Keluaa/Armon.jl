
abstract type TwoStateTestCase <: TestCase end
struct Sod{A}    <: TwoStateTestCase end
struct Sod_circ  <: TwoStateTestCase end
struct Bizarrium <: TwoStateTestCase end
struct Sedov{T}  <: TwoStateTestCase
    r::T
    center_cells_count::Int
end

create_test(::ArmonParameters, ::Type{Test}) where {Test <: TestCase} = Test()

function create_test(params::ArmonParameters{T}, ::Type{Sedov}) where {T}
    # We want to place the energy at the center of the domain.
    # `/2` to only include the center of cells at the center of the domain.
    Δx = params.domain_size ./ params.N ./ 2
    r_Sedov = T(hypot(Δx...))

    # Adjustment to make sure rounding errors don't exclude a cell from the center
    r_Sedov *= one(T) + 100*eps(T)

    # The energy will be distributed to 2 cells in directions where there is an even number
    # of cells in the mesh, up to `2^dim` cells. This way we don't have to compute
    # intersections of cell regions with the center sphere of radius `r_Sedov`.
    center_cells_count = 2^count(iseven, params.N)

    return Sedov{T}(r_Sedov, center_cells_count)
end

test_from_name(::Val{:Sod})       = Sod{Axis.X}
test_from_name(::Val{:Sod_x})     = Sod{Axis.X}
test_from_name(::Val{:Sod_y})     = Sod{Axis.Y}
test_from_name(::Val{:Sod_z})     = Sod{Axis.Z}
test_from_name(::Val{:Sod_circ})  = Sod_circ
test_from_name(::Val{:Bizarrium}) = Bizarrium
test_from_name(::Val{:Sedov})     = Sedov

test_from_name(::Val{s}) where s = solver_error(:config, "Unknown test case: '$s'")
test_from_name(s::Symbol) = test_from_name(Val(s))

test_name(::Test) where {Test <: TestCase} = nameof(Test)

default_domain_size(::Type{<:TestCase}, D) = ntuple(Returns(1), D)
default_domain_size(::Type{Sedov}, D) = ntuple(Returns(2), D)  # 2×2 domain centered at (0,0)

default_domain_origin(::Type{<:TestCase}, D) = ntuple(Returns(0), D)
default_domain_origin(::Type{Sedov}, D) = ntuple(Returns(1), D)

default_CFL(::Union{Sod, Sod_circ}) = 0.95
default_CFL(::Bizarrium) = 0.6
default_CFL(::Sedov) = 0.7

default_max_time(::Union{Sod, Sod_circ}) = 0.20
default_max_time(::Bizarrium) = 80e-6
default_max_time(::Sedov) = 1.0

specific_heat_ratio(::TestCase) = 7/5  # Di-atomic perfect gas

is_conservative(::TestCase) = true
is_conservative(::Bizarrium) = false

has_source_term(::TestCase) = false

Base.show(io::IO, ::Sod{A}) where {A} = print(io, "Sod shock tube (along the $(Symbol(A)) axis)")
Base.show(io::IO, ::Sod_circ)  = print(io, "Sod shock tube (high-region in centered circle)")
Base.show(io::IO, ::Bizarrium) = print(io, "Bizarrium")
Base.show(io::IO, ::Sedov)     = print(io, "Sedov")

test_region_high(x::Tuple{Vararg{T}}, ::Sod{A})    where {T, A} = x[A] ≤ 0.5
test_region_high(x::Tuple{Vararg{T}}, ::Sod_circ)  where {T} = sum((x .- T(0.5)).^2) ≤ T(0.3)^2
test_region_high(x::Tuple{Vararg{T}}, ::Bizarrium) where {T} = x[1] ≤ 0.5
test_region_high(x::Tuple{Vararg{T}}, s::Sedov{T}) where {T} = sum(x.^2) ≤ s.r^2


struct InitTestParamsTwoState{T, D}
    high_ρ::T
    low_ρ::T
    high_E::T
    low_E::T
    high_u::NTuple{D, T}
    low_u::NTuple{D, T}

    function InitTestParamsTwoState(;
        high_ρ::T, low_ρ::T, high_E::T, low_E::T, high_u::NTuple{D, T}, low_u::NTuple{D, T}
    ) where {T, D}
        new{T, D}(high_ρ, low_ρ, high_E, low_E, high_u, low_u)
    end
end


function init_test_params(::Union{Sod, Sod_circ}, ::Type{T}, D) where {T}
    return InitTestParamsTwoState(
        high_ρ = T(1.),
         low_ρ = T(0.125),
        high_E = T(2.5),
         low_E = T(2.0),
        high_u = ntuple(Returns(zero(T)), D),
         low_u = ntuple(Returns(zero(T)), D),
    )
end

function init_test_params(::Bizarrium, ::Type{T}, D) where {T}
    return InitTestParamsTwoState(
        high_ρ = T(1.42857142857e+4),
         low_ρ = T(10000.),
        high_E = T(4.48657821135e+6),
         low_E = T(0.5 * 250^2),
        high_u = ntuple(Returns(zero(0)), D),
         low_u = ntuple(d -> d == 1 ? T(250.) : zero(T), D),
    )
end

function init_test_params(p::Sedov, ::Type{T}, D) where {T}
    return InitTestParamsTwoState(
        high_ρ = T(1.),
         low_ρ = T(1.),
        # E so that the blast wave reaches r=1 at t=1, and spread among the cells at the center of the mesh
        # See https://en.wikipedia.org/wiki/Taylor%E2%80%93von_Neumann%E2%80%93Sedov_blast_wave#Mathematical_description
        high_E = T((1/1.033)^5 / p.center_cells_count),
         low_E = T(2.5e-14),
        high_u = ntuple(Returns(zero(T)), D),
         low_u = ntuple(Returns(zero(T)), D),
    )
end


@enum BC FreeFlow Dirichlet


function boundary_condition(test, side::Side.T, D::Val{Dim}, ::Type{T}) where {Dim, T}
    condition = boundary_condition(test, D)[side]
    if condition == FreeFlow
        return ntuple(Returns(one(T)), Dim)
    else  # if condition == Dirichlet
        # mirror along the axis of side
        return ifelse.(axis_of(side) .== axes_of(Dim), -one(T), one(T))
    end
end


function boundary_condition(::Sod{A}, dim) where {A}
    return Neighbours(dim) do axis, _
        return axis == A ? Dirichlet : FreeFlow
    end
end


boundary_condition(::Sod_circ, dim) = Neighbours(Returns(Dirichlet), dim)


function boundary_condition(::Bizarrium, dim)
    return Neighbours(dim) do _, side
        return side == Side.Right ? FreeFlow : Dirichlet
    end
end


boundary_condition(::Sedov, dim) = Neighbours(Returns(FreeFlow), dim)

#
# DebugIndexes
#

struct DebugIndexes <: TestCase end

test_from_name(::Val{:DebugIndexes}) = DebugIndexes

default_CFL(::DebugIndexes) = 0
default_max_time(::DebugIndexes) = 0

Base.show(io::IO, ::DebugIndexes) = print(io, "DebugIndexes")

boundary_condition(::DebugIndexes, dim) = Neighbours(Returns(Dirichlet), dim)
