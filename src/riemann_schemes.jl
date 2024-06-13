
struct RiemannGodunov <: RiemannScheme end
struct RiemannGAD     <: RiemannScheme end

scheme_from_name(::Val{:Godunov}) = RiemannGodunov()
scheme_from_name(::Val{:GAD})     = RiemannGAD()

scheme_from_name(::Val{s}) where s = solver_error(:config, "Unknown scheme: '$s'")
scheme_from_name(s::Symbol) = scheme_from_name(Val(s))

Base.show(io::IO, ::RiemannGodunov) = print(io, "Godunov")
Base.show(io::IO, ::RiemannGAD)     = print(io, "GAD")

uses_limiter(::RiemannGodunov) = false
uses_limiter(::RiemannGAD)     = true

stencil_width(::RiemannGodunov) = 1
stencil_width(::RiemannGAD)     = 2


@kernel_function function acoustic_Godunov(
    ρᵢ::T, ρᵢ₋₁::T, cᵢ::T, cᵢ₋₁::T,
    uᵢ::T, uᵢ₋₁::T, pᵢ::T, pᵢ₋₁::T
) where T
    rc_l = ρᵢ₋₁ * cᵢ₋₁
    rc_r = ρᵢ   * cᵢ
    uˢᵢ = (rc_l * uᵢ₋₁ + rc_r * uᵢ +               (pᵢ₋₁ - pᵢ)) / (rc_l + rc_r)
    pˢᵢ = (rc_r * pᵢ₋₁ + rc_l * pᵢ + rc_l * rc_r * (uᵢ₋₁ - uᵢ)) / (rc_l + rc_r)
    return uˢᵢ, pˢᵢ
end


@generic_kernel function acoustic!(
    s::Int, uˢ_::V, pˢ_::V,
    ρ::V, uₐ::V, p::V, c::V
) where V
    u = uₐ  # `u` or `v` depending on the current axis
    i = @index_2D_lin()
    uˢ_[i], pˢ_[i] = acoustic_Godunov(
        ρ[i], ρ[i-s], c[i], c[i-s],
        u[i], u[i-s], p[i], p[i-s]
    )
end


function numerical_fluxes!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, ::RiemannGodunov)
    range = block_domain_range(blk.size, state.steps_ranges.fluxes)
    blk_data = block_device_data(blk)
    u = state.axis == Axis.X ? blk_data.u : blk_data.v
    s = stride_along(blk.size, state.axis)
    return acoustic!(params, blk_data, range, s, blk_data.uˢ, blk_data.pˢ, u)
end


@generic_kernel function acoustic_GAD!(
    s::Int, dt::T, dx::T,
    uˢ::V, pˢ::V, ρ::V, uₐ::V, p::V, c::V,
    lim::LimiterType
) where {T, V <: AbstractArray{T}, LimiterType <: Limiter}
    i = @index_2D_lin()

    u = uₐ  # `u` or `v` depending on the current axis

    # First order acoustic solver on the left cell
    uˢ_i₋, pˢ_i₋ = acoustic_Godunov(
        ρ[i-s], ρ[i-2s], c[i-s], c[i-2s],
        u[i-s], u[i-2s], p[i-s], p[i-2s]
    )

    # First order acoustic solver on the current cell
    uˢ_i, pˢ_i = acoustic_Godunov(
        ρ[i], ρ[i-s], c[i], c[i-s],
        u[i], u[i-s], p[i], p[i-s]
    )

    # First order acoustic solver on the right cell
    uˢ_i₊, pˢ_i₊ = acoustic_Godunov(
        ρ[i+s], ρ[i], c[i+s], c[i],
        u[i+s], u[i], p[i+s], p[i]
    )

    # Second order GAD acoustic solver on the current cell

    r_u₋ = (uˢ_i₊  -  u[i]) / (uˢ_i - u[i-s] + T(1e-6))
    r_p₋ = (pˢ_i₊  -  p[i]) / (pˢ_i - p[i-s] + T(1e-6))
    r_u₊ = (u[i-s] - uˢ_i₋) / (u[i] - uˢ_i   + T(1e-6))
    r_p₊ = (p[i-s] - pˢ_i₋) / (p[i] - pˢ_i   + T(1e-6))

    r_u₋ = limiter(r_u₋, lim)
    r_p₋ = limiter(r_p₋, lim)
    r_u₊ = limiter(r_u₊, lim)
    r_p₊ = limiter(r_p₊, lim)

    dm_l = ρ[i-s] * dx
    dm_r = ρ[i]   * dx
    Dm   = (dm_l + dm_r) / 2

    rc_l = ρ[i-s] * c[i-s]
    rc_r = ρ[i]   * c[i]
    θ    = T(0.5) * (1 - (rc_l + rc_r) / 2 * (dt / Dm))

    uˢ[i] = uˢ_i + θ * (r_u₊ * (u[i] - uˢ_i) - r_u₋ * (uˢ_i - u[i-s]))
    pˢ[i] = pˢ_i + θ * (r_p₊ * (p[i] - pˢ_i) - r_p₋ * (pˢ_i - p[i-s]))
end


function numerical_fluxes!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, ::RiemannGAD)
    range = block_domain_range(blk.size, state.steps_ranges.fluxes)
    blk_data = block_device_data(blk)
    u = state.axis == Axis.X ? blk_data.u : blk_data.v
    s = stride_along(blk.size, state.axis)
    return acoustic_GAD!(params, blk_data, range, s, state.dt, state.dx, u, state.riemann_limiter)
end


numerical_fluxes!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock) =
    numerical_fluxes!(params, state, blk, state.riemann_scheme)

function numerical_fluxes!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    @iter_blocks for blk in grid
        numerical_fluxes!(params, state, blk)
    end
end
