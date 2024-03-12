
struct EulerProjection    <: ProjectionScheme end
struct Euler2ndProjection <: ProjectionScheme end

scheme_from_name(::Val{:euler})     = EulerProjection()
scheme_from_name(::Val{:euler_2nd}) = Euler2ndProjection()

Base.show(io::IO, ::EulerProjection)    = print(io, "1ˢᵗ order projection")
Base.show(io::IO, ::Euler2ndProjection) = print(io, "2ⁿᵈ order projection")

stencil_width(::EulerProjection)    = 1
stencil_width(::Euler2ndProjection) = 2


@kernel_function function slope_minmod(uᵢ₋::T, uᵢ::T, uᵢ₊::T, r₋::T, r₊::T) where T
    Δu₊ = r₊ * (uᵢ₊ - uᵢ )
    Δu₋ = r₋ * (uᵢ  - uᵢ₋)
    s = sign(Δu₊)
    return s * max(0, min(s * Δu₊, s * Δu₋))
end


@generic_kernel function euler_projection!(
    s::Int, dx::T, dt::T,
    uˢ::V, ρ::V, u::V, v::V, E::V,
    advection_ρ::V, advection_uρ::V, advection_vρ::V, advection_Eρ::V
) where {T, V <: AbstractArray{T}}
    i = @index_2D_lin()

    dX = dx + dt * (uˢ[i+s] - uˢ[i])

    tmp_ρ  = (dX * ρ[i]        - (advection_ρ[i+s]  - advection_ρ[i] )) / dx
    tmp_uρ = (dX * ρ[i] * u[i] - (advection_uρ[i+s] - advection_uρ[i])) / dx
    tmp_vρ = (dX * ρ[i] * v[i] - (advection_vρ[i+s] - advection_vρ[i])) / dx
    tmp_Eρ = (dX * ρ[i] * E[i] - (advection_Eρ[i+s] - advection_Eρ[i])) / dx

    ρ[i] = tmp_ρ
    u[i] = tmp_uρ / tmp_ρ
    v[i] = tmp_vρ / tmp_ρ
    E[i] = tmp_Eρ / tmp_ρ
end


function euler_projection!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    projection_range = block_domain_range(blk.size, state.steps_ranges.projection)
    s = stride_along(blk.size, state.axis)
    blk_data = block_device_data(blk)
    euler_projection!(
        params, blk_data, projection_range, s, state.dx, state.dt,
        blk_data.work_1, blk_data.work_2, blk_data.work_3, blk_data.work_4
    )
end


function euler_projection!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    @section "Projection" @iter_blocks for blk in all_blocks(grid)
        euler_projection!(params, state, blk)
    end
end


@generic_kernel function advection_first_order!(
    s::Int, dt::T,
    uˢ::V, ρ::V, u::V, v::V, E::V,
    advection_ρ::V, advection_uρ::V, advection_vρ::V, advection_Eρ::V
) where {T, V <: AbstractArray{T}}
    i = @index_2D_lin()
    is = i
    disp = dt * uˢ[i]
    if disp > 0
        i = i - s
    end

    advection_ρ[is]  = disp * (ρ[i]       )
    advection_uρ[is] = disp * (ρ[i] * u[i])
    advection_vρ[is] = disp * (ρ[i] * v[i])
    advection_Eρ[is] = disp * (ρ[i] * E[i])
end


function advection_fluxes!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, ::EulerProjection)
    advection_range = block_domain_range(blk.size, state.steps_ranges.advection)
    s = stride_along(blk.size, state.axis)
    blk_data = block_device_data(blk)
    advection_first_order!(
        params, blk_data, advection_range, s, state.dt,
        blk_data.work_1, blk_data.work_2, blk_data.work_3, blk_data.work_4
    )
end


@generic_kernel function advection_second_order!(
    s::Int, dx::T, dt::T,
    uˢ::V, ρ::V, u::V, v::V, E::V,
    advection_ρ::V, advection_uρ::V, advection_vρ::V, advection_Eρ::V
) where {T, V <: AbstractArray{T}}
    i = @index_2D_lin()
    is = i
    disp = dt * uˢ[i]
    if disp > 0
        Δxₑ = -(dx - dt * uˢ[i-s])
        i = i - s
    else
        Δxₑ = dx + dt * uˢ[i+s]
    end

    Δxₗ₋  = dx + dt * (uˢ[i]    - uˢ[i-s])
    Δxₗ   = dx + dt * (uˢ[i+s]  - uˢ[i]  )
    Δxₗ₊  = dx + dt * (uˢ[i+2s] - uˢ[i+s])

    r₋  = (2 * Δxₗ) / (Δxₗ + Δxₗ₋)
    r₊  = (2 * Δxₗ) / (Δxₗ + Δxₗ₊)

    slopes_ρ  = slope_minmod(ρ[i-s]         , ρ[i]       , ρ[i+s]         , r₋, r₊)
    slopes_uρ = slope_minmod(ρ[i-s] * u[i-s], ρ[i] * u[i], ρ[i+s] * u[i+s], r₋, r₊)
    slopes_vρ = slope_minmod(ρ[i-s] * v[i-s], ρ[i] * v[i], ρ[i+s] * v[i+s], r₋, r₊)
    slopes_Eρ = slope_minmod(ρ[i-s] * E[i-s], ρ[i] * E[i], ρ[i+s] * E[i+s], r₋, r₊)

    length_factor = Δxₑ / (2 * Δxₗ)
    advection_ρ[is]  = disp * (ρ[i]        - slopes_ρ  * length_factor)
    advection_uρ[is] = disp * (ρ[i] * u[i] - slopes_uρ * length_factor)
    advection_vρ[is] = disp * (ρ[i] * v[i] - slopes_vρ * length_factor)
    advection_Eρ[is] = disp * (ρ[i] * E[i] - slopes_Eρ * length_factor)
end


function advection_fluxes!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, ::Euler2ndProjection)
    advection_range = block_domain_range(blk.size, state.steps_ranges.advection)
    s = stride_along(blk.size, state.axis)
    blk_data = block_device_data(blk)
    advection_second_order!(
        params, blk_data, advection_range, s, state.dx, state.dt,
        blk_data.work_1, blk_data.work_2, blk_data.work_3, blk_data.work_4
    )
end


advection_fluxes!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock) =
    advection_fluxes!(params::ArmonParameters, state, blk::LocalTaskBlock, state.projection_scheme)

function advection_fluxes!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    @section "Advection" @iter_blocks for blk in all_blocks(grid)
        advection_fluxes!(params, state, blk)
    end
end


function projection_remap!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    advection_fluxes!(params, state, blk)
    euler_projection!(params, state, blk)
end


function projection_remap!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    advection_fluxes!(params, state, grid)
    euler_projection!(params, state, grid)
end
