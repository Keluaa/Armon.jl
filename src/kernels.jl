
# TODO: make the stride deductible at compile-time (i.e. one kernel instanciation per axis)

@generic_kernel function perfect_gas_EOS!(
    γ::T,
    ρ::V, E::V, u::V, v::V, p::V, c::V, g::V
) where {T, V <: AbstractArray{T}}
    i = @index_2D_lin()
    e = E[i] - 0.5 * (u[i]^2 + v[i]^2)
    p[i] = (γ - 1.) * ρ[i] * e
    c[i] = sqrt(γ * p[i] / ρ[i])
    g[i] = (1. + γ) / 2
end


@generic_kernel function bizarrium_EOS!(
    ρ::V, u::V, v::V, E::V, p::V, c::V, g::V
) where {T, V <: AbstractArray{T}}
    i = @index_2D_lin()

    # O. Heuzé, S. Jaouen, H. Jourdren, 
    # "Dissipative issue of high-order shock capturing schemes with non-convex equations of state"
    # JCP 2009

    @kernel_init begin
        rho0::T = 10000.
        K0::T   = 1e+11
        Cv0::T  = 1000.
        T0::T   = 300.
        eps0::T = 0.
        G0::T   = 1.5
        s::T    = 1.5
        q::T    = -42080895/14941154
        r::T    = 727668333/149411540
    end

    x = ρ[i] / rho0 - 1
    G = G0 * (1-rho0 / ρ[i])

    f0 = (1+(s/3-2)*x+q*x^2+r*x^3)/(1-s*x)
    f1 = (s/3-2+2*q*x+3*r*x^2+s*f0)/(1-s*x)
    f2 = (2*q+6*r*x+2*s*f1)/(1-s*x)
    f3 = (6*r+3*s*f2)/(1-s*x)

    epsk0     = eps0 - Cv0*T0*(1+G) + 0.5*(K0/rho0)*x^2*f0
    pk0       = -Cv0*T0*G0*rho0 + 0.5*K0*x*(1+x)^2*(2*f0+x*f1)
    pk0prime  = -0.5*K0*(1+x)^3*rho0 * (2*(1+3x)*f0 + 2*x*(2+3x)*f1 + x^2*(1+x)*f2)
    pk0second = 0.5*K0*(1+x)^4*rho0^2 * (12*(1+2x)*f0 + 6*(1+6x+6*x^2)*f1 + 
                                                    6*x*(1+x)*(1+2x)*f2 + x^2*(1+x)^2*f3)

    e = E[i] - 0.5 * (u[i]^2 + v[i]^2)
    p[i] = pk0 + G0 * rho0 * (e - epsk0)
    c[i] = sqrt(G0 * rho0 * (p[i] - pk0) - pk0prime) / ρ[i]
    g[i] = 0.5 / (ρ[i]^3 * c[i]^2) * (pk0second + (G0 * rho0)^2 * (p[i] - pk0))
end


@generic_kernel function cell_update!(
    s::Int, dx::T, dt::T, 
    uˢ::V, pˢ::V, ρ::V, uₐ::V, E::V
) where {T, V <: AbstractArray{T}}
    i = @index_2D_lin()
    u = uₐ  # `u` or `v` depending on the current axis
    dm = ρ[i] * dx
    ρ[i]  = dm / (dx + dt * (uˢ[i+s] - uˢ[i]))
    u[i] += dt / dm * (pˢ[i]         - pˢ[i+s]          )
    E[i] += dt / dm * (pˢ[i] * uˢ[i] - pˢ[i+s] * uˢ[i+s])
end


@kernel_function function init_vars(
    test_case::TwoStateTestCase, test_params::InitTestParamsTwoState, X::NTuple,
    i, ρ::V, E::V, u::V, v::V, p::V, c::V, g::V
) where {V}
    if test_region_high(X, test_case)
        ρ[i] = test_params.high_ρ
        E[i] = test_params.high_E
        u[i] = test_params.high_u
        v[i] = test_params.high_v
    else
        ρ[i] = test_params.low_ρ
        E[i] = test_params.low_E
        u[i] = test_params.low_u
        v[i] = test_params.low_v
    end

    p[i] = zero(eltype(V))
    c[i] = zero(eltype(V))
    g[i] = zero(eltype(V))
end


@kernel_function function init_vars(
    ::DebugIndexes, i, global_i, ρ::V, E::V, u::V, v::V, p::V, c::V, g::V
) where {V}
    ρ[i] = global_i
    E[i] = global_i
    u[i] = global_i
    v[i] = global_i
    p[i] = global_i
    c[i] = global_i
    g[i] = global_i
end


@generic_kernel function init_test(
    global_pos::NTuple{2, Int}, N::NTuple{2, Int}, bsize::BSize,
    origin::NTuple{2, T}, ΔX::NTuple{2, T},
    x::V, y::V, mask::V, ρ::V, E::V, u::V, v::V, p::V, c::V, g::V, vars_to_zero::Tuple{Vararg{V}},
    test_case::Test
) where {T, V <: AbstractArray{T}, Test <: TestCase, BSize <: BlockSize}
    @kernel_init begin
        if Test <: TwoStateTestCase
            test_init_params = init_test_params(test_case, T)
        end
    end

    i = @index_2D_lin()
    I = position(bsize, i)  # Position in the block's real cells

    # Index in the global grid (0-indexed)
    gI = I .+ global_pos .- 1

    # Position in the global grid
    (x[i], y[i]) = gI .* ΔX .+ origin

    # Set the domain mask to 1 if the cell is real or 0 otherwise
    mask[i] = is_ghost(bsize, i) ? 0 : 1

    # Middle point of the cell
    mid = (x[i], y[i]) .+ ΔX ./ 2

    if Test <: TwoStateTestCase
        init_vars(test_case, test_init_params, mid, i, ρ, E, u, v, p, c, g)
    elseif Test <: DebugIndexes
        global_i = sum(gI .* Base.size_to_strides(1, N...)) + 1
        init_vars(test_case, i, global_i, ρ, E, u, v, p, c, g)
    else
        init_vars(test_case, mid, i, ρ, E, u, v, p, c, g)
    end

    for var in vars_to_zero
        var[i] = zero(T)
    end
end

#
# Wrappers
#

function update_EOS!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, tc::TestCase)
    range = block_domain_range(blk.size, state.steps_ranges.EOS)
    gamma = eltype(blk)(specific_heat_ratio(tc))
    return perfect_gas_EOS!(params, block_device_data(blk), range, gamma)
end


function update_EOS!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, ::Bizarrium)
    range = block_domain_range(blk.size, state.steps_ranges.EOS)
    return bizarrium_EOS!(params, block_device_data(blk), range)
end


function update_EOS!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    return update_EOS!(params, state, blk, state.test_case)
end


function update_EOS!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    @iter_blocks for blk in all_blocks(grid)
        update_EOS!(params, state, blk)
    end
end


function init_test(params::ArmonParameters, blk::LocalTaskBlock)
    blk_domain = block_domain_range(blk.size, first(params.steps_ranges).full_domain)

    # Position of the origin of this block
    real_static_bsize = params.block_size .- 2*params.nghost
    blk_global_pos = params.cart_coords .* params.N .+ (Tuple(blk.pos) .- 1) .* real_static_bsize

    # Cell dimensions
    ΔX = params.domain_size ./ params.global_grid

    # Make sure all variables to save are initialized
    vars_names_to_zero = setdiff(saved_vars(), (:x, :y, :ρ, :E, :u, :v, :p, :c, :g, :mask))
    vars_to_zero = Tuple(get_vars(blk, vars_names_to_zero))

    init_test(params, block_device_data(blk), blk_domain, blk_global_pos, blk.size, ΔX, vars_to_zero, params.test)
end


function init_test(params::ArmonParameters, grid::BlockGrid)
    @iter_blocks for blk in all_blocks(grid)
        init_test(params, blk)
    end
end


function cell_update!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    blk_domain = block_domain_range(blk.size, state.steps_ranges.cell_update)
    blk_data = block_device_data(blk)
    u = state.axis == Axis.X ? blk_data.u : blk_data.v
    s = stride_along(blk.size, state.axis)
    cell_update!(params, blk_data, blk_domain, s, state.dx, state.dt, u)
end


function cell_update!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    @iter_blocks for blk in all_blocks(grid)
        cell_update!(params, state, blk)
    end
end
