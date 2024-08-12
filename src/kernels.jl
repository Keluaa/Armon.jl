
# TODO: make the stride deductible at compile-time (i.e. one kernel instanciation per axis)

@generic_kernel function perfect_gas_EOS!(
    γ::T, ρ::V, E::V, p::V, c::V, g::V, U::NTuple{D, V}
) where {T, V <: AbstractArray{T}, D}
    i = @kt_i()
    e = E[i] - 0.5 * sum(get_tuple(U, i).^2)
    p[i] = (γ - 1.) * ρ[i] * e
    c[i] = sqrt(γ * p[i] / ρ[i])
    g[i] = (1. + γ) / 2
end


@generic_kernel function bizarrium_EOS!(
    ρ::V, E::V, p::V, c::V, g::V, U::NTuple{D, V}
) where {T, V <: AbstractArray{T}, D}
    i = @kt_i()

    # O. Heuzé, S. Jaouen, H. Jourdren, 
    # "Dissipative issue of high-order shock capturing schemes with non-convex equations of state"
    # JCP 2009

    # TODO: @kernel_init begin
        rho0::T = 10000.
        K0::T   = 1e+11
        Cv0::T  = 1000.
        T0::T   = 300.
        eps0::T = 0.
        G0::T   = 1.5
        s::T    = 1.5
        q::T    = -42080895/14941154
        r::T    = 727668333/149411540
    # end

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

    e = E[i] - 0.5 * sum(get_tuple(U, i).^2)
    p[i] = pk0 + G0 * rho0 * (e - epsk0)
    c[i] = sqrt(G0 * rho0 * (p[i] - pk0) - pk0prime) / ρ[i]
    g[i] = 0.5 / (ρ[i]^3 * c[i]^2) * (pk0second + (G0 * rho0)^2 * (p[i] - pk0))
end


@generic_kernel function cell_update!(
    s::Int, dx::T, dt::T,
    uˢ::V, pˢ::V, ρ::V, uₐ::V, E::V
) where {T, V <: AbstractArray{T}}
    i = @kt_i()
    u = uₐ  # `u` or `v` depending on the current axis
    dm = ρ[i] * dx
    ρ[i]  = dm / (dx + dt * (uˢ[i+s] - uˢ[i]))
    u[i] += dt / dm * (pˢ[i]         - pˢ[i+s]          )
    E[i] += dt / dm * (pˢ[i] * uˢ[i] - pˢ[i+s] * uˢ[i+s])
end


@kernel_function function init_vars(
    test_case::TwoStateTestCase, test_params::InitTestParamsTwoState, X::NTuple{D},
    i, ρ::V, E::V, ::NTuple{D, V}, p::V, c::V, g::V
) where {V, D}
    if test_region_high(X, test_case)
        ρ[i] = test_params.high_ρ
        E[i] = test_params.high_E
        set_tuple!(U, test_params.high_u, i)
    else
        ρ[i] = test_params.low_ρ
        E[i] = test_params.low_E
        set_tuple!(U, test_params.low_u, i)
    end

    p[i] = zero(eltype(V))
    c[i] = zero(eltype(V))
    g[i] = zero(eltype(V))
end


@kernel_function function init_vars(
    ::DebugIndexes, i, global_i, ρ::V, E::V, U::NTuple{D, V}, p::V, c::V, g::V
) where {V, D}
    ρ[i] = global_i
    E[i] = global_i
    set_tuple!(U, global_i, i)
    p[i] = global_i
    c[i] = global_i
    g[i] = global_i
end


@generic_kernel function init_test(
    global_pos::NTuple{D, Int}, N::NTuple{D, Int}, bsize::BSize,
    origin::NTuple{D, T}, ΔX::NTuple{D, T},
    X::NTuple{D, V}, mask::V, ρ::V, E::V, U::NTuple{D, V}, p::V, c::V, g::V, vars_to_zero::Tuple{Vararg{V}},
    test_case::Test
) where {T, V <: AbstractArray{T}, D, Test <: TestCase, BSize <: BlockSize{D}}
    #= TODO: @kernel_init =# begin
        if Test <: TwoStateTestCase
            test_init_params = init_test_params(test_case, T, D)
        end
    end

    i = @kt_i()
    # TODO: use `@kt_I` + an offset instead
    I = position(bsize, i)  # Position in the block's real cells

    # Index in the global grid (0-indexed)
    gI = I .+ global_pos .- 1

    # Position in the global grid
    pos = set_tuple!(X, gI .* ΔX .- origin, i)

    # Set the domain mask to 1 if the cell is real or 0 otherwise
    mask[i] = is_ghost(bsize, I) ? 0 : 1

    # Middle point of the cell
    mid = pos .+ ΔX ./ 2

    if Test <: TwoStateTestCase
        init_vars(test_case, test_init_params, mid, i, ρ, E, U, p, c, g)
    elseif Test <: DebugIndexes
        global_i = sum(gI .* Base.size_to_strides(1, N...)) + 1
        init_vars(test_case, i, global_i, ρ, E, U, p, c, g)
    else
        init_vars(test_case, mid, i, ρ, E, U, p, c, g)
    end

    for var in vars_to_zero
        var[i] = zero(T)
    end
end

#
# Wrappers
#

function update_EOS!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, tc::TestCase)
    domain = block_domain_range(blk.size, state.steps_ranges.EOS)
    γ = eltype(blk)(specific_heat_ratio(tc))
    data = block_data(blk)
    u = data.dim_vars.u
    (; ρ, E, p, c, g) = data.scalar_vars
    return perfect_gas_EOS!(γ, ρ, E, p, c, g, u; ctx=params.kernel_ctx, domain)
end


function update_EOS!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock, ::Bizarrium)
    domain = block_domain_range(blk.size, state.steps_ranges.EOS)
    data = block_data(blk)
    u = data.dim_vars.u
    (; ρ, E, p, c, g) = data.scalar_vars
    return bizarrium_EOS!(ρ, E, p, c, g, u; ctx=params.kernel_ctx, domain)
end


function update_EOS!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    return update_EOS!(params, state, blk, state.test_case)
end


function update_EOS!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    @iter_blocks for blk in grid
        update_EOS!(params, state, blk)
    end
end


function init_test(params::ArmonParameters, blk::LocalTaskBlock)
    domain = block_domain_range(blk.size, first(params.steps_ranges).full_domain)

    # Position of the origin of this block
    real_static_bsize = params.block_size .- 2*params.nghost
    blk_global_pos = params.N_origin .- 1 .+ (Tuple(blk.pos) .- 1) .* real_static_bsize

    # Cell dimensions
    ΔX = params.domain_size ./ params.global_grid

    # Make sure all variables are initialized. This is also to make sure to touch every memory page
    # of the block's variables, ensuring that they are stored in the right NUMA group.
    vars_names_to_zero = setdiff(block_vars(), (:x, :ρ, :E, :u, :p, :c, :g, :mask))
    vars_to_zero = var_arrays(blk, vars_names_to_zero)

    data = block_data(blk)
    (; u) = data.dim_vars
    (; X, mask, ρ, E, p, c, g) = data.scalar_vars
    init_test(
        blk_global_pos, params.N, blk.size, params.origin, ΔX,
        X, mask, ρ, E, u, p, c, g, vars_to_zero, params.test;
        ctx=params.kernel_ctx, domain
    )

    if params.numa_aware
        # Now is the best time to move the pages, as we know they exist physically and that they are
        # currently cannot be used by other threads.
        target_numa_node = NUMA.current_numa_node()
        move_pages(blk, target_numa_node)
        params.lock_memory && lock_pages(blk)

        # Do the exact same with the MPI buffers associated with the block
        for neighbour in blk.neighbours
            !(neighbour isa RemoteTaskBlock) && continue
            move_pages(neighbour, target_numa_node)
            params.lock_memory && lock_pages(neighbour)
        end
    end
end


function init_test(params::ArmonParameters, grid::BlockGrid)
    @iter_blocks for blk in grid
        init_test(params, blk)
    end
end


function cell_update!(params::ArmonParameters, state::SolverState, blk::LocalTaskBlock)
    domain = block_domain_range(blk.size, state.steps_ranges.cell_update)
    s = stride_along(blk.size, state.axis)
    data = block_device_data(blk)
    uₐ = data.dim_vars.u[state.axis]
    (; ρ, E, uˢ, pˢ) = data.scalar_vars
    cell_update!(s, state.dx, state.dt, uˢ, pˢ, ρ, uₐ, E; ctx=params.kernel_ctx, domain)
end


function cell_update!(params::ArmonParameters, state::SolverState, grid::BlockGrid)
    @iter_blocks for blk in grid
        cell_update!(params, state, blk)
    end
end
