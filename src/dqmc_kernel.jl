using Random
using LinearAlgebra
import JDQMCFramework as jdqmcf
import JDQMCMeasurements as jdqmcm
import LatticeUtilities as lu

function local_update_G00!(
    G00::AbstractMatrix{T},
    Gτ0::AbstractMatrix{T},
    G0τ::AbstractMatrix{T},
    R,
    Δ,
    i::Int,
    u::AbstractVector{T},
    v::AbstractVector{T}
)::Nothing where {T<:Number}

    copyto!(u, @view G0τ[:, i])

    # vᵀ = G(τ,0)[i, :]
    vt = adjoint(v)
    copyto!(vt, @view Gτ0[i, :])

    # G(0,0) ← G(0,0) + (Δ/R) G(0,τ)[:,i] G(τ,0)[i,:]
    BLAS.ger!(Δ / R, u, v, G00)

    return nothing
end

function local_update_Gτ0!(
    Gτ0::AbstractMatrix{T},
    Gττ::AbstractMatrix{T},
    R,
    Δ,
    i::Int,
    u::AbstractVector{T},
    v::AbstractVector{T}
)::Nothing where {T<:Number}

    copyto!(u, @view Gττ[:, i])

    # vᵀ = G(τ,0)[i, :]
    vt = adjoint(v)
    copyto!(vt, @view Gτ0[i, :])

    BLAS.ger!(Δ / R, u, v, Gτ0)

    return nothing
end

function local_update_G0τ!(
    G0τ::AbstractMatrix{T},
    Gττ::AbstractMatrix{T},
    R,
    Δ,
    i::Int,
    u::AbstractVector{T},
    v::AbstractVector{T}
)::Nothing where {T<:Number}

    copyto!(u, @view G0τ[:, i])

    # vᵀ = G(τ,τ)[i,:] - I[i,:]
    vt = adjoint(v)
    copyto!(vt, @view Gττ[i, :])
    v[i] -= one(T)

    BLAS.ger!(Δ / R, u, v, G0τ)

    return nothing
end

function grover_matrix(G1::AbstractMatrix, G2::AbstractMatrix, C::Vector{Int})
    T = promote_type(eltype(G1), eltype(G2))
    nC = length(C)

    IC = Matrix{T}(I, nC, nC)

    if nC == 0
        return IC
    end

    G1C = Matrix{T}(G1[C, C])
    G2C = Matrix{T}(G2[C, C])

    gC = G1C * G2C + (IC - G1C) * (IC - G2C)

    return gC
end

function logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)
    gCup = grover_matrix(Gup_00, G2up_00, C)
    gCdn = grover_matrix(Gdn_00, G2dn_00, C)

    X = det(gCup) * det(gCdn)
    aX = abs(X)

    return aX == 0 ? -Inf : log(aX)
end

function fused_ratio_from_logs(logX_new, logX_old, λ)
    if λ == 0
        return 1.0
    elseif logX_new == -Inf && logX_old == -Inf
        return 1.0
    elseif logX_new == -Inf
        return 0.0
    elseif logX_old == -Inf
        return Inf
    else
        logR = λ * (logX_new - logX_old)
        return exp(clamp(logR, -700.0, 700.0))
    end
end

function normalize_boundary_condition(boundary_condition)
    normalized =
        boundary_condition isa Symbol ? boundary_condition :
        boundary_condition isa AbstractString ? Symbol(boundary_condition) :
        throw(ArgumentError("boundary_condition must be \"pbc\" or \"apbc\""))

    normalized in (:pbc, :apbc) ||
        throw(ArgumentError("boundary_condition must be \"pbc\" or \"apbc\""))
    return normalized
end

function hubbard_chain_kinetic_matrix(t, L::Integer, boundary_condition = :pbc)
    L >= 2 || throw(ArgumentError("L must be at least 2"))
    boundary_condition = normalize_boundary_condition(boundary_condition)

    K = zeros(typeof(t), L, L)
    for site in 1:(L - 1)
        K[site, site + 1] = -t
        K[site + 1, site] = -conj(t)
    end

    boundary_hopping = boundary_condition == :pbc ? -t : t
    K[L, 1] = boundary_hopping
    K[1, L] = conj(boundary_hopping)
    return K
end

function initializer(
    t,
    U,
    μ,
    β,
    L,
    Δτ,
    seed;
    n_stab = 10,
    boundary_condition = :pbc,
)
    boundary_condition = normalize_boundary_condition(boundary_condition)
    Lτ=jdqmcf.eval_length_imaginary_axis(β,Δτ)
    rng = Xoshiro(seed)
    symmetric = false
    checkerboard = false
    unit_cell = lu.UnitCell(lattice_vecs = [[1.]], basis_vecs = [[0.]])
    lattice = lu.Lattice(L = [L], periodic = [true])

    bond_x = lu.Bond(orbitals = (1,1), displacement = [1])

    neighbor_table = lu.build_neighbor_table([bond_x], unit_cell, lattice)

    N = lu.nsites(unit_cell, lattice)

    N_bonds = size(neighbor_table, 2)

    Δτ′ = symmetric ? Δτ/2 : Δτ

    if checkerboard
        bond_hoppings = fill(t, N_bonds)
        if boundary_condition == :apbc
            seam = findfirst(
                bond -> neighbor_table[1, bond] == L && neighbor_table[2, bond] == 1,
                axes(neighbor_table, 2),
            )
            seam === nothing && error("Periodic seam bond L -> 1 was not found")
            bond_hoppings[seam] = -t
        end
        expnΔτ′K = cb.CheckerboardMatrix(neighbor_table, bond_hoppings, Δτ′)
    else
        K = hubbard_chain_kinetic_matrix(t, L, boundary_condition)
        expnΔτ′K = exp(-Δτ′*K)
        exppΔτ′K = exp(+Δτ′*K)
    end;

    α = acosh(exp(Δτ*U/2))
    s = rand(rng, -1:2:1, N, Lτ)
    T_expnΔτK = eltype(t)
    T_expnΔτV = typeof(α)

    if checkerboard && symmetric

        Bup = jdqmcf.SymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.SymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]

    elseif checkerboard && !symmetric

        Bup = jdqmcf.AsymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.AsymChkbrdPropagator{T_expnΔτK, T_expnΔτV}[]

    elseif !checkerboard && symmetric

        Bup = jdqmcf.SymExactPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.SymExactPropagator{T_expnΔτK, T_expnΔτV}[]

    elseif !checkerboard && !symmetric

        Bup = jdqmcf.AsymExactPropagator{T_expnΔτK, T_expnΔτV}[]
        Bdn = jdqmcf.AsymExactPropagator{T_expnΔτK, T_expnΔτV}[]
    end;

    for l in 1:Lτ

        s_l = @view s[:,l]

        expnΔτVup = zeros(T_expnΔτV, N)
        @. expnΔτVup = exp(+α * s_l + Δτ*μ)

        expnΔτVdn = zeros(T_expnΔτV, N)
        @. expnΔτVdn = exp(-α * s_l + Δτ*μ)

        if checkerboard && symmetric

            push!(Bup, jdqmcf.SymChkbrdPropagator(expnΔτVup, expnΔτ′K))
            push!(Bdn, jdqmcf.SymChkbrdPropagator(expnΔτVdn, expnΔτ′K))

        elseif checkerboard && !symmetric

            push!(Bup, jdqmcf.AsymChkbrdPropagator(expnΔτVup, expnΔτ′K))
            push!(Bdn, jdqmcf.AsymChkbrdPropagator(expnΔτVdn, expnΔτ′K))

        elseif !checkerboard && symmetric

            push!(Bup, jdqmcf.SymExactPropagator(expnΔτVup, expnΔτ′K, exppΔτ′K))
            push!(Bdn, jdqmcf.SymExactPropagator(expnΔτVdn, expnΔτ′K, exppΔτ′K))

        elseif !checkerboard && !symmetric

            push!(Bup, jdqmcf.AsymExactPropagator(expnΔτVup, expnΔτ′K, exppΔτ′K))
            push!(Bdn, jdqmcf.AsymExactPropagator(expnΔτVdn, expnΔτ′K, exppΔτ′K))
        end
    end

    fermion_greens_calc_up = jdqmcf.FermionGreensCalculator(Bup, β, Δτ, n_stab)
    fermion_greens_calc_dn = jdqmcf.FermionGreensCalculator(Bdn, β, Δτ, n_stab);

    Gup = zeros(typeof(t), N, N)
    logdetGup, sgndetGup = jdqmcf.calculate_equaltime_greens!(Gup, fermion_greens_calc_up)

    Gdn = zeros(typeof(t), N, N)
    logdetGdn, sgndetGdn = jdqmcf.calculate_equaltime_greens!(Gdn, fermion_greens_calc_dn);

    Gup_τ0 = zero(Gup)
    Gup_0τ = zero(Gup)
    Gup_00 = copy(Gup)
    Gdn_τ0 = zero(Gdn)
    Gdn_0τ = zero(Gdn)
    Gdn_00 = copy(Gdn);

    jdqmcf.initialize_unequaltime_greens!(Gup_τ0,Gup_0τ,Gup,Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_τ0,Gdn_0τ,Gdn,Gdn_00)

return (
    t = t,
    U = U,
    μ = μ,
    β = β,
    L = L,
    Δτ = Δτ,
    Lτ = Lτ,
    n_stab = n_stab,
    boundary_condition = boundary_condition,
    α=α,

    rng = rng,
    s = s,

    unit_cell = unit_cell,
    lattice = lattice,
    neighbor_table = neighbor_table,

    Bup = Bup,
    Bdn = Bdn,
    fermion_greens_calc_up = fermion_greens_calc_up,
    fermion_greens_calc_dn = fermion_greens_calc_dn,

    Gup = Gup,
    logdetGup = logdetGup,
    sgndetGup = sgndetGup,
    Gup_τ0 = Gup_τ0,
    Gup_0τ = Gup_0τ,
    Gup_00 = Gup_00,

    Gdn = Gdn,
    logdetGdn = logdetGdn,
    sgndetGdn = sgndetGdn,
    Gdn_τ0 = Gdn_τ0,
    Gdn_0τ = Gdn_0τ,
    Gdn_00 = Gdn_00
)

end

function _local_update_replica!(
    Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup,
    Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn,
    s, μ, α, Δτ, rng, perm, u, v
)

    shuffle!(rng, perm)

    accepted = 0
    N = length(perm)

    for i in perm
        ΔVup = -α / Δτ * (-2 * s[i])
        ΔVdn = +α / Δτ * (-2 * s[i])

        Rup, Δup = jdqmcf.local_update_det_ratio(Gup, ΔVup, i, Δτ)
        Rdn, Δdn = jdqmcf.local_update_det_ratio(Gdn, ΔVdn, i, Δτ)

        P = min(1.0, abs(Rup * Rdn))

        if rand(rng) < P
            accepted += 1
            s[i] = -s[i]

            local_update_G00!(Gup_00, Gup_τ0, Gup_0τ, Rup, Δup, i, u, v)
            local_update_Gτ0!(Gup_τ0, Gup, Rup, Δup, i, u, v)
            local_update_G0τ!(Gup_0τ, Gup, Rup, Δup, i, u, v)

            logdetGup, sgndetGup = jdqmcf.local_update_greens!(
                Gup, logdetGup, sgndetGup, Bup, Rup, Δup, i, u, v
            )

            local_update_G00!(Gdn_00, Gdn_τ0, Gdn_0τ, Rdn, Δdn, i, u, v)
            local_update_Gτ0!(Gdn_τ0, Gdn, Rdn, Δdn, i, u, v)
            local_update_G0τ!(Gdn_0τ, Gdn, Rdn, Δdn, i, u, v)

            logdetGdn, sgndetGdn = jdqmcf.local_update_greens!(
                Gdn, logdetGdn, sgndetGdn, Bdn, Rdn, Δdn, i, u, v
            )
        end
    end

    acceptance_rate = accepted / N

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, acceptance_rate
end

function _local_update_replica_fused!(
    Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup,
    Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn,
    G2up_00, G2dn_00,
    s, μ, α, Δτ, rng, perm, u, v,
    C, λ
)

    shuffle!(rng, perm)

    accepted = 0
    N = length(perm)

    # Cache the current Grover log weight for this time slice.
    # It changes only when the tentative G00 update is accepted.
    logX_current = logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)

    for i in perm
        ΔVup = -α / Δτ * (-2 * s[i])
        ΔVdn = +α / Δτ * (-2 * s[i])

        Rup, Δup = jdqmcf.local_update_det_ratio(Gup, ΔVup, i, Δτ)
        Rdn, Δdn = jdqmcf.local_update_det_ratio(Gdn, ΔVdn, i, Δτ)

        Gup_00_old = copy(Gup_00)
        Gdn_00_old = copy(Gdn_00)

        logX_old = logX_current

        # Tentatively update G00 only.
        local_update_G00!(Gup_00, Gup_τ0, Gup_0τ, Rup, Δup, i, u, v)
        local_update_G00!(Gdn_00, Gdn_τ0, Gdn_0τ, Rdn, Δdn, i, u, v)

        logX_new = logabsX_from_G00(Gup_00, Gdn_00, G2up_00, G2dn_00, C)

        Rλ = fused_ratio_from_logs(logX_new, logX_old, λ)

        P = min(1.0, abs(Rup * Rdn) * Rλ)

        if rand(rng) < P
            accepted += 1
            s[i] = -s[i]

            local_update_Gτ0!(Gup_τ0, Gup, Rup, Δup, i, u, v)
            local_update_G0τ!(Gup_0τ, Gup, Rup, Δup, i, u, v)

            logdetGup, sgndetGup = jdqmcf.local_update_greens!(
                Gup, logdetGup, sgndetGup, Bup, Rup, Δup, i, u, v
            )

            local_update_Gτ0!(Gdn_τ0, Gdn, Rdn, Δdn, i, u, v)
            local_update_G0τ!(Gdn_0τ, Gdn, Rdn, Δdn, i, u, v)

            logdetGdn, sgndetGdn = jdqmcf.local_update_greens!(
                Gdn, logdetGdn, sgndetGdn, Bdn, Rdn, Δdn, i, u, v
            )

            logX_current = logX_new
        else
            copyto!(Gup_00, Gup_00_old)
            copyto!(Gdn_00, Gdn_00_old)
        end
    end

    acceptance_rate = accepted / N

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, acceptance_rate
end

function local_update_replica!(
    Gup::Matrix{T}, logdetGup, sgndetGup, Bup, fermion_greens_calc_up,
    Gdn::Matrix{T}, logdetGdn, sgndetGdn, Bdn, fermion_greens_calc_dn,
    Gup_00, Gup_τ0, Gup_0τ,
    Gdn_00, Gdn_τ0, Gdn_0τ,
    s, μ, α, Δτ, δG, rng
) where {T<:Number}

    Lτ = length(Bup)
    N = size(Gup, 1)

    u = zeros(T, N)
    v = zeros(T, N)
    perm = collect(1:N)

    acceptance_rate = 0.0

    jdqmcf.initialize_unequaltime_greens!(Gup_τ0, Gup_0τ, Gup, Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_τ0, Gdn_0τ, Gdn, Gdn_00)

    for l in fermion_greens_calc_up
        jdqmcf.propagate_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup,
            fermion_greens_calc_up, Bup
        )

        jdqmcf.propagate_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn,
            fermion_greens_calc_dn, Bdn
        )

        s_l = @view s[:, l]

        logdetGup, sgndetGup,
        logdetGdn, sgndetGdn,
        acceptance_rate_l = _local_update_replica!(
            Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup[l],
            Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn[l],
            s_l, μ, α, Δτ, rng, perm, u, v
        )

        acceptance_rate += acceptance_rate_l / Lτ

        logdetGup, sgndetGup, δGup, δθup = jdqmcf.stabilize_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup, Gup_00,
            logdetGup, sgndetGup,
            fermion_greens_calc_up, Bup,
            update_B̄ = true
        )

        logdetGdn, sgndetGdn, δGdn, δθdn = jdqmcf.stabilize_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn, Gdn_00,
            logdetGdn, sgndetGdn,
            fermion_greens_calc_dn, Bdn,
            update_B̄ = true
        )

        δG = max(δG, δGup, δGdn)

        iterate(fermion_greens_calc_dn, fermion_greens_calc_up)
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, acceptance_rate
end

function local_update_replica_fused!(
    Gup::Matrix{T}, logdetGup, sgndetGup, Bup, fermion_greens_calc_up,
    Gdn::Matrix{T}, logdetGdn, sgndetGdn, Bdn, fermion_greens_calc_dn,
    Gup_00, Gup_τ0, Gup_0τ,
    Gdn_00, Gdn_τ0, Gdn_0τ,
    G2up_00, G2dn_00,
    s, μ, α, Δτ, δG, rng,
    C, λ
) where {T<:Number}

    Lτ = length(Bup)
    N = size(Gup, 1)

    u = zeros(T, N)
    v = zeros(T, N)
    perm = collect(1:N)

    acceptance_rate = 0.0

    jdqmcf.initialize_unequaltime_greens!(Gup_τ0, Gup_0τ, Gup, Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_τ0, Gdn_0τ, Gdn, Gdn_00)

    for l in fermion_greens_calc_up
        jdqmcf.propagate_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup,
            fermion_greens_calc_up, Bup
        )

        jdqmcf.propagate_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn,
            fermion_greens_calc_dn, Bdn
        )

        s_l = @view s[:, l]

        logdetGup, sgndetGup,
        logdetGdn, sgndetGdn,
        acceptance_rate_l = _local_update_replica_fused!(
            Gup, Gup_00, Gup_τ0, Gup_0τ, logdetGup, sgndetGup, Bup[l],
            Gdn, Gdn_00, Gdn_τ0, Gdn_0τ, logdetGdn, sgndetGdn, Bdn[l],
            G2up_00, G2dn_00,
            s_l, μ, α, Δτ, rng, perm, u, v,
            C, λ
        )

        acceptance_rate += acceptance_rate_l / Lτ

        logdetGup, sgndetGup, δGup, δθup = jdqmcf.stabilize_unequaltime_greens!(
            Gup_τ0, Gup_0τ, Gup, Gup_00,
            logdetGup, sgndetGup,
            fermion_greens_calc_up, Bup,
            update_B̄ = true
        )

        logdetGdn, sgndetGdn, δGdn, δθdn = jdqmcf.stabilize_unequaltime_greens!(
            Gdn_τ0, Gdn_0τ, Gdn, Gdn_00,
            logdetGdn, sgndetGdn,
            fermion_greens_calc_dn, Bdn,
            update_B̄ = true
        )

        δG = max(δG, δGup, δGdn)

        iterate(fermion_greens_calc_dn, fermion_greens_calc_up)
    end

    return logdetGup, sgndetGup, logdetGdn, sgndetGdn, δG, acceptance_rate
end

function update_replica!(replica, δG)
    logdetGup, sgndetGup,
    logdetGdn, sgndetGdn,
    δG, acceptance_rate = local_update_replica!(
        replica.Gup, replica.logdetGup, replica.sgndetGup,
        replica.Bup, replica.fermion_greens_calc_up,

        replica.Gdn, replica.logdetGdn, replica.sgndetGdn,
        replica.Bdn, replica.fermion_greens_calc_dn,

        replica.Gup_00, replica.Gup_τ0, replica.Gup_0τ,
        replica.Gdn_00, replica.Gdn_τ0, replica.Gdn_0τ,

        replica.s, replica.μ, replica.α, replica.Δτ,
        δG, replica.rng
    )

    replica = merge(replica, (
        logdetGup = logdetGup,
        sgndetGup = sgndetGup,
        logdetGdn = logdetGdn,
        sgndetGdn = sgndetGdn,
    ))

    return replica, δG, acceptance_rate
end

function update_replica_fused!(active, passive, δG, C, λ)

    # Exact fallback to the already-validated independent estimator.
    if λ == 0
        return update_replica!(active, δG)
    end

    logdetGup, sgndetGup,
    logdetGdn, sgndetGdn,
    δG, acceptance_rate = local_update_replica_fused!(
        active.Gup, active.logdetGup, active.sgndetGup,
        active.Bup, active.fermion_greens_calc_up,

        active.Gdn, active.logdetGdn, active.sgndetGdn,
        active.Bdn, active.fermion_greens_calc_dn,

        active.Gup_00, active.Gup_τ0, active.Gup_0τ,
        active.Gdn_00, active.Gdn_τ0, active.Gdn_0τ,

        passive.Gup_00, passive.Gdn_00,

        active.s, active.μ, active.α, active.Δτ,
        δG, active.rng,

        C, λ
    )

    active = merge(active, (
        logdetGup = logdetGup,
        sgndetGup = sgndetGup,
        logdetGdn = logdetGdn,
        sgndetGdn = sgndetGdn,
    ))

    return active, δG, acceptance_rate
end

function power_increment_statistics(
    avg_phase,
    avg_den,
    avg_num
)
    phase, phase_err =
        jdqmcm.jackknife(
            identity,
            avg_phase
        )

    den, den_err =
        jdqmcm.jackknife(
            identity,
            avg_den
        )

    num, num_err =
        jdqmcm.jackknife(
            identity,
            avg_num
        )

    # Z₁ / Z_λ
    Z1_over_Zlambda, Z1_over_Zlambda_err =
        jdqmcm.jackknife(
            (num_bin, phase_bin) ->
                real(num_bin / phase_bin),
            avg_num,
            avg_phase
        )

    # Z_λ / Z₀
    Zlambda_over_Z0, Zlambda_over_Z0_err =
        jdqmcm.jackknife(
            (phase_bin, den_bin) ->
                real(phase_bin / den_bin),
            avg_phase,
            avg_den
        )

    # Z₁ / Z₀
    Z1_over_Z0, Z1_over_Z0_err =
        jdqmcm.jackknife(
            (num_bin, den_bin) ->
                real(num_bin / den_bin),
            avg_num,
            avg_den
        )

    # S₂ = -log|Z₁/Z₀|
    S2, S2_err =
        jdqmcm.jackknife(
            (num_bin, den_bin) ->
                -log(abs(num_bin / den_bin)),
            avg_num,
            avg_den
        )

    return (
        phase = real(phase),
        phase_err = phase_err,
        den = real(den),
        den_err = den_err,
        num = real(num),
        num_err = num_err,
        Z1_over_Zlambda = real(Z1_over_Zlambda),
        Z1_over_Zlambda_err = Z1_over_Zlambda_err,
        Zlambda_over_Z0 = real(Zlambda_over_Z0),
        Zlambda_over_Z0_err = Zlambda_over_Z0_err,
        Z1_over_Z0 = real(Z1_over_Z0),
        Z1_over_Z0_err = Z1_over_Z0_err,
        S2 = real(S2),
        S2_err = S2_err
    )
end
