"""
    square_lattice_kinetic_matrix(Lx, Ly; t=1.0, periodic=(true, true))

Nearest-neighbor hopping matrix for a rectangular square lattice. Site `(x,y)`
is `x + (y-1)*Lx`. Each undirected bond occurs once, including when a periodic
axis has length two; an axis of length one adds no self hopping.
The Hamiltonian contains `sum(K[i,j] * c†[i] * c[j])` for each spin.
"""
function square_lattice_kinetic_matrix(Lx::Integer, Ly::Integer;
    t::Real = 1.0, periodic = (true, true))
    Lx > 0 && Ly > 0 || throw(ArgumentError("Lx and Ly must be positive"))
    length(periodic) == 2 && all(x -> x isa Bool, periodic) ||
        throw(ArgumentError("periodic must contain two Booleans"))
    isfinite(t) || throw(ArgumentError("t must be finite"))
    unit_cell = lu.UnitCell(lattice_vecs = [[1.0, 0.0], [0.0, 1.0]],
        basis_vecs = [[0.0, 0.0]])
    lattice = lu.Lattice(L = [Lx, Ly], periodic = collect(periodic))
    bonds = [lu.Bond(orbitals = (1, 1), displacement = [1, 0]),
             lu.Bond(orbitals = (1, 1), displacement = [0, 1])]
    neighbors = lu.build_neighbor_table(bonds, unit_cell, lattice)
    K = zeros(Float64, Lx * Ly, Lx * Ly)
    for bond in axes(neighbors, 2)
        i, j = neighbors[:, bond]
        i == j && continue
        K[i, j] = K[j, i] = -Float64(t)
    end
    return K
end

"""
    initializer_from_kinetic(K, U, mu, beta, dtau, seed; n_stab=10)

Initialize the existing stabilized JDQMCFramework kernel on a real symmetric
hopping matrix `K`. The interaction is `U*(nup-1/2)*(ndn-1/2)` and the chemical
potential term is `-mu*(nup+ndn)`. Repulsive `U >= 0` is supported.
"""
function initializer_from_kinetic(K::AbstractMatrix{<:Real}, U::Real, mu::Real,
    beta::Real, dtau::Real, seed::Integer; n_stab::Integer = 10)
    size(K, 1) == size(K, 2) > 0 || throw(ArgumentError("K must be square"))
    issymmetric(K) && all(isfinite, K) ||
        throw(ArgumentError("K must be finite and real symmetric"))
    all(isfinite, (U, mu, beta, dtau)) || throw(ArgumentError("parameters must be finite"))
    U, mu, beta, dtau = Float64.((U, mu, beta, dtau))
    n_stab = Int(n_stab)
    U >= 0 || throw(ArgumentError("the spin Hubbard-Stratonovich transform requires U >= 0"))
    beta > 0 && dtau > 0 || throw(ArgumentError("beta and dtau must be positive"))
    Ltau = round(Int, beta / dtau)
    Ltau > 0 && isapprox(Ltau * dtau, beta; rtol = 1e-10, atol = 1e-12) ||
        throw(ArgumentError("beta must be an integer multiple of dtau"))
    1 <= n_stab <= Ltau || throw(ArgumentError("n_stab must be between 1 and beta/dtau"))
    N = size(K, 1)
    rng = Xoshiro(seed)
    alpha = acosh(exp(dtau * U / 2))
    s = rand(rng, -1:2:1, N, Ltau)
    K = Matrix{Float64}(K)
    expminusK, expplusK = exp(-dtau * K), exp(dtau * K)
    Bup = jdqmcf.AsymExactPropagator{Float64, Float64}[]
    Bdn = jdqmcf.AsymExactPropagator{Float64, Float64}[]
    for l in 1:Ltau
        potential_up = exp.(alpha .* s[:, l] .+ dtau * mu)
        potential_dn = exp.(-alpha .* s[:, l] .+ dtau * mu)
        push!(Bup, jdqmcf.AsymExactPropagator(potential_up, expminusK, expplusK))
        push!(Bdn, jdqmcf.AsymExactPropagator(potential_dn, expminusK, expplusK))
    end
    calc_up = jdqmcf.FermionGreensCalculator(Bup, beta, dtau, n_stab)
    calc_dn = jdqmcf.FermionGreensCalculator(Bdn, beta, dtau, n_stab)
    Gup, Gdn = zeros(N, N), zeros(N, N)
    logdetGup, sgndetGup = jdqmcf.calculate_equaltime_greens!(Gup, calc_up)
    logdetGdn, sgndetGdn = jdqmcf.calculate_equaltime_greens!(Gdn, calc_dn)
    Gup_00, Gdn_00 = copy(Gup), copy(Gdn)
    Gup_tau0, Gup_0tau, Gdn_tau0, Gdn_0tau = zero(Gup), zero(Gup), zero(Gdn), zero(Gdn)
    jdqmcf.initialize_unequaltime_greens!(Gup_tau0, Gup_0tau, Gup, Gup_00)
    jdqmcf.initialize_unequaltime_greens!(Gdn_tau0, Gdn_0tau, Gdn, Gdn_00)
    return (; K, t = 1.0, U = Float64(U), μ = Float64(mu), β = Float64(beta),
        L = N, Δτ = Float64(dtau), Lτ = Ltau, n_stab,
        boundary_condition = :custom, α = alpha, rng, s, Bup, Bdn,
        fermion_greens_calc_up = calc_up, fermion_greens_calc_dn = calc_dn,
        Gup, logdetGup, sgndetGup, Gup_00, Gup_τ0 = Gup_tau0, Gup_0τ = Gup_0tau,
        Gdn, logdetGdn, sgndetGdn, Gdn_00, Gdn_τ0 = Gdn_tau0, Gdn_0τ = Gdn_0tau)
end

function dqmc_observables(state)
    N = size(state.Gup, 1)
    nup = 1 .- diag(state.Gup)
    ndn = 1 .- diag(state.Gdn)
    density = sum(nup + ndn) / N
    double_occupancy = sum(nup .* ndn) / N
    kinetic = sum(state.K .* (2I - transpose(state.Gup) - transpose(state.Gdn)))
    interaction = state.U * sum(nup .* ndn - (nup + ndn) / 2 .+ 0.25)
    energy = kinetic + interaction - state.μ * N * density
    return (; energy = real(energy), density = real(density),
        double_occupancy = real(double_occupancy))
end

"""
    run_square_dqmc(; Lx=4, Ly=4, t=1, U=4, mu=0, beta=4, dtau=0.1, ...)

Run finite-temperature, grand-canonical Hubbard DQMC with the established
stabilized local-update kernel. Return sign-reweighted energy, density and
double occupancy, their bin jackknife errors, and the bins. `root_seed` is the
only user seed. Errors measure sampling uncertainty, not time-step bias.
"""
function run_square_dqmc(; Lx::Integer = 4, Ly::Integer = 4, t::Real = 1.0,
    U::Real = 4.0, mu::Real = 0.0, beta::Real = 4.0, dtau::Real = 0.1,
    periodic = (true, true), n_stab::Integer = 10, root_seed::Integer = 1234,
    burnin::Integer = 200, measurements::Integer = 2000, bins::Integer = 20,
    sweeps_per_measurement::Integer = 1)
    burnin >= 0 || throw(ArgumentError("burnin must be nonnegative"))
    2 <= bins <= measurements && measurements % bins == 0 ||
        throw(ArgumentError("measurements must divide evenly into at least two bins"))
    sweeps_per_measurement > 0 || throw(ArgumentError("sweeps_per_measurement must be positive"))
    K = square_lattice_kinetic_matrix(Lx, Ly; t, periodic)
    state = initializer_from_kinetic(K, U, mu, beta, dtau, root_seed; n_stab)
    deltaG, acceptance_sum, nsweeps = 0.0, 0.0, 0
    for _ in 1:burnin
        state, deltaG, acceptance = update_replica!(state, deltaG)
        acceptance_sum += acceptance
        nsweeps += 1
    end
    signs, energies, densities, doubles = zeros(bins), zeros(bins), zeros(bins), zeros(bins)
    binsize = measurements ÷ bins
    for sample in 1:measurements
        for _ in 1:sweeps_per_measurement
            state, deltaG, acceptance = update_replica!(state, deltaG)
            acceptance_sum += acceptance
            nsweeps += 1
        end
        bin = (sample - 1) ÷ binsize + 1
        weight_sign = sign(real(state.sgndetGup * state.sgndetGdn))
        obs = dqmc_observables(state)
        signs[bin] += weight_sign / binsize
        energies[bin] += weight_sign * obs.energy / binsize
        densities[bin] += weight_sign * obs.density / binsize
        doubles[bin] += weight_sign * obs.double_occupancy / binsize
    end
    abs(sum(signs)) > eps() || error("average sign is zero; increase sampling or change parameters")
    ratio_stats(values) = jdqmcm.jackknife((x, s) -> x / s, values, signs)
    energy, energy_error = ratio_stats(energies)
    density, density_error = ratio_stats(densities)
    double_occupancy, double_occupancy_error = ratio_stats(doubles)
    return (; energy, energy_error, energy_per_site = energy / (Lx * Ly),
        energy_per_site_error = energy_error / (Lx * Ly), density, density_error,
        double_occupancy, double_occupancy_error, average_sign = sum(signs) / bins,
        acceptance = acceptance_sum / nsweeps, deltaG_max = deltaG,
        root_seed, burnin, measurements, bins, sweeps_per_measurement,
        bin_data = (; sign = signs, energy = energies, density = densities,
            double_occupancy = doubles))
end
