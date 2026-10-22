using LinearAlgebra

function _ed_parameters(K, U, mu, max_dimension)
    n = size(K, 1)
    size(K) == (n, n) || throw(ArgumentError("K must be square"))
    1 <= n <= 16 || throw(ArgumentError("ED supports 1 to 16 sites"))
    eltype(K) <: Real || throw(ArgumentError("K must be real"))
    all(isfinite, K) && issymmetric(K) ||
        throw(ArgumentError("K must be finite and symmetric"))
    isfinite(U) && isfinite(mu) || throw(ArgumentError("U and mu must be finite"))
    max_dimension >= 1 || throw(ArgumentError("max_dimension must be positive"))
    return n, Matrix{Float64}(K)
end

function _ed_masks(n, particles)
    return [mask for mask in UInt(0):((UInt(1) << n) - UInt(1))
            if count_ones(mask) == particles]
end

_ed_occupied(mask, site) = !iszero(mask & (UInt(1) << (site - 1)))

function _ed_hop_sign(mask, source, destination)
    lo, hi = minmax(source, destination)
    crossings = count(site -> _ed_occupied(mask, site), (lo + 1):(hi - 1))
    return isodd(crossings) ? -1.0 : 1.0
end

function _ed_spin_hops(K, masks)
    index = Dict(mask => i for (i, mask) in enumerate(masks))
    hops = [Tuple{Int, Float64}[] for _ in masks]
    n = size(K, 1)
    for (column, mask) in enumerate(masks), source in 1:n, destination in 1:n
        source == destination && continue
        iszero(K[destination, source]) && continue
        if _ed_occupied(mask, source) && !_ed_occupied(mask, destination)
            newmask = xor(mask, (UInt(1) << (source - 1)) | (UInt(1) << (destination - 1)))
            amplitude = K[destination, source] * _ed_hop_sign(mask, source, destination)
            push!(hops[column], (index[newmask], amplitude))
        end
    end
    return hops
end

function _ed_sector(K, U, mu, Nup, Ndn, max_dimension)
    n, hopping = _ed_parameters(K, U, mu, max_dimension)
    Nup isa Integer && Ndn isa Integer && 0 <= Nup <= n && 0 <= Ndn <= n ||
        throw(ArgumentError("Nup and Ndn must be integers between 0 and the site count"))
    dimension = binomial(n, Nup) * binomial(n, Ndn)
    dimension <= max_dimension || throw(ArgumentError(
        "This sector has $dimension states; the dense ED limit is $max_dimension. " *
        "Use a smaller lattice or raise max_dimension after checking memory."))
    up_masks, down_masks = _ed_masks(n, Nup), _ed_masks(n, Ndn)
    up_hops = _ed_spin_hops(hopping, up_masks)
    down_hops = _ed_spin_hops(hopping, down_masks)
    nd = length(down_masks)
    H = zeros(dimension, dimension)
    double_counts = zeros(dimension)
    for (iu, up) in enumerate(up_masks), (id, down) in enumerate(down_masks)
        column = (iu - 1) * nd + id
        doubles = count_ones(up & down)
        double_counts[column] = doubles
        H[column, column] = U * (doubles - (Nup + Ndn) / 2 + n / 4) - mu * (Nup + Ndn)
        for site in 1:n
            H[column, column] += hopping[site, site] *
                (_ed_occupied(up, site) + _ed_occupied(down, site))
        end
        for (row_up, amplitude) in up_hops[iu]
            H[(row_up - 1) * nd + id, column] += amplitude
        end
        for (row_down, amplitude) in down_hops[id]
            H[(iu - 1) * nd + row_down, column] += amplitude
        end
    end
    basis = (; nsites=n, Nup, Ndn, up_masks, down_masks)
    return H, double_counts, basis
end

"""
    ed_ground_state(K; U=4.0, mu=0.0, Nup=cld(size(K,1),2),
                    Ndn=fld(size(K,1),2), max_dimension=4096)

Solve a small Hubbard model in a fixed `(Nup, Ndn)` sector. `K[i,j]` is the
one-body hopping matrix, including any site energies. The Hamiltonian is
`sum(K[i,j] * c†[i,spin] * c[j,spin]) + U*sum((nup-.5)*(ndn-.5)) - mu*N`.
At `mu=0`, this interaction convention gives half filling on a bipartite lattice.

Return the total `energy`, `energy_per_site`, `density`, mean per-site
`double_occupancy`, normalized `state`, full sector `spectrum`, `basis`, and
eigenvector `residual`. A degenerate ground state returns one eigenvector;
its entanglement can depend on that choice. Dense ED is intended for small
checks. A sector of dimension 4096 already needs 128 MiB per dense matrix.
"""
function ed_ground_state(K; U=4.0, mu=0.0, Nup=cld(size(K, 1), 2),
                         Ndn=fld(size(K, 1), 2), max_dimension=4096)
    H, doubles, basis = _ed_sector(K, U, mu, Nup, Ndn, max_dimension)
    solution = eigen(Symmetric(H))
    energy = first(solution.values)
    state = solution.vectors[:, 1]
    n = basis.nsites
    return (; energy, energy_per_site=energy / n,
            density=(Nup + Ndn) / n,
            double_occupancy=dot(abs2.(state), doubles) / n,
            state, spectrum=solution.values, basis,
            residual=norm(H * state - energy * state))
end

"""
    ed_thermal(K; U=4.0, mu=0.0, beta=2.0, max_sites=6, max_dimension=4096)

Compute grand-canonical thermal observables by diagonalizing every particle
sector. Uses the same Hamiltonian as [`ed_ground_state`](@ref), without a
Trotter approximation. `beta` is the inverse temperature in units where
`k_B=1`; `beta=0` is allowed. Return total `energy`, `energy_per_site`,
`density`, mean per-site `double_occupancy`, `logZ`, and `beta`.
The default six-site limit prevents accidentally starting a large dense job.
"""
function ed_thermal(K; U=4.0, mu=0.0, beta=2.0, max_sites=6, max_dimension=4096)
    n, _ = _ed_parameters(K, U, mu, max_dimension)
    isfinite(beta) && beta >= 0 || throw(ArgumentError("beta must be finite and nonnegative"))
    n <= max_sites || throw(ArgumentError(
        "Thermal ED is limited to $max_sites sites; use a smaller lattice."))
    energies, particles, doubles = Float64[], Float64[], Float64[]
    for Nup in 0:n, Ndn in 0:n
        H, double_counts, _ = _ed_sector(K, U, mu, Nup, Ndn, max_dimension)
        solution = eigen(Symmetric(H))
        append!(energies, solution.values)
        append!(particles, fill(Nup + Ndn, length(solution.values)))
        append!(doubles, transpose(abs2.(solution.vectors)) * double_counts)
    end
    ground_energy = minimum(energies)
    weights = exp.(-beta .* (energies .- ground_energy))
    normalization = sum(weights)
    probabilities = weights ./ normalization
    energy = dot(probabilities, energies)
    return (; energy, energy_per_site=energy / n,
            density=dot(probabilities, particles) / n,
            double_occupancy=dot(probabilities, doubles) / n,
            logZ=log(normalization) - beta * ground_energy, beta)
end

function _ed_compress(mask, sites)
    result = UInt(0)
    for (position, site) in enumerate(sites)
        _ed_occupied(mask, site) && (result |= UInt(1) << (position - 1))
    end
    return result
end

function _ed_partition_sign(up, down, sitesA, sitesB, n)
    # Reorder occupied modes from [up(all), down(all)] to [up(A), down(A),
    # up(B), down(B)]. Count all inversions, including noncontiguous cuts.
    order = vcat(sitesA, sitesA .+ n, sitesB, sitesB .+ n)
    occupied_modes = filter(mode -> mode <= n ? _ed_occupied(up, mode) :
                            _ed_occupied(down, mode - n), order)
    crossings = 0
    for i in eachindex(occupied_modes), j in (i + 1):length(occupied_modes)
        crossings += occupied_modes[i] > occupied_modes[j]
    end
    return isodd(crossings) ? -1.0 : 1.0
end

"""
    ed_renyi2(ground_state, sites; max_subsystem_dimension=4096)

Return `S2 = -log(tr(rho_A^2))` for a ground state returned by
[`ed_ground_state`](@ref). `sites` contains distinct, one-based site numbers;
for example, `[1, 3]` selects the left column of an x-fast 2×2 square.
The spatial partial trace includes the fermionic mode-reordering signs.
"""
function ed_renyi2(ground_state, sites; max_subsystem_dimension=4096)
    basis = ground_state.basis
    n = basis.nsites
    sitesA = collect(sites)
    all(site -> site isa Integer && 1 <= site <= n, sitesA) ||
        throw(ArgumentError("sites must contain integers between 1 and $n"))
    length(unique(sitesA)) == length(sitesA) || throw(ArgumentError("sites must be distinct"))
    sort!(sitesA)
    sitesB = setdiff(1:n, sitesA)
    length(sitesA) > length(sitesB) && ((sitesA, sitesB) = (sitesB, sitesA))
    nA, nB = length(sitesA), length(sitesB)
    dimension_A = 1 << (2 * nA)
    dimension_A <= max_subsystem_dimension || throw(ArgumentError(
        "Reduced density matrix dimension $dimension_A exceeds limit $max_subsystem_dimension"))
    groups = Dict{Int, Vector{Tuple{Int, Float64}}}()
    nd = length(basis.down_masks)
    for (iu, up) in enumerate(basis.up_masks), (id, down) in enumerate(basis.down_masks)
        a = Int(_ed_compress(up, sitesA) | (_ed_compress(down, sitesA) << nA)) + 1
        b = Int(_ed_compress(up, sitesB) | (_ed_compress(down, sitesB) << nB)) + 1
        coefficient = ground_state.state[(iu - 1) * nd + id] *
            _ed_partition_sign(up, down, sitesA, sitesB, n)
        push!(get!(groups, b, Tuple{Int, Float64}[]), (a, coefficient))
    end
    rhoA = zeros(dimension_A, dimension_A)
    for entries in values(groups), (a1, c1) in entries, (a2, c2) in entries
        rhoA[a1, a2] += c1 * c2
    end
    return -log(sum(abs2, rhoA))
end
