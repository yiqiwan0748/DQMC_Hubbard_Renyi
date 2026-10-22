"""Optional ITensorMPS backend, loaded by `dmrg_ground_state` on first use."""
module HubbardDMRG

using ITensors
using ITensorMPS
using LinearAlgebra
using Random

function validate_inputs(K, U, mu, Nup, Ndn, cut, nsweeps, maxdim, cutoff)
    N = size(K, 1)
    size(K, 2) == N && N >= 2 || throw(ArgumentError("K must be square with at least two sites"))
    all(isfinite, K) && issymmetric(K) || throw(ArgumentError("K must be finite, real, and symmetric"))
    isfinite(U) && isfinite(mu) || throw(ArgumentError("U and mu must be finite"))
    0 <= Nup <= N && 0 <= Ndn <= N || throw(ArgumentError("Nup and Ndn must lie between 0 and the number of sites"))
    1 <= cut < N || throw(ArgumentError("cut must lie between 1 and N-1"))
    nsweeps >= 2 || throw(ArgumentError("Use at least two sweeps to monitor convergence"))
    !isempty(maxdim) && all(>(0), maxdim) || throw(ArgumentError("maxdim must contain positive integers"))
    isfinite(cutoff) && cutoff > 0 || throw(ArgumentError("cutoff must be finite and positive"))
    return N
end

function initial_occupations(N, Nup, Ndn)
    # Spread particles over the chain and start half filling in a Neel state.
    up_order = vcat(collect(1:2:N), collect(2:2:N))
    dn_order = vcat(collect(2:2:N), collect(1:2:N))
    up = falses(N)
    dn = falses(N)
    up[up_order[1:Nup]] .= true
    dn[dn_order[1:Ndn]] .= true
    return [up[i] ? (dn[i] ? "UpDn" : "Up") : (dn[i] ? "Dn" : "Emp") for i in 1:N]
end

function hubbard_mpo(K, sites, U, mu)
    N = size(K, 1)
    terms = OpSum()
    for i in 1:N, j in (i + 1):N
        iszero(K[i, j]) && continue
        for (create, destroy) in (("Cdagup", "Cup"), ("Cdagdn", "Cdn"))
            terms += K[i, j], create, i, destroy, j
            terms += K[i, j], create, j, destroy, i
        end
    end
    for i in 1:N
        # U(nup-1/2)(ndn-1/2) = U nup ndn - U Ntot/2 + U/4.
        # Keeping the constant makes the reported energy match ED and DQMC.
        terms += U, "Nupdn", i
        terms += K[i, i] - mu - U / 2, "Ntot", i
        terms += U / 4, "Id", i
    end
    return MPO(terms, sites)
end

function prefix_renyi2(psi, cut)
    orthogonalize!(psi, cut)
    sites = siteinds(psi)
    left = cut == 1 ? (sites[cut],) : (linkind(psi, cut - 1), sites[cut])
    _, singular, _ = svd(psi[cut], left)
    weights = [real(singular[i, i])^2 for i in 1:dim(singular, 1)]
    weights ./= sum(weights)
    return -log(sum(abs2, weights))
end

"""
    dmrg_ground_state(K; U=4, mu=0, Nup=cld(N,2), Ndn=fld(N,2), cut=N÷2,
                      nsweeps=10, maxdim=[16,32,64,128], cutoff=1e-12,
                      root_seed=1234, outputlevel=0)

Find a fixed-particle-number Hubbard ground state with ITensorMPS DMRG.
The default sector has one electron per site, including for odd site counts.
`K` is a real symmetric single-particle matrix, with hopping `-t` on a bond.
The interaction is centered: `U*(nup-1/2)*(ndn-1/2) - mu*(nup+ndn)`.
`cut` gives the prefix subsystem `1:cut`, in the input site order. Square
lattices use the same x-fast order as the ED and DQMC examples.

Increase the bond dimensions and sweeps, and compare energies and entropies,
before using a large-system result. A small energy variance alone does not
establish that the state is the ground state. This wrapper uses ITensorMPS;
it does not implement a separate DMRG engine.
"""
function dmrg_ground_state(K::AbstractMatrix{<:Real}; U::Real=4.0, mu::Real=0.0,
        Nup::Integer=cld(size(K, 1), 2), Ndn::Integer=fld(size(K, 1), 2),
        cut::Integer=size(K, 1)÷2, nsweeps::Integer=10,
        maxdim::AbstractVector{<:Integer}=[16, 32, 64, 128],
        cutoff::Real=1e-12, root_seed::Integer=1234, outputlevel::Integer=0)
    N = validate_inputs(K, U, mu, Nup, Ndn, cut, nsweeps, maxdim, cutoff)
    sites = siteinds("Electron", N; conserve_qns=true)
    H = hubbard_mpo(K, sites, U, mu)
    rng = Xoshiro(root_seed)
    psi0 = random_mps(rng, sites, initial_occupations(N, Nup, Ndn); linkdims=min(8, first(maxdim)))
    observer = DMRGObserver()
    noise = nsweeps == 2 ? [1e-8, 0.0] : [1e-7, 1e-8, 0.0]
    elapsed = @elapsed _, psi = dmrg(H, psi0; nsweeps=Int(nsweeps),
        maxdim=Int.(maxdim), cutoff=Float64(cutoff), noise,
        eigsolve_krylovdim=8, eigsolve_tol=1e-12, observer,
        outputlevel=Int(outputlevel))
    normalize!(psi)
    energy = real(inner(psi', H, psi))
    raw_variance = real(inner(H, psi, H, psi)) - energy^2
    sweep_energies = Float64.(energies(observer))
    truncation_errors = Float64.(truncerrors(observer))
    return (
        energy=energy, S2=prefix_renyi2(psi, Int(cut)), state=psi,
        Nup=Int(Nup), Ndn=Int(Ndn), cut=Int(cut),
        energy_variance=max(0.0, raw_variance),
        sweep_energy_change=abs(sweep_energies[end] - sweep_energies[end - 1]),
        maxlinkdim=maxlinkdim(psi),
        final_truncation_error=last(truncation_errors),
        sweep_energies=sweep_energies, truncation_errors=truncation_errors,
        nsweeps=Int(nsweeps), maxdim=Int.(maxdim), cutoff=Float64(cutoff),
        root_seed=Int(root_seed), elapsed_seconds=elapsed,
        density_up=expect(psi, "Nup"), density_dn=expect(psi, "Ndn"),
    )
end

end # module HubbardDMRG
