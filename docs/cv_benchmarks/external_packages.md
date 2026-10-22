# Comparisons with established packages

All three comparisons below actually ran the external package. They tested
public source commit `6f8b949275a648f9592410085fd8263b1286f2ba` on October 22,
2026. The new benchmark scripts and evidence do not change the numerical core.

| Interface | External reference | Tested scope | Observed agreement |
|---|---|---|---|
| DQMC | SmoQyDQMC.jl 2.0.14 | Open 2x2 square, U/t=0 and 2, beta*t=2, dtau*t=0.05; energy, density and double occupancy | All 6 comparisons passed; largest interacting difference 0.31 combined standard errors |
| ED | QuSpin 1.0.1 | Seven particle-sector spectra (546 eigenvalues), 17 spatial second Renyi entropies, six grand-canonical thermal cases | All 108 grouped checks passed; largest spectrum difference 2.85e-14 t and S2 difference 3.36e-14 |
| ITensorMPS-based DMRG | TeNPy 1.1.1 | Open 2x2, half-filled and one-hole 3x2, and 6x1; two convergence settings each | All 32 observable checks passed; largest total-energy difference 1.25e-14 t and S2 difference 3.20e-13 |

These are finite-cluster checks of specified observables, rather than a claim
of convergence on larger lattices. Each external code built its own model and
state. Inputs, dependency versions, execution logs, source hashes, tolerances
and per-observable results are saved with the reports below.

## DQMC: official SmoQyDQMC and an independent exact reference

The official SmoQyDQMC sampler ran in its own Julia environment, using its
high-level local updates and native energy measurements. Each engine used
four disjoint random seeds with 1,000 discarded sweeps and 20,000 measurements
per seed. The comparison uses sign-reweighted bin-jackknife standard errors;
doubling the bin length gives similar errors. Both average signs were one.

At U/t=2 the public package gave energy/site
`-1.07348056 +/- 0.00082467`; SmoQyDQMC gave
`-1.07312977 +/- 0.00079045`. Their difference is 0.307 combined standard
errors. Double occupancy differs by 0.129 combined standard errors. The free
case agrees with the analytic answer to roundoff.

The engines share JDQMCFramework 1.2.6, JDQMCMeasurements 1.9.2 and
LatticeUtilities 2.0.7. To check against an independent implementation too,
QuSpin constructed the full grand-canonical operators and SciPy evaluated
the exact finite-step trace of
`[exp(-dtau V) exp(-dtau K)]^(beta/dtau)`. Its interacting reference is
energy/site `-1.0732735407761982`, double occupancy
`0.16113539181310976`, and density one. All 12 comparisons of the two
Monte Carlo engines with this reference passed; the largest interacting
residual is 0.36 sampling standard errors.

This finite-step reference matches the sampled discretization. Comparing to
the continuous-time Hamiltonian instead would also include Trotter bias.
These tests cover thermal observables, and do not validate the replicated
Renyi estimator or measure performance against SmoQyDQMC.

[DQMC settings, evidence and reproduction](../../results/external_packages_20261022/dqmc/README.md).

## ED: independent QuSpin bases and operators

QuSpin built its own spinful-fermion bases, centered Hubbard operators,
eigensystems and spatial partial traces. Cases include open squares,
a periodic rectangle, custom non-nearest-neighbor hopping and site energies,
repulsive and attractive interactions, a nonzero chemical potential, and
vacuum/full sectors. Noncontiguous cuts test fermionic subsystem ordering.
Ground-state observables are excluded when the ground state is degenerate.

Thermal checks compare energy, density, double occupancy and log partition
function, including an infinite-temperature case. The 108 grouped checks
also include matched hopping matrices, Hermiticity and eigenvector residuals.
The largest residual or comparison error is `1.18e-13`, below the `2e-10`
general tolerance. The full-spectrum maximum error is `2.85e-14 t`.

[ED cases, evidence and reproduction](../../results/external_packages_20261022/ed/REPORT.md).

## DMRG: independent TeNPy optimization

TeNPy used its stock FermiHubbardModel and TwoSiteDMRGEngine, independently
of the public ITensorMPS optimizer. Site order, every hopping bond, particle
counts, and the entropy cut were explicitly matched. The centered interaction
was obtained from TeNPy's uncentered convention by setting
`mu_TeNPy=mu+U/2` and adding `Nsites*U/4` to its energy. This is the identity
`U*(nup-1/2)*(ndn-1/2)=U*nup*ndn-U*(nup+ndn)/2+U/4`;
the state and entanglement are unchanged.

Both packages ran standard and tighter bond-dimension/cutoff schedules from
different initial states. Energy and entropy remained stable. These runs
check small open systems with at most six sites; they do not establish
larger-system or periodic-boundary convergence.

[DMRG settings, evidence and reproduction](../../results/external_packages_20261022/dmrg/REPORT.md).

## CV wording and timing scope

"Cross-checked square-lattice DQMC thermodynamics against SmoQyDQMC, ED
against QuSpin, and ITensorMPS ground states and Renyi entropy against TeNPy
on small Hubbard lattices."

The earlier **4.59x** timing result remains a measured local rank-one Grover
update-kernel comparison with direct Grover-matrix rebuilds. It is not an
end-to-end speedup or a speed comparison with an external package. Its saved
inputs and numerical source hashes are unchanged; see the
[original timing report](dqmc.md).

## Reference methods

- [SmoQyDQMC: Cohen-Stead et al., SciPost Physics Codebases 29 (2024)](https://scipost.org/10.21468/SciPostPhysCodeb.29).
- [QuSpin: Weinberg and Bukov, SciPost Physics 7, 020 (2019)](https://scipost.org/10.21468/SciPostPhys.7.2.020).
- [TeNPy: Hauschild et al., Tensor Network Python version 1 (2024)](https://arxiv.org/abs/2408.02010).
