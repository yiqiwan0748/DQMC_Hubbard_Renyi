# Comparison with official SmoQyDQMC

Run date: October 22, 2026. This comparison actually ran the registered
[SmoQyDQMC.jl](https://github.com/SmoQySuite/SmoQyDQMC.jl) package, version
2.0.14, in a separate Julia environment and process. Its high-level
`local_updates!` sampled its own Hubbard-Stratonovich fields; its native hopping
and Hubbard energy measurements supplied the external energies. It did not call
the local package sampler.

The public package tested was clean commit
`6f8b949275a648f9592410085fd8263b1286f2ba` of `DQMC_Hubbard_Renyi`.
The SmoQyDQMC registered source tree hash was
`0d1f705a9a20cdabe2ccb1e7601b55dfe8425e75`. The published method is
[Cohen-Stead et al., SciPost Physics Codebases 29 (2024)](https://scipost.org/10.21468/SciPostPhysCodeb.29).

## Matched model and sampling

Both engines used the centered Hubbard Hamiltonian
`H = sum(K_ij cdag_i c_j) + U sum((nup-1/2)(ndn-1/2)) - mu N`, in the
grand-canonical ensemble, with `t=1`, `mu=0`, `beta=2`, `dtau=0.05`, and
`U=0` or `U=2`. The open 2-by-2 square has exactly four undirected hopping
bonds, `(1,2)`, `(1,3)`, `(2,4)`, `(3,4)`:

```text
K = [ 0 -1 -1  0
     -1  0  0 -1
     -1  0  0 -1
      0 -1 -1  0 ]
```

Both scripts assert this matrix explicitly. SmoQyDQMC requires periodic cell
coordinates, so its model uses one four-orbital cell with only the four internal
zero-displacement bonds. There are no bonds between cells. This yields exactly
the open square Hamiltonian above, with no duplicate periodic length-two bonds.
Its metadata records both the computational cell flags and physical boundaries.

Both engines use `AsymExactPropagator{Float64,Float64}`, with
`symmetric=false`, `checkerboard=false`, and the same finite-step definition
`B = exp(-dtau V) exp(-dtau K)`. Thus this comparison does not mix different
Trotter or checkerboard approximations. Stabilization occurs every 10 slices.

Each engine used four independent, disjoint seeds, each with 1,000 discarded
sweeps followed by 20,000 measurements, one local sweep per measurement.
Each stream stored 40 bins of 500 measurements. The table pools 160 bins per
engine and interaction strength; errors are sign-reweighted delete-one-bin
jackknife standard errors. The largest stabilization correction was
`2.04e-12` for the local package and `1.05e-12` for SmoQyDQMC.

## Observed agreement

| U | Observable per site | Public package | SmoQyDQMC 2.0.14 | Difference / combined error |
|---:|---|---:|---:|---:|
| 0 | Energy | -0.964027580075812 | -0.964027580075812 | Roundoff only |
| 0 | Density | 1 | 1 | Exact to roundoff |
| 0 | Double occupancy | 0.25 | 0.25 | Exact to roundoff |
| 2 | Energy | -1.07348056 +/- 0.00082467 | -1.07312977 +/- 0.00079045 | -0.3071 |
| 2 | Density | 1 | 1 | Exact to roundoff |
| 2 | Double occupancy | 0.16116259 +/- 0.00015229 | 0.16119075 +/- 0.00015549 | -0.1294 |

The combined error is `sqrt(error_local^2 + error_SmoQy^2)`, since the two
engines used independent streams. Average sign was one in every case.
At `U=0`, the analytic values are `E/site=-tanh(2)`, density `1`, and double
occupancy `1/4`; both engines agree within `6e-15` for energy and roundoff for
the other observables. `comparison.csv` retains the complete precision,
roundoff floor, thresholds, and pass flags for all six comparisons.

Doubling the bin length to 1,000 sweeps changes the interacting energy error
from `0.00082467` to `0.00082640` locally and from `0.00079045` to `0.00081726`
in SmoQyDQMC. The corresponding double-occupancy errors change from
`0.00015229` to `0.00015393`, and from `0.00015549` to `0.00016248`.
Per-seed estimates and seed-to-seed errors are also retained.

## Independent QuSpin observable reference

We also used [QuSpin 1.0.1](https://quspin.github.io/QuSpin/) to construct the
full 256-dimensional many-body Hamiltonian in a separate Python process. Its
kinetic and interaction matrices give the exact finite-step transfer matrix
`T=exp(-dtau*V) exp(-dtau*K)`. Direct traces of `T^40` provide an independent
reference at the same `beta=2` and `dtau=0.05`, without either JDQMC backend.
The reference output is [`../ed/quspin_trotter.json`](../ed/quspin_trotter.json).

| U=2 observable per site | QuSpin exact finite step | Public package residual / sampling error | SmoQyDQMC residual / sampling error |
|---|---:|---:|---:|
| Energy | -1.073273540776198 | -0.2510 | +0.1819 |
| Density | 1 | Exact to roundoff | Exact to roundoff |
| Double occupancy | 0.161135391813110 | +0.1786 | +0.3561 |

Both engines also match the `U=0` QuSpin values within `5e-15`. All 12
engine/interaction/observable comparisons passed. `quspin_comparison.csv`
retains every result and `quspin_comparison_metadata.json` records input hashes,
the physical hopping matrix and its SHA-256, package versions and thresholds.
The comparator checks matching geometry, Hamiltonian parameters, ensemble and
propagator settings before comparing values. Deterministic limits use a
`1e-10` absolute roundoff floor.

QuSpin also gives the continuous-time energy/site `-1.072447929941280` and
double occupancy `0.161384971037859` at `U=2`. Their differences from the
finite-step values demonstrate why matching the finite-step trace matters.
This independent reference corroborates these specific thermal observables;
it does not establish correctness of every shared backend routine.

## Scope

This is an external package comparison of small-square thermal observables.
The engines share `JDQMCFramework 1.2.6`, `JDQMCMeasurements 1.9.2`, and
`LatticeUtilities 2.0.7`; it is therefore not an independent validation of those
shared low-level libraries. The external model construction, high-level update
implementation, native energy measurements, and random streams are separate.
It does not validate the replicated Renyi estimator or imply any speed advantage
over SmoQyDQMC. Finite-step bias remains when comparing either engine to the
continuous-time Hubbard model.

The run used Julia 1.12.1 on Linux/WSL, an AMD Ryzen 7 7700X, one Julia thread,
and one OpenBLAS thread. Each engine's metadata records package versions,
source hashes, lattice definitions, hardware, seeds, and numerical settings.
`smoqy_environment/` contains the exact external Project and Manifest.

## Reproduction

Run these PowerShell-compatible commands from the repository root, using a new
output directory.

```powershell
# Restore the isolated official package environment without changing the root project.
julia --startup-file=no --project=.\results\external_packages_20261022\dqmc\smoqy_environment -e 'using Pkg; Pkg.instantiate()'

# Use --project=. when running from the published package checkout.
julia --startup-file=no --project=. .\analysis\external_packages\dqmc_local.jl .\results\external_packages_new\dqmc
julia --startup-file=no --project=.\results\external_packages_20261022\dqmc\smoqy_environment .\analysis\external_packages\dqmc_smoqy.jl .\results\external_packages_new\dqmc
julia --startup-file=no --project=.\results\external_packages_20261022\dqmc\smoqy_environment .\analysis\external_packages\summarize_dqmc.jl .\results\external_packages_new\dqmc

# After generating the QuSpin reference, compare both engines against its exact trace.
# This comparator needs Python 3.11 or later and uses only the standard library.
python .\analysis\external_packages\dqmc_quspin_compare.py .\results\external_packages_20261022\ed\quspin_trotter.json .\results\external_packages_new\dqmc .\results\external_packages_new\dqmc\quspin_comparison.csv
```

For this local research checkout, the tested public package project was
`.\output\github_release`, so its local-engine command used
`--project=.\output\github_release` instead of `--project=.`.
The scripts refuse to replace existing evidence.
