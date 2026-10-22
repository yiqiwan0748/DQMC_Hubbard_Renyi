# ED comparison with QuSpin

On 2026-10-22, the public `HubbardRenyiDQMC` ED API at commit
`6f8b949275a648f9592410085fd8263b1286f2ba` passed **108 of 108 grouped checks
against QuSpin 1.0.1**. The comparison covers seven fixed-particle-sector cases
(546 eigenvalues), six grand-canonical thermal cases, and 17 spatial Renyi-2 cuts.
QuSpin supplied the fermion basis, Hamiltonian and observable operators,
diagonalization interface, and fermionic partial trace.

| Quantity | Largest absolute difference |
| --- | ---: |
| Ground or thermal total energy | 1.4211e-14 |
| Full fixed-sector spectrum | 2.8422e-14 |
| Density per site | 9.9921e-16 |
| Double occupancy per site | 4.0246e-16 |
| Spatial S2 | 3.3529e-14 |
| Thermal log Z | 1.0659e-14 |
| Kinetic matrix entries | 0 |

The absolute comparison tolerance was 2e-10; kinetic matrices used 1e-14.
All external Hamiltonians were Hermitian to the 1e-13 check threshold.
The largest external ground-state residual was 2.9766e-14. The largest full
thermal eigensystem Frobenius residual was 1.1763e-13. These residuals measure
the external solver's numerical accuracy and are listed separately from
cross-package differences.

## Matched physical inputs

Both implementations used

`H = sum(K[i,j] cdag[i,s] c[j,s]) + U sum((nup-.5)(ndn-.5)) - mu*N`.

`energy` includes the centered interaction constant and the chemical-potential
term. Density and double occupancy are per physical site. QuSpin's `z|z`
operator directly implements the centered interaction. Square and rectangle
sites use x-fast order, `i = x + (y-1)*Lx`, with each undirected bond included
once. Full square-lattice bond lists were supplied to the site-indexed
`spinful_fermion_basis_1d` Fock basis without symmetry reduction.

The complete inputs, including every cut, are in
[ed_cases.toml](../../../analysis/external_packages/ed_cases.toml).

| Fixed-sector case | Boundary / hopping | U | mu | (Nup,Ndn) |
| --- | --- | ---: | ---: | --- |
| 2x2 square | Open, t=1 | 4 | 0 | (2,2) |
| 3x2 rectangle | Open, t=1 | 4 | 0 | (3,3) |
| 3x2 rectangle | Periodic both axes, t=0.7 | 2 | 0.13 | (2,1) |
| Custom three-site K | Matrix below | 3.1 | 0.23 | (2,1) |
| Custom three-site K | Matrix below | -2 | 0.23 | (1,1) |
| Custom three-site K, vacuum | Matrix below | 3.1 | 0.23 | (0,0) |
| Custom three-site K, full | Matrix below | 3.1 | 0.23 | (3,3) |

```text
K = [ 0.1  -0.7   0.2
     -0.7  -0.2  -1.1
      0.2  -1.1   0.3 ]
```

The thermal cases sum all particle sectors: open 2x2 at `(U,mu,beta)` equal
to `(2,0,2)`, `(0,0,2)`, `(4,0.37,1.3)`, and `(4,0.37,0)`; and the custom
matrix at `(3.1,0.23,0.8)` and `(-2,0.23,0.8)`.

Spatial entropies use QuSpin's `alpha=2`, `density=False`, and the same
physical sites for both spin species. Cuts include contiguous, noncontiguous,
single-site, and complementary subsystems. All nontrivial ground states in
this run were nondegenerate. The script excludes eigenvector-dependent
ground-state observables if the gap is at most 1e-9; that exclusion was not
needed in this run. No entropy cut was requested for the periodic doped case.

## Independent finite-step reference for DQMC

The separate [Trotter helper](../../../analysis/external_packages/ed_quspin_trotter.py)
uses QuSpin's full 256-state matrices and SciPy matrix exponentials to evaluate
`T = exp(-dtau*V) exp(-dtau*K)`, `rho = T^40`, and
`<O> = trace(O*rho)/trace(rho)`. This matches the asymmetric finite-step
DQMC propagator at beta=2 and dtau=0.05 on the open four-bond 2x2 square.

| U | Finite-step energy/site | Finite-step double occupancy | Continuum energy/site | Continuum double occupancy |
| --- | ---: | ---: | ---: | ---: |
| 0 | -0.9640275800758169 | 0.25 | -0.9640275800758162 | 0.2500000000000002 |
| 2 | -1.0732735407761982 | 0.16113539181310976 | -1.0724479299412795 | 0.16138497103785854 |

Both finite-step densities are exactly 1 to the reported precision. The U=2
finite-step log Z is 11.563703768818065. The U=0 finite-step and continuum
energies agree within 3e-15. All values are in [quspin_trotter.json](quspin_trotter.json).

## Environment and reproduction

The successful run used isolated WSL Linux Python 3.12.15, QuSpin 1.0.1,
quspin-extensions 0.1.8, parallel-sparse-tools 0.2.5, NumPy 2.5.4, and
SciPy 1.18.1. Julia was 1.12.1 and the tested package version was 0.1.0.
All BLAS, OpenMP, and Numba thread settings were 1. Dependencies are pinned
in [requirements_linux_frozen.txt](requirements_linux_frozen.txt).
The earlier Windows wheel installation could not import because
`libomp140.x86_64.dll` was absent. The validated results below came from the
successful Linux environment, which required no DLL or library-path workaround.

These commands can be launched from PowerShell. They use Julia and Python
3.12 inside WSL, an isolated Python environment, and a new evidence directory.
Replace `/path/to/DQMC_Hubbard_Renyi` with the Linux path of a public-repository
clone. Instantiate the root Julia project first, as described in the README.

```powershell
wsl.exe -d Ubuntu -- bash -lc 'cd /path/to/DQMC_Hubbard_Renyi && python3.12 -m venv output/external_packages_repeat/quspin_venv && output/external_packages_repeat/quspin_venv/bin/python -m pip install -r results/external_packages_20261022/ed/requirements_linux_frozen.txt'
wsl.exe -d Ubuntu -- bash -lc 'cd /path/to/DQMC_Hubbard_Renyi && OPENBLAS_NUM_THREADS=1 julia --project=. analysis/external_packages/ed_julia.jl analysis/external_packages/ed_cases.toml results/ed_repeat/julia_results.toml'
wsl.exe -d Ubuntu -- bash -lc 'cd /path/to/DQMC_Hubbard_Renyi && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 output/external_packages_repeat/quspin_venv/bin/python analysis/external_packages/ed_quspin.py analysis/external_packages/ed_cases.toml results/ed_repeat/julia_results.toml results/ed_repeat'
wsl.exe -d Ubuntu -- bash -lc 'cd /path/to/DQMC_Hubbard_Renyi && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 output/external_packages_repeat/quspin_venv/bin/python analysis/external_packages/ed_quspin_trotter.py results/ed_repeat/quspin_trotter.json'
```

The recorded run used the research checkout's `output/github_release` Julia
project and workspace-local Python runtime. The commands above restore its
Python dependencies without requiring those local directories. The comparison
scripts refuse to replace existing result files.

Evidence: [summary](summary.json), [all comparisons](comparison.csv),
[Julia output](julia_results_final.toml), [QuSpin output](quspin_results.json),
[Linux run log](quspin_linux_run.log), and
[source and script SHA-256 hashes](provenance.json). The finite-step JSON also
records its script hashes and dependency versions. These are small-cluster
validation results for real hopping matrices; they do not establish accuracy
for untested model sizes or parameters.

Reference package: [QuSpin documentation](https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html)
and [Weinberg and Bukov, SciPost Physics 7, 020 (2019)](https://scipost.org/10.21468/SciPostPhys.7.2.020).
