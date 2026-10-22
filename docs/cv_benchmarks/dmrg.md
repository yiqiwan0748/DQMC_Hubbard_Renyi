# Hubbard ground states with DMRG

The DMRG example finds a ground state and its second Renyi entropy with
[ITensorMPS.jl](https://github.com/ITensor/ITensorMPS.jl). It wraps that library;
the tensor-network optimizer is ITensorMPS's implementation of DMRG. The new
code supplies a common Hubbard Hamiltonian, editable lattice settings, and
checks against exact solutions. It needs no Monte Carlo settings or experience.

## Run a small square

From the repository directory in PowerShell:

```powershell
julia --project=. .\examples\dmrg_square.jl
julia --project=. .\examples\dmrg_square.jl .\examples\dmrg_square.toml --output .\results\my_dmrg_square.csv
```

The first command displays the result. The second writes a new CSV and refuses
to replace an existing file. Install the project dependencies with
`julia --project=. -e 'using Pkg; Pkg.instantiate()'` if needed.

Edit `examples/dmrg_square.toml` to set `Lx`, `Ly`, hopping `t`, interaction `U`,
and the two periodic-boundary flags. The supplied example is an open 2 by 2
square, with two up and two down electrons. For an open 3 by 2 rectangle at half
filling, set `Lx = 3`, `Ly = 2`, `Nup = 3`, `Ndn = 3`, and `cut = 3`.
Both methods use site number `x + (y-1)Lx`; entropy refers to sites `1:cut`.
Use a non-square lattice or custom hopping by passing its real symmetric
single-particle matrix directly to `dmrg_ground_state(K; ...)`.

This is a zero-temperature calculation in the specified particle sector.
For larger lattices, raise the bond dimensions and number of sweeps and compare
both energy and entropy. Narrow open strips are practical starting points for
DMRG. Two-dimensional width increases the required bond dimensions rapidly;
large fully periodic squares need particular care. The interface accepting
such a lattice does not establish that a chosen calculation has converged.
This follows the official
[ITensorMPS convergence guidance](https://itensor.github.io/ITensorMPS.jl/stable/faq/DMRG.html).

## Hamiltonian and entropy

All new DQMC, ED, and DMRG examples use

\[
H=\sum_{ij\sigma}K_{ij}c^\dagger_{i\sigma}c_{j\sigma}
+U\sum_i(n_{i\uparrow}-\tfrac12)(n_{i\downarrow}-\tfrac12)
-\mu\sum_i(n_{i\uparrow}+n_{i\downarrow}).
\]

Each undirected lattice bond contributes hopping `-t` in both directions.
Length-two periodic axes use one bond between the two sites. The centered
interaction is implemented using the exact identity
`U*nup*ndn - (U/2)*Ntot + U/4` at every site. Relative to the conventional
`U*nup*ndn` Hamiltonian in a fixed-particle sector,
`E_centered = E_conventional - U*(Nup+Ndn)/2 + U*Nsites/4`.
No estimator or entropy definition changes: `S2 = -log(tr(rho_A^2))`.
DMRG obtains it from normalized squared Schmidt singular values.

## Exact reference and reproducibility

The named analytic reference is the **exact Hubbard dimer solution** described
by Carrascal, Ferrer, Smith, and Burke,
[*The Hubbard dimer: a density functional case study of a many-body problem*](https://arxiv.org/abs/1502.02194),
J. Phys.: Condens. Matter **27**, 393001 (2015).
For the symmetric half-filled dimer,

\[
E_{\mathrm{centered}}=-\tfrac12\sqrt{U^2+16t^2},\qquad
q=\tfrac12\left(1-\frac{U}{\sqrt{U^2+16t^2}}\right),\qquad
S_2=-\log\frac{q^2+(1-q)^2}{2}.
\]

At `U/t = 4`, these give energy `-2.8284271247461903 t` and one-site entropy
`0.9808292530117262`. The dimer checks the interaction normalization and entropy
against formulas independent of either numerical solver.
The square and rectangular clusters are then checked against a separately
constructed full Fock-space Hamiltonian and Julia/LAPACK exact diagonalization.
Their ED Hilbert-space dimensions are 36 and 400 respectively. This small-system
comparison is specifically recommended in the official ITensorMPS convergence
guidance. Increasing sweep counts and bond caps and changing the starting seed
checks the stability of both energy and entropy.

```powershell
julia --project=. .\examples\dmrg_benchmark.jl .\results\dmrg_validation_new.csv
julia --project=. .\test\test_dmrg_package.jl
```

The benchmark runner writes all reference values, absolute errors, particle
counts, sweep counts, bond caps, seeds, truncation errors, and energy variances
to the requested new CSV. Its elapsed times include compilation when applicable
and are diagnostic times, not a DMRG performance claim. A small energy variance
checks an eigenstate; agreement with ED checks that it is the ground state on
these clusters.

Measured results from 21 October 2026 are saved in
`results/cv_benchmarks_20261021/dmrg_benchmark_v1.csv`. All cases used `t=1`,
`U=4`, open boundaries, equal up/down counts at half filling, and cutoff
`1e-12`. Julia was version 1.12.1, with one Julia thread and one OpenBLAS thread.

| Lattice | Reference energy | Largest energy error | Largest S2 error |
| --- | ---: | ---: | ---: |
| 2-site dimer | -2.8284271247461903 | 0 at printed precision | 0 at printed precision |
| 2 by 2 square | -6.102748483462078 | 6.22e-15 | 3.32e-13 |
| 3 by 2 rectangle | -9.61932132395754 | 1.60e-14 | 5.20e-14 |

The last two rows include both the initial run (8 sweeps, seed 1234) and a
convergence run (12 sweeps, seed 5678) with a higher bond cap. The largest
actual MPS bond dimensions were 16 and 64. These are small clusters where the
MPS can represent the full state. This validates the wrapper and entropy
measurement; it is not evidence of large-system DMRG convergence.
The focused package tests passed 14 of 14 checks, including analytic dimer
energy and entropy, particle counts away from half filling, square-cluster ED
agreement, and input validation. Their log is saved alongside the CSV as
`dmrg_tests.log`.

The algorithm reference is White,
[*Density matrix formulation for quantum renormalization groups*](https://doi.org/10.1103/PhysRevLett.69.2863),
Phys. Rev. Lett. **69**, 2863 (1992). This paper identifies the algorithm; the
exact dimer and full diagonalization are the references actually tested here.

Suggested CV wording, supported by the saved validation results:
“Built an ITensorMPS workflow for Hubbard-model ground states and Renyi entropy;
validated against the exact Hubbard dimer and full exact diagonalization on
small square-lattice clusters, with energy errors below `2e-14 t`.”
