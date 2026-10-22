# HubbardRenyiDQMC

Julia tools for the Hubbard model: finite-temperature determinantal quantum
Monte Carlo (DQMC), small-system exact diagonalization (ED), and ground-state
DMRG through ITensorMPS. Start with a square lattice, edit a few settings, and
compare with an exact answer before increasing the size.

The DQMC code uses the SmoQySuite libraries **JDQMCFramework**,
**JDQMCMeasurements**, and **LatticeUtilities**, following their explicit
lattice, propagator, and stabilized Green-function conventions. The DMRG
optimizer is **ITensorMPS**; this package supplies the model and workflow.

## Install

Install Julia 1.12.1, then run these commands in PowerShell. The first setup
may take several minutes while Julia installs and compiles dependencies.

```powershell
git clone https://github.com/yiqiwan0748/DQMC_Hubbard_Renyi.git
Set-Location .\DQMC_Hubbard_Renyi
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Try an exact answer first

```powershell
julia --project=. .\examples\ed_square.jl
```

This solves an open 2×2 square. It displays the ground-state energy and
second Rényi entropy, then computes the thermal energy, density, and double
occupancy. Edit the size, hopping, boundaries, and interaction in that short
example. ED has no sampling error, but its memory cost grows quickly.
Thermal ED defaults to at most six sites.

## Run square-lattice DQMC

```powershell
julia --project=. .\run\run_square_hubbard.jl .\parameters\square_4x4.toml
```

The supplied example is a half-filled 4×4 square. It writes `summary.csv`,
`bins.csv`, and a copy of the input settings to a new directory. For another
run, choose a new directory:

```powershell
julia --project=. .\run\run_square_hubbard.jl .\parameters\square_4x4.toml --output-dir .\results\my_square_run
```

The output directory is reserved before sampling starts. If a run fails, its
input settings remain there; choose a new directory when retrying. A directory
with parameters but no summary is an incomplete run.

Edit `parameters/square_4x4.toml`:

| Setting | Meaning |
| --- | --- |
| `Lx`, `Ly` | Number of sites along each direction. |
| `periodic` | Wrap each direction: `[true, true]` is a torus; `[false, false]` has open edges. |
| `t`, `U`, `mu` | Hopping, repulsion, and chemical potential. Use `t=1` as the energy unit. |
| `beta` | Inverse temperature. Larger values mean lower temperature. |
| `dtau` | Imaginary-time step. Halve it to check discretization bias; `beta/dtau` must be an integer. |
| `root_seed` | Seed for reproducing the random stream. Use independent seeds to check stability. |
| `burnin` | Initial sweeps discarded before collecting measurements. |
| `measurements`, `bins` | Samples and groups used for jackknife errors. Measurements must divide evenly into at least two bins. |

The errors describe sampling uncertainty. Check them by increasing burn-in,
measurements, and bin size. Inspect `average_sign`: a small value makes the
result much harder to resolve. Even bipartite squares at `mu=0` provide a
useful starting point without a sign problem. Fully periodic odd lengths can
break bipartiteness.

This square-lattice command measures energy, density, and double occupancy.
The established Rényi lambda-grid runner remains the separate 1D workflow in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md). These square-lattice thermodynamic
checks do not validate a square-lattice entropy estimator.

## Run DMRG on the same small lattice

```powershell
julia --project=. .\examples\dmrg_square.jl
```

Edit `examples/dmrg_square.toml` to change the lattice, particle counts, or
solver settings. DMRG finds a zero-temperature state in a fixed particle
sector. Raise the bond dimensions and sweeps until both energy and entropy
stop changing. Start with a small square or narrow open strip; wider 2D
systems need much larger bond dimensions.

## Model conventions

All three new interfaces use

```text
H = sum(K[i,j] * c†[i,spin] * c[j,spin])
    + U * sum((n_up - 1/2) * (n_down - 1/2))
    - mu * sum(n_up + n_down).
```

Nearest-neighbor entries of `K` are `-t`. Each undirected bond occurs once,
including periodic axes of length two; a length-one axis has no self hopping.
Sites are numbered `x + (y-1)*Lx`. A custom real symmetric `K` can be passed
to `initializer_from_kinetic`, `ed_ground_state`, `ed_thermal`, or
`dmrg_ground_state`.

DQMC and thermal ED are grand canonical at finite temperature. Ground-state
ED and DMRG use specified particle counts. Match the Hamiltonian, boundaries,
temperature, and ensemble before comparing results. The centered interaction
above differs by known number and constant terms from `U*n_up*n_down`.

## Benchmarks and tests

Actual external-package runs cross-check **DQMC against SmoQyDQMC.jl 2.0.14**,
**ED against QuSpin 1.0.1**, and **ITensorMPS-based DMRG against TeNPy 1.1.1**.
All comparisons passed on the small Hubbard systems tested. QuSpin also
provides an independent exact finite-step reference for both DQMC engines.
See the [external-package report](docs/cv_benchmarks/external_packages.md)
for numerical differences, tested scope, shared DQMC dependencies and saved
evidence.

Additional references are the **analytic Hubbard dimer solution** of Carrascal
et al. (2015), the **noninteracting Fermi-Dirac limit**, and **full exact
diagonalization** of small square clusters. Details and reproduction commands:

- [ED validation](docs/cv_benchmarks/ed.md)
- [DMRG validation](docs/cv_benchmarks/dmrg.md)
- [DQMC validation and rank-one timings](docs/cv_benchmarks/dqmc.md)
- [Regression and independent full-Fock checks](docs/cv_benchmarks/crosscheck.md)

Run the regression suite and optional DMRG checks:

```powershell
julia --project=. --startup-file=no .\test\runtests.jl
julia --project=. --startup-file=no .\test\test_dmrg_package.jl
```

The earlier 1D Rényi study is described in
[docs/dqmc_renyi_ee.md](docs/dqmc_renyi_ee.md). Historical benchmark snapshots
are separate from these new package checks. New benchmark commands refuse to
overwrite existing outputs.
