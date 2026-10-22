# Square-lattice DQMC and rank-one benchmark

Run date: October 21, 2026. New evidence is in
[`results/cv_benchmarks_20261021`](../../results/cv_benchmarks_20261021/).

## What the package provides

`run_square_dqmc` measures finite-temperature energy, density and double occupancy
for a rectangular square lattice. Edit `parameters/square_4x4.toml` to change its
size, hopping, interaction, temperature, boundaries or sample count. The existing
chain Rényi workflow remains available through `run/run_lambda_node.jl`.

The square-lattice entry point uses the same JDQMCFramework propagators and
stabilized update routines as the existing code. LatticeUtilities defines the
unit cell, finite lattice and nearest-neighbor bonds. This follows the
[JDQMCFramework square Hubbard tutorial](https://smoqysuite.github.io/JDQMCFramework.jl/stable/examples/square_hubbard/)
and the Hamiltonian convention in the
[SmoQyDQMC square Hubbard tutorial](https://smoqysuite.github.io/SmoQyDQMC.jl/stable/tutorials/hubbard_square/).
This package uses those framework components; it is not the SmoQyDQMC package.

The shared Hamiltonian is

```text
H = sum_ij,spin K[i,j] cdag[i,spin] c[j,spin]
    + U sum_i (nup[i]-1/2)(ndn[i]-1/2) - mu sum_i (nup[i]+ndn[i]).
```

`mu=0` is half filling on the even nearest-neighbor square lattices used here.
Sites are numbered with x increasing first: `site=x+(y-1)*Lx`.
Each undirected bond occurs once. A periodic axis of length two therefore has
one bond between its two sites, rather than two parallel bonds. The 2-by-2
benchmark has four bonds; its open and periodic hopping matrices coincide.
`initializer_from_kinetic` accepts a custom finite real symmetric hopping
matrix for researchers who need another geometry.

From PowerShell:

```powershell
julia --project=. --startup-file=no --history-file=no .\run\run_square_hubbard.jl .\parameters\square_4x4.toml
```

The command saves a summary, signed measurement bins and the input parameters.
It refuses to replace an existing output directory; use `--output-dir` to choose
a new one. Energy is the total centered Hamiltonian energy, including `-mu*N`;
the summary also gives energy per site. Density and double occupancy are per
site. Errors are bin jackknife standard errors. Increase burn-in, measurements
and bin size for a new research problem, and compare smaller `dtau` values.
Low average sign makes estimates harder to converge. The default half-filled,
even square lattice has no sign problem.

## Matched finite-temperature ED check

The reference is full-spectrum exact diagonalization, summed over all spin and
particle-number sectors. Both methods use the same 2-by-2 four-bond matrix,
`t=1`, `U=2`, `mu=0`, `beta=2` and the grand-canonical ensemble. This comparison
does not substitute a ground-state or fixed-sector result for thermal ED.

Four independent root seeds were used for each time step. Each stream discarded
1,000 sweeps, then measured 20,000 configurations in 40 bins of 500 sweeps.
The table pools the 160 bins per time step. Average sign was one.

| Observable | Full-spectrum ED | DQMC, dtau=0.10 | DQMC, dtau=0.05 |
|---|---:|---:|---:|
| Total energy | -4.2897917198 | -4.3020885651 +/- 0.0029684075 | -4.2958807288 +/- 0.0033254466 |
| Density per site | 1.0000000000 | 1.0000000000 | 1.0000000000 |
| Double occupancy per site | 0.1613849710 | 0.1604164936 +/- 0.0001560272 | 0.1613043445 +/- 0.0001526590 |

The coarse time step has visible discretization bias. At `dtau=0.05`, energy
differs from ED by 1.83 sampling standard errors and double occupancy by 0.53.
Sampling errors exclude time-step bias. This is a small-system validation,
not a claim of a converged continuum-time result for every model.

An additional independent exact calculation evaluates the same finite-step
transfer matrix, `rho=(exp(-dtau*K_manybody)*exp(-dtau*V_manybody))^(beta/dtau)`,
by dense many-body matrix exponentiation and traces over every particle sector.
This avoids sampling and separates discretization error from Monte Carlo error.

| Time step | Exact finite-step energy | DQMC minus exact, in sampling errors | Exact finite-step double occupancy | DQMC minus exact, in sampling errors |
|---|---:|---:|---:|---:|
| 0.10 | -4.3029049590 | +0.28 | 0.1603914854 | +0.16 |
| 0.05 | -4.2930941631 | -0.84 | 0.1611353918 | +1.11 |

The exact finite-step results approach continuous-time ED when the time step is
halved. The Monte Carlo data agree with their matched finite-step references.
See `examples/benchmark_ed_trotter.jl` and `ed_trotter_square.csv` for all values,
including an exact `dtau=0.025` reference.

The separate 4-by-4 example completed with 500 burn-in and 4,000 measurements:
energy per site `-1.81945048 +/- 0.00381879`, density `1`, double occupancy
`0.12337615 +/- 0.00046678`, average sign `1`. It demonstrates the editable
runner; it is not an ED benchmark.

Reproduce in a new output directory:

```powershell
julia --project=. --startup-file=no --history-file=no .\analysis\cv_benchmark_dqmc.jl thermal .\results\cv_benchmarks_new
```

## Rank-one timing protocol

The baseline rebuilds and factors the Grover matrix for each proposed change.
The optimized path updates its determinant and inverse using the determinant
lemma and Sherman-Morrison formula. Both paths already use rank-one Green
function updates, so the comparison measures the additional Grover optimization.

For an active subsystem Green block `G1` and passive block `G2`,

```text
g = G1*G2 + (I-G1)*(I-G2)
G1_new = G1 + p*q^T
g_new = g + p*r^T,  r^T = q^T*(2G2-I)
det(g_new)/det(g) = 1 + r^T*inv(g)*p.
```

This preserves the sampled lambda weight. The existing guard and explicit
fallback paths remain active. All 2,318 existing identity, tracker, fallback
and matched-seed checks passed again. The new square tests also compare the
noninteracting thermal result with the Fermi-Dirac formula.

Timing uses Julia 1.12.1 on an AMD Ryzen 7 7700X under WSL/Linux, with one Julia
thread and one OpenBLAS thread. Each of seven repeats uses matched seeds, two
untimed warm-up rounds, and a garbage collection before 100 timed two-replica
update rounds. Direct and optimized execution order alternates. Initialization,
file output and observable measurements are excluded. A separate three-round
check compares fields, all Green functions, determinant metadata and RNG stream
positions for each seed pair. The source hashes, timings, allocations, GC time
and maximum stabilization corrections are saved.

| Workload | Direct mean time | Rank-one mean time | Ratio of mean times | Direct / rank-one allocations |
|---|---:|---:|---:|---:|
| 4-by-4 square, U=4, beta=4, dtau=0.1, nA=8 | 0.787904 s | 0.171716 s | **4.5884x** | 2244.01 / 11.18 MiB |
| 8-site chain, U=2, beta=32, dtau=0.05, nA=4 | 3.037387 s | 1.187841 s | **2.5571x** | 6394.61 / 56.95 MiB |

Both workloads use `lambda=0.5` and stabilization every ten time slices.
The square uses `C=1:8` and the chain `C=1:4`. Reported times cover 100
two-replica rounds. Square paired-repeat ratios range from 3.06x to 5.19x,
with median 4.81x; chain ratios range from 2.43x to 2.78x, with median 2.50x.
Garbage collection inside the timed region is included. The optimized path
also removes allocation-heavy direct rebuilds, so that saving is part of the
measured improvement. Maximum corrected Green-function discrepancies were
`2.00e-7` for the square and `5.08e-11` for the chain. The three-round matched
checks had identical final Green functions and auxiliary fields for every seed.

A concise CV statement supported by this run is:

> Accelerated square-lattice Rényi DQMC update kernels by 4.6x using rank-one
> Grover determinant updates, compared with direct rebuilds on a 4-by-4 Hubbard
> benchmark; validated thermal observables against exact diagonalization.

This is a measured kernel comparison on the stated hardware and workload.
It does not measure total application runtime or compare against SmoQyDQMC.

```powershell
julia --project=. --startup-file=no --history-file=no .\analysis\cv_benchmark_dqmc.jl timing .\results\cv_benchmarks_new
```

Keep other computational jobs idle while collecting performance numbers.
See `docs/grover_rank1_speedup.md` for the older timing record; its historical
3.8236x comparison involved separate production runs and is distinct from this
paired kernel benchmark.
