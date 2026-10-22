# Exact diagonalization: use and validation

The ED tools solve a small Hubbard lattice without Monte Carlo sampling.
Use `ed_ground_state` for a fixed particle sector and `ed_thermal` for a
grand-canonical thermal reference. Both use the particle-hole symmetric
interaction `U(n_up-1/2)(n_down-1/2) - mu(n_up+n_down)`.

```powershell
julia --project=. .\examples\ed_square.jl
julia --project=. .\examples\benchmark_ed.jl .\results\my_new_ed_check
```

Edit the lattice size, hopping, and boundaries in `examples/ed_square.jl`.
The same real symmetric hopping matrix can also describe a custom graph or
site potentials. The default dense sector limit is 4096 states; thermal ED
defaults to at most six sites. These limits make ED suitable for validating
small systems. DMRG or DQMC handles larger lattices.

## Named reference

**Analytic Hubbard dimer solution**, as presented by Carrascal, Ferrer,
Smith, and Burke, *J. Phys.: Condens. Matter* **27**, 393001 (2015),
[Eq. (36)](https://arxiv.org/html/1502.02194v1#S2.E36).
Its conventional half-filled ground-state energy is
`(U - sqrt(U^2 + 16t^2))/2`. Our centered convention subtracts `U/2`
and the chemical potential subtracts `2mu`. The benchmark separately
constructs the analytic many-particle levels to check thermal observables.
It also checks the noninteracting square lattice against Fermi-Dirac
occupations of the one-particle levels.

## Measured results

Run on October 21, 2026 (America/New_York), using Julia 1.12.1 on Linux
(WSL), an AMD Ryzen 7 7700X, one Julia thread, and one OpenBLAS thread.
Each case passed the absolute tolerance `2e-12`.

| Independent check | Cases | Largest absolute difference |
| --- | ---: | ---: |
| Four dimer eigenvalues in the `(Nup,Ndn)=(1,1)` sector | 24 | `5.33e-15` |
| Dimer thermal energy, density, double occupancy, and log partition function | 48 | `9.77e-15` |
| Free 2×2 square thermal observables | 9 | `4.62e-14` |

The dimer grid used `t = 0.5, 1.0, 1.7`, `U = 0, 2, 4, 8`,
and `mu = -0.3, 0`. Thermal cases used `beta = 0.2, 2.0`.
The free square used `mu = -0.3, 0, 0.2` and `beta = 0.2, 1.3, 4.0`.
The unit tests additionally check dimer double occupancy and spatial
Rényi entropy, infinite-temperature limits, and input errors.

Evidence is in `results/cv_benchmarks_20261021/ed_dimer_ground.csv`,
`ed_dimer_thermal.csv`, `ed_free_square.csv`, and `ed_metadata.toml`.
`ed_square_references.csv` supplies matched interacting references for
DQMC and DMRG on open 2×1, 2×2, and 3×2 lattices. Benchmark scripts
refuse to overwrite these files.

Suggested CV wording: “Developed configurable Hubbard-model exact
diagonalization tools and validated them against the analytic Hubbard
dimer solution and noninteracting square-lattice limits.”

The tested reference is an analytic finite-system solution. These results
do not claim validation against a thermodynamic Bethe-ansatz result.

## Exact finite-step reference for DQMC

`examples/benchmark_ed_trotter.jl` computes the exact many-body transfer
matrix at each imaginary-time step, then sums traces across every particle
sector. This checks the discretized DQMC target as well as the continuous
Hamiltonian. For the open 2×2 square at `U=2`, `mu=0`, `beta=2`:

| `dtau` | Total centered energy | Double occupancy per site |
| --- | ---: | ---: |
| 0 (continuous Hamiltonian) | -4.289791719765 | 0.161384971038 |
| 0.1 | -4.302904959004 | 0.160391485401 |
| 0.05 | -4.293094163105 | 0.161135391813 |
| 0.025 | -4.290618858270 | 0.161322500183 |

These values are saved in `ed_trotter_square.csv`. The transfer matrix uses
`exp(-dtau*K) exp(-dtau*V)`. Reversing that asymmetric ordering preserves
the partition function and the kinetic, potential, and double-occupancy
expectation traces by cyclicity and the relevant commutation relations.
This permits a direct comparison to the DQMC slice measurements. The
finite-step energy and double occupancy approach the continuous values
quadratically in `dtau`.

```powershell
julia --project=. .\examples\benchmark_ed_trotter.jl .\results\my_new_trotter_check.csv
```
