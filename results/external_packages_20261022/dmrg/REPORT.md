# DMRG comparison with TeNPy

The published ITensorMPS-based Hubbard DMRG wrapper agrees with an independent
TeNPy 1.1.1 optimization on four small open systems. All 32 cross-package
observable checks passed. The largest differences were **1.25e-14 t in total
energy** and **3.20e-13 in second Renyi entropy**.

The public source tested was commit
`6f8b949275a648f9592410085fd8263b1286f2ba`, with Julia 1.12.1,
ITensorMPS 0.4.1 and ITensors 0.9.30. The independent reference used Python
3.12.14, TeNPy 1.1.1, NumPy 2.5.4 and SciPy 1.18.1. TeNPy built its own stock
`FermiHubbardModel` and optimized its own MPS with `TwoSiteDMRGEngine`; it did
not consume a Julia Hamiltonian, MPO, wavefunction or ED eigenvector.

## Matched systems

All cases have nearest-neighbor hopping t=1 and open boundaries on both axes.
The subsystem contains the first `cut` sites in the public x-fast site order.

| System | U/t | mu/t | Nup, Ndn | cut | Public energy/t, tight run | Public S2, tight run |
|---|---:|---:|---:|---:|---:|---:|
| 2x2 square | 4 | 0 | 2, 2 | 2 | -6.102748483462076 | 0.752122036001990 |
| 3x2 rectangle | 4 | 0 | 3, 3 | 3 | -9.619321323957557 | 1.237787082849253 |
| 3x2 rectangle, one hole | 4 | 0.35 | 3, 2 | 3 | -10.257817706153952 | 1.902105277138321 |
| 6x1 chain | 2 | 0 | 3, 3 | 3 | -7.546313794436001 | 1.230607720070336 |

TeNPy's default square-lattice order advances its y coordinate fastest.
The reference transposes lattice axes, `(Lx, Ly)_TeNPy = (Ly, Lx)_public`,
so every MPS site and the prefix subsystem match exactly. The script asserts
the full coordinate order and all undirected nearest-neighbor bonds before
each optimization. No snake ordering is used.

The public package uses the centered interaction

`H = hopping + U*sum((nup-1/2)*(ndn-1/2)) - mu*N`.

TeNPy's stock model uses `U*nup*ndn - mu_TeNPy*N`. The comparison sets
`mu_TeNPy = mu + U/2` and adds the constant `Nsites*U/4` to its total energy.
This is exactly the identity
`U*(nup-1/2)*(ndn-1/2) = U*nup*ndn - U*(nup+ndn)/2 + U/4`.
It preserves the wavefunction, entanglement and measured occupations.

## Accuracy and convergence

Each package ran two independent initial-state settings. The standard setting
used a bond-dimension schedule 16, 32, 64, 128 and cutoff 1e-13; the tight
setting used 32, 64, 128, 256 and cutoff 1e-14. The public wrapper performed
10 and 14 sweeps. TeNPy performed 11 and 15 sweeps under its stopping-counter
convention; actual sweep counts are recorded in the CSV. Initial-state seed
labels were 1234 and 4321, with different particle arrangements in TeNPy.

| Observable | Largest absolute package difference | Acceptance tolerance |
|---|---:|---:|
| Total energy/t | 1.24345e-14 | 1e-8 |
| Prefix S2 | 3.19745e-13 | 1e-7 |
| Mean density | 4.44090e-16 | 1e-10 |
| Mean double occupancy | 7.88259e-15 | 1e-8 |

Across the two settings, the largest public energy and S2 changes were
1.07e-14 t and 6.41e-14. The corresponding TeNPy changes were 1.96e-14 t and
1.56e-15. The final bond dimensions were 16 for the square, 64 for the
half-filled rectangle and chain, and 36 for the one-hole rectangle.

As an additional check, each public result was compared with the package's
full sector ED. The largest energy difference was 1.78e-14 t. All four sector
ground states were nondegenerate; the smallest measured gap was 0.2955 t.
The independent package comparison is TeNPy, rather than this internal ED
check. These finite-cluster tests do not establish convergence for larger
two-dimensional systems or periodic boundaries.

## Evidence and reproduction

- `public.csv`: all eight public wrapper runs and supporting ED values.
- `tenpy/tenpy.csv`: all eight independent TeNPy optimizations.
- `tenpy/comparison.csv`: every measured difference, tolerance and pass flag.
- `tenpy/environment.json`: exact Python dependency versions and conventions.
- `julia_environment_verified.txt`: Julia and ITensor versions.
- `convergence_summary.json`: changes between standard and tight settings.
- `public.log` and `tenpy.log`: captured execution logs.

From the published repository root, use a new output directory on each repeat.
These commands are PowerShell compatible:

```powershell
julia --project=. .\analysis\external_packages\dmrg_public.jl .\results\external_packages_new\dmrg\public.csv
python -m venv .\output\external_packages_new\tenpy_venv
& .\output\external_packages_new\tenpy_venv\Scripts\python.exe -m pip install physics-tenpy==1.1.1 numpy==2.5.4 scipy==1.18.1
& .\output\external_packages_new\tenpy_venv\Scripts\python.exe .\analysis\external_packages\dmrg_tenpy.py --public .\results\external_packages_new\dmrg\public.csv --output-dir .\results\external_packages_new\dmrg\tenpy
```

The DMRG benchmark input is `analysis/external_packages/dmrg_cases.toml`.
The scripts refuse to replace existing result files. Julia can also be run
through WSL with the corresponding Linux paths.

Suggested CV wording: **"Cross-validated ITensorMPS-based Hubbard DMRG against
TeNPy on small open square-lattice clusters."**

## Reference package

The [official TeNPy Hubbard-model documentation](https://tenpy.readthedocs.io/en/latest/reference/tenpy.models.hubbard.FermiHubbardModel.html)
defines its uncentered interaction and chemical potential. The independent
optimizer is the documented
[TwoSiteDMRGEngine](https://tenpy.readthedocs.io/en/latest/reference/tenpy.algorithms.dmrg.TwoSiteDMRGEngine.html).
The package publication is
[Hauschild et al., *Tensor Network Python (TeNPy) version 1*, arXiv:2408.02010](https://arxiv.org/abs/2408.02010).
