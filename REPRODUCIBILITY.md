# Reproduce the package checks

Use Julia 1.12.1 from this repository directory. All commands below work in
PowerShell. The benchmark commands refuse existing output paths; choose a
fresh path for each run.

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --startup-file=no .\test\runtests.jl
julia --project=. --startup-file=no .\test\test_dmrg_package.jl
julia --project=. --startup-file=no .\analysis\validate_grover_rank1_update.jl
```

Run the exact references, thermal DQMC comparison, and DMRG comparison:

```powershell
julia --project=. .\examples\benchmark_ed.jl .\results\new_ed
julia --project=. .\examples\benchmark_ed_trotter.jl .\results\new_trotter.csv
julia --project=. .\analysis\cv_benchmark_dqmc.jl thermal .\results\new_thermal
julia --project=. .\examples\dmrg_benchmark.jl .\results\new_dmrg.csv
```

For the speed comparison, keep other computational jobs idle:

```powershell
julia --project=. .\analysis\cv_benchmark_dqmc.jl timing .\results\new_timing
```

Reports are in `docs/cv_benchmarks/`; the measured CSVs and settings are in
`results/cv_benchmarks_20261021/`. See the README for the editable square
examples and the shared Hamiltonian convention.

The retained 1D Renyi runner also supports a new temporary output directory:

```powershell
$raw = Join-Path $env:TEMP ('HubbardRenyiDQMC_' + [guid]::NewGuid().ToString('N'))
julia --project=. .\run\run_grid_local.jl .\parameters\L10_U2_beta40_Ninc8.toml --all-seeds --max-jobs 4 --output-dir $raw
julia --project=. .\analysis\summarize_grid.jl $raw --parameter-file .\parameters\L10_U2_beta40_Ninc8.toml --output (Join-Path $raw 'summary.csv')
```

This last command pair is a larger Monte Carlo job. The older study notes
record earlier results; the research archive, accepted raw data, and old
processed snapshots are not part of this source release.
