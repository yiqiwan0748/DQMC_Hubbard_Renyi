# U=4 Mott crossover pilot

## Question

The U=2 APBC exact data at L=8 and L=12 give `c_eff = 1.9886`, even though
the asymptotic half-filled repulsive Hubbard chain has only one gapless spin
mode and `c = 1`. Test whether the same small systems move toward `c = 1` when
the interaction is increased to U=4 and the charge correlation length becomes
shorter.

The two-point definition used here is

```text
c_eff = 4 [S2(L=12) - S2(L=8)] / log(12/8).
```

Both sizes use closed-shell APBC and a half-chain subsystem.

## Exact result

| U | S2 ED, L=8 | S2 ED, L=12 | c_eff |
|---:|---:|---:|---:|
| 2 | 1.438139939952 | 1.639720441050 | 1.9886 |
| 4 | 1.145870324242 | 1.258158791329 | 1.1077 |

The exact small-system slope moves from approximately 2 to approximately 1.
This strongly supports the expected interpretation that the U=2 result lies on
the short-distance side of the Mott crossover.

## DQMC screening run

The U=4 screening run used `dtau=0.05`, `beta=4L`, `Ninc=8`, 500 burn-in
sweeps, 1,000 measurements, 20 bins, and four roots. Raw data were written to a
scratch directory outside the repository.

| L | S2 DQMC | S2 ED | difference |
|---:|---:|---:|---:|
| 8 | 1.185789049266 +/- 0.020607329183 | 1.145870324242 | +1.94 sigma |
| 12 | 1.349354013793 +/- 0.031463441265 | 1.258158791329 | +2.90 sigma |

The corresponding DQMC slope is `c_eff = 1.61 +/- 0.37`, which is consistent
with the exact slope but is not precise. All 72 nodes completed, every phase was
exactly 1, adjacent-ratio consistency checks were within 1.64 sigma, and resume
validation accepted every node.

The first five-bin block was unusually high at L=12, and the four individual
L=12 root estimates were more dispersed than at U=2. Therefore these short
DQMC runs are screening data, not accepted benchmarks. A longer-thermalization
or time-step study would be required to isolate equilibration, autocorrelation,
finite-temperature, and Trotter effects. That follow-up was deliberately not
run here.

## Conclusion

The exact two-size contrast strongly supports the Mott-crossover explanation:
accessible U=2 chains look nearly `c=2`, while at U=4 the same sizes already
give `c_eff` close to 1. The short DQMC pilot is qualitatively compatible with
this crossover but is not an accepted quantitative U=4 benchmark.

## Reproduce in PowerShell

```powershell
$scratch = Join-Path $env:TEMP ("HubbardRenyiDQMC_U4_" + [guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Force $scratch | Out-Null

julia .\analysis\compute_s2_ed.jl --L 8 --U 4 --nA 4 --boundary apbc --output (Join-Path $scratch "L8_ED.csv")
julia .\analysis\compute_s2_ed.jl --L 12 --U 4 --nA 6 --boundary apbc --output (Join-Path $scratch "L12_ED.csv")

julia .\run\run_grid_local.jl .\parameters\L8_U4_beta32_nA4_APBC_Ninc8_pilot.toml --all-seeds --max-jobs 4 --output-dir (Join-Path $scratch "L8_raw")
julia .\analysis\summarize_grid.jl (Join-Path $scratch "L8_raw") --parameter-file .\parameters\L8_U4_beta32_nA4_APBC_Ninc8_pilot.toml --output (Join-Path $scratch "L8_summary.csv")

julia .\run\run_grid_local.jl .\parameters\L12_U4_beta48_nA6_APBC_Ninc8_pilot.toml --all-seeds --max-jobs 4 --output-dir (Join-Path $scratch "L12_raw")
julia .\analysis\summarize_grid.jl (Join-Path $scratch "L12_raw") --parameter-file .\parameters\L12_U4_beta48_nA6_APBC_Ninc8_pilot.toml --output (Join-Path $scratch "L12_summary.csv")
```
