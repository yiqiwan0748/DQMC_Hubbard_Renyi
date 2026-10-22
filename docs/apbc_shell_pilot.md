# APBC closed-shell pilot

## Question

For half-filled chains, use PBC for `L = 4m + 2` and APBC for `L = 4m` so the
noninteracting Fermi sea is closed shell. Does this remove the strong finite-size
branching in the half-chain Renyi entropy?

APBC changes only the wrap hopping `L <-> 1` from `-t` to `+t`. The lattice is
still a ring, and the Renyi estimator is unchanged.

## Screening setup

- `U = 2`, `mu = 0`, `dtau = 0.05`, `beta = 4L`
- half-chain subsystem `nA = L/2`
- `Ninc = 8`, 500 burn-in sweeps, 1,000 measurements, 20 bins
- roots `1234`, `11234`, `21234`, `31234`
- 36 nodes per size

## Result

| L | boundary | exact ED S2 | DQMC S2 | DQMC - ED |
|---:|:---:|---:|---:|---:|
| 8 | PBC | 1.309407045626 | 1.318734225375 +/- 0.007790202478 | +0.009327179749 |
| 8 | APBC | 1.438139939952 | 1.446183297421 +/- 0.008290348556 | +0.008043357469 |
| 12 | PBC | 1.542354369358 | 1.561005400404 +/- 0.007601610412 | +0.018651031046 |
| 12 | APBC | 1.639720441050 | 1.657932218556 +/- 0.010169092474 | +0.018211777506 |

The APBC DQMC values agree with ED within `0.97 sigma` at `L=8` and `1.79 sigma`
at `L=12`. An independent L8 APBC DMRG calculation gave
`S2 = 1.438139939904`, agreeing with ED.

For a two-point fit `S2 = (c_eff/4) log(L) + constant` using `L=8,12`:

- PBC ED: `c_eff = 2.2981`
- APBC ED: `c_eff = 1.9886`
- PBC DQMC: `c_eff = 2.3901`
- APBC DQMC: `c_eff = 2.0889`

## Conclusion

Closed-shell APBC moves both even-`nA` points upward and reduces the excessive
logarithmic slope, so shell filling is a real part of the finite-size effect. It
does not recover the expected `c = 1` scaling at these sizes. Strong subsystem
parity and other subleading corrections remain. A high-stat APBC production run
is therefore not justified by this pilot.

All 72 APBC nodes completed without failure. Every raw bin had phase exactly 1,
all adjacent-ratio consistency checks were within `1.8 sigma`, and a full resume
pass validated and skipped every node. Raw pilot data were kept outside the
repository. A same-seed PBC regression node reproduced its accepted raw bins
byte-for-byte, confirming that the default PBC path is unchanged.

## Reproduce in PowerShell

```powershell
$scratch = Join-Path $env:TEMP "HubbardRenyiDQMC_APBC"
New-Item -ItemType Directory -Force $scratch | Out-Null

julia .\run\run_grid_local.jl .\parameters\L8_U2_beta32_nA4_APBC_Ninc8_pilot.toml --all-seeds --max-jobs 4 --output-dir (Join-Path $scratch "L8_raw")
julia .\analysis\summarize_grid.jl (Join-Path $scratch "L8_raw") --parameter-file .\parameters\L8_U2_beta32_nA4_APBC_Ninc8_pilot.toml --output (Join-Path $scratch "L8_summary.csv")

julia .\run\run_grid_local.jl .\parameters\L12_U2_beta48_nA6_APBC_Ninc8_pilot.toml --all-seeds --max-jobs 4 --output-dir (Join-Path $scratch "L12_raw")
julia .\analysis\summarize_grid.jl (Join-Path $scratch "L12_raw") --parameter-file .\parameters\L12_U2_beta48_nA6_APBC_Ninc8_pilot.toml --output (Join-Path $scratch "L12_summary.csv")
```
