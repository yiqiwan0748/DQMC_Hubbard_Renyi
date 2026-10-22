# L10 U/t=2 beta=40 block-size benchmark

> Historical validation record. See [`dqmc_renyi_ee.md`](dqmc_renyi_ee.md)
> for the current scientific summary and conclusions.

## Scope and interpretation

This benchmark compares the Rényi-2 entropy of contiguous subsystems `nA=1:5`
in the periodic, half-filled, one-dimensional Hubbard model at `L=10`, `t=1`,
and `U/t=2`.

The DQMC values use `mu=0`, finite inverse temperature `beta=40`, and a
lambda grid with `Ninc=16`. Their reported errors are one Monte Carlo standard
error. The ED values are ground-state results in the fixed half-filled sector
`Nup=Ndn=5`. Finite-beta DQMC and ground-state ED therefore do not represent
identical state ensembles. Agreement within the Monte Carlo errors is a
benchmark of the finite-beta calculation, not proof that the two ensembles are
identical.

## Benchmark values

The plotting input is the tracked block-size summary at commit
`88c21ae3345ad7414d86b37e777e10c35e6451f2`.

| nA | S2 DQMC | MC 1 sigma | ground-state ED | DQMC - ED | z |
|---:|---:|---:|---:|---:|---:|
| 1 | 1.308138170759 | 0.001209659236 | 1.308461130863 | -0.000322960104 | -0.267 |
| 2 | 1.274574334430 | 0.002188963049 | 1.273890131213 | +0.000684203217 | +0.313 |
| 3 | 1.634894520051 | 0.002334522219 | 1.633561808649 | +0.001332711401 | +0.571 |
| 4 | 1.524569914598 | 0.003041891872 | 1.524217285675 | +0.000352628924 | +0.116 |
| 5 | 1.716441955736 | 0.003040377094 | 1.716551661433 | -0.000109705697 | -0.036 |

Here `z=(S2_DQMC-S2_ED)/S2_error`. Recalculation from the displayed source
values gives a largest absolute value of `0.5708711578` at `nA=3`; all five
DQMC points are within `0.6` Monte Carlo standard errors of ED.

## Half-chain ED precision

The repository-authoritative half-chain ED checkpoint is
`S2=1.716551661433`, as displayed in the block-size summary. The tracked
detailed ED table records `1.716551661443483`; the difference is
`+1.0483e-11`. This extra precision is retained as a reference check and does
not replace the authoritative checkpoint. The difference is negligible
relative to the DQMC half-chain error `0.003040377094`.

The parameter file separately stores `1.7165516614333372`. These are distinct
serialization precisions of the same repository benchmark and are not treated
as independent measurements.

## Provenance

All paths below are repository-relative and all hashes are SHA-256 values at
commit `88c21ae3345ad7414d86b37e777e10c35e6451f2`.

| Role | Path | SHA-256 |
|---|---|---|
| Authoritative plotted values | `results/processed/BLOCKSIZE_U2_beta40_Ninc16_summary.csv` | `1f799d83d74b6ad31f21122f5bb5a40935846e449603096aa2e08f8e987fc3db` |
| Higher-precision ED check | `results/processed/ED_L10_U2_beta40_nA1to5.csv` | `ce957e089fdabeff940230cc9d84b40ddda2b3adb99a8a27e39eb1c76dcd01de` |
| Half-chain run parameters | `parameters/L10_U2_beta40_Ninc16_highstat.toml` | `15e6baaf36109add9787145ec74c2888e7a54f4fa12f7e2fc69a8e86b513a338` |
| Active lambda-node runner | `run/run_lambda_node.jl` | `04b26dc03926ff8b17bb9aabccfcf6b79446ab9e867ada1e04737bc24464e2b6` |
| Validated periodic DQMC kernel | `archive/validated_old_pipeline/run_lambda_grid_step_validated.jl` | `f0868a2d015939bbf9569abe4a45a6e5717d7de446e5402335aaba9882b26a51` |
| Periodic ground-state ED implementation | `analysis/compute_s2_ed.jl` | `a8969d398fbef4144835e9d3cf730606d4d0302d21a93aa51e47c9faff266a87` |

The active runner includes the validated kernel, whose one-dimensional lattice
is explicitly periodic. The ED implementation constructs the periodic
nearest-neighbor chain and the detailed ED table records `Nup=Ndn=5`.

## Deterministic figure generation

The commands in this section require the original research scripts and
processed snapshot, which are omitted from this source release. For checks
that run with the included files, see [REPRODUCIBILITY.md](../REPRODUCIBILITY.md).

The SVG generator uses Julia Base and `Printf` only. It reads the tracked
summary, refuses to overwrite an existing path, and emits no timestamp,
hostname, source path, or package-environment metadata.

From the repository root in PowerShell:

```powershell
$env:JULIA_LOAD_PATH = "@stdlib"
$out = Join-Path $env:TEMP ("L10_U2_beta40_blocksize_" + [guid]::NewGuid().ToString("N") + ".svg")
julia --startup-file=no --history-file=no --compiled-modules=no .\analysis\plot_L10_U2_beta40_blocksize.jl --input-summary .\results\processed\BLOCKSIZE_U2_beta40_Ninc16_summary.csv --output $out
Get-FileHash -Algorithm SHA256 $out
```

The command creates only the new path stored in `$out`; it performs no Julia
package operation and does not read or modify raw simulation data.
