# One-dimensional Hubbard U/t=2 FSS extension through L=22

> Historical validation record. See [`dqmc_renyi_ee.md`](dqmc_renyi_ee.md)
> for the current scientific summary and conclusions.

> Source-release note: the historical data, figure scripts, and separate FSS
> environment referenced below are omitted. Reproducing those artifacts requires
> the original research workspace. Included package checks are listed in
> [REPRODUCIBILITY.md](../REPRODUCIBILITY.md).

## Result

This note extends the accepted finite-size dataset with validated DQMC points
at L=20 and L=22. The earlier 20-observation input, 273-fit summary, figures,
and detailed diagnostic note remain unchanged.

The new points agree with the corresponding ground-state DMRG references and
make the persistent `L mod 4` oscillation clearer. They do **not** make the
simple finite-window fits adequate, so no asymptotic central-charge value is
claimed. For the repulsive half-filled chain the expected asymptote has one
gapless sector (`c=1`), but U/t=2 and L<=22 can remain in a long finite-size
crossover. The fitted `c_eff_diagnostic` values below are not evidence for a
different asymptotic central charge.

## New validated DQMC points

Both production grids completed all 68 lambda-node/root pairs and all 2,720
bin rows with zero failed nodes. The displayed errors are the accepted
production absolute DQMC one-sigma estimates; the canonical fits retain these
values and treat cross-size errors as independent. The targeted post-production
audit below gives a separate conservative interpretation without changing the
fit input.

| L | nA | beta | S2 DQMC | DQMC 1 sigma | S2 DMRG | DQMC-DMRG | difference / DQMC sigma |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 10 | 80 | 1.805179361272 | 0.008209270540 | 1.808766655046 | -0.003587293774 | -0.437 |
| 22 | 11 | 88 | 2.025158749526 | 0.004730458606 | 2.021914792180 | +0.003243957346 | +0.686 |

Across all five common sizes L=14,16,18,20,22, the descriptive sum of squared
DQMC-only z-scores is 4.4033. The inverse-DQMC-variance weighted offset is
`DQMC-DMRG = 0.0015832 +/- 0.0022139`, or 0.715 nominal combined standard
errors under the diagonal-covariance assumption. This is a consistency check,
not a formal combined-method goodness-of-fit test: the DMRG rows have no quoted
uncertainty and cross-size DQMC covariance is unavailable.

## Provenance

Production used commit
`c4c873ffa16908a6641f851cb6c9ad41a9a22c29`. The exact promoted summary bytes
were committed separately at
`c766667ee7ee327bbf0744bfacbe01b82fd4cea0` before constructing the extended
fit input.

| L | promoted snapshot | SHA-256 | Git blob OID | combined raw-manifest SHA-256 |
|---:|---|---|---|---|
| 20 | `results/processed/FSS_1D_U2_beta4L_L20_Ninc16.csv` | `d6c734257a969f6ccc465972fc7c079fa2772b0c9d9794820538a902db6b048b` | `47cec81f8eae2d5196113b23f54b270671383d7f` | `1d9e9fd097352a9cccf404a43d973abdca1aadc40d618dd5727c9024708d1609` |
| 22 | `results/processed/FSS_1D_U2_beta4L_L22_Ninc16.csv` | `bf82b619ae420912ec3795731b34c2791f3e4c37258a9195cadb8179cbe3dafd` | `9cfd7a90cb7962e970937e0f0019f3a38dc7e14b` | `368ec6bcce670690c3105192c4afb58926a4c38a377378611335b83926c09eca` |

The combined raw-manifest hashes cover the 68 bin files and 68 node-summary
files for each size. They are SHA-256 hashes of 136 LF-terminated GNU
`sha256sum` records, ordered lexicographically by basename and serialized as
`<file SHA-256><two spaces>./<basename>`. Raw production data remain in the
external run directory and were not copied into or modified inside the
repository.

The production-code hashes were:

| Role | SHA-256 |
|---|---|
| `run/run_lambda_node.jl` | `dae2051cfee2a443084890f40ee35cab628eb4b3393e6c7c427b4c82711462e6` |
| `run/grover_rank1_update.jl` | `443f3be1b6d39ad06909f20a9c9a54834625dcf5b6e79d84a2f9389484fe168c` |
| `analysis/summarize_grid.jl` | `6319928e853d0840cc4a0c4c334122c829b7ab0ab9f277a25ce9cac5ebc55ce3` |

The L20/L22 parameters use `dtau=0.05`, `Ninc=16`, four root seeds,
`N_burnin=500`, 4,000 measurements, and 40 bins per node.

## Targeted long-burn validation

Targeted checks raised `N_burnin` from 500 to 8,000 while holding `dtau`,
`Ninc`, the 4,000-measurement count, 40-bin partition, estimator, and four root
seeds fixed. Their outputs were isolated from the accepted production grids.
L20 tested three nodes selected from the raw audit's localized correlation and
root-agreement watches (`k=6,11,13`); L22 tested `k=14`, its largest production
early/late node discrepancy. The corrected long-burn launches completed all 16
node/root jobs with zero failures.

Old-versus-check comparisons used each estimate's larger bin-jackknife or
within-root contiguous block-10 jackknife error. The shift diagnostic replaces
only the tested forward entropy increments; it does not define a new full-grid
central value.

| L | tested nodes | completed | largest old/check difference | tested-forward shift diagnostic | max early/late | max block-10 inflation | result |
|---:|---|---:|---:|---:|---:|---:|---|
| 20 | `k=6,11,13` | 12/12 | 1.142 sigma | +0.002995418018 (0.485 sigma) | 2.564 sigma | 1.271 | PASS |
| 22 | `k=14` | 4/4 | 0.455 sigma | -0.000763112398 (0.161 full-result sigma) | 1.157 sigma | 0.940 | PASS |

These checks cover 3 of 16 L20 forward nodes and 1 of 16 L22 forward nodes.
They support no detected burn-in bias at the audited nodes, not a global
equilibration proof.

For L20, all six plus/minus comparisons passed and the minimum root-agreement
`p` value was 0.5535. For L22, the separate plus/minus differences were 0.455
and 0.403 standard errors, and the minimum root-agreement `p` value was 0.6022.
Both validation corpora passed exact coverage, schema, finite-value, and
bin-to-summary reproduction checks.

| L | validation parameter file | parameter SHA-256 | validation raw-CSV manifest SHA-256 |
|---:|---|---|---|
| 20 | `parameters/L20_U2_beta80_nA10_Ninc16_burnin8000.toml` | `06c42ecb5f178221fa7448827a865de132ee93d5768838f3a8c9407c5c27152c` | `f4f6c068b487348a319ec77fda548e38bd96aa130076c4acf4fbfad2bc39efa3` |
| 22 | `parameters/L22_U2_beta88_nA11_Ninc16_burnin8000.toml` | `7da5c1fc525b4ca7c029b1ba3ceb64ffabf3563e3ca74bd3d2e4a22e321cac59` | `e3986effae133e8bd6f41d372c7dc961a95bd193d598a7f7c04678759cbe7269` |

The validation manifests cover 24 L20 CSVs (12 bin files and 12 node
summaries) and 8 L22 CSVs (4 bin files and 4 node summaries), using the same
sorted GNU `sha256sum`-record serialization defined above. Logs are excluded.
The raw checks remain external under
`C:\Users\YQW\HubbardRenyiDQMC_runs\fss_L20_burnin8000_check_20260821\raw`
and
`C:\Users\YQW\HubbardRenyiDQMC_runs\fss_L22_burnin8000_check_20260821\raw`.

The full L20 raw audit found a root-cluster jackknife error of 0.0094715 and a
block-10 error of 0.0092614, compared with the production 0.0082093. Therefore
the conservative standalone report is `S2(L20)=1.80518 +/- 0.0095`; the
accepted snapshot and canonical fits retain their nominal
`sigma=0.008209270540`. L22 blocking was stable, so its standalone report
remains `S2(L22)=2.025159 +/- 0.004730`. No full-grid rerun is warranted by
these diagnostics.

## Extended fits

The new canonical input contains 10 DQMC, 5 ED, and 7 DMRG observations. It
keeps the previously approved Lmin choices unchanged. Only the automatic
contiguous-window and leave-one-out suites grow when L20/L22 are added.

| Fit role | rows |
|---|---:|
| principal | 15 |
| approved Lmin sensitivity | 58 |
| contiguous windows | 84 |
| leave one out | 132 |
| total | 289 |

The five principal DQMC diagnostics are:

| Model | branch / p | c_eff diagnostic | nominal WLS SE | chi2/dof |
|---|---|---:|---:|---:|
| total | all | 1.490492211 | 0.011227150 | 477.892675 |
| branch | L mod 4 = 0 | 1.811460749 | 0.028552501 | 29.278026 |
| branch | L mod 4 = 2 | 1.504059961 | 0.012270437 | 3.899204 |
| oscillation | p=1 | 1.791547725 | 0.012403920 | 80.555709 |
| oscillation | p=2 | 1.754515144 | 0.012610513 | 244.178174 |

All five full-window `chi2/dof` values are too large for a precision
interpretation. The nominal WLS errors only propagate the reported diagonal
Monte Carlo errors through a chosen regression; they do not include model,
thermal, Trotter, equilibration, or cross-size covariance uncertainty.

The DQMC total-fit sensitivity ranges are 1.41595-1.79060 for the approved
Lmin suite, 1.32269-1.85677 for contiguous five-size windows, and
1.41595-1.92449 for full-fit leave-one-out checks. Across all sensitivity
roles, branch fits span 1.45699-2.24589 and oscillatory fits span
1.58979-2.51667. The spread is much larger than any principal nominal SE.

## Artifacts and validation

| Artifact | SHA-256 |
|---|---|
| `results/processed/oneD_hubbard_fss_fit_input_L4_to_L22_20260820.csv` | `755be72882eaf42dd4e382446fb0e56be5e5885a309040b7a3d23a68f8b6eeee` |
| `results/processed/oneD_hubbard_fss_fit_summary_L4_to_L22_20260820.csv` | `09381d05b0155e1ba28b2b6181d28adbb0090f744f11daddb7b4e79c5127738a` |
| `results/figures/oneD_hubbard_fss_L4_to_L22_20260820/` | four PNG/SVG pairs |

The summary has 289 data rows. Every fit has positive degrees of freedom and
full SVD rank. There are no conditioning warnings; the maximum normalized
condition number is 31.0218, below the hard limit of `1e6`. The largest
QR-versus-SVD coefficient difference is `1.38556e-13`, within the fitter's
absolute-plus-relative tolerance. A clean input rebuild matched the committed
input byte for byte. An independent repeated fit reproduced the summary and
all eight figures byte for byte. Visual inspection found the four PNG figures
complete and unclipped; SVG counterparts were also generated.

The final artifacts were generated with already-installed WSL Julia 1.12.1
and Plots 1.41.6; no package operation was performed. An independent earlier
Windows Julia 1.12.1 / Plots 1.41.2 run agreed at every reported digit and
passed the same gates. The checked-in FSS project pins Plots 1.41.2 as the
locked reproduction target. The command below assumes that environment is
already instantiated and performs no package operation; raw serialized floats
and graphics bytes can differ across BLAS and rendering platforms.

```powershell
$run = Join-Path $env:TEMP ('HubbardRenyiDQMC_fss_L22_' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $run | Out-Null
$input = Join-Path $run 'fit_input.csv'
$figures = Join-Path $run 'figures'
$summary = Join-Path $run 'fit_summary.csv'
julia --startup-file=no --history-file=no .\analysis\build_oneD_hubbard_fss_L4_to_L22_input.jl --output $input --no-clobber
julia --startup-file=no --history-file=no --project=.\environments\fss .\analysis\fit_oneD_hubbard_fss_with_oscillation.jl --input $input --summary-output $summary --figure-dir $figures --no-clobber
```

## Limits and next step

- DQMC is finite-beta and grand-canonical; DMRG is a fixed-sector ground-state
  reference.
- L20/L22 DMRG values are first-pass, maxlinkdim-limited references with no
  reported uncertainty.
- Targeted `N_burnin=8000` checks passed at L20 `k=6,11,13` and L22 `k=14`,
  but they are not full-grid long-burn convergence tests, do not directly
  observe the original pre-measurement trajectories, and do not cover untested
  nodes.
- L20/L22 do not yet have independent beta or `dtau` convergence checks.
- Agreement with DMRG is reassuring but cannot separately bound thermal,
  Trotter, or equilibration bias.

The present FSS milestone is therefore: validated DQMC through L=22,
cross-method consistency at large L, and confirmed mod-4 shell oscillations.
The next scientific step is beta/`dtau` convergence control and/or larger sizes
with a theory-motivated correction model, not selecting one of the current
finite-window slopes as the asymptotic answer.
