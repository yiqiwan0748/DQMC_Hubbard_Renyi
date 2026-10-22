# One-dimensional Hubbard U/t=2 finite-size diagnostics

> Historical validation record. See [`dqmc_renyi_ee.md`](dqmc_renyi_ee.md)
> for the current scientific summary and conclusions.

## Scope and scientific status

This report records finite-size diagnostic regressions for the half-chain
second Rényi entropy of the periodic one-dimensional Hubbard model at `t=1`
and `U/t=2`. The canonical input is tied to source commit
`88c21ae3345ad7414d86b37e777e10c35e6451f2` and contains exactly 20
observations: 8 finite-temperature DQMC points, 5 ground-state ED references,
and 7 ground-state DMRG references.

The fitted quantity is named `c_eff_diagnostic`. It is a finite-size slope
diagnostic. It is **not** an asymptotic central-charge determination, and its
nominal regression standard error is not a physical precision estimate. DQMC
and ED/DMRG are never pooled into one regression because they use different
ensembles and state conventions.

Two ground-state reference variants are always retained:

- `reference_without_L26` is the baseline ED+DMRG diagnostic through L=24.
- `reference_with_L26_candidate` adds the convergence-qualified but provisional
  L=26 DMRG candidate as a sensitivity diagnostic.

## Reproducible analysis environment

The root `Project.toml` and ignored root `Manifest.toml` remain frozen. Their
SHA-256 values are, respectively,
`236d1edb91a54227e6ee68f1b3959d8a1b19f24b3b72e4591204322d24e1073f`
and
`edbc59b60f4c017c2627cec9c94d92264f549d861e88ff8ec800b21f977b2d6c`.
The FSS analysis uses only the dedicated tracked environment:

| Item | Recorded value |
|---|---|
| Canonical platform | Ubuntu under WSL2, x86_64 |
| Julia | 1.12.1 |
| Julia launcher | `/home/yqw/.juliaup/bin/julia` |
| Launcher SHA-256 | `8ac1a8ee6370064715c37d1a117baf15aeab88a63c65babf069963c17241cb40` |
| Julia binary | `/home/yqw/.julia/juliaup/julia-1.12.1+0.x64.linux.gnu/bin/julia` |
| Julia binary SHA-256 | `4ae85cc900bff34224392ac6e41328642ea306429ecc5805518090e732cbfadd` |
| Plots | 1.41.2 |
| GR / GR_jll | 0.73.26 / 0.73.26+0 |
| FSS `Project.toml` SHA-256 | `e4319ccc1893f9a5ecd5e3860166d1a71f06ad97fe92c42f8b1279a206732c71` |
| FSS `Manifest-v1.12.toml` SHA-256 | `d09a5593637c67beb20218db12e32d078ff1ea818d24aa3cb473138dbc376f5d` |
| Julia threads during resolution/load test | 1 |

The project declares `LinearAlgebra`, `Printf`, and `Plots`, with
`julia = "=1.12.1"` and `Plots = "~1.41"`. Here `~1.41` means
`>=1.41.0,<1.42.0`; it is a patch-series restriction, unlike caret-compatible
`1.41`. The exact manifest locks Plots 1.41.2 and the Julia/JLL wrapper graph.
It cannot by itself guarantee byte-identical graphics across operating
systems, installed fonts, or native rendering stacks. The canonical rendering
environment is the platform above.

Package resolution used unique native `/tmp` locations:

```text
project=/tmp/HubbardRenyiDQMC_fss_project_corrected_p6JBfYYh
depot=/tmp/HubbardRenyiDQMC_fss_depot_corrected_58A63YND
TMP=TEMP=TMPDIR=/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1
```

The exact successful resolution command was:

```bash
TMP=/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1 TEMP=/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1 TMPDIR=/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1 JULIA_DEPOT_PATH=/tmp/HubbardRenyiDQMC_fss_depot_corrected_58A63YND JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/yqw/.juliaup/bin/julia --startup-file=no --history-file=no --project=/tmp/HubbardRenyiDQMC_fss_project_corrected_p6JBfYYh -e 'using Pkg; println("julia_executable=", joinpath(Sys.BINDIR, Base.julia_exename())); println("julia_version=", VERSION); println("tmpdir=", tempdir()); println("depot=", DEPOT_PATH); Pkg.Registry.add("General"); Pkg.add(Pkg.PackageSpec(name="Plots", version=v"1.41.2")); Pkg.resolve(); Pkg.instantiate(); deps=Pkg.dependencies(); plots=only(filter(x -> x.second.name == "Plots", collect(deps))); @assert plots.second.version == v"1.41.2"; println("plots_version=", plots.second.version); Pkg.status(mode=Pkg.PKGMODE_MANIFEST)'
```

The first attempt inherited a Windows temporary directory through WSL and was
interrupted before the environment lock completed after package tree-hash
checks exposed DrvFS executable-mode mismatches. It performed resolver/download
activity but produced no manifest and did not update the temporary project.
The one permitted correction pass fixed
`TMP`, `TEMP`, and `TMPDIR` to the native directory above; it did not select a
different package version or mutate the root/global environment.

The exact successful load-and-render test was:

```bash
TMP=/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1 TEMP=/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1 TMPDIR=/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1 JULIA_DEPOT_PATH=/tmp/HubbardRenyiDQMC_fss_depot_corrected_58A63YND JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=1 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 GKSwstype=100 /home/yqw/.juliaup/bin/julia --startup-file=no --history-file=no --project=/tmp/HubbardRenyiDQMC_phaseB_candidates_Rn5aI9m5/environments/fss -e 'using LinearAlgebra, Printf, Plots; import Pkg; deps=Pkg.dependencies(); plots=only(filter(x -> x.second.name == "Plots", collect(deps))); gr=only(filter(x -> x.second.name == "GR", collect(deps))); @assert VERSION == v"1.12.1"; @assert plots.second.version == v"1.41.2"; println("julia_executable=", joinpath(Sys.BINDIR, Base.julia_exename())); println("julia_version=", VERSION); println("plots_version=", plots.second.version); println("gr_version=", gr.second.version); println("threads=", Threads.nthreads()); p=plot([1.0,2.0],[1.0,4.0],label="load test",size=(320,240)); savefig(p,"/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1/load_test.png"); savefig(p,"/tmp/HubbardRenyiDQMC_fss_tmp_corrected_48GIShJ1/load_test.svg"); println("load_test=PASS")'
```

## Canonical input and provenance

The committed input is
`results/processed/oneD_hubbard_fss_fit_input.csv`, SHA-256
`713cb6de4f5be017debf3d3a39701b47307d333d3908f2ca6106d914b74d1427`
before Git promotion. Source numeric tokens are copied verbatim; generated
numeric output uses round-trip-safe serialization rather than `%.12g`.

For tracked sources, every `measurement_source_sha256` and
`parameter_source_sha256` is the SHA-256 of the canonical Git blob bytes from
`88c21ae3345ad7414d86b37e777e10c35e6451f2`, equivalent to
`git show <commit>:<path>`. The adjacent Git blob OID disambiguates identity.
This matters for CRLF-marked files: a checkout's working-byte SHA-256 can
differ even though the committed blob is identical. For the ignored L26 log,
the SHA-256 is over raw file bytes and the Git OID field explicitly reads
`not_available_untracked_source`.

For DQMC rows, `sweeps_text=1` is the parameter file's model-level
`N_sweeps` token; it is not the total burn-in or measurement count. The latter
sampling fields and four root seeds remain traceable in each row's hashed
parameter TOML.

The convention implementations are pinned by these canonical blob hashes:

| Role | Path | Git-blob SHA-256 |
|---|---|---|
| Active DQMC entry point | `run/run_lambda_node.jl` | `04b26dc03926ff8b17bb9aabccfcf6b79446ab9e867ada1e04737bc24464e2b6` |
| DQMC validated kernel | `archive/validated_old_pipeline/run_lambda_grid_step_validated.jl` | `f0868a2d015939bbf9569abe4a45a6e5717d7de446e5402335aaba9882b26a51` |
| ED generator | `analysis/compute_s2_ed.jl` | `a8969d398fbef4144835e9d3cf730606d4d0302d21a93aa51e47c9faff266a87` |
| DMRG generator | `analysis/compute_s2_itensor_dmrg.jl` | `2b5d35294751c049eb5bc379b5c0a6e090404580847b7552b5b0fcc15ebd5019` |
| DMRG L14–24 table | `results/processed/DMRG_S2_U2_periodic.csv` | `a66bd493e4fadf975061f67c5fa879c99a837e30278291d834360ddaa0acf21c` |

The 20 canonical measurements are:

| Method | L | nA | beta/state | S2 source token | absolute MC 1σ token | status |
|---|---:|---:|---|---:|---:|---|
| DQMC | 4 | 2 | 16.0 | `1.113651203984e+00` | `9.975514074400e-03` | accepted processed snapshot |
| DQMC | 6 | 3 | 24.0 | `1.533082262570e+00` | `2.249188431683e-03` | accepted processed snapshot |
| DQMC | 8 | 4 | 32.0 | `1.318734225375e+00` | `7.790202477854e-03` | accepted processed snapshot |
| DQMC | 10 | 5 | 40.0 | `1.716441955736e+00` | `3.040377094374e-03` | accepted processed snapshot |
| DQMC | 12 | 6 | 48.0 | `1.561005400404e+00` | `7.601610411883e-03` | accepted processed snapshot |
| DQMC | 14 | 7 | 56.0 | `1.843169617829e+00` | `3.662701883375e-03` | accepted processed snapshot |
| DQMC | 16 | 8 | 64.0 | `1.705923245750e+00` | `7.900417155705e-03` | accepted processed snapshot |
| DQMC | 18 | 9 | 72 | `1.949714833089e+00` | `4.306010012762e-03` | accepted snapshot; burn-in 8000 |
| ED | 4 | 2 | ground state | `0.880358089685708` | — | deterministic reference |
| ED | 6 | 3 | ground state | `1.535020904892800` | — | deterministic reference |
| ED | 8 | 4 | ground state | `1.309407045625562` | — | deterministic reference |
| ED | 10 | 5 | ground state | `1.716551661433` | — | repository-authoritative checkpoint |
| ED | 12 | 6 | ground state | `1.542354369357964` | — | deterministic reference |
| DMRG | 14 | 7 | ground state | `1.8460982712452576e+00` | — | first-pass reference |
| DMRG | 16 | 8 | ground state | `1.6962908755077759e+00` | — | first-pass reference |
| DMRG | 18 | 9 | ground state | `1.9442405246094483e+00` | — | first-pass reference |
| DMRG | 20 | 10 | ground state | `1.8087666550463910e+00` | — | first-pass reference |
| DMRG | 22 | 11 | ground state | `2.0219147921796217e+00` | — | first-pass reference |
| DMRG | 24 | 12 | ground state | `1.8960694735263326e+00` | — | first-pass reference |
| DMRG | 26 | 13 | ground state | `2.0852574848816734` | — | provisional convergence-qualified candidate |

The L=10 canonical ED value deliberately remains the repository-authoritative
`1.716551661433`. The tracked journal-check table contains
`1.716551661443483`, which is `+1.0483e-11` higher. It is a precision
cross-check, not a replacement observation, and the difference has no effect
at the DQMC uncertainty scale. The journal-check source is
`results/processed/ED_L10_U1_U2_U4_betaInf_nA1to9_journal_check.csv`, canonical
Git-blob SHA-256
`0386905163d79d2690afc2d847cfa45d94d95daf2511dd1dc429f44668bbbb75`.

## Physical conventions

All observations use `t=1`, `U=2`, `U/t=2`, half filling, a periodic
nearest-neighbor chain, `nA=L/2`, and the contiguous subsystem sites `1:nA`.
They represent the order-2 Rényi definition
`S2=-log(Tr(rho_A^2))`. The DQMC estimator obtains the same quantity through
the replica identity `S2=-log(Z1/Z0)`.

The state conventions differ and are kept explicit:

| Method | Ensemble | State convention | Filling convention |
|---|---|---|---|
| DQMC | grand canonical | finite beta, beta=4L | mu=0; half filling by particle-hole symmetry |
| ED | fixed particle number | ground state | Nup=Ndn=L/2 |
| DMRG | fixed conserved quantum numbers | ground state | Nup=Ndn=L/2 and Sz=0 |

Thus the DQMC-versus-reference comparison tests finite-temperature
convergence and method consistency; it does not assert that the ensembles are
identical.

## L26 eligibility evidence

The canonical L26 observation is the seed-1234, 22-sweep, configured
`maxdim=5600` run. The seed-11234 run is an independent convergence validation
only; it is not averaged with the canonical observation and is not a Monte
Carlo replica.

| Resource level | Seed | sweeps / maxdim | S2 | energy density | normcheck | final truncation error |
|---|---:|---:|---:|---:|---:|---:|
| level 1 | 1234 | 18 / 5200 | `2.0852568099861148` | `-0.84679839136275659` | `0.99999999999999922` | `5.00e-09` |
| level 2 canonical | 1234 | 22 / 5600 | `2.0852574848816734` | `-0.84679839487120911` | `1.0000000000000018` | `3.56E-09` |
| level 2 repeat | 11234 | 22 / 5600 | `2.0852574834044932` | `-0.84679839485401831` | `0.99999999999999944` | `3.55E-09` |

The absolute level1–level2 changes are `6.748955585500482e-07` in S2 and
`3.5084525285000723e-09` in energy density. The level2–repeat changes are
`1.4771801559732012e-09` and `1.7190804335598386e-11`, respectively. The
maximum norm deviation is `1.7763568394002505e-15`, and every final
truncation error is below `1e-8`. These values pass the accepted L26 thresholds
of `1e-4`, `1e-8`, `1e-10`, and `1e-8`.

The canonical log is
`results/raw/dmrg_queue/L26_U2_periodic_nA13_convergence_20260806T154756Z_oBDewq/level2_maxdim5600_seed1234.log`,
raw-byte SHA-256
`6eabf88493dfb05aa0dd0d5adc6edeb2fef9feb9dda446369291b772a207ff20`.
The repeat log SHA-256 is
`cd63028ae843e191d6314b86fa0e4da3d129378fffbc8f9fa3e6b325695469a7`.
The preserved maxdim-4800 first-pass log did not record `sum_p2`, so its S2
cannot be reconstructed independently from raw Schmidt-probability
diagnostics. No common-size ED–DMRG S2 overlap artifact exists at this commit;
`overlap_status=not_available_at_HEAD`.

The campaign directory is ignored and is not part of the commit. A fresh clone
therefore has the committed input's paths, selectors, diagnostics, and hashes,
but not the raw L26 evidence. That is a material archival limitation. L26 is
eligible only as a provisional Phase B candidate, and all scientific reading
must compare fits with and without it.

## Regression models and uncertainty semantics

All fits use `x=log(L/pi)` and write the finite-size slope as
`c_eff_diagnostic/4`:

```text
total:       S2 = (c_eff_diagnostic/4) log(L/pi) + b
branch:      same form, separately for L mod 4 = 0 and L mod 4 = 2
oscillation: S2 = (c_eff_diagnostic/4) log(L/pi) + b
                  + A (-1)^(L/2) / L^p, for fixed p=1 and p=2
```

For DQMC, each reported Monte Carlo error is treated as an independent
absolute one-sigma value because no cross-size covariance estimate is
available. The analysis uses
`Sigma=diag(sigma_i^2)`, `W=Sigma^-1`, and covariance
`(X' W X)^-1` without reduced-chi-square rescaling. Coefficient errors are
named `nominal_wls_se`; the reported statistic is `chi2_per_dof`. A large
value diagnoses model inadequacy rather than increased precision.

ED and DMRG have no fabricated statistical error. Their OLS scatter estimate
uses `s^2 (X'X)^-1` with `s^2=RSS/dof`, names coefficient errors
`nominal_ols_se`, and reports `rss_per_dof`. These are nominal regression
scatter values, not ED precision, DMRG convergence uncertainty, or physical
uncertainty.

Production coefficients are solved by QR and independently checked with SVD.
Every attempted fit records dataset variant, analysis role, model, branch,
fixed `p`, L range, omitted point, input hash, point and parameter counts,
degrees of freedom, SVD rank, raw and column-normalized condition numbers,
coefficients, SE semantics, statistic semantics, and L26 status. A fit is
accepted only when `n_points>n_parameters`, `dof>0`, the SVD rank is full,
all DQMC sigma values are finite and positive, and the column-normalized
condition number is below `1e6` (warning at `1e4`). Independent numeric checks
use

```text
abs(expected-actual) <= 1e-10 + 1e-9*max(abs(expected),abs(actual)).
```

## Principal diagnostic fits

The generated summary contains 273 fit rows and has SHA-256
`a8914f19c8843aeadd1e5c0417ed8ad3b7205b67b8ac639d4db05920b0e7eef4`.
Its 15 full-window principal fits are:

| Dataset variant | Model | branch / p | c_eff_diagnostic | nominal SE | statistic |
|---|---|---|---:|---:|---:|
| DQMC | total | all | 1.49398424 | 0.01350933 WLS | chi2/dof = 573.688127 |
| DQMC | branch | L mod 4 = 0 | 1.76798678 | 0.03469582 WLS | chi2/dof = 41.484966 |
| DQMC | branch | L mod 4 = 2 | 1.49098672 | 0.01468874 WLS | chi2/dof = 4.538057 |
| DQMC | oscillation | p=1 | 1.78809368 | 0.01454855 WLS | chi2/dof = 95.113124 |
| DQMC | oscillation | p=2 | 1.77041742 | 0.01483990 WLS | chi2/dof = 283.264225 |
| reference_without_L26 | total | all | 2.07499959 | 0.34011262 OLS | RSS/dof = 0.02355283 |
| reference_without_L26 | branch | L mod 4 = 0 | 2.27199423 | 0.05302021 OLS | RSS/dof = 0.00038574 |
| reference_without_L26 | branch | L mod 4 = 2 | 1.50238357 | 0.01936392 OLS | RSS/dof = 0.00002466 |
| reference_without_L26 | oscillation | p=1 | 1.78714012 | 0.06070736 OLS | RSS/dof = 0.00069368 |
| reference_without_L26 | oscillation | p=2 | 1.65606847 | 0.20458483 OLS | RSS/dof = 0.00693769 |
| reference_with_L26_candidate | total | all | 2.10626711 | 0.30095100 OLS | RSS/dof = 0.02134391 |
| reference_with_L26_candidate | branch | L mod 4 = 0 | 2.27199423 | 0.05302021 OLS | RSS/dof = 0.00038574 |
| reference_with_L26_candidate | branch | L mod 4 = 2 | 1.50865683 | 0.01498194 OLS | RSS/dof = 0.00002074 |
| reference_with_L26_candidate | oscillation | p=1 | 1.80699885 | 0.05660797 OLS | RSS/dof = 0.00068537 |
| reference_with_L26_candidate | oscillation | p=2 | 1.72629703 | 0.18817772 OLS | RSS/dof = 0.00681913 |

The very large DQMC `chi2_per_dof` values reject a precision interpretation of
the small nominal WLS errors and expose substantial finite-size model
inadequacy. The reference results also depend strongly on whether the total,
branch, or oscillatory model is used.

These quantities must be read as finite-window diagnostics. In particular,
neither apparent stability nor a small nominal regression SE establishes an
asymptotic central charge.

## Window and leave-one-out sensitivity

The analysis evaluates the approved Lmin sets, every contiguous five-size
all-size window, every contiguous three-point window within each `L mod 4`
branch, and leave-one-out variants of every full-window baseline model. It
does not apply leave-one-out to a minimal window when omission would
underdetermine the fit.

For the all-size total model, the observed `c_eff_diagnostic` ranges are:

| Dataset variant | Lmin suite | contiguous five-size windows | full-fit leave one out |
|---|---:|---:|---:|
| DQMC | 1.39553–1.94512 | 1.32269–1.81829 | 1.38656–2.08040 |
| reference_without_L26 | 1.44735–2.07500 | 1.50163–2.38403 | 1.56970–2.37278 |
| reference_with_L26_candidate | 1.66056–2.15461 | 1.50163–2.38403 | 1.67741–2.38703 |

The smaller branch and oscillatory windows are at least as model-sensitive.
Across all requested Lmin, contiguous-window, and leave-one-out roles, DQMC
branch fits span 1.45699–2.24589 and DQMC oscillation fits span
1.58979–2.51667. The corresponding reference-without-L26 spans are
1.46472–2.41737 and 1.32810–2.35488; the with-L26 set retains those broad
extrema because windows that end before L26 are intentionally still present.

For otherwise identical principal reference fits, adding the provisional L26
candidate shifts `c_eff_diagnostic` as follows:

| Model | with-L26 minus without-L26 |
|---|---:|
| total | +0.03126752 |
| L mod 4 = 0 branch | 0.00000000 |
| L mod 4 = 2 branch | +0.00627326 |
| oscillation p=1 | +0.01985873 |
| oscillation p=2 | +0.07022856 |

The zero shift for `L mod 4 = 0` is structural: L=26 belongs to
`L mod 4 = 2`. The range widths and model-to-model differences are much larger
than a precision-level reading of any one nominal SE, which is why the
analysis retains the full sensitivity rows rather than selecting a preferred
window after inspection.

## Generated artifacts and validation

The analysis CLI accepts `--input`, `--summary-output`, `--figure-dir`, and
`--no-clobber`. A PowerShell-compatible invocation for fresh output paths is:

```powershell
$env:GKSwstype = '100'
$run = Join-Path $env:TEMP ('HubbardRenyiDQMC_fss_' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $run | Out-Null
$figures = Join-Path $run 'figures'
New-Item -ItemType Directory -Path $figures | Out-Null
$summary = Join-Path $run 'oneD_hubbard_fss_fit_summary.csv'
julia --startup-file=no --history-file=no --project=.\environments\fss .\analysis\fit_oneD_hubbard_fss_with_oscillation.jl --input .\results\processed\oneD_hubbard_fss_fit_input.csv --summary-output $summary --figure-dir $figures --no-clobber
```

`--no-clobber` refuses any existing requested output. Candidate outputs were
generated in a unique `/tmp` directory before promotion. The validated result
set consists of the summary plus four PNG/SVG pairs: canonical input,
all-size total fit, `L mod 4` branch fits, and fixed-p oscillation fits.

The candidate analysis script has SHA-256
`aab2c712e83f91d575d02eb5ec657e0fafe4afa57aea3e850f4acf2412ba3063`.
CLI help and negative tests for missing/unknown arguments, malformed input,
schema/count/convention failures, and no-clobber behavior passed.

An independent repeated run under
`/tmp/HubbardRenyiDQMC_phaseB_fit_canon1_82m8A8pa` reproduced all nine
generated outputs (summary plus eight figure files) byte for byte. The summary
contains 273 accepted fits: 15 principal, 58 Lmin, 76 contiguous-window, and
124 leave-one-out rows. Every row has positive degrees of freedom and full SVD
rank. The maximum column-normalized condition number is
`31.02182418544026`, far below `1e6`, and the largest QR-versus-SVD
coefficient difference is `8.171241461241152e-14`, within the required
absolute-plus-relative tolerance.

All four PNGs have valid signatures and CRCs. The input, total, and
oscillation PNGs are 2160×1440; the branch PNG is 2700×1368. All four SVGs
parse as XML. Their corresponding dimensions are 1200×800 and 1500×760.
Raster and vector versions use the same panel structure, labels, legends, and
data selections; explicit structural checks find the provisional L26 marker
in the canonical-input view and every required with-L26 fit view. Visual
review found no clipped axes, titles, annotations, or legends. Candidate
script, summary, report, and figure modes are 0644.

## Limitations

- DQMC is finite-temperature and grand-canonical, while ED/DMRG is a
  fixed-particle-number ground-state reference family.
- Cross-size DQMC covariances are unavailable, so diagonal independent-sigma
  weighting is an explicit approximation.
- The tracked L14–L24 DMRG rows have no reported statistical error and were not
  subjected to the new systematic two-level convergence campaign. Their
  configured schedule cap was 4800; the input separately records observed
  `maxlinkdim`.
- The L26 campaign passes its numerical gate but remains provisional because
  no common-size ED–DMRG overlap artifact is available and its raw logs are
  ignored rather than committed.
- Finite-size model inadequacy and oscillatory corrections dominate the
  interpretation; nominal regression standard errors do not capture those
  systematics.
- Graphics reproducibility is canonical only for Ubuntu WSL2 x86_64 with
  Julia 1.12.1 and the checked-in manifest. Numerical results are validated by
  tolerances rather than by cross-platform byte identity.

Accordingly, this result set supports method and finite-size sensitivity
diagnostics only. It does not claim an asymptotic central-charge measurement.
