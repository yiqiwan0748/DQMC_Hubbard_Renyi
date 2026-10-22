# Grover Rank-One Speedup

> Historical validation record. See [`dqmc_renyi_ee.md`](dqmc_renyi_ee.md)
> for the current scientific summary and conclusions.

Date: 2026-08-19

## Goal

Accelerate nonzero-lambda Renyi updates by 2.5x without changing the sampled
distribution, estimator definition, RNG convention, or measured observables.

The accepted reference implementation under `archive/validated_old_pipeline/`
remains unchanged. The optimized implementation is isolated in
`src/grover_rank1_update.jl` and is selected by the active lambda-node runner.

## Preserved Identity

For subsystem Green's-function blocks `G1` and `G2`, the Grover matrix is

```text
g = G1*G2 + (I - G1)*(I - G2).
```

A local auxiliary-field proposal changes the active block by

```text
G1_new = G1 + p*q^T.
```

Therefore

```text
g_new = g + p*r^T,
r^T   = q^T*(2*G2 - I),
det(g_new)/det(g) = 1 + r^T*inv(g)*p.
```

The determinant ratio is evaluated in `O(nA^2)` instead of rebuilding and
factoring a Grover matrix in `O(nA^3)` for every proposal. The lambda weight
is still exactly

```text
exp(lambda * (logabsdet(g_new) - logabsdet(g))).
```

## Redundant Work Removed

- The current Grover determinant and inverse are carried between proposals.
- Accepted proposals use rank-one updates; rejected proposals do not copy or
  restore full matrices.
- Between numerical-stabilization boundaries, only `G00[C,C]` is needed for
  the Grover factor. This block is maintained by the tracker, so redundant
  full `G00` rank-one updates are skipped.
- JDQMCFramework reconstructs the full `G00` from stabilized propagators at
  every required boundary and at sweep completion. The tracker is
  resynchronized at those boundaries, leaving the measured full state intact.

## Numerical Safeguards

- Scale-aware Sherman-Morrison denominator check.
- Condition estimate refreshed after every accepted update.
- Two-sided inverse residual checks at regular intervals.
- Explicit Grover rebuild and LU factorization when a guard fails.
- Ill-conditioned matrices remain on the explicit fallback path.
- The explicit formula uses the same unexpanded operation order as the
  validated implementation to avoid additional cancellation.

## Validation

The validation harness is `analysis/validate_grover_rank1_update.jl`.

It covers:

- 400 random determinant-lemma checks.
- 1,401 ordinary tracker proposal checks, including a required fast-path count.
- Sequential accepted updates and periodic inverse checks.
- Rejected proposals and forced explicit fallbacks.
- Ill-conditioned fallback acceptance and the following proposal.
- Twelve same-seed direct-versus-optimized DQMC rounds, including auxiliary
  fields, all Green's functions, determinant metadata, Grover measurements,
  acceptance rates, and RNG-stream position.

All 2,318 checks pass. A separate active-runner smoke node completed and wrote
valid summary and bin CSV files only to a temporary directory.

## Performance Result

Benchmark parameters:

```text
L=8, beta=32, dtau=0.05, nA=4, lambda=0.5, n_stab=10
20 two-replica update rounds, 5 repeats, BLAS threads=1
```

Each timed section begins after warm-up and an isolated garbage collection.

| Metric | Validated direct kernel | Rank-one kernel |
|---|---:|---:|
| Mean time | 0.550769 s | 0.204545 s |
| Mean allocation | 1278.946 MiB | 11.415 MiB |
| Mean GC time | 0.073517 s | 0.000000 s |

Measured mean speedup:

```text
2.692658x
```

This exceeds the 2.5x target. Applied to the prior 44.6-second matched
workload, the projected runtime is about 16.6 seconds; a long production-sized
timing should be used to replace that projection before publishing a final
performance claim.

## Pre-Commit Physics Checks

All fresh outputs below were written under
`/tmp/hubbard_rank1_ed_dmrg_20260819`; repository raw and processed snapshots
were not modified.

### L=10 ED checkpoint

The complete `L10_U2_beta40_Ninc8` grid was rerun with all four standard root
seeds: 36 lambda-seed nodes, 100 burn-in sweeps, and 1000 measurement sweeps.

```text
S2_DQMC = 1.709038001302 +/- 0.008325520030
S2_ED   = 1.716551661433
DeltaS2 = -0.007513660131
z       = -0.902
```

All 36 optimized bin files match the historical validated bin files exactly
after excluding the intentionally different run-name field. Mean runtime for
the 28 interior-lambda nodes changed from 95.9854 seconds historically to
25.1035 seconds, a production-workload speedup of 3.8236x.

### L=14 DMRG checkpoint

A fresh one-seed `L14_U2_beta56_Ninc16` grid used 100 burn-in sweeps and 1000
measurement sweeps per node.

```text
S2_DQMC = 1.843461576701 +/- 0.015942703821
S2_DMRG = 1.8460982712452576
DeltaS2 = -0.002636694544
z_DQMC  = -0.165
```

The z value uses only the DQMC Monte Carlo error because the repository's DMRG
table does not quote an uncertainty. This comparison is a finite-temperature,
grand-canonical DQMC result against a fixed-sector ground-state DMRG result, so
it is a consistency check rather than an equality test between identical
ensembles.

## Reproduction

From PowerShell, using the tracked repository environment:

```powershell
julia --project=. --startup-file=no --history-file=no .\analysis\validate_grover_rank1_update.jl --benchmark
```
