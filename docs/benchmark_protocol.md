# Benchmark Protocol

## Baseline

Use `parameters/L10_U2_beta40_Ninc8.toml` as the default smoke benchmark.

## Required Output Metadata

Each raw output row should include:

- `run_name`
- `L`
- `t`
- `U`
- `mu`
- `beta`
- `dtau`
- `n_stab`
- `N_sweeps`
- `nA`
- `Ninc`
- `delta_lambda`
- `lambda_index`
- `lambda`
- `root_seed`
- `direction`
- `lambda_from`
- `lambda_to`
- `ratio`
- `ratio_err`
- `entropy_increment`
- `entropy_increment_err`
- `acceptance`
- `deltaG`
- `runtime_seconds`

## Required Checks

- Seed reproducibility: same parameter file and root seed should reproduce output.
- Replica independence: derived replica RNG streams must not be identical.
- Endpoint rules: `lambda = 0` only has forward ratio; `lambda = 1` only has backward ratio.
- Bidirectional consistency: `R_plus(lambda_k) * R_minus(lambda_{k+1}) ~= 1`.
- ED comparison: final `S2_DQMC` should agree with `S2_ED` within uncertainty for the default benchmark.
