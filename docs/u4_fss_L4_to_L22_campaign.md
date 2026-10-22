# U=4 PBC FSS campaign through L=22

## Objective

Repeat the accepted U=2 half-chain finite-size series at U=4 and compare DQMC
with ground-state ED/DMRG references. This is a new PBC series; the earlier
L=8,12 APBC screening data are not mixed into it.

## Fixed matrix

- `L = 4,6,...,22`, PBC, `nA=L/2`, `beta=4L`
- `U=4`, `dtau=0.05`, `n_stab=10`, one sweep per measurement
- `Ninc=16`, 4,000 measurements, 40 bins
- roots `1234`, `11234`, `21234`, `31234`
- 8,000 burn-in sweeps and four simultaneous DQMC workers

The longer burn-in is intentional. The U=4 APBC screening pilot with 500
burn-in sweeps had an anomalously high early L=12 block and finished 2.90 DQMC
standard errors above ED. Starting with 8,000 avoids making that known risk the
reason for a full production rerun.

Expected DQMC wall time is approximately 40--48 hours. The sequential ED/DMRG
references are expected to add roughly 4--8 hours. Raw DQMC storage is
approximately 11 MiB.

## References

- ED: PBC at L=4,6,8,10,12.
- DMRG: PBC at L=12,14,16,18,20,22, with L=12 as the ED overlap gate.
- ED stops at L=12: the current full-reorthogonalization implementation would
  exceed this machine's memory at larger sizes.
- ED is accepted only when its reported Lanczos residual is below `1e-11`.
- Every primary DMRG reference uses 22 sweeps, cutoff `1e-12`, seed 1234, and
  the explicit schedule ending with four maxdim-5600 sweeps.
- L=20 and L=22 additionally run a maxdim-5200 resource check and an independent
  maxdim-5600 seed-11234 repeat. Across all three, the S2 span must be at most
  `1e-4`, the energy-density span at most `1e-8`, and each final-sweep maximum
  truncation error at most `1e-8`.

DQMC here is finite-temperature (`beta=4L`), grand-canonical, and uses finite
`dtau=0.05` (`U*dtau=0.2`). ED/DMRG are fixed-particle-number ground states.
The comparison column `z_DQMC_only` divides by the reported DQMC Monte Carlo
error only; it does not include DMRG convergence uncertainty or beta/dtau
systematics.

## Detached and recoverable execution

The authoritative owner is one Windows scheduled task,
`HubbardRenyiDQMC-U4-FSS`. It runs the idempotent WSL campaign driver every five
minutes for up to 14 days with overlapping starts disabled. AC sleep is already
disabled on the desktop.

The task survives laptop/SSH disconnection. It uses the same non-admin
Interactive principal as the previously successful DQMC tasks on this desktop,
so the desktop account must remain signed in. After a WSL restart it resumes on
the next attempt; after a desktop restart it resumes after that account signs
in. The desktop still has to be powered for computation to proceed.

The launcher returns success only after the task is running and a fresh
heartbeat appears. Production sources must also be tracked, committed, and
unchanged from the recorded source manifest.

Recovery granularity is explicit:

- DQMC resumes from validated lambda/root nodes; at most four in-flight nodes
  are repeated.
- ED and DMRG resume from completed sizes; an interrupted current DMRG size is
  repeated because the DMRG implementation has no within-sweep checkpoint.
- Interrupted temporary files, incomplete CSV pairs, and attempt logs are moved
  to `quarantine/`; nothing is silently deleted and `--force` is never used.
- Three consecutive handled failures stop and disable the task for inspection.

The launch command from Windows PowerShell is:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File '\\wsl.localhost\Ubuntu\home\yqw\code\HubbardRenyiDQMC\run\start_u4_fss_detached.ps1'
```

Status is read-only:

```powershell
wsl.exe -d Ubuntu -u yqw --cd /home/yqw/code/HubbardRenyiDQMC -- /bin/bash run/run_u4_fss_campaign.sh --status /mnt/c/Users/YQW/HubbardRenyiDQMC_runs/u4_fss_L4_L22_20260822
Get-ScheduledTask -TaskName HubbardRenyiDQMC-U4-FSS | Get-ScheduledTaskInfo
```

## Completion gates

Each DQMC size must have exactly 68 summary/bin pairs, pass the strict existing
summarizer, and pass its SHA-256 manifest. Reference CSVs are schema-checked;
L=12 ED and DMRG must agree within `1e-6`, and the DMRG endpoint convergence
checks above must pass. Only then does the driver create:

- `processed/u4_fss_L4_to_L22_comparison.csv`
- `processed/u4_fss_L4_to_L22_fits.csv`
- `processed/u4_fss_L4_to_L22_s2_vs_L.svg`
- `state/final_artifacts.sha256`
- the overall `DONE` marker

The fit table contains DQMC weighted and reference unweighted total, mod-4
branch, and fixed `1/L` and `1/L^2` oscillatory fits. Every fitted `c_eff` is
explicitly a finite-window diagnostic, not an asymptotic central-charge claim.

The task disables itself only after those gates pass.
