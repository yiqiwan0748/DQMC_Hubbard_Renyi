# Regression and independent full-Fock checks

This report covers analytic, full-Fock and internal ED/DMRG checks. Actual
runs against **SmoQyDQMC, QuSpin and TeNPy** were completed separately on
October 22, 2026; see the [external-package comparison](external_packages.md).

The public package was checked from a fresh HTTPS clone of commit
`f8c6386dd6083a4a51c5ea5455b78306be444945`, then checked again after the
corrections below. Julia was 1.12.1. The complete installation and numerical
checks ran under Linux/WSL with one OpenBLAS thread, using the existing Julia
dependency cache. The output-file reservation helper was also tested with
native Windows Julia. This does not establish a fresh Windows dependency
installation or large-system convergence.

## Findings and corrections

- **Odd-site DMRG filling:** the previous defaults used `Nup=Ndn=N÷2`, which
  leaves one hole when the site count is odd. They now use `Nup=cld(N,2)` and
  `Ndn=fld(N,2)`, matching ED and preserving `Nup+Ndn=N`. Explicit particle
  counts still select the requested sector. The centered Hamiltonian and
  `S2=-log(tr(rho_A^2))` are unchanged.
- **Concurrent output writes:** the square DQMC runner now reserves its new
  directory before sampling. The DMRG runner reserves its CSV with exclusive
  creation before solving and holds that file open. Competing runs cannot
  replace an output created after an earlier existence check. A failed DQMC
  run retains its parameters without a summary; a failed DMRG run retains an
  empty reserved CSV. Retrying requires a new path.
- **DMRG provenance:** its CSV now also records hopping, both boundary flags,
  truncation cutoff, and the full bond-dimension schedule.
- **Historical documentation:** the public notes now identify omitted
  original data and scripts. A broken historical figure embed was replaced
  with its provenance path. The active workflow links resolve.

## Verification

| Check | Result |
| --- | --- |
| Fresh public clone: `Pkg.instantiate(); Pkg.test(; allow_reresolve=false)` | Passed, 458 original regression checks |
| Corrected release: standard `Pkg.test` | Passed, 494 checks including 36 output-safety checks |
| Independent full-Fock spectra and thermal observables | Passed, 60 checks |
| Interacting custom-hopping ED/DMRG comparison in five sectors | Passed, 20 checks |
| DMRG regression, including odd filling and empty/full sectors | Passed, 26 checks |
| Native Windows exclusive file creation and sentinel preservation | Passed, 5 checks |
| Real ED, DMRG, and default 4×4 DQMC commands, saved outputs, and refusal to overwrite | Passed, 34 checks |

The independent reference explicitly constructs fermion creation and
annihilation matrices on the full three-site Fock space. It tests all particle
sectors, non-nearest-neighbor hopping, diagonal site energies, a nonzero
chemical potential, and attractive, free, and repulsive interactions. It does
not reuse the package's spin-mask hopping construction.

Reproduce the numerical and output-safety checks from PowerShell:

```powershell
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
julia --project=. --startup-file=no .\analysis\crosscheck_numerics.jl
```

The independent numerical script includes the 26 DMRG regression checks.
It is separate from the default suite because ITensor compilation and DMRG
calculations take longer. Logs are in `results/crosscheck_20261021/`.

The real DMRG command wrote all 20 CSV columns and agreed with ED to
`9.77e-15` in energy and `2.45e-13` in entropy. The default square DQMC
summary, all 40 bins, and copied parameters were byte-identical to the
previous saved example. Both runners rejected an existing output without
changing its saved bytes.

## Benchmark claims

The saved timing CSV gives `mean(direct)/mean(rank1)=4.5884097528` on the 4×4
square. Thus the CV's **4.59×** is the rounded measured ratio. The seven paired
ratios range from 3.06× to 5.19×. This remains an update-kernel comparison with
direct Grover rebuilds, excluding setup and measurements.

All four numerical-source hashes recorded in the timing metadata match the
released benchmark sources. None of the corrections above changes those
sources. The ED and DMRG error claims also match their saved CSVs. No new
performance claim is made by this cross-check. See the
[DQMC](dqmc.md), [ED](ed.md), and [DMRG](dmrg.md) reports for the original
benchmark settings and limitations.
