"""Compare the public Julia DMRG wrapper with independent TeNPy optimization.

Uses TeNPy's stock FermiHubbardModel and TwoSiteDMRGEngine. The only energy
conversion is the documented centered-interaction constant; no Julia MPO or
wavefunction is supplied to TeNPy. Run with --help for paths.
"""
from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import logging
import os
import platform
import sys
import time
import tomllib
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import scipy
import tenpy
from tenpy.algorithms import dmrg
from tenpy.models.hubbard import FermiHubbardModel
from tenpy.networks.mps import MPS


def stock_model(case):
    # TeNPy default Square ordering has its y coordinate fastest. Transposing
    # coordinates makes site i be public (x=i%Lx, y=i//Lx), without snake turns.
    model = FermiHubbardModel({
        "lattice": "Square", "Lx": case["Ly"], "Ly": case["Lx"],
        "bc_x": "open", "bc_y": "open", "bc_MPS": "finite",
        "order": "default", "t": case["t"], "U": case["U"],
        "mu": case["mu"] + case["U"] / 2.0,
        "cons_N": "N", "cons_Sz": "Sz",
    })
    n = case["Lx"] * case["Ly"]
    expected_order = np.array([[i // case["Lx"], i % case["Lx"], 0] for i in range(n)])
    np.testing.assert_array_equal(model.lat.order, expected_order)
    expected_bonds = set()
    for y in range(case["Ly"]):
        for x in range(case["Lx"]):
            i = x + y * case["Lx"]
            if x + 1 < case["Lx"]:
                expected_bonds.add((i, i + 1))
            if y + 1 < case["Ly"]:
                expected_bonds.add((i, i + case["Lx"]))
    actual_bonds = set()
    for u1, u2, dx in model.lat.pairs["nearest_neighbors"]:
        left, right, _, _ = model.lat.possible_couplings(u1, u2, dx)
        actual_bonds.update(tuple(sorted((int(i), int(j)))) for i, j in zip(left, right))
    assert actual_bonds == expected_bonds, (actual_bonds, expected_bonds)
    return model


def initial_state(case, seed):
    n = case["Lx"] * case["Ly"]
    # Independent initial states in the same conserved sector. Reversing the
    # odd/even occupation priority changes the start while preserving counts.
    up_order = list(range(0, n, 2)) + list(range(1, n, 2))
    down_order = list(range(1, n, 2)) + list(range(0, n, 2))
    if seed != 1234:
        up_order, down_order = up_order[::-1], down_order[::-1]
    up, down = set(up_order[:case["Nup"]]), set(down_order[:case["Ndn"]])
    return ["full" if i in up and i in down else "up" if i in up else "down" if i in down else "empty" for i in range(n)]


def run_case(case, settings):
    np.random.seed(settings["root_seed"])
    model = stock_model(case)
    psi = MPS.from_product_state(model.lat.mps_sites(), initial_state(case, settings["root_seed"]),
                                 bc="finite", unit_cell_width=model.lat.mps_unit_cell_width)
    options = {
        "mixer": True,
        "mixer_params": {"amplitude": 1e-7, "decay": 10.0, "disable_after": 4},
        "trunc_params": {"chi_max": settings["maxdim"][-1], "svd_min": settings["cutoff"]},
        "chi_list": {i: chi for i, chi in enumerate(settings["maxdim"])},
        "min_sweeps": settings["nsweeps"], "max_sweeps": settings["nsweeps"],
        "max_E_err": 1e-12, "max_S_err": 1e-10,
        "N_sweeps_check": 1,
        "lanczos_params": {"N_min": 3, "N_max": 100, "P_tol": 1e-14, "E_tol": 1e-14},
        "combine": True,
    }
    engine = dmrg.TwoSiteDMRGEngine(psi, model, options)
    start = time.perf_counter()
    energy_uncentered, psi = engine.run()
    elapsed = time.perf_counter() - start
    psi.canonical_form()
    n = case["Lx"] * case["Ly"]
    energy = float(np.real(model.H_MPO.expectation_value(psi))) + n * case["U"] / 4.0
    entropy = float(psi.entanglement_entropy(n=2)[case["cut"] - 1])
    energies = [float(v) for v in engine.sweep_stats["E"]]
    return {
        "case": case["name"], "setting": settings["name"], "package": "TeNPy",
        "package_version": tenpy.__version__, "Lx": case["Lx"], "Ly": case["Ly"],
        "t": case["t"], "U": case["U"], "mu": case["mu"],
        "Nup": case["Nup"], "Ndn": case["Ndn"], "cut": case["cut"],
        "root_seed": settings["root_seed"], "nsweeps": engine.sweeps,
        "maxdim_cap": settings["maxdim"][-1], "cutoff": settings["cutoff"],
        "energy": energy, "S2": entropy,
        "density": float(np.mean(psi.expectation_value("Ntot"))),
        "double_occupancy": float(np.mean(psi.expectation_value("NuNd"))),
        "maxlinkdim": max(psi.chi),
        "sweep_energy_change": abs(energies[-1] - energies[-2]),
        "final_truncation_error": float(engine.sweep_stats["max_trunc_err"][-1]),
        "elapsed_seconds": elapsed,
        "reported_engine_energy_uncentered": float(energy_uncentered),
        "energy_offset": n * case["U"] / 4.0,
    }


def compare(public_path, tenpy_rows, output):
    with public_path.open(newline="", encoding="utf-8") as stream:
        public = {(r["case"], r["setting"]): r for r in csv.DictReader(stream)}
    rows = []
    tolerances = {"energy": 1e-8, "S2": 1e-7, "density": 1e-10, "double_occupancy": 1e-8}
    for reference in tenpy_rows:
        own = public[(reference["case"], reference["setting"])]
        for observable, tolerance in tolerances.items():
            own_value, reference_value = float(own[observable]), float(reference[observable])
            error = abs(own_value - reference_value)
            rows.append({"case": reference["case"], "setting": reference["setting"],
                         "observable": observable, "public_value": own_value,
                         "tenpy_value": reference_value, "absolute_error": error,
                         "tolerance": tolerance, "pass": error <= tolerance})
    with output.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    failed = [row for row in rows if not row["pass"]]
    print(json.dumps({"checks": len(rows), "passed": len(rows) - len(failed),
                      "failed": failed, "max_absolute_errors": {
                          key: max(r["absolute_error"] for r in rows if r["observable"] == key)
                          for key in tolerances}}, indent=2), flush=True)
    return not failed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public", type=Path, required=True, help="CSV generated by dmrg_public.jl")
    parser.add_argument("--output-dir", type=Path, required=True, help="new results directory")
    parser.add_argument("--compare-only", action="store_true", help="compare existing tenpy.csv after the public run completes")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.compare_only:
        with (args.output_dir / "tenpy.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        if not compare(args.public, rows, args.output_dir / "comparison.csv"):
            raise SystemExit(1)
        return
    for filename in ("tenpy.csv", "comparison.csv", "environment.json"):
        if (args.output_dir / filename).exists():
            raise FileExistsError(args.output_dir / filename)
    with (Path(__file__).parent / "dmrg_cases.toml").open("rb") as stream:
        cfg = tomllib.load(stream)
    logging.basicConfig(level=logging.WARNING)
    environment = {"python": sys.version, "platform": platform.platform(), "tenpy": tenpy.__version__,
                   "numpy": np.__version__, "scipy": scipy.__version__, "engine": "TwoSiteDMRGEngine",
                   "reference_model": "stock FermiHubbardModel on open Square lattice",
                   "interaction_conversion": "TeNPy mu = public mu + U/2; add Nsites*U/4 to TeNPy energy",
                   "site_order": "TeNPy axes transposed; exact match to public x-fast physical prefix",
                   "packages": {p.metadata["Name"]: p.version for p in importlib.metadata.distributions()}}
    (args.output_dir / "environment.json").write_text(json.dumps(environment, indent=2), encoding="utf-8")
    results = []
    with (args.output_dir / "tenpy.csv").open("x", newline="", encoding="utf-8") as stream:
        writer = None
        for case in cfg["cases"]:
            for settings in cfg["settings"]:
                print(f"TeNPy DMRG: {case['name']} / {settings['name']}", flush=True)
                row = run_case(case, settings)
                if writer is None:
                    writer = csv.DictWriter(stream, row.keys())
                    writer.writeheader()
                writer.writerow(row)
                stream.flush()
                results.append(row)
    if not compare(args.public, results, args.output_dir / "comparison.csv"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
