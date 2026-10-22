"""Compare the public Julia ED API with independently built QuSpin operators.

Run ed_julia.jl first. QuSpin constructs its own fermion basis, Hamiltonian,
eigensystems, observable operators, and fermionic spatial partial traces.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
from pathlib import Path
import sys
import tomllib

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import numpy as np
from quspin.basis import spinful_fermion_basis_1d
from quspin.operators import hamiltonian


def kinetic(case):
    if "K" in case:
        return np.array(case["K"], dtype=float)
    lx, ly = case["shape"]
    k = np.zeros((lx * ly, lx * ly))
    for y in range(ly):
        for x in range(lx):
            i = x + y * lx
            for axis, size in enumerate((lx, ly)):
                if size == 1:
                    continue
                point = [x, y]
                point[axis] += 1
                if point[axis] == size:
                    if not case["periodic"][axis]:
                        continue
                    point[axis] = 0
                j = point[0] + point[1] * lx
                k[i, j] = k[j, i] = -case["t"]
    return k


def make_operators(case, sector=None):
    k = kinetic(case)
    n = len(k)
    basis = spinful_fermion_basis_1d(n, Nf=sector)
    hop = [[float(k[i, j]), i, j] for i in range(n) for j in range(n)
           if i != j and k[i, j] != 0]
    onsite = [[float(k[i, i] - case["mu"]), i] for i in range(n)]
    interaction = [[float(case["U"]), i, i] for i in range(n)]
    # QuSpin z = n - 1/2, so this exactly matches the centered interaction.
    static = [["+-|", hop], ["|+-", hop], ["n|", onsite],
              ["|n", onsite], ["z|z", interaction]]
    options = dict(basis=basis, dtype=np.float64, check_herm=False,
                   check_symm=False, check_pcon=False)
    h = hamiltonian(static, [], **options)
    d = hamiltonian([["n|n", [[1.0 / n, i, i] for i in range(n)]]], [], **options)
    particles = [[1.0 / n, i] for i in range(n)]
    number = hamiltonian([["n|", particles], ["|n", particles]], [], **options)
    return k, basis, h, d, number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", type=Path)
    parser.add_argument("julia_results", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    cases = tomllib.loads(args.cases.read_text(encoding="utf-8-sig"))
    julia = tomllib.loads(args.julia_results.read_text(encoding="utf-8-sig"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename in ("quspin_results.json", "summary.json", "comparison.csv", "provenance.json"):
        if (args.output_dir / filename).exists():
            raise SystemExit("Output already exists; choose a new evidence directory: " + str(args.output_dir))
    metrics = []
    external = {"ground": [], "thermal": []}
    tolerance = 2e-10

    def check(case_id, metric, reference, actual, tol=tolerance):
        expected_array, actual_array = np.asarray(reference), np.asarray(actual)
        if expected_array.shape != actual_array.shape:
            raise ValueError(f"{case_id}: {metric} shape mismatch")
        delta = float(np.max(np.abs(expected_array - actual_array)))
        metrics.append(dict(case=case_id, metric=metric, max_abs_difference=delta,
                            tolerance=tol, passed=delta <= tol))

    julia_ground = {row["id"]: row for row in julia["ground"]}
    for case in cases["ground"]:
        k, basis, h, d, number = make_operators(case, (case["Nup"], case["Ndn"]))
        energies, vectors = h.eigh()
        psi = vectors[:, 0]
        gap = float(energies[1] - energies[0]) if len(energies) > 1 else None
        renyi2 = [float(basis.ent_entropy(psi,
            sub_sys_A=(tuple(site - 1 for site in cut), tuple(site - 1 for site in cut)),
            density=False, alpha=2.0)["Sent_A"]) for cut in case["cuts"]]
        row = dict(id=case["id"], energy=float(energies[0]),
                   dimension=int(basis.Ns), gap=gap, spectrum=energies.tolist(),
                   density=float(number.expt_value(psi)),
                   double_occupancy=float(d.expt_value(psi)), renyi2=renyi2,
                   residual=float(np.linalg.norm(h.dot(psi) - energies[0] * psi)))
        reference = julia_ground[case["id"]]
        if len(reference["renyi2"]) != len(renyi2):
            raise ValueError(case["id"] + ": Renyi cut count mismatch")
        check(case["id"], "kinetic_matrix", reference["K"], k, tol=1e-14)
        check(case["id"], "external_hermiticity", h.toarray(), h.toarray().T, tol=1e-13)
        check(case["id"], "external_residual", 0.0, row["residual"])
        for metric in ("energy", "spectrum", "density"):
            check(case["id"], metric, reference[metric], row[metric])
        # A degenerate eigenvector is not a unique observable benchmark.
        if gap is None or gap > 1e-9:
            check(case["id"], "double_occupancy", reference["double_occupancy"], row["double_occupancy"])
            for i, (expected, actual) in enumerate(zip(reference["renyi2"], renyi2)):
                check(case["id"], "renyi2_cut_" + str(i + 1), expected, actual)
        else:
            row["ground_observables_excluded"] = "degenerate ground state"
        external["ground"].append(row)
        print(case["id"], "E=", row["energy"], "gap=", gap, flush=True)

    julia_thermal = {row["id"]: row for row in julia["thermal"]}
    for case in cases["thermal"]:
        k, basis, h, d, number = make_operators(case)
        energies, vectors = h.eigh()
        check(case["id"], "external_hermiticity", h.toarray(), h.toarray().T, tol=1e-13)
        residual = np.linalg.norm(h.dot(vectors) - vectors * energies, ord="fro")
        check(case["id"], "external_eigensystem_residual", 0.0, residual)
        weights = np.exp(-case["beta"] * (energies - energies[0]))
        probabilities = weights / weights.sum()
        double_diag = np.asarray(d.diagonal()).real
        number_diag = np.asarray(number.diagonal()).real
        row = dict(id=case["id"], dimension=int(basis.Ns),
            energy=float(probabilities @ energies),
            density=float(probabilities @ (np.abs(vectors).T**2 @ number_diag)),
            double_occupancy=float(probabilities @ (np.abs(vectors).T**2 @ double_diag)),
            logZ=float(np.log(weights.sum()) - case["beta"] * energies[0]))
        reference = julia_thermal[case["id"]]
        check(case["id"], "kinetic_matrix", reference["K"], k, tol=1e-14)
        for metric in ("energy", "density", "double_occupancy", "logZ"):
            check(case["id"], metric, reference[metric], row[metric])
        external["thermal"].append(row)
        print(case["id"], "thermal E=", row["energy"], flush=True)

    versions = {name: importlib.metadata.version(name) for name in
                ("quspin", "quspin-extensions", "parallel-sparse-tools", "numpy", "scipy")}
    summary = dict(passed=all(row["passed"] for row in metrics),
        comparison_metrics=len(metrics), ground_cases=len(cases["ground"]),
        thermal_cases=len(cases["thermal"]),
        spectrum_eigenvalues=sum(row["dimension"] for row in external["ground"]),
        renyi2_cuts=sum(len(row["renyi2"]) for row in external["ground"]
                       if "ground_observables_excluded" not in row),
        max_abs_difference=max(row["max_abs_difference"] for row in metrics),
        tolerance=tolerance, python=sys.version, platform=platform.platform(), versions=versions,
        julia_version=julia["julia_version"], julia_package_version=julia["package_version"],
        julia_source_path=julia["source_path"],
        julia_source_commit=julia.get("source_commit", "not recorded"),
        documentation="https://quspin.github.io/QuSpin/generated/quspin.basis.spinful_fermion_basis_1d.html",
        paper="https://scipost.org/10.21468/SciPostPhys.7.2.020")
    (args.output_dir / "quspin_results.json").write_text(json.dumps(external, indent=2) + "\n")
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    source_files = [args.cases, args.julia_results, Path(__file__), Path(__file__).with_name("ed_julia.jl")]
    package_source = Path(julia["source_path"]).parent
    if not package_source.exists() and Path("output/github_release/src").exists():
        # The Julia run may use a WSL path while this comparison runs on Windows.
        package_source = Path("output/github_release/src")
    source_files.extend(package_source / name for name in
                        ("HubbardRenyiDQMC.jl", "exact_diagonalization.jl", "square_lattice.jl"))
    provenance = {"sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files},
                  "installed_distributions": {dist.metadata["Name"]: dist.version
                      for dist in importlib.metadata.distributions()}}
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    with (args.output_dir / "comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise SystemExit("QuSpin comparison failed; inspect comparison.csv")


if __name__ == "__main__":
    main()
