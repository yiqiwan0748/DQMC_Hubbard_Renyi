"""Compare each DQMC engine with an independent QuSpin finite-step trace.

Only reads existing evidence. Both output files must be new.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import tomllib


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_matrix(matrix):
    assert len(matrix) == 4 and all(len(row) == 4 for row in matrix)
    return [[float(value) for value in row] for row in matrix]


def matrix_hash(matrix):
    payload = json.dumps(canonical_matrix(matrix), separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode()).hexdigest()


def require_equal(value, expected, name):
    if not math.isclose(float(value), float(expected), rel_tol=0, abs_tol=1e-12):
        raise ValueError(f"Unmatched {name}: {value} != {expected}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quspin_reference", type=Path)
    parser.add_argument("dqmc_directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    metadata_path = args.output.with_name(args.output.stem + "_metadata.json")
    if args.output.exists() or metadata_path.exists():
        raise SystemExit("Output exists; choose new evidence paths")

    expected_k = [[0.,-1.,-1.,0.],[-1.,0.,0.,-1.],[-1.,0.,0.,-1.],[0.,-1.,-1.,0.]]
    expected_hash = matrix_hash(expected_k)
    reference = json.loads(args.quspin_reference.read_text())
    assert reference["transfer"] == "exp(-dtau*V) @ exp(-dtau*K)"
    assert reference["versions"]["quspin"] == "1.0.1"
    refs = {}
    for case in reference["cases"]:
        assert case["shape"] == [2,2] and case["periodic"] == [False,False]
        for key, value in (("t",1), ("mu",0), ("beta",2), ("dtau",0.05), ("steps",40), ("dimension",256)):
            require_equal(case[key], value, "QuSpin " + key)
        assert canonical_matrix(case["K"]) == expected_k
        assert matrix_hash(case["K"]) == expected_hash
        require_equal(case["energy"], 4 * case["energy_per_site"], "energy normalization")
        interaction = float(case["U"])
        assert interaction not in refs
        refs[interaction] = case
    assert set(refs) == {0.,2.}

    source_paths = [args.quspin_reference]
    for prefix in ("local", "smoqy"):
        path = args.dqmc_directory / (prefix + "_metadata.toml")
        source_paths.append(path)
        metadata = tomllib.loads(path.read_text())
        settings = metadata["settings"]
        assert settings["periodic"] == [False,False]
        assert settings["ensemble"] == "grand canonical"
        assert settings["U_values"] == [0.,2.]
        assert settings["symmetric"] is False and settings["checkerboard"] is False
        assert settings["propagator"] == "exp(-dtau V) exp(-dtau K)"
        for key, value in (("Lx",2), ("Ly",2), ("sites",4), ("t",1), ("mu",0), ("beta",2), ("dtau",0.05)):
            require_equal(settings[key], value, prefix + " " + key)
        assert canonical_matrix(metadata["K_rows"]) == expected_k
        assert matrix_hash(metadata["K_rows"]) == expected_hash

    comparison_path = args.dqmc_directory / "comparison.csv"
    source_paths.append(comparison_path)
    with comparison_path.open(newline="") as handle:
        comparisons = list(csv.DictReader(handle))
    expected_pairs = {(u, observable) for u in (0.,2.) for observable in
        ("energy_per_site", "density", "double_occupancy")}
    observed_pairs = {(float(row["U"]), row["observable"]) for row in comparisons}
    assert observed_pairs == expected_pairs and len(comparisons) == 6

    rows = []
    for comparison in comparisons:
        interaction = float(comparison["U"])
        observable = comparison["observable"]
        exact = refs[interaction][observable]
        assert math.isfinite(exact)
        for prefix, engine in (("local", "HubbardRenyiDQMC"), ("smoqy", "SmoQyDQMC")):
            mean = float(comparison[prefix + "_mean"])
            error = float(comparison[prefix + "_standard_error"])
            assert math.isfinite(mean) and math.isfinite(error) and error >= 0
            residual = mean - exact
            tolerance = max(5 * error, 1e-10)
            z = residual / error if error > 1e-12 else ""
            continuum = refs[interaction].get("continuum_" + observable, "")
            row = dict(engine=engine, U=interaction, observable=observable,
                mean=mean, standard_error=error, quspin_finite_step=exact,
                residual=residual, z=z, absolute_tolerance=tolerance,
                passed=abs(residual) <= tolerance, quspin_continuum=continuum)
            rows.append(row)
            print(json.dumps(row), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    provenance = dict(purpose="Independent finite-step observable reference; no Renyi or performance claim",
        python=sys.version, platform=platform.platform(), K=expected_k,
        K_sha256=expected_hash, quspin_versions=reference["versions"],
        roundoff_floor=1e-10, sampling_error_multiplier=5,
        sources={str(path): sha256(path) for path in source_paths},
        script_sha256=sha256(Path(__file__)), comparisons=len(rows),
        passed=sum(row["passed"] for row in rows))
    with metadata_path.open("x") as handle:
        json.dump(provenance, handle, indent=2)
        handle.write("\n")
    if not all(row["passed"] for row in rows):
        raise SystemExit("A DQMC engine disagrees with the QuSpin finite-step reference")


if __name__ == "__main__":
    main()
