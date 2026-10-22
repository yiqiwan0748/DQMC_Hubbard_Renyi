"""Independent finite-step grand-canonical references from QuSpin matrices.

The exact finite-step trace uses T = exp(-dtau*V) exp(-dtau*K), not a
continuum approximation. Observables are traced against T**(beta/dtau).
"""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys

from ed_quspin import make_operators
import numpy as np
from scipy.linalg import expm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Output exists; choose a new path: " + str(args.output))
    beta, dtau = 2.0, 0.05
    steps = round(beta / dtau)
    rows = []
    for u in (0.0, 2.0):
        case = dict(shape=[2, 2], periodic=[False, False], t=1.0, U=u, mu=0.0)
        k, basis, h, doubles, number = make_operators(case)
        _, _, kinetic, _, _ = make_operators(dict(case, U=0.0))
        h_matrix = h.toarray()
        k_matrix = kinetic.toarray()
        v_matrix = h_matrix - k_matrix
        transfer = expm(-dtau * v_matrix) @ expm(-dtau * k_matrix)
        rho = np.linalg.matrix_power(transfer, steps)
        normalization = np.trace(rho)
        assert np.isfinite(normalization) and normalization > 0
        assert np.allclose(h_matrix, h_matrix.T, atol=1e-14, rtol=0)
        expectation = lambda operator: float(np.trace(operator @ rho) / normalization)
        eigenvalues, vectors = h.eigh()
        weights = np.exp(-beta * (eigenvalues - eigenvalues[0]))
        weights /= weights.sum()
        continuum_energy = float(weights @ eigenvalues)
        continuum_double = float(weights @ (np.abs(vectors).T**2 @ doubles.diagonal()))
        row = dict(U=u, t=1.0, mu=0.0, shape=[2, 2], periodic=[False, False],
            beta=beta, dtau=dtau, steps=steps, dimension=int(basis.Ns),
            K=k.tolist(), energy=expectation(h_matrix),
            energy_per_site=expectation(h_matrix) / len(k),
            density=expectation(number.toarray()),
            double_occupancy=expectation(doubles.toarray()),
            logZ=float(np.log(normalization)),
            continuum_energy=continuum_energy,
            continuum_energy_per_site=continuum_energy / len(k),
            continuum_double_occupancy=continuum_double)
        rows.append(row)
        print(json.dumps(row), flush=True)
    assert abs(rows[0]["energy"] - rows[0]["continuum_energy"]) < 1e-11
    files = [Path(__file__), Path(__file__).with_name("ed_quspin.py")]
    result = dict(hamiltonian="sum(Kij cdag_i c_j)+U sum((nup-.5)(ndn-.5))-mu*N",
        transfer="exp(-dtau*V) @ exp(-dtau*K)", observables="trace(O @ transfer**steps)/trace(transfer**steps)",
        python=sys.version, platform=platform.platform(),
        versions={name: importlib.metadata.version(name) for name in
            ("quspin", "quspin-extensions", "parallel-sparse-tools", "numpy", "scipy")},
        sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
        cases=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
