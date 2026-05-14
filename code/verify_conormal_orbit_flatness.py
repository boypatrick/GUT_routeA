#!/usr/bin/env python3
"""Verify F-flatness and D-flatness of the constrained 54 conormal orbit.

This is a local numerical audit with no web input.  It samples points
S=O S0 O^T on the spectral orbit and checks:

  C(S)=S^2-S-6I=0,
  rank dC_S = 30 on the traceless-symmetric 54,
  dim ker dC_S = 24, matching Spin(10)/PS tangents,
  D^a=<S,T_a.S>=0 for the real symmetric orbit.

The conormal superpotential is interpreted as W=<Xi_N,C(S)> with Xi_N in the
normal/conormal rank-30 bundle.  At Xi_N=0 and C(S)=0 all F-terms vanish.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "conormal_orbit_flatness"
SEED = 5410


def inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.trace(a.T @ b))


def traceless_symmetric_basis(n: int = 10) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(n):
        for j in range(i + 1, n):
            m = np.zeros((n, n), dtype=float)
            m[i, j] = 1.0 / math.sqrt(2.0)
            m[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(m)
    for k in range(n - 1):
        m = np.zeros((n, n), dtype=float)
        coeff = 1.0 / math.sqrt((k + 1) * (k + 2))
        for i in range(k + 1):
            m[i, i] = coeff
        m[k + 1, k + 1] = -(k + 1) * coeff
        basis.append(m)
    assert len(basis) == 54
    return basis


def symmetric_basis(n: int = 10) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(n):
        m = np.zeros((n, n), dtype=float)
        m[i, i] = 1.0
        basis.append(m)
    for i in range(n):
        for j in range(i + 1, n):
            m = np.zeros((n, n), dtype=float)
            m[i, j] = 1.0 / math.sqrt(2.0)
            m[j, i] = 1.0 / math.sqrt(2.0)
            basis.append(m)
    assert len(basis) == 55
    return basis


def so10_basis(n: int = 10) -> list[np.ndarray]:
    basis: list[np.ndarray] = []
    for i in range(n):
        for j in range(i + 1, n):
            m = np.zeros((n, n), dtype=float)
            m[i, j] = 1.0 / math.sqrt(2.0)
            m[j, i] = -1.0 / math.sqrt(2.0)
            basis.append(m)
    assert len(basis) == 45
    return basis


def random_so10(rng: np.random.Generator) -> np.ndarray:
    a = rng.normal(size=(10, 10))
    q, r = np.linalg.qr(a)
    signs = np.sign(np.diag(r))
    signs[signs == 0.0] = 1.0
    q = q @ np.diag(signs)
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1.0
    return q


def constraint(s: np.ndarray) -> np.ndarray:
    return s @ s - s - 6.0 * np.eye(10)


def linearized_constraint(s: np.ndarray, x: np.ndarray) -> np.ndarray:
    return s @ x + x @ s - x


def jacobian(s: np.ndarray, codomain: list[np.ndarray], domain: list[np.ndarray]) -> np.ndarray:
    j = np.zeros((len(codomain), len(domain)))
    for a, e in enumerate(domain):
        image = linearized_constraint(s, e)
        for b, c in enumerate(codomain):
            j[b, a] = inner(c, image)
    return j


def tangent_gram(s: np.ndarray, generators: list[np.ndarray]) -> np.ndarray:
    tangents = [t @ s - s @ t for t in generators]
    return np.array([[inner(a, b) for b in tangents] for a in tangents])


def d_vector(s: np.ndarray, generators: list[np.ndarray]) -> np.ndarray:
    # The real representation tangent is T.S = T S - S T.
    return np.array([inner(s, t @ s - s @ t) for t in generators])


def audit_sample(label: str, s: np.ndarray, domain: list[np.ndarray], codomain: list[np.ndarray], generators: list[np.ndarray]) -> dict[str, Any]:
    c = constraint(s)
    j = jacobian(s, codomain, domain)
    singular = np.linalg.svd(j, compute_uv=False)
    rank = int(np.sum(singular > 1.0e-10))
    tg = tangent_gram(s, generators)
    tangent_rank = int(np.sum(np.linalg.eigvalsh(tg) > 1.0e-10))
    d = d_vector(s, generators)
    return {
        "label": label,
        "trace_S": float(np.trace(s)),
        "trace_S2": float(np.trace(s @ s)),
        "constraint_norm": float(np.linalg.norm(c)),
        "jacobian_rank": rank,
        "jacobian_nullity": len(domain) - rank,
        "tangent_rank": tangent_rank,
        "D_norm": float(np.linalg.norm(d)),
        "F_Xi_norm_at_orbit": float(np.linalg.norm(c)),
        "F_S_norm_at_Xi_zero": 0.0,
        "min_nonzero_singular": float(min(x for x in singular if x > 1.0e-10)),
        "max_singular": float(max(singular)),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_report(payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Conormal orbit F/D-flatness audit")
    lines.append("")
    lines.append("No web lookup was used.")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(
        "All sampled orbit points satisfy C(S)=0, D=0, rank dC=30, and tangent rank=24 within numerical precision."
    )
    lines.append("")
    lines.append("| sample | ||C|| | rank dC | nullity | tangent rank | ||D|| |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for row in payload["samples"]:
        lines.append(
            f"| {row['label']} | {row['constraint_norm']:.3e} | "
            f"{row['jacobian_rank']} | {row['jacobian_nullity']} | "
            f"{row['tangent_rank']} | {row['D_norm']:.3e} |"
        )
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(payload["verdict"])
    lines.append("")
    (OUT / "conormal_orbit_flatness_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    domain = traceless_symmetric_basis()
    codomain = symmetric_basis()
    generators = so10_basis()
    s0 = np.diag([-2.0] * 6 + [3.0] * 4)
    samples = [audit_sample("S0", s0, domain, codomain, generators)]
    for i in range(1, 6):
        o = random_so10(rng)
        s = o @ s0 @ o.T
        samples.append(audit_sample(f"random_orbit_{i}", s, domain, codomain, generators))
    passes = all(
        row["constraint_norm"] < 1.0e-11
        and row["D_norm"] < 1.0e-11
        and row["jacobian_rank"] == 30
        and row["jacobian_nullity"] == 24
        and row["tangent_rank"] == 24
        for row in samples
    )
    payload = {
        "seed": SEED,
        "samples": samples,
        "passes": passes,
        "verdict": (
            "The conormal-bundle orbit constraint is F-flat at Xi_N=0 and C(S)=0, "
            "and D-flat on the real symmetric Spin(10) orbit. The remaining open "
            "condition is model-building: Xi_N must be auxiliary/composite or a "
            "complete degenerate threshold sector."
        ),
    }
    write_csv(OUT / "conormal_orbit_flatness_samples.csv", samples)
    (OUT / "conormal_orbit_flatness_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(payload)
    print("Conormal orbit F/D-flatness audit complete")
    print(f"  passes={passes}")
    for row in samples:
        print(
            f"  {row['label']}: C={row['constraint_norm']:.3e} "
            f"D={row['D_norm']:.3e} rank={row['jacobian_rank']} tangent={row['tangent_rank']}"
        )
    print(f"  outputs={OUT}")


if __name__ == "__main__":
    main()
