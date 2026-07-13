#!/usr/bin/env python3
"""Audit whether the PS adjoint cubic can be a primitive Spin(10) 45^3 term.

The Pati-Salam Goldstone-locking derivation used

    W_C = m_C/2 Tr Sigma_C^2 + lambda_C/3 Tr Sigma_C^3

for an SU(4)_C adjoint.  This script checks the next embedding question:
can the cubic be the restriction of a renormalizable Spin(10) invariant built
only from one 45?

Answer: no.  The adjoint invariant polynomials of so(10) have degrees
2,4,5,6,8, so there is no symmetric cubic invariant.  In vector-index form,
the only metric contraction for one antisymmetric matrix A is Tr A^3, and this
vanishes identically because A^T=-A.  The SU(4) cubic is a different invariant,
the fundamental trace/d-symbol of SU(4), and must come either from a larger
Spin(10) field with a cubic invariant (for example a 210-sector contraction) or
from matching to a Pati-Salam EFT below the first GUT-breaking step.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PS_JSON = ROOT / "output" / "ps_goldstone_locking" / "ps_goldstone_locking_summary.json"
OUT = ROOT / "output" / "spin10_cubic_embedding"


def random_antisymmetric_trace_test(samples: int = 200, size: int = 10, seed: int = 20260506) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(samples):
        raw = rng.normal(size=(size, size))
        a = raw - raw.T
        values.append(float(np.trace(a @ a @ a)))
    return {
        "samples": samples,
        "size": size,
        "seed": seed,
        "max_abs_trace_A3": float(np.max(np.abs(values))),
        "mean_abs_trace_A3": float(np.mean(np.abs(values))),
        "sample_values_first_5": values[:5],
    }


def su4_cubic_examples() -> list[dict[str, object]]:
    examples = []
    matrices = {
        "breaking_direction_diag_1_1_1_minus3": np.diag([1.0, 1.0, 1.0, -3.0]),
        "rank_two_diag_1_minus1_0_0": np.diag([1.0, -1.0, 0.0, 0.0]),
        "generic_traceless_diag_4_1_minus2_minus3": np.diag([4.0, 1.0, -2.0, -3.0]),
    }
    for name, h in matrices.items():
        norm = np.sqrt(np.trace(h @ h).real)
        hn = h / norm
        examples.append(
            {
                "name": name,
                "trace_H": float(np.trace(h).real),
                "trace_H2": float(np.trace(h @ h).real),
                "trace_H3": float(np.trace(h @ h @ h).real),
                "trace_normalized_H3": float(np.trace(hn @ hn @ hn).real),
            }
        )
    return examples


def load_ps_lock() -> dict[str, object]:
    payload = json.loads(PS_JSON.read_text(encoding="utf-8"))
    return {
        "eigenvalue_summary": payload["eigenvalue_summary"],
        "all_eigen_residuals_small": payload["all_eigen_residuals_small"],
        "all_locking_bounds_pass": payload["all_locking_bounds_pass"],
    }


def write_report(payload: dict[str, object]) -> None:
    lines: list[str] = []
    lines.append("# Spin(10) cubic embedding audit")
    lines.append("")
    lines.append("No web lookup was used.  This audit tests whether the Pati-Salam")
    lines.append("adjoint cubic used for Goldstone locking can be a primitive renormalizable")
    lines.append("`45^3` invariant of Spin(10).")
    lines.append("")
    lines.append("## No-go for one Spin(10) adjoint")
    lines.append("")
    lines.append("Represent the `45` by a real antisymmetric matrix `A`, so `A^T=-A`.")
    lines.append("Then")
    lines.append("")
    lines.append("```text")
    lines.append("Tr A^3 = Tr[(A^3)^T] = Tr[(A^T)^3] = Tr[-A^3] = - Tr A^3,")
    lines.append("```")
    lines.append("")
    lines.append("hence `Tr A^3=0`.  The antisymmetric structure-constant invariant")
    lines.append("`f_abc Phi^a Phi^b Phi^c` also vanishes for a single commuting chiral")
    lines.append("superfield.  Equivalently, the primitive invariant degrees of `D5=so(10)`")
    lines.append("are `2,4,5,6,8`, with no degree-three symmetric invariant.")
    lines.append("")
    tr = payload["so10_trace_test"]
    lines.append("Numerical antisymmetric-matrix check:")
    lines.append("")
    lines.append("```text")
    lines.append(f"samples = {tr['samples']}")
    lines.append(f"max |Tr A^3| = {tr['max_abs_trace_A3']:.3e}")
    lines.append(f"mean |Tr A^3| = {tr['mean_abs_trace_A3']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("## Why the Pati-Salam cubic is different")
    lines.append("")
    lines.append("The SU(4)_C adjoint has a fundamental-trace cubic/d-symbol.  For example:")
    lines.append("")
    lines.append("| example | Tr H | Tr H^2 | Tr H^3 | normalized Tr H^3 |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in payload["su4_cubic_examples"]:
        lines.append(
            f"| `{row['name']}` | {row['trace_H']:.1f} | {row['trace_H2']:.1f} | "
            f"{row['trace_H3']:.1f} | {row['trace_normalized_H3']:.6f} |"
        )
    lines.append("")
    lines.append("Thus the PS cubic used in")
    lines.append("`W_C = m_C Tr Sigma_C^2/2 + lambda_C Tr Sigma_C^3/3` is not the")
    lines.append("restriction of a primitive one-field `45^3` Spin(10) invariant.")
    lines.append("")
    lines.append("## Consequence for the framework")
    lines.append("")
    lines.append("The PS-stage Goldstone lock remains valid: its Hessian gives")
    ps = payload["ps_goldstone_lock"]
    for sector, row in ps["eigenvalue_summary"].items():
        vals = ",".join(str(v) for v in row["unique_eigenvalues_over_m"])
        lines.append(f"- `{sector}`: multiplicity {row['multiplicity']}, eigenvalue/m_C `{vals}`.")
    lines.append("")
    lines.append("But a full Spin(10) paper has to choose one of two honest completions:")
    lines.append("")
    lines.append("1. derive the same PS cubic from a field with a genuine Spin(10) cubic")
    lines.append("   invariant, with the `210_H^3` contraction the natural next target; or")
    lines.append("2. state the cubic as a renormalizable Pati-Salam EFT term obtained after")
    lines.append("   the first Spin(10)-breaking step.")
    lines.append("")
    lines.append("This converts the previous embedding caveat into a precise no-go and a")
    lines.append("well-defined next calculation, rather than leaving it as vague model-building.")
    (OUT / "spin10_cubic_embedding_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "note": "No web lookup used. 45^3 is not a primitive one-field Spin(10) source of the PS cubic.",
        "input_files": {"ps_goldstone_locking_summary": str(PS_JSON)},
        "so10_adjoint_primitive_invariant_degrees": [2, 4, 5, 6, 8],
        "so10_trace_test": random_antisymmetric_trace_test(),
        "su4_cubic_examples": su4_cubic_examples(),
        "ps_goldstone_lock": load_ps_lock(),
        "verdict": {
            "one_field_45_cubic_available": False,
            "ps_cubic_is_restriction_of_45_cubic": False,
            "next_full_spin10_target": "Compute explicit 210_H^3 or 210-mediated contraction and match its (15,1,1) branch to the PS cubic.",
            "fallback": "State the SU(4)_C cubic as the renormalizable Pati-Salam EFT below the first Spin(10)-breaking step.",
        },
    }
    (OUT / "spin10_cubic_embedding_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_report(payload)

    tr = payload["so10_trace_test"]
    su4 = payload["su4_cubic_examples"][0]
    print("Spin(10) cubic embedding audit")
    print(f"  max |Tr A^3| for SO(10) antisymmetric A: {tr['max_abs_trace_A3']:.3e}")
    print(f"  SU(4) breaking-direction Tr H^3: {su4['trace_H3']:.3f}")
    print("  one-field 45^3 available: False")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
