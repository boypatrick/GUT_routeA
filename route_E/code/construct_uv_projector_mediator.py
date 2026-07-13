#!/usr/bin/env python3
"""Renormalizable mediator realization of the 54_H spectral projector.

This script upgrades the effective X_(6,2,2) projector lift to an explicit
renormalizable Spin(10) superpotential with two heavy adjoint mediators B and C:

    W/MG = 1/2 A K0 A
         + R B C
         + eta A (F54 - 2) B
         + eta A (F54 + 4/3) C,

where K0 = mu + a54 F54 + b210 D210.  All terms are bilinear masses or cubic
45-45-54 / 45-45-210 couplings after inserting the 54_H and 210_H vevs.

Tree-level elimination of B,C gives the Schur-complement operator

    K_eff = K0 - (2 eta^2/R) (F54 - 2)(F54 + 4/3)
          = K0 + [50 eta^2/(9R)] P_X.

For finite mediator mass the physical light eigenvalue differs from the Schur
limit because the Kähler metric is also corrected.  We therefore solve the
full 3-by-3 mass matrix in each Clebsch sector and fit the light eigenvalues
directly.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import root


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "output" / "clebsch_210" / "clebsch_210_lifting_summary.json"
OUT = ROOT / "output" / "uv_projector_mediator"

STATE_DATA = {
    "Sigma_L": {
        "label": "(1,3,1)",
        "F54": 2.0,
        "D210": 1.0,
        "b_full45_block": [0.0, 2.0, 0.0],
    },
    "Sigma_R": {
        "label": "(1,1,3)",
        "F54": 2.0,
        "D210": -1.0,
        "b_full45_block": [6.0 / 5.0, 0.0, 0.0],
    },
    "Sigma8_block": {
        "label": "(15,1,1)",
        "F54": -4.0 / 3.0,
        "D210": 0.0,
        # Full (15,1,1) chiral block: octet + color triplet pair + singlet.
        # The octet-only threshold used in the previous minimal spectrum is
        # recovered after the non-octet Goldstone/lifting sector is specified.
        "b_full45_block": [8.0 / 5.0, 0.0, 4.0],
    },
    "X_622": {
        "label": "(6,2,2)",
        "F54": 1.0 / 3.0,
        "D210": 0.0,
        "b_full45_block": [26.0 / 5.0, 6.0, 4.0],
    },
}

PROJECTOR = np.eye(3) - np.ones((3, 3)) / 3.0


def mass_matrix(mu: float, a54: float, b210: float, eta: float, r_med: float, f54: float, d210: float) -> np.ndarray:
    k0 = mu + a54 * f54 + b210 * d210
    u = eta * (f54 - 2.0)
    v = eta * (f54 + 4.0 / 3.0)
    return np.array(
        [
            [k0, u, v],
            [u, 0.0, r_med],
            [v, r_med, 0.0],
        ],
        dtype=float,
    )


def light_eigenvalue(evals: np.ndarray, target: float) -> float:
    return float(evals[np.argmin(np.abs(evals - target))])


def solve_finite_r(targets: dict[str, float], r_med: float, initial: dict[str, float]) -> dict[str, object]:
    eta0 = math.sqrt(r_med * 9.0 * initial["delta_X_projector"] / 50.0)
    x0 = np.array([initial["mu"], initial["a54"], initial["b210"], eta0], dtype=float)

    def residual(x: np.ndarray) -> np.ndarray:
        mu, a54, b210, eta = x
        vals = []
        for name, target in targets.items():
            data = STATE_DATA[name]
            evals = np.linalg.eigvalsh(
                mass_matrix(mu, a54, b210, eta, r_med, data["F54"], data["D210"])
            )
            vals.append(light_eigenvalue(evals, target) - target)
        return np.array(vals, dtype=float)

    sol = root(residual, x0, method="hybr")
    mu, a54, b210, eta = [float(v) for v in sol.x]
    state_rows = {}
    for name, target in targets.items():
        data = STATE_DATA[name]
        mat = mass_matrix(mu, a54, b210, eta, r_med, data["F54"], data["D210"])
        evals = np.linalg.eigvalsh(mat)
        light = light_eigenvalue(evals, target)
        heavy = [float(abs(v)) for v in evals if abs(v - light) > 1.0e-7]
        state_rows[name] = {
            "label": data["label"],
            "F54": data["F54"],
            "D210": data["D210"],
            "target_kappa": target,
            "matrix": mat.tolist(),
            "eigenvalues_signed": [float(v) for v in evals],
            "light_kappa": light,
            "heavy_kappas_abs": heavy,
            "light_abs_error": abs(light - target),
        }

    return {
        "r_med": r_med,
        "success": bool(sol.success),
        "message": sol.message,
        "parameters": {
            "mu": mu,
            "a54": a54,
            "b210": b210,
            "eta": eta,
            "eta_squared": eta * eta,
            "schur_delta_X": 50.0 * eta * eta / (9.0 * r_med),
        },
        "residual_vector": residual(sol.x).tolist(),
        "residual_l2": float(np.linalg.norm(residual(sol.x))),
        "states": state_rows,
    }


def heavy_threshold_residual(solution: dict[str, object]) -> dict[str, object]:
    delta = np.zeros(3)
    rows = {}
    for name, state in solution["states"].items():
        b = np.array(STATE_DATA[name]["b_full45_block"], dtype=float)
        logs = [-math.log(kappa) for kappa in state["heavy_kappas_abs"]]
        contrib = b * sum(logs) / (2.0 * math.pi)
        delta += contrib
        rows[name] = {
            "b_full45_block": b.tolist(),
            "heavy_kappas_abs": state["heavy_kappas_abs"],
            "log_MG_over_M_sum": sum(logs),
            "threshold_contribution": contrib.tolist(),
        }
    projected = PROJECTOR @ delta
    return {
        "block_rows": rows,
        "delta_full": delta.tolist(),
        "projected_delta": projected.tolist(),
        "projected_l2": float(np.linalg.norm(projected)),
    }


def build_payload() -> dict[str, object]:
    source = json.loads(IN.read_text(encoding="utf-8"))
    target = source["scenarios"]["target"]
    initial = source["scenarios"]["projector_lifted_solution"]["coefficients"]
    mg = target["MG_GeV"]
    targets = {
        "Sigma_L": target["kappa_Sigma3"],
        "Sigma_R": 1.0,
        "Sigma8_block": target["kappa_Sigma8"],
        "X_622": 1.0,
    }

    benchmark_r = 50.0
    benchmark = solve_finite_r(targets, benchmark_r, initial)
    benchmark["heavy_threshold"] = heavy_threshold_residual(benchmark)
    benchmark["masses_GeV"] = {
        name: {
            "light_mass_GeV": row["light_kappa"] * mg,
            "heavy_masses_GeV": [k * mg for k in row["heavy_kappas_abs"]],
        }
        for name, row in benchmark["states"].items()
    }

    decoupling_scan = []
    for r_med in [5.0, 10.0, 20.0, 50.0, 100.0]:
        sol = solve_finite_r(targets, r_med, initial)
        thresh = heavy_threshold_residual(sol)
        decoupling_scan.append(
            {
                "r_med": r_med,
                "eta": sol["parameters"]["eta"],
                "schur_delta_X": sol["parameters"]["schur_delta_X"],
                "light_fit_residual_l2": sol["residual_l2"],
                "heavy_threshold_projected_l2": thresh["projected_l2"],
            }
        )

    return {
        "input_summary": str(IN),
        "note": "No web lookup used.  The UV mediator sector uses only renormalizable mass and cubic Spin(10)-invariant couplings.",
        "superpotential": {
            "formula": "W/MG = 1/2 A K0 A + R B C + eta A(F54-2)B + eta A(F54+4/3)C",
            "K0": "K0 = mu + a54 F54 + b210 D210",
            "fields": {
                "A": "light 45_Sigma adjoint",
                "B": "heavy 45 mediator",
                "C": "heavy 45 mediator",
            },
            "schur_complement": "K_eff = K0 - (2 eta^2/R)(F54-2)(F54+4/3) = K0 + [50 eta^2/(9R)] P_X",
        },
        "MG_GeV": mg,
        "targets": targets,
        "benchmark": benchmark,
        "decoupling_scan": decoupling_scan,
    }


def write_report(payload: dict[str, object]) -> None:
    bench = payload["benchmark"]
    pars = bench["parameters"]
    thresh = bench["heavy_threshold"]
    lines: list[str] = []
    lines.append("# Renormalizable 45-mediator realization of the projector lift")
    lines.append("")
    lines.append("No web lookup was used.  The construction uses two heavy adjoint mediators")
    lines.append("`B_45` and `C_45` and only renormalizable mass/cubic Spin(10) invariants.")
    lines.append("")
    lines.append("## Superpotential")
    lines.append("")
    lines.append("```text")
    lines.append("W/MG = 1/2 A K0 A")
    lines.append("     + R B C")
    lines.append("     + eta A (F54 - 2) B")
    lines.append("     + eta A (F54 + 4/3) C")
    lines.append("K0 = mu + a54 F54 + b210 D210")
    lines.append("```")
    lines.append("")
    lines.append("Here `A` is the light `45_Sigma`, while `B,C` are heavy `45` mediators.")
    lines.append("The shifted factors are ordinary bilinear masses plus `45-45-54` cubic")
    lines.append("couplings after the `54_H` vev is inserted.")
    lines.append("")
    lines.append("Tree-level elimination gives")
    lines.append("")
    lines.append("```text")
    lines.append("K_eff = K0 - (2 eta^2/R)(F54 - 2)(F54 + 4/3)")
    lines.append("      = K0 + [50 eta^2/(9R)] P_X.")
    lines.append("```")
    lines.append("")
    lines.append("## Finite-mediator benchmark")
    lines.append("")
    lines.append("The benchmark solves the full 3-by-3 physical mass matrix, not only the")
    lines.append("Schur-complement limit.")
    lines.append("")
    lines.append("```text")
    lines.append(f"R = M_med/MG = {bench['r_med']:.6g}")
    lines.append(f"mu   = {pars['mu']:.9f}")
    lines.append(f"a54  = {pars['a54']:.9f}")
    lines.append(f"b210 = {pars['b210']:.9f}")
    lines.append(f"eta  = {pars['eta']:.9f}")
    lines.append(f"Schur delta_X = {pars['schur_delta_X']:.9f}")
    lines.append(f"light fit residual_l2 = {bench['residual_l2']:.3e}")
    lines.append("```")
    lines.append("")
    lines.append("| sector | light kappa | target kappa | heavy kappas |")
    lines.append("|---|---:|---:|---:|")
    for name, state in bench["states"].items():
        heavy = ", ".join(f"{v:.6f}" for v in state["heavy_kappas_abs"])
        lines.append(
            f"| `{name}` | {state['light_kappa']:.9f} | {state['target_kappa']:.9f} | {heavy} |"
        )
    lines.append("")
    lines.append("The light spectrum is therefore")
    lines.append("")
    lines.append("```text")
    for name, masses in bench["masses_GeV"].items():
        lines.append(f"{name}: {masses['light_mass_GeV']:.6e} GeV")
    lines.append("```")
    lines.append("")
    lines.append("## Heavy-threshold check")
    lines.append("")
    lines.append("Using complete `45` block beta vectors for the two heavy mediator eigenstates")
    lines.append("gives")
    lines.append("")
    lines.append("```text")
    lines.append(f"projected heavy-threshold residual_l2 = {thresh['projected_l2']:.6e}")
    lines.append(f"projected vector = {thresh['projected_delta']}")
    lines.append("```")
    lines.append("")
    lines.append("The residual decreases in the decoupling limit:")
    lines.append("")
    lines.append("| R | eta | Schur delta_X | heavy residual |")
    lines.append("|---:|---:|---:|---:|")
    for row in payload["decoupling_scan"]:
        lines.append(
            f"| {row['r_med']:.0f} | {row['eta']:.6f} | {row['schur_delta_X']:.6f} | {row['heavy_threshold_projected_l2']:.6e} |"
        )
    lines.append("")
    lines.append("Thus the projector sector has a renormalizable UV realization.  The remaining")
    lines.append("paper-level task is to include these heavy mediator thresholds in the full")
    lines.append("two-loop fit, or to push `R` high enough that their non-universal contribution")
    lines.append("is below the desired matching tolerance.")
    (OUT / "uv_projector_mediator_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    (OUT / "uv_projector_mediator_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    write_report(payload)
    bench = payload["benchmark"]
    pars = bench["parameters"]
    print("Renormalizable 45-mediator projector realization")
    print(f"  R={bench['r_med']:.6g} eta={pars['eta']:.9f}")
    print(f"  mu={pars['mu']:.9f} a54={pars['a54']:.9f} b210={pars['b210']:.9f}")
    print(f"  light residual={bench['residual_l2']:.3e}")
    print(f"  heavy threshold residual={bench['heavy_threshold']['projected_l2']:.6e}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
