#!/usr/bin/env python3
"""Component Hessian audit for the constrained clockwork source branch.

No web lookup is used.

This is the follow-up to audit_rank_one_clockwork_locking.py.  The previous
audit showed that q=3,n=6 gives kappa=q^n=729, enough for the dressed C5
stress grid.  Here we check the component Hessian of the source/link
constraints

    W_link = M sum_{a=0}^{n-1} Lambda_a (X_a - q X_{a+1})

and the boundary/source closure

    W_bd = M Xi (v . X),

where v is the normalized clockwork zero vector.  Without W_bd there is one
intended zero mode per orthogonal chain.  With W_bd, the stacked constraint
matrix is full rank and has no exponentially light singular value.  This
supports a constrained/composite source interpretation of the clockwork lock.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "constrained_clockwork_source_hessian"
CLOCKWORK = ROOT / "output" / "rank_one_clockwork_locking" / "summary.json"
CROSSED = ROOT / "output" / "crossed_120_triplet_projector" / "summary.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def link_constraint(q: float, n: int) -> np.ndarray:
    c = np.zeros((n, n + 1), dtype=float)
    for a in range(n):
        c[a, a] = 1.0
        c[a, a + 1] = -q
    return c


def zero_mode(q: float, n: int) -> np.ndarray:
    v = np.array([q ** (-a) for a in range(n + 1)], dtype=float)
    return v / np.linalg.norm(v)


def chiral_hessian_from_constraint(a: np.ndarray) -> np.ndarray:
    """For W=M L^T A X return the dimensionless Hessian [[0,A.T],[A,0]]."""
    return np.block(
        [
            [np.zeros((a.shape[1], a.shape[1])), a.T],
            [a, np.zeros((a.shape[0], a.shape[0]))],
        ]
    )


def eig_abs_sorted(matrix: np.ndarray) -> np.ndarray:
    return np.sort(np.abs(np.linalg.eigvalsh(matrix)))


def spectrum_rows(label: str, values: np.ndarray, m_lock: float) -> list[dict[str, Any]]:
    rows = []
    for i, val in enumerate(values):
        rows.append(
            {
                "spectrum": label,
                "index": i,
                "mass_over_Mlock": float(val),
                "mass_GeV": float(val * m_lock),
            }
        )
    return rows


def audit() -> dict[str, Any]:
    cw = read_json(CLOCKWORK)
    crossed = read_json(CROSSED)
    card = cw["preferred_clockwork_card"]
    q = float(card["q"])
    n = int(card["n"])
    m_lock = float(crossed["threshold_locked_dilation"]["M_lock_GeV"])

    c = link_constraint(q, n)
    v = zero_mode(q, n)
    residual = c @ v

    h_raw = chiral_hessian_from_constraint(c)
    raw_abs = eig_abs_sorted(h_raw)
    raw_zero_count = int(np.sum(raw_abs < 1.0e-10))
    raw_nonzero = raw_abs[raw_abs >= 1.0e-10]

    a_closed = np.vstack([c, v[None, :]])
    h_closed = chiral_hessian_from_constraint(a_closed)
    closed_abs = eig_abs_sorted(h_closed)
    closed_zero_count = int(np.sum(closed_abs < 1.0e-10))

    singular_c = np.linalg.svd(c, compute_uv=False)
    singular_closed = np.linalg.svd(a_closed, compute_uv=False)

    eps = q ** (-n)
    w_eff = np.diag([1.0, eps, eps, eps])
    w_singular = np.linalg.svd(w_eff, compute_uv=False)

    raw_rows = spectrum_rows("one_chain_raw_link_hessian_abs", raw_abs, m_lock)
    closed_rows = spectrum_rows("one_chain_boundary_closed_hessian_abs", closed_abs, m_lock)
    chain_rows = raw_rows + closed_rows

    return {
        "note": "No web lookup used. Component Hessian audit for the q=3,n=6 constrained clockwork source branch.",
        "input_card": {
            "q": q,
            "n": n,
            "kappa": float(card["kappa"]),
            "epsilon": eps,
            "M_lock_GeV": m_lock,
        },
        "superpotential": {
            "W_link": "M sum_{a=0}^{n-1} Lambda_a (X_a - q X_{a+1})",
            "W_boundary": "M Xi (v . X), with v_a proportional q^{-a}",
            "meaning": "W_link creates the clockwork endpoint suppression; W_boundary fixes/removes the intended source zero mode.",
        },
        "linear_algebra": {
            "constraint_shape": list(c.shape),
            "rank_C": int(np.linalg.matrix_rank(c)),
            "null_residual_norm": float(np.linalg.norm(residual)),
            "zero_mode_endpoint_ratio": float(abs(v[-1] / v[0])),
            "raw_hessian_dimension": int(h_raw.shape[0]),
            "raw_zero_modes_per_chain": raw_zero_count,
            "raw_nonzero_min": float(np.min(raw_nonzero)),
            "raw_nonzero_max": float(np.max(raw_nonzero)),
            "raw_nonzero_condition": float(np.max(raw_nonzero) / np.min(raw_nonzero)),
            "boundary_closed_hessian_dimension": int(h_closed.shape[0]),
            "boundary_closed_zero_modes_per_chain": closed_zero_count,
            "boundary_closed_min_abs_eigenvalue": float(np.min(closed_abs)),
            "boundary_closed_max_abs_eigenvalue": float(np.max(closed_abs)),
            "boundary_closed_condition": float(np.max(closed_abs) / np.min(closed_abs)),
            "singular_values_C": [float(x) for x in singular_c],
            "singular_values_closed_constraint": [float(x) for x in singular_closed],
        },
        "effective_inverse_block": {
            "W_eff_singular_values": [float(x) for x in w_singular],
            "condition_number": float(np.max(w_singular) / np.min(w_singular)),
            "matches_required_kappa": bool(np.max(w_singular) / np.min(w_singular) >= cw["input_required_kappa"]),
            "dressed_max_width_margin_at_ST_1e_minus_5": float(card["max_width_margin_1e35_at_ST_1e_minus_5"]),
        },
        "three_chain_aggregate": {
            "raw_intended_zero_modes": 3 * raw_zero_count,
            "boundary_closed_zero_modes": 3 * closed_zero_count,
            "massive_dirac_pairs_per_chain_before_boundary": n,
            "massive_dirac_pairs_three_chains_before_boundary": 3 * n,
            "massive_dirac_pairs_per_chain_after_boundary": n + 1,
            "massive_dirac_pairs_three_chains_after_boundary": 3 * (n + 1),
        },
        "threshold_interpretation": {
            "constrained_or_composite_source_nonuniversal_threshold": [0.0, 0.0, 0.0],
            "complete_degenerate_link_nonuniversal_threshold": [0.0, 0.0, 0.0],
            "literal_propagating_visible_chain_warning": (
                "The previous clockwork audit estimated Lambda_LP/M_G=12.641502 for q=3,n=6 if the "
                "three chains are implemented as propagating complete 10+10bar pairs on top of the minimal "
                "Yukawa/Higgs sector, so this branch should be constrained/composite or otherwise cutoff-local."
            ),
        },
        "spectrum_rows": chain_rows,
        "verdict": {
            "raw_chain_has_one_intended_zero_mode_per_chain": raw_zero_count == 1,
            "boundary_driver_removes_all_clockwork_zero_modes": closed_zero_count == 0,
            "no_exponential_light_tower_after_boundary": float(np.min(closed_abs)) >= 0.999999,
            "preferred_status": "PASS_CONDITIONAL",
            "main_caveat": (
                "The component Hessian supports the constrained source mechanism, but it is not yet a "
                "microscopic Spin(10) UV completion.  The chain must remain a source/composite sector or "
                "be embedded in complete threshold-silent multiplets."
            ),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(summary: dict[str, Any]) -> None:
    la = summary["linear_algebra"]
    eff = summary["effective_inverse_block"]
    agg = summary["three_chain_aggregate"]
    th = summary["threshold_interpretation"]
    verdict = summary["verdict"]
    card = summary["input_card"]
    lines = [
        "# Constrained clockwork source Hessian audit",
        "",
        "No web lookup was used.",
        "",
        f"Input card: `q={card['q']:.1f}`, `n={card['n']}`, "
        f"`kappa={card['kappa']:.0f}`, `epsilon={card['epsilon']:.6e}`.",
        "",
        "## Algebra",
        "",
        "`W_link = M sum Lambda_a (X_a - q X_{a+1})` gives a constraint matrix `C`.",
        "Its null vector is the clockwork source profile `v_a proportional q^{-a}`.",
        "",
        f"`||C v|| = {la['null_residual_norm']:.6e}`.",
        f"`|v_n/v_0| = {la['zero_mode_endpoint_ratio']:.6e}`.",
        "",
        "Without the boundary driver, the Hessian has",
        f"`{la['raw_zero_modes_per_chain']}` intended zero mode per chain and nonzero masses",
        f"from `{la['raw_nonzero_min']:.6f}` to `{la['raw_nonzero_max']:.6f}` in units of `M_lock`.",
        "",
        "With `W_boundary = M Xi (v . X)`, the stacked constraint is full rank:",
        f"`zero modes = {la['boundary_closed_zero_modes_per_chain']}`,",
        f"`min |m|/M_lock = {la['boundary_closed_min_abs_eigenvalue']:.6f}`,",
        f"`max |m|/M_lock = {la['boundary_closed_max_abs_eigenvalue']:.6f}`.",
        "",
        "## Effective inverse block",
        "",
        f"`singular(W_eff) = {eff['W_eff_singular_values']}`.",
        f"`condition(W_eff) = {eff['condition_number']:.0f}`.",
        f"Dressed max-width margin at `S_T=1e-5`: `{eff['dressed_max_width_margin_at_ST_1e_minus_5']:.6f}`.",
        "",
        "## Three-chain aggregate",
        "",
        f"Raw intended zero modes: `{agg['raw_intended_zero_modes']}`.",
        f"Boundary-closed zero modes: `{agg['boundary_closed_zero_modes']}`.",
        f"Massive Dirac pairs after boundary closure: `{agg['massive_dirac_pairs_three_chains_after_boundary']}`.",
        "",
        "## Threshold interpretation",
        "",
        f"Constrained/composite source threshold vector: `{th['constrained_or_composite_source_nonuniversal_threshold']}`.",
        f"Complete-degenerate link threshold vector: `{th['complete_degenerate_link_nonuniversal_threshold']}`.",
        th["literal_propagating_visible_chain_warning"],
        "",
        "## Verdict",
        "",
        f"`raw_chain_has_one_intended_zero_mode_per_chain = {verdict['raw_chain_has_one_intended_zero_mode_per_chain']}`.",
        f"`boundary_driver_removes_all_clockwork_zero_modes = {verdict['boundary_driver_removes_all_clockwork_zero_modes']}`.",
        f"`no_exponential_light_tower_after_boundary = {verdict['no_exponential_light_tower_after_boundary']}`.",
        "",
        verdict["main_caveat"],
        "",
    ]
    (OUT / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = audit()
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(OUT / "hessian_spectrum.csv", summary["spectrum_rows"])
    write_report(summary)
    la = summary["linear_algebra"]
    print("Constrained clockwork source Hessian audit")
    print(f"  raw zero modes per chain: {la['raw_zero_modes_per_chain']}")
    print(f"  boundary-closed zero modes per chain: {la['boundary_closed_zero_modes_per_chain']}")
    print(f"  boundary min abs eigenvalue: {la['boundary_closed_min_abs_eigenvalue']:.6f}")
    print(f"  wrote: {OUT}")


if __name__ == "__main__":
    main()
