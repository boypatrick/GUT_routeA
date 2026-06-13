#!/usr/bin/env python3
"""Audit 1b covariant-rank contract.

This is not a flavor fit.  It is the design check that prevents Audit 1b from
becoming a hollow test.  The local model is the shared-spurion Veronese-sector
map

    T_a = alpha_a Phi + beta_a H(Phi),  a=1,...,k,

where Phi is a binary quartic and H(Phi) is its Hessian, both represented in
the coefficient convention

    a x^4 + 4 b x^3 y + 6 c x^2 y^2 + 4 d x y^3 + e y^4.

For k tensors, the ambient Veronese-sector space has dimension 5k over C.  A
finite-difference complex Jacobian check at generic points verifies

    rank = 2k + 3,
    constraints = 5k - (2k+3) = 3(k-1) over C.

This validates the corrected Audit-1b design: a literature two-tensor (h,f)
test with unknown U(3) basis is hollow, while an internal companion test with
three or more aligned tensors has positive net constraints after the contact
direction fixes the basis up to O(3).
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit1b"

SOURCES = {
    "audit1b_builder": ROOT / "code" / "audit1b_covariant_rank_contract.py",
}


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_digest(obj: Any) -> str:
    encoded = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def manifest() -> list[dict[str, Any]]:
    rows = []
    for label, path in SOURCES.items():
        rows.append(
            {
                "label": label,
                "path": str(path.relative_to(ROOT)),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else None,
                "sha256": sha256(path),
            }
        )
    return rows


def convolve(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    out = np.zeros(len(p) + len(q) - 1, dtype=complex)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def dx(poly: np.ndarray, degree: int) -> np.ndarray:
    return np.array([(degree - m) * poly[m] for m in range(degree)], dtype=complex)


def dy(poly: np.ndarray, degree: int) -> np.ndarray:
    return np.array([(m + 1) * poly[m + 1] for m in range(degree)], dtype=complex)


def hessian_quartic(phi: np.ndarray) -> np.ndarray:
    # Convert invariant coefficients (a,b,c,d,e) to ordinary monomial
    # coefficients for x^(4-m)y^m.
    a, b, c, d, e = phi
    f = np.array([a, 4.0 * b, 6.0 * c, 4.0 * d, e], dtype=complex)
    f_xx = dx(dx(f, 4), 3)
    f_xy = dy(dx(f, 4), 3)
    f_yy = dy(dy(f, 4), 3)
    ordinary = convolve(f_xx, f_yy) - convolve(f_xy, f_xy)
    # Return to the same binary-quartic convention.
    return np.array(
        [
            ordinary[0],
            ordinary[1] / 4.0,
            ordinary[2] / 6.0,
            ordinary[3] / 4.0,
            ordinary[4],
        ],
        dtype=complex,
    )


def shared_spurion_map(params: np.ndarray, k: int) -> np.ndarray:
    phi = params[:5]
    h_phi = hessian_quartic(phi)
    rows = []
    offset = 5
    for _ in range(k):
        alpha = params[offset]
        beta = params[offset + 1]
        offset += 2
        rows.extend(alpha * phi + beta * h_phi)
    return np.array(rows, dtype=complex)


def jacobian_rank(k: int, seed: int, eps: float = 1.0e-7) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    params = (rng.normal(size=5 + 2 * k) + 1j * rng.normal(size=5 + 2 * k)) * 0.2
    base = shared_spurion_map(params, k)
    jac = np.zeros((5 * k, 5 + 2 * k), dtype=complex)
    for col in range(len(params)):
        step = np.zeros_like(params)
        step[col] = eps
        jac[:, col] = (shared_spurion_map(params + step, k) - base) / eps
    singular_values = np.linalg.svd(jac, compute_uv=False)
    tolerance = 1.0e-7 * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    expected_rank = 2 * k + 3
    ambient_dim = 5 * k
    raw_parameter_dim = 5 + 2 * k
    return {
        "k_tensor_count": k,
        "ambient_complex_dimension": ambient_dim,
        "raw_parameter_complex_dimension": raw_parameter_dim,
        "rank_complex": rank,
        "expected_rank_complex": expected_rank,
        "fiber_complex_dimension": raw_parameter_dim - rank,
        "expected_fiber_complex_dimension": 2,
        "constraint_complex_dimension": ambient_dim - rank,
        "expected_constraint_complex_dimension": 3 * (k - 1),
        "rank_matches": rank == expected_rank,
        "constraint_matches": ambient_dim - rank == 3 * (k - 1),
        "singular_values": [float(x) for x in singular_values],
        "rank_tolerance": float(tolerance),
        "seed": seed,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [jacobian_rank(k, seed=20260613 + k) for k in range(1, 6)]

    literature_hf = {
        "tensor_count_k": 2,
        "visible_complex_constraints": 3,
        "visible_real_constraints": 6,
        "unknown_basis_real_dimension": 9,
        "net_real_constraints_after_unknown_U3_basis": -3,
        "verdict": "hollow_without_fixed_basis",
        "interpretation": "Do not use a literature (h,f) two-tensor test with unknown U(3) orientation as a pass/fail result.",
    }
    internal_tests = [
        {
            "tensor_count_k": 3,
            "visible_complex_constraints": 6,
            "visible_real_constraints": 12,
            "residual_basis_real_dimension_after_K_contact": 3,
            "net_real_constraints": 9,
        },
        {
            "tensor_count_k": 4,
            "visible_complex_constraints": 9,
            "visible_real_constraints": 18,
            "residual_basis_real_dimension_after_K_contact": 3,
            "net_real_constraints": 15,
        },
    ]

    card: dict[str, Any] = {
        "audit": "audit1b_covariant_rank_contract",
        "status": "rank/design contract only; not a flavor fit",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": manifest(),
        "map_definition": {
            "ambient": "k Veronese-sector binary quartics, dimension 5k over C",
            "shared_spurion_map": "T_a = alpha_a Phi + beta_a H(Phi)",
            "coefficient_convention": "a x^4 + 4 b x^3 y + 6 c x^2 y^2 + 4 d x y^3 + e y^4",
            "hessian_convention": "H(Phi) = f_xx f_yy - f_xy^2, returned in the same binary-quartic coefficient convention",
        },
        "rank_rows": rows,
        "all_rank_checks_pass": all(row["rank_matches"] and row["constraint_matches"] for row in rows),
        "design_implications": {
            "empty_test_to_avoid": literature_hf,
            "internal_companion_tests_with_teeth": internal_tests,
            "required_reporting": [
                "Report residuals tensor-by-tensor rather than only a pooled residual.",
                "Each added aligned tensor should add 3 complex / 6 real constraints before residual basis subtraction.",
                "Only use literature matrices as a non-binding warm start unless their family basis is fixed relative to K_tr.",
            ],
        },
        "publication_boundary": {
            "not_claimed": [
                "global flavor fit",
                "projection residual for real benchmark matrices",
                "RGE or target-table freshness",
            ],
            "claim": "Audit 1b rank contract fixes which covariant-subspace projection tests are meaningful.",
        },
    }
    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)

    json_path = OUT / "covariant_rank_contract.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Audit 1b Covariant-Rank Contract",
        "",
        "This is a design/rank check, not a flavor fit.",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- all rank checks pass: `{card['all_rank_checks_pass']}`",
        "",
        "## Rank Rows",
        "",
        "| k | ambient C dim | rank | expected | constraints C | expected | fiber C |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['k_tensor_count']} | {row['ambient_complex_dimension']} | "
            f"{row['rank_complex']} | {row['expected_rank_complex']} | "
            f"{row['constraint_complex_dimension']} | {row['expected_constraint_complex_dimension']} | "
            f"{row['fiber_complex_dimension']} |"
        )
    lines += [
        "",
        "## Design Implications",
        "",
        "- Do not run the literature `(h,f)` two-tensor test as a pass/fail audit when the family basis is unknown; its net real constraints are negative after scanning a full `U(3)` orientation.",
        "- Use the internal companion test with at least three aligned family tensors after the `K_tr` contact direction fixes the basis up to a residual `O(3)`.",
        "- Report tensor-by-tensor residual increments; a pooled residual can hide whether the expected 6 real constraints per added tensor are doing work.",
        "",
    ]
    md_path = OUT / "covariant_rank_contract.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {md_path.relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")
    print(f"all_rank_checks_pass = {card['all_rank_checks_pass']}")


if __name__ == "__main__":
    main()
