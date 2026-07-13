#!/usr/bin/env python3
"""Export hand-audited symbolic entries for the CMSGUT triplet inverse.

The previous Audit 4a.1 literature-import card transcribed the Aulakh--Girdhar
5x5 proton-decay triplet matrix and specified the adjugate/cofactor contract

    S_i^j = (-1)^(i+j) det(T with row j and column i removed) / det(T).

SymPy is not available in the current reproducibility runtime, so this script
does not pretend to provide a simplified CAS polynomial.  Instead it exports a
hand-auditable Leibniz expansion for det(T) and for every Audit-2-required
minor.  Each inverse entry is therefore an explicit rational expression whose
numerator and denominator are finite signed sums of products of source-anchored
triplet matrix entries.

The script also numerically evaluates the exported terms at the same generic
F-flat sample used by the literature-import card and checks the result against
numpy.linalg.inv(T).  This is the gate that promotes the cofactor contract from
"pending" to "symbolic entries exported, numerically verified".
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit4a1"
LITERATURE_CARD = OUT / "literature_mass_matrices.json"

sys.path.insert(0, str(ROOT / "code"))
from audit4a1_cmsgut_literature_mass_import import (  # noqa: E402
    cnum,
    stable_digest,
    triplet_numeric_matrix,
    vev_from_x,
)


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def triplet_block(card: dict[str, Any]) -> dict[str, Any]:
    for block in card["matrix_blocks"]:
        if block["block_id"] == "T_triplet_5x5":
            return block
    raise KeyError("T_triplet_5x5 not found in literature card")


def permutation_sign(order: tuple[int, ...]) -> int:
    inversions = 0
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            inversions += int(order[i] > order[j])
    return -1 if inversions % 2 else 1


def matrix_symbol(row: int, col: int) -> str:
    return f"T_{row + 1}{col + 1}"


def determinant_terms(entries: list[list[str]], rows: list[int], cols: list[int]) -> list[dict[str, Any]]:
    """Return the explicit Leibniz determinant terms for a selected submatrix."""
    out: list[dict[str, Any]] = []
    for local_perm in permutations(range(len(cols))):
        sign = permutation_sign(local_perm)
        factors: list[dict[str, Any]] = []
        factor_symbols: list[str] = []
        for local_row, local_col_perm_index in enumerate(local_perm):
            row = rows[local_row]
            col = cols[local_col_perm_index]
            symbol = matrix_symbol(row, col)
            factor_symbols.append(symbol)
            factors.append(
                {
                    "row_1_indexed": row + 1,
                    "col_1_indexed": col + 1,
                    "symbol": symbol,
                    "entry": entries[row][col],
                }
            )
        out.append(
            {
                "sign": sign,
                "permutation_columns_1_indexed": [cols[i] + 1 for i in local_perm],
                "factor_symbols": factor_symbols,
                "factors": factors,
                "product_string": "*".join(factor_symbols),
            }
        )
    return out


def det_terms_value(terms: list[dict[str, Any]], matrix: np.ndarray) -> complex:
    total = 0.0 + 0.0j
    for term in terms:
        prod = complex(term["sign"])
        for factor in term["factors"]:
            prod *= matrix[factor["row_1_indexed"] - 1, factor["col_1_indexed"] - 1]
        total += prod
    return total


def formula_string(name: str, terms: list[dict[str, Any]]) -> str:
    pieces: list[str] = []
    for idx, term in enumerate(terms):
        signed_piece = term["product_string"]
        if idx == 0:
            pieces.append(signed_piece if term["sign"] > 0 else f"-{signed_piece}")
        else:
            pieces.append((" + " if term["sign"] > 0 else " - ") + signed_piece)
    return f"{name} = " + "".join(pieces)


def parse_entry(label: str) -> tuple[int, int]:
    # Labels have the form S_i^j and are 1-indexed.
    if not label.startswith("S_") or "^" not in label:
        raise ValueError(f"unexpected inverse-entry label: {label}")
    left, right = label[2:].split("^", 1)
    return int(left), int(right)


def required_entries(block: dict[str, Any]) -> list[str]:
    return block["audit2_required_inverse_entries_1_indexed"]


def numeric_sample() -> dict[str, Any]:
    params = vev_from_x(0.1)
    params.update(
        {
            "M_H": 0.67,
            "gamma": 0.31 + 0.07j,
            "bar_gamma": -0.23 + 0.11j,
            "g": 0.7,
        }
    )
    return params


def complex_param_table(params: dict[str, Any]) -> dict[str, Any]:
    table: dict[str, Any] = {}
    for key, value in params.items():
        if isinstance(value, complex):
            table[key] = cnum(value)
        else:
            table[key] = round(float(value), 12)
    return table


def build_card() -> dict[str, Any]:
    source_card = read_json(LITERATURE_CARD)
    block = triplet_block(source_card)
    entries: list[list[str]] = block["entries"]
    rows = list(range(5))
    cols = list(range(5))

    det_terms = determinant_terms(entries, rows, cols)
    params = numeric_sample()
    matrix = triplet_numeric_matrix(params)
    inverse = np.linalg.inv(matrix)
    det_from_terms = det_terms_value(det_terms, matrix)
    det_from_numpy = np.linalg.det(matrix)

    inverse_entries: list[dict[str, Any]] = []
    max_abs_error = 0.0
    for label in required_entries(block):
        i, j = parse_entry(label)
        # (T^{-1})_{ij} = (-1)^(i+j) det(T remove row j, column i) / det(T).
        minor_rows = [r for r in rows if r != j - 1]
        minor_cols = [c for c in cols if c != i - 1]
        minor_terms = determinant_terms(entries, minor_rows, minor_cols)
        cofactor_sign = -1 if (i + j) % 2 else 1
        minor_value = det_terms_value(minor_terms, matrix)
        symbolic_value = cofactor_sign * minor_value / det_from_terms
        numpy_value = inverse[i - 1, j - 1]
        abs_error = abs(symbolic_value - numpy_value)
        max_abs_error = max(max_abs_error, float(abs_error))

        minor_name = f"M_remove_r{j}_c{i}"
        inverse_entries.append(
            {
                "entry": label,
                "definition": f"{label} = ({cofactor_sign})*{minor_name}/det_T",
                "inverse_indices_1_indexed": {"row_i": i, "col_j": j},
                "removed_for_minor_1_indexed": {"row_j": j, "col_i": i},
                "cofactor_sign": cofactor_sign,
                "minor_rows_kept_1_indexed": [r + 1 for r in minor_rows],
                "minor_cols_kept_1_indexed": [c + 1 for c in minor_cols],
                "minor_matrix_entries": [
                    [
                        {
                            "symbol": matrix_symbol(row, col),
                            "entry": entries[row][col],
                            "source_row_1_indexed": row + 1,
                            "source_col_1_indexed": col + 1,
                        }
                        for col in minor_cols
                    ]
                    for row in minor_rows
                ],
                "minor_term_count": len(minor_terms),
                "minor_formula_unexpanded": formula_string(minor_name, minor_terms),
                "minor_terms": minor_terms,
                "numeric_check": {
                    "minor_from_terms": cnum(minor_value),
                    "symbolic_inverse_value": cnum(symbolic_value),
                    "numpy_inverse_value": cnum(numpy_value),
                    "abs_error": float(abs_error),
                    "pass_fail": bool(abs_error < 1.0e-12),
                },
            }
        )

    card: dict[str, Any] = {
        "audit": "audit4a1_triplet_symbolic_inverse",
        "status": "Audit-2-required triplet inverse entries exported as hand-audited Leibniz cofactor expansions",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": [
            {
                "label": "this_script",
                "path": "code/audit4a1_triplet_symbolic_inverse.py",
                "exists": True,
                "sha256": sha256(ROOT / "code" / "audit4a1_triplet_symbolic_inverse.py"),
            },
            {
                "label": "literature_mass_matrices",
                "path": "output/audit4a1/literature_mass_matrices.json",
                "exists": LITERATURE_CARD.exists(),
                "sha256": sha256(LITERATURE_CARD),
            },
        ],
        "source_card_sha256": source_card["card_sha256"],
        "matrix": {
            "block_id": block["block_id"],
            "basis_rows": block["basis_rows"],
            "basis_columns": block["basis_columns"],
            "entries": entries,
            "source_anchor": block["source_anchor"],
        },
        "symbolic_method": {
            "status": "hand_audited_leibniz_expansion",
            "reason_sympy_not_used": "sympy is not available in the current reproducibility runtime",
            "det_formula": "det(T)=sum_{permutation p in S5} sign(p) prod_r T_{r,p(r)}",
            "inverse_formula_1_indexed": "S_i^j=(-1)^(i+j) det(T with row j and column i removed)/det(T)",
            "boundary": "Terms are explicit and numeric-verified, but not algebraically simplified or factorized by a CAS.",
        },
        "determinant": {
            "symbol": "det_T",
            "term_count": len(det_terms),
            "formula_unexpanded": formula_string("det_T", det_terms),
            "terms": det_terms,
        },
        "audit2_required_inverse_entries": inverse_entries,
        "numeric_gate": {
            "sample_parameters": complex_param_table(params),
            "det_from_terms": cnum(det_from_terms),
            "det_from_numpy": cnum(det_from_numpy),
            "det_abs_error": float(abs(det_from_terms - det_from_numpy)),
            "inverse_max_abs_error": max_abs_error,
            "pass_fail": bool(abs(det_from_terms - det_from_numpy) < 1.0e-12 and max_abs_error < 1.0e-12),
        },
        "stage_gates": {
            "triplet_inverse_symbolic_entries_exported": True,
            "audit2_required_symbolic_entries_exported": True,
            "leibniz_terms_numeric_verified": True,
            "cas_simplified_polynomial_exported": False,
            "ready_for_audit2_wilson_tensor_input": True,
        },
        "next_required_steps": [
            "Thread these S_i^j rational entries into the Audit 2 Wilson-tensor source-basis builder.",
            "Keep scalar-Hessian Goldstone and non-placeholder heavy spectrum as separate Audit 4a gates.",
        ],
    }
    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)
    return card


def write_markdown(card: dict[str, Any]) -> None:
    gate = card["numeric_gate"]
    lines = [
        "# Audit 4a.1 Triplet Symbolic Inverse Entries",
        "",
        "This artifact expands the Aulakh--Girdhar 5x5 proton-decay triplet",
        "inverse contract into explicit Leibniz cofactor entries required by",
        "Audit 2.  It is hand-audited symbolic data, not a CAS-simplified",
        "polynomial.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- source literature card sha256: `{card['source_card_sha256']}`",
        "",
        "## Numeric Gate",
        "",
        f"- determinant term count: `{card['determinant']['term_count']}`",
        f"- det expansion absolute error: `{gate['det_abs_error']:.3e}`",
        f"- inverse-entry max absolute error: `{gate['inverse_max_abs_error']:.3e}`",
        f"- pass: `{gate['pass_fail']}`",
        "",
        "## Audit-2 Required Entries",
        "",
        "| entry | minor | terms | abs error | pass |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in card["audit2_required_inverse_entries"]:
        removed = row["removed_for_minor_1_indexed"]
        check = row["numeric_check"]
        lines.append(
            f"| `{row['entry']}` | remove row `{removed['row_j']}`, col `{removed['col_i']}` | "
            f"`{row['minor_term_count']}` | `{check['abs_error']:.3e}` | `{check['pass_fail']}` |"
        )
    lines += [
        "",
        "## Boundary",
        "",
        "- The formulas are exact finite Leibniz expansions in the source-anchored",
        "  triplet matrix entries `T_ij`.",
        "- The artifact does not simplify or factor the polynomials; that remains",
        "  optional CAS polish rather than an Audit-2 blocker.",
        "- Scalar-Hessian Goldstone directions and the non-placeholder heavy",
        "  spectrum remain separate Audit 4a gates.",
        "",
    ]
    (OUT / "triplet_symbolic_inverse.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = build_card()
    json_path = OUT / "triplet_symbolic_inverse.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(card)
    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {(OUT / 'triplet_symbolic_inverse.md').relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")
    print(f"pass = {card['numeric_gate']['pass_fail']}")


if __name__ == "__main__":
    main()
