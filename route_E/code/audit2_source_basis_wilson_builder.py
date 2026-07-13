#!/usr/bin/env python3
"""Build the Audit 2 source-basis Wilson-tensor contract.

This is the first Audit 2 artifact that consumes the CMSGUT triplet inverse
entries exported by Audit 4a.1.  It deliberately stops at the source basis:
no physical flavor rotations, no SUSY dressing, and no proton-channel widths
are computed here.

The source-basis formulas are

    C5L_{ab,cd} = sum_{i,j} S_i^j (Y_QQ^i)_{ab} (Y_QL^j)_{cd},
    C5R_{ab,cd} = sum_{i,j} S_i^j (Y_UE^i)_{ab} (Y_UD^j)_{cd},

where S=T^{-1} is the Aulakh--Girdhar triplet inverse.  The index orientation
is inherited from the symbolic inverse card:

    S_i^j = (T^{-1})_{i,j}
          = (-1)^(i+j) det(T with row j and column i removed)/det(T).

The output is therefore a formal Wilson-builder contract plus a numerical
random-tensor gate proving that the hand-audited symbolic S entries produce the
same source-basis tensors as the direct inverse entries on the generic F-flat
sample.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit2"
SYMBOLIC_INVERSE = ROOT / "output" / "audit4a1" / "triplet_symbolic_inverse.json"
LITERATURE_MASSES = ROOT / "output" / "audit4a1" / "literature_mass_matrices.json"


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_digest(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def cnum(z: complex, ndigits: int = 12) -> dict[str, float]:
    return {"re": round(float(np.real(z)), ndigits), "im": round(float(np.imag(z)), ndigits)}


def complex_from_cell(cell: dict[str, float]) -> complex:
    return complex(cell["re"], cell["im"])


def find_triplet_block(card: dict[str, Any]) -> dict[str, Any]:
    for block in card["matrix_blocks"]:
        if block["block_id"] == "T_triplet_5x5":
            return block
    raise KeyError("T_triplet_5x5 not found")


def parse_entry(label: str) -> tuple[int, int]:
    left, right = label[2:].split("^", 1)
    return int(left), int(right)


def source_label(block: dict[str, Any], index_1: int, side: str) -> str:
    if side == "row":
        return block["basis_rows"][index_1 - 1]
    if side == "column":
        return block["basis_columns"][index_1 - 1]
    raise ValueError(side)


def random_symmetric_matrix(rng: np.random.Generator) -> np.ndarray:
    raw = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    return 0.5 * (raw + raw.T)


def random_general_matrix(rng: np.random.Generator) -> np.ndarray:
    return rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))


def source_basis_tensor(
    entries: list[dict[str, Any]],
    s_values: dict[str, complex],
    left_tensors: dict[int, np.ndarray],
    right_tensors: dict[int, np.ndarray],
) -> np.ndarray:
    out = np.zeros((3, 3, 3, 3), dtype=complex)
    for entry in entries:
        i = entry["inverse_indices_1_indexed"]["row_i"]
        j = entry["inverse_indices_1_indexed"]["col_j"]
        label = entry["entry"]
        out += s_values[label] * np.einsum("ab,cd->abcd", left_tensors[i], right_tensors[j], optimize=True)
    return out


def build_card() -> dict[str, Any]:
    symbolic = read_json(SYMBOLIC_INVERSE)
    literature = read_json(LITERATURE_MASSES)
    triplet = find_triplet_block(literature)
    entries = symbolic["audit2_required_inverse_entries"]

    source_entries: list[dict[str, Any]] = []
    for row in entries:
        label = row["entry"]
        i, j = parse_entry(label)
        source_entries.append(
            {
                "entry": label,
                "inverse_indices_1_indexed": {"row_i": i, "col_j": j},
                "left_triplet_source_column": source_label(triplet, i, "column"),
                "right_triplet_source_row": source_label(triplet, j, "row"),
                "cofactor_definition": row["definition"],
                "minor_formula_unexpanded": row["minor_formula_unexpanded"],
                "symbolic_numeric_value": row["numeric_check"]["symbolic_inverse_value"],
                "direct_inverse_numeric_value": row["numeric_check"]["numpy_inverse_value"],
                "entry_abs_error": row["numeric_check"]["abs_error"],
            }
        )

    symbolic_s = {
        row["entry"]: complex_from_cell(row["numeric_check"]["symbolic_inverse_value"])
        for row in entries
    }
    direct_s = {
        row["entry"]: complex_from_cell(row["numeric_check"]["numpy_inverse_value"])
        for row in entries
    }

    rng = np.random.default_rng(20260616)
    left_indices = sorted({row["inverse_indices_1_indexed"]["row_i"] for row in entries})
    right_indices = sorted({row["inverse_indices_1_indexed"]["col_j"] for row in entries})

    y_qq = {idx: random_symmetric_matrix(rng) for idx in left_indices}
    y_ue = {idx: random_symmetric_matrix(rng) for idx in left_indices}
    y_ql = {idx: random_general_matrix(rng) for idx in right_indices}
    y_ud = {idx: random_general_matrix(rng) for idx in right_indices}

    c5l_symbolic = source_basis_tensor(entries, symbolic_s, y_qq, y_ql)
    c5l_direct = source_basis_tensor(entries, direct_s, y_qq, y_ql)
    c5r_symbolic = source_basis_tensor(entries, symbolic_s, y_ue, y_ud)
    c5r_direct = source_basis_tensor(entries, direct_s, y_ue, y_ud)

    c5l_error = float(np.max(np.abs(c5l_symbolic - c5l_direct)))
    c5r_error = float(np.max(np.abs(c5r_symbolic - c5r_direct)))

    card: dict[str, Any] = {
        "audit": "audit2_source_basis_wilson_builder",
        "status": "source-basis C5L/C5R Wilson-tensor contract built from Audit 4a.1 symbolic triplet inverse entries",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "source_manifest": [
            {
                "label": "this_script",
                "path": "code/audit2_source_basis_wilson_builder.py",
                "exists": True,
                "sha256": sha256(ROOT / "code" / "audit2_source_basis_wilson_builder.py"),
            },
            {
                "label": "triplet_symbolic_inverse",
                "path": "output/audit4a1/triplet_symbolic_inverse.json",
                "exists": SYMBOLIC_INVERSE.exists(),
                "sha256": sha256(SYMBOLIC_INVERSE),
            },
            {
                "label": "literature_mass_matrices",
                "path": "output/audit4a1/literature_mass_matrices.json",
                "exists": LITERATURE_MASSES.exists(),
                "sha256": sha256(LITERATURE_MASSES),
            },
        ],
        "input_digests": {
            "triplet_symbolic_inverse_card_sha256": symbolic["card_sha256"],
            "literature_mass_matrix_card_sha256": literature["card_sha256"],
        },
        "triplet_index_convention": {
            "matrix": "T_triplet_5x5",
            "row_basis_barred_sources": triplet["basis_rows"],
            "column_basis_unbarred_sources": triplet["basis_columns"],
            "inverse_orientation": "S_i^j=(T^{-1})_{i,j}=(-1)^(i+j) det(T remove row j, column i)/det(T)",
            "left_source_indices_used": left_indices,
            "right_source_indices_used": right_indices,
            "orientation_boundary": (
                "This card is a source-basis contract.  The final assignment of "
                "source labels to physical QQ, QL, UE, and UD Yukawa tensors is "
                "deferred to Audit 1 physical flavor rotations and the triplet "
                "Yukawa map."
            ),
        },
        "source_entries": source_entries,
        "formal_wilson_contract": {
            "C5L": "C5L_{ab,cd}=sum_{(i,j)} S_i^j (Y_QQ^i)_{ab} (Y_QL^j)_{cd}",
            "C5R": "C5R_{ab,cd}=sum_{(i,j)} S_i^j (Y_UE^i)_{ab} (Y_UD^j)_{cd}",
            "family_index_range": "a,b,c,d=1..3",
            "summed_source_entries": [row["entry"] for row in source_entries],
            "LLLL_tensor_symmetry_assumption": "Y_QQ^i is treated as symmetric in the random numerical gate only; the formal contract keeps the tensor label explicit.",
            "RRRR_tensor_symmetry_assumption": "Y_UE^i is treated as symmetric in the random numerical gate only; the formal contract keeps the tensor label explicit.",
        },
        "numeric_gate": {
            "random_seed": 20260616,
            "random_tensor_role": "checks symbolic S entries against direct inverse entries in source-basis C5 contractions; not a phenomenological flavor point",
            "C5L_symbolic_vs_direct_max_abs_error": c5l_error,
            "C5R_symbolic_vs_direct_max_abs_error": c5r_error,
            "S_entry_max_abs_error": float(max(row["entry_abs_error"] for row in source_entries)),
            "pass_fail": bool(c5l_error < 1.0e-12 and c5r_error < 1.0e-12),
        },
        "stage_gates": {
            "triplet_symbolic_inverse_consumed": True,
            "source_basis_C5L_contract_exported": True,
            "source_basis_C5R_contract_exported": True,
            "random_tensor_numeric_gate_passed": bool(c5l_error < 1.0e-12 and c5r_error < 1.0e-12),
            "physical_flavor_rotations_attached": False,
            "channel_dressing_or_widths_computed": False,
            "publication_level_d5_proton_decay": False,
        },
        "next_required_steps": [
            "Attach Audit 1 physical flavor rotations and triplet Yukawa maps.",
            "Promote the source-basis C5L/C5R tensors to mass-basis Wilson tensors.",
            "Only after the mass-basis tensors exist, add SUSY dressing and channel widths.",
        ],
    }

    digest_payload = {k: v for k, v in card.items() if k not in {"card_sha256", "created_utc"}}
    card["card_sha256"] = stable_digest(digest_payload)
    return card


def write_markdown(card: dict[str, Any]) -> None:
    gate = card["numeric_gate"]
    lines = [
        "# Audit 2 Source-Basis Wilson Contract",
        "",
        "This artifact consumes the Audit 4a.1 symbolic triplet inverse entries",
        "and exports formal source-basis `C5L` and `C5R` Wilson-tensor contracts.",
        "It is not yet a mass-basis proton-decay calculation.",
        "",
        "## Digest",
        "",
        f"- card sha256: `{card['card_sha256']}`",
        f"- symbolic inverse card: `{card['input_digests']['triplet_symbolic_inverse_card_sha256']}`",
        "",
        "## Formal Contract",
        "",
        f"- `C5L`: `{card['formal_wilson_contract']['C5L']}`",
        f"- `C5R`: `{card['formal_wilson_contract']['C5R']}`",
        "",
        "## Source Entries",
        "",
        "| entry | left column source | right row source | source error |",
        "| --- | --- | --- | --- |",
    ]
    for row in card["source_entries"]:
        lines.append(
            f"| `{row['entry']}` | `{row['left_triplet_source_column']}` | "
            f"`{row['right_triplet_source_row']}` | `{row['entry_abs_error']:.3e}` |"
        )
    lines += [
        "",
        "## Numeric Gate",
        "",
        f"- random seed: `{gate['random_seed']}`",
        f"- `C5L` symbolic/direct max error: `{gate['C5L_symbolic_vs_direct_max_abs_error']:.3e}`",
        f"- `C5R` symbolic/direct max error: `{gate['C5R_symbolic_vs_direct_max_abs_error']:.3e}`",
        f"- pass: `{gate['pass_fail']}`",
        "",
        "## Boundary",
        "",
        "- This is a source-basis Wilson contract only.",
        "- Physical flavor rotations, triplet Yukawa maps, SUSY dressing, lattice",
        "  matrix elements, and channel widths remain deferred.",
        "",
    ]
    (OUT / "source_basis_wilson_contract.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    card = build_card()
    json_path = OUT / "source_basis_wilson_contract.json"
    json_path.write_text(json.dumps(card, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(card)
    print(f"wrote {json_path.relative_to(ROOT)}")
    print(f"wrote {(OUT / 'source_basis_wilson_contract.md').relative_to(ROOT)}")
    print(f"card_sha256 = {card['card_sha256']}")
    print(f"pass = {card['numeric_gate']['pass_fail']}")


if __name__ == "__main__":
    main()
