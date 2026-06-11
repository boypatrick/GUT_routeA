#!/usr/bin/env python3
"""Verify the small algebra/index claims used in Route-D D1-F.

This is intentionally modest.  It checks only the local representation
dimension identity and the P1 line-bundle cohomology arithmetic.  It does not
construct a global F-theory compactification or verify hypercharge flux
masslessness.
"""

from __future__ import annotations

import json
from pathlib import Path


def p1_h0_h1(degree: int) -> tuple[int, int]:
    """Return h0, h1 for O(degree) on P1."""
    h0 = max(degree + 1, 0)
    h1 = max(-degree - 1, 0)
    return h0, h1


def main() -> None:
    e6_branch = {
        "45_0": 45,
        "1_0": 1,
        "16_-3": 16,
        "bar16_+3": 16,
    }
    e6_total = sum(e6_branch.values())

    canonical_half_degree = -1  # K_{P1}^{1/2} = O(-1)
    flux_degree = 3
    carrier_degree = canonical_half_degree + flux_degree
    h0, h1 = p1_h0_h1(carrier_degree)

    report = {
        "e6_to_so10_u1_branching": e6_branch,
        "e6_adjoint_dimension_sum": e6_total,
        "e6_dimension_check_pass": e6_total == 78,
        "p1_spin_half_degree": canonical_half_degree,
        "flux_degree": flux_degree,
        "effective_carrier": f"O({carrier_degree})",
        "h0_P1_effective_carrier": h0,
        "h1_P1_effective_carrier": h1,
        "three_family_index_check_pass": (h0, h1) == (3, 0),
        "global_f_theory_claim": "not_checked",
        "hypercharge_flux_masslessness": "conditional_global_topology_required",
    }

    out_dir = Path(__file__).resolve().parents[1] / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "d1_f_theory_local_placement.json"
    md_path = out_dir / "d1_f_theory_local_placement.md"

    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    md_path.write_text(
        "# D1-F F-theory Local Placement Verification\n\n"
        f"E6 branching dimension sum: `{e6_total}`.\n"
        f"Branching dimension check passes: `{report['e6_dimension_check_pass']}`.\n"
        f"Effective P1 carrier: `{report['effective_carrier']}` from "
        f"`O({canonical_half_degree}) tensor O({flux_degree})`.\n"
        f"`h0 = {h0}`, `h1 = {h1}`.\n"
        f"Three-family index check passes: `{report['three_family_index_check_pass']}`.\n\n"
        "Boundary: this verification does not construct a global F-theory "
        "compactification and does not check the global topological condition "
        "for massless hypercharge flux.\n"
    )

    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    print(
        "D1-F checks: "
        f"E6 dimension={e6_total}, carrier=O({carrier_degree}), h0={h0}, h1={h1}"
    )


if __name__ == "__main__":
    main()
