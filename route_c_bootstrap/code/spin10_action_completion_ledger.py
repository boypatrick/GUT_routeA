#!/usr/bin/env python3
"""Route-C P9 Spin(10)-locked action-completion ledger.

P9 stops broadening the comparison set and locks the Route-C branch to the
P8-leading finite field-theory candidate, Spin(10).  It supplies the first
machine-checkable action-level data requested by P4--P7:

* a 16-state Pati--Salam/SM-face basis;
* sparse broken-generator transition maps for the Spin(10) half-spinor;
* Goldstone/Higgs completion formulae for broken-vector Ward identities;
* symbolic mediator masses/couplings;
* physical-flavor rotation and proton-bound input interfaces.

This is still a completion ledger, not a physical proton-lifetime calculation.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
P4_JSON = OUT / "ward_identity_ledger.json"
P6_JSON = OUT / "high_energy_growth_audit.json"
P7_JSON = OUT / "low_energy_matching_proton_report.json"
P8_JSON = OUT / "theory_comparison_scorecard.json"
OUT_JSON = OUT / "spin10_action_completion_ledger.json"
OUT_MD = OUT / "spin10_action_completion_ledger.md"


COLORS = ["r", "g", "b"]
SU4 = ["r", "g", "b", "ell"]
SU2L = ["u", "d"]
SU2R = ["minus", "plus"]


def load_json(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"missing {path}; run earlier Route-C stages first")
    return json.loads(path.read_text(encoding="utf-8"))


def frac(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fjson(value: Fraction) -> dict[str, int | float | str]:
    return {
        "fraction": frac(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
        "float": float(value),
    }


def basis_states() -> list[dict]:
    states: list[dict] = []
    index = 0

    # (4,2,1): Q + L.
    for A in SU4:
        for alpha in SU2L:
            if A in COLORS:
                multiplet = "Q"
                label = f"Q_{A}_{alpha}"
                y = Fraction(1, 6)
                bl = Fraction(1, 3)
                su3 = "3"
                color = A
            else:
                multiplet = "L"
                label = f"L_{alpha}"
                y = Fraction(-1, 2)
                bl = Fraction(-1, 1)
                su3 = "1"
                color = "1"
            t3l = Fraction(1, 2) if alpha == "u" else Fraction(-1, 2)
            states.append(
                {
                    "index": index,
                    "label": label,
                    "multiplet": multiplet,
                    "ps_irrep": "(4,2,1)",
                    "su4_index": A,
                    "su2l_index": alpha,
                    "su2r_index": "singlet",
                    "SU3_C": su3,
                    "color": color,
                    "Y": fjson(y),
                    "B_minus_L": fjson(bl),
                    "T3L": fjson(t3l),
                    "T3R": fjson(Fraction(0)),
                }
            )
            index += 1

    # (bar4,1,2): u^c + d^c + nu^c + e^c.
    for A in SU4:
        for rho in SU2R:
            if A in COLORS and rho == "minus":
                multiplet = "u^c"
                label = f"u^c_{A}"
                y = Fraction(-2, 3)
                bl = Fraction(-1, 3)
                su3 = "bar3"
                color = A
            elif A in COLORS and rho == "plus":
                multiplet = "d^c"
                label = f"d^c_{A}"
                y = Fraction(1, 3)
                bl = Fraction(-1, 3)
                su3 = "bar3"
                color = A
            elif A == "ell" and rho == "minus":
                multiplet = "nu^c"
                label = "nu^c"
                y = Fraction(0)
                bl = Fraction(1)
                su3 = "1"
                color = "1"
            else:
                multiplet = "e^c"
                label = "e^c"
                y = Fraction(1)
                bl = Fraction(1)
                su3 = "1"
                color = "1"
            t3r = Fraction(-1, 2) if rho == "minus" else Fraction(1, 2)
            states.append(
                {
                    "index": index,
                    "label": label,
                    "multiplet": multiplet,
                    "ps_irrep": "(bar4,1,2)",
                    "su4_index": A,
                    "su2l_index": "singlet",
                    "su2r_index": rho,
                    "SU3_C": su3,
                    "color": color,
                    "Y": fjson(y),
                    "B_minus_L": fjson(bl),
                    "T3L": fjson(Fraction(0)),
                    "T3R": fjson(t3r),
                }
            )
            index += 1

    return states


def idx_by_quantum_numbers(states: list[dict]) -> dict[tuple[str, str, str], int]:
    out = {}
    for st in states:
        out[(st["ps_irrep"], st["su4_index"], st["su2l_index"])] = st["index"]
        out[(st["ps_irrep"], st["su4_index"], st["su2r_index"])] = st["index"]
    return out


def state_by_index(states: list[dict]) -> dict[int, dict]:
    return {st["index"]: st for st in states}


def transition(
    states_by_i: dict[int, dict], source: int, target: int, coefficient: str = "1"
) -> dict:
    src = states_by_i[source]
    tgt = states_by_i[target]
    dy = Fraction(tgt["Y"]["numerator"], tgt["Y"]["denominator"]) - Fraction(
        src["Y"]["numerator"], src["Y"]["denominator"]
    )
    dbl = Fraction(
        tgt["B_minus_L"]["numerator"], tgt["B_minus_L"]["denominator"]
    ) - Fraction(src["B_minus_L"]["numerator"], src["B_minus_L"]["denominator"])
    return {
        "source_index": source,
        "source_label": src["label"],
        "source_multiplet": src["multiplet"],
        "target_index": target,
        "target_label": tgt["label"],
        "target_multiplet": tgt["multiplet"],
        "coefficient": coefficient,
        "delta_Y": fjson(dy),
        "delta_B_minus_L": fjson(dbl),
    }


def sparse_map(label: str, family: str, entries: list[dict], note: str) -> dict:
    return {
        "label": label,
        "generator_family": family,
        "sparse_transition_entries": entries,
        "entry_count": len(entries),
        "normalization": "sparse Clebsch map; hermitian vector generators are obtained by adding the conjugate map and its i-rotated partner when appropriate",
        "note": note,
    }


def leptoquark_maps(states: list[dict]) -> list[dict]:
    idx = idx_by_quantum_numbers(states)
    by_i = state_by_index(states)
    maps = []
    for color in COLORS:
        # E_{ell,color}: Q -> L and nu/e -> u/d in the conjugate block.
        entries_down = []
        for alpha in SU2L:
            entries_down.append(
                transition(
                    by_i,
                    idx[("(4,2,1)", color, alpha)],
                    idx[("(4,2,1)", "ell", alpha)],
                )
            )
        for rho in SU2R:
            entries_down.append(
                transition(
                    by_i,
                    idx[("(bar4,1,2)", "ell", rho)],
                    idx[("(bar4,1,2)", color, rho)],
                    "-1",
                )
            )
        maps.append(
            sparse_map(
                f"SU4/E_SM leptoquark E_ell{color}",
                "Pati-Salam SU(4)_C leptoquark",
                entries_down,
                "One complex broken SU(4) root; its conjugate is E_color,ell.",
            )
        )

        entries_up = []
        for alpha in SU2L:
            entries_up.append(
                transition(
                    by_i,
                    idx[("(4,2,1)", "ell", alpha)],
                    idx[("(4,2,1)", color, alpha)],
                )
            )
        for rho in SU2R:
            entries_up.append(
                transition(
                    by_i,
                    idx[("(bar4,1,2)", color, rho)],
                    idx[("(bar4,1,2)", "ell", rho)],
                    "-1",
                )
            )
        maps.append(
            sparse_map(
                f"SU4/E_SM leptoquark E_{color}ell",
                "Pati-Salam SU(4)_C leptoquark",
                entries_up,
                "Conjugate complex broken SU(4) root.",
            )
        )
    return maps


def su2r_maps(states: list[dict]) -> list[dict]:
    idx = idx_by_quantum_numbers(states)
    by_i = state_by_index(states)
    entries_plus = []
    entries_minus = []
    for A in SU4:
        entries_plus.append(
            transition(
                by_i,
                idx[("(bar4,1,2)", A, "minus")],
                idx[("(bar4,1,2)", A, "plus")],
            )
        )
        entries_minus.append(
            transition(
                by_i,
                idx[("(bar4,1,2)", A, "plus")],
                idx[("(bar4,1,2)", A, "minus")],
            )
        )
    return [
        sparse_map(
            "SU2_R T_R^+",
            "broken SU(2)_R charged current",
            entries_plus,
            "Raises T3R and carries SM hypercharge +1 after SU(2)_R breaking.",
        ),
        sparse_map(
            "SU2_R T_R^-",
            "broken SU(2)_R charged current",
            entries_minus,
            "Lowers T3R and carries SM hypercharge -1 after SU(2)_R breaking.",
        ),
    ]


def off_face_maps(states: list[dict]) -> list[dict]:
    idx = idx_by_quantum_numbers(states)
    by_i = state_by_index(states)
    maps = []
    for i, A in enumerate(SU4):
        for B in SU4[i + 1 :]:
            for alpha in SU2L:
                for rho in SU2R:
                    entries = [
                        transition(
                            by_i,
                            idx[("(4,2,1)", A, alpha)],
                            idx[("(bar4,1,2)", B, rho)],
                            "1/sqrt(2)",
                        ),
                        transition(
                            by_i,
                            idx[("(4,2,1)", B, alpha)],
                            idx[("(bar4,1,2)", A, rho)],
                            "-1/sqrt(2)",
                        ),
                    ]
                    maps.append(
                        sparse_map(
                            f"(6,2,2) C_[{A}{B}],{alpha},{rho}",
                            "Spin(10)/Pati-Salam off-face generator",
                            entries,
                            (
                                "Antisymmetric SU(4) Clebsch map from (4,2,1) "
                                "to (bar4,1,2).  The conjugate map gives the "
                                "reverse chirality-block transition."
                            ),
                        )
                    )
    return maps


def multiplet_transition_set(maps: Iterable[dict], include_reverse: bool = True) -> set[tuple[str, str]]:
    transitions = set()
    for mp in maps:
        for entry in mp["sparse_transition_entries"]:
            pair = (entry["source_multiplet"], entry["target_multiplet"])
            transitions.add(pair)
            if include_reverse:
                transitions.add((pair[1], pair[0]))
    return transitions


def p4_broken_transition_set(p4: dict) -> set[tuple[str, str]]:
    pairs = set()
    for row in p4["transition_checks"]:
        if row["requires_broken_sector_completion"]:
            pairs.add((row["source"], row["target"]))
    return pairs


def generator_checks(p4: dict, maps: list[dict]) -> dict:
    p4_pairs = p4_broken_transition_set(p4)
    generated = multiplet_transition_set(maps)
    missing = sorted(p4_pairs - generated)
    extra = sorted(generated - p4_pairs)
    return {
        "p4_broken_transition_pairs": len(p4_pairs),
        "transition_pairs_generated_by_spin10_maps": len(generated),
        "missing_p4_pairs": [{"source": a, "target": b} for a, b in missing],
        "missing_p4_pair_interpretation": (
            "These P4 bookkeeping transitions are not generated by the Spin(10) "
            "adjoint gauge-vector maps in the 16.  They remain possible "
            "scalar/source/auxiliary completion channels if a later branch "
            "activates them, but P9 must not misidentify them as D5 gauge "
            "generators."
        ),
        "extra_generated_pairs": [{"source": a, "target": b} for a, b in extra],
        "covers_all_p4_broken_pairs": not missing,
        "generated_pairs_are_subset_of_p4_broken_pairs": not extra,
        "spin10_adjoint_vector_pair_count": len(generated),
    }


def completion_templates() -> list[dict]:
    return [
        {
            "sector": "Spin(10) -> Pati-Salam",
            "broken_generators": "(6,2,2)",
            "candidate_order_parameters": ["45_H", "210_H", "adjoint/source spurion"],
            "mass_formula": "M_X^2 = g_10^2 <v| T_X T_X |v> after diagonalizing the broken-generator mass matrix",
            "goldstone_formula": "phi_X proportional to T_X v",
            "ward_target": "p_mu M^mu(V_X)=M_X M(phi_X)",
            "status": "template_only_order_parameter_not_chosen",
        },
        {
            "sector": "SU(4)_C -> SU(3)_C x U(1)_{B-L}",
            "broken_generators": "SU(4) leptoquark roots E_{ell a}, E_{a ell}",
            "candidate_order_parameters": ["45_H/210_H direction", "Pati-Salam breaking source"],
            "mass_formula": "M_LQ^2 = g_4^2 <v_4| T_LQ T_LQ |v_4>",
            "goldstone_formula": "phi_LQ proportional to T_LQ v_4",
            "ward_target": "p_mu M^mu(V_LQ)=M_LQ M(phi_LQ)",
            "status": "template_only_order_parameter_not_chosen",
        },
        {
            "sector": "SU(2)_R -> U(1)_{T3R}",
            "broken_generators": "T_R^+, T_R^-",
            "candidate_order_parameters": ["126_H or 16_H B-L breaking branch", "right-handed doublet/triplet source"],
            "mass_formula": "M_WR^2 = g_R^2 <v_R| T_R^- T_R^+ |v_R>",
            "goldstone_formula": "phi_R^+ proportional to T_R^+ v_R and conjugate",
            "ward_target": "p_mu M^mu(W_R)=M_WR M(phi_R)",
            "status": "template_only_order_parameter_not_chosen",
        },
    ]


def mediator_inputs() -> list[dict]:
    return [
        {
            "symbol": "g_10",
            "meaning": "Spin(10) unified gauge coupling at matching scale",
            "required_for": ["all gauge-vector Wilson coefficients", "threshold matching"],
            "status": "symbolic_required",
        },
        {
            "symbol": "M_{(6,2,2)}",
            "meaning": "mass eigenvalues for Spin(10)/Pati-Salam broken vectors",
            "required_for": ["P6 high-energy completion", "P7 dimension-six matching"],
            "status": "symbolic_required",
        },
        {
            "symbol": "M_{LQ}",
            "meaning": "Pati-Salam SU(4) leptoquark vector masses",
            "required_for": ["proton/leptoquark matching", "threshold matching"],
            "status": "symbolic_required",
        },
        {
            "symbol": "M_{W_R}",
            "meaning": "broken SU(2)_R vector mass",
            "required_for": ["right-current matching", "Ward completion"],
            "status": "symbolic_required",
        },
        {
            "symbol": "M_{T}, Y_T",
            "meaning": "colored Higgs/scalar triplet masses and Yukawa couplings if scalar branches are activated",
            "required_for": ["dimension-five proton audit", "scalar/auxiliary P6 branches"],
            "status": "not_supplied_in_p9",
        },
    ]


def flavor_rotation_inputs() -> dict:
    return {
        "gauge_to_mass_basis": "f_gauge_i = (U_f)_{i a} f_mass_a for f in {Q,L,u^c,d^c,nu^c,e^c}",
        "unitary_matrices_required": ["U_Q", "U_L", "U_uc", "U_dc", "U_nuc", "U_ec"],
        "observable_constraints": {
            "CKM": "built from the left quark rotations after electroweak breaking",
            "PMNS": "built from charged-lepton and light-neutrino rotations after seesaw/majorana matching",
        },
        "wilson_rotation_template": (
            "C_phys[a b c d] = U_f1[i a] U_f2[j b] U_f3[k c] "
            "U_f4[l d] C_gauge[i j k l]"
        ),
        "status": "template_only_no_numeric_flavor_fit",
    }


def proton_bound_inputs(p7: dict) -> dict:
    return {
        "p7_current_status": p7["summary"],
        "required_external_inputs": [
            "current experimental lower lifetime limit tau_exp for each proton channel",
            "completed P6 branch and mediator spectrum",
            "physical Wilson basis and chiral contractions",
            "flavor rotations from gauge basis to mass basis",
            "RG evolution factors from matching scale to hadronic scale",
            "lattice/chiral hadronic matrix H_ij^(channel)",
        ],
        "bound_template": (
            "Gamma(p -> a)=C_phys^dagger H^(a) C_phys, "
            "tau_exp(a) > 1/Gamma(p -> a)"
        ),
        "p9_boundary": (
            "P9 supplies the interface and Spin(10) generator data.  It does "
            "not insert numerical current experimental limits or compute lifetimes."
        ),
    }


def build_payload() -> dict:
    p4 = load_json(P4_JSON)
    p6 = load_json(P6_JSON)
    p7 = load_json(P7_JSON)
    p8 = load_json(P8_JSON)
    if p8["summary"]["leading_conditional_field_theory_branch"] != "Spin(10)":
        raise SystemExit("P9 expects P8 to select Spin(10) as the leading branch")

    states = basis_states()
    maps = leptoquark_maps(states) + su2r_maps(states) + off_face_maps(states)
    checks = generator_checks(p4, maps)
    return {
        "description": "Route-C P9 Spin(10)-locked action-completion ledger",
        "boundary": (
            "P9 starts action-level completion for the P8-leading Spin(10) "
            "branch.  It supplies generator and interface data, but it does not "
            "yet choose a full Higgs potential, mediator spectrum, flavor fit, "
            "or numerical proton-bound evaluation."
        ),
        "spin10_lock": {
            "source": "P8 comparative scorecard",
            "leading_conditional_field_theory_branch": "Spin(10)",
            "broad_comparison_closed_for_this_pass": True,
        },
        "basis": {
            "dimension": len(states),
            "decomposition": "16 -> (4,2,1) + (bar4,1,2) under Pati-Salam",
            "states": states,
        },
        "broken_generator_maps": {
            "map_count": len(maps),
            "su4_leptoquark_root_maps": 6,
            "su2r_charged_root_maps": 2,
            "spin10_off_face_clebsch_maps": 24,
            "maps": maps,
        },
        "verification": {
            "basis_dimension_is_16": len(states) == 16,
            "expected_map_count_is_32": len(maps) == 32,
            "p4_transition_coverage": checks,
            "p6_tower_like_uv_forced_before_p9": p6["summary"]["tower_like_uv_forced_now"],
            "p7_physical_proton_bounds_evaluable_before_p9": p7["summary"][
                "physical_proton_bounds_evaluable_now"
            ],
            "p9_generator_layer_passes": len(states) == 16 and len(maps) == 32 and checks[
                "generated_pairs_are_subset_of_p4_broken_pairs"
            ],
        },
        "goldstone_higgs_completion_templates": completion_templates(),
        "mediator_mass_and_coupling_inputs": mediator_inputs(),
        "physical_flavor_rotation_inputs": flavor_rotation_inputs(),
        "proton_bound_inputs": proton_bound_inputs(p7),
        "p9_status": {
            "broken_generator_matrices": "sparse transition maps supplied and checked against P4",
            "goldstone_higgs_completion": "formulae supplied; order parameter branch not chosen",
            "mediator_masses_couplings": "symbolic inputs listed; no spectrum fixed",
            "physical_flavor_rotations": "interface supplied; no flavor benchmark inserted",
            "proton_bound_inputs": "input interface supplied; no numerical lifetime claim",
        },
        "next_stage": {
            "recommended": "P10_choose_spin10_breaking_branch_and_mass_spectrum",
            "required_decision": (
                "Choose a concrete Spin(10)-breaking order parameter/source "
                "branch, then diagonalize broken-vector masses and replay P6/P7 "
                "with explicit masses, couplings, and flavor rotations."
            ),
        },
    }


def write_markdown(payload: dict) -> str:
    ver = payload["verification"]
    coverage = ver["p4_transition_coverage"]
    lines = [
        "# Route-C P9 Spin(10)-Locked Action-Completion Ledger",
        "",
        payload["boundary"],
        "",
        "## Summary",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| basis dimension | {payload['basis']['dimension']} |",
        f"| broken generator sparse maps | {payload['broken_generator_maps']['map_count']} |",
        f"| SU(4) leptoquark root maps | {payload['broken_generator_maps']['su4_leptoquark_root_maps']} |",
        f"| SU(2)_R charged root maps | {payload['broken_generator_maps']['su2r_charged_root_maps']} |",
        f"| Spin(10)/Pati-Salam off-face Clebsch maps | {payload['broken_generator_maps']['spin10_off_face_clebsch_maps']} |",
        f"| P4 broken transition pairs | {coverage['p4_broken_transition_pairs']} |",
        f"| Spin(10) adjoint vector transition pairs | {coverage['spin10_adjoint_vector_pair_count']} |",
        f"| non-adjoint P4 completion pairs | {len(coverage['missing_p4_pairs'])} |",
        f"| generated pairs are subset of P4 broken pairs | {coverage['generated_pairs_are_subset_of_p4_broken_pairs']} |",
        f"| P9 generator layer passes | {ver['p9_generator_layer_passes']} |",
        "",
        "## Generator Families",
        "",
        "| family | maps | role |",
        "| --- | ---: | --- |",
        "| Pati-Salam SU(4)_C leptoquark | 6 | completes Q/L and conjugate-sector leptoquark transitions |",
        "| broken SU(2)_R charged current | 2 | completes u^c-d^c and nu^c-e^c transitions |",
        "| Spin(10)/Pati-Salam off-face (6,2,2) | 24 | connects (4,2,1) to (bar4,1,2) via antisymmetric SU(4) Clebsch maps |",
        "",
        "P9 intentionally does not force every P4 charge-bookkeeping transition to be",
        "a Spin(10) adjoint gauge vector.  The non-adjoint P4 pairs remain",
        "scalar/source/auxiliary completion candidates if a later branch activates",
        "them.",
        "",
        "## Goldstone/Higgs Completion Templates",
        "",
        "| sector | broken generators | status |",
        "| --- | --- | --- |",
    ]
    for row in payload["goldstone_higgs_completion_templates"]:
        lines.append(
            f"| {row['sector']} | {row['broken_generators']} | {row['status']} |"
        )
    lines.extend(
        [
            "",
            "The Ward target for every broken vector remains",
            "",
            "```text",
            "p_mu M^mu(V_X) = M_X M(phi_X)",
            "```",
            "",
            "## Required Symbolic Inputs",
            "",
            "| input | status |",
            "| --- | --- |",
        ]
    )
    for row in payload["mediator_mass_and_coupling_inputs"]:
        lines.append(f"| `{row['symbol']}` | {row['status']} |")
    lines.extend(
        [
            "",
            "## Physical Flavor and Proton-Bound Boundary",
            "",
            f"- flavor status: `{payload['physical_flavor_rotation_inputs']['status']}`",
            f"- proton-bound status: `{payload['proton_bound_inputs']['p9_boundary']}`",
            "",
            "P9 supplies the Spin(10) action-completion interface.  It does not claim",
            "that proton lifetimes or full flavor observables are now evaluable.",
            "",
            "## Next Stage",
            "",
            f"`{payload['next_stage']['recommended']}`",
            "",
            payload["next_stage"]["required_decision"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "basis_dimension": payload["basis"]["dimension"],
                "broken_generator_maps": payload["broken_generator_maps"]["map_count"],
                "covers_all_p4_broken_pairs": payload["verification"][
                    "p4_transition_coverage"
                ]["covers_all_p4_broken_pairs"],
                "p9_generator_layer_passes": payload["verification"][
                    "p9_generator_layer_passes"
                ],
                "next_stage": payload["next_stage"]["recommended"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
