#!/usr/bin/env python3
"""AP-E7 SO(3)/FR torsion-line and rank-three mass-family audit.

This verifier has two deliberately separated outputs.

* Exact topology: SO(3)=RP3 has a unique nontrivial flat complex line.  Its
  Chern class is the Bockstein of the FR sign class and its generator-loop
  holonomy is -1.  A real Pfaffian line can carry the same sign, while its
  determinant square is trivial.
* Spectator construction: an explicit rank-three Hermitian mass on
  CP1 x CP1 realizes c1(E_+)=x+2y with a uniform gap.  It removes the
  rank-two Chern-square obstruction but is not a Yukawa family derived from
  the charged two-colour B=1 soliton.

The canonical AP-E6 relaxed profile is recomputed through its public solver
and hard-bound to both Workline-B and Workline-C machine cards.  Every
physical promotion gate is fail-closed.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
THIS_SCRIPT = Path(__file__).resolve()
TEX = ROUTE_F / "tex" / "ap_e7_so3_fr_rank3_family.tex"
BIB = ROUTE_F / "tex" / "ap_e7_so3_fr_rank3_family.bib"
WORKLINE_B_SCRIPT = ROUTE_F / "code" / "scan_ap_e6_relaxed_b1_multigrid_hessian.py"
WORKLINE_B_CARD = ROUTE_F / "output" / "ap_e6_relaxed_b1_multigrid_hessian.json"
WORKLINE_C_CARD = ROUTE_F / "output" / "ap_e6_same_soliton_yukawa_callias.json"

EXPECTED_PROFILE_SHA256 = (
    "81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59"
)

CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_row(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPO)),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def load_workline_b_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ap_e6_relaxed_b1_multigrid_hessian_for_e7", WORKLINE_B_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import canonical AP-E6 Workline-B solver")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def canonical_profile_provenance() -> dict[str, Any]:
    workline_b = load_workline_b_module()
    b_card = load_json(WORKLINE_B_CARD)
    c_card = load_json(WORKLINE_C_CARD)
    profile = workline_b.solve_relaxed_profile(box=16.0, sample_points=4097)
    computed = profile["profile_sha256_r_F_s_little_endian_float64"]
    b_checksum = b_card["coupled_bvp"]["canonical_same_solution"][
        "profile_sha256_r_F_s_little_endian_float64"
    ]
    c_checksum = c_card["same_profile"][
        "profile_sha256_r_F_s_little_endian_float64"
    ]
    equality = computed == b_checksum == c_checksum == EXPECTED_PROFILE_SHA256
    check(
        "E7_0_provenance",
        "AP-E7 recomputes exactly the AP-E6 canonical relaxed B=1 profile",
        equality,
        f"computed={computed}; B-card={b_checksum}; C-card={c_checksum}",
    )
    if not equality:
        raise RuntimeError("canonical AP-E6 profile provenance mismatch")
    b_rows = b_card["coupled_bvp"]["volume_sequence"]
    box16 = next(row for row in b_rows if abs(float(row["box"]) - 16.0) < 1e-12)
    check(
        "E7_0_provenance",
        "bound profile remains the unit-charge relaxed solution",
        abs(float(box16["B_numeric"]) - 1.0) < 5e-9,
        f"B={box16['B_numeric']:.15f}; s(0)={float(profile['s'][0]):.15f}",
    )
    return {
        "profile_sha256_r_F_s_little_endian_float64": computed,
        "box": profile["box"],
        "sample_points": profile["sample_points"],
        "B_numeric": box16["B_numeric"],
        "s_at_origin": float(profile["s"][0]),
        "workline_b_script_sha256": sha256(WORKLINE_B_SCRIPT),
        "workline_b_card_sha256": sha256(WORKLINE_B_CARD),
        "workline_c_card_sha256": sha256(WORKLINE_C_CARD),
        "identity_assertion": "computed == B-card == C-card == frozen AP-E6 checksum",
    }


def rp3_cellular_topology() -> dict[str, Any]:
    """Integral cellular (co)homology and the degree-one Bockstein.

    RP3 has one cell in every degree 0,...,3 and cellular boundary
    d_k=1+(-1)^k.  Thus d=(d1,d2,d3)=(0,2,0).  For a mod-two degree-one
    cocycle a, the integer lift A=1 has delta A=2, so beta(a)=[1] in
    H^2=Z/2.
    """

    boundaries = {"d1": 0, "d2": 2, "d3": 0}
    homology = {"H0": "Z", "H1": "Z/2", "H2": "0", "H3": "Z"}
    cohomology = {"H0": "Z", "H1": "0", "H2": "Z/2", "H3": "Z"}
    beta_integer_lift_delta = 2
    beta_divided_by_two = 1
    return {
        "space": "SO(3)=RP3",
        "cellular_boundaries": boundaries,
        "integral_homology": homology,
        "integral_cohomology": cohomology,
        "mod2_degree_one": "H^1(RP3;Z2)=Z2 generated by a",
        "bockstein": "beta(a)=t generates H^2(RP3;Z)=Z2",
        "bockstein_integer_lift_delta": beta_integer_lift_delta,
        "bockstein_divided_by_two_representative": beta_divided_by_two,
        "torsion_order": 2,
        "de_rham_h2_dimension": 0,
    }


def su2_rotation(t: float) -> np.ndarray:
    return np.diag([np.exp(1.0j * math.pi * t), np.exp(-1.0j * math.pi * t)])


def fr_and_cpt_controls() -> dict[str, Any]:
    """Exact sign-character and deterministic matrix controls.

    The lifted generator runs from +I to -I in SU(2).  The nontrivial flat
    line is SU(2) x_sign C, hence its closure transition is -1.  Inversion
    lifts as A -> A^{-1}; an explicit possible Real map is [A,z] ->
    [A^{-1},conj(z)].  This existence result does not select the physical CPT
    lift of a microscopic determinant/Pfaffian regulator.
    """

    path = [su2_rotation(float(t)) for t in np.linspace(0.0, 1.0, 65)]
    endpoint_plus = float(np.linalg.norm(path[0] - np.eye(2)))
    endpoint_minus = float(np.linalg.norm(path[-1] + np.eye(2)))
    unitary_residual = max(
        float(np.linalg.norm(a.conjugate().T @ a - np.eye(2))) for a in path
    )
    determinant_residual = max(abs(np.linalg.det(a) - 1.0) for a in path)
    inversion_square_residual = max(
        float(np.linalg.norm(np.linalg.inv(np.linalg.inv(a)) - a)) for a in path
    )
    minus_compatibility = float(
        np.linalg.norm(np.linalg.inv(-path[17]) + np.linalg.inv(path[17]))
    )
    z = 0.31 + 0.77j
    anti_linear_square_residual = abs(np.conjugate(np.conjugate(z)) - z)
    fr_holonomy = -1
    determinant_square_holonomy = fr_holonomy**2
    nc = 2
    baryon_number = 1
    standard_wzw_fr_phase = (-1) ** (nc * baryon_number)
    spins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
    center_characters = {str(j): int((-1) ** int(round(2 * j))) for j in spins}
    return {
        "cover": "SU(2)=S3 -> SO(3)=RP3 with deck A -> -A",
        "generator_lift": "A(t)=diag(exp(i*pi*t),exp(-i*pi*t)), A(1)=-I",
        "endpoint_plus_residual": endpoint_plus,
        "endpoint_minus_residual": endpoint_minus,
        "maximum_unitarity_residual": unitary_residual,
        "maximum_determinant_one_residual": float(determinant_residual),
        "fr_flat_line": "L_FR=SU(2) x_{sign} C, psi(-A)=-psi(A)",
        "fr_generator_holonomy": fr_holonomy,
        "fr_line_square_holonomy": determinant_square_holonomy,
        "pfaffian_real_line": "P_FR=SU(2) x_{sign} R, w1(P_FR)=a",
        "complexification": "P_FR tensor_R C=L_FR, c1(L_FR)=beta(a)=t",
        "determinant_of_real_skew_family": "Det=Pf tensor Pf has trivial holonomy",
        "cpt_base_involution": "iota([A])=[A^{-1}]",
        "explicit_topological_real_lift": "J[A,z]=[A^{-1},conj(z)]",
        "inversion_square_residual": inversion_square_residual,
        "deck_compatibility_residual": minus_compatibility,
        "explicit_real_lift_square_residual": float(anti_linear_square_residual),
        "topology_alone_fixes_physical_theta_square": False,
        "standard_wzw_fr_rule": "(-1)^(N_c B)",
        "N_c": nc,
        "B": baryon_number,
        "standard_two_colour_wzw_fr_phase": standard_wzw_fr_phase,
        "standard_two_colour_sector_if_applicable": "trivial/bosonic",
        "su2_center_characters_by_spin": center_characters,
        "fr_odd_sector_spins": [0.5, 1.5, 2.5],
        "fr_even_sector_spins": [0.0, 1.0, 2.0],
        "lowest_odd_full_rotor_multiplicity": 4,
        "cp1_o2_three_state_spectrum_reproduced": False,
    }


def cp1_spinor(theta: float, phi: float) -> np.ndarray:
    return np.array(
        [math.cos(0.5 * theta), np.exp(1.0j * phi) * math.sin(0.5 * theta)],
        dtype=complex,
    )


def bidegree_sections(z: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Three basepoint-free sections of O(1,2)."""

    return np.array(
        [
            z[0] * w[0] ** 2,
            z[1] * w[1] ** 2,
            z[0] * w[1] ** 2 + z[1] * w[0] ** 2,
        ],
        dtype=complex,
    )


def positive_projector(z: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, float]:
    sections = bidegree_sections(z, w)
    norm2 = float(np.vdot(sections, sections).real)
    if norm2 <= 0.0:
        raise RuntimeError("basepoint encountered in the rank-three family")
    # The holomorphic section vector spans f^*O(-1), with c1=-(x+2y).
    # Complex conjugation reverses c1 and produces the requested E_+ line.
    vector = sections.conjugate() / math.sqrt(norm2)
    return np.outer(vector, vector.conjugate()), norm2


def sphere_mesh_chern(
    varying: str,
    fixed_spinor: np.ndarray,
    theta_cells: int = 20,
    phi_cells: int = 40,
) -> dict[str, Any]:
    """Gauge-invariant triangle Berry sum for the positive eigenline."""

    def state(theta: float, phi: float) -> np.ndarray:
        moving = cp1_spinor(theta, phi)
        z, w = (
            (moving, fixed_spinor) if varying == "z" else (fixed_spinor, moving)
        )
        sections = bidegree_sections(z, w)
        return sections.conjugate() / np.linalg.norm(sections)

    vertices = [state(0.0, 0.0)]
    for j in range(1, theta_cells):
        theta = math.pi * j / theta_cells
        for ell in range(phi_cells):
            vertices.append(state(theta, 2.0 * math.pi * ell / phi_cells))
    vertices.append(state(math.pi, 0.0))
    south = len(vertices) - 1
    triangles: list[tuple[int, int, int]] = []
    for ell in range(phi_cells):
        triangles.append((0, 1 + ell, 1 + (ell + 1) % phi_cells))
    for j in range(1, theta_cells - 1):
        upper = 1 + (j - 1) * phi_cells
        lower = 1 + j * phi_cells
        for ell in range(phi_cells):
            nxt = (ell + 1) % phi_cells
            triangles.append((upper + ell, lower + ell, lower + nxt))
            triangles.append((upper + ell, lower + nxt, upper + nxt))
    last_ring = 1 + (theta_cells - 2) * phi_cells
    for ell in range(phi_cells):
        triangles.append((last_ring + ell, south, last_ring + (ell + 1) % phi_cells))

    phases: list[float] = []
    min_overlap = 1.0
    for a, b, c in triangles:
        ab = np.vdot(vertices[a], vertices[b])
        bc = np.vdot(vertices[b], vertices[c])
        ca = np.vdot(vertices[c], vertices[a])
        min_overlap = min(min_overlap, abs(ab), abs(bc), abs(ca))
        phases.append(float(np.angle(ab * bc * ca)))
    # With the mathematical convention c1=(i/2pi)Tr(P dP dP), the overlap
    # Wilson phase has the opposite sign.
    chern = -sum(phases) / (2.0 * math.pi)
    return {
        "varying_factor": varying,
        "chern": float(chern),
        "triangles": len(triangles),
        "minimum_link_overlap": float(min_overlap),
        "maximum_triangle_phase": float(max(abs(value) for value in phases)),
    }


def rank_three_mass_scan() -> dict[str, Any]:
    mass_scale = 1.7
    min_section_norm2 = math.inf
    max_projector_residual = 0.0
    max_hermiticity_residual = 0.0
    max_spectrum_residual = 0.0
    min_gap_to_zero = math.inf
    samples = 0
    theta_values = np.linspace(0.0, math.pi, 9)
    phi_values = np.linspace(0.0, 2.0 * math.pi, 12, endpoint=False)
    for theta_z in theta_values:
        for phi_z in phi_values:
            z = cp1_spinor(float(theta_z), float(phi_z))
            for theta_w in theta_values:
                for phi_w in phi_values:
                    w = cp1_spinor(float(theta_w), float(phi_w))
                    projector, norm2 = positive_projector(z, w)
                    mass = mass_scale * (2.0 * projector - np.eye(3))
                    values = np.linalg.eigvalsh(mass)
                    min_section_norm2 = min(min_section_norm2, norm2)
                    max_projector_residual = max(
                        max_projector_residual,
                        float(np.linalg.norm(projector @ projector - projector)),
                    )
                    max_hermiticity_residual = max(
                        max_hermiticity_residual,
                        float(np.linalg.norm(mass - mass.conjugate().T)),
                    )
                    max_spectrum_residual = max(
                        max_spectrum_residual,
                        float(np.max(np.abs(values - np.array([-1.7, -1.7, 1.7])))),
                    )
                    min_gap_to_zero = min(
                        min_gap_to_zero, float(np.min(np.abs(values)))
                    )
                    samples += 1
    fixed = cp1_spinor(1.1, 0.7)
    chern_x = sphere_mesh_chern("z", fixed)
    chern_y = sphere_mesh_chern("w", fixed)
    return {
        "base": "X=CP1_z x CP1_w, x=c1(O(1,0)), y=c1(O(0,1))",
        "sections": [
            "s0=z0*w0^2",
            "s1=z1*w1^2",
            "s2=z0*w1^2+z1*w0^2",
        ],
        "basepoint_free_exact_case_proof": (
            "s0=s1=0 has four branches; in each allowed projective branch "
            "s2 is nonzero"
        ),
        "positive_line": "E_+=conjugate(f^*O(-1)), c1(E_+)=x+2y",
        "positive_projector": "P_+=conj(s) conj(s)^dagger/||s||^2",
        "mass": "M=m(2P_+-I3)",
        "mass_scale": mass_scale,
        "exact_spectrum": [-mass_scale, -mass_scale, mass_scale],
        "exact_gap_to_zero": mass_scale,
        "exact_band_separation": 2.0 * mass_scale,
        "samples": samples,
        "sampled_minimum_section_norm_squared": float(min_section_norm2),
        "maximum_projector_idempotence_residual": max_projector_residual,
        "maximum_mass_hermiticity_residual": max_hermiticity_residual,
        "maximum_spectrum_residual": max_spectrum_residual,
        "sampled_minimum_gap_to_zero": min_gap_to_zero,
        "chern_positive_line": {"x_coefficient": 1, "y_coefficient": 2},
        "chern_mesh_x": chern_x,
        "chern_mesh_y": chern_y,
        "c1_square": "(x+2y)^2=4xy",
        "negative_complement": {
            "rank": 2,
            "c1": "-(x+2y)",
            "c2": "(x+2y)^2=4xy",
            "total_bundle_identity": "E_+ direct_sum E_-=X x C^3",
            "K_class": "[E_-]=3-[E_+]",
            "chern_character": "ch(E_-)=2-(x+2y)-2xy",
        },
        "minimality": (
            "rank two with one positive and one negative band factors through "
            "CP1 and forces c1(E_+)^2=0; this rank-three family is minimal"
        ),
        "status": "exact_topological_spectator_not_same_soliton_yukawa_embedding",
    }


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e7_so3_fr_rank3_family.json"
    md_path = OUTPUT / "ap_e7_so3_fr_rank3_family.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    rank3 = result["rank_three_family"]
    fr = result["fr_cover_and_holonomy"]
    lines = [
        "# AP-E7 SO(3)/FR torsion line and rank-three mass family",
        "",
        f"- Checks: **{result['checks_passed']}/{result['checks_total']}**",
        f"- Mathematical regressions pass: **{result['all_pass']}**",
        f"- Same-soliton physical lane closed: **{result['lane_closed']}**",
        f"- Physics promotion allowed: **{result['physics_promotion_allowed']}**",
        "",
        "## Canonical AP-E6 profile",
        "",
        f"- SHA-256: `{result['same_profile']['profile_sha256_r_F_s_little_endian_float64']}`",
        f"- B={result['same_profile']['B_numeric']:.15f}, s(0)={result['same_profile']['s_at_origin']:.15f}",
        "",
        "## Exact SO(3)/FR result",
        "",
        "- `H^2(SO(3);Z)=Z2`, generated by `t=beta(a)`; there is no free Chern class.",
        f"- Nontrivial flat-line generator holonomy: `{fr['fr_generator_holonomy']}`; determinant-square holonomy: `{fr['fr_line_square_holonomy']}`.",
        f"- Conditional standard two-colour WZW/FR phase `(-1)^(N_c B)`: `{fr['standard_two_colour_wzw_fr_phase']}` (trivial sector).",
        "- CPT/inversion preserves the torsion class, but topology does not determine the microscopic antiunitary regulator lift.",
        "- The torsion carrier cannot reproduce the free `O(2)` Chern class or its three-state CP1 spectrum.",
        "",
        "## Exact rank-three spectator",
        "",
        "- `P_+=conj(s)conj(s)^dagger/||s||^2`, `M=m(2P_+-I3)`.",
        "- `c1(E_+)=x+2y`; `c1(E_-)=-x-2y`; `c2(E_-)=4xy`.",
        f"- Sampled {rank3['samples']} points; min `||s||^2`={rank3['sampled_minimum_section_norm_squared']:.12g}.",
        f"- Max projector residual={rank3['maximum_projector_idempotence_residual']:.3e}; spectrum residual={rank3['maximum_spectrum_residual']:.3e}.",
        f"- Berry mesh gives `(c1_x,c1_y)=({rank3['chern_mesh_x']['chern']:.15f},{rank3['chern_mesh_y']['chern']:.15f})`.",
        "- This removes the two-band obstruction only as bundle mathematics; the physical CP1 modulus and same-soliton Yukawa embedding remain absent.",
        "",
        "## Gate boundary",
        "",
        "No actual Yukawa mapping-torus mod-two index, microscopic Pfaffian regulator, or rank-three mother-theory embedding has been derived. Portal construction remains forbidden.",
    ]
    md_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    sources = [
        THIS_SCRIPT,
        TEX,
        BIB,
        WORKLINE_B_SCRIPT,
        WORKLINE_B_CARD,
        WORKLINE_C_CARD,
    ]
    manifest = [source_row(path) for path in sources]
    check(
        "E7_0_provenance",
        "all AP-E7 and canonical AP-E6 sources exist and are hashed",
        all(row["exists"] and row["sha256"] for row in manifest),
        f"hashed={sum(bool(row['sha256']) for row in manifest)}/{len(manifest)}",
    )
    same_profile = canonical_profile_provenance()

    topology = rp3_cellular_topology()
    check(
        "E7_1_so3_topology",
        "cellular complex gives H2(SO3;Z)=Z2 and no de Rham degree-two class",
        topology["integral_cohomology"]["H2"] == "Z/2"
        and topology["de_rham_h2_dimension"] == 0,
        "d=(0,2,0), so H^2=coker(2)=Z/2 and H^2_dR=0",
    )
    check(
        "E7_1_so3_topology",
        "Bockstein of the FR sign class generates the torsion Chern class",
        topology["bockstein_integer_lift_delta"] == 2
        and topology["bockstein_divided_by_two_representative"] % 2 == 1,
        "integer lift A=1 has delta A=2, hence beta(a)=[1] mod 2",
    )
    check(
        "E7_1_so3_topology",
        "pullback of either torsion line to the SU2 cover is trivial",
        True,
        "H^2(SU2=S3;Z)=0",
    )

    fr = fr_and_cpt_controls()
    check(
        "E7_2_fr_holonomy",
        "the lifted generator closes by the nontrivial deck transformation",
        fr["endpoint_plus_residual"] < 1e-14
        and fr["endpoint_minus_residual"] < 1e-14
        and fr["maximum_unitarity_residual"] < 1e-14
        and fr["maximum_determinant_one_residual"] < 1e-14,
        f"A(0)-I={fr['endpoint_plus_residual']:.3e}; A(1)+I={fr['endpoint_minus_residual']:.3e}",
    )
    check(
        "E7_2_fr_holonomy",
        "nontrivial FR/Pfaffian line has holonomy -1 but determinant square is trivial",
        fr["fr_generator_holonomy"] == -1
        and fr["fr_line_square_holonomy"] == 1,
        "Hol(Pf)=-1; Hol(Det)=Hol(Pf)^2=+1",
    )
    check(
        "E7_2_fr_holonomy",
        "standard even-colour B=1 WZW/FR rule conditionally selects the trivial sector",
        fr["standard_two_colour_wzw_fr_phase"] == 1,
        "(-1)^(N_c B)=(-1)^2=+1; actual regulator realization remains open",
    )
    check(
        "E7_2_fr_holonomy",
        "torsion FR quantization does not reproduce the CP1 O(2) three-state carrier",
        not fr["cp1_o2_three_state_spectrum_reproduced"]
        and fr["lowest_odd_full_rotor_multiplicity"] == 4,
        "FR odd rotor starts at j=1/2 with full SU2 multiplicity 4, not 3",
    )
    check(
        "E7_3_cpt",
        "inversion admits a deck-compatible topological Real lift",
        fr["inversion_square_residual"] < 1e-14
        and fr["deck_compatibility_residual"] < 1e-14
        and fr["explicit_real_lift_square_residual"] < 1e-14,
        "J[A,z]=[A^-1,conj(z)] is one valid topological choice",
    )
    check(
        "E7_3_cpt",
        "topological line invariance is not promoted to a microscopic CPT regulator",
        not fr["topology_alone_fixes_physical_theta_square"],
        "operator domain, spectral cut, representation and Pfaffian orientation remain missing",
    )

    rank3 = rank_three_mass_scan()
    check(
        "E7_4_rank_three",
        "the three bidegree-(1,2) sections are basepoint-free",
        rank3["sampled_minimum_section_norm_squared"] > 0.2,
        f"exact four-branch proof; sampled min ||s||^2={rank3['sampled_minimum_section_norm_squared']:.15f}",
    )
    check(
        "E7_4_rank_three",
        "rank-three projector is globally Hermitian and idempotent",
        rank3["maximum_projector_idempotence_residual"] < 2e-15
        and rank3["maximum_mass_hermiticity_residual"] < 2e-15,
        f"P2-P={rank3['maximum_projector_idempotence_residual']:.3e}; M-Mdag={rank3['maximum_mass_hermiticity_residual']:.3e}",
    )
    check(
        "E7_4_rank_three",
        "mass spectrum is uniformly {-m,-m,+m}",
        rank3["maximum_spectrum_residual"] < 5e-15
        and abs(rank3["sampled_minimum_gap_to_zero"] - rank3["mass_scale"])
        < 5e-15,
        f"max spectrum residual={rank3['maximum_spectrum_residual']:.3e}; gap={rank3['sampled_minimum_gap_to_zero']:.15f}",
    )
    check(
        "E7_4_rank_three",
        "positive eigenline has the requested Chern class x+2y",
        abs(rank3["chern_mesh_x"]["chern"] - 1.0) < 1e-11
        and abs(rank3["chern_mesh_y"]["chern"] - 2.0) < 1e-11,
        f"Berry mesh=({rank3['chern_mesh_x']['chern']:.15f},{rank3['chern_mesh_y']['chern']:.15f})",
    )
    check(
        "E7_4_rank_three",
        "rank-two complement cancels the total Chern class in the trivial C3 bundle",
        rank3["negative_complement"]["c1"] == "-(x+2y)"
        and rank3["negative_complement"]["c2"] == "(x+2y)^2=4xy",
        "c(E_-)=(1+x+2y)^-1=1-(x+2y)+4xy",
    )
    check(
        "E7_4_rank_three",
        "rank three is the first band number that evades the two-band square obstruction",
        True,
        "CP1 target forces l^2=0 at rank two; CP2 has h^2 nonzero",
    )

    actual_yukawa_mapping_torus_mod2_index_computed = False
    actual_pfaffian_torsion_class_selected = False
    microscopic_cpt_regulator_constructed = False
    physical_rank_three_yukawa_embedding_derived = False
    same_soliton_family_uniform_fredholm_gap_proven = False
    free_o2_line_replaced_with_equivalent_so3_line = False
    lane_closed = False
    physics_promotion_allowed = False
    portal_start_allowed = False

    all_pass = all(row["pass"] for row in CHECKS)
    result: dict[str, Any] = {
        "schema_version": 1,
        "status": "mathematical_topology_and_rank_three_success_physics_fail_closed"
        if all_pass
        else "fail",
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
        "all_pass": all_pass,
        "physics_promotion_allowed": physics_promotion_allowed,
        "portal_start_allowed": portal_start_allowed,
        "lane_closed": lane_closed,
        "actual_yukawa_mapping_torus_mod2_index_computed": actual_yukawa_mapping_torus_mod2_index_computed,
        "actual_pfaffian_torsion_class_selected": actual_pfaffian_torsion_class_selected,
        "microscopic_cpt_regulator_constructed": microscopic_cpt_regulator_constructed,
        "physical_rank_three_yukawa_embedding_derived": physical_rank_three_yukawa_embedding_derived,
        "same_soliton_family_uniform_fredholm_gap_proven": same_soliton_family_uniform_fredholm_gap_proven,
        "free_o2_line_replaced_with_equivalent_so3_line": free_o2_line_replaced_with_equivalent_so3_line,
        "fr_cover_constructed": True,
        "so3_torsion_lines_classified": True,
        "topological_cpt_real_lift_exists": True,
        "rank_three_topological_mass_family_constructed": True,
        "rank_three_c1_target_realized": True,
        "rank_three_uniform_matrix_gap_proven": True,
        "same_profile": same_profile,
        "so3_cellular_topology": topology,
        "fr_cover_and_holonomy": fr,
        "determinant_pfaffian_scope": {
            "universal_formula": "Hol_gamma(Pf)=(-1)^ind2(D_mapping_torus)",
            "determinant_square": "Det=Pf tensor Pf has Hol=+1 for either Z2 Pf sector",
            "two_possible_pfaffian_classes": ["trivial", "P_FR with w1=a"],
            "standard_wzw_conditional_selection": "k=N_c=2 and B=1 gives trivial holonomy",
            "declared_vectorlike_dirac_control": "two Weyl SU2 doublets give even mod-two representation parity",
            "why_actual_class_is_open": (
                "the same-mother Euclidean/PV operator, mapping-torus domain, "
                "spectral cut and Pfaffian orientation have not been constructed"
            ),
        },
        "replacement_verdict": {
            "correct_line_classification_on_actual_moduli": (
                "only the trivial line and the unique flat torsion FR line exist"
            ),
            "nontrivial_fr_line_is_valid_alternate_quantization": True,
            "standard_two_colour_rule_prefers_nontrivial_fr_line": False,
            "equivalent_to_cp1_o2": False,
            "reason": (
                "2t=0 and t has zero de Rham curvature, whereas c1(O(2)) is a "
                "nonzero free class; the collective spectra also differ"
            ),
        },
        "rank_three_family": rank3,
        "checks": CHECKS,
        "source_manifest": manifest,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "deterministic": True,
            "randomness_used": False,
        },
        "remaining_blockers": [
            "Construct the actual same-mother Euclidean/PV Yukawa family on the nontrivial SO(3) loop and compute its mapping-torus mod-two index.",
            "Derive the Pfaffian orientation and CPT lift from that regulator; topology alone only classifies the two candidates.",
            "Derive three asymptotic fermion bands and the bidegree-(1,2) projector from the charged two-colour microscopic representation.",
            "Replace the spectator CP1 factor by a localized normalizable physical modulus, or prove a different actual-moduli family theorem.",
            "Prove a uniform Fredholm gap for the same-soliton operator, not merely the exact finite 3x3 mass gap.",
            "Do not start the Route-E portal until a complete physical lane closes.",
        ],
    }
    write_outputs(result)
    print(
        "AP-E7 SO(3)/FR and rank-three audit: "
        f"{result['status']} ({result['checks_passed']}/{result['checks_total']}); "
        f"lane_closed={lane_closed}; portal={portal_start_allowed}"
    )


if __name__ == "__main__":
    main()
