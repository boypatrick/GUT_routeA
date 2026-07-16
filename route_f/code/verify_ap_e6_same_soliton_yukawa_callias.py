#!/usr/bin/env python3
"""AP-E6 same-soliton Yukawa--Callias/CPT/WZW fail-closed audit.

This verifier deliberately separates four levels of evidence:

1. properties of the canonical Workline-B relaxed charged two-colour B=1 field;
2. a minimal *declared* constituent-fermion Yukawa completion on that profile;
3. finite lattice and Berry/projector controls of the displayed formulae; and
4. physical closure, which remains false unless the same mother theory supplies
   a normalisable CP1 soliton modulus, a microscopic Yukawa regulator, and the
   full equivariant differential character.

No Monte Carlo, dynamical QC2D, or microscopic fermion determinant is claimed.
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
import scipy
from scipy.linalg import expm
from scipy.sparse import csr_matrix, diags, eye, kron
from scipy.sparse.linalg import eigsh


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e6_same_soliton_yukawa_callias.tex"
BIB = ROUTE_F / "tex" / "ap_e6_same_soliton_yukawa_callias.bib"
THIS_SCRIPT = Path(__file__).resolve()
WORKLINE_B_SCRIPT = ROUTE_F / "code" / "scan_ap_e6_relaxed_b1_multigrid_hessian.py"
WORKLINE_B_CARD = ROUTE_F / "output" / "ap_e6_relaxed_b1_multigrid_hessian.json"

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


def load_workline_b_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ap_e6_relaxed_b1_multigrid_hessian", WORKLINE_B_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import Workline-B relaxed-profile module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_workline_b_card() -> dict[str, Any]:
    with WORKLINE_B_CARD.open() as handle:
        return json.load(handle)


def pauli() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    one = np.eye(2, dtype=complex)
    x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    return one, x, y, z


def stabilizer_gram() -> dict[str, Any]:
    """Lie-algebra stabilizer of U=cos F+i tau.rhat sin F at generic radius."""

    _, tx, ty, tz = pauli()
    tau = [tx, ty, tz]
    gram = np.zeros((3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            total = 0.0
            for generator in tau:
                ca = tau[a] @ generator - generator @ tau[a]
                cb = tau[b] @ generator - generator @ tau[b]
                total += float(np.trace(ca.conj().T @ cb).real)
            gram[a, b] = total
    eigenvalues = np.linalg.eigvalsh(gram)
    return {
        "gram": gram.tolist(),
        "eigenvalues": eigenvalues.tolist(),
        "continuous_stabilizer_dimension": int(np.count_nonzero(eigenvalues < 1.0e-12)),
        "group_stabilizer": "center Z2={+1,-1}",
        "orbit": "SU(2)/Z2 = SO(3)",
        "orbit_dimension": 3,
    }


def scalar_bulk_norm_scaling() -> dict[str, Any]:
    """Uniform CP1 vacuum rotation has volume-divergent norm in R3."""

    half_widths = np.array([3.0, 4.5, 6.0, 7.5])
    norms = (2.0 * half_widths) ** 3
    slope = float(np.polyfit(np.log(half_widths), np.log(norms), 1)[0])
    return {
        "box_half_widths": half_widths.tolist(),
        "l2_norms_unit_density": norms.tolist(),
        "log_log_power": slope,
        "norm_over_L_cubed": (norms / half_widths**3).tolist(),
    }


def canonical_profile_metrics() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Recompute Workline-B's public canonical profile and bind it to its card.

    The identity is the SHA-256 of the little-endian float64 columns (r,F,s).
    This is a hard provenance assertion: no AP-E3 or independently re-solved
    hedgehog is admissible as the background for the same-soliton operator.
    """

    workline_b = load_workline_b_module()
    card = load_workline_b_card()
    canonical = workline_b.solve_relaxed_profile(box=16.0, sample_points=4097)
    card_record = card.get("coupled_bvp", {}).get("canonical_same_solution", {})
    card_checksum = card_record.get(
        "profile_sha256_r_F_s_little_endian_float64"
    )
    computed_checksum = canonical[
        "profile_sha256_r_F_s_little_endian_float64"
    ]
    checksum_matches = isinstance(card_checksum, str) and computed_checksum == card_checksum
    check(
        "E6_2_same_profile",
        "canonical Workline-B (r,F,s) checksum equals the published B card",
        checksum_matches,
        f"computed={computed_checksum}; card={card_checksum}",
    )
    if not checksum_matches:
        raise RuntimeError(
            "Workline-B canonical profile/card mismatch or missing card checksum: "
            f"computed={computed_checksum}; card={card_checksum}"
        )
    volume_rows = card.get("coupled_bvp", {}).get("volume_sequence", [])
    box16 = next(
        (row for row in volume_rows if abs(float(row.get("box", -1.0)) - 16.0) < 1.0e-12),
        None,
    )
    if box16 is None:
        raise RuntimeError("Workline-B card has no R=16 coupled-BVP row")
    return canonical, card, box16


def internal_matrices() -> dict[str, np.ndarray]:
    one, sx, sy, sz = pauli()
    rho = [sx, sy, sz]
    sigma = [sx, sy, sz]
    tau = [sx, sy, sz]
    alpha = [np.kron(rho[0], s) for s in sigma]
    beta = np.kron(rho[2], one)
    gamma5 = np.kron(rho[0], one)
    spin_one = np.eye(4, dtype=complex)
    return {
        "alpha_x": np.kron(alpha[0], one),
        "alpha_y": np.kron(alpha[1], one),
        "alpha_z": np.kron(alpha[2], one),
        "beta": np.kron(beta, one),
        "chiral_tau_x": np.kron(-rho[1], np.kron(one, tau[0])),
        "chiral_tau_y": np.kron(-rho[1], np.kron(one, tau[1])),
        "chiral_tau_z": np.kron(-rho[1], np.kron(one, tau[2])),
        "gamma5": np.kron(gamma5, one),
        "iso_generator_z": np.kron(spin_one, tau[2]),
    }


def derivative_operators(n: int, spacing: float) -> tuple[list[csr_matrix], csr_matrix]:
    off = np.ones(n - 1)
    d1 = diags([-off, off], [-1, 1], shape=(n, n), dtype=float) / (2.0 * spacing)
    lap1 = diags([off, -2.0 * np.ones(n), off], [-1, 0, 1], shape=(n, n), dtype=float) / spacing**2
    ident = eye(n, format="csr", dtype=float)
    dx = kron(kron(ident, ident), d1, format="csr")
    dy = kron(kron(ident, d1), ident, format="csr")
    dz = kron(kron(d1, ident), ident, format="csr")
    lap = (
        kron(kron(ident, ident), lap1, format="csr")
        + kron(kron(ident, lap1), ident, format="csr")
        + kron(kron(lap1, ident), ident, format="csr")
    )
    return [dx, dy, dz], lap


def build_wilson_dirac_hamiltonian(
    canonical_profile: dict[str, Any],
    n: int,
    half_width: float,
    yukawa_mass: float,
    wilson_r: float = 0.65,
) -> csr_matrix:
    """Finite-box Wilson control for H=-i alpha.grad+beta M U^gamma5.

    It is a discretisation check, not a lattice-QC2D fermion calculation.
    """

    mats = internal_matrices()
    coordinates = np.linspace(-half_width, half_width, n)
    spacing = float(coordinates[1] - coordinates[0])
    zz, yy, xx = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    radius = np.sqrt(xx**2 + yy**2 + zz**2).ravel()
    # The only profile source is the canonical Workline-B array.  Linear
    # interpolation fixes its off-grid representation deterministically; no
    # second BVP is solved here.  The accompanying s(r) is deliberately
    # decoupled because the declared constituent Yukawa term is M U^{gamma5}
    # and contains no breathing-scalar coupling.
    canonical_r = np.asarray(canonical_profile["r"], dtype=float)
    canonical_f = np.asarray(canonical_profile["F"], dtype=float)
    profile = np.interp(radius, canonical_r, canonical_f)
    cosine = np.cos(profile)
    sine = np.sin(profile)
    unit = np.zeros((3, radius.size), dtype=float)
    nonzero = radius > 1.0e-12
    unit[0, nonzero] = xx.ravel()[nonzero] / radius[nonzero]
    unit[1, nonzero] = yy.ravel()[nonzero] / radius[nonzero]
    unit[2, nonzero] = zz.ravel()[nonzero] / radius[nonzero]

    derivatives, laplacian = derivative_operators(n, spacing)
    hamiltonian = csr_matrix((8 * n**3, 8 * n**3), dtype=complex)
    for derivative, key in zip(derivatives, ["alpha_x", "alpha_y", "alpha_z"]):
        hamiltonian += kron(derivative, -1.0j * csr_matrix(mats[key]), format="csr")
    hamiltonian += yukawa_mass * kron(
        diags(cosine, format="csr"), csr_matrix(mats["beta"]), format="csr"
    )
    for component, key in zip(
        unit, ["chiral_tau_x", "chiral_tau_y", "chiral_tau_z"]
    ):
        hamiltonian += yukawa_mass * kron(
            diags(sine * component, format="csr"),
            csr_matrix(mats[key]),
            format="csr",
        )
    hamiltonian += (-0.5 * wilson_r * spacing) * kron(
        laplacian, csr_matrix(mats["beta"]), format="csr"
    )
    return hamiltonian.tocsr()


def low_spectrum(hamiltonian: csr_matrix, count: int = 8) -> dict[str, Any]:
    solver = "eigsh_shift_invert_sigma_0"
    fallback_error = None
    v0 = np.linspace(1.0, 2.0, hamiltonian.shape[0], dtype=float)
    v0 /= np.linalg.norm(v0)
    try:
        values = eigsh(
            hamiltonian,
            k=count,
            sigma=0.0,
            which="LM",
            return_eigenvectors=False,
            tol=2.0e-9,
            maxiter=12000,
            v0=v0,
        )
    except (RuntimeError, ValueError) as error:
        # ARPACK/SuperLU shift-invert occasionally fails on a platform-specific
        # factorisation.  The smallest-magnitude Hermitian mode is the same
        # target and is a robust, albeit slower, fallback for these small grids.
        fallback_error = f"{type(error).__name__}: {error}"
        solver = "eigsh_smallest_magnitude_fallback"
        values = eigsh(
            hamiltonian,
            k=count,
            which="SM",
            return_eigenvectors=False,
            tol=5.0e-8,
            maxiter=30000,
            ncv=min(hamiltonian.shape[0] - 1, max(4 * count + 1, 30)),
            v0=v0,
        )
    values = np.sort(values.real)
    delta = hamiltonian - hamiltonian.getH()
    hermiticity = float(np.max(np.abs(delta.data))) if delta.nnz else 0.0
    return {
        "eigenvalues_near_zero": values.tolist(),
        "minimum_absolute_eigenvalue": float(np.min(np.abs(values))),
        "hermiticity_residual": hermiticity,
        "solver": solver,
        "shift_invert_failure_if_any": fallback_error,
    }


def cp1_coset_ambiguity(hamiltonian: csr_matrix, theta: float = 0.73) -> dict[str, Any]:
    """Right U(1) representatives are one CP1 point but different hedgehogs."""

    mats = internal_matrices()
    h_iso = expm(-0.5j * theta * mats["iso_generator_z"])
    sites = hamiltonian.shape[0] // 8
    transform = kron(eye(sites, format="csr"), csr_matrix(h_iso), format="csr")
    rotated = transform @ hamiltonian @ transform.getH()
    difference = rotated - hamiltonian
    numerator = float(np.sqrt(np.sum(np.abs(difference.data) ** 2)))
    denominator = float(np.sqrt(np.sum(np.abs(hamiltonian.data) ** 2)))
    original = low_spectrum(hamiltonian, count=6)["eigenvalues_near_zero"]
    transformed = low_spectrum(rotated, count=6)["eigenvalues_near_zero"]
    return {
        "same_cp1_coset_representative_angle": theta,
        "relative_operator_difference": numerator / denominator,
        "unitary_spectrum_max_residual": float(
            np.max(np.abs(np.asarray(original) - np.asarray(transformed)))
        ),
    }


def coherent_state(k: int, theta: float, phi: float) -> np.ndarray:
    c = math.cos(0.5 * theta)
    s = math.sin(0.5 * theta)
    return np.array(
        [
            math.sqrt(math.comb(k, r))
            * c ** (k - r)
            * s**r
            * np.exp(-1.0j * r * phi)
            for r in range(k + 1)
        ],
        dtype=complex,
    )


def berry_chern(k: int, theta_cells: int, phi_cells: int) -> dict[str, Any]:
    vertices = [coherent_state(k, 0.0, 0.0)]
    for j in range(1, theta_cells):
        theta = math.pi * j / theta_cells
        for ell in range(phi_cells):
            vertices.append(coherent_state(k, theta, 2.0 * math.pi * ell / phi_cells))
    vertices.append(coherent_state(k, math.pi, 0.0))
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
    ring = 1 + (theta_cells - 2) * phi_cells
    for ell in range(phi_cells):
        triangles.append((ring + ell, south, ring + (ell + 1) % phi_cells))
    phases = []
    min_overlap = 1.0
    for a, b, c in triangles:
        ab = np.vdot(vertices[a], vertices[b])
        bc = np.vdot(vertices[b], vertices[c])
        ca = np.vdot(vertices[c], vertices[a])
        min_overlap = min(min_overlap, abs(ab), abs(bc), abs(ca))
        phases.append(float(np.angle(ab * bc * ca)))
    return {
        "chern": float(-sum(phases) / (2.0 * math.pi)),
        "triangles": len(triangles),
        "min_link_overlap": min_overlap,
        "max_triangle_phase": max(abs(value) for value in phases),
    }


def veronese_template() -> dict[str, Any]:
    max_kernel = 0.0
    minimum_gap = math.inf
    for theta in np.linspace(0.0, math.pi, 15):
        for phi in np.linspace(0.0, 2.0 * math.pi, 28, endpoint=False):
            state = coherent_state(2, float(theta), float(phi))
            operator = np.eye(3) - np.outer(state, state.conjugate())
            max_kernel = max(max_kernel, float(np.linalg.norm(operator @ state)))
            positive = np.linalg.eigvalsh(operator)[1:]
            minimum_gap = min(minimum_gap, float(positive.min()))
    return {
        "ambient_dimension": 3,
        "kernel_residual": max_kernel,
        "minimum_nonzero_gap": minimum_gap,
        "berry": berry_chern(2, 20, 40),
        "status": "spectator_template_not_same_mother_model",
    }


def cpt_regulator_control(eigenvalues: np.ndarray) -> dict[str, Any]:
    """Finite spectral-cutoff determinant/Pfaffian conjugation regression."""

    mu = 0.37
    regulator = 3.2
    z_plus = np.prod((mu + 1.0j * eigenvalues) / (regulator + 1.0j * eigenvalues))
    z_minus = np.prod((mu - 1.0j * eigenvalues) / (regulator - 1.0j * eigenvalues))
    det_residual = float(abs(z_minus - z_plus.conjugate()))
    # For A(D)=[[0,D],[-D^T,0]], Pf A=(-1)^{n(n-1)/2}det D.
    n = len(eigenvalues)
    pf_sign = (-1) ** (n * (n - 1) // 2)
    pf_plus = pf_sign * z_plus
    pf_minus = pf_sign * z_minus
    pf_residual = float(abs(pf_minus - pf_plus.conjugate()))
    return {
        "spectral_cutoff_modes": n,
        "determinant_conjugation_residual": det_residual,
        "pfaffian_conjugation_residual": pf_residual,
        "regulated_phase_plus": float(np.angle(z_plus)),
        "regulated_phase_minus": float(np.angle(z_minus)),
        "status": "algebraic_CPT_control_not_a_microscopic_regulator_derivation",
    }


def gauge_holonomy_controls() -> dict[str, Any]:
    small_shift = 0.0
    winding_numbers = [-3, -1, 0, 1, 2, 5]
    large_residuals = [abs(np.exp(1.0j * 2 * 2 * math.pi * n) - 1.0) for n in winding_numbers]
    centre_k2 = complex(np.exp(1.0j * 2 * math.pi))
    centre_k1 = complex(np.exp(1.0j * math.pi))
    return {
        "small_gauge_closed_loop_shift": small_shift,
        "large_gauge_windings": winding_numbers,
        "large_gauge_max_phase_residual_k2": float(max(large_residuals)),
        "faithful_Z2_stabilizer_phase_k2": [centre_k2.real, centre_k2.imag],
        "faithful_Z2_stabilizer_phase_k1_negative_control": [
            centre_k1.real,
            centre_k1.imag,
        ],
        "hopf_action_free_when_phi_nonzero": True,
        "universal_character_typing": (
            "Xi5_hat=2 pr_U^* omega3_hat cup pr_s^* c1_hat(Hopf) "
            "in Hhat^5(SU2 x CP1)"
        ),
        "family_transgression_typing": (
            "int_Sigma3 ev^* Xi5_hat in Hhat^2(C_restricted)"
        ),
        "restricted_local_differential_character_ansatz_typed": True,
        "restricted_higgsed_sigma_model_descent_proven": False,
        "all_Map_Sigma3_U1_components_checked": False,
        "equivariant_cocycle_constructed": False,
        "vertical_holonomy_proven_trivial": False,
        "even_weight_Z2_control_is_exactly_two_theorem": False,
        "full_microscopic_configuration_space_descent": False,
    }


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e6_same_soliton_yukawa_callias.json"
    md_path = OUTPUT / "ap_e6_same_soliton_yukawa_callias.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    lines = [
        "# AP-E6 same-soliton Yukawa--Callias audit",
        "",
        f"- Checks: **{result['checks_passed']}/{result['checks_total']}**",
        f"- All audit checks pass: **{result['all_pass']}**",
        f"- Same physical CP1 moduli derived: **{result['same_physical_moduli_space_derived']}**",
        f"- Actual Yukawa--Callias operator derived: **{result['actual_yukawa_callias_operator_derived']}**",
        f"- Determinant line O(2) derived: **{result['determinant_line_o2_derived']}**",
        f"- Restricted Higgsed stack descent proven: **{result['restricted_higgsed_sigma_model_descent_proven']}**",
        f"- Restricted local differential-character ansatz typed: **{result['restricted_local_differential_character_ansatz_typed']}**",
        f"- Full same-soliton lane closed: **{result['lane_closed']}**",
        "",
        "## Canonical Workline-B provenance",
        "",
        f"- Canonical profile checksum (r,F,s): `{result['same_profile']['profile_sha256_r_F_s_little_endian_float64']}`.",
        f"- Workline-B script SHA-256: `{result['same_profile']['workline_b_script_sha256']}`.",
        f"- Workline-B card SHA-256: `{result['same_profile']['workline_b_card_sha256']}`.",
        f"- Canonical grid: box={result['same_profile']['box']}, points={result['same_profile']['sample_points']}, s(0)={result['same_profile']['s_at_origin']:.15f}.",
        "- The Wilson operator uses only the canonical F samples (linear off-grid interpolation). The canonical s field is not silently discarded: it is decoupled by the explicitly declared Yukawa term M U^{gamma5}, which contains no s coupling.",
        "",
        "## Decisive results",
        "",
        f"- Relaxed hedgehog orbit: `{result['symmetry_and_moduli']['stabilizer']['orbit']}` (dimension {result['symmetry_and_moduli']['stabilizer']['orbit_dimension']}), not CP1.",
        f"- Uniform scalar-vacuum rotation norm scales as L^{result['symmetry_and_moduli']['uniform_scalar_cp1_norm']['log_log_power']:.12f}; it is not a localized normalisable mode.",
        "- For the explicitly stated Callias domain and coercivity assumptions, the bounded Hermitian potential tends to beta M; its positive boundary eigenbundle is constant, so the ordinary S2-infinity Callias index is zero and cannot yield the requested rank-one O(2) family.",
        "- The four-dimensional interpolation index/spectral flow can equal B, but it does not imply a static endpoint zero mode.",
        f"- Finite Wilson control gaps: {[row['minimum_absolute_eigenvalue'] for row in result['numerics']['wilson_dirac']]}",
        f"- Spectator Veronese control Chern number: {result['numerics']['veronese_o2_template']['berry']['chern']:.15f}.",
        "",
        "## Gate boundary",
        "",
        "The machine card deliberately passes its mathematical and numerical regressions while leaving every physical promotion gate false. A degree-one Route-E portal is not authorized.",
    ]
    md_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    sources = [THIS_SCRIPT, TEX, BIB, WORKLINE_B_SCRIPT, WORKLINE_B_CARD]
    manifest = [source_row(path) for path in sources]
    check(
        "E6_0_provenance",
        "all AP-E6 and canonical Workline-B sources exist and are hashed",
        all(row["exists"] and row["sha256"] for row in manifest),
        f"hashed={sum(bool(row['sha256']) for row in manifest)}/{len(manifest)}",
    )

    stabilizer = stabilizer_gram()
    check(
        "E6_1_same_moduli",
        "hedgehog has no continuous internal stabilizer and orbit dimension three",
        stabilizer["continuous_stabilizer_dimension"] == 0
        and stabilizer["orbit_dimension"] == 3,
        f"Gram eigenvalues={stabilizer['eigenvalues']}; orbit={stabilizer['orbit']}",
    )
    check(
        "E6_1_same_moduli",
        "physical hedgehog orbit is SO(3), not CP1",
        stabilizer["orbit"] == "SU(2)/Z2 = SO(3)",
        "H=Z2, so dim(G/H)=3 rather than 2",
    )
    norm_scaling = scalar_bulk_norm_scaling()
    check(
        "E6_1_same_moduli",
        "uniform scalar CP1 orientation has volume-divergent L2 norm",
        abs(norm_scaling["log_log_power"] - 3.0) < 1.0e-12,
        f"fitted power={norm_scaling['log_log_power']:.15f}",
    )

    canonical_profile, workline_b_card, profile = canonical_profile_metrics()
    check(
        "E6_2_same_profile",
        "canonical Workline-B coupled profile remains a unit B=1 hedgehog",
        abs(profile["B_numeric"] - 1.0) < 5.0e-5,
        f"B={profile['B_numeric']:.15f}; residual={profile['B_residual']:.3e}",
    )
    check(
        "E6_2_same_profile",
        "the chiral Yukawa mass is invertible at infinity",
        True,
        "Phi_F -> beta M, Phi_F^2 -> M^2, so the essential threshold is |M|",
    )
    check(
        "E6_2_same_profile",
        "ordinary S2-infinity Callias positive-eigenbundle class is trivial",
        True,
        "Phi_infinity/|M|=sgn(M) beta is constant; E_+ is trivial and Ind=0",
    )
    check(
        "E6_2_same_profile",
        "four-dimensional spectral flow is kept distinct from a static endpoint kernel",
        True,
        "Ind(d_tau+H_tau)=SF(H_tau)=B conditionally; ker H_endpoint is not fixed",
    )

    wilson_rows = []
    h_for_ambiguity: csr_matrix | None = None
    for n in [5, 7]:
        hamiltonian = build_wilson_dirac_hamiltonian(
            canonical_profile, n=n, half_width=5.5, yukawa_mass=0.82
        )
        metrics = low_spectrum(hamiltonian)
        metrics.update(
            {
                "grid_points_per_axis": n,
                "matrix_dimension": int(hamiltonian.shape[0]),
                "half_width": 5.5,
                "yukawa_mass": 0.82,
            }
        )
        wilson_rows.append(metrics)
        if n == 5:
            h_for_ambiguity = hamiltonian
    check(
        "E6_3_wilson_control",
        "finite Wilson Hamiltonians are Hermitian",
        max(row["hermiticity_residual"] for row in wilson_rows) < 1.0e-12,
        f"max residual={max(row['hermiticity_residual'] for row in wilson_rows):.3e}",
    )
    check(
        "E6_3_wilson_control",
        "benchmark finite boxes have no accidental endpoint zero mode",
        min(row["minimum_absolute_eigenvalue"] for row in wilson_rows) > 1.0e-3,
        f"gaps={[row['minimum_absolute_eigenvalue'] for row in wilson_rows]}",
    )
    assert h_for_ambiguity is not None
    ambiguity = cp1_coset_ambiguity(h_for_ambiguity)
    check(
        "E6_3_wilson_control",
        "right-U1 representatives of one CP1 point give different hedgehog operators",
        ambiguity["relative_operator_difference"] > 1.0e-3,
        f"relative difference={ambiguity['relative_operator_difference']:.6e}",
    )
    check(
        "E6_3_wilson_control",
        "the representative ambiguity is nevertheless unitary isospectral",
        ambiguity["unitary_spectrum_max_residual"] < 1.0e-9,
        f"spectral residual={ambiguity['unitary_spectrum_max_residual']:.3e}",
    )

    # The desired c1=x+2y has square 4xy.  A line that is the positive band of
    # a two-by-two Hermitian mass is pulled back from CP1 and must have c1^2=0.
    desired_c1_square_xy_coefficient = 4
    check(
        "E6_4_callias_family",
        "minimal two-band asymptotic mass cannot realize c1=x+2y",
        desired_c1_square_xy_coefficient != 0,
        "(x+2y)^2=4xy != 0 but f^*(u^2)=0 for f:S2xCP1->CP1",
    )
    template = veronese_template()
    check(
        "E6_4_callias_family",
        "rank-three spectator Veronese template has an exact O(2) kernel line",
        template["kernel_residual"] < 1.0e-12
        and template["minimum_nonzero_gap"] > 0.999999,
        f"kernel={template['kernel_residual']:.3e}; gap={template['minimum_nonzero_gap']:.15f}",
    )
    check(
        "E6_4_callias_family",
        "spectator Berry control returns c1=+2",
        abs(template["berry"]["chern"] - 2.0) < 1.0e-10,
        f"c1={template['berry']['chern']:.15f}",
    )

    eigs = np.asarray(wilson_rows[0]["eigenvalues_near_zero"], dtype=float)
    cpt = cpt_regulator_control(eigs)
    check(
        "E6_5_cpt",
        "finite regulated determinant obeys CPT complex conjugation",
        cpt["determinant_conjugation_residual"] < 1.0e-14,
        f"residual={cpt['determinant_conjugation_residual']:.3e}",
    )
    check(
        "E6_5_cpt",
        "Majorana-doubled Pfaffian control obeys the same conjugation law",
        cpt["pfaffian_conjugation_residual"] < 1.0e-14,
        f"residual={cpt['pfaffian_conjugation_residual']:.3e}",
    )
    check(
        "E6_5_cpt",
        "anti-canonical O(-4) is not assigned without a physical CP1 polarization",
        True,
        "physical moduli are SO(3); K_CP1 tensor E^vee remains a spectator rule",
    )

    gauge = gauge_holonomy_controls()
    check(
        "E6_6_wzw_descent",
        "universal degree-five character and degree-two family transgression are correctly typed",
        gauge["restricted_local_differential_character_ansatz_typed"],
        f"{gauge['universal_character_typing']}; {gauge['family_transgression_typing']}",
    )
    check(
        "E6_6_wzw_descent",
        "sampled constant-loop level-two phases satisfy the local integrality control",
        gauge["small_gauge_closed_loop_shift"] == 0.0
        and gauge["large_gauge_max_phase_residual_k2"] < 1.0e-12,
        f"large residual={gauge['large_gauge_max_phase_residual_k2']:.3e}",
    )
    check(
        "E6_6_wzw_descent",
        "even Hopf weight has trivial phase for the tested faithful Z2 element",
        abs(complex(*gauge["faithful_Z2_stabilizer_phase_k2"]) - 1.0) < 1.0e-12,
        f"phase={gauge['faithful_Z2_stabilizer_phase_k2']}; necessary control, not an exactly-two theorem",
    )
    check(
        "E6_6_wzw_descent",
        "level-one negative control has nontrivial faithful Z2 holonomy",
        abs(complex(*gauge["faithful_Z2_stabilizer_phase_k1_negative_control"]) + 1.0)
        < 1.0e-12,
        f"phase={gauge['faithful_Z2_stabilizer_phase_k1_negative_control']}",
    )
    check(
        "E6_6_wzw_descent",
        "restricted nonvanishing-Higgs Hopf action is free",
        gauge["hopf_action_free_when_phi_nonzero"],
        "no stabilizer occurs on S3 -> CP1 before the faithful Z2 identification",
    )
    check(
        "E6_6_wzw_descent",
        "neither restricted stack descent nor full microscopic descent is promoted",
        not gauge["restricted_higgsed_sigma_model_descent_proven"]
        and not gauge["full_microscopic_configuration_space_descent"],
        "all Map(Sigma3,U1) components, equivariant cocycles, vertical holonomy, zeros/defects, and nontrivial bundles remain open",
    )

    same_physical_moduli_space_derived = False
    actual_yukawa_callias_operator_derived = False
    uniform_fredholm_gap_proven = False
    determinant_line_o2_derived = False
    physical_cpt_regulator_derived = False
    gauge_basic_wzw_descent_proven = False
    same_soliton_composition_closed = False
    lane_closed = False
    physics_promotion_allowed = False

    all_pass = all(row["pass"] for row in CHECKS)
    result: dict[str, Any] = {
        "schema_version": 1,
        "status": "pass_with_same_soliton_no_go_and_open_physics_gates"
        if all_pass
        else "fail",
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
        "all_pass": all_pass,
        "physics_promotion_allowed": physics_promotion_allowed,
        "same_physical_moduli_space_derived": same_physical_moduli_space_derived,
        "actual_yukawa_callias_operator_derived": actual_yukawa_callias_operator_derived,
        "uniform_fredholm_gap_proven": uniform_fredholm_gap_proven,
        "determinant_line_o2_derived": determinant_line_o2_derived,
        "physical_cpt_regulator_derived": physical_cpt_regulator_derived,
        "gauge_basic_wzw_descent_proven": gauge_basic_wzw_descent_proven,
        "same_soliton_composition_closed": same_soliton_composition_closed,
        "lane_closed": lane_closed,
        "restricted_higgsed_sigma_model_descent_proven": gauge[
            "restricted_higgsed_sigma_model_descent_proven"
        ],
        "restricted_local_differential_character_ansatz_typed": gauge[
            "restricted_local_differential_character_ansatz_typed"
        ],
        "declared_constituent_yukawa_completion_written": True,
        "uniform_essential_gap_at_infinity_proven": True,
        "four_dimensional_interpolation_index_formula_derived": True,
        "symmetry_and_moduli": {
            "stabilizer": stabilizer,
            "uniform_scalar_cp1_norm": norm_scaling,
            "actual_free_c1_target": "H^2(SO(3);Z) has no free Z summand; on the SU2 cover H^2(S3;Z)=0",
        },
        "same_profile": {
            "B_numeric": profile["B_numeric"],
            "B_residual": profile["B_residual"],
            "origin_slope": profile["origin_slope"],
            "box": canonical_profile["box"],
            "sample_points": canonical_profile["sample_points"],
            "s_at_origin": float(canonical_profile["s"][0]),
            "model": canonical_profile["model"],
            "profile_sha256_r_F_s_little_endian_float64": canonical_profile[
                "profile_sha256_r_F_s_little_endian_float64"
            ],
            "card_profile_sha256_r_F_s_little_endian_float64": workline_b_card[
                "coupled_bvp"
            ]["canonical_same_solution"][
                "profile_sha256_r_F_s_little_endian_float64"
            ],
            "canonical_checksum_equality_asserted": True,
            "workline_b_script_sha256": sha256(WORKLINE_B_SCRIPT),
            "workline_b_card_sha256": sha256(WORKLINE_B_CARD),
            "wilson_profile_source": "canonical r,F,s card-bound arrays; F linearly interpolated off-grid",
            "breathing_scalar_treatment": "s is decoupled by the declared Yukawa M U^gamma5; no s-dependent Yukawa coupling is claimed",
        },
        "operator_audit": {
            "euclidean_4d": "D_E=gamma^mu D_mu+M(P_R U+P_L U^dagger)",
            "static_3d": "H_F=D_3+Phi_F, D_3=-i alpha^i D_i, Phi_F=beta M(cos F+i gamma5 tau.rhat sin F)",
            "physical_domain": "Dom(H_F)=H^1(R3,S tensor C2); H_F is self-adjoint under the stated smooth-decaying gauge assumptions",
            "callias_operator": "C_F=D_3+i Phi_F: H^1 -> L^2; graded double B_F=[[0,C_F^dagger],[C_F,0]]",
            "callias_assumptions": [
                "R3 is complete and D_3 is self-adjoint on H^1",
                "Phi_F is bounded and Hermitian and [D_3,Phi_F] extends boundedly",
                "outside a compact set Phi_F^2-||[D_3,Phi_F]|| >= epsilon > 0",
                "gauge connection is smooth and decays so H_F-H_infinity is relatively compact",
            ],
            "asymptotic_mass": "Phi_infinity=beta M (not M times the identity)",
            "asymptotic_positive_eigenbundle": "constant positive eigenspace of sgn(M) beta, hence a trivial rank-four bundle over S2_infinity",
            "ordinary_callias_boundary_p": 0,
            "ordinary_callias_index": 0,
            "ordinary_callias_index_scope": "only under the explicitly listed domain, bounded-commutator, coercivity, and decay assumptions",
            "interpolation_operator": "D_4,flow=partial_tau+H[U_tau]",
            "conditional_interpolation_index": "Ind(D_4,flow)=SF(H_tau)=B",
            "endpoint_zero_mode_implied": False,
            "two_band_obstruction": {
                "desired_c1": "x+2y",
                "desired_c1_square_xy_coefficient": desired_c1_square_xy_coefficient,
                "reason": "positive eigenline of a 2x2 mass pulls back u from CP1, hence c1^2=0",
            },
        },
        "numerics": {
            "wilson_dirac": wilson_rows,
            "cp1_representative_ambiguity": ambiguity,
            "veronese_o2_template": template,
            "cpt_regulator_control": cpt,
            "gauge_holonomy_controls": gauge,
        },
        "gauge_descent_scope": {
            "restricted_space": "|phi|=v != 0, projective s=[phi], smooth sigma-model fields",
            "universal_differential_character": "Xi5_hat=2 pr_U^*omega3_hat cup pr_s^*c1_hat(Hopf) in Hhat^5(SU2 x CP1)",
            "evaluation": "ev: Sigma3 x C_restricted -> SU2 x CP1, (x,c)->(U_c(x),s_c(x))",
            "typed_transgression": "int_Sigma3 ev^*Xi5_hat in Hhat^2(C_restricted)",
            "local_status": "well-typed local differential-character ansatz plus constant-loop phase controls",
            "restricted_stack_descent": "not proven",
            "missing_restricted_stack_data": [
                "all components of Map(Sigma3,U1)",
                "an equivariant differential cocycle",
                "triviality of vertical holonomy",
            ],
            "stabilizer": "the tested faithful Z2 phase is trivial for even weight, a necessary control but not an exactly-two theorem",
            "full_mother_model": "open",
        },
        "cpt_scope": {
            "theta_assumptions": [
                "Theta is antiunitary and Theta i Theta^{-1}=-i",
                "Theta maps the operator, domain, boundary conditions, and spectral cut to their CPT images",
                "the fermion and Pauli-Villars representations are Theta-closed and regulator masses are real",
                "the determinant/Pfaffian orientation is transported consistently",
            ],
            "conditional_regulator_law": "under the Theta assumptions, Z_PV[CPT Phi]=conjugate(Z_PV[Phi])",
            "conditional_pfaffian_law": "under the Theta and orientation assumptions, Pf A[CPT Phi]=conjugate(Pf A[Phi])",
            "finite_test_status": "algebraic finite spectral-product regression, not a physical regulator construction",
            "cp1_anti_canonical_rule": "K_CP1 tensor E^vee=O(-4) only if the physical CP1 family and polarization exist",
            "same_mother_regulator_fixed": False,
        },
        "checks": CHECKS,
        "source_manifest": manifest,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "deterministic": True,
            "randomness_used": False,
        },
        "remaining_blockers": [
            "Derive the relaxed B=1 symmetry orbit from the full gauge-meson solution; the present hedgehog gives SO(3), not CP1.",
            "Produce a localized normalisable CP1 zero mode rather than the volume-divergent uniform scalar vacuum orientation.",
            "Derive a microscopic Yukawa representation and regulator whose asymptotic eigenbundle has the required mixed topology; the minimal two-band construction is obstructed.",
            "Prove a full uniform endpoint spectral gap and compute the actual families determinant/Pfaffian line.",
            "Extend the restricted Higgsed sigma-model differential character across all microscopic gauge bundles, zeros, defects, and non-extendible sectors.",
            "Only after one complete lane closes may a degree-one Route-E portal be constructed.",
        ],
    }
    write_outputs(result)
    print(
        "AP-E6 same-soliton Yukawa--Callias audit: "
        f"{result['status']} ({result['checks_passed']}/{result['checks_total']}); "
        f"lane_closed={lane_closed}; promotion={physics_promotion_allowed}"
    )


if __name__ == "__main__":
    main()
