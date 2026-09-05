#!/usr/bin/env python3
"""Fixed-VEV Majorana tadpole and four-copy matching interface.

All Yukawa matrices here are deterministic synthetic complex matrices, not a
Spin(10) flavor fit. The bosonic matrix is the existing same-action result;
combining it with synthetic fermions tests assembly only, not a physical point.
The local variable f denotes the Majorana-normalized f_M, not the
unit-bidoublet coupling f_D; the companion Clifford audit gives f_M=2sqrt(6)f_D.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"
AUDIT = RF / "code/verify_p54_theory_audit.py"
BOSON = RF / "output/p54_full_doublet_cw.json"
PS = RF / "output/p54_ps_active_census.json"
SPINOR = RF / "output/p54_spinor_intertwiners.json"
JSON_OUT = RF / "output/p54_fermion_tadpole.json"
MD_OUT = RF / "output/p54_fermion_tadpole.md"


def helpers():
    spec = importlib.util.spec_from_file_location("p54_audit_fermion_helpers", AUDIT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def herm(a):
    return (a + a.conj().T) / 2


def realification(a):
    return np.block([[a.real, -a.imag], [a.imag, a.real]])


def decode(a):
    return np.array(a["real"]) + 1j * np.array(a["imag"])


def cj(a):
    return {"real": np.asarray(a).real.tolist(), "imag": np.asarray(a).imag.tolist()}


def relative(a, b):
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)) / max(np.linalg.norm(b), 1e-15))


def fixed_vev_fermion(yukawas, f, sigma, mu, p126, audit):
    """Here f is f_M. Hold mu fixed in every field derivative."""
    if sigma <= 0 or mu <= 0 or not np.allclose(f, f.T, atol=1e-13):
        raise ValueError("positive sigma, mu and complex symmetric f required")
    mr = sigma * f
    lam = np.linalg.eigvalsh(mr.conj().T @ mr)
    if lam.min() <= 0:
        raise ValueError("invertible hard Majorana matrix required")
    tadpole = -float(np.sum(lam**2 * (np.log(lam / mu**2) - 1))) / (8 * math.pi**2 * sigma)
    dnu2 = tadpole / sigma
    pi = audit.complex_cw(yukawas, mr, mu)
    ct = -dnu2 * p126 / 2
    return {"mr": mr, "tadpole": tadpole, "delta_nu2": dnu2, "pi": pi, "ct": ct, "fixed": herm(pi + ct)}


def schur_retune(m0, correction, p10):
    tree_values, tree_vectors = np.linalg.eigh(herm(m0))
    light = tree_vectors[:, 0]
    heavy = tree_vectors[:, 1:]
    pre = herm(m0 + correction)
    low = lambda dx: float(np.linalg.eigvalsh(pre - dx * p10)[0])
    if not (low(-1) > 0 and low(1) < 0):
        return {"success": False, "reason": "synthetic assembly does not bracket a stable first zero"}
    dx = brentq(low, -1, 1, xtol=1e-14)
    tuned = herm(pre - dx * p10)
    values, vectors = np.linalg.eigh(tuned)
    c = vectors[:, 0]
    cblock = herm(heavy.conj().T @ tuned @ heavy)
    b = heavy.conj().T @ tuned @ light
    a = float(np.vdot(light, tuned @ light).real)
    rhs = float(np.vdot(b, np.linalg.solve(cblock, b)).real)
    return {
        "success": True,
        "delta_xi02": dx,
        "matrix": cj(tuned),
        "eigenvalues": values.tolist(),
        "light_coefficients": cj(c),
        "nearest_heavy_gap": float(values[1]),
        "eigenpair_residual": float(np.linalg.norm(tuned @ c)),
        "Schur_A": a,
        "Schur_bCinvb": rhs,
        "Schur_residual": abs(a - rhs),
        "tree_light_rotation_radians": math.acos(min(1.0, abs(np.vdot(light, c)))),
        "tree_heavy_values": tree_values[1:].tolist(),
    }


def run():
    audit = helpers()
    boson = json.loads(BOSON.read_text())
    ps = json.loads(PS.read_text())
    sigma = float(boson["vacuum"]["sigma"])
    mu = float(boson["scheme"]["mu_over_omega"])
    p126 = np.diag([0.0, 0.0, 1.0, 1.0]).astype(complex)
    p10 = np.eye(4) - p126
    rng = np.random.default_rng(20260906)

    def unitary(n):
        u, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        return u

    def symmetric(scale):
        z = scale * (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
        return (z + z.T) / 2

    cases = []
    first = None
    for index, masses in enumerate(([0.3, 0.6, 0.9], [0.45, 0.45, 0.8], [0.12, 0.4, 1.0])):
        u = unitary(3)
        f = u.conj() @ np.diag(masses) @ u.conj().T
        ys = [symmetric(scale) for scale in (0.07, 0.09, 0.08, 0.05)]
        result = fixed_vev_fermion(ys, f, sigma, mu, p126, audit)
        radial_step = sigma * 2e-5
        zero = np.zeros(8)
        radial = lambda s: audit.potential(zero, ys, s * f, mu)
        fd_t = (radial(sigma + radial_step) - radial(sigma - radial_step)) / (2 * radial_step)
        # CT is fixed at the matching background while the field varies.
        fixed_radial = lambda s: radial(s) - result["delta_nu2"] * s**2 / 2
        fd_fixed = (fixed_radial(sigma + radial_step) - fixed_radial(sigma - radial_step)) / (2 * radial_step)
        corrected_potential = lambda q: (
            audit.potential(q, ys, result["mr"], mu)
            - result["delta_nu2"] * float(q @ realification(p126) @ q) / 4
        )
        analytic_real = realification(result["fixed"])
        fd_h = audit.numerical_hessian(corrected_potential, 8, 1e-4)
        fd_h2 = audit.numerical_hessian(corrected_potential, 8, 3e-4)
        # Independent active/sterile basis transformations are permitted for
        # MD. Use a common family congruence here to preserve symmetric Ys.
        uf = unitary(3)
        f_family = uf.T @ f @ uf
        ys_family = [uf.T @ y @ uf for y in ys]
        family = fixed_vev_fermion(ys_family, f_family, sigma, mu, p126, audit)
        # If z=V z_new, Ys_new[b]=sum_a Ys[a] V[a,b]. Transform P126 as well.
        vc = unitary(4)
        ys_copy = [sum(ys[a] * vc[a, b] for a in range(4)) for b in range(4)]
        p_copy = vc.conj().T @ p126 @ vc
        copy = fixed_vev_fermion(ys_copy, f, sigma, mu, p_copy, audit)
        qnew = rng.normal(size=8) * 0.003
        qold = realification(vc) @ qnew
        potential_copy_error = abs(
            audit.potential(qold, ys, result["mr"], mu)
            - audit.potential(qnew, ys_copy, result["mr"], mu)
        )
        row = {
            "case": index,
            "f_Takagi_values": list(masses),
            "radial_tadpole": result["tadpole"],
            "radial_tadpole_finite_difference": fd_t,
            "radial_relative_error": relative(fd_t, result["tadpole"]),
            "finite_CT_radial_relative_residual": abs(fd_fixed) / max(abs(result["tadpole"]), 1e-15),
            "full_8_real_Hessian_relative_error_step_1e4": relative(fd_h, analytic_real),
            "full_8_real_Hessian_relative_error_step_3e4": relative(fd_h2, analytic_real),
            "family_tadpole_relative_error": relative(family["tadpole"], result["tadpole"]),
            "family_fixed_matrix_relative_error": relative(family["fixed"], result["fixed"]),
            "copy_fixed_matrix_relative_error": relative(copy["fixed"], vc.conj().T @ result["fixed"] @ vc),
            "copy_potential_absolute_error": potential_copy_error,
            "projector_covariance_error": relative(copy["ct"], vc.conj().T @ result["ct"] @ vc),
        }
        cases.append(row)
        if first is None:
            first = (ys, f, result)
    ys, f, result = first
    m0 = decode(boson["tree_doublet_matrix_over_omega2"])
    boson_fixed = decode(boson["bosonic_CW"]["total_matrix_over_omega2"]) + decode(boson["tadpoles"]["doublet_CT_over_omega2"])
    assembled = schur_retune(m0, boson_fixed + result["fixed"], p10)
    v_assembly = unitary(4)
    rotated = schur_retune(
        v_assembly.conj().T @ m0 @ v_assembly,
        v_assembly.conj().T @ (boson_fixed + result["fixed"]) @ v_assembly,
        v_assembly.conj().T @ p10 @ v_assembly,
    )
    checks = [
        {"name": "fixed-mu radial derivative on three nontrivial complex families", "pass": max(c["radial_relative_error"] for c in cases) < 2e-8},
        {"name": "fixed finite mass CT cancels the radial tadpole", "pass": max(c["finite_CT_radial_relative_residual"] for c in cases) < 2e-8},
        {"name": "all eight real doublet directions agree with exact complex fixed-VEV matrix", "pass": max(c["full_8_real_Hessian_relative_error_step_1e4"] for c in cases) < 3e-5},
        {"name": "second finite-difference step confirms eight-direction Hessian", "pass": max(c["full_8_real_Hessian_relative_error_step_3e4"] for c in cases) < 3e-5},
        {"name": "family-unitary covariance of tadpole and fixed matrix", "pass": max(max(c["family_tadpole_relative_error"], c["family_fixed_matrix_relative_error"]) for c in cases) < 1e-11},
        {"name": "copy-unitary covariance includes P126 counterterm", "pass": max(max(c["copy_fixed_matrix_relative_error"], c["projector_covariance_error"]) for c in cases) < 1e-11},
        {"name": "finite-field potential respects copy basis change", "pass": max(c["copy_potential_absolute_error"] for c in cases) < 1e-14},
        {"name": "degenerate heavy Takagi masses need no eigenvector derivative", "pass": cases[1]["full_8_real_Hessian_relative_error_step_1e4"] < 3e-5},
        {"name": "synthetic fermion plus actual bosonic matrix has a stable Schur-retuned zero", "pass": assembled.get("success", False) and assembled.get("nearest_heavy_gap", -1) > 0 and assembled.get("Schur_residual", 1) < 1e-10 and assembled.get("eigenpair_residual", 1) < 1e-10},
        {"name": "assembled Schur retuning is copy-basis covariant", "pass": rotated.get("success", False) and abs(rotated.get("delta_xi02", 99) - assembled.get("delta_xi02", 0)) < 1e-11},
        {"name": "PS logarithmic repair audit is available and passed", "pass": ps["summary"]["all_pass"]},
    ]
    for c in checks:
        c["pass"] = bool(c["pass"])
    return {
        "schema": "route-f-p54-fermion-fixed-vev-interface-v1",
        "date": "2026-09-05",
        "scope": "deterministic synthetic Yukawa interface tests; not a Spin(10) fit, not actual Clebsch matching",
        "family_matrix_convention": "f in the synthetic code is f_M; MR=sigma*f_M and f_M=2*sqrt(6)*f_D for unit canonical 126 doublet coupling f_D",
        "scheme": {"mu_over_omega": mu, "mu_held_fixed_under_field_derivatives": True, "sigma_over_omega": sigma, "Majorana_mode_count": 3, "complex_copy_count": 4},
        "formulae": {
            "X": "MR^dagger MR, MR=sigma*f_M (local variable f)",
            "t_sigma_F": "-Tr[X^2(log(X/mu^2)-1)]/(8*pi^2*sigma)",
            "delta_nu2_F": "t_sigma_F/sigma",
            "doublet_fixed_vev_delta": "Pi_F - delta_nu2_F*P126/2",
            "radial_CT_potential": "-delta_nu2_F*sigma_field^2/2 with delta_nu2_F fixed at matching background",
        },
        "synthetic_example": {"f": cj(f), "Y_copies": [cj(y) for y in ys], "P126": cj(p126), "MR": cj(result["mr"]), "t_sigma_F": result["tadpole"], "delta_nu2_F": result["delta_nu2"], "Pi_F": cj(result["pi"]), "doublet_CT_F": cj(result["ct"]), "fixed_vev_delta_F": cj(result["fixed"])},
        "independent_regressions": cases,
        "synthetic_boson_fermion_assembly": {"scope": "actual historical-background bosonic matrix plus arbitrary synthetic Yukawa vertices in its declared four-copy basis; assembly regression only", "retuning": assembled, "unitary_rotated_retuning": rotated},
        "fit_readiness": {
            "fermion_fixed_vev_interface": "implemented_and_synthetic_tested",
            "historical_PS_census_and_MI_matching": "failed",
            "four_parent_PS_repair_logarithmic_identities": "passed_all_six_gauge_directions",
            "finite_PS_threshold_matching": "open",
            "P54_complete_PS_Yukawa_beta_system": "open",
            "spinor_absolute_Clebsch_normalization": "verified_in_companion_20_check_audit; f_M=2sqrt(6)f_D",
            "spinor_four_copy_relative_phases": "open",
            "actual_quantum_eigenvector_flavor_fit": "open",
            "type_I_plus_type_II_threshold_matching": "open",
            "physical_global_optimizer_run": False,
            "physical_full_loop_point_certified": False,
        },
        "checks": checks,
        "summary": {"passed": sum(c["pass"] for c in checks), "total": len(checks), "all_pass": all(c["pass"] for c in checks)},
        "sources": [{"path": str(p.relative_to(ROOT)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in (AUDIT, BOSON, PS, SPINOR, Path(__file__).resolve())],
    }


def markdown(r):
    x = r["synthetic_example"]
    s = r["synthetic_boson_fermion_assembly"]["retuning"]
    return "\n".join([
        "# Fixed-VEV fermionic tadpole and doublet interface", "",
        "Status: executable synthetic-matrix regression. No flavor fit or actual Spin(10) Clebsch assignment is claimed.", "",
        r"For $M_R=\sigma f_M$, $X=M_R^\dagger M_R$, at fixed renormalization scale $\mu$:", "",
        r"The local code variable `f` denotes $f_M$, not the unit-bidoublet $f_D$. The independently verified canonical translation is $f_M=2\sqrt6 f_D$. The arbitrary synthetic $Y_a$ are not assigned actual Spin(10) Clebsch tensors.", "",
        r"$$V_F(0)=-\frac{\mathrm{Tr}[X^2(\log(X/\mu^2)-3/2)]}{32\pi^2},\quad \frac{dX}{d\sigma}=\frac{2X}{\sigma}.$$", "",
        r"Differentiating the spectral trace gives $t_{\sigma,F}=-\mathrm{Tr}[X^2(\log(X/\mu^2)-1)]/(8\pi^2\sigma)$.", "",
        r"The invariant finite counterterm has radial potential $-\delta\nu_F^2\sigma_{\rm field}^2/2$, hence $\delta\nu_F^2=t_{\sigma,F}/\sigma$. Its canonical doublet Hessian is $-\delta\nu_F^2P_{126}/2$.", "",
        r"$$\Delta M_{D,F}^{2,\mathrm{fixed\ VEV}}=\Pi_F-\frac{\delta\nu_F^2}{2}P_{126}.$$", "",
        "The same counterterm value is held fixed during the finite-difference field variation. The scale is not varied with the radial field.", "",
        f"Synthetic example: t_sigma,F={x['t_sigma_F']:.12g}, delta_nu2,F={x['delta_nu2_F']:.12g} in omega=1 units.", "",
        "Three nontrivial complex 3x3 families, including a degenerate Takagi pair, test the radial derivative, finite-counterterm cancellation, all eight real doublet directions, family unitaries, and copy unitaries with P126 rotated consistently.", "",
        f"Actual bosonic matrix plus synthetic fermion assembly: delta_xi02={s.get('delta_xi02', float('nan')):.12g}; heavy gap={s.get('nearest_heavy_gap', float('nan')):.12g}; Schur residual={s.get('Schur_residual', float('nan')):.3e}. This is an assembly regression, not a physical fit point.", "",
        "Physical fit gates: historical PS census/matching failed; repaired logarithmic identities pass; absolute spinor normalization is independently checked. Finite PS thresholds, complete P54 PS Yukawa betas, common four-copy phases, and type-I plus type-II matching remain open.", "",
        f"Checks: **{r['summary']['passed']}/{r['summary']['total']}**. Full matrices and residuals are in `p54_fermion_tadpole.json`.", "",
    ])


if __name__ == "__main__":
    report = run()
    JSON_OUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    MD_OUT.write_text(markdown(report))
    print(f"Fermion fixed-VEV interface: {report['summary']['passed']}/{report['summary']['total']} checks")
    for c in report["checks"]:
        if not c["pass"]:
            print("FAIL:", c["name"])
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)
