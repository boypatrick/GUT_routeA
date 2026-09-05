#!/usr/bin/env python3
"""Finite complex Yukawa matching algebra, not uncomputed PS diagrams.

Uses the actual P54 copy dictionary with synthetic family matrices/finite
self-energies to test the canonical matching operation.  No missing vertex,
kinetic or EFT term is silently set to zero in the physical-fit interface.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

RF = Path(__file__).resolve().parents[1]
PHASE = RF / "output/p54_common_yukawa_phase.json"
BOSON = RF / "output/p54_full_doublet_cw.json"
OUT = RF / "output/p54_finite_yukawa_interface"
LOOP = 16 * math.pi**2


def decode(x):
    return np.asarray(x["real"]) + 1j * np.asarray(x["imag"])


def cjson(x):
    x = np.asarray(x)
    return {"real": x.real.tolist(), "imag": x.imag.tolist()}


def herm(x):
    return (x + x.conj().T) / 2


def hermitian_power(x, power):
    if not np.allclose(x, x.conj().T, atol=1e-12, rtol=1e-12):
        raise ValueError("kinetic matrix must be Hermitian")
    w, u = np.linalg.eigh(x)
    if w.min() <= 0:
        raise ValueError("canonical matching requires positive kinetic metric")
    return (u * w**power) @ u.conj().T


def loop_f_h(x, y, mu):
    """Finite MSbar kernels from nonsingular Feynman-parameter integrals.

    x,y are squared masses. Both zero is deliberately rejected: a hard
    matching kernel must not absorb a purely light IR diagram.
    """
    if x < 0 or y < 0 or mu <= 0 or x + y == 0:
        raise ValueError("nonnegative squared masses, one hard mass and mu>0 required")
    logmass = lambda t: math.log(((1-t)*x + t*y) / mu**2)
    f = -quad(logmass, 0, 1, epsabs=1e-12, epsrel=1e-12)[0] / LOOP
    h = quad(lambda t: (1-t)*logmass(t), 0, 1,
             epsabs=1e-12, epsrel=1e-12)[0] / LOOP
    return np.array([f, h])


def finite_match(vertices, k_left, k_right, k_scalar, curvature=None):
    """Weyl convention psi_L^T Y_a psi_R z_a, with all inputs supplied.

    Z=I+K. Exact algebraic square roots provide a covariance regression;
    finite one-loop inputs do not make their higher powers a two-loop result.
    """
    if any(x is None for x in (vertices, k_left, k_right, k_scalar)):
        raise ValueError("finite vertices and all three hard kinetic corrections are required")
    zl = np.eye(k_left.shape[0]) + k_left
    zr = np.eye(k_right.shape[0]) + k_right
    zh = np.eye(k_scalar.shape[0]) + k_scalar
    sl, sr, sh = [hermitian_power(z, -.5) for z in (zl, zr, zh)]
    yc = np.einsum("ij,ajk,kl,ab->bil", sl.T, vertices, sr, sh, optimize=True)
    dc = None if curvature is None else herm(sh.conj().T @ curvature @ sh)
    return {"vertices": yc, "curvature": dc, "scalar_map": sh,
            "left_map": sl, "right_map": sr, "scalar_metric": zh}


def first_order_match(tree, delta_vertex, kl, kr, kh):
    leg = np.einsum("ij,ajk->aik", kl.T, tree)
    leg += np.einsum("aij,jk->aik", tree, kr)
    leg += np.einsum("aij,ab->bij", tree, kh)
    return tree + delta_vertex - leg / 2


def require_physical_inputs(payload):
    """Fail closed rather than promote a synthetic assembly to a fit."""
    required = ("upper_finite_gauge", "lower_finite_gauge",
                "upper_finite_Yukawa", "lower_finite_Yukawa",
                "all_active_PS_flow", "lower_EFT_flow",
                "sequential_Weinberg_matching", "same_action_scalar_feedback")
    missing = [name for name in required if payload.get(name) is not True]
    if missing:
        raise ValueError("physical fit map incomplete: " + ", ".join(missing))


def run():
    phase = json.loads(PHASE.read_text())
    boson = json.loads(BOSON.read_text())
    rng = np.random.default_rng(20260905)
    unitary = lambda n: np.linalg.qr(rng.normal(size=(n, n)) +
                                    1j*rng.normal(size=(n, n)))[0]
    def random_matrix(n, scale):
        return scale*(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))
    h = random_matrix(3, .12); h = (h+h.T)/2
    f = random_matrix(3, .025); f = (f+f.T)/2
    coefficients = phase["raw_copy_coefficients"]
    tree = {}
    for species, rows in coefficients.items():
        tree[species] = (decode(rows["h_raw"])[:, None, None]*h +
                         decode(rows["f_raw"])[:, None, None]*f)
    kl, kr, kh = [herm(random_matrix(n, .007)) for n in (3, 3, 4)]
    dy = np.array([random_matrix(3, .001) for _ in range(4)])
    # These finite arrays are explicitly synthetic; actual P54 loop state
    # sums must replace them before any physical-fit call is allowed.
    y = tree["u"]
    c = decode(phase["bosonic_light_coefficients"])
    d0 = decode(boson["tree_doublet_matrix_over_omega2"])
    pi = decode(boson["bosonic_CW"]["total_matrix_over_omega2"])
    ct = decode(boson["tadpoles"]["doublet_CT_over_omega2"])
    dx = boson["retuned_bosonic_eigenpair"]["delta_xi02"]
    d = herm(d0+pi+ct-dx*np.diag([1, 1, 0, 0]))
    matched = finite_match(y+dy, kl, kr, kh, d)
    zh = matched["scalar_metric"]
    norm = math.sqrt(np.vdot(c, zh@c).real)
    cc = hermitian_power(zh, .5)@c/norm
    physical_direction = matched["scalar_map"]@cc
    physical_y = np.einsum("a,aij->ij", cc, matched["vertices"])
    direct_y = matched["left_map"].T @ np.einsum("a,aij->ij", c/norm, y+dy) @ matched["right_map"]

    ul, ur, v = unitary(3), unitary(3), unitary(4)
    yrot = np.einsum("ij,ajk,kl,ab->bil", ul.T, y+dy, ur, v)
    rot = finite_match(yrot, ul.conj().T@kl@ul, ur.conj().T@kr@ur,
                       v.conj().T@kh@v, v.conj().T@d@v)
    predicted_yrot = np.einsum("ij,ajk,kl,ab->bil", ul.T,
                              matched["vertices"], ur, v)
    # Down/e multiply z*, so their copy metric is Z_H*, not Z_H.
    down = finite_match(tree["d"], kl, kr, kh.conj())
    down_rot_vertices = np.einsum("aij,ab->bij", tree["d"], v.conj())
    down_rot = finite_match(down_rot_vertices, kl, kr,
                            (v.conj().T@kh@v).conj())
    down_predicted = np.einsum("aij,ab->bij", down["vertices"], v.conj())
    down_contraction = np.einsum("a,aij->ij", cc.conj(), down["vertices"])
    down_direct = down["left_map"].T @ np.einsum("a,aij->ij", c.conj()/norm, tree["d"]) @ down["right_map"]

    expansion_errors = []
    for epsilon in (1., .5, .25):
        exact = finite_match(y+epsilon*dy, epsilon*kl, epsilon*kr,
                             epsilon*kh)["vertices"]
        linear = first_order_match(y, epsilon*dy, epsilon*kl,
                                   epsilon*kr, epsilon*kh)
        expansion_errors.append(float(np.linalg.norm(exact-linear)))
    mu = .73
    kernels = []
    kernel_error, slope_error = 0., 0.
    for x, yy in ((1., 0.), (0., 1.), (1., 1.), (1., 1+1e-12),
                  (3., .07), (.07, 3.)):
        got = loop_f_h(x, yy, mu)
        if yy == 0:
            expected = np.array([1-math.log(x/mu**2), .5*math.log(x/mu**2)-.25])/LOOP
        elif x == 0:
            expected = np.array([1-math.log(yy/mu**2), .5*math.log(yy/mu**2)-.75])/LOOP
        elif abs(yy/x-1) < 1e-10:
            expected = np.array([-math.log(x/mu**2), .5*math.log(x/mu**2)])/LOOP
        else:
            r = yy/x
            expected = np.array([1-(x*math.log(x/mu**2)-yy*math.log(yy/mu**2))/(x-yy),
                .5*math.log(x/mu**2)+(.5*r*r*math.log(r)-.75*r*r+r-.25)/(1-r)**2])/LOOP
        kernel_error = max(kernel_error, float(np.max(abs(got-expected))))
        derivative = (loop_f_h(x, yy, mu*math.exp(1e-4))-
                      loop_f_h(x, yy, mu*math.exp(-1e-4)))/(2e-4)
        slope_error = max(slope_error, float(np.max(abs(derivative-np.array([2., -1.])/LOOP))))
        kernels.append({"squared_masses": [x, yy], "f_h": got.tolist()})
    rejected = {}
    for name, function in (
        ("pure_light_kernel", lambda: loop_f_h(0, 0, mu)),
        ("negative_kinetic_metric", lambda: finite_match(y, kl, kr, -2*np.eye(4))),
        ("missing_finite_vertex", lambda: finite_match(None, kl, kr, kh)),
        ("incomplete_physical_fit", lambda: require_physical_inputs({"upper_finite_gauge": True})),
    ):
        try:
            function(); rejected[name] = False
        except ValueError:
            rejected[name] = True
    checks = []
    def check(name, error, tolerance=1e-11):
        checks.append({"name": name, "residual": float(error), "tolerance": tolerance,
                       "pass": bool(error < tolerance)})
    check("actual common phase prerequisite", 0 if phase["summary"]["all_pass"] else 1)
    check("actual bosonic null vector", np.linalg.norm(d@c))
    check("finite scalar normalization preserves null direction", np.linalg.norm(matched["curvature"]@cc))
    check("canonical light vector has unit norm", abs(np.vdot(cc, cc)-1))
    check("zero-mode original direction changes only by kinetic normalization", np.linalg.norm(physical_direction-c/norm))
    check("canonical Yukawa and original-coordinate projection agree", np.linalg.norm(physical_y-direct_y))
    check("independent left right and complex copy unitary covariance", np.linalg.norm(rot["vertices"]-predicted_yrot))
    check("same finite scalar map transports curvature", np.linalg.norm(rot["curvature"]-v.conj().T@matched["curvature"]@v))
    check("conjugate-copy down matching covariance", np.linalg.norm(down_rot["vertices"]-down_predicted))
    check("down contraction uses conjugate null vector and conjugate kinetic metric", np.linalg.norm(down_contraction-down_direct))
    check("finite kernel limits and closed forms", kernel_error)
    check("finite kernel matching-scale derivatives", slope_error)
    for name, passed in rejected.items():
        check("reject "+name, 0 if passed else 1)
    ratios = [expansion_errors[i]/expansion_errors[i+1] for i in range(2)]
    check("linear matching agrees through first order", max(abs(r-4) for r in ratios), .25)
    report = {
        "schema": "p54-finite-yukawa-canonical-interface-v1", "date": "2026-09-05",
        "scope": "actual scalar/copy dictionary and curvature; synthetic family matrices and synthetic finite diagrams; not physical PS thresholds or a fit",
        "formulae": {
            "canonical": "Yc[a]=sum_b ZL^(-T/2) (Y[b]+deltaY[b]) ZR^(-1/2) [ZH^(-1/2)]ba",
            "first_order": "Yc[a]=Y[a]+deltaY[a]-(KL^T Y[a]+Y[a] KR+sum_b Y[b] KH[ba])/2",
            "curvature": "Dc=ZH^(-1/2) D ZH^(-1/2)",
            "zero_mode": "Dc cc=0 with cc=ZH^(1/2)c/sqrt(c^dagger ZH c); original physical direction=c/sqrt(c^dagger ZH c)",
            "down": "Down/e copy vertices multiply z*; transform with conjugate ZH and conjugate c, not the up convention",
            "f": "-integral_0^1 log(((1-t)x+t*y)/mu^2) dt/(16*pi^2)",
            "h": "integral_0^1 (1-t) log(((1-t)x+t*y)/mu^2) dt/(16*pi^2)",
        },
        "synthetic_regression": {"h_raw": cjson(h), "f_raw": cjson(f), "K_left": cjson(kl),
            "K_right": cjson(kr), "K_scalar": cjson(kh), "delta_vertex_up": cjson(dy),
            "kinetic_light_norm": norm, "light_canonical": cjson(cc),
            "expansion_errors": expansion_errors, "error_ratios": ratios},
        "loop_kernel_regressions": kernels,
        "missing_physical_terms": ["P54-specific heavy and mixed vertex graphs",
            "all fermion/scalar kinetic state sums, including heavy vectors in the declared gauge",
            "UV minus EFT subtraction and both-site finite Yukawa coefficients",
            "sequential type-I/type-II Wilson matching and all same-action scalar feedback"],
        "physical_fit_performed": False,
        "checks": checks, "summary": {"passed": sum(c["pass"] for c in checks),
            "total": len(checks), "all_pass": all(c["pass"] for c in checks)},
        "references": [{"url": "https://arxiv.org/html/2310.16563v2", "use": "Appendices A/B: canonical finite matching and scalar loop kernels only; no SU(5) group factors imported"}],
        "sources": [{"path": str(p.relative_to(RF.parent)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                    for p in (Path(__file__), PHASE, BOSON)],
    }
    return report


def main():
    report = run()
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n")
    OUT.with_suffix(".md").write_text("\n".join([
        "# Finite complex Yukawa matching interface", "", report["scope"], "",
        f"Checks: {report['summary']['passed']}/{report['summary']['total']}.", "",
        "Canonical normalization must transform vertices and curvature with the same scalar map. At an exact zero, the original-coordinate light direction changes only by its kinetic normalization, not an independent wavefunction-induced rotation.", "",
        "The conjugate Higgs vertices use the conjugate scalar kinetic matrix and light vector. All finite corrections are mandatory arguments; no absent PS threshold is silently set to zero.", "",
        "The f/h kernels include unequal, equal and one-zero masses, with matching-scale derivatives checked. Pure-light integrals are rejected by this hard-kernel API.", "",
        "These tests use synthetic finite diagrams, not calculated P54 finite Yukawa thresholds. The physical fit remains disabled until the required EFT maps exist.", ""]))
    print(json.dumps(report["summary"]))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
