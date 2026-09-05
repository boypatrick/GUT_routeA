#!/usr/bin/env python3
"""Same-action tree light quartic after relaxing every massive scalar.

This tests a local stability condition missed by a semidefinite Hessian.
No parameter or background is fitted in this calculation.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def main():
    p1path = RF / "code/verify_p54_p1_hessian_spectrum.py"
    searchpath = RF / "code/search_p54_hierarchical_p1.py"
    p2path = RF / "output/p54_p2_two_site_matching.json"
    p1 = module("p54_quartic_p1", p1path)
    search = module("p54_quartic_search", searchpath)
    p2 = json.loads(p2path.read_text())
    hierarchy = json.loads((RF / "output/p54_hierarchical_p1_search.json").read_text())
    sigma, vs = float(p2["hierarchical_vacuum_ratio"]), float(hierarchy["vevs"]["vs"])
    params = search.parameters_at(p1, sigma, vs)
    params["xi02"] = float(p2["doublet_tuning"]["xi02"])
    x0 = p1.vacuum_vector(1., sigma, vs)
    potential = p1.potential_factory(params)
    vf = p1.jax.jit(potential)
    hf = p1.jax.jit(p1.hessian(potential))

    def hess(x):
        h = np.asarray(hf(p1.anp.asarray(x)), dtype=float)
        return (h + h.T) / 2

    h0 = hess(x0)
    lam, u = np.linalg.eigh(h0)
    tol = float(p2["doublet_tuning"]["zero_tolerance"])
    heavy = lam > tol
    inv = (u[:, heavy] / lam[heavy]) @ u[:, heavy].T
    zero = u[:, np.abs(lam) < tol]
    sym, _, _ = search.symmetry_complement(p1, x0)
    symproj = sym @ np.linalg.pinv(sym, rcond=1e-11)
    light, sv, _ = np.linalg.svd(zero - symproj @ zero, full_matrices=False)
    light = light[:, sv > 2e-8]
    v0 = float(vf(p1.anp.asarray(x0)))
    rows = []
    for k in (0, 1):
        q = light[:, k]
        jets = []
        for eps in (.02, .01):
            hp, hm = hess(x0 + eps*q), hess(x0 - eps*q)
            j = ((hp-hm)/(2*eps)) @ q
            fourth = float(q @ ((hp+hm-2*h0)/eps**2) @ q)
            direct = fourth / 6
            exchange = float(j @ inv @ j) / 2
            jets.append({"epsilon": eps, "direct_lambda": direct,
                         "heavy_exchange_subtraction": exchange,
                         "relaxed_lambda": direct-exchange,
                         "source_zero_space_norm": float(np.linalg.norm(zero.T @ j))})
        z2 = -0.5 * inv @ j
        series = []
        for t in (.04, .02, .01):
            vp = float(vf(p1.anp.asarray(x0+t*q+t*t*z2)))
            vm = float(vf(p1.anp.asarray(x0-t*q+t*t*z2)))
            # Subtract the tiny residual quadratic eigenvalue before dividing.
            quartic = 4*((vp+vm)/2-v0-0.5*t*t*float(q@h0@q))/t**4
            series.append({"t": t, "relaxed_path_quartic_estimate": quartic})
        rows.append({"light_direction": k, "jets": jets, "direct_path": series})
    target = rows[0]["jets"][-1]["relaxed_lambda"]
    checks = [
        {"name": "four real light directions", "pass": light.shape[1] == 4},
        {"name": "290 positive scalar modes", "pass": int(heavy.sum()) == 290},
        {"name": "no negative quadratic modes", "pass": float(lam.min()) > -tol},
        {"name": "source orthogonal to all zero modes", "pass": max(j["source_zero_space_norm"] for r in rows for j in r["jets"]) < 1e-8},
        {"name": "quartic stable under derivative step change", "pass": max(abs(r["jets"][0]["relaxed_lambda"]-r["jets"][1]["relaxed_lambda"]) for r in rows) < 1e-7},
        {"name": "SM-equivalent directions agree", "pass": abs(target-rows[1]["jets"][-1]["relaxed_lambda"]) < 1e-7},
        {"name": "direct relaxed-path check", "pass": max(abs(r["direct_path"][-1]["relaxed_path_quartic_estimate"]-target) for r in rows) < 1e-4},
    ]
    passed = sum(c["pass"] for c in checks)
    result = {"date": "2026-09-05", "model": "P54PQ-v2", "scope": "tree background; full massive-scalar relaxation at quartic order",
              "formula": "lambda_EFT=V4[q,q,q,q]/6-J^T H_heavy^{-1} J/2; J=V3[.,q,q]",
              "sigma_over_omega": sigma, "vs_over_omega": vs, "xi02": params["xi02"],
              "relaxed_lambda": target, "local_light_quartic_positive": target > 0,
              "global_stability_proved": False, "loop_stationarity_proved": False,
              "directions": rows, "checks": checks,
              "summary": {"passed": passed, "total": len(checks), "all_pass": passed == len(checks)},
              "sources": [{"path": str(p.relative_to(ROOT)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in (p1path, searchpath, p2path, RF/"output/p54_hierarchical_p1_search.json", Path(__file__))]}
    out = RF / "output/p54_relaxed_light_quartic.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(f"Same-action relaxed quartic: lambda={target:.12g}; checks {passed}/{len(checks)}")
    print(f"Positive local quartic: {target > 0}; global/loop stability remains unproved")
    if passed != len(checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
