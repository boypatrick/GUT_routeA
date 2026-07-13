#!/usr/bin/env python3
"""DYN-4a: benchmark replay and a connected local zeta profile likelihood.

Provenance.  The benchmark flavor artifacts (the publication closure card
with the fitted Yukawa sectors, and the source-Majorana texture-rank
summary) were pruned from the lean tree; the author restored the archival
output set under route_E/output/.  This audit reads them read-only, records
their sha256, and REGENERATES the paper's reproducibility anchors from the
recipe preserved in HEAD code:

    zeta = 0.1076472949 + 0.0736514853 i ,   contact fraction 1.304275e-1.

Pipeline:
  1. Replay: Y_e, Y_nu (closure card) + the fixed light-nu benchmark
     (hardcoded in HEAD audit_source_majorana_texture_rank.py) -> inverse
     type-I seesaw M_R -> Veronese + contact decomposition (c_hat = -K_tr
     orientation, M_* = sigma_max).  Gates: zeta and the contact fraction
     reproduce the archived card AND the paper digits; the forward seesaw
     closes; the PMNS replay returns the old benchmark targets exactly.
  2. Rank-contract compliance (Audit 1b): the decomposition is reported
     tensor-by-tensor over an aligned 6-tensor basis (5 Veronese + contact)
     with a residual cascade, not a hollow two-tensor pass/fail.
  3. Target refresh: NuFit-6.0 (JHEP 12 (2024) 216, arXiv:2410.05380)
     IC24-with-SK-atm as primary, IC19-without as variant; shift table vs
     the old hardcoded targets (the sin^2 theta_23 octant flip is ~ -6.9
     sigma).
  4. zeta profile likelihood, CONDITIONAL on the benchmark Dirac sector: fit
     M_R(zeta') = M_*(M_V + zeta' c_hat) over the complex zeta' plane,
     chi^2 against the refreshed targets (5 observables), report the joint
     best candidate and connected-local one-parameter profile Delta-chi^2
     intervals;
     delta_CP is reported as a prediction, not fitted.  This is not called a
     normalized Bayesian posterior because no prior or evidence integral is
     defined.
  5. Sensitivity-window refresh: the old loose windows (Delta_s =
     3.208798e-5, Delta_delta = 5.424119e-5 rad) are re-derived with the
     historical tolerance definition as a regression, then recomputed as
     Delta-chi^2 <= 1 windows with the refreshed targets.

Boundary: the profile likelihood is conditional on the benchmark Y_nu, Y_e
(the full covariant refit of all Yukawa sectors is DYN-4b).  The likelihood
contains a needle-like basin: the connected component around the reproducible
best candidate is converged, but a global-minimum or disconnected-component
completeness theorem is not claimed.  This is a fit, not a prediction;
normal ordering is assumed, m1 is fixed at 1e-3 eV, and v_u = 100 GeV follows
the archival convention.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import bisect, minimize

from route_e_paths import ARCHIVE_OUTPUT, AUDIT_OUTPUT, REPO_ROOT

ROOT = REPO_ROOT
ARCH = ARCHIVE_OUTPUT
OUT = AUDIT_OUTPUT / "audit1"

CHECKS: list[tuple[str, bool]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cmat(raw) -> np.ndarray:
    return np.array([[c["re"] + 1j * c["im"] for c in row] for row in raw],
                    dtype=complex)


# --------------------------------------------------------------- section 1
print("== DYN-4a section 1: benchmark replay from archival inputs ==")

card_path = ARCH / "publication_closure_card" / "publication_closure_card.json"
rank_path = ARCH / "source_majorana_texture_rank" / "summary.json"
card = json.loads(card_path.read_text())
arch_rank = json.loads(rank_path.read_text())
yuk = {k: cmat(v) for k, v in card["selected_row"]["Yukawa_fit"].items()}
Y_e, Y_nu = yuk["charged_lepton"], yuk["neutrino_dirac"]

# fixed light-nu benchmark hardcoded in HEAD audit_source_majorana_texture_rank.py
OLD = {"m1_eV": 1.0e-3, "dm21": 7.42e-5, "dm31": 2.517e-3,
       "s12": 0.304, "s13": 0.0222, "s23": 0.573,
       "delta": 1.20 * math.pi, "a21": 0.35 * math.pi, "a31": 1.10 * math.pi}


def standard_pmns(s12, s13, s23, delta, a21, a31):
    sq = math.sqrt
    s12_, s13_, s23_ = sq(s12), sq(s13), sq(s23)
    c12, c13, c23 = sq(1 - s12), sq(1 - s13), sq(1 - s23)
    eid, emid = np.exp(1j * delta), np.exp(-1j * delta)
    u = np.array([
        [c12 * c13, s12_ * c13, s13_ * emid],
        [-s12_ * c23 - c12 * s23_ * s13_ * eid,
         c12 * c23 - s12_ * s23_ * s13_ * eid, s23_ * c13],
        [s12_ * s23_ - c12 * c23 * s13_ * eid,
         -c12 * s23_ - s12_ * c23 * s13_ * eid, c23 * c13]], dtype=complex)
    return u @ np.diag([1.0, np.exp(0.5j * a21), np.exp(0.5j * a31)])


def left_rotation(y):
    vals, vecs = np.linalg.eigh(y @ y.conjugate().T)
    order = np.argsort(vals)
    return vecs[:, order]


def takagi_masses_vectors(m):
    vals, vecs = np.linalg.eigh(m.conjugate().T @ m)
    order = np.argsort(vals)
    return np.sqrt(np.maximum(vals[order], 0.0)), vecs[:, order]


def mixing_angles(u):
    p2 = np.abs(u) ** 2
    s13 = float(p2[0, 2])
    return {"s12": float(p2[0, 1] / max(1 - s13, 1e-30)), "s13": s13,
            "s23": float(p2[1, 2] / max(1 - s13, 1e-30))}


U_e = left_rotation(Y_e)
m_diag = np.array([OLD["m1_eV"],
                   math.sqrt(OLD["m1_eV"] ** 2 + OLD["dm21"]),
                   math.sqrt(OLD["m1_eV"] ** 2 + OLD["dm31"])])
U_pmns_old = standard_pmns(OLD["s12"], OLD["s13"], OLD["s23"],
                           OLD["delta"], OLD["a21"], OLD["a31"])
U_nu = U_e @ U_pmns_old
m_light = U_nu.conjugate() @ np.diag(m_diag) @ U_nu.conjugate().T
mD_eV = Y_nu * 100.0e9                       # v_u = 100 GeV, in eV
MR_eV = -(mD_eV.T @ np.linalg.inv(m_light) @ mD_eV)
MR_GeV = MR_eV / 1.0e9

sing = np.linalg.svd(MR_GeV, compute_uv=False)
M_star = float(sing[0])
Mn = MR_GeV / M_star
c0 = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], dtype=complex)
c_hat = c0 / math.sqrt(3.0)
zeta = complex(np.vdot(c_hat, Mn))
M_V = Mn - zeta * c_hat
contact_fraction = float(np.linalg.norm(zeta * c_hat) / np.linalg.norm(Mn))

arch_zeta = complex(arch_rank["decomposition"]["contact_coefficient"]["re"],
                    arch_rank["decomposition"]["contact_coefficient"]["im"])
check("replayed zeta equals the archived card value",
      abs(zeta - arch_zeta) < 1e-12 * abs(arch_zeta), f"zeta = {zeta!r}")
check("zeta reproduces the PAPER reproducibility anchor 0.1076472949 + 0.0736514853 i",
      abs(zeta - complex(0.1076472949, 0.0736514853)) < 1e-9)
check("contact fraction reproduces the paper anchor 1.304275e-1",
      abs(contact_fraction - 0.1304275166688152) < 1e-12,
      f"{contact_fraction:.12f}")
check("Veronese + contact residual is machine-zero",
      float(np.linalg.norm(Mn - M_V - zeta * c_hat) / np.linalg.norm(Mn)) < 1e-14)
check("orientation note RECORDED: c_hat = -K_tr(paper); M_* = sigma_max(M_R)",
      np.allclose(c_hat, -np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]]) / math.sqrt(3)),
      f"M_* = {M_star:.6e} GeV")

# forward closure: seesaw with the replayed M_R returns the OLD targets exactly
def forward_observables(MR_eV_matrix):
    m_l = -(mD_eV @ np.linalg.inv(MR_eV_matrix) @ mD_eV.T)
    masses, vecs = takagi_masses_vectors(m_l)
    u = U_e.conjugate().T @ vecs
    ang = mixing_angles(u)
    return {**ang, "dm21": float(masses[1] ** 2 - masses[0] ** 2),
            "dm31": float(masses[2] ** 2 - masses[0] ** 2),
            "masses_eV": [float(x) for x in masses], "U": u}


obs_bench = forward_observables(MR_eV)
old_close = (abs(obs_bench["s12"] - OLD["s12"]) < 1e-9
             and abs(obs_bench["s13"] - OLD["s13"]) < 1e-9
             and abs(obs_bench["s23"] - OLD["s23"]) < 1e-9
             and abs(obs_bench["dm21"] / OLD["dm21"] - 1) < 1e-9
             and abs(obs_bench["dm31"] / OLD["dm31"] - 1) < 1e-9)
check("forward seesaw closes: the replayed M_R returns the OLD benchmark "
      "targets to 1e-9 (inverse <-> forward loop shut)", old_close)

heavy_masses_GeV = [float(s) for s in np.linalg.svd(MR_GeV, compute_uv=False)]
check("heavy Majorana spectrum exported for DYN-7 (leptogenesis feed)",
      all(m > 0 for m in heavy_masses_GeV),
      f"M_R singular values = {[f'{m:.3e}' for m in heavy_masses_GeV]} GeV")

# --------------------------------------------------------------- section 2
print("== DYN-4a section 2: Audit-1b rank-contract compliance ==")

E = {}
E["v1"] = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=complex)
E["v2"] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex) / math.sqrt(2)
E["v3"] = np.array([[0, 0, 1 / math.sqrt(2)], [0, math.sqrt(2), 0],
                    [1 / math.sqrt(2), 0, 0]], dtype=complex) / math.sqrt(3)
E["v4"] = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex) / math.sqrt(2)
E["v5"] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=complex)
E["contact"] = c_hat
gram_ok = all(abs(np.vdot(E[a], E[b]) - (a == b)) < 1e-12
              for a in E for b in E)
cascade = []
resid = Mn.copy()
for name in ["contact", "v3", "v1", "v2", "v4", "v5"]:
    coef = complex(np.vdot(E[name], resid))
    resid = resid - coef * E[name]
    cascade.append({"tensor": name, "coefficient": [coef.real, coef.imag],
                    "residual_after": float(np.linalg.norm(resid) / np.linalg.norm(Mn))})
check("aligned 6-tensor basis is orthonormal (5 Veronese + contact); "
      "Audit-1b requires >= 3 aligned tensors with per-tensor residuals",
      gram_ok and len(cascade) >= 3)
check("residual cascade terminates at machine zero (complete decomposition, "
      "not a hollow two-tensor test)", cascade[-1]["residual_after"] < 1e-14,
      f"cascade tail {cascade[-1]['residual_after']:.2e}")

# --------------------------------------------------------------- section 3
print("== DYN-4a section 3: target-table refresh (NuFit-6.0) ==")

NUFIT = {
    "IC24_with_SK": {"s12": (0.308, 0.0115), "s13": (0.02215, 0.00057),
                     "s23": (0.470, 0.015), "dm21": (7.49e-5, 0.19e-5),
                     "dm31": (2.513e-3, 0.020e-3), "delta_deg": (212, 33.5)},
    "IC19_without_SK": {"s12": (0.307, 0.0115), "s13": (0.02195, 0.00056),
                        "s23": (0.561, 0.0135), "dm21": (7.49e-5, 0.19e-5),
                        "dm31": (2.534e-3, 0.024e-3), "delta_deg": (177, 19.5)},
    "provenance": "NuFit-6.0, JHEP 12 (2024) 216, arXiv:2410.05380, Table 1 "
                  "(Normal Ordering); 1-sigma symmetrized",
}
PRIMARY = NUFIT["IC24_with_SK"]
shift_sigma = {k: (OLD_v - PRIMARY[k][0]) / PRIMARY[k][1]
               for k, OLD_v in [("s12", OLD["s12"]), ("s13", OLD["s13"]),
                                ("s23", OLD["s23"]), ("dm21", OLD["dm21"]),
                                ("dm31", OLD["dm31"])]}
check("target shift table computed; the sin^2 theta_23 octant flip is the "
      "dominant refresh (old benchmark sits several sigma high)",
      abs(shift_sigma["s23"]) > 3,
      f"shifts (old - new)/sigma: "
      f"{ {k: round(v, 2) for k, v in shift_sigma.items()} }")

# --------------------------------------------------------------- section 4
print("== DYN-4a section 4: zeta profile likelihood (conditional Dirac sector) ==")

OBS = ["s12", "s13", "s23", "dm21", "dm31"]


def chi2_of(obs, targets):
    return float(sum(((obs[k] - targets[k][0]) / targets[k][1]) ** 2 for k in OBS))


def observables_at_zeta(z):
    MR = M_star * 1.0e9 * (M_V + z * c_hat)      # back to eV
    return forward_observables(MR)


chi2_bench_old = chi2_of(obs_bench, {k: (OLD[{"s12": "s12", "s13": "s13",
                                              "s23": "s23", "dm21": "dm21",
                                              "dm31": "dm31"}[k]], PRIMARY[k][1])
                                     for k in OBS})
check("pipeline sanity: chi^2 of the benchmark against the OLD targets is ~ 0",
      chi2_bench_old < 1e-12, f"{chi2_bench_old:.2e}")

chi2_bench_new = chi2_of(obs_bench, PRIMARY)
pulls_bench = {k: (obs_bench[k] - PRIMARY[k][0]) / PRIMARY[k][1] for k in OBS}
check("benchmark pull table vs NuFit-6.0 DISCLOSED (theta_23 drives it)", True,
      f"chi^2 = {chi2_bench_new:.1f}; pulls "
      f"{ {k: round(v, 1) for k, v in pulls_bench.items()} }")

# Work in q = log(|zeta'|/|zeta_ref|) and the wrapped relative phase phi.
# The logarithmic radial coordinate keeps |zeta'| positive and avoids the
# coordinate singularity of polar variables at zero.  A deterministic coarse
# search discovers several basins; the needle-like benchmark basin is polished
# independently with Nelder-Mead and Powell.  This is evidence for a stable
# connected local candidate, not a proof of global basin completeness.
Q_BOUNDS = (-math.log(1000.0), math.log(2.0))
PHI_BOUNDS = (-math.pi, math.pi)


def zeta_at_qphi(q, phi):
    return (math.exp(float(q)) * abs(zeta)
            * np.exp(1j * (np.angle(zeta) + float(phi))))


def chi2_qphi(x):
    q, phi = float(x[0]), float(x[1])
    if not (Q_BOUNDS[0] <= q <= Q_BOUNDS[1]
            and PHI_BOUNDS[0] <= phi <= PHI_BOUNDS[1]):
        return 1.0e300
    try:
        value = chi2_of(observables_at_zeta(zeta_at_qphi(q, phi)), PRIMARY)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return 1.0e300
    return value if np.isfinite(value) else 1.0e300


JOINT_SEEDS = [
    (0.0, 0.0),
    (math.log(0.5), 0.0), (math.log(1.5), 0.0),
    (0.0, -0.5 * math.pi), (0.0, 0.5 * math.pi),
    (math.log(0.5), -0.5 * math.pi),
    (math.log(0.5), 0.5 * math.pi),
    (math.log(1.5), -0.5 * math.pi),
    (math.log(1.5), 0.5 * math.pi),
]
JOINT_METHODS = {
    "Nelder-Mead": {"xatol": 1.0e-12, "fatol": 1.0e-11,
                    "maxiter": 5000},
    "Powell": {"xtol": 1.0e-12, "ftol": 1.0e-12,
               "maxiter": 5000},
}
joint_runs = []
for method, options in JOINT_METHODS.items():
    for seed in JOINT_SEEDS:
        result = minimize(chi2_qphi, np.array(seed), method=method,
                          options=options)
        joint_runs.append({
            "method": method, "seed": [float(seed[0]), float(seed[1])],
            "success": bool(result.success),
            "x": [float(result.x[0]), float(result.x[1])],
            "chi2": float(result.fun), "nfev": int(result.nfev),
        })

finite_joint = [r for r in joint_runs if r["success"] and np.isfinite(r["chi2"])]
if not finite_joint:
    raise RuntimeError("joint zeta minimization produced no finite solution")
seed_best = min(finite_joint, key=lambda r: r["chi2"])
polish = {}
for method, options in JOINT_METHODS.items():
    result = minimize(chi2_qphi, np.array(seed_best["x"]), method=method,
                      options=options)
    polish[method] = result

method_minima = {
    method: min((r for r in finite_joint if r["method"] == method),
                key=lambda r: r["chi2"])
    for method in JOINT_METHODS
}
nm, pw = polish["Nelder-Mead"], polish["Powell"]
optimizer_dx = float(np.max(np.abs(nm.x - pw.x)))
optimizer_dchi2 = float(abs(nm.fun - pw.fun))
optimizer_agreement = (nm.success and pw.success
                       and optimizer_dx < 1.0e-8
                       and optimizer_dchi2 < 1.0e-6)
best_coarse_candidate_agreement = (
    abs(method_minima["Nelder-Mead"]["chi2"]
        - method_minima["Powell"]["chi2"]) < 1.0e-6
    and np.max(np.abs(np.array(method_minima["Nelder-Mead"]["x"])
                      - np.array(method_minima["Powell"]["x"]))) < 1.0e-8)
check("best coarse-search candidate agrees between two optimizers",
      optimizer_agreement and best_coarse_candidate_agreement,
      f"|Delta x|_inf = {optimizer_dx:.2e}, Delta chi2 = "
      f"{optimizer_dchi2:.2e}; {len(finite_joint)}/{len(joint_runs)} finite runs")

best_result = nm if nm.fun <= pw.fun else pw
q_best, phi_best = (float(best_result.x[0]), float(best_result.x[1]))
chi2_best = float(best_result.fun)
z_best = complex(zeta_at_qphi(q_best, phi_best))
check("adaptive joint zeta fit located a finite interior local candidate",
      np.isfinite(chi2_best)
      and Q_BOUNDS[0] < q_best < Q_BOUNDS[1]
      and PHI_BOUNDS[0] < phi_best < PHI_BOUNDS[1],
      f"q = {q_best:.8e}, phi = {phi_best:.8e}, "
      f"zeta_best = {z_best:.9f}, chi^2_min = {chi2_best:.6f} "
      f"(benchmark chi^2 = {chi2_bench_new:.6f})")

# Resolve the discovered narrow basin from a stencil whose offsets are of the
# same order as its eventual confidence widths.  Coarse starts far outside
# this basin converge to much higher stationary points, so their failure to
# hit it is recorded rather than hidden behind the word "multi-start".
LOCAL_STENCIL = [
    (dq, dp)
    for dq in (-5.0e-6, 0.0, 5.0e-6)
    for dp in (-1.0e-5, 0.0, 1.0e-5)
]
local_runs = []
for method, options in JOINT_METHODS.items():
    for dq, dp in LOCAL_STENCIL:
        result = minimize(
            chi2_qphi, np.array([q_best + dq, phi_best + dp]),
            method=method, options=options)
        local_runs.append({
            "method": method,
            "seed": [q_best + dq, phi_best + dp],
            "success": bool(result.success),
            "x": [float(result.x[0]), float(result.x[1])],
            "chi2": float(result.fun),
        })
local_dx = max(
    float(np.max(np.abs(np.array(run["x"]) - np.array([q_best, phi_best]))))
    for run in local_runs)
local_dchi2 = max(abs(run["chi2"] - chi2_best) for run in local_runs)
local_basin_stable = (
    all(run["success"] for run in local_runs)
    and local_dx < 1.0e-8 and local_dchi2 < 1.0e-6
)
coarse_basin_hits = sum(
    abs(run["chi2"] - chi2_best) < 1.0e-6 for run in finite_joint
)
check("needle-basin local stencil converges to the same candidate",
      local_basin_stable,
      f"{len(local_runs)}/{len(local_runs)} local runs; "
      f"max |Delta x|={local_dx:.2e}, max Delta chi2={local_dchi2:.2e}; "
      f"coarse basin hits={coarse_basin_hits}/{len(finite_joint)}")

obs_best = observables_at_zeta(z_best)
pulls_best = {k: (obs_best[k] - PRIMARY[k][0]) / PRIMARY[k][1] for k in OBS}
check("residual tension after the zeta-only refit DISCLOSED: with the Dirac "
      "sector frozen, zeta alone cannot absorb the theta_23 octant flip if "
      "chi^2_min stays large (a finding either way)", True,
      f"pulls at best { {k: round(v, 1) for k, v in pulls_best.items()} }")

# Genuine local profiles: at every fixed coordinate, minimize over the other
# coordinate.  The bracket is deliberately local to the connected confidence
# component around the verified local candidate.  No global-minimum claim is
# made.  A 1e-5 initial simplex resolves the narrow local
# nuisance minimum without asking a coarse bounded scalar search to find it.
PROFILE_SIMPLEX_STEP = 1.0e-5


def profile_q(q):
    result = minimize(
        lambda x: chi2_qphi((q, float(x[0]))), np.array([phi_best]),
        method="Nelder-Mead",
        options={"initial_simplex": np.array([
                     [phi_best], [phi_best + PROFILE_SIMPLEX_STEP]]),
                 "xatol": 1.0e-13, "fatol": 1.0e-11,
                 "maxiter": 1000})
    return float(result.fun), float(result.x[0])


def profile_phi(phi):
    result = minimize(
        lambda x: chi2_qphi((float(x[0]), phi)), np.array([q_best]),
        method="Nelder-Mead",
        options={"initial_simplex": np.array([
                     [q_best], [q_best + PROFILE_SIMPLEX_STEP]]),
                 "xatol": 1.0e-13, "fatol": 1.0e-11,
                 "maxiter": 1000})
    return float(result.fun), float(result.x[0])


def profile_side(profile, center, sign, xtol):
    def root_function(width):
        return profile(center + sign * width)[0] - chi2_best - 1.0

    hi = 1.0e-8
    while root_function(hi) < 0.0 and hi < 1.0e-2:
        hi *= 2.0
    if hi >= 1.0e-2 and root_function(hi) < 0.0:
        raise RuntimeError("failed to bracket a Delta-chi2=1 profile root")
    return float(bisect(root_function, 0.0, hi, xtol=xtol,
                        rtol=1.0e-12, maxiter=200))


profile_roots_coarse = {
    "q_minus": profile_side(profile_q, q_best, -1.0, 5.0e-11),
    "q_plus": profile_side(profile_q, q_best, +1.0, 5.0e-11),
    "phi_minus": profile_side(profile_phi, phi_best, -1.0, 5.0e-11),
    "phi_plus": profile_side(profile_phi, phi_best, +1.0, 5.0e-11),
}
profile_roots = {
    "q_minus": profile_side(profile_q, q_best, -1.0, 1.0e-13),
    "q_plus": profile_side(profile_q, q_best, +1.0, 1.0e-13),
    "phi_minus": profile_side(profile_phi, phi_best, -1.0, 1.0e-13),
    "phi_plus": profile_side(profile_phi, phi_best, +1.0, 1.0e-13),
}
root_refinement_rel = max(
    abs(profile_roots[k] - profile_roots_coarse[k]) / profile_roots[k]
    for k in profile_roots)
root_residual = max(
    abs(profile_q(q_best - profile_roots["q_minus"])[0] - chi2_best - 1.0),
    abs(profile_q(q_best + profile_roots["q_plus"])[0] - chi2_best - 1.0),
    abs(profile_phi(phi_best - profile_roots["phi_minus"])[0] - chi2_best - 1.0),
    abs(profile_phi(phi_best + profile_roots["phi_plus"])[0] - chi2_best - 1.0),
)
check("profile roots pass tolerance refinement and Delta-chi2 residual gates",
      root_refinement_rel < 0.01 and root_residual < 1.0e-4,
      f"max refinement shift = {root_refinement_rel:.2e}, "
      f"max root residual = {root_residual:.2e}")

# Local Hessian cross-check.  For Delta chi2 = (1/2) dx^T H dx, profiling
# over the nuisance coordinate gives the Schur-complement curvature.
def numerical_hessian(step):
    x0 = np.array([q_best, phi_best])
    h = float(step)
    value0 = chi2_qphi(x0)
    H = np.zeros((2, 2))
    for i in range(2):
        ei = np.zeros(2)
        ei[i] = h
        H[i, i] = (chi2_qphi(x0 + ei) - 2.0 * value0
                   + chi2_qphi(x0 - ei)) / h**2
    e0, e1 = np.array([h, 0.0]), np.array([0.0, h])
    H[0, 1] = H[1, 0] = (
        chi2_qphi(x0 + e0 + e1) - chi2_qphi(x0 + e0 - e1)
        - chi2_qphi(x0 - e0 + e1) + chi2_qphi(x0 - e0 - e1)
    ) / (4.0 * h**2)
    return H


def hessian_profile_scales(H):
    cq = H[0, 0] - H[0, 1] ** 2 / H[1, 1]
    cp = H[1, 1] - H[0, 1] ** 2 / H[0, 0]
    return np.array([math.sqrt(2.0 / cq), math.sqrt(2.0 / cp)])


H_FINE = numerical_hessian(1.0e-7)
H_COARSE = numerical_hessian(3.0e-7)
hessian_eigenvalues = np.linalg.eigvalsh(H_FINE)
hessian_scales = hessian_profile_scales(H_FINE)
hessian_scales_coarse = hessian_profile_scales(H_COARSE)
profile_mean_scales = np.array([
    0.5 * (profile_roots["q_minus"] + profile_roots["q_plus"]),
    0.5 * (profile_roots["phi_minus"] + profile_roots["phi_plus"]),
])
hessian_profile_rel = float(np.max(
    np.abs(hessian_scales - profile_mean_scales) / profile_mean_scales))
hessian_step_rel = float(np.max(
    np.abs(hessian_scales - hessian_scales_coarse) / hessian_scales))
check("joint minimum has a positive, step-stable Hessian consistent with profiles",
      np.all(hessian_eigenvalues > 0.0)
      and hessian_profile_rel < 0.10 and hessian_step_rel < 0.05,
      f"eig(H) = {[f'{x:.3e}' for x in hessian_eigenvalues]}, "
      f"profile-vs-Hessian = {hessian_profile_rel:.3f}, "
      f"step refinement = {hessian_step_rel:.3f}")

q_interval = [q_best - profile_roots["q_minus"],
              q_best + profile_roots["q_plus"]]
phi_interval = [phi_best - profile_roots["phi_minus"],
                phi_best + profile_roots["phi_plus"]]
abs_zeta_interval = [abs(zeta) * math.exp(q_interval[0]),
                     abs(zeta) * math.exp(q_interval[1])]
arg_zeta_interval = [np.angle(zeta) + phi_interval[0],
                     np.angle(zeta) + phi_interval[1]]
post = {
    "zeta_best": [z_best.real, z_best.imag],
    "abs_zeta_best": abs(z_best), "arg_zeta_best": float(np.angle(z_best)),
    "chi2_min": chi2_best, "chi2_benchmark": chi2_bench_new,
    "q_best": q_best, "relative_phase_best": phi_best,
    "q_68_profile": q_interval,
    "relative_phase_68_profile": phi_interval,
    "abs_zeta_68": abs_zeta_interval,
    "arg_zeta_68": arg_zeta_interval,
    "distance_from_benchmark": abs(z_best - zeta),
    "conditional_on": "benchmark Y_nu, Y_e frozen; NO; m1 = 1e-3 eV",
    "statistical_object": "connected local conditional profile likelihood; "
                          "diagonal-Gaussian NuFit summary; not a normalized "
                          "Bayesian posterior; no global basin-completeness claim",
    "thresholds": {"one_parameter_68p27": 1.0,
                   "joint_two_parameter_68p27": 2.30,
                   "joint_two_parameter_95": 5.99},
}
check("connected local zeta profile published with nuisance minimization "
      "and Delta-chi2 <= 1 one-parameter intervals",
      all(profile_roots[k] > 0.0 for k in profile_roots),
      f"|zeta| in [{post['abs_zeta_68'][0]:.5f}, {post['abs_zeta_68'][1]:.5f}], "
      f"arg in [{post['arg_zeta_68'][0]:.5f}, {post['arg_zeta_68'][1]:.5f}] rad")

# delta_CP prediction at best point (not fitted)
def jarlskog_delta(u):
    J = float(np.imag(u[0, 1] * np.conj(u[0, 2]) * np.conj(u[1, 1]) * u[1, 2]))
    a = mixing_angles(u)
    s12, s13, s23 = (math.sqrt(a["s12"]), math.sqrt(a["s13"]), math.sqrt(a["s23"]))
    c12, c13, c23 = (math.sqrt(1 - a["s12"]), math.sqrt(1 - a["s13"]),
                     math.sqrt(1 - a["s23"]))
    denom = s12 * c12 * s23 * c23 * s13 * c13**2
    return math.degrees(math.asin(max(-1.0, min(1.0, J / denom)))) if denom > 0 else float("nan")


delta_pred = jarlskog_delta(obs_best["U"])
branch_a = delta_pred % 360
branch_b = (180 - delta_pred) % 360
pull_a = min(abs(branch_a - PRIMARY["delta_deg"][0]),
             360 - abs(branch_a - PRIMARY["delta_deg"][0])) / PRIMARY["delta_deg"][1]
pull_b = min(abs(branch_b - PRIMARY["delta_deg"][0]),
             360 - abs(branch_b - PRIMARY["delta_deg"][0])) / PRIMARY["delta_deg"][1]
check("delta_CP at the best point reported as a PREDICTION vs NuFit "
      "(both sin-branches; not fitted)", True,
      f"delta in {{{branch_a:.0f}, {branch_b:.0f}}} deg vs NuFit "
      f"{PRIMARY['delta_deg'][0]} +/- {PRIMARY['delta_deg'][1]} deg "
      f"(pulls {pull_a:.1f} / {pull_b:.1f} sigma)")

# --------------------------------------------------------------- section 5
print("== DYN-4a section 5: sensitivity-window regression and refresh ==")

LOOSE_ANGLE, LOOSE_DM = 3.0e-2, 3.0e-1


def passes_loose(obs):
    return (abs(obs["s12"] - OLD["s12"]) <= LOOSE_ANGLE
            and abs(obs["s13"] - OLD["s13"]) <= LOOSE_ANGLE
            and abs(obs["s23"] - OLD["s23"]) <= LOOSE_ANGLE
            and abs(obs["dm21"] / OLD["dm21"] - 1) <= LOOSE_DM
            and abs(obs["dm31"] / OLD["dm31"] - 1) <= LOOSE_DM)


def passes_tight(obs):
    return (abs(obs["s12"] - OLD["s12"]) <= 1.0e-3
            and abs(obs["s13"] - OLD["s13"]) <= 1.0e-3
            and abs(obs["s23"] - OLD["s23"]) <= 1.0e-3
            and abs(obs["dm21"] / OLD["dm21"] - 1) <= 1.0e-2
            and abs(obs["dm31"] / OLD["dm31"] - 1) <= 1.0e-2)


def bisect_side(direction, sign, criterion):
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        z = ((1 + sign * mid) * zeta if direction == "scale"
             else zeta * np.exp(1j * sign * mid))
        try:
            ok = criterion(observables_at_zeta(z))
        except np.linalg.LinAlgError:
            ok = False
        if ok:
            lo = mid
        else:
            hi = mid
    return lo


def half_width(direction, criterion):
    up = bisect_side(direction, +1, criterion)
    dn = bisect_side(direction, -1, criterion)
    return (up + dn) / 2, up, dn


win_s, up_s, dn_s = half_width("scale", passes_loose)
win_d, up_d, dn_d = half_width("phase", passes_loose)
win_s_t, _, _ = half_width("scale", passes_tight)
win_d_t, _, _ = half_width("phase", passes_tight)
check("historical LOOSE windows REGRESS as half-widths of the two-sided "
      "interval: paper anchors 3.208798e-5 and 5.424119e-5 rad",
      abs(win_s - 3.2087975865957574e-5) < 1e-4 * 3.2e-5
      and abs(win_d - 5.424119231661329e-5) < 1e-4 * 5.4e-5,
      f"Delta_s = {win_s:.6e} (up {up_s:.4e}, down {dn_s:.4e}), "
      f"Delta_delta = {win_d:.6e}")
check("historical TIGHT windows also regress against the archived card "
      "(9.9379e-7 scale, 1.7993e-6 rad phase)",
      abs(win_s_t - 9.93792673542604e-7) < 1e-3 * 1e-6
      and abs(win_d_t - 1.799256517286863e-6) < 1e-3 * 1.8e-6,
      f"tight Delta_s = {win_s_t:.6e}, tight Delta_delta = {win_d_t:.6e}")


# Backward-compatible scalar summaries are the means of the asymmetric genuine
# profile half-widths.  Keep the two sides in the ledger; downstream consumers
# that need a conservative tolerance should use max(minus, plus), not the mean.
s_best = math.exp(q_best)
s_interval = [math.exp(q_interval[0]), math.exp(q_interval[1])]
s_minus = s_best - s_interval[0]
s_plus = s_interval[1] - s_best
phi_minus = profile_roots["phi_minus"]
phi_plus = profile_roots["phi_plus"]
win_s_new = 0.5 * (s_minus + s_plus)
win_d_new = 0.5 * (phi_minus + phi_plus)
check("refreshed connected-local profile Delta-chi2 <= 1 windows summarized without "
      "discarding their asymmetric sides",
      s_minus > 0 and s_plus > 0 and phi_minus > 0 and phi_plus > 0,
      f"Delta_s = -{s_minus:.3e}/+{s_plus:.3e}; "
      f"Delta_phi = -{phi_minus:.3e}/+{phi_plus:.3e} rad")

# --------------------------------------------------------------- ledger
npass = sum(1 for _, ok in CHECKS if ok)
report = {
    "audit": "audit1_dyn4a_seesaw_replay_zeta_posterior",
    "dyn_item": "DYN-4a",
    "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "checks_total": len(CHECKS), "checks_passed": npass,
    "all_pass": npass == len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "provenance": {
        "closure_card": {"path": str(card_path.relative_to(ROOT)),
                         "sha256": sha256(card_path),
                         "note": "user-restored archival artifact"},
        "texture_rank_summary": {"path": str(rank_path.relative_to(ROOT)),
                                 "sha256": sha256(rank_path)},
        "recipe_source": "HEAD code/audit_source_majorana_texture_rank.py + "
                         "code/verify_seesaw_item3.py (conventions "
                         "transcribed verbatim)",
        "nufit": NUFIT["provenance"],
    },
    "derivation_log": {
        "step_1_replay": "Y_e -> U_e (eigh ascending); U_nu = U_e U_PMNS(old "
                         "benchmark incl. Majorana phases 0.35pi, 1.10pi); "
                         "m_light = U_nu* D U_nu*^T; mD = Y_nu * 100 GeV; "
                         "M_R = -mD^T m_light^-1 mD (inverse type-I); "
                         "M_* = sigma_max; zeta = <c_hat, M_R/M_*>, "
                         "c_hat = -K_tr(paper)",
        "step_2_rank_contract": "orthonormal 6-tensor cascade (contact + 5 "
                                "Veronese), per-tensor coefficients and "
                                "residuals recorded",
        "step_3_refresh": "old hardcoded targets vs NuFit-6.0 NO (both "
                          "variants); theta_23 octant flip dominates",
        "step_4_profile_likelihood": "M_R(zeta') = M_*(M_V + zeta' c_hat), "
                            "forward seesaw, diagonal-Gaussian chi^2 over 5 "
                            "observables (IC24 1-sigma); q=log(|zeta'|/"
                            "|zeta_ref|), phi=arg(zeta'/zeta_ref); dual-"
                            "optimizer deterministic multi-start joint fit; "
                            "connected-local nuisance-minimized Delta-chi^2=1 "
                            "one-parameter profile intervals; global basin "
                            "completeness not claimed",
        "step_5_windows": "historical loose windows re-derived with the "
                          "archival tolerance definition (regression), then "
                          "recomputed as Delta-chi^2 windows at the new "
                          "best point",
    },
    "benchmark": {
        "zeta": [zeta.real, zeta.imag],
        "contact_fraction": contact_fraction,
        "M_star_GeV": M_star,
        "heavy_MR_singular_values_GeV": heavy_masses_GeV,
        "old_targets": OLD,
        "pulls_vs_nufit60_ic24": {k: float(v) for k, v in pulls_bench.items()},
        "chi2_vs_nufit60_ic24": chi2_bench_new,
    },
    "rank_contract_cascade": cascade,
    "target_shift_sigma": {k: float(v) for k, v in shift_sigma.items()},
    # Retain the historical key for downstream schema compatibility.  Its
    # statistical_object field explicitly says this is not a Bayesian posterior.
    "zeta_posterior": post,
    "zeta_profile_likelihood": post,
    "joint_fit_diagnostics": {
        "coordinates": ["q=log(|zeta|/|zeta_ref|)",
                        "phi=arg(zeta/zeta_ref)"],
        "bounds": {"q": list(Q_BOUNDS), "phi": list(PHI_BOUNDS)},
        "runs": joint_runs,
        "local_stencil_runs": local_runs,
        "coarse_basin_hits": coarse_basin_hits,
        "coarse_runs_total": len(finite_joint),
        "local_basin_stable": local_basin_stable,
        "optimizer_polish": {
            "Nelder-Mead": {"success": bool(nm.success),
                            "x": [float(x) for x in nm.x],
                            "chi2": float(nm.fun), "nfev": int(nm.nfev)},
            "Powell": {"success": bool(pw.success),
                       "x": [float(x) for x in pw.x],
                       "chi2": float(pw.fun), "nfev": int(pw.nfev)},
        },
        "agreement": {"max_abs_coordinate": optimizer_dx,
                      "abs_chi2": optimizer_dchi2,
                      "pass": bool(optimizer_agreement
                                   and best_coarse_candidate_agreement)},
    },
    "profile_diagnostics": {
        "definition": {"q": "min_phi chi2(q,phi)",
                       "phi": "min_q chi2(q,phi)"},
        "delta_chi2": 1.0,
        "root_half_widths": profile_roots,
        "root_half_widths_coarse": profile_roots_coarse,
        "root_refinement_max_relative": root_refinement_rel,
        "root_residual_max_abs_delta_chi2": root_residual,
        "hessian": H_FINE.tolist(),
        "hessian_eigenvalues": [float(x) for x in hessian_eigenvalues],
        "hessian_profile_scales": [float(x) for x in hessian_scales],
        "profile_mean_scales": [float(x) for x in profile_mean_scales],
        "hessian_vs_profile_max_relative": hessian_profile_rel,
        "hessian_step_refinement_max_relative": hessian_step_rel,
    },
    "pulls_at_best": {k: float(v) for k, v in pulls_best.items()},
    "delta_cp_prediction_deg_branches": [branch_a, branch_b],
    "delta_cp_branch_pulls_sigma": [float(pull_a), float(pull_b)],
    "windows": {
        "loose_regression": {"Delta_s": win_s, "Delta_delta": win_d},
        "chi2_refreshed": {
            "definition": "one-parameter nuisance-minimized profile, "
                          "Delta-chi2=1",
            "Delta_s": win_s_new, "Delta_s_minus": s_minus,
            "Delta_s_plus": s_plus, "scale_interval": s_interval,
            "Delta_delta": win_d_new, "Delta_delta_minus": phi_minus,
            "Delta_delta_plus": phi_plus,
            "relative_phase_interval": phi_interval,
        },
    },
    # negative-boundary flags
    "fit_not_prediction": True,
    "profile_likelihood_conditional_on_benchmark_dirac_sector": True,
    "profile_is_connected_local_component": True,
    "global_minimum_or_disconnected_component_completeness_claimed": False,
    "normalized_bayesian_posterior_claimed": False,
    "full_covariant_yukawa_refit_deferred_to_DYN4b": True,
    "normal_ordering_assumed_m1_fixed": True,
    "zeta_value_derived": False,
}

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "dyn4a_seesaw_zeta_posterior.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n")
(OUT / "dyn4a_seesaw_zeta_posterior.md").write_text(f"""# DYN-4a: Seesaw Replay, Target Refresh, and the zeta Profile Likelihood

`{npass}/{len(CHECKS)}` checks passed.

## Provenance closed

The pruned benchmark artifacts were restored by the author under
route_E/output/ (sha256 recorded).  Using ONLY the recipe preserved in HEAD
code, this audit regenerates the paper's reproducibility anchors from the
archival Yukawas: zeta = {zeta.real:.10f} + {zeta.imag:.10f} i and contact
fraction {contact_fraction:.7f} -- both match the printed digits.  The
forward seesaw closes on the old targets to 1e-9, and the heavy Majorana
spectrum ({', '.join(f'{m:.2e}' for m in heavy_masses_GeV)} GeV) is exported
for DYN-7.

## Target refresh (NuFit-6.0, arXiv:2410.05380, NO)

Old hardcoded benchmark vs IC24-with-SK-atm: the sin^2 theta_23 octant flip
dominates ({shift_sigma['s23']:+.1f} sigma); other shifts are mild.  Both
NuFit variants are recorded in the ledger.

## Connected local zeta profile (conditional on the benchmark Dirac sector)

Benchmark chi^2 vs NuFit-6.0 = {chi2_bench_new:.1f} (pull table in ledger,
theta_23-driven).  Jointly minimizing M_R(zeta') = M_*(M_V + zeta' c_hat):
best zeta' = {z_best.real:.6f} + {z_best.imag:.6f} i with
chi^2_min = {chi2_best:.2f}; connected-local one-parameter profile 68.27%
intervals
|zeta| in [{post['abs_zeta_68'][0]:.5f}, {post['abs_zeta_68'][1]:.5f}],
arg zeta in [{post['arg_zeta_68'][0]:.5f}, {post['arg_zeta_68'][1]:.5f}] rad;
distance from the benchmark {post['distance_from_benchmark']:.5f}.
The fit uses a diagonal-Gaussian NuFit summary.  The local stencil and profile
roots converge, but the basin is needle-like and coarse starts land on much
higher stationary points.  This is not a proof of the global minimum or of
disconnected-component completeness, and it is not a normalized Bayesian
posterior.
Residual pulls at the best point are in the ledger: with the Dirac sector
frozen, the zeta direction can only partially absorb the octant flip --
the remainder is DYN-4b's covariant refit problem.
delta_CP (sin branch) at the best point: {delta_pred:.1f} deg vs NuFit
{PRIMARY['delta_deg'][0]} +/- {PRIMARY['delta_deg'][1]} deg (prediction,
not fitted).

## Windows

Historical loose windows REGRESS to the paper anchors:
Delta_s = {win_s:.6e} (anchor 3.208798e-5),
Delta_delta = {win_d:.6e} rad (anchor 5.424119e-5).
Refreshed nuisance-minimized Delta-chi^2 <= 1 profile half-widths:
Delta_s = -{s_minus:.3e}/+{s_plus:.3e},
Delta_phi = -{phi_minus:.3e}/+{phi_plus:.3e} rad.

## Boundary

Connected local profile conditional on benchmark Y_nu, Y_e (full covariant refit =
DYN-4b); NO assumed, m1 = 1e-3 eV fixed; v_u = 100 GeV archival convention;
a fit, not a prediction -- Theorem VII.4 stands.
""")
print(f"Wrote {OUT / 'dyn4a_seesaw_zeta_posterior.json'} (+ .md)")
print(f"DYN-4a: {npass}/{len(CHECKS)} checks passed; zeta anchors regenerated; "
      f"connected local profile: chi2 {chi2_bench_new:.1f} -> {chi2_best:.1f}, "
      f"zeta_best = {z_best:.4f}.")
if npass != len(CHECKS):
    raise SystemExit(1)
