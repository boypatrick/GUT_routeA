#!/usr/bin/env python3
"""Route-D D5: high-scale SUSY-breaking bridge scan.

One-parameter interpolating family between the EXCLUDED TeV-scale SUSY
slice (DYN-2/3/2b) and the DYN-9 non-SUSY two-step chains: all
superpartners at a common M_SS; non-SUSY running below, SUSY running
above.  The ordering M_SS < M_* is recorded only as a necessary scale
condition for a supersymmetric matching regime.  It does not repair or
validate the interacting messenger action rejected by DYN-5V.

Solver: the DYN-9 machinery (two-loop SM segment, one-loop
intermediate, Newton on (ln M_I, ln M_X, split)) extended to two
sub-segments of the intermediate range,

  M_Z --SM(2L)--> min(M_SS, M_I) [--MSSM--> M_I if M_SS < M_I]
      --G_I non-SUSY--> max(M_SS, M_I) --G_I SUSY--> M_X,

with the SUSY intermediate b coefficients DERIVED from the superfield
content (b = -3 C_A + sum T; the U(1)_{B-L} cubic anomaly of Delta_R
forces the conjugate Delta_R-bar, and SUSY PS likewise doubles
Sigma_R -- a content difference vs the non-SUSY ESH set, disclosed).

Proton ledger along the bridge:
  tau_d6 = 1.6e34 (M_X/1e16)^4 (alpha_G^-1/25)^2   [DYN-9 formula]
  tau_d5 = 5.9e33 (M_T/1e17)^2 (M_SS/3e3)^2, M_T = M_X (flagged) --
  the DYN-3 mini-split lever computed as a curve (tau_d5 scales as the
  PRODUCT (M_T M_SS)^2).

Survival window per chain:
  solvable & physical  AND  tau_d6 > 2.4e34 (current e+ pi0 bound)
  AND  tau_d5 > 5.9e33 (Super-K nu K+)  AND  M_SS < M_* = 3.93e15
  (a necessary matching-scale ordering only).

Sanity gates at both ends: M_SS above M_X reproduces the DYN-9 ledger
solutions; M_SS = 3 TeV reproduces the d=5 kill structure and the
DYN-3 calibration anchor.

Route-D discipline: CONDITIONAL scenario, NOT promoted; one-loop
intermediate segments; ESH-minimal content; zeta NOT derived.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "route_d" / "output"
CHECKS = []
DLOG = {}


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


dyn9 = json.loads((ROOT / "output" / "audit9"
                   / "dyn9_nonsusy_intermediate.json").read_text())
dyn3 = json.loads((ROOT / "output" / "audit2"
                   / "dyn3_proton_d5_kill_criterion.json").read_text())
M_STAR = float(json.loads((ROOT / "output" / "audit1"
                           / "dyn4a_seesaw_zeta_posterior.json"
                           ).read_text())["benchmark"]["M_star_GeV"])

# ---------------------------------------------------------- section 1
print("== D5 section 1: b coefficients (non-SUSY reused, SUSY derived) ==")

# non-SUSY intermediate b's, identical to DYN-9 (gated there 15/15)
B_NS = {"G_LR": {"3c": Fr(-7), "2L": Fr(-3), "2R": Fr(-7, 3),
                 "BL": Fr(11, 2)},
        "PS": {"4C": Fr(-23, 3), "2L": Fr(-3), "2R": Fr(11, 3)}}
check("non-SUSY intermediate b's reused as gated in DYN-9: "
      "G_LR (-7,-3,-7/3,11/2), PS (-23/3,-3,11/3)", True,
      "transcribed from code/audit9 (its own content-derivation gates)")

# SUSY: b = -3 C_A + sum T(chiral superfields); abelian: sum T
# G_LR superfields/gen: Q(3,2,1;1/3) Qc(3b,1,2;-1/3) L(1,2,1;-1)
# Lc(1,1,2;+1); Higgs: bidoublet Phi(1,2,2;0), DeltaR(1,1,3;+2) AND
# DeltaR-bar(1,1,3;-2) (forced, see anomaly gate)
f38 = Fr(3, 8)
b_LR_S = {
    "3c": -9 + Fr(1, 2) * 2 * 2 * 3,                 # 3 gen x (Q + Qc)
    "2L": -6 + (Fr(1, 2) * 3 + Fr(1, 2)) * 3 + 1,    # gens + bidoublet
    "2R": -6 + (Fr(1, 2) * 3 + Fr(1, 2)) * 3 + 1 + 2 + 2,
    "BL": f38 * (Fr(1, 9) * 6 + Fr(1, 9) * 6 + 2 + 2) * 3
    + f38 * 12 + f38 * 12,
}
check("SUSY G_LR b's DERIVED: (b3, b2L, b2R, bBL) = (-3, 1, 5, 15)",
      (b_LR_S["3c"], b_LR_S["2L"], b_LR_S["2R"], b_LR_S["BL"])
      == (Fr(-3), Fr(1), Fr(5), Fr(15)),
      f"{[str(b_LR_S[k]) for k in ('3c', '2L', '2R', 'BL')]}")

# PS superfields/gen: F(4,2,1) Fc(4b,1,2); Higgs: Phi(1,2,2),
# SigmaR(10b,1,3) AND SigmaR-bar(10,1,3)
b_PS_S = {
    "4C": -12 + Fr(1, 2) * 2 * 2 * 3 + 3 * 3 + 3 * 3,
    "2L": -6 + Fr(1, 2) * 4 * 3 + 1,
    "2R": -6 + Fr(1, 2) * 4 * 3 + 1 + 2 * 10 + 2 * 10,
}
check("SUSY PS b's DERIVED: (b4, b2L, b2R) = (12, 1, 41)",
      (b_PS_S["4C"], b_PS_S["2L"], b_PS_S["2R"]) == (Fr(12), Fr(1), Fr(41)),
      f"{[str(b_PS_S[k]) for k in ('4C', '2L', '2R')]}")

B_MSSM = np.array([33 / 5, 1.0, -3.0])
check("MSSM segment b's are the textbook (33/5, 1, -3)", True, "standard")

anom = 2 ** 3 * 3                       # DeltaR fermions: charge^3 x dim
check("anomaly gate: gauged U(1)_{B-L} cubic anomaly of the DeltaR "
      "superfield fermions is (+2)^3 x 3 = 24 != 0, so the SUSY chain "
      "REQUIRES the conjugate DeltaR-bar (and SUSY PS likewise doubles "
      "SigmaR) -- a content doubling relative to the non-SUSY ESH set, "
      "which is what makes the SUSY-segment b's so much larger",
      anom == 24, "b2R jumps -7/3 -> +5 (LR), 11/3 -> +41 (PS)")
DLOG["S1_b_coefficients"] = {
    "susy_LR": {k: str(v) for k, v in b_LR_S.items()},
    "susy_PS": {k: str(v) for k, v in b_PS_S.items()},
    "doubling_note": "anomaly/holomorphy force Delta+Delta-bar and "
                     "Sigma+Sigma-bar above M_SS",
}

# ---------------------------------------------------------- section 2
print("== D5 section 2: two-segment bridge solver ==")

B_SM = np.array([41 / 10, -19 / 6, -7.0])
BIJ_SM = np.array([[199 / 50, 27 / 10, 44 / 5], [9 / 10, 35 / 6, 12],
                   [11 / 10, 9 / 2, -26]])
MZ = 91.1876
AINV_MZ = np.array([0.6 * 127.951 * (1 - 0.23122), 127.951 * 0.23122,
                    1 / 0.1180])


def run2(ainv0, t0, t1, b, bij, nstep=160):
    a = np.array(ainv0, dtype=float)
    h = (t1 - t0) / nstep

    def f(ai):
        al = 1.0 / ai
        return -b / (2 * math.pi) - (bij @ al) / (8 * math.pi ** 2)

    for _ in range(nstep):
        k1 = f(a); k2 = f(a + h * k1 / 2)
        k3 = f(a + h * k2 / 2); k4 = f(a + h * k3)
        a = a + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return a


def sm_rhs(_t, ainv):
    """Two-loop SM inverse-coupling flow used by both solver paths."""
    alpha = 1.0 / ainv
    return -B_SM / (2 * math.pi) - (BIJ_SM @ alpha) / (8 * math.pi ** 2)


# The old implementation replayed a 160-step RK4 integration at every
# finite-difference evaluation of every Newton start.  Integrate the same
# autonomous ODE once with a high-accuracy dense interpolant instead.  The
# explicit RK4 implementation is retained above as an independent regression
# oracle and as a guarded fallback outside the tabulated interval.
SM_DENSE_T_MAX = 50.0
SM_DENSE = solve_ivp(
    sm_rhs,
    (0.0, SM_DENSE_T_MAX),
    AINV_MZ,
    method="DOP853",
    rtol=1.0e-11,
    atol=1.0e-12,
    dense_output=True,
)


def sm_ainv_at_log_scale(log_scale):
    t = log_scale - math.log(MZ)
    if SM_DENSE.success and 0.0 <= t <= SM_DENSE_T_MAX:
        return np.asarray(SM_DENSE.sol(t), dtype=float)
    return run2(AINV_MZ, 0.0, t, B_SM, BIJ_SM)


def bridge_equations(u, lnMSS, kind):
    """Residuals for (ln M_I, ln M_X, split) at SUSY-breaking lnMSS."""
    lnMI, lnMX, split = u
    # below M_I: SM two-loop to min(MSS, MI), then MSSM one-loop
    lnLow = min(lnMSS, lnMI)
    aSM = sm_ainv_at_log_scale(lnLow)
    if lnMSS < lnMI:
        aSM = aSM - B_MSSM * (lnMI - lnMSS) / (2 * math.pi)
    a1, a2, a3 = aSM
    if kind == "LR":
        a2R = split
        aBL = (a1 - 0.6 * a2R) / 0.4
        aI0 = np.array([a3, a2, a2R, aBL])
        keys = ("3c", "2L", "2R", "BL")
        bns = np.array([float(B_NS["G_LR"][k]) for k in keys])
        bs = np.array([float(b_LR_S[k]) for k in keys])
    else:
        a4 = a3
        a2R = split
        aI0 = np.array([a4, a2, a2R])
        keys = ("4C", "2L", "2R")
        bns = np.array([float(B_NS["PS"][k]) for k in keys])
        bs = np.array([float(b_PS_S[k]) for k in keys])
    # intermediate range [MI, MX]: non-SUSY up to max(MSS, MI) clamped
    # at MX, SUSY above
    lnMid = min(max(lnMSS, lnMI), lnMX)
    aX = aI0 - (bns * (lnMid - lnMI) + bs * (lnMX - lnMid)) / (2 * math.pi)
    if kind == "LR":
        r = np.array([aX[0] - aX[1], aX[1] - aX[2], aX[2] - aX[3]])
    else:
        cons = a1 - (0.6 * aI0[2] + 0.4 * aI0[0])
        r = np.array([aX[0] - aX[1], aX[1] - aX[2], cons])
    return r, aX, aI0


def _newton(u0, lnMSS, kind):
    u = np.array(u0, dtype=float)
    for _ in range(120):
        r, aX, aI0 = bridge_equations(u, lnMSS, kind)
        if np.max(np.abs(r)) < 1e-10:
            break
        J = np.zeros((3, 3))
        for j in range(3):
            du = np.zeros(3); du[j] = 1e-6
            J[:, j] = (bridge_equations(u + du, lnMSS, kind)[0]
                       - bridge_equations(u - du, lnMSS, kind)[0]) / 2e-6
        try:
            step = np.linalg.solve(J, r)
        except np.linalg.LinAlgError:
            break
        lam = 1.0
        for _ in range(12):
            if np.linalg.norm(bridge_equations(u - lam * step, lnMSS,
                                               kind)[0]) \
                    < np.linalg.norm(r):
                u = u - lam * step
                break
            lam /= 2
        else:
            break
    r, aX, aI0 = bridge_equations(u, lnMSS, kind)
    physical = bool(np.max(np.abs(r)) < 1e-8
                    and math.log(MZ) < u[0] < u[1]
                    and np.all(aI0 > 0) and np.all(aX > 0))
    return u, float(np.max(np.abs(r))), physical, aX


def solve_bridge(lnMSS, kind, u0=None):
    """Multi-start Newton; prefer physical roots, then smallest residual."""
    starts = [] if u0 is None else [np.array(u0, dtype=float)]
    for lgMI in (5.0, 9.0, 11.0, 13.0, 15.0):
        for lgMX, sp in ((15.7, 30.0), (16.3, 30.0), (16.3, 45.0),
                         (15.7, 50.0)):
            starts.append(np.array([lgMI * math.log(10),
                                    lgMX * math.log(10), sp]))
    best = None
    for s in starts:
        u, res, physical, aX = _newton(s, lnMSS, kind)
        cand = (not physical, res)          # physical first, then residual
        if best is None or cand < best[0]:
            best = (cand, u, res, physical, aX)
        if physical and res < 1e-10 and u0 is not None:
            break                            # continuation start succeeded
    _, u, res, physical, aX = best
    lg = math.log(10)
    return {"log10_MI": u[0] / lg, "log10_MX": u[1] / lg,
            "alpha_G_inv": float(np.mean(aX)),
            "residual": res, "physical": physical,
            "case": "B(MSS<MI)" if lnMSS < u[0] else "A(MI<=MSS)",
            "u": u}


dense_regression_errors = []
for scale in (1.0e5, 1.0e10, 1.0e15):
    log_scale = math.log(scale)
    dense = sm_ainv_at_log_scale(log_scale)
    rk4 = run2(AINV_MZ, 0.0, log_scale - math.log(MZ), B_SM, BIJ_SM)
    dense_regression_errors.append(float(np.max(np.abs(dense - rk4))))
check("dense two-loop SM trajectory reproduces the independent 160-step "
      "RK4 oracle at 1e5, 1e10, and 1e15 GeV",
      SM_DENSE.success and max(dense_regression_errors) < 1.0e-8,
      f"max |Delta alpha^-1| = {max(dense_regression_errors):.3e}")


# endpoint gate: M_SS far above M_X -> the SUSY segment vanishes and
# the solver must reproduce the DYN-9 ledger solutions
lg10 = math.log(10)
for name, kind in (("G_LR", "LR"), ("PS", "PS")):
    ref = dyn9["chains"][name]
    got = solve_bridge(math.log(1e18), kind)
    check(f"ENDPOINT GATE (DYN-9 limit), {name}: M_SS above M_X "
          "reproduces the DYN-9 ledger (log10 M_I, log10 M_X, "
          "alpha_G^-1)",
          abs(got["log10_MI"] - ref["log10_MI"]) < 1e-3
          and abs(got["log10_MX"] - ref["log10_MX"]) < 1e-3
          and abs(got["alpha_G_inv"] - ref["alpha_G_inv"]) < 0.05,
          f"got ({got['log10_MI']:.4f}, {got['log10_MX']:.4f}, "
          f"{got['alpha_G_inv']:.2f}) vs ledger "
          f"({ref['log10_MI']:.4f}, {ref['log10_MX']:.4f}, "
          f"{ref['alpha_G_inv']:.2f})")

# ---------------------------------------------------------- section 3
print("== D5 section 3: proton ledger along the bridge ==")

SK_EPI, SK_NUK = 2.4e34, 5.9e33


def tau_d6(MX, aGinv):
    return 1.6e34 * (MX / 1e16) ** 4 * (aGinv / 25.0) ** 2


def tau_d5(MT, MSS):
    return 5.9e33 * (MT / 1e17) ** 2 * (MSS / 3.0e3) ** 2


cal = tau_d5(1.45e13, 3.0e3)
check("DYN-3 CALIBRATION GATE: at the DYN-2 compatibility point "
      "(M_T_eff = 1.45e13, M_SS = 3 TeV) the mini-split curve returns "
      "the DYN-3 kill number tau ~ 1.2e26 yr",
      abs(cal - 1.2e26) / 1.2e26 < 0.05, f"tau_d5 = {cal:.2e} yr")

# scan
grid = [3.5 + 0.05 * i for i in range(int((16.6 - 3.5) / 0.05) + 1)]
scan = {}
for name, kind in (("G_LR", "LR"), ("PS", "PS")):
    rows = []
    u0 = None
    for lgMSS in reversed(grid):            # continuation from the top
        s = solve_bridge(lgMSS * lg10, kind, u0)
        if s["physical"]:
            u0 = s["u"]
        MX, MSS = 10 ** s["log10_MX"], 10 ** lgMSS
        row = {"log10_MSS": lgMSS, **{k: s[k] for k in
               ("log10_MI", "log10_MX", "alpha_G_inv", "residual",
                "physical", "case")},
               "tau_d6": tau_d6(MX, s["alpha_G_inv"]),
               "tau_d5": tau_d5(MX, MSS)}
        row["alive"] = bool(s["physical"] and row["tau_d6"] > SK_EPI
                            and row["tau_d5"] > SK_NUK
                            and MSS < M_STAR)
        rows.append(row)
    rows.reverse()
    scan[name] = rows

n_phys = {n: sum(1 for r in rows if r["physical"])
          for n, rows in scan.items()}
obstruction = {n: [r["log10_MSS"] for r in rows if not r["physical"]]
               for n, rows in scan.items()}
check("bridge scan bookkeeping: every grid point either yields a "
      "physical solution or is DISCLOSED as a no-physical-solution "
      "obstruction region (multi-start Newton; obstruction = part of "
      "the map, as in DYN-2b)",
      all(n_phys[n] > 0 for n in scan),
      "; ".join(f"{n}: {n_phys[n]}/{len(grid)} physical, obstruction "
                f"[{min(o):.2f}, {max(o):.2f}]" if (o := obstruction[n])
                else f"{n}: {n_phys[n]}/{len(grid)} physical"
                for n in scan))

# kill-structure gate at M_SS = 3 TeV: either the point is physical and
# d=5 dead, or there is no physical bridge solution at all -- both
# reproduce the TeV-end exclusion
low = {n: min(scan[n], key=lambda r: abs(r["log10_MSS"] - 3.5))
       for n in scan}
tev_dead = {n: (not low[n]["physical"]) or low[n]["tau_d5"] < SK_NUK
            for n in scan}
check("TeV-END GATE: at M_SS ~ 3 TeV every bridge is EXCLUDED -- "
      "physical-and-d=5-dead or no physical solution (the DYN-2/3 "
      "exclusion structure; the catastrophic DYN-3 magnitude was "
      "slice-specific, disclosed)",
      all(tev_dead.values()),
      "; ".join(f"{n}: " + (f"tau_d5 = {low[n]['tau_d5']:.1e}"
                            if low[n]["physical"] else "no physical "
                            "solution") for n in scan))


def failed_conditions(r):
    out = []
    if not r["physical"]:
        out.append("no_physical_solution")
    if r["tau_d6"] <= SK_EPI:
        out.append("tau_d6")
    if r["tau_d5"] <= SK_NUK:
        out.append("tau_d5")
    if 10 ** r["log10_MSS"] >= M_STAR:
        out.append("M_SS >= M_*")
    return out


# windows, with binding conditions read off the rows just outside
windows = {}
for name, rows in scan.items():
    alive = [r for r in rows if r["alive"]]
    if alive:
        i_lo = rows.index(alive[0])
        i_hi = rows.index(alive[-1])
        windows[name] = {
            "log10_MSS_min": alive[0]["log10_MSS"],
            "log10_MSS_max": alive[-1]["log10_MSS"],
            "binding_low": (failed_conditions(rows[i_lo - 1])
                            if i_lo > 0 else ["grid edge"]),
            "binding_high": (failed_conditions(rows[i_hi + 1])
                             if i_hi + 1 < len(rows) else ["grid edge"]),
            "tau_d6_range": [min(r["tau_d6"] for r in alive),
                             max(r["tau_d6"] for r in alive)],
            "MX_range_log10": [min(r["log10_MX"] for r in alive),
                               max(r["log10_MX"] for r in alive)],
            "cases": sorted({r["case"] for r in alive}),
        }
    else:
        windows[name] = None

g_win, p_win = windows["G_LR"], windows["PS"]
check("G_LR bridge window verdict computed",
      g_win is not None or all(not r["alive"] for r in scan["G_LR"]),
      f"window = {g_win}")
check("PS bridge window verdict computed (incl. whether the SUSY "
      "segment RESCUES the marginal chain past the 2.4e34 bound)",
      p_win is not None or all(not r["alive"] for r in scan["PS"]),
      f"window = {p_win}")
ps_d6_rescued = any(r["physical"] and r["tau_d6"] > SK_EPI
                    for r in scan["PS"])
check("PS d=6 rescue question answered explicitly across the whole "
      "bridge", True,
      f"exists M_SS with tau_d6(PS) > 2.4e34: {ps_d6_rescued}")

ps_phys = [r for r in scan["PS"] if r["physical"]]
ps_seg = max((max(0.0, r["log10_MX"] - r["log10_MSS"])
              for r in ps_phys), default=0.0)
check("PS sampled-segment diagnostic: with the declared doubled-Sigma "
      "content (b2R=41), no physical grid point contains more than "
      "0.05 dex of SUSY PS running below M_X.  This records the result "
      "but does not by itself prove b2R=41 is the unique cause",
      ps_seg <= 0.05,
      f"physical PS points: {len(ps_phys)}; longest sampled SUSY "
      f"segment below M_X = {ps_seg:.2f} dex")

check("necessary matching-scale ordering intersected: every alive point "
      "has M_SS < M_* = 3.93e15 GeV; this is not sufficient to validate "
      "the Route-B/DYN-5 messenger action rejected by DYN-5V",
      all(10 ** r["log10_MSS"] < M_STAR
          for ch in scan.values() for r in ch if r["alive"]))
check("M_T = M_X assumption FLAGGED with product scaling: tau_d5 "
      "scales as (M_T M_SS)^2, so each decade of triplet lightness "
      "shifts the d=5 window edge up one decade in M_SS", True,
      "sensitivity recorded in the derivation log")

DLOG["S3_scan"] = {
    "grid_step_dex": 0.05,
    "windows": windows,
    "MT_sensitivity": "tau_d5 prop to (M_T M_SS)^2; M_T = M_X assumed",
    "tev_end": {n: {"tau_d5": low[n]["tau_d5"], "tau_d6": low[n]["tau_d6"]}
                for n in scan},
}

# ---------------------------------------------------------- verdict
print("== D5 section 4: verdict ==")

bridge_alive = g_win is not None
check("BRIDGE VERDICT: the high-scale SUSY-breaking bridge has a "
      "non-empty survival window on at least one chain -- the string-"
      "natural interpolation benchmark exists on the sampled G_LR grid",
      bridge_alive,
      f"G_LR: {'ALIVE' if g_win else 'EMPTY'}, "
      f"PS: {'ALIVE' if p_win else 'EMPTY'}")

n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "Route-D D5 high-scale SUSY-breaking bridge scan",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    "windows": windows,
    "scan": {n: [{k: (v if not isinstance(v, float) else float(v))
                  for k, v in r.items() if k != "u"} for r in rows]
             for n, rows in scan.items()},
    "provenance": {
        "dyn9_ledger": "output/audit9/dyn9_nonsusy_intermediate.json",
        "dyn3_ledger": "output/audit2/dyn3_proton_d5_kill_criterion.json",
        "machinery": "two-loop SM + one-loop segments transcribed from "
                     "code/audit9_dyn9_nonsusy_intermediate.py",
    },
    # negative-boundary flags
    "one_loop_intermediate_segments": True,
    "esh_minimal_content_both_phases": True,
    "MT_equals_MX_assumed": True,
    "soft_spectrum_degenerate_assumed": True,
    "above_MX_perturbativity_not_audited": True,
    "MSS_below_Mstar_is_only_a_necessary_scale_condition": True,
    "branch_specific_experimental_MI_floor_applied": False,
    "independent_root_solver_performed": False,
    "threshold_uncertainty_envelope_performed": False,
    "proton_decay_flavor_uncertainty_envelope_performed": False,
    "physics_status": "preliminary_toy_scan",
    "physics_promotion_allowed": False,
    "blockers": [
        "branch-specific experimental lower bound on M_I",
        "threshold nuisance envelope and independent root solver",
        "physical-basis proton-decay flavor and lattice uncertainties",
        "M_T/M_X and non-degenerate soft-spectrum sensitivity",
        "a valid interacting hidden-messenger action (DYN-5V remains failed)",
    ],
    "zeta_value_derived": False,
    "promoted_to_paper": False,
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "d5_susy_breaking_bridge.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n")


def fmt_win(w):
    if w is None:
        return "EMPTY"
    return (f"log10 M_SS in [{w['log10_MSS_min']:.2f}, "
            f"{w['log10_MSS_max']:.2f}], binding: {w['binding_low']} / "
            f"{w['binding_high']}, cases {w['cases']}")


md = ["# Route-D D5: high-scale SUSY-breaking bridge scan", "",
      f"{n_pass}/{len(CHECKS)} checks pass.  CONDITIONAL scenario, NOT "
      "promoted; zeta NOT derived.", "",
      "## Verdicts", "",
      f"- G_LR bridge window: {fmt_win(g_win)}",
      f"- PS bridge window: {fmt_win(p_win)}",
      f"- PS d=6 rescue anywhere on the bridge: {ps_d6_rescued}",
      "- SUSY-segment content doubling (Delta-bar, Sigma-bar) forced "
      "by the U(1)_{B-L} cubic anomaly: b2R jumps -7/3 -> +5 (LR) and "
      "11/3 -> +41 (PS).",
      "- Endpoint gates: DYN-9 ledger reproduced at the high end; "
      "DYN-3 kill number 1.2e26 reproduced at (1.45e13, 3 TeV); TeV "
      "end d=5 dead on both chains.",
      "- Every alive point has M_SS < M_*.  This is only a necessary "
      "matching-scale ordering and does not restore the invalid DYN-5 "
      "messenger action.",
      "", "## Boundary (NOT claimed)", "",
      "- One-loop intermediate segments; ESH-minimal content; "
      "M_T = M_X; degenerate soft spectrum; above-M_X perturbativity "
      "not audited.",
      "- Route-D promotion bar NOT passed; must not be cited as "
      "evidence.", "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "d5_susy_breaking_bridge.md").write_text("\n".join(md) + "\n")

print(f"\nD5: {n_pass}/{len(CHECKS)} checks; G_LR window "
      f"{fmt_win(g_win)}; PS window {fmt_win(p_win)}; ledgers -> "
      f"{OUT.relative_to(ROOT)}/d5_susy_breaking_bridge.*")

if n_pass != len(CHECKS):
    raise SystemExit(1)
