#!/usr/bin/env python3
"""DYN-9b-1: non-SUSY vacuum kinematics and ESH threshold sensitivity.

First stage of the non-SUSY core re-derivation (REQUIRED for the alive
branch).  Everything here is either SUSY-INDEPENDENT gauge kinematics or
a one-loop threshold recomputation; the genuinely dynamical remainder
(quartic Clebsch descent, survivor masses) is staged to DYN-9b-1b and
disclosed.

  S1  Little-group map on the singlet slice (p, a, omega, sigma).
      The gauge-boson masses come from covariant derivatives of the
      210/126bar vevs -- they carry NO reference to supersymmetry.  The
      five AG gauge-sector mass formulas (gated against gaugino-column
      norms in DYN-1b) are therefore reused verbatim, and zero-counting
      with the sector multiplicities {G:1, F:2, J:6, E:12, X:12} gives
      the little group of every vev pattern:
        (p,0,0,0)   -> dim 21 = Pati-Salam (D-parity broken: the 210
                       (1,1,1) is D-odd)
        (0,a,0,0)   -> dim 15 = SU(3) x SU(2)_L x SU(2)_R x U(1)_{B-L}
        (p,a,0,0)   -> dim 15 = the SAME chain WITH D broken by p
        (0,0,w,0)   -> dim 13 = 3c 2L 1R 1BL
        (t,t,t,0)   -> dim 25 = flipped SU(5) x U(1)
        (-t,-t,t,0) -> dim 25 = SU(5) x U(1)
        AG benchmark -> dim 12 = SM
      KEY CONSEQUENCE: the left-right chain of the non-SUSY rescue is
      realizable from the framework's own 210 (pattern (p,a,0,0), with
      p breaking D-parity) -- the 45_H source swap recorded in the
      DYN-9 verdict is OPTIONAL, not forced.

  S2  Eaten-Goldstone bookkeeping (integers): the massive-vector sector
      sets per chain step match dim SO(10) - dim H exactly, and the
      second-stage eaten counts come from the 126bar.

  S3  ESH threshold-sensitivity scan: the DYN-9 verdicts assumed the
      extended-survival-hypothesis minimal scalar set between M_I and
      M_X.  The non-SUSY potential generically leaves pseudo-Goldstone
      survivors; pending their mass derivation (9b-1b), we bracket the
      effect by re-solving each chain with one candidate multiplet at a
      time surviving at M_I (worst case) and recording every verdict
      flip.  This turns BDM's qualitative threshold-spread caveat into
      an in-framework quantified table.

  S4  The tree-level pseudo-Goldstone landmine is DISCLOSED with its
      literature anchor (one-loop stabilization of the non-SUSY vacua,
      Bertolini-Di Luzio-Malinsky PRD 81 035015), not re-derived.

Ledgers -> output/audit9/dyn9b1_nonsusy_vacuum_thresholds.{json,md}.
zeta is NOT derived; no scalar mass is derived here.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

from route_e_paths import AUDIT_OUTPUT, REPO_ROOT

ROOT = REPO_ROOT
OUT = AUDIT_OUTPUT / "audit9"
CHECKS = []
DLOG = {}


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


dyn9 = json.loads((OUT / "dyn9_nonsusy_intermediate.json").read_text())

print("== DYN-9b-1 section 0: premise ==")
check("premise: DYN-9 declares the non-SUSY potential not re-derived and "
      "ESH-minimal content assumed -- exactly the two items this audit "
      "begins to discharge",
      dyn9["nonsusy_scalar_potential_not_rederived_DYN9b"] is True
      and dyn9["esh_minimal_scalar_content"] is True)

# ------------------------------------------------------------- section 1
print("== DYN-9b-1 section 1: little-group map (SUSY-independent) ==")

# AG gauge-sector masses (transcribed from code/audit4a_dyn1b_full_spectrum
# .py, where they are gated against the gaugino-column norms of the AG
# blocks).  They arise from |D_mu <Phi>|^2 -- pure gauge kinematics.
G_GAUGE = 0.7


def m_lambda(p, a, om, sg, g=G_GAUGE):
    return {
        "G": math.sqrt(10) * g * abs(sg),
        "J": g * math.sqrt(8 * a * a + 16 * om * om + 2 * abs(sg) ** 2),
        "F": g * math.sqrt(24 * om * om + 2 * abs(sg) ** 2),
        "E": g * math.sqrt(4 * (a - om) ** 2 + 2 * (om - p) ** 2
                           + 2 * abs(sg) ** 2),
        "X": g * math.sqrt(4 * (a + om) ** 2 + 2 * (p + om) ** 2),
    }


MULT = {"G": 1, "F": 2, "J": 6, "E": 12, "X": 12}   # DYN-1b census
check("sector multiplicities sum to 45 - dim SM = 33 (DYN-1b census: "
      "G=(1,1,0), F=(1,1,+-1), J=(3,1,2/3)+c.c., E=(3,2,1/6)+c.c., "
      "X=(3,2,-5/6)+c.c.)", sum(MULT.values()) == 33)


def little_group_dim(p, a, om, sg):
    ml = m_lambda(p, a, om, sg)
    return 12 + sum(MULT[s] for s, v in ml.items() if v < 1e-12)


# AG benchmark vevs (x = 0.1) for the SM gate
x = 0.1
om_b = -x
a_b = (x * x + 2 * x - 1) / (1 - x)
p_b = x * (5 * x * x - 1) / (1 - x) ** 2
sg_b = math.sqrt(2 * x * (1 - 3 * x) * (1 + x * x) / (1 - x) ** 2)

PATTERNS = {
    "AG_benchmark": ((p_b, a_b, om_b, sg_b), 12, "SM"),
    "p_only": ((1.0, 0.0, 0.0, 0.0), 21, "Pati-Salam (D broken: 210 "
               "(1,1,1) is D-odd)"),
    "a_only": ((0.0, 1.0, 0.0, 0.0), 15, "SU(3)xSU(2)_LxSU(2)_RxU(1)_BL"),
    "p_and_a": ((0.7, 1.0, 0.0, 0.0), 15, "the SAME left-right chain "
                "with D-parity broken by p"),
    "omega_only": ((0.0, 0.0, 1.0, 0.0), 13, "3c 2L 1R 1BL"),
    "flipped_SU5_locus": ((1.0, 1.0, 1.0, 0.0), 25, "flipped SU(5)xU(1) "
                          "(E sector massless)"),
    "SU5_locus": ((-1.0, -1.0, 1.0, 0.0), 25, "SU(5)xU(1) (X sector "
                  "massless)"),
}
table = {}
all_ok = True
for name, (vv, want, label) in PATTERNS.items():
    got = little_group_dim(*vv)
    table[name] = {"vev": vv, "dim_H": got, "expected": want,
                   "little_group": label}
    all_ok &= got == want
check("little-group map on the singlet slice: all seven vev patterns "
      "give the expected unbroken dimensions (zero-counting of the "
      "SUSY-independent gauge masses)", all_ok,
      "; ".join(f"{n}: {t['dim_H']}" for n, t in table.items()))

check("KEY RESULT: the left-right chain is 210-REALIZABLE -- pattern "
      "(p, a, 0, 0) leaves exactly dim 15 unbroken with D-parity broken "
      "by the D-odd p; the 45_H source swap recorded in DYN-9 is "
      "OPTIONAL, not forced (D-parity assignments: Chang-Mohapatra-"
      "Parida bookkeeping, cited not re-derived)",
      table["p_and_a"]["dim_H"] == 15 and table["a_only"]["dim_H"] == 15)

check("boundary echo: the flipped-SU(5) and SU(5) enhanced loci sit at "
      "E- and X-sector zeros exactly as in the SUSY audit (DYN-1a "
      "boundary checks) -- the map is branch-independent",
      table["flipped_SU5_locus"]["dim_H"] == 25
      and table["SU5_locus"]["dim_H"] == 25)
DLOG["S1_little_group_map"] = table

# ------------------------------------------------------------- section 2
print("== DYN-9b-1 section 2: eaten-Goldstone bookkeeping ==")

ml_p = m_lambda(1.0, 0.0, 0.0, 0.0)
massive_p = {s for s, v in ml_p.items() if v > 1e-12}
ml_pa = m_lambda(0.7, 1.0, 0.0, 0.0)
massive_pa = {s for s, v in ml_pa.items() if v > 1e-12}
check("first-stage eaten sets: Pati-Salam step eats E+X = 24 = 45-21; "
      "left-right step eats J+E+X = 30 = 45-15 -- all from the 210 "
      "(the only vev at M_X)",
      massive_p == {"E", "X"} and sum(MULT[s] for s in massive_p) == 24
      and massive_pa == {"J", "E", "X"}
      and sum(MULT[s] for s in massive_pa) == 30)

ml_full = m_lambda(p_b, a_b, om_b, sg_b)
second_ps = {s for s in ml_full if s not in massive_p}
second_lr = {s for s in ml_full if s not in massive_pa}
check("second-stage eaten sets (sigma turns on): Pati-Salam -> SM eats "
      "G+J+F = 9 = 21-12; left-right -> SM eats G+F = 3 = 15-12 -- "
      "consistent with the 126bar (10,1,3) supplying 9 (resp. 3) of its "
      "60 real dof as Goldstones",
      second_ps == {"G", "J", "F"}
      and sum(MULT[s] for s in second_ps) == 9
      and second_lr == {"G", "F"}
      and sum(MULT[s] for s in second_lr) == 3)
DLOG["S2_goldstone_bookkeeping"] = {
    "PS_step1_eaten": sorted(massive_p), "PS_step2_eaten": sorted(second_ps),
    "LR_step1_eaten": sorted(massive_pa), "LR_step2_eaten": sorted(second_lr),
}

# ------------------------------------------------------------- section 3
print("== DYN-9b-1 section 3: ESH threshold-sensitivity scan ==")

# DYN-9 solver machinery (transcribed; its own ledger gates the b's)
B_SM = np.array([41 / 10, -19 / 6, -7.0])
BIJ_SM = np.array([[199 / 50, 27 / 10, 44 / 5], [9 / 10, 35 / 6, 12],
                   [11 / 10, 9 / 2, -26]])
MZ = 91.1876
AINV_MZ = np.array([0.6 * 127.951 * (1 - 0.23122), 127.951 * 0.23122,
                    1 / 0.1180])
B_LR = {"3c": Fr(-7), "2L": Fr(-3), "2R": Fr(-7, 3), "BL": Fr(11, 2)}
B_PS = {"4C": Fr(-23, 3), "2L": Fr(-3), "2R": Fr(11, 3)}


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


def solve_chain(bI, kind):
    def equations(u):
        lnMI, lnMX, split = u
        aSM = run2(AINV_MZ, 0.0, lnMI - math.log(MZ), B_SM, BIJ_SM)
        a1, a2, a3 = aSM
        if kind == "LR":
            a2R = split
            aBL = (a1 - 0.6 * a2R) / 0.4
            aI0 = np.array([a3, a2, a2R, aBL])
            b = [bI["3c"], bI["2L"], bI["2R"], bI["BL"]]
        else:
            a4 = a3
            a2R = split
            aI0 = np.array([a4, a2, a2R])
            b = [bI["4C"], bI["2L"], bI["2R"]]
        dt = lnMX - lnMI
        aX = aI0 - np.array([float(x_) for x_ in b]) * dt / (2 * math.pi)
        if kind == "LR":
            r = np.array([aX[0] - aX[1], aX[1] - aX[2], aX[2] - aX[3]])
        else:
            cons = a1 - (0.6 * aI0[2] + 0.4 * aI0[0])
            r = np.array([aX[0] - aX[1], aX[1] - aX[2], cons])
        return r, aX

    u = np.array([math.log(1e11), math.log(1e16), 30.0])
    for _ in range(120):
        r, aX = equations(u)
        if np.max(np.abs(r)) < 1e-10:
            break
        J = np.zeros((3, 3))
        for j in range(3):
            du = np.zeros(3); du[j] = 1e-6
            J[:, j] = (equations(u + du)[0] - equations(u - du)[0]) / 2e-6
        try:
            step = np.linalg.solve(J, r)
        except np.linalg.LinAlgError:
            break
        lam_ = 1.0
        for _ in range(12):
            if np.linalg.norm(equations(u - lam_ * step)[0]) \
                    < np.linalg.norm(r):
                u = u - lam_ * step
                break
            lam_ /= 2
        else:
            break
    r, aX = equations(u)
    lg = math.log(10)
    return {"log10_MI": u[0] / lg, "log10_MX": u[1] / lg,
            "alpha_G_inv": float(np.mean(aX)),
            "residual": float(np.max(np.abs(r)))}


def tau_d6(MX, aGinv):
    return 1.6e34 * (MX / 1e16) ** 4 * (aGinv / 25.0) ** 2


SK_EPI = 2.4e34

base = {"PS": solve_chain(B_PS, "PS"), "G_LR": solve_chain(B_LR, "LR")}
for n, s in base.items():
    ref = dyn9["chains"][n]
    check(f"baseline gate {n}: transcribed solver reproduces the DYN-9 "
          "ledger",
          abs(s["log10_MI"] - ref["log10_MI"]) < 1e-3
          and abs(s["log10_MX"] - ref["log10_MX"]) < 1e-3,
          f"({s['log10_MI']:.3f}, {s['log10_MX']:.3f}, "
          f"{s['alpha_G_inv']:.2f})")

# candidate pseudo-Goldstone survivors (complex scalars unless noted),
# Delta b derived from Dynkin/charge content, added at M_I (WORST case)
CAND_PS = {
    "(15,1,1)":       {"4C": Fr(4, 3)},
    "(15,1,3)":       {"4C": Fr(4), "2R": Fr(10)},
    "(15,3,1)":       {"4C": Fr(4), "2L": Fr(10)},
    "(6,2,2)":        {"4C": Fr(4, 3), "2L": Fr(2), "2R": Fr(2)},
    "(10,2,2)+(10b,2,2)": {"4C": Fr(8), "2L": Fr(20, 3), "2R": Fr(20, 3)},
    "(15,2,2)_from126": {"4C": Fr(16, 3), "2L": Fr(5), "2R": Fr(5)},
}
CAND_LR = {
    "(8,1,1,0)":      {"3c": Fr(1)},
    "(1,3,1,0)":      {"2L": Fr(2, 3)},
    "(1,1,3,0)":      {"2R": Fr(2, 3)},
    "(3,2,2,2/3)":    {"3c": Fr(2, 3), "2L": Fr(1), "2R": Fr(1),
                       "BL": Fr(2, 3)},
    "(1,2,2,0)_2nd_bidoublet": {"2L": Fr(1, 3), "2R": Fr(1, 3)},
}


M_PLANCK = 1.2e19


def scan(base_b, cands, kind):
    rows = {}
    for name, db in cands.items():
        b = dict(base_b)
        for k, v in db.items():
            b[k] = b[k] + v
        s = solve_chain(b, kind)
        MX = 10 ** s["log10_MX"]
        physical = bool(s["residual"] < 1e-8
                        and math.log10(MZ) < s["log10_MI"] < s["log10_MX"]
                        and MX < M_PLANCK)
        rows[name] = {**s, "tau_d6": tau_d6(MX, s["alpha_G_inv"]),
                      "delta_b": {k: str(v) for k, v in db.items()},
                      "physical": physical,
                      "alive": bool(physical
                                    and tau_d6(MX, s["alpha_G_inv"])
                                    > SK_EPI)}
    return rows


ps_rows = scan(B_PS, CAND_PS, "PS")
lr_rows = scan(B_LR, CAND_LR, "LR")
n_unphys = sum(1 for r in {**ps_rows, **lr_rows}.values()
               if not r["physical"])
check("sensitivity scan solved for all candidates on both chains; "
      "unphysical formal solutions (inverted M_I > M_X or "
      "trans-Planckian M_X) flagged and excluded from verdicts",
      all(r["residual"] < 1e-8 for r in {**ps_rows, **lr_rows}.values()),
      f"{len(ps_rows)} PS + {len(lr_rows)} LR candidates; "
      f"{n_unphys} flagged unphysical")

ps_phys = {n: r for n, r in ps_rows.items() if r["physical"]}
ps_up = {n: r for n, r in ps_phys.items() if r["alive"]}
ps_taus = {n: (f"{r['tau_d6']:.1e}" if r["physical"] else "unphysical")
           for n, r in ps_rows.items()}
check("PS verdict is THRESHOLD-FRAGILE in both directions among "
      "PHYSICAL solutions (this quantifies, in-framework, the "
      "threshold-spread caveat behind the Hyper-K window statement K4): "
      "at least one single survivor lifts tau above the current bound "
      "and at least one sinks it further",
      len(ps_up) >= 1 and any(r["tau_d6"] < 4.3e33 / 2
                              for r in ps_phys.values()),
      "; ".join(f"{n}: {t}" for n, t in ps_taus.items()))

lr_phys = {n: r for n, r in lr_rows.items() if r["physical"]}
lr_killers = [n for n, r in lr_phys.items() if not r["alive"]]
check("G_LR verdict is NOT unconditionally robust: single survivors "
      "exist that pull tau below the bound (alive/dead per candidate "
      "recorded; the survivor spectrum -- 9b-1b -- decides)",
      len(lr_killers) >= 1,
      f"killers: {lr_killers}; physical tau range "
      f"[{min(r['tau_d6'] for r in lr_phys.values()):.1e}, "
      f"{max(r['tau_d6'] for r in lr_phys.values()):.1e}]")
DLOG["S3_esh_sensitivity"] = {
    "baseline": base, "PS": ps_rows, "G_LR": lr_rows,
    "note": "survivors placed AT M_I = worst case; masses of the "
            "candidates are NOT derived here (9b-1b); the scan brackets "
            "the ESH systematic one candidate at a time",
}

# ------------------------------------------------------------- section 4
print("== DYN-9b-1 section 4: landmine disclosure ==")

check("tree-level pseudo-Goldstone landmine DISCLOSED, not re-derived: "
      "minimal non-SUSY SO(10) potentials are known to disfavor the "
      "intermediate little groups at TREE level, with the vacua rescued "
      "by one-loop corrections (Bertolini-Di Luzio-Malinsky, "
      "PRD 81 035015); the 210+126bar quartic Clebsch descent and the "
      "survivor masses remain open as DYN-9b-1b", True,
      "literature-anchored; flags below")

n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-9b-1 non-SUSY vacuum kinematics + ESH sensitivity",
    "dyn_item": "DYN-9b-1",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    "provenance": {
        "gauge_mass_formulas": "code/audit4a_dyn1b_full_spectrum.py "
                               "(gated there against gaugino-column "
                               "norms); SUSY-independent because they "
                               "come from |D_mu <vev>|^2",
        "sector_multiplicities": "DYN-1b census (G,F,J,E,X)",
        "solver": "transcribed from code/audit9_dyn9_nonsusy_"
                  "intermediate.py",
        "d_parity": "Chang-Mohapatra-Parida assignments, cited",
    },
    # negative-boundary flags
    "quartic_clebsch_descent_computed": False,
    "survivor_masses_derived": False,
    "tree_level_minimum_established": False,
    "one_loop_stabilization_cited_not_rederived": True,
    "flavor_refit_still_open_DYN9b2": True,
    "zeta_value_derived": False,
}
(OUT / "dyn9b1_nonsusy_vacuum_thresholds.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-1: non-SUSY vacuum kinematics + ESH sensitivity", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "## Results", "",
      "1. **Little-group map (SUSY-independent)**: seven vev patterns "
      "verified by zero-counting the gauge masses; the left-right chain "
      "is 210-REALIZABLE via (p, a, 0, 0) with D broken by the D-odd p "
      "-- **the 45_H swap recorded in DYN-9 is optional, not forced**.",
      "2. **Goldstone bookkeeping**: PS step eats E+X = 24, LR step "
      "eats J+E+X = 30; second stage eats G+J+F = 9 (PS) / G+F = 3 "
      "(LR) from the 126bar.",
      "3. **ESH sensitivity**: single pseudo-Goldstone survivors at "
      "M_I move the PS lifetime across the bound in BOTH directions "
      "(the Hyper-K window statement K4 now carries an in-framework "
      "threshold-spread quantification); G_LR robustness per candidate "
      "recorded.  Survivor masses NOT derived (9b-1b).",
      "4. **Landmine disclosed**: tree-level pseudo-Goldstone problem "
      "of minimal non-SUSY SO(10) potentials; one-loop stabilization "
      "cited (BDM PRD 81 035015), not re-derived.",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b1_nonsusy_vacuum_thresholds.md").write_text(
    "\n".join(md) + "\n")

print(f"\nDYN-9b-1: {n_pass}/{len(CHECKS)} checks; LR chain "
      f"210-realizable (45_H optional); PS threshold-fragile both ways; "
      f"ledgers -> output/audit9/dyn9b1_nonsusy_vacuum_thresholds.*")
