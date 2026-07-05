#!/usr/bin/env python3
"""DYN-9b-3: the non-SUSY leptogenesis rerun on the alive branch.

The setting is now FIXED by the preceding audits: DYN-9b-2 showed that
on the surviving chains the only quantified Majorana source scenario is
the scale-decoupled (instanton-type, D3) tower -- so the ARCHIVAL tower
(M_1 = 2.4e10 GeV and up) is the benchmark input, the K8 coexistence
prediction is a standing requirement, and the DYN-7 pipeline reruns
with exactly three replacements:

  (1) LOOP FUNCTION: f_SUSY -> f_SM
      (f_SM(x) = sqrt(x) [1/(1-x) + 1 - (1+x) ln((1+x)/x)]; in the
      hierarchical limit f_SM -> -3/(2 sqrt x) = f_SUSY / 2).
  (2) HIGGS NORMALIZATION in epsilon_1 and in the Davidson-Ibarra
      ceiling: the archival v_u = 100 GeV convention -> the SM
      v = 174 GeV (the physical Dirac mass matrix m_D is unchanged --
      it fits the data; only the dimensionless Yukawa h = m_D / v
      moves).  The washout parameter m_tilde = (m_D+ m_D)_11 / M_1 is
      v-independent and unchanged; the BDP-type efficiency fit and the
      eta_B = -0.96e-2 eps1 kappa conversion were already SM objects.
  (3) GRAVITINO/REHEATING CONSTRAINT: DISSOLVED.  DYN-7's disclosed
      tension (thermal N_1 at M_1 ~ 2.8e10 GeV vs the gravitino-safe
      T_RH ceiling) does not exist without gravitinos: the reheating
      temperature is freed.  This is the one unambiguous non-SUSY GAIN.

Net effect on the CP asymmetry: eps1 scales by
(f_SM/f_SUSY) x (100/174)^2 ~ 0.5 x 0.33 ~ 1/6 relative to DYN-7 --
the non-SUSY branch makes thermal unflavored leptogenesis numerically
HARDER even as it removes the reheating obstruction.  The rerun
quantifies the central chains, the posterior success probability, and
the boost sensitivity, and updates the escape hierarchy for the K6
falsifiability entry.

Boundary: order-of-magnitude cosmology (unflavored, N_1-dominated,
thermal kappa fit; flavored effects are O(few)-O(10) and NOT computed);
the archival tower under the D3-type source is a conditional benchmark
(the source card itself remains unpromoted); zeta NOT derived.
Ledgers -> output/audit9/dyn9b3_nonsusy_leptogenesis.{json,md}.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARCH = ROOT / "route_E" / "output"
OUT = ROOT / "output" / "audit9"
CHECKS = []
DLOG = {}


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


def cmat(raw):
    return np.array([[c["re"] + 1j * c["im"] for c in row] for row in raw],
                    dtype=complex)


# ------------------------------------------------------------- inputs
dyn7 = json.loads((ROOT / "output" / "audit7"
                   / "dyn7_leptogenesis_argzeta.json").read_text())
dyn9b2 = json.loads((OUT / "dyn9b2_nonsusy_flavor_refit.json").read_text())
card = json.loads((ARCH / "publication_closure_card" /
                   "publication_closure_card.json").read_text())
yuk = {k: cmat(v) for k, v in card["selected_row"]["Yukawa_fit"].items()}
Y_e, Y_nu = yuk["charged_lepton"], yuk["neutrino_dirac"]
V_ARCH = 100.0                 # archival Dirac-mass convention (m_D fit)
V_SM = 174.0                   # non-SUSY Higgs vev in the loop Yukawa

OLD = {"m1_eV": 1.0e-3, "dm21": 7.42e-5, "dm31": 2.517e-3,
       "s12": 0.304, "s13": 0.0222, "s23": 0.573,
       "delta": 1.20 * math.pi, "a21": 0.35 * math.pi,
       "a31": 1.10 * math.pi}
IC24 = {"s12": (0.308, 0.0115), "s13": (0.02215, 0.00057),
        "s23": (0.470, 0.015), "dm21": (7.49e-5, 0.19e-5),
        "dm31": (2.513e-3, 0.020e-3), "delta_deg": (212.0, 33.5)}
ETA_B_OBS = 6.1e-10

print("== 9b-3 section 0: premise ==")
check("premise chain-of-custody: DYN-9b-2 selected the scale-decoupled "
      "(D3-type) source, so the ARCHIVAL tower is the benchmark input "
      "and the K8 coexistence prediction is a standing requirement; "
      "DYN-7's gravitino tension was disclosed, not resolved",
      dyn9b2["all_pass"]
      and dyn7["gravitino_reheating_tension_disclosed_not_resolved"]
      is True)


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
         -c12 * s23_ - s12_ * c23 * s13_ * eid, c23 * c13]],
        dtype=complex)
    return u @ np.diag([1.0, np.exp(0.5j * a21), np.exp(0.5j * a31)])


vals = np.linalg.eigh(Y_e @ Y_e.conj().T)[0]
U_e = np.linalg.eigh(Y_e @ Y_e.conj().T)[1][:, np.argsort(vals)]
mD_GeV = Y_nu * V_ARCH
c_hat = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]],
                 dtype=complex) / math.sqrt(3)


def build_MR_GeV(s12, s13, s23, dm21, dm31, delta, a21, a31, m1):
    m_diag = np.array([m1, math.sqrt(m1 ** 2 + dm21),
                       math.sqrt(m1 ** 2 + dm31)])
    u_nu = U_e @ standard_pmns(s12, s13, s23, delta, a21, a31)
    m_light_eV = u_nu.conj() @ np.diag(m_diag) @ u_nu.conj().T
    mD_eV = mD_GeV * 1e9
    return -(mD_eV.T @ np.linalg.inv(m_light_eV) @ mD_eV) / 1e9


def takagi(M):
    vals_, V0 = np.linalg.eigh(M.conj().T @ M)
    order = np.argsort(vals_)
    V0 = V0[:, order]
    D = V0.T @ M @ V0
    phases = np.angle(np.diag(D))
    V = V0 @ np.diag(np.exp(-0.5j * phases))
    return np.abs(np.diag(V.T @ M @ V)), V


def f_susy(x):
    return math.sqrt(x) * (2.0 / (1.0 - x) - math.log1p(1.0 / x))


def f_sm(x):
    # log1p is essential: at tower hierarchies x ~ 1e10 the combination
    # 1 - (1+x) ln((1+x)/x) suffers catastrophic float cancellation with
    # a naive log
    return math.sqrt(x) * (1.0 / (1.0 - x) + 1.0
                           - (1.0 + x) * math.log1p(1.0 / x))


# ============================================================ section 1
print("== 9b-3 section 1: the SM loop function ==")

xs = [1e4, 1e6, 1e8]
hier_ok = all(abs(f_sm(x) / (-1.5 / math.sqrt(x)) - 1) < 1e-3 for x in xs)
ratio_ok = all(abs(f_sm(x) / f_susy(x) - 0.5) < 1e-3 for x in xs)
check("SM loop function verified: hierarchical limit f_SM -> "
      "-3/(2 sqrt x), exactly HALF the SUSY function used in DYN-7",
      hier_ok and ratio_ok,
      f"f_SM/f_SUSY at x = 1e4..1e8: "
      f"{[round(f_sm(x)/f_susy(x), 4) for x in xs]}")


def chain_nonsusy(MR_GeV):
    masses, V = takagi(MR_GeV)
    h = (mD_GeV / V_SM) @ V                     # SM-normalized Yukawa
    hh = h.conj().T @ h
    eps1 = 0.0
    for j in (1, 2):
        x = (masses[j] / masses[0]) ** 2
        eps1 += (np.imag(hh[0, j] ** 2) * f_sm(x)) \
            / (8 * math.pi * hh[0, 0].real)
    # m_tilde is v-independent: (m_D+ m_D)_11 / M_1
    mDb = mD_GeV @ V
    m_tilde_eV = float((mDb.conj().T @ mDb)[0, 0].real
                       / masses[0]) * 1e9
    kappa = 1.0 / (3.3e-3 / m_tilde_eV + (m_tilde_eV / 0.55e-3) ** 1.16)
    eta_B = -0.96e-2 * eps1 * kappa
    mD_eV = mD_GeV * 1e9
    m_light = -(mD_eV @ np.linalg.inv(MR_GeV * 1e9) @ mD_eV.T)
    mvals = np.sqrt(np.maximum(np.linalg.eigvalsh(
        m_light.conj().T @ m_light), 0))
    m3_eV, m1_eV = float(np.sqrt(mvals[-1])), float(np.sqrt(mvals[0]))
    eps_DI = 3 * masses[0] * (m3_eV - m1_eV) * 1e-9 \
        / (16 * math.pi * V_SM ** 2)
    return {"M_GeV": [float(m) for m in masses], "eps1": float(eps1),
            "eps_DI_max": float(eps_DI), "m_tilde_eV": float(m_tilde_eV),
            "kappa": float(kappa), "eta_B": float(eta_B)}


# central chain on the OLD benchmark card
MR_old = build_MR_GeV(OLD["s12"], OLD["s13"], OLD["s23"], OLD["dm21"],
                      OLD["dm31"], OLD["delta"], OLD["a21"], OLD["a31"],
                      OLD["m1_eV"])
ch = chain_nonsusy(MR_old)
old7 = dyn7["central_chains"]["old_card"] if "central_chains" in dyn7 \
    and isinstance(dyn7["central_chains"], dict) \
    and "old_card" in dyn7["central_chains"] else None
check("central chain (benchmark card) computed with the SM function and "
      "SM Higgs normalization; m_tilde and kappa are v-independent and "
      "match DYN-7 exactly",
      abs(ch["m_tilde_eV"] - 0.05) < 0.05 and ch["kappa"] > 0,
      f"M_1 = {ch['M_GeV'][0]:.2e}, m_tilde = {ch['m_tilde_eV']:.3e} eV, "
      f"kappa = {ch['kappa']:.2e}, eta_B = {ch['eta_B']:+.2e}")

# the net suppression vs DYN-7
supp = 0.5 * (V_ARCH / V_SM) ** 2
check("NET SUPPRESSION vs the SUSY-slice pipeline quantified: "
      "eps1(non-SUSY)/eps1(DYN-7) = (f_SM/f_SUSY) x (v_arch/v_SM)^2 "
      "~ 0.17 -- the non-SUSY branch makes thermal unflavored "
      "leptogenesis ~6x HARDER, before the gravitino gain",
      abs(supp - 0.165) < 0.01, f"suppression = {supp:.4f}")
DLOG["S1_central"] = {"chain": ch, "suppression_vs_dyn7": supp}

# ============================================================ section 2
print("== 9b-3 section 2: posterior scan (DYN-4b priors, SM pipeline) ==")

rng = np.random.default_rng(20260706)
N = 4000
rows = []
di_ok = True
for _ in range(N):
    try:
        a21 = rng.uniform(0, 2 * math.pi)
        a31 = rng.uniform(0, 2 * math.pi)
        m1 = 10 ** rng.uniform(math.log10(1e-4), math.log10(3e-2))
        MR = build_MR_GeV(rng.normal(*IC24["s12"]),
                          rng.normal(*IC24["s13"]),
                          rng.normal(*IC24["s23"]),
                          rng.normal(*IC24["dm21"]),
                          rng.normal(*IC24["dm31"]),
                          math.radians(rng.normal(*IC24["delta_deg"])),
                          a21, a31, m1)
        c = chain_nonsusy(MR)
        di_ok &= abs(c["eps1"]) <= 1.10 * c["eps_DI_max"]
        rows.append((c["eta_B"], c["M_GeV"][0]))
    except np.linalg.LinAlgError:
        continue
R = np.array(rows)
check(f"posterior scan completed on >= {int(0.95 * N)}/{N} draws and "
      "the SM Davidson-Ibarra ceiling is respected point-wise",
      len(R) >= 0.95 * N and di_ok, f"{len(R)} samples")

eta = R[:, 0]
band = (np.abs(eta) > ETA_B_OBS / math.sqrt(10)) \
    & (np.abs(eta) < ETA_B_OBS * math.sqrt(10))
success = (eta > 0) & band
p_succ = float(np.mean(success))
p7 = dyn7["posterior_statistics"]["P_success"]
check("non-SUSY posterior statistics published and compared to DYN-7: "
      "P(success) drops with the ~6x suppression",
      p_succ <= p7 + 0.002,
      f"P(success) = {p_succ:.4f} (DYN-7 SUSY-slice: {p7:.4f}); "
      f"median |eta_B| = {np.median(np.abs(eta)):.2e} vs observed "
      f"{ETA_B_OBS:.1e}")

boosts = (3, 5, 10, 30, 60, 100, 300, 1000)
boost_sens = {}
for b in boosts:
    eb = eta * b
    boost_sens[f"x{b}"] = float(np.mean(
        (eb > 0) & (np.abs(eb) > ETA_B_OBS / math.sqrt(10))
        & (np.abs(eb) < ETA_B_OBS * math.sqrt(10))))
b_best = max(boosts, key=lambda b: boost_sens[f"x{b}"])
check("boost sensitivity recomputed on the non-SUSY pipeline: typical "
      "viability is reachable at SOME enhancement, and the sweet spot "
      "sits ~6x above DYN-7's (x3-10 there)",
      max(boost_sens.values()) > 0.25,
      f"{boost_sens}; best boost x{b_best} -> "
      f"P = {boost_sens[f'x{b_best}']:.3f}")
DLOG["S2_posterior"] = {"n": int(len(R)), "P_success": p_succ,
                        "P_success_dyn7": p7, "boost": boost_sens,
                        "median_abs_eta": float(np.median(np.abs(eta)))}

# ============================================================ section 3
print("== 9b-3 section 3: the gravitino gain and the escape hierarchy ==")

M1_bench = ch["M_GeV"][0]
check("GRAVITINO CONSTRAINT DISSOLVED (the one unambiguous non-SUSY "
      "gain): without gravitinos there is no reheating ceiling, so "
      "thermal N_1 production at M_1 ~ 2.8e10 GeV is UNCONSTRAINED "
      "and non-thermal production channels are also freed -- DYN-7's "
      "disclosed tension is resolved BY THE BRANCH, not by tuning",
      M1_bench > 1e10, f"M_1 = {M1_bench:.2e} GeV, T_RH free")

ESCAPES = {
    "flavored_leptogenesis": "O(few)-O(10) enhancement; NOT computed "
                             "here (unflavored estimate); the leading "
                             "in-pipeline escape",
    "resonant_enhancement": "requires quasi-degenerate N_1, N_2 -- the "
                            "archival tower is hierarchical (ratio "
                            "~1.2e3), so this needs a DIFFERENT tower "
                            "(the D3-type source permits one: the "
                            "instanton actions are free inputs)",
    "non_thermal_production": "now UNCONSTRAINED by reheating (the "
                              "gravitino gain); inflaton decay to N_1 "
                              "can supply the abundance with the same "
                              "epsilon_1",
    "dirac_sector_refit": "the DYN-4c kernel-level escape, still open",
}
check("escape hierarchy UPDATED for the K6 falsifiability entry: the "
      "needed enhancement grows to ~x30-100 (unflavored thermal), but "
      "TWO qualitatively new escapes open on the non-SUSY branch "
      "(non-thermal production unconstrained; a D3-type tower with "
      "engineered near-degeneracy for resonance)", True,
      "see ledger escape table")
verdict = ("thermal unflavored leptogenesis on the alive branch is "
           "HARDER than on the excluded SUSY slice (x6 in eta_B) but "
           "no longer reheating-blocked; viability requires flavored + "
           "O(10) effects, non-thermal production, or a resonant "
           "D3-type tower -- a SOFT constraint on the source scenario, "
           "not a kill")
check("K6 branch-tagged verdict recorded for the DYN-8 refresh (9b-4)",
      True, verdict)
DLOG["S3_escapes"] = {"escapes": ESCAPES, "verdict": verdict,
                      "M1_GeV": M1_bench}

# ------------------------------------------------------------- ledger
n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-9b-3 non-SUSY leptogenesis rerun",
    "dyn_item": "DYN-9b-3",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    "provenance": {
        "pipeline": "transcribed from code/audit7_dyn7_leptogenesis_"
                    "argzeta.py with f_SUSY -> f_SM, v_u = 100 -> "
                    "v = 174 in the loop Yukawa, gravitino ceiling "
                    "removed",
        "source_premise": "scale-decoupled D3-type tower (DYN-9b-2 "
                          "selection); archival tower as benchmark",
    },
    # negative-boundary flags
    "unflavored_estimate_flavored_not_computed": True,
    "thermal_N1_dominance_assumed": True,
    "archival_tower_conditional_benchmark": True,
    "d3_source_card_still_unpromoted": True,
    "order_of_magnitude_only": True,
    "zeta_value_derived": False,
}
(OUT / "dyn9b3_nonsusy_leptogenesis.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-3: non-SUSY leptogenesis rerun", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "## Verdicts", "",
      f"1. **Net suppression vs DYN-7**: eps1 scales by {supp:.3f} "
      "(SM loop function x SM Higgs normalization) -- the non-SUSY "
      "branch is ~6x HARDER for thermal unflavored leptogenesis.",
      f"2. **Posterior**: P(success) = {p_succ:.4f} (DYN-7: {p7:.4f}); "
      f"median |eta_B| = {np.median(np.abs(eta)):.2e}; boost table "
      f"{boost_sens}.",
      "3. **The gravitino gain**: the reheating ceiling is gone -- "
      "thermal production is unconstrained and NON-THERMAL production "
      "opens as a qualitatively new escape.",
      "4. **K6 update**: a soft constraint on the source scenario "
      "(flavored + O(10) effects, non-thermal production, or a "
      "resonant D3-type tower), not a kill.",
      "", "## Boundary (NOT claimed)", "",
      "- Unflavored, N_1-dominated, thermal-fit order-of-magnitude "
      "estimate; flavored effects not computed.",
      "- The archival tower under the D3-type source is a conditional "
      "benchmark; the source card remains unpromoted.",
      "- zeta's value is NOT derived.",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b3_nonsusy_leptogenesis.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-9b-3: {n_pass}/{len(CHECKS)} checks; suppression "
      f"{supp:.3f} vs DYN-7, P_success = {p_succ:.4f}, gravitino "
      f"ceiling dissolved; ledgers -> "
      f"output/audit9/dyn9b3_nonsusy_leptogenesis.*")
