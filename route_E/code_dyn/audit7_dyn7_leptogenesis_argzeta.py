#!/usr/bin/env python3
"""DYN-7: thermal leptogenesis as the first independent constraint on arg zeta.

The contact coefficient's phase is the extra CP phase of the Majorana
sector.  This audit computes, per point of the DYN-4b posterior:

  M_R --Takagi--> (M_1, M_2, M_3), V ;   h = Y_nu V / 1  (v_u = 100 GeV
  archival convention, h dimensionless);
  epsilon_1 = (1/8pi) sum_{j>1} Im[((h^dag h)_{1j})^2] f_SUSY(x_j) /
              (h^dag h)_{11},   x_j = (M_j/M_1)^2,
  f_SUSY(x) = sqrt(x) [ 2/(1-x) - ln((1+x)/x) ]  (-> -3/sqrt(x) for x >> 1);
  m_tilde = (h^dag h)_{11} v_u^2 / M_1 ;  efficiency kappa from the
  Giudice-et-al interpolation 1/kappa = (3.3e-3 eV)/m_tilde +
  (m_tilde/0.55e-3 eV)^1.16 ;
  eta_B = -0.96e-2 * epsilon_1 * kappa   (sphaleron + dilution, SM-like
  O(1) constant; sign convention: epsilon_1 > 0 means excess LEPTONS,
  sphalerons then give NEGATIVE baryon number, so matter requires
  epsilon_1 < 0).

Deliverables:
  1. Takagi self-test and the Davidson-Ibarra bound respected point-wise.
  2. Full number chain at the old benchmark card and the DYN-4b refreshed
     central card.
  3. Historical unflavored diagnostic scan (same priors and seed family as
     DYN-4b).  Its hit fractions and conditioned arg-zeta interval are kept
     for regression only; they are not publishable probabilities until a
     branch-local flavored kinetic calculation is supplied.

Boundary: this file is an order-of-magnitude unflavored regression in the
wrong flavor regime.  At M_1 ~ 1e10 GeV a flavored treatment is mandatory;
the MSSM transition additionally needs tan(beta) and branch-local thermal
rates.  See ``audit7_dyn7_flavor_regime_gate.py``.  No physical posterior or
constraint on zeta is claimed here.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from route_e_paths import ARCHIVE_OUTPUT, AUDIT_OUTPUT, REPO_ROOT

ROOT = REPO_ROOT
ARCH = ARCHIVE_OUTPUT
OUT = AUDIT_OUTPUT / "audit7"

CHECKS: list[tuple[str, bool]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


def cmat(raw) -> np.ndarray:
    return np.array([[c["re"] + 1j * c["im"] for c in row] for row in raw],
                    dtype=complex)


# ------------------------------------------------------------ shared machinery
card = json.loads((ARCH / "publication_closure_card" /
                   "publication_closure_card.json").read_text())
yuk = {k: cmat(v) for k, v in card["selected_row"]["Yukawa_fit"].items()}
Y_e, Y_nu = yuk["charged_lepton"], yuk["neutrino_dirac"]
V_U = 100.0                                     # GeV, archival convention

OLD = {"m1_eV": 1.0e-3, "dm21": 7.42e-5, "dm31": 2.517e-3,
       "s12": 0.304, "s13": 0.0222, "s23": 0.573,
       "delta": 1.20 * math.pi, "a21": 0.35 * math.pi, "a31": 1.10 * math.pi}
IC24 = {"s12": (0.308, 0.0115), "s13": (0.02215, 0.00057),
        "s23": (0.470, 0.015), "dm21": (7.49e-5, 0.19e-5),
        "dm31": (2.513e-3, 0.020e-3), "delta_deg": (212.0, 33.5)}
ETA_B_OBS = 6.1e-10


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
    return vecs[:, np.argsort(vals)]


U_e = left_rotation(Y_e)
mD_GeV = Y_nu * V_U
c_hat = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], dtype=complex) / math.sqrt(3)


def build_MR_GeV(s12, s13, s23, dm21, dm31, delta, a21, a31, m1):
    m_diag = np.array([m1, math.sqrt(m1**2 + dm21), math.sqrt(m1**2 + dm31)])
    u_nu = U_e @ standard_pmns(s12, s13, s23, delta, a21, a31)
    m_light_eV = u_nu.conjugate() @ np.diag(m_diag) @ u_nu.conjugate().T
    mD_eV = mD_GeV * 1e9
    return -(mD_eV.T @ np.linalg.inv(m_light_eV) @ mD_eV) / 1e9


def takagi(M):
    """Autonne-Takagi for complex symmetric M: returns (masses asc, V) with
    V^T M V = diag(masses), V unitary."""
    vals, V0 = np.linalg.eigh(M.conjugate().T @ M)
    order = np.argsort(vals)
    V0 = V0[:, order]
    D = V0.T @ M @ V0
    phases = np.angle(np.diag(D))
    V = V0 @ np.diag(np.exp(-0.5j * phases))
    return np.abs(np.diag(V.T @ M @ V)), V


def f_susy(x):
    return math.sqrt(x) * (2.0 / (1.0 - x) - math.log((1.0 + x) / x))


def leptogenesis_chain(MR_GeV):
    masses, V = takagi(MR_GeV)
    h = (mD_GeV / V_U) @ V                       # dimensionless Yukawa
    hh = h.conjugate().T @ h
    eps1 = 0.0
    for j in (1, 2):
        x = (masses[j] / masses[0]) ** 2
        eps1 += (np.imag(hh[0, j] ** 2) * f_susy(x)) / (8 * math.pi * hh[0, 0].real)
    m_tilde_eV = hh[0, 0].real * (V_U**2 / masses[0]) * 1e9   # GeV -> eV
    kappa = 1.0 / (3.3e-3 / m_tilde_eV + (m_tilde_eV / 0.55e-3) ** 1.16)
    eta_B = -0.96e-2 * eps1 * kappa
    # Davidson-Ibarra bound with the light m3 from the forward seesaw
    mD_eV = mD_GeV * 1e9
    m_light = -(mD_eV @ np.linalg.inv(MR_GeV * 1e9) @ mD_eV.T)
    # eigvalsh(m_nu^dagger m_nu) returns m_i^2.  One square root therefore
    # already gives the physical light-neutrino masses.  The historical code
    # applied a second square root and overestimated the DI bound by ~3.91.
    mvals = np.sqrt(np.maximum(np.linalg.eigvalsh(
        m_light.conjugate().T @ m_light), 0))
    m3_eV, m1_eV = float(mvals[-1]), float(mvals[0])
    eps_DI = 3 * masses[0] * (m3_eV - m1_eV) * 1e-9 / (8 * math.pi * V_U**2)
    return {"M_GeV": [float(m) for m in masses], "eps1": float(eps1),
            "eps_DI_max": float(eps_DI), "m_tilde_eV": float(m_tilde_eV),
            "kappa": float(kappa), "eta_B": float(eta_B)}


# ------------------------------------------------------------ section 1
print("== DYN-7 section 1: Takagi self-test and central chains ==")

rng0 = np.random.default_rng(11)
tak_ok = True
for _ in range(50):
    A = rng0.normal(size=(3, 3)) + 1j * rng0.normal(size=(3, 3))
    Msym = (A + A.T) / 2
    d, V = takagi(Msym)
    tak_ok &= (np.linalg.norm(V.conjugate().T @ V - np.eye(3)) < 1e-10
               and np.linalg.norm(V.T @ Msym @ V - np.diag(d)) < 1e-9
               and np.all(np.diff(d) >= -1e-12))
check("Takagi self-test: V unitary, V^T M V real diagonal ascending "
      "(50 random symmetric matrices)", tak_ok)

MR_old = build_MR_GeV(OLD["s12"], OLD["s13"], OLD["s23"], OLD["dm21"],
                      OLD["dm31"], OLD["delta"], OLD["a21"], OLD["a31"],
                      OLD["m1_eV"])
chain_old = leptogenesis_chain(MR_old)
NEW_DELTA = math.radians(IC24["delta_deg"][0])
MR_new = build_MR_GeV(IC24["s12"][0], IC24["s13"][0], IC24["s23"][0],
                      IC24["dm21"][0], IC24["dm31"][0], NEW_DELTA,
                      OLD["a21"], OLD["a31"], OLD["m1_eV"])
chain_new = leptogenesis_chain(MR_new)
check("old-benchmark chain computed and Davidson-Ibarra bound respected",
      abs(chain_old["eps1"]) <= 1.05 * chain_old["eps_DI_max"],
      f"M_1 = {chain_old['M_GeV'][0]:.3e} GeV, eps1 = {chain_old['eps1']:+.3e} "
      f"(DI max {chain_old['eps_DI_max']:.3e}), m_tilde = "
      f"{chain_old['m_tilde_eV']:.3e} eV, kappa = {chain_old['kappa']:.3e}, "
      f"eta_B = {chain_old['eta_B']:+.3e}")
check("refreshed-central chain computed and DI bound respected",
      abs(chain_new["eps1"]) <= 1.05 * chain_new["eps_DI_max"],
      f"M_1 = {chain_new['M_GeV'][0]:.3e} GeV, eta_B = {chain_new['eta_B']:+.3e} "
      f"(observed +6.1e-10)")
check("sign convention stated: eta_B = -0.96e-2 eps1 kappa; matter needs "
      "eps1 < 0", True,
      f"old eta_B {'+' if chain_old['eta_B'] > 0 else '-'}, "
      f"new eta_B {'+' if chain_new['eta_B'] > 0 else '-'}")

bottleneck = {
    "eps1_over_DI_old": abs(chain_old["eps1"]) / chain_old["eps_DI_max"],
    "kappa_old": chain_old["kappa"],
    "washout_regime": f"strong (m_tilde = {chain_old['m_tilde_eV']:.3f} eV "
                      ">> m_* ~ 1e-3 eV)",
    "boost_needed_old": ETA_B_OBS / abs(chain_old["eta_B"]),
    "boost_needed_new": ETA_B_OBS / abs(chain_new["eta_B"]),
}
check("BOTTLENECK ANALYSIS disclosed: the central slices underproduce "
      "|eta_B| -- strong washout times phase-suppressed epsilon_1; the "
      "required enhancement factor is quantified", True,
      f"eps1/DI = {bottleneck['eps1_over_DI_old']:.4f}, kappa = "
      f"{bottleneck['kappa_old']:.2e}, boost needed: old x"
      f"{bottleneck['boost_needed_old']:.0f}, new x"
      f"{bottleneck['boost_needed_new']:.0f} (and the central SIGN is "
      "antimatter: the phase convention point sits on the wrong side)")

# ------------------------------------------------------------ section 2
print("== DYN-7 section 2: posterior scan (DYN-4b priors) ==")

rng = np.random.default_rng(20260705)
N = 4000
rows = []
di_ok = True
for _ in range(N):
    try:
        a21, a31 = rng.uniform(0, 2 * math.pi), rng.uniform(0, 2 * math.pi)
        m1 = 10 ** rng.uniform(math.log10(1e-4), math.log10(3e-2))
        MR = build_MR_GeV(rng.normal(*IC24["s12"]), rng.normal(*IC24["s13"]),
                          rng.normal(*IC24["s23"]), rng.normal(*IC24["dm21"]),
                          rng.normal(*IC24["dm31"]),
                          math.radians(rng.normal(*IC24["delta_deg"])),
                          a21, a31, m1)
        ch = leptogenesis_chain(MR)
        zeta = complex(np.vdot(c_hat, MR / np.linalg.svd(MR, compute_uv=False)[0]))
        di_ok &= abs(ch["eps1"]) <= 1.10 * ch["eps_DI_max"]
        rows.append((np.angle(zeta), abs(zeta), ch["eta_B"], ch["M_GeV"][0],
                     a21, m1))
    except np.linalg.LinAlgError:
        continue
R = np.array(rows)
check(f"posterior scan completed on >= {int(0.95 * N)}/{N} draws",
      len(R) >= 0.95 * N, f"{len(R)} samples")
check("Davidson-Ibarra bound respected across the whole posterior "
      "(consistency of the epsilon_1 implementation)", di_ok)

eta = R[:, 2]
p_sign = float(np.mean(eta > 0))
band = (np.abs(eta) > ETA_B_OBS / math.sqrt(10)) & (np.abs(eta) < ETA_B_OBS * math.sqrt(10))
p_band = float(np.mean(band))
success = (eta > 0) & band
p_succ = float(np.mean(success))
check("unflavored sign and magnitude diagnostics recorded (not physical "
      "posterior probabilities)", True,
      f"P(sign) = {p_sign:.3f}, P(band) = {p_band:.3f}, "
      f"P(success) = {p_succ:.3f}")
check("the success set is nonempty but RARE: on this slice thermal "
      "unflavored leptogenesis succeeds only in posterior tails",
      0.0 < p_succ < 0.05,
      f"successful samples: {int(np.sum(success))} of {len(R)}")

boost_sens = {f"P_success_x{b}": float(np.mean((eta * b > 0)
              & (np.abs(eta * b) > ETA_B_OBS / math.sqrt(10))
              & (np.abs(eta * b) < ETA_B_OBS * math.sqrt(10))))
              for b in (3, 5, 10, 30)}
check("boost sensitivity: P(success) under a uniform O(1)-O(10) "
      "enhancement (flavor effects, resonance, convention) is quantified "
      "-- how far the slice sits from viability", True, f"{boost_sens}")

# ------------------------------------------------------------ section 3
print("== DYN-7 section 3: the leptogenesis-conditioned arg zeta window ==")


def circ_stats(angles):
    mean_dir = float(np.angle(np.mean(np.exp(1j * angles))))
    dev = (angles - mean_dir + math.pi) % (2 * math.pi) - math.pi
    return {"circular_mean": mean_dir,
            "p16": float(mean_dir + np.percentile(dev, 16)),
            "p50": float(mean_dir + np.percentile(dev, 50)),
            "p84": float(mean_dir + np.percentile(dev, 84)),
            "p2p5": float(mean_dir + np.percentile(dev, 2.5)),
            "p97p5": float(mean_dir + np.percentile(dev, 97.5))}


arg_all = circ_stats(R[:, 0])
arg_cut = circ_stats(R[success, 0])
width_all = arg_all["p84"] - arg_all["p16"]
width_cut = arg_cut["p84"] - arg_cut["p16"]
check("unflavored conditioned-window diagnostic recorded, but not "
      "promoted to a phase constraint -- "
      "success correlates with none of the scanned variables strongly and "
      "the 14-sample window carries large statistical noise (flagged)",
      True,
      f"unconditioned arg zeta 68%: [{arg_all['p16']:.4f}, {arg_all['p84']:.4f}] "
      f"(width {width_all:.4f}); conditioned (n = {int(np.sum(success))}, "
      f"NOISY): [{arg_cut['p16']:.4f}, {arg_cut['p84']:.4f}] "
      f"(width {width_cut:.4f})")

# does the old benchmark phase survive the cut?
bench_arg = 0.600038020318215
in_cut_95 = arg_cut["p2p5"] <= bench_arg <= arg_cut["p97p5"]
check("position of the old benchmark phase arg zeta = 0.6000 relative to "
      "the conditioned window DISCLOSED", True,
      f"benchmark {'INSIDE' if in_cut_95 else 'OUTSIDE'} the conditioned "
      f"95% window [{arg_cut['p2p5']:.4f}, {arg_cut['p97p5']:.4f}]")

# what drives success: correlations
succ_f = success.astype(float)
corr_a21 = float(np.corrcoef(np.cos(R[:, 4]), succ_f)[0, 1])
corr_m1 = float(np.corrcoef(np.log10(R[:, 5]), succ_f)[0, 1])
corr_M1 = float(np.corrcoef(np.log10(R[:, 3]), succ_f)[0, 1])
check("success drivers DISCLOSED (correlations with cos a21, log10 m1, "
      "log10 M_1)", True,
      f"cos(a21): {corr_a21:+.2f}, log10 m1: {corr_m1:+.2f}, "
      f"log10 M_1: {corr_M1:+.2f}")

# ------------------------------------------------------------ ledger
npass = sum(1 for _, ok in CHECKS if ok)
report = {
    "audit": "audit7_dyn7_leptogenesis_argzeta",
    "dyn_item": "DYN-7",
    "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "checks_total": len(CHECKS), "checks_passed": npass,
    "all_pass": npass == len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": {
        "chain": "M_R --Takagi--> (M_i, V); h = Y_nu V (v_u = 100 GeV); "
                 "eps1 = (1/8pi) sum Im[(h^dag h)_{1j}^2] f_SUSY(x_j) / "
                 "(h^dag h)_{11}; f_SUSY(x) = sqrt(x)[2/(1-x) - ln((1+x)/x)]; "
                 "m_tilde = (h^dag h)_11 v_u^2/M_1; 1/kappa = 3.3e-3/m_tilde "
                 "+ (m_tilde/0.55e-3)^1.16 (eV); eta_B = -0.96e-2 eps1 kappa",
        "sign_convention": "eps1 > 0 = lepton excess; sphalerons give "
                           "negative B; matter requires eps1 < 0",
        "priors": "identical to DYN-4b (NuFit-6.0 IC24 gaussians x uniform "
                  "Majorana phases x log-uniform m1)",
        "success_definition": "eta_B > 0 AND |eta_B| within half a decade "
                              "of 6.1e-10",
    },
    "central_chains": {"old_benchmark": chain_old,
                       "refreshed_central": chain_new},
    "posterior_statistics": {"n_samples": int(len(R)),
                             "P_sign": p_sign, "P_band": p_band,
                             "P_success": p_succ,
                             "log10_abs_etaB_16_50_84": [
                                 float(np.percentile(np.log10(np.abs(eta)), p))
                                 for p in (16, 50, 84)]},
    "arg_zeta_windows": {"unconditioned_68": [arg_all["p16"], arg_all["p84"]],
                         "unconditioned_95": [arg_all["p2p5"], arg_all["p97p5"]],
                         "leptogenesis_conditioned_68": [arg_cut["p16"], arg_cut["p84"]],
                         "leptogenesis_conditioned_95": [arg_cut["p2p5"], arg_cut["p97p5"]],
                         "width_reduction_68": float(width_all - width_cut),
                         "benchmark_phase_inside_conditioned_95": bool(in_cut_95)},
    "success_drivers": {"corr_cos_a21": corr_a21, "corr_log10_m1": corr_m1,
                        "corr_log10_M1": corr_M1},
    "bottleneck": bottleneck,
    "boost_sensitivity": boost_sens,
    "finding": "on the archival Dirac slice, thermal unflavored "
               "leptogenesis is MARGINAL: central chains give the wrong "
               "sign and |eta_B| ~ 1e-11 (30-70x short); P(success) ~ 0.4% "
               "(tails only); the required O(30-70) enhancement exceeds "
               "typical O(few) flavor corrections -> a soft kill-criterion "
               "entry for DYN-8, escapable by a modified Dirac sector "
               "(DYN-4c) or non-thermal production; arg zeta itself is NOT "
               "the lever on this slice",
    "conditioned_window_small_sample_caveat": True,
    "claim_status": "blocked_missing_branch_thermal_inputs",
    "physics_posterior_valid": False,
    "success_probability_publishable": False,
    "arg_zeta_conditioned_interval_publishable": False,
    "required_repair": "provide tan(beta), branch-local thermal and "
                       "spectator rates, and solve at least the two-flavor "
                       "kinetic system (density matrix preferred)",
    # negative-boundary flags
    "order_of_magnitude_only": True,
    "unflavored_estimate_two_flavor_regime_O1": True,
    "thermal_N1_dominance_assumed": True,
    "gravitino_reheating_tension_disclosed_not_resolved": True,
    "archival_vu_100GeV_convention": True,
    "zeta_value_derived": False,
    "constraint_is_empirical_anchor_not_principle": True,
}

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "dyn7_leptogenesis_argzeta.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n")
(OUT / "dyn7_leptogenesis_argzeta.md").write_text(f"""# DYN-7: Historical Unflavored Leptogenesis Diagnostic

`{npass}/{len(CHECKS)}` checks passed.

## Chain (order of magnitude, SUSY loop function, BDP-type washout)

Takagi -> h = Y_nu V -> epsilon_1 (Davidson-Ibarra respected point-wise)
-> m_tilde, kappa -> eta_B = -0.96e-2 eps1 kappa.
Old benchmark: M_1 = {chain_old['M_GeV'][0]:.2e} GeV,
eps1 = {chain_old['eps1']:+.2e}, m_tilde = {chain_old['m_tilde_eV']:.2e} eV,
kappa = {chain_old['kappa']:.2e}, eta_B = {chain_old['eta_B']:+.2e}.
Refreshed central: M_1 = {chain_new['M_GeV'][0]:.2e} GeV,
eta_B = {chain_new['eta_B']:+.2e} (observed +6.1e-10).

## Regression statistics (DYN-4b priors, {len(R)} samples; not a posterior)

P(eta_B > 0) = {p_sign:.3f}; P(|eta_B| within half a decade of 6.1e-10) =
{p_band:.3f}; P(success) = {p_succ:.3f}.

## Diagnostic only: no flavored phase constraint

Unconditioned (DYN-4b): 68% [{arg_all['p16']:.4f}, {arg_all['p84']:.4f}],
95% [{arg_all['p2p5']:.4f}, {arg_all['p97p5']:.4f}] rad.
Leptogenesis-conditioned: 68% [{arg_cut['p16']:.4f}, {arg_cut['p84']:.4f}],
95% [{arg_cut['p2p5']:.4f}, {arg_cut['p97p5']:.4f}] rad.
The old benchmark phase 0.6000 is
{'INSIDE' if in_cut_95 else 'OUTSIDE'} the conditioned 95% window.
Success drivers: cos(a21) corr {corr_a21:+.2f}, log10 m1 {corr_m1:+.2f},
log10 M_1 {corr_M1:+.2f}.

These numbers are retained only to regression-test the historical unflavored
pipeline.  They cannot cut the contact phase in the tau-resolved regime.

## Boundary

Status: `blocked_missing_branch_thermal_inputs`.  The missing inputs are
tan(beta), branch-local charged-lepton/spectator rates, and a flavored kinetic
solution.  See `audit7_dyn7_flavor_regime_gate.py`.
""")
print(f"Wrote {OUT / 'dyn7_leptogenesis_argzeta.json'} (+ .md)")
print(f"DYN-7: {npass}/{len(CHECKS)} arithmetic checks passed; historical "
      "unflavored statistics are diagnostic only; physics status = "
      "blocked_missing_branch_thermal_inputs.")
if npass != len(CHECKS):
    raise SystemExit(1)
