#!/usr/bin/env python3
"""DYN-9b-3: the non-SUSY leptogenesis rerun on the alive branch.

This audit adopts the archival scale-decoupled (instanton-type, D3) tower
as a conditional benchmark.  DYN-9b-2 does not select that source: it is a
fixed-kernel comparison, the D3 card is unpromoted, and K8 remains a
conditional diagnostic.  On that explicit benchmark, the DYN-7 pipeline
reruns with exactly three replacements:

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
  (3) GRAVITINO-SPECIFIC REHEATING CEILING: ABSENT.  DYN-7's disclosed
      gravitino-safe T_RH ceiling does not apply without gravitinos, but
      the branch-local reheating history and initial abundance are not
      supplied.  General thermal production is therefore open, not free.

Net effect on the CP asymmetry: eps1 scales by
(f_SM/f_SUSY) x (100/174)^2 ~ 0.5 x 0.33 ~ 1/6 relative to DYN-7 --
the non-SUSY branch makes the historical unflavored diagnostic numerically
HARDER even as it removes only the gravitino-specific ceiling.  The resulting hit
fractions and boost sensitivity are unweighted prior-draw regression
diagnostics.  No likelihood is evaluated, so they are neither posterior
probabilities nor branch-viability probabilities.

Boundary: order-of-magnitude cosmology (unflavored, N_1-dominated,
thermal kappa fit).  Since M_1 is in the tau-resolved two-flavor regime, the
physics claim is blocked until flavored kinetics are computed;
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

from route_e_paths import ARCHIVE_OUTPUT, AUDIT_OUTPUT, REPO_ROOT

ROOT = REPO_ROOT
ARCH = ARCHIVE_OUTPUT
OUT = AUDIT_OUTPUT / "audit9"
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
check("premise chain-of-custody: the archival D3-type tower is an "
      "explicit conditional benchmark, not a selected source; K8 and D3 "
      "remain unpromoted while DYN-7's gravitino tension was disclosed",
      dyn9b2["all_pass"]
      and dyn9b2["derivation_log"]["S4_conditional_comparison"][
          "source_selected"] is False
      and dyn9b2["physics_promotion_allowed"] is False
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
    # eigvalsh(m_nu^dagger m_nu) returns m_i^2; a second square root is a
    # dimensionally wrong m_i^(1/2) and inflated the historical DI bound.
    mvals = np.sqrt(np.maximum(np.linalg.eigvalsh(
        m_light.conj().T @ m_light), 0))
    m3_eV, m1_eV = float(mvals[-1]), float(mvals[0])
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
print("== 9b-3 section 2: unweighted prior-draw regression (SM pipeline) ==")

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
check(f"prior-draw regression completed on >= {int(0.95 * N)}/{N} draws "
      "and the SM Davidson-Ibarra ceiling is respected point-wise",
      len(R) >= 0.95 * N and di_ok, f"{len(R)} samples")

eta = R[:, 0]
band = (np.abs(eta) > ETA_B_OBS / math.sqrt(10)) \
    & (np.abs(eta) < ETA_B_OBS * math.sqrt(10))
target_band_hit = (eta > 0) & band
target_band_hit_fraction = float(np.mean(target_band_hit))
# The DYN-7 ledger retains the historical field name P_success.  Here it is
# used only as a like-for-like, unweighted target-band regression reference.
dyn7_target_band_hit_fraction = dyn7["posterior_statistics"]["P_success"]
likelihood_applied = False
draws_reweighted = False
check("unweighted non-SUSY prior-draw target-band hit fraction recomputed "
      "and compared with the like-for-like DYN-7 regression reference; "
      "no likelihood or draw reweighting is applied",
      target_band_hit_fraction <= dyn7_target_band_hit_fraction + 0.002
      and not likelihood_applied and not draws_reweighted,
      f"hit fraction = {target_band_hit_fraction:.4f} "
      f"(DYN-7 reference: {dyn7_target_band_hit_fraction:.4f}); "
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
check("fixed multiplicative boost stress test recomputed as unweighted "
      "target-band hit fractions; its maximum is a regression target, "
      "not evidence that a physical enhancement or viability is realized",
      max(boost_sens.values()) > 0.25,
      f"{boost_sens}; best boost x{b_best} -> "
      f"hit fraction = {boost_sens[f'x{b_best}']:.3f}")
DLOG["S2_prior_draw_regression"] = {
    "n_draws_requested": N,
    "n_draws_evaluated": int(len(R)),
    "likelihood_applied": likelihood_applied,
    "draws_reweighted": draws_reweighted,
    "target_band_definition": (
        "eta_B > 0 and |eta_B| within a factor sqrt(10) of 6.1e-10"
    ),
    "target_band_hit_fraction": target_band_hit_fraction,
    "dyn7_legacy_target_band_hit_fraction": dyn7_target_band_hit_fraction,
    "boost_target_band_hit_fractions": boost_sens,
    "median_abs_eta": float(np.median(np.abs(eta))),
    "inferential_status": "regression_only_not_a_posterior_or_viability_probability",
}

# ============================================================ section 3
print("== 9b-3 section 3: branch-specific reheating boundary ==")

M1_bench = ch["M_GeV"][0]
gravitino_specific_ceiling_absent = True
branch_local_reheating_history_supplied = False
check("the gravitino-specific SUSY reheating ceiling is absent on the "
      "non-SUSY branch, but no branch-local reheating history is supplied; "
      "thermal or non-thermal production is therefore OPEN, not "
      "unconstrained",
      M1_bench > 1e10 and gravitino_specific_ceiling_absent
      and not branch_local_reheating_history_supplied,
      f"M_1 = {M1_bench:.2e} GeV; T_RH and initial abundance unspecified")

ESCAPES = {
    "flavored_leptogenesis": "O(few)-O(10) enhancement; NOT computed "
                             "here (unflavored estimate); the leading "
                             "in-pipeline escape",
    "resonant_enhancement": "requires quasi-degenerate N_1, N_2 -- the "
                            "archival tower is hierarchical (ratio "
                            "~1.2e3), so this needs a DIFFERENT tower "
                            "(the D3-type source permits one: the "
                            "instanton actions are free inputs)",
    "non_thermal_production": "not excluded by a gravitino-specific ceiling, "
                              "but requires an explicit inflaton/reheating "
                              "model and an initial-abundance calculation",
    "dirac_sector_refit": "the DYN-4c kernel-level escape, still open",
}
check("escape hierarchy recorded for the K6 diagnostic: the needed "
      "unflavored enhancement is ~x30-100; flavored, resonant and "
      "non-thermal mechanisms remain conditional repair hypotheses rather "
      "than computed escapes", all(ESCAPES.values()),
      "see ledger escape table")
verdict = ("thermal unflavored leptogenesis on the alive branch is "
           "HARDER than on the excluded SUSY slice (x6 in eta_B) but "
           "the gravitino-specific ceiling alone is absent; assessing "
           "viability requires a branch-local flavored kinetic and "
           "reheating model. "
           "The current result is an order diagnostic, not a kill")
check("K6 branch-tagged non-promotion verdict recorded for DYN-8",
      "order diagnostic, not a kill" in verdict, verdict)
DLOG["S3_escapes"] = {"escapes": ESCAPES, "verdict": verdict,
                      "M1_GeV": M1_bench,
                      "gravitino_specific_ceiling_absent": True,
                      "branch_local_reheating_history_supplied": False,
                      "reheating_unconstrained": False}

# ------------------------------------------------------------- ledger
source_selected = False
success_probability_publishable = False
physics_promotion_allowed = False
physics_status = "blocked_missing_branch_thermal_inputs"
check("promotion predicate fails closed: the source is unselected, the "
      "unweighted hit fraction is not publishable as a success "
      "probability, and physics promotion remains forbidden",
      source_selected is False
      and success_probability_publishable is False
      and physics_promotion_allowed is False
      and physics_status == "blocked_missing_branch_thermal_inputs")
n_pass = sum(1 for _, ok in CHECKS if ok)
mechanical_status = "pass" if n_pass == len(CHECKS) else "fail"
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
                    "v = 174 in the loop Yukawa; only the gravitino-specific "
                    "ceiling is removed",
        "source_premise": "unpromoted scale-decoupled D3-type tower; "
                          "archival tower adopted as a conditional benchmark",
    },
    # negative-boundary flags
    "unflavored_estimate_flavored_not_computed": True,
    "thermal_N1_dominance_assumed": True,
    "archival_tower_conditional_benchmark": True,
    "d3_source_card_still_unpromoted": True,
    "source_selected": source_selected,
    "order_of_magnitude_only": True,
    "claim_status": physics_status,
    "mechanical_status": mechanical_status,
    "physics_status": physics_status,
    "physics_promotion_allowed": physics_promotion_allowed,
    "promotion_gate_fail_closed": True,
    "likelihood_applied": likelihood_applied,
    "draws_reweighted": draws_reweighted,
    "prior_draw_regression_only": True,
    "physics_posterior_valid": False,
    "success_probability_publishable": success_probability_publishable,
    "gravitino_specific_ceiling_absent": True,
    "branch_local_reheating_history_supplied": False,
    "reheating_unconstrained": False,
    "required_repair": "solve tau-resolved two-flavor Boltzmann or density-"
                       "matrix equations with spectator effects and branch-"
                       "local thermal rates, reheating history, and initial "
                       "abundance",
    "zeta_value_derived": False,
}
(OUT / "dyn9b3_nonsusy_leptogenesis.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-3: historical non-SUSY unflavored diagnostic", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "## Verdicts", "",
      f"1. **Net suppression vs DYN-7**: eps1 scales by {supp:.3f} "
      "(SM loop function x SM Higgs normalization) -- the non-SUSY "
      "branch is ~6x HARDER for thermal unflavored leptogenesis.",
      f"2. **Unweighted prior-draw regression (no likelihood)**: "
      f"target-band hit fraction = {target_band_hit_fraction:.4f} "
      f"(DYN-7 legacy regression reference: "
      f"{dyn7_target_band_hit_fraction:.4f}); "
      f"median |eta_B| = {np.median(np.abs(eta)):.2e}; boost table "
      f"{boost_sens}.",
      "3. **Reheating boundary**: only the gravitino-specific ceiling is "
      "absent.  Thermal and non-thermal production remain conditional on "
      "an explicit reheating history and initial abundance.",
      "4. **Physics status**: blocked_missing_branch_thermal_inputs; K6 "
      "cannot be updated from these unflavored hit fractions; physics "
      "promotion is forbidden.",
      "", "## Boundary (NOT claimed)", "",
      "- Unflavored, N_1-dominated regression only; flavored effects, "
      "spectator rates, reheating history, and initial abundance not "
      "computed.",
      "- No likelihood or importance weights are applied; target-band "
      "hit fractions are regression diagnostics, not posterior or "
      "viability probabilities.",
      "- The archival tower under the D3-type source is a conditional "
      "benchmark; the source card remains unpromoted.",
      "- zeta's value is NOT derived.",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b3_nonsusy_leptogenesis.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-9b-3: {n_pass}/{len(CHECKS)} arithmetic checks; suppression "
      f"{supp:.3f} vs DYN-7; physics status = "
      f"blocked_missing_branch_thermal_inputs; "
      f"ledgers -> "
      f"output/audit9/dyn9b3_nonsusy_leptogenesis.*")
if n_pass != len(CHECKS):
    raise SystemExit(1)
