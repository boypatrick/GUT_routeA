#!/usr/bin/env python3
"""DYN-4b: the refreshed Majorana card and the UNCONDITIONAL zeta posterior.

DYN-4a showed that with the full benchmark M_R frozen except for zeta, the
zeta direction is orthogonal to the NuFit-6.0 theta_23 octant flip.  DYN-4b
therefore frees the ENTIRE Majorana sector: the inverse type-I seesaw with
the archival Dirac sector absorbs the refreshed data exactly, and the
questions become (i) where do zeta, the Veronese part, and the contact
fraction move, and (ii) what is the data-driven posterior of zeta once the
unknowns nobody measures (Majorana phases, m1, delta_CP) are marginalized.

Pipeline:
  1. Regression: the new code path reproduces the old benchmark card
     (zeta anchors, I/J invariants of the old Veronese part against the
     Audit-0 card digits).
  2. Refreshed CENTRAL card: NuFit-6.0 IC24 centrals (delta = 212 deg; the
     old Majorana-phase convention 0.35 pi / 1.10 pi and m1 = 1e-3 eV are
     kept for the central card and DISCLOSED as convention): M_R', zeta',
     contact fraction', I', J', heavy spectrum'; theta_23 = 0.470 absorbed
     exactly; movement diagnostics in units of the old sensitivity windows.
  3. UNCONDITIONAL zeta posterior: Monte Carlo over the NuFit-6.0
     likelihood (5 observables + delta_CP) x uniform Majorana phases x
     log-uniform m1 in [1e-4, 3e-2] eV; per-sample inverse seesaw and
     decomposition; circular statistics for arg zeta; contact-essentiality
     robustness P(contact fraction > thresholds); correlations with the
     nuisances.
  4. DYN-7 feed: the posterior distribution of the lightest heavy Majorana
     mass M_1 and P(M_1 > 1e9 GeV) (Davidson-Ibarra viability indicator).

Boundary: the Dirac kernels stay archival (kernel-level refit against
refreshed charged-fermion/CKM data = DYN-4c); normal ordering; v_u =
100 GeV archival convention; a fit/posterior, not a prediction -- the
boundary theorem stands.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

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


def cmat(raw) -> np.ndarray:
    return np.array([[c["re"] + 1j * c["im"] for c in row] for row in raw],
                    dtype=complex)


# ------------------------------------------------------------ machinery
card = json.loads((ARCH / "publication_closure_card" /
                   "publication_closure_card.json").read_text())
yuk = {k: cmat(v) for k, v in card["selected_row"]["Yukawa_fit"].items()}
Y_e, Y_nu = yuk["charged_lepton"], yuk["neutrino_dirac"]

OLD = {"m1_eV": 1.0e-3, "dm21": 7.42e-5, "dm31": 2.517e-3,
       "s12": 0.304, "s13": 0.0222, "s23": 0.573,
       "delta": 1.20 * math.pi, "a21": 0.35 * math.pi, "a31": 1.10 * math.pi}
IC24 = {"s12": (0.308, 0.0115), "s13": (0.02215, 0.00057),
        "s23": (0.470, 0.015), "dm21": (7.49e-5, 0.19e-5),
        "dm31": (2.513e-3, 0.020e-3),
        "delta_deg": (212.0, 33.5)}


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
mD_eV = Y_nu * 100.0e9
c_hat = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], dtype=complex) / math.sqrt(3)


def build_card(s12, s13, s23, dm21, dm31, delta, a21, a31, m1):
    m_diag = np.array([m1, math.sqrt(m1**2 + dm21), math.sqrt(m1**2 + dm31)])
    u_nu = U_e @ standard_pmns(s12, s13, s23, delta, a21, a31)
    m_light = u_nu.conjugate() @ np.diag(m_diag) @ u_nu.conjugate().T
    MR_eV = -(mD_eV.T @ np.linalg.inv(m_light) @ mD_eV)
    MR_GeV = MR_eV / 1e9
    sing = np.linalg.svd(MR_GeV, compute_uv=False)
    M_star = float(sing[0])
    Mn = MR_GeV / M_star
    zeta = complex(np.vdot(c_hat, Mn))
    return {"MR_GeV": MR_GeV, "M_star": M_star, "Mn": Mn, "zeta": zeta,
            "M_V": Mn - zeta * c_hat,
            "contact_fraction": float(abs(zeta) / np.linalg.norm(Mn)),
            "heavy_GeV": [float(x) for x in sing]}


def quartic_invariants(mv):
    rt2 = math.sqrt(2.0)
    a, b, c = mv[0, 0], mv[0, 1] / rt2, mv[0, 2]
    d, e = mv[1, 2] / rt2, mv[2, 2]
    I = a * e - 4.0 * b * d + 3.0 * c * c
    J = a * c * e + 2.0 * b * c * d - a * d * d - b * b * e - c**3
    return complex(I), complex(J)


def takagi_angles(MR_eV_matrix):
    m_l = -(mD_eV @ np.linalg.inv(MR_eV_matrix) @ mD_eV.T)
    vals, vecs = np.linalg.eigh(m_l.conjugate().T @ m_l)
    order = np.argsort(vals)
    masses = np.sqrt(np.maximum(vals[order], 0))
    u = U_e.conjugate().T @ vecs[:, order]
    p2 = np.abs(u) ** 2
    s13 = float(p2[0, 2])
    return {"s12": float(p2[0, 1] / (1 - s13)), "s13": s13,
            "s23": float(p2[1, 2] / (1 - s13)),
            "dm21": float(masses[1] ** 2 - masses[0] ** 2),
            "dm31": float(masses[2] ** 2 - masses[0] ** 2)}


# ------------------------------------------------------------ section 1
print("== DYN-4b section 1: regression through the new code path ==")

old_card = build_card(OLD["s12"], OLD["s13"], OLD["s23"], OLD["dm21"],
                      OLD["dm31"], OLD["delta"], OLD["a21"], OLD["a31"],
                      OLD["m1_eV"])
check("old benchmark reproduces through the DYN-4b path: zeta anchor bit-level",
      abs(old_card["zeta"] - complex(0.1076472949, 0.0736514853)) < 1e-9,
      f"zeta = {old_card['zeta']:.12f}")
I_old, J_old = quartic_invariants(old_card["M_V"])
check("Audit-0 invariants REGRESS on the old Veronese part: "
      "I = 0.008692433865+0.018779423508i, J = -8.0628e-5+5.6043e-4i",
      abs(I_old - complex(0.008692433865, 0.018779423508)) < 1e-10
      and abs(J_old - complex(-8.062794523e-5, 5.604303361e-4)) < 1e-11,
      f"I = {I_old:.12f}, J = {J_old:.6e}")

# ------------------------------------------------------------ section 2
print("== DYN-4b section 2: refreshed central card (NuFit-6.0 IC24) ==")

NEW_DELTA = math.radians(IC24["delta_deg"][0])
new_card = build_card(IC24["s12"][0], IC24["s13"][0], IC24["s23"][0],
                      IC24["dm21"][0], IC24["dm31"][0], NEW_DELTA,
                      OLD["a21"], OLD["a31"], OLD["m1_eV"])
obs_new = takagi_angles(new_card["MR_GeV"] * 1e9)
check("theta_23 octant flip ABSORBED exactly: forward seesaw returns the "
      "IC24 targets to 1e-9 (s23 = 0.470)",
      abs(obs_new["s23"] - 0.470) < 1e-9 and abs(obs_new["s12"] - 0.308) < 1e-9
      and abs(obs_new["dm31"] / 2.513e-3 - 1) < 1e-9)
zeta_new = new_card["zeta"]
I_new, J_new = quartic_invariants(new_card["M_V"])
move_phase = abs(np.angle(zeta_new) - np.angle(old_card["zeta"]))
move_scale = abs(abs(zeta_new) / abs(old_card["zeta"]) - 1)
check("refreshed card published: zeta', contact fraction', I', J', M_*'",
      True,
      f"zeta' = {zeta_new:.6f} (|.| {abs(zeta_new):.6f}, arg {np.angle(zeta_new):.6f}), "
      f"cf' = {new_card['contact_fraction']:.6f}, M_*' = {new_card['M_star']:.4e} GeV")
check("movement diagnostics vs the old sensitivity windows DISCLOSED: the "
      "data refresh moves zeta by orders of magnitude more than the old "
      "loose windows (the old anchor digits are decisively superseded)",
      move_phase > 100 * 5.424119e-5 or move_scale > 100 * 3.208798e-5,
      f"phase moved {move_phase:.4f} rad = {move_phase / 5.424119e-5:.0f} loose "
      f"windows; scale moved {move_scale:.4f} = {move_scale / 3.208798e-5:.0f} windows")
check("refreshed invariants I', J' published (Audit-0.5 refresh input)", True,
      f"I' = {I_new:.9f}, J' = {J_new:.6e}")

# Audit-0.5-style finite phase-transfer table vs the refreshed arg zeta
cand = {"arg I'": np.angle(I_new), "2 arg I'": 2 * np.angle(I_new),
        "arg J'": np.angle(J_new), "arg I' + arg J'": np.angle(I_new) + np.angle(J_new)}
best_name, best_miss = min(
    ((k, min(abs((v - np.angle(zeta_new) + math.pi) % (2 * math.pi) - math.pi),
             abs((v - np.angle(zeta_new)) % math.pi)))
     for k, v in cand.items()), key=lambda t: t[1])
check("Audit-0.5 phase-transfer re-test on the refreshed card: no candidate "
      "lands near arg zeta' (no-hit verdict persists)", best_miss > 0.05,
      f"closest = {best_name}, miss = {best_miss:.4f} rad")

# ------------------------------------------------------------ section 3
print("== DYN-4b section 3: UNCONDITIONAL zeta posterior (MC marginalization) ==")

rng = np.random.default_rng(20260705)
N = 4000
samples = []
for _ in range(N):
    try:
        cardi = build_card(
            rng.normal(*IC24["s12"]), rng.normal(*IC24["s13"]),
            rng.normal(*IC24["s23"]), rng.normal(*IC24["dm21"]),
            rng.normal(*IC24["dm31"]),
            math.radians(rng.normal(*IC24["delta_deg"])),
            rng.uniform(0, 2 * math.pi), rng.uniform(0, 2 * math.pi),
            10 ** rng.uniform(math.log10(1e-4), math.log10(3e-2)))
        samples.append((abs(cardi["zeta"]), np.angle(cardi["zeta"]),
                        cardi["contact_fraction"], cardi["heavy_GeV"][-1],
                        cardi["M_star"]))
    except np.linalg.LinAlgError:
        continue
S = np.array(samples)
check(f"Monte Carlo completed on >= {int(0.95 * N)}/{N} draws", len(S) >= 0.95 * N,
      f"{len(S)} samples")

q = lambda col, p: float(np.percentile(S[:, col], p))
mean_dir = np.angle(np.mean(np.exp(1j * S[:, 1])))
dev = (S[:, 1] - mean_dir + math.pi) % (2 * math.pi) - math.pi
arg_q = lambda p: float(mean_dir + np.percentile(dev, p))
post = {
    "abs_zeta_percentiles_2p5_16_50_84_97p5": [q(0, 2.5), q(0, 16), q(0, 50),
                                               q(0, 84), q(0, 97.5)],
    "arg_zeta_circular_mean": float(mean_dir),
    "arg_zeta_percentiles_2p5_16_50_84_97p5": [arg_q(2.5), arg_q(16), arg_q(50),
                                               arg_q(84), arg_q(97.5)],
    "contact_fraction_percentiles_16_50_84": [q(2, 16), q(2, 50), q(2, 84)],
    "old_benchmark_abs_zeta": abs(old_card["zeta"]),
    "old_benchmark_arg_zeta": float(np.angle(old_card["zeta"])),
}
check("UNCONDITIONAL zeta posterior published (percentiles of |zeta|, "
      "circular arg zeta, contact fraction)", True,
      f"|zeta| = {post['abs_zeta_percentiles_2p5_16_50_84_97p5'][2]:.4f} "
      f"[{post['abs_zeta_percentiles_2p5_16_50_84_97p5'][1]:.4f}, "
      f"{post['abs_zeta_percentiles_2p5_16_50_84_97p5'][3]:.4f}], "
      f"arg = {post['arg_zeta_circular_mean']:.4f} "
      f"[{post['arg_zeta_percentiles_2p5_16_50_84_97p5'][1]:.4f}, "
      f"{post['arg_zeta_percentiles_2p5_16_50_84_97p5'][3]:.4f}] rad")

ess = {f"P(cf > {t})": float(np.mean(S[:, 2] > t)) for t in (0.01, 0.05, 0.13)}
check("contact-essentiality robustness: the Veronese-only branch stays "
      "excluded across the posterior (P(contact fraction > 0.01) reported)",
      ess["P(cf > 0.01)"] > 0.95, f"{ess}")

in_band = float(np.mean((S[:, 0] > 0.9 * abs(old_card["zeta"]))
                        & (S[:, 0] < 1.1 * abs(old_card["zeta"]))))
check("posterior position of the OLD benchmark zeta DISCLOSED "
      "(fraction of posterior within +/-10% of the old |zeta|)", True,
      f"{in_band:.3f} of posterior mass; old anchor is a point in a broad "
      "posterior, not a prediction")

# ------------------------------------------------------------ section 4
print("== DYN-4b section 4: DYN-7 feed -- lightest heavy Majorana mass ==")

M1_stats = {"log10_M1_GeV_16_50_84": [math.log10(q(3, 16)),
                                      math.log10(q(3, 50)),
                                      math.log10(q(3, 84))],
            "P_M1_above_1e9_GeV": float(np.mean(S[:, 3] > 1e9)),
            "central_card_heavy_GeV": new_card["heavy_GeV"]}
check("M_1 posterior published; Davidson-Ibarra viability indicator "
      "P(M_1 > 1e9 GeV) reported for DYN-7", True,
      f"log10 M_1 = {M1_stats['log10_M1_GeV_16_50_84'][1]:.2f} "
      f"[{M1_stats['log10_M1_GeV_16_50_84'][0]:.2f}, "
      f"{M1_stats['log10_M1_GeV_16_50_84'][2]:.2f}]; "
      f"P(M_1 > 1e9) = {M1_stats['P_M1_above_1e9_GeV']:.3f}")

# correlations with nuisances (rerun small MC recording nuisances)
rng2 = np.random.default_rng(7)
rows = []
for _ in range(1500):
    m1 = 10 ** rng2.uniform(-4, math.log10(3e-2))
    a21, a31 = rng2.uniform(0, 2 * math.pi), rng2.uniform(0, 2 * math.pi)
    try:
        cardi = build_card(IC24["s12"][0], IC24["s13"][0], IC24["s23"][0],
                           IC24["dm21"][0], IC24["dm31"][0], NEW_DELTA,
                           a21, a31, m1)
        rows.append((math.log10(m1), a21, a31, abs(cardi["zeta"]),
                     cardi["contact_fraction"]))
    except np.linalg.LinAlgError:
        continue
R = np.array(rows)
corr_m1 = float(np.corrcoef(R[:, 0], R[:, 3])[0, 1])
check("nuisance attribution: correlation of |zeta| with log10 m1 and the "
      "Majorana phases DISCLOSED (which unknowns drive the posterior width)",
      True,
      f"corr(|zeta|, log10 m1) = {corr_m1:+.2f}, "
      f"corr(|zeta|, a21) = {float(np.corrcoef(R[:, 1], R[:, 3])[0, 1]):+.2f}, "
      f"corr(|zeta|, a31) = {float(np.corrcoef(R[:, 2], R[:, 3])[0, 1]):+.2f}")

# ------------------------------------------------------------ ledger
npass = sum(1 for _, ok in CHECKS if ok)
report = {
    "audit": "audit1_dyn4b_refreshed_card_unconditional_zeta",
    "dyn_item": "DYN-4b",
    "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "checks_total": len(CHECKS), "checks_passed": npass,
    "all_pass": npass == len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": {
        "design": "DYN-4a proved zeta (M_V frozen) is orthogonal to the "
                  "theta_23 flip; DYN-4b frees the full Majorana sector: "
                  "the inverse seesaw absorbs the refreshed data exactly "
                  "and the posterior of the DECOMPOSITION (zeta, contact "
                  "fraction) is computed by marginalizing the unmeasured "
                  "inputs",
        "priors": "5 observables + delta_CP ~ NuFit-6.0 IC24 gaussians; "
                  "Majorana phases uniform [0, 2pi); m1 log-uniform "
                  "[1e-4, 3e-2] eV (cosmology-safe); NO",
        "central_card_convention": "delta = 212 deg, old phase convention "
                                   "0.35pi/1.10pi, m1 = 1e-3 eV -- disclosed",
        "invariants": "classical binary-quartic I, J in the audit-0 "
                      "(x^2, sqrt2 xy, y^2) convention; old card regressed "
                      "to the paper digits",
    },
    "refreshed_central_card": {
        "zeta": [zeta_new.real, zeta_new.imag],
        "abs_zeta": abs(zeta_new), "arg_zeta": float(np.angle(zeta_new)),
        "contact_fraction": new_card["contact_fraction"],
        "M_star_GeV": new_card["M_star"],
        "heavy_GeV": new_card["heavy_GeV"],
        "I": [I_new.real, I_new.imag], "J": [J_new.real, J_new.imag],
        "movement_vs_old_windows": {
            "phase_rad": move_phase,
            "phase_in_loose_windows": move_phase / 5.424119e-5,
            "scale_rel": move_scale,
            "scale_in_loose_windows": move_scale / 3.208798e-5},
        "audit05_retest_closest_miss_rad": [best_name, best_miss],
    },
    "unconditional_zeta_posterior": post,
    "contact_essentiality": ess,
    "old_anchor_posterior_band_fraction": in_band,
    "dyn7_feed": M1_stats,
    # negative-boundary flags
    "fit_not_prediction": True,
    "dirac_kernels_archival_refit_deferred_to_DYN4c": True,
    "normal_ordering_assumed": True,
    "central_card_phase_convention_disclosed": True,
    "zeta_value_derived": False,
}

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "dyn4b_unconditional_zeta.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n")
(OUT / "dyn4b_unconditional_zeta.md").write_text(f"""# DYN-4b: Refreshed Majorana Card and the Unconditional zeta Posterior

`{npass}/{len(CHECKS)}` checks passed.

## Design

DYN-4a proved the zeta direction (M_V frozen) is orthogonal to the NuFit-6.0
theta_23 octant flip.  Here the FULL Majorana sector is freed: the inverse
type-I seesaw with the archival Dirac sector absorbs the refreshed data
exactly, and the posterior of the covariant decomposition is computed.

## Refreshed central card (IC24 centrals; conventions disclosed)

zeta' = {zeta_new.real:.6f} + {zeta_new.imag:.6f} i
(|zeta'| = {abs(zeta_new):.6f}, arg = {np.angle(zeta_new):.6f} rad),
contact fraction {new_card['contact_fraction']:.6f},
M_*' = {new_card['M_star']:.4e} GeV,
I' = {I_new:.9f}, J' = {J_new:.6e}.
The refresh moves arg zeta by {move_phase / 5.424119e-5:.0f} old loose
windows -- the printed benchmark digits are decisively superseded by data,
exactly as the boundary theorem says data (not principles) must do.
Audit-0.5 phase-transfer re-test: still no hit
(closest {best_name}, miss {best_miss:.3f} rad).

## Unconditional zeta posterior (NuFit-6.0 x phases x m1)

|zeta|: {post['abs_zeta_percentiles_2p5_16_50_84_97p5'][2]:.4f}
[{post['abs_zeta_percentiles_2p5_16_50_84_97p5'][1]:.4f},
{post['abs_zeta_percentiles_2p5_16_50_84_97p5'][3]:.4f}] (68%),
[{post['abs_zeta_percentiles_2p5_16_50_84_97p5'][0]:.4f},
{post['abs_zeta_percentiles_2p5_16_50_84_97p5'][4]:.4f}] (95%).
arg zeta: {post['arg_zeta_circular_mean']:.4f}
[{post['arg_zeta_percentiles_2p5_16_50_84_97p5'][1]:.4f},
{post['arg_zeta_percentiles_2p5_16_50_84_97p5'][3]:.4f}] rad (68%, circular).
Contact essentiality: {ess}.
Old anchor band (+/-10% in |zeta|): {in_band:.3f} of posterior mass.

## DYN-7 feed

log10 M_1/GeV = {M1_stats['log10_M1_GeV_16_50_84'][1]:.2f}
[{M1_stats['log10_M1_GeV_16_50_84'][0]:.2f},
{M1_stats['log10_M1_GeV_16_50_84'][2]:.2f}];
P(M_1 > 1e9 GeV) = {M1_stats['P_M1_above_1e9_GeV']:.3f}.

## Boundary

Dirac kernels archival (DYN-4c open); NO assumed; central-card phase
convention disclosed; this is a posterior, not a prediction.
""")
print(f"Wrote {OUT / 'dyn4b_unconditional_zeta.json'} (+ .md)")
print(f"DYN-4b: {npass}/{len(CHECKS)} checks passed; refreshed card zeta' = "
      f"{zeta_new:.4f}; posterior and DYN-7 feed published.")
if npass != len(CHECKS):
    raise SystemExit(1)
