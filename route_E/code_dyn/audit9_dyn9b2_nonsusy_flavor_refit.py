#!/usr/bin/env python3
"""DYN-9b-2: fixed-kernel non-SUSY rescaling and zeta invariance.

The perturbativity gate (DYN-9 patch) showed the archival Majorana
tower cannot be sourced renormalizably on the surviving chains
(f = M_R/v_R needs up to 379x more than 4 pi).  This audit does the
fixed-kernel work that follows from the LIGHT-neutrino side, without
pretending to a full non-SUSY Yukawa fit (disclosed):

  S1  Archival replay + fixed-kernel ceiling diagnostics.  The type-I seesaw
      with the oscillation data FIXES the required tower given the
      Dirac Yukawa scale: M_R ~ (Y v)^2 / m_nu.  The relations
      Y_nu = h - 3f and Y_u = h + f do not by themselves force
      Y_nu ~ Y_u (h=3f is an immediate counterexample).  The
      perturbative ceiling 4 pi M_I caps the allowed fixed-shape Dirac
      scale at y_max per chain.  We report both the actual archival-kernel
      suppression y_arch/y_max and, separately, the tension obtained only
      after adding a top-like third-generation ansatz y_nu,3 ~ y_top(M_X).

  S2  Type-II escape ESTIMATE (order of magnitude, flagged): the
      induced triplet vev v_L ~ lambda v^2 v_R / M_Delta^2 with the
      Delta_L-type block at the gauge scale (DYN-9b-1c) falls short of
      m_nu/f by many orders on both chains: type-II does not rescue
      the renormalizable source.

  S3  POSITIVE-REAL RESCALING THEOREM.  For y > 0, rescaling the Dirac
      kernel Y_nu -> y Y_nu gives M_R -> y^2 M_R and
      M_* -> y^2 M_*, so zeta, the contact direction, and the contact
      fraction are exactly invariant.  For complex y, instead,
      zeta -> (y^2/|y|^2) zeta: its phase is covariant, not invariant.
      Both statements are machine-verified.  The re-extracted M_* windows
      at v_R ~ M_I use the positive-real fixed-shape family.

  S4  CONDITIONAL LEPTOGENESIS SOURCE COMPARISON.  Under the rescaled
      renormalizable source, M_1 -> y^2 M_1 falls BELOW the
      Davidson-Ibarra floor on both chains (PS: ~1e8, LR: ~2e5 GeV);
      under a scale-decoupled (instanton-type, D3) source the archival
      M_1 ~ 2.8e10 remains above that order-estimate floor.  This is a
      fixed-kernel comparison, not a source-selection theorem: a global
      branch-local Yukawa fit and flavored thermal kinetics have not
      been performed.  The D3 coexistence property K8 (N_2, N_3 above
      the intermediate gauge scale) therefore remains an unpromoted
      conditional diagnostic.

Boundary: the Dirac KERNEL SHAPE is retained from the archival card
(conditional input; a full non-SUSY Yukawa fit is NOT performed); the
top-like comparison is an additional ansatz, not an SO(10) theorem; the
type-II and Davidson-Ibarra numbers are order estimates (flagged); zeta's
VALUE is NOT derived.
Ledgers -> output/audit9/dyn9b2_nonsusy_flavor_refit.{json,md}.
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
OUT = AUDIT_OUTPUT / "audit9"
CHECKS = []
DLOG = {}


def check(name, ok, detail=""):
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


# ------------------------------------------------------------- inputs
dyn4a = json.loads((ROOT / "output" / "audit1"
                    / "dyn4a_seesaw_zeta_posterior.json").read_text())
dyn4b = json.loads((ROOT / "output" / "audit1"
                    / "dyn4b_unconditional_zeta.json").read_text())
dyn9 = json.loads((ROOT / "output" / "audit9"
                   / "dyn9_nonsusy_intermediate.json").read_text())
d3 = json.loads((ROOT / "route_d" / "output"
                 / "d3_instanton_majorana_pricing.json").read_text())
card_path = ARCHIVE_OUTPUT / "publication_closure_card" \
    / "publication_closure_card.json"
card_sha = sha256(card_path)

print("== 9b-2 section 0: premise chain-of-custody ==")
check("archival closure card sha256 matches the DYN-4a provenance; the "
      "DYN-9 ledger declares the card slice-local (the premise this "
      "audit sharpens)",
      card_sha == dyn4a["provenance"]["closure_card"]["sha256"]
      and dyn9["archival_zeta_card_transplantable_to_nonsusy_chains"]
      is False)

# ------------------------------------------------------------ section 1
print("== 9b-2 section 1: archival replay + fixed-kernel ceiling ==")

card = json.loads(card_path.read_text())
yuk = {k: cmat(v) for k, v in card["selected_row"]["Yukawa_fit"].items()}
Y_e, Y_nu = yuk["charged_lepton"], yuk["neutrino_dirac"]
OLDT = {"m1_eV": 1.0e-3, "dm21": 7.42e-5, "dm31": 2.517e-3,
        "s12": 0.304, "s13": 0.0222, "s23": 0.573,
        "delta": 1.20 * math.pi, "a21": 0.35 * math.pi,
        "a31": 1.10 * math.pi}


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
m_diag = np.array([OLDT["m1_eV"],
                   math.sqrt(OLDT["m1_eV"] ** 2 + OLDT["dm21"]),
                   math.sqrt(OLDT["m1_eV"] ** 2 + OLDT["dm31"])])
U_nu = U_e @ standard_pmns(OLDT["s12"], OLDT["s13"], OLDT["s23"],
                           OLDT["delta"], OLDT["a21"], OLDT["a31"])
m_light = U_nu.conj() @ np.diag(m_diag) @ U_nu.conj().T


def tower(Ynu):
    mD = Ynu * 100.0e9                       # v_u = 100 GeV, in eV
    MR = -(mD.T @ np.linalg.inv(m_light) @ mD) / 1e9     # GeV
    sv = np.linalg.svd(MR, compute_uv=False)
    Mst = float(sv[0])
    c_hat = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]],
                     dtype=complex) / math.sqrt(3)
    z = complex(np.vdot(c_hat, MR / Mst))
    cf = float(np.linalg.norm(z * c_hat) / np.linalg.norm(MR / Mst))
    return MR, sv, Mst, z, cf


MR0, sv0, Mst0, zeta0, cf0 = tower(Y_nu)
check("archival replay reproduces the paper anchors (zeta, M_*, tower)",
      abs(zeta0 - complex(0.1076472949, 0.0736514853)) < 1e-9
      and abs(Mst0 - 3.9265e15) / 3.9265e15 < 1e-3,
      f"zeta = {zeta0.real:.6f}+{zeta0.imag:.6f}i, "
      f"M_* = {Mst0:.4e} GeV")

# y_top at M_X from the one-loop SM running (optional top-like reference).
# Here b1=41/10 is the SU(5)-normalised convention, and 0.462 is already
# g1=sqrt(5/3) g_Y at the electroweak matching point.  Multiplying 0.462 by
# sqrt(5/3) a second time double-normalises hypercharge and biases y_t(M_X).
g1, g2, g3, yt = 0.462, 0.648, 1.166, 0.936


def run_yt(lnratio, n=4000):
    global g1, g2, g3, yt
    h = lnratio / n
    a1, a2, a3, y = g1, g2, g3, yt
    for _ in range(n):
        b1 = 41 / 10 * a1 ** 3 / (16 * math.pi ** 2)
        b2 = -19 / 6 * a2 ** 3 / (16 * math.pi ** 2)
        b3 = -7 * a3 ** 3 / (16 * math.pi ** 2)
        by = y * (4.5 * y * y - 8 * a3 * a3 - 2.25 * a2 * a2
                  - 0.85 * a1 * a1) / (16 * math.pi ** 2)
        a1, a2, a3, y = a1 + h * b1, a2 + h * b2, a3 + h * b3, y + h * by
    return y


yt_MX = run_yt(math.log(2e16 / 173.0))
check("y_top(M_X) machine-derived from the one-loop SM running for the "
      "optional top-like third-generation benchmark",
      0.3 < yt_MX < 0.6, f"y_t(2e16) = {yt_MX:.3f}")

# The displayed minimal 10_H + 126bar_H relations do not imply a lock:
# h=3f gives Y_nu=h-3f=0 while Y_u=h+f=4f.  A top-like Y_nu therefore
# has to be declared as an additional benchmark ansatz.
f_counter = 1.0
h_counter = 3.0 * f_counter
ynu_counter = h_counter - 3.0 * f_counter
yu_counter = h_counter + f_counter
check("algebraic counterexample: Y_nu=h-3f and Y_u=h+f do not force "
      "Y_nu ~ Y_u; h=3f gives Y_nu=0 and Y_u=4f",
      ynu_counter == 0.0 and yu_counter == 4.0 * f_counter)

y_arch = float(np.linalg.svd(Y_nu, compute_uv=False)[0])
check("archival Dirac kernel scale recorded: sigma_max(Y_nu) happens to "
      "be top-like within an order-one factor in this card; this is a "
      "benchmark property, not a consequence of the SO(10) relations",
      0.2 < y_arch / yt_MX < 5.0,
      f"sigma_max(Y_nu) = {y_arch:.3f} = {y_arch/yt_MX:.2f} x y_t(M_X)")

FOUR_PI = 4 * math.pi
MI = {n: 10 ** v["log10_MI"] for n, v in dyn9["chains"].items()
      if n in ("PS", "G_LR")}
lock = {}
for n, mi in MI.items():
    MR_ceiling = FOUR_PI * mi
    y_max = y_arch * math.sqrt(MR_ceiling / Mst0)
    lock[n] = {
        "M_I": mi,
        "MR_ceiling_GeV": MR_ceiling,
        "y_max": y_max,
        "archival_kernel_suppression": y_arch / y_max,
        "top_like_ansatz_tension": yt_MX / y_max,
    }
check("FIXED-ARCHIVAL-KERNEL diagnostic: the perturbative 4 pi M_I "
      "ceiling requires a uniform norm suppression y_arch/y_max of about "
      "19 (PS) / 342 (G_LR)",
      lock["PS"]["archival_kernel_suppression"] > 10
      and lock["G_LR"]["archival_kernel_suppression"] > 300,
      "; ".join(f"{n}: y_max = {v['y_max']:.4f}, archival suppression = "
                f"{v['archival_kernel_suppression']:.1f}x"
                for n, v in lock.items()))
check("OPTIONAL TOP-LIKE ANSATZ diagnostic: if one additionally imposes "
      "y_nu,3 ~ y_t(M_X), the same ceilings correspond to tensions of "
      "about 10 (PS) / 169 (G_LR); this is not an SO(10) theorem",
      lock["PS"]["top_like_ansatz_tension"] > 5
      and lock["G_LR"]["top_like_ansatz_tension"] > 50,
      "; ".join(f"{n}: top-like tension = "
                f"{v['top_like_ansatz_tension']:.1f}x"
                for n, v in lock.items()))
DLOG["S1_lock"] = {
    "yt_MX": yt_MX,
    "y_arch": y_arch,
    "so10_relations_force_top_like_Ynu": False,
    "counterexample": {"h_over_f": 3.0, "Ynu_over_f": 0.0,
                       "Yu_over_f": 4.0},
    "chains": lock,
}

# ------------------------------------------------------------ section 2
print("== 9b-2 section 2: the type-II escape estimate ==")

v_EW = 174.0
m_nu3 = 0.05e-9                                   # GeV
typeII = {}
for n, mi in MI.items():
    M_Delta = 10 ** dyn9["chains"][n]["log10_MX"]  # gauge-scale block
    v_L = v_EW ** 2 * mi / M_Delta ** 2            # lambda ~ 1
    m_typeII = FOUR_PI * v_L                       # f <= 4 pi
    typeII[n] = {"v_L_GeV": v_L, "m_max_eV": m_typeII * 1e9,
                 "deficit_orders": math.log10(m_nu3 / m_typeII)}
check("TYPE-II ESCAPE FAILS (order estimate, flagged; lambda ~ 1 "
      "and f at 4 pi, i.e. GENEROUS): with the Delta_L-type block at "
      "the gauge scale (DYN-9b-1c) the induced triplet vev falls short "
      "of the atmospheric scale by 3.7 (PS) / 7.3 (G_LR) orders",
      all(v["deficit_orders"] > 3 for v in typeII.values()),
      "; ".join(f"{n}: max m_nu = {v['m_max_eV']:.1e} eV, deficit "
                f"{v['deficit_orders']:.1f} orders"
                for n, v in typeII.items()))
DLOG["S2_typeII"] = typeII

# ------------------------------------------------------------ section 3
print("== 9b-2 section 3: positive-real invariance / complex covariance ==")

inv_ok = True
worst = 0.0
for y in (0.5, 0.1, lock["PS"]["y_max"] / y_arch,
          lock["G_LR"]["y_max"] / y_arch):
    MRy, svy, Msty, zy, cfy = tower(y * Y_nu)
    dz = abs(zy - zeta0)
    dM = abs(Msty / Mst0 - y * y)
    dcf = abs(cfy - cf0)
    worst = max(worst, dz, dM, dcf)
    inv_ok &= dz < 1e-12 and dM < 1e-12 and dcf < 1e-12
check("POSITIVE-REAL RESCALING THEOREM (machine-verified): for y>0, "
      "Y_nu -> y Y_nu gives M_R -> y^2 M_R and M_* -> y^2 M_*, so "
      "zeta, the projective contact direction, and contact fraction are "
      "exactly invariant within this fixed-kernel one-parameter family",
      inv_ok, f"worst deviation over four rescales: {worst:.2e}")

theta = 0.37
y_complex = 0.4 * np.exp(1j * theta)
_, _, Mst_complex, zeta_complex, cf_complex = tower(y_complex * Y_nu)
phase_factor = y_complex ** 2 / abs(y_complex) ** 2
complex_covariance_error = max(
    abs(Mst_complex / Mst0 - abs(y_complex) ** 2),
    abs(zeta_complex - phase_factor * zeta0),
    abs(cf_complex - cf0),
)
check("COMPLEX-RESCALING BOUNDARY (machine-verified): for complex y, "
      "M_* -> |y|^2 M_* and zeta -> (y^2/|y|^2) zeta; the contact "
      "fraction is invariant but the complex zeta phase is not",
      complex_covariance_error < 1e-12
      and abs(zeta_complex - zeta0) > 1e-3,
      f"theta={theta:.2f}, covariance error={complex_covariance_error:.2e}")

Mstar_win = {n: {"M_star_at_ceiling_GeV": (v["y_max"] / y_arch) ** 2 * Mst0}
             for n, v in lock.items()}
check("re-extracted M_* windows at v_R ~ M_I recorded per chain "
      "(renormalizable-source scenario at the perturbative ceiling): "
      "M_*' = 4 pi M_I by construction",
      all(abs(Mstar_win[n]["M_star_at_ceiling_GeV"]
              / (FOUR_PI * MI[n]) - 1) < 1e-9 for n in MI),
      "; ".join(f"{n}: M_*' <= {Mstar_win[n]['M_star_at_ceiling_GeV']:.2e}"
                for n in MI))
check("contact essentiality carries over to the positive-real rescaled "
      "card (cf invariant; the DYN-4b posterior statement "
      "P(cf > 0.01) = 1.000 remains conditional on the Dirac kernel "
      "SHAPE only, which the rescale does not touch)",
      dyn4b["contact_essentiality"]["P(cf > 0.01)"] == 1.0
      and abs(cf0 - 0.1304275166688152) < 1e-9,
      f"cf = {cf0:.10f} (invariant)")
DLOG["S3_invariance"] = {
    "domain_of_exact_zeta_invariance": "uniform positive-real y",
    "zeta": [zeta0.real, zeta0.imag],
    "contact_fraction": cf0,
    "positive_real_worst_error": worst,
    "complex_rescaling": {
        "theta_rad": theta,
        "phase_factor": [phase_factor.real, phase_factor.imag],
        "covariance_error": complex_covariance_error,
        "zeta_phase_invariant": False,
    },
    "M_star_windows": Mstar_win,
}

# ------------------------------------------------------------ section 4
print("== 9b-2 section 4: conditional source comparison (unpromoted) ==")

M1_arch = 10 ** dyn4b["dyn7_feed"]["log10_M1_GeV_16_50_84"][1]
DI_FLOOR = 5.0e8                        # Davidson-Ibarra, order estimate
sel = {}
for n, v in lock.items():
    M1_resc = (v["y_max"] / y_arch) ** 2 * M1_arch
    sel[n] = {"M1_rescaled_GeV": M1_resc,
              "below_DI_floor": bool(M1_resc < DI_FLOOR)}
check("CONDITIONAL SOURCE COMPARISON: under the fixed-kernel rescaled "
      "source M_1 falls BELOW the Davidson-Ibarra floor on both chains "
      "(PS ~ 7e7, G_LR ~ 2e5 GeV); under a scale-decoupled "
      "(instanton-type) benchmark the archival M_1 ~ 2.8e10 stays "
      "above the order-estimate floor.  This does not select a source "
      "without a global flavor fit and flavored thermal kinetics",
      all(v["below_DI_floor"] for v in sel.values())
      and M1_arch > DI_FLOOR,
      "; ".join(f"{n}: M_1' = {v['M1_rescaled_GeV']:.1e}"
                for n, v in sel.items())
      + f"; archival M_1 = {M1_arch:.1e}")

coex = d3["price_card"]["beyond_zeta_consequence"]
check("K8 NON-PROMOTION: the fixed archival tower used by the D3 "
      "benchmark has N_2 and N_3 above the intermediate gauge scale; "
      "the instanton mechanism permits but does not generically require "
      "this ordering, and the D3 card is explicitly unpromoted",
      d3["all_pass"] and "heavier than the intermediate" in coex
      and not d3["promoted_to_paper"],
      "global flavor fit and flavored/density-matrix kinetics remain open")
DLOG["S4_conditional_comparison"] = {
    "M1_archival": M1_arch, "DI_floor": DI_FLOOR,
    "rescaled": sel, "source_selected": False,
    "reason_not_selected": (
        "fixed Dirac-kernel shape, order-estimate DI floor, and missing "
        "flavored thermal kinetics")}

# ------------------------------------------------------------- ledger
n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-9b-2 fixed-kernel non-SUSY rescaling audit",
    "dyn_item": "DYN-9b-2",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    "provenance": {"closure_card_sha256": card_sha,
                   "dyn4a": "output/audit1/dyn4a_seesaw_zeta_posterior"
                            ".json",
                   "dyn4b": "output/audit1/dyn4b_unconditional_zeta"
                            ".json",
                   "dyn9": "output/audit9/dyn9_nonsusy_intermediate"
                           ".json",
                   "d3": "route_d/output/d3_instanton_majorana_pricing"
                         ".json"},
    # negative-boundary flags
    "full_nonsusy_yukawa_fit_performed": False,
    "dirac_kernel_shape_retained_conditional": True,
    "so10_relations_force_top_like_Ynu": False,
    "top_like_Ynu_is_additional_ansatz": True,
    "fixed_kernel_ceiling_quantified_not_refit": True,
    "exact_zeta_invariance_domain": "uniform positive-real rescaling",
    "complex_rescaling_changes_zeta_phase": True,
    "typeII_and_DI_are_order_estimates": True,
    "zeta_value_derived": False,
    "mechanical_status": "pass" if n_pass == len(CHECKS) else "fail",
    "physics_status": "preliminary_missing_global_fit_and_flavored_kinetics",
    "physics_promotion_allowed": False,
    "blockers": [
        "full branch-local non-SUSY charged-fermion and neutrino Yukawa fit",
        "chain-specific threshold matching and intermediate-scale RGEs",
        "flavored or density-matrix leptogenesis with thermal initial data",
        "valid UV source construction replacing unpromoted D3 pricing",
    ],
}
(OUT / "dyn9b2_nonsusy_flavor_refit.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-2: fixed-kernel non-SUSY rescaling audit", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "## Verdicts", "",
      f"1. **Fixed-kernel ceiling**: y_t(M_X) = {yt_MX:.3f} "
      "is an optional top-like reference, not an SO(10) theorem.  The "
      "perturbative ceiling caps the nu-Dirac "
      f"scale at y_max = {lock['PS']['y_max']:.4f} (PS) / "
      f"{lock['G_LR']['y_max']:.4f} (G_LR).  The actual archival-kernel "
      f"suppressions are {lock['PS']['archival_kernel_suppression']:.1f}x / "
      f"{lock['G_LR']['archival_kernel_suppression']:.1f}x; the additional "
      f"top-like ansatz gives {lock['PS']['top_like_ansatz_tension']:.1f}x / "
      f"{lock['G_LR']['top_like_ansatz_tension']:.1f}x.",
      "2. **Type-II fails** by "
      f"{typeII['PS']['deficit_orders']:.1f} / "
      f"{typeII['G_LR']['deficit_orders']:.1f} orders (PS / G_LR).",
      "3. **Positive-real rescaling theorem**: for uniform y>0, zeta, "
      "projective contact direction, and contact fraction are exactly "
      "invariant while M_* moves (M_*' = 4 pi M_I at the ceiling: "
      f"{Mstar_win['PS']['M_star_at_ceiling_GeV']:.2e} PS, "
      f"{Mstar_win['G_LR']['M_star_at_ceiling_GeV']:.2e} G_LR).  For "
      "complex y, zeta acquires the phase y^2/|y|^2.",
      "4. **Conditional source comparison (not a selection)**: "
      "rescaled M_1 falls below the Davidson-Ibarra floor on both "
      "chains, while the archival M_1 lies above the adopted floor in "
      "the instanton-type benchmark.  K8 remains unpromoted because "
      "the comparison retains a fixed Dirac kernel and omits flavored "
      "thermal kinetics.",
      "", "## Boundary (NOT claimed)", "",
      "- No full non-SUSY Yukawa fit; the Dirac kernel shape is a "
      "retained conditional input.",
      "- Y_nu=h-3f and Y_u=h+f do not force a top-like Y_nu; the "
      "top-like comparison is an additional benchmark ansatz.",
      "- Exact complex-zeta invariance is restricted to uniform "
      "positive-real rescaling; complex rescaling is phase-covariant.",
      "- Type-II and Davidson-Ibarra numbers are order estimates.",
      "- zeta's value is NOT derived (boundary theorem intact).",
      "- Physics promotion is blocked pending the global fit, threshold "
      "matching, flavored kinetics, and a valid UV source construction.",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b2_nonsusy_flavor_refit.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-9b-2: {n_pass}/{len(CHECKS)} checks; archival suppression "
      f"{lock['PS']['archival_kernel_suppression']:.1f}x/"
      f"{lock['G_LR']['archival_kernel_suppression']:.1f}x; zeta invariant "
      "for uniform positive-real rescaling, M_* "
      f"branch-local; source comparison remains unpromoted; "
      f"ledgers -> output/audit9/dyn9b2_nonsusy_flavor_refit.*")

if n_pass != len(CHECKS):
    raise SystemExit(1)
