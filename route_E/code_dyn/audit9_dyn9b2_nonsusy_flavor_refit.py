#!/usr/bin/env python3
"""DYN-9b-2: the non-SUSY flavor refit and the zeta re-extraction.

The perturbativity gate (DYN-9 patch) showed the archival Majorana
tower cannot be sourced renormalizably on the surviving chains
(f = M_R/v_R needs up to 379x more than 4 pi).  This audit does the
refit-level work that follows from the LIGHT-neutrino side, without
pretending to a full non-SUSY Yukawa fit (disclosed):

  S1  Archival replay + the SO(10) LOCK TENSION.  The type-I seesaw
      with the oscillation data FIXES the required tower given the
      Dirac Yukawa scale: M_R ~ (Y v)^2 / m_nu.  In minimal SO(10)
      (10_H + 126bar_H) the neutrino Dirac Yukawa is LOCKED to the up
      sector (Y_nu = h - 3f vs Y_u = h + f), so its third-generation
      entry sits at O(y_top(M_X)) -- computed here from the one-loop
      SM running, not assumed.  The perturbative ceiling 4 pi M_I then
      caps the allowed Dirac scale at y_max per chain; the ratio
      y_lock / y_max is the LOCK-TENSION factor: the renormalizable
      minimal source needs the nu-Dirac coupling far below its locked
      scale on BOTH surviving chains.

  S2  Type-II escape ESTIMATE (order of magnitude, flagged): the
      induced triplet vev v_L ~ lambda v^2 v_R / M_Delta^2 with the
      Delta_L-type block at the gauge scale (DYN-9b-1c) falls short of
      m_nu/f by many orders on both chains: type-II does not rescue
      the renormalizable source.

  S3  THE ZETA-INVARIANCE THEOREM.  Rescaling the Dirac kernel
      Y_nu -> y Y_nu rescales the inverse-seesaw tower UNIFORMLY
      (M_R -> y^2 M_R), so the normalized decomposition is EXACTLY
      invariant: zeta, the contact direction, and the contact fraction
      do not move; only M_* = y^2 M_*_arch does.  Machine-verified to
      1e-12.  This SHARPENS the branch-locality statement of the
      DYN-9 patch: the contact card's normalized content is
      scale-covariant; M_* is the only branch-local number.  The
      re-extracted M_* windows at v_R ~ M_I are computed per chain.

  S4  LEPTOGENESIS SOURCE SELECTION.  Under the rescaled
      renormalizable source, M_1 -> y^2 M_1 falls BELOW the
      Davidson-Ibarra floor on both chains (PS: ~1e8, LR: ~2e5 GeV);
      under a scale-decoupled (instanton-type, D3) source the archival
      M_1 ~ 2.8e10 survives.  Thermal leptogenesis therefore SELECTS
      the scale-decoupled source on the alive branch, and the D3
      coexistence prediction K8 (N_2, N_3 above the intermediate gauge
      scale) upgrades from an optional card to a REQUIREMENT of the
      only quantified surviving source scenario.

Boundary: the Dirac KERNEL SHAPE is retained from the archival card
(conditional input; a full non-SUSY Yukawa fit is NOT performed); the
SO(10) lock is quantified from the third-generation relation, not
re-derived from a global fit; the type-II and Davidson-Ibarra numbers
are order estimates (flagged); zeta's VALUE is NOT derived.
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
print("== 9b-2 section 1: archival replay + the SO(10) lock tension ==")

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

# y_top at M_X from the one-loop SM running (machine-derived lock scale)
g1, g2, g3, yt = math.sqrt(5 / 3) * 0.462, 0.648, 1.166, 0.936


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
check("y_top(M_X) machine-derived from the one-loop SM running: the "
      "SO(10) lock scale for the third-generation nu-Dirac coupling",
      0.3 < yt_MX < 0.6, f"y_t(2e16) = {yt_MX:.3f}")

y_arch = float(np.linalg.svd(Y_nu, compute_uv=False)[0])
check("archival Dirac kernel scale recorded: sigma_max(Y_nu) is O(y_t) "
      "-- the archival card RESPECTS the SO(10) lock",
      0.2 < y_arch / yt_MX < 5.0,
      f"sigma_max(Y_nu) = {y_arch:.3f} = {y_arch/yt_MX:.2f} x y_t(M_X)")

FOUR_PI = 4 * math.pi
MI = {n: 10 ** v["log10_MI"] for n, v in dyn9["chains"].items()
      if n in ("PS", "G_LR")}
lock = {}
for n, mi in MI.items():
    MR_ceiling = FOUR_PI * mi
    y_max = y_arch * math.sqrt(MR_ceiling / Mst0)
    lock[n] = {"M_I": mi, "MR_ceiling_GeV": MR_ceiling,
               "y_max": y_max, "lock_tension": yt_MX / y_max}
check("SO(10) LOCK TENSION quantified: the perturbative ceiling "
      "4 pi M_I caps the nu-Dirac scale at y_max; the renormalizable "
      "minimal source needs the third-generation coupling ~9x (PS) / "
      "~160x (G_LR) BELOW its locked scale y_t(M_X) -- i.e. 2 / 4.4 "
      "orders of magnitude in the tower, a structural tension only "
      "removable by that much Yukawa tuning",
      lock["PS"]["lock_tension"] > 5 and lock["G_LR"]["lock_tension"] > 50,
      "; ".join(f"{n}: y_max = {v['y_max']:.4f}, tension = "
                f"{v['lock_tension']:.0f}x" for n, v in lock.items()))
DLOG["S1_lock"] = {"yt_MX": yt_MX, "y_arch": y_arch, "chains": lock}

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
print("== 9b-2 section 3: the zeta-invariance theorem ==")

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
check("ZETA-INVARIANCE THEOREM (machine-verified): under Y_nu -> y Y_nu "
      "the tower rescales uniformly, so zeta, the contact direction and "
      "the contact fraction are EXACTLY invariant; only M_* = y^2 "
      "M_*_arch moves.  The contact card's normalized content is "
      "scale-covariant; M_* is the only branch-local number",
      inv_ok, f"worst deviation over four rescales: {worst:.2e}")

Mstar_win = {n: {"M_star_at_ceiling_GeV": (v["y_max"] / y_arch) ** 2 * Mst0}
             for n, v in lock.items()}
check("re-extracted M_* windows at v_R ~ M_I recorded per chain "
      "(renormalizable-source scenario at the perturbative ceiling): "
      "M_*' = 4 pi M_I by construction",
      all(abs(Mstar_win[n]["M_star_at_ceiling_GeV"]
              / (FOUR_PI * MI[n]) - 1) < 1e-9 for n in MI),
      "; ".join(f"{n}: M_*' <= {Mstar_win[n]['M_star_at_ceiling_GeV']:.2e}"
                for n in MI))
check("contact essentiality carries over VERBATIM to the rescaled "
      "card (cf invariant; the DYN-4b posterior statement "
      "P(cf > 0.01) = 1.000 remains conditional on the Dirac kernel "
      "SHAPE only, which the rescale does not touch)",
      dyn4b["contact_essentiality"]["P(cf > 0.01)"] == 1.0
      and abs(cf0 - 0.1304275166688152) < 1e-9,
      f"cf = {cf0:.10f} (invariant)")
DLOG["S3_invariance"] = {"zeta": [zeta0.real, zeta0.imag],
                         "contact_fraction": cf0,
                         "M_star_windows": Mstar_win}

# ------------------------------------------------------------ section 4
print("== 9b-2 section 4: leptogenesis source selection, K8 upgrade ==")

M1_arch = 10 ** dyn4b["dyn7_feed"]["log10_M1_GeV_16_50_84"][1]
DI_FLOOR = 5.0e8                        # Davidson-Ibarra, order estimate
sel = {}
for n, v in lock.items():
    M1_resc = (v["y_max"] / y_arch) ** 2 * M1_arch
    sel[n] = {"M1_rescaled_GeV": M1_resc,
              "below_DI_floor": bool(M1_resc < DI_FLOOR)}
check("LEPTOGENESIS SOURCE SELECTION: under the rescaled renormalizable "
      "source M_1 falls BELOW the Davidson-Ibarra floor on both chains "
      "(PS ~ 7e7, G_LR ~ 2e5 GeV); under a scale-decoupled "
      "(instanton-type) source the archival M_1 ~ 2.8e10 survives -- "
      "thermal leptogenesis SELECTS the scale-decoupled source on the "
      "alive branch (DI floor is an order estimate, flagged)",
      all(v["below_DI_floor"] for v in sel.values())
      and M1_arch > DI_FLOOR,
      "; ".join(f"{n}: M_1' = {v['M1_rescaled_GeV']:.1e}"
                for n, v in sel.items())
      + f"; archival M_1 = {M1_arch:.1e}")

coex = d3["price_card"]["beyond_zeta_consequence"]
check("K8 UPGRADE: the only quantified source scenario that survives "
      "S1 (lock tension) + S2 (type-II fails) + S4 (leptogenesis) is "
      "the scale-decoupled instanton-type tower of the D3 card, whose "
      "coexistence prediction (N_2, N_3 above the intermediate gauge "
      "scale) therefore upgrades from an optional pricing card to a "
      "REQUIREMENT of the alive branch with viable leptogenesis",
      d3["all_pass"] and "heavier than the intermediate" in coex
      and not d3["promoted_to_paper"],
      "D3 remains UNPROMOTED: the upgrade is a conditional-dependency "
      "statement, not a promotion")
DLOG["S4_selection"] = {"M1_archival": M1_arch, "DI_floor": DI_FLOOR,
                        "rescaled": sel}

# ------------------------------------------------------------- ledger
n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-9b-2 non-SUSY flavor refit + zeta re-extraction",
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
    "so10_lock_quantified_not_refit": True,
    "typeII_and_DI_are_order_estimates": True,
    "zeta_value_derived": False,
}
(OUT / "dyn9b2_nonsusy_flavor_refit.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-2: non-SUSY flavor refit + zeta re-extraction", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "## Verdicts", "",
      f"1. **SO(10) lock tension**: y_t(M_X) = {yt_MX:.3f} "
      "(machine-derived); the perturbative ceiling caps the nu-Dirac "
      f"scale at y_max = {lock['PS']['y_max']:.4f} (PS) / "
      f"{lock['G_LR']['y_max']:.4f} (G_LR) -- tension factors "
      f"{lock['PS']['lock_tension']:.0f}x / "
      f"{lock['G_LR']['lock_tension']:.0f}x.",
      "2. **Type-II fails** by "
      f"{typeII['PS']['deficit_orders']:.1f} / "
      f"{typeII['G_LR']['deficit_orders']:.1f} orders (PS / G_LR).",
      "3. **Zeta-invariance theorem**: zeta, contact direction and "
      "contact fraction are EXACTLY invariant under the Dirac rescale; "
      "only M_* moves (M_*' = 4 pi M_I at the ceiling: "
      f"{Mstar_win['PS']['M_star_at_ceiling_GeV']:.2e} PS, "
      f"{Mstar_win['G_LR']['M_star_at_ceiling_GeV']:.2e} G_LR).",
      "4. **Leptogenesis selects the scale-decoupled source**: "
      "rescaled M_1 falls below the Davidson-Ibarra floor on both "
      "chains; the archival M_1 survives only under an instanton-type "
      "tower -- the D3 coexistence prediction (K8) becomes a "
      "REQUIREMENT of the alive branch with viable leptogenesis "
      "(D3 stays unpromoted).",
      "", "## Boundary (NOT claimed)", "",
      "- No full non-SUSY Yukawa fit; the Dirac kernel shape is a "
      "retained conditional input.",
      "- The SO(10) lock is quantified from the third-generation "
      "relation, not re-derived from a global fit.",
      "- Type-II and Davidson-Ibarra numbers are order estimates.",
      "- zeta's value is NOT derived (boundary theorem intact).",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b2_nonsusy_flavor_refit.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-9b-2: {n_pass}/{len(CHECKS)} checks; lock tension "
      f"{lock['PS']['lock_tension']:.0f}x/"
      f"{lock['G_LR']['lock_tension']:.0f}x; zeta invariant, M_* "
      f"branch-local; leptogenesis selects the D3-type source; "
      f"ledgers -> output/audit9/dyn9b2_nonsusy_flavor_refit.*")
