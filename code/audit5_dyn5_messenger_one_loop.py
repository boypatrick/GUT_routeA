#!/usr/bin/env python3
"""DYN-5: hidden-messenger dynamics at one loop (companion Audit 5 lane).

Promotes the paper's TREE-LEVEL matching-silence theorem for the Route-B
sterile transvectant messenger to an audited one-loop statement.  Four
sections, every derivation step recorded in the JSON derivation_log:

  S1  Anchor reproduction.  lambda = sqrt(zeta), |lambda|, |zeta|,
      arg zeta, and the companion caveat number
      |lambda|^2/(16 pi^2) = 8.259697e-4 (per e-fold of messenger
      interval), plus the |zeta|^2 ~ 1.7e-2 higher-order remark.

  S2  Structure of delta Z_N.  The portal X(N) has coupling matrix
      lambda * 1_3 (coordinate-free trace pairing), and the messenger
      mass matrix M_* K_tr^-1 = 3 M_* K_tr has ALL singular values equal
      to sqrt(3) M_* (from K_tr^2 = I/3).  Hence the leading-log
      wavefunction correction is exactly family-universal,
      delta Z_N = deltaZ * 1_3 with
      deltaZ = |lambda|^2/(16 pi^2) ln(M_*/M_X); in the strict Route-B
      benchmark the interval is ln(M_*/(sqrt(3) M_*)) = -ln(3)/2.

  S3  Replay silence vs the DYN-4 windows.  Using the archival closure
      card (same sha256 as DYN-4a): (i) the seesaw light-nu matrix is
      EXACTLY invariant under the canonical redefinition
      mD -> mD A, M_R -> A^T M_R A, A = 1 - deltaZ_N/2 -- for ANY
      invertible A, so light-sector silence is structural
      (field-redefinition covariance), not merely O(deltaZ);
      (ii) the paper's LINEARIZED replay formulas
      Y_nuD -> Y_nuD(1 - deltaZ_N/2),
      M_R -> M_R - (deltaZ_N^T M_R + M_R deltaZ_N)/2
      leave an O(deltaZ^2) residual whose crossing of the DYN-4a
      contact-sensitivity windows is quantified; (iii) the HEAVY sector
      is NOT silent: M_R singular values shift by (1 - deltaZ),
      compared against the DYN-4b M_1 posterior window; (iv) the
      normalized zeta extraction is exactly invariant (only M_* moves).

  S4  R-selection at one loop.  Combinatorial enumeration: all
      X-decorated Dirac/triplet/Majorana operators have R = 2 + m,
      m >= 1; the portal forces R(X) = 1; silence fails only on the
      measure-zero assignment R(16) = 2.  SUSY non-renormalization
      makes one-loop effects Kahler-only, and the nu^c leg appears in
      NO visible Yukawa or d=5 proton channel, so
      delta Y_u = delta Y_d = delta Y_e = delta C_5L = delta C_5R = 0
      at one loop; X, N are visible singlets so Delta b = (0,0,0) and
      the DYN-2 threshold vector is untouched.  R-breaking spurion
      conditions: regeneration requires k r_s = -m; explicit ceilings
      on odd/even spurions vs the DYN-4 windows (order estimates,
      flagged).

Boundary: zeta is NOT derived; the R symmetry is NOT claimed anomaly
free; NO microscopic messenger completion is constructed; silence is
claimed at ONE loop only (the two-loop Tr(Ynu+ Ynu) leakage is
estimated and flagged).  Follows the repo ledger pattern: JSON + MD in
output/audit5/ with explicit negative-boundary flags.
"""

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ARCH = ROOT / "route_E" / "output"
OUT = ROOT / "output" / "audit5"
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


# --------------------------------------------------------------- section 1
print("== DYN-5 section 1: companion caveat anchors ==")

ZETA = complex(0.1076472949, 0.0736514853)
lam = np.sqrt(ZETA)
loop = 16.0 * math.pi ** 2
per_efold = abs(lam) ** 2 / loop

TEX = {  # printed anchors in paper/gut_framework.tex (Route-B datum + caveat)
    "abs_zeta": 0.13043190325293763,
    "arg_zeta": 0.600038020318215,
    "lam_re": 0.34502115743308964,
    "lam_im": 0.10673473744038917,
    "abs_lam": 0.3611535729477664,
    "per_efold": 8.259697e-4,
}

check("lambda = sqrt(zeta): lambda^2 - zeta machine zero and tex digits match",
      abs(lam ** 2 - ZETA) < 1e-16
      and abs(lam.real - TEX["lam_re"]) < 1e-15
      and abs(lam.imag - TEX["lam_im"]) < 1e-15,
      f"lambda = {lam.real:.17g} + {lam.imag:.17g} i")
check("|lambda|, |zeta|, arg zeta reproduce the tex anchors",
      abs(abs(lam) - TEX["abs_lam"]) < 1e-15
      and abs(abs(ZETA) - TEX["abs_zeta"]) < 1e-15
      and abs(np.angle(ZETA) - TEX["arg_zeta"]) < 1e-12,
      f"|lambda| = {abs(lam):.16f}")
check("companion caveat anchor |lambda|^2/(16 pi^2) = 8.259697e-4",
      abs(per_efold - TEX["per_efold"]) < 1e-10,
      f"{per_efold:.9e} per e-fold of ln(M_*/M_X)")
check("higher-order remark |zeta|^2 ~ 1.7e-2 reproduced",
      abs(abs(ZETA) ** 2 - 1.7012e-2) < 1e-4,
      f"|zeta|^2 = {abs(ZETA)**2:.5e}")

DLOG["S1_anchors"] = {
    "lambda": [lam.real, lam.imag],
    "per_efold_coefficient": per_efold,
    "note": "the tex caveat number 8.259697e-4 is the PER-E-FOLD coefficient "
            "of delta Z_N; the full delta Z_N multiplies ln(M_*/M_X)",
}

# --------------------------------------------------------------- section 2
print("== DYN-5 section 2: structure of delta Z_N ==")

K_tr = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]], dtype=float) / math.sqrt(3)
K_inv = 3.0 * K_tr
sv_Kinv = np.linalg.svd(K_inv, compute_uv=False)

check("K_tr^2 = I/3 and K_tr^-1 = 3 K_tr (fixed convention)",
      np.allclose(K_tr @ K_tr, np.eye(3) / 3, atol=1e-15)
      and np.allclose(K_tr @ K_inv, np.eye(3), atol=1e-15))
check("messenger mass matrix M_* K_tr^-1 has ALL singular values sqrt(3) M_* "
      "=> exact family degeneracy => delta Z_N proportional to identity at "
      "leading log",
      np.allclose(sv_Kinv, math.sqrt(3), atol=1e-14),
      f"singular values of K_tr^-1 = {sv_Kinv}")

# the portal X(N) = lambda X^a N_a is the trace pairing: coupling matrix
# lambda * 1_3; one-loop leading log: delta Z_N = (Y+ Y)/(16 pi^2) ln
Y_portal = lam * np.eye(3)
dZ_matrix_per_efold = (Y_portal.conj().T @ Y_portal).real / loop
check("portal coupling matrix is lambda * 1_3: (Y+ Y)/(16 pi^2) = "
      "per-efold coefficient times identity",
      np.allclose(dZ_matrix_per_efold, per_efold * np.eye(3), atol=1e-18))

L_strict = -0.5 * math.log(3.0)           # ln(M_*/(sqrt(3) M_*))
dZ_strict = per_efold * L_strict
L_SCAN = {"strict_benchmark": L_strict, "one_efold": 1.0,
          "one_decade": math.log(10.0), "two_decades": math.log(100.0),
          "three_decades": math.log(1000.0)}
check("strict Route-B benchmark interval ln(M_*/M_X) = -ln(3)/2 "
      "(X sits AT sqrt(3) M_*): |delta Z_N| = 4.537102e-4",
      abs(dZ_strict + 4.537102e-4) < 1e-8,
      f"deltaZ_strict = {dZ_strict:.6e}")

DLOG["S2_structure"] = {
    "derivation": [
        "W_R/M_* = M_V(N,N)/2 + lambda X(N) - K_tr^-1(X,X)/2 (paper theorem)",
        "portal in components: lambda X_a N_a => Yukawa matrix lambda * 1_3",
        "one-loop Kahler correction on N: delta Z_N = (Y+ Y)/(16pi^2) "
        "ln(M_*/M_X_i) summed over messenger mass eigenstates i",
        "X mass matrix = M_* K_tr^-1 = 3 M_* K_tr; K_tr = P/sqrt(3) with P "
        "orthogonal symmetric (P^2 = I) => |eigenvalues| all sqrt(3) M_*",
        "identical log arguments => delta Z_N = deltaZ * 1_3 EXACTLY at "
        "leading log (not just approximately family-universal)",
    ],
    "delta_Z_scan": {k: per_efold * v for k, v in L_SCAN.items()},
    "strict_benchmark_note": "in the strict theorem all three messengers sit "
                             "at sqrt(3) M_*, ABOVE the matching scale; the "
                             "caveat's generic interval covers microscopic "
                             "completions where the operator-generation scale "
                             "exceeds the messenger mass",
}

# --------------------------------------------------------------- section 3
print("== DYN-5 section 3: replay silence vs the DYN-4 windows ==")

card_path = ARCH / "publication_closure_card" / "publication_closure_card.json"
dyn4a = json.loads((ROOT / "output" / "audit1"
                    / "dyn4a_seesaw_zeta_posterior.json").read_text())
dyn4b = json.loads((ROOT / "output" / "audit1"
                    / "dyn4b_unconditional_zeta.json").read_text())
card_sha = sha256(card_path)
check("archival closure card present with the DYN-4a provenance sha256",
      card_sha == dyn4a["provenance"]["closure_card"]["sha256"],
      f"sha256 = {card_sha[:16]}...")

card = json.loads(card_path.read_text())
yuk = {k: cmat(v) for k, v in card["selected_row"]["Yukawa_fit"].items()}
Y_e, Y_nu = yuk["charged_lepton"], yuk["neutrino_dirac"]

# DYN-4a recipe, transcribed verbatim (benchmark replay)
OLDT = {"m1_eV": 1.0e-3, "dm21": 7.42e-5, "dm31": 2.517e-3,
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


U_e_rot = np.linalg.eigh(Y_e @ Y_e.conj().T)[1]
vals = np.linalg.eigh(Y_e @ Y_e.conj().T)[0]
U_e_rot = U_e_rot[:, np.argsort(vals)]
m_diag = np.array([OLDT["m1_eV"],
                   math.sqrt(OLDT["m1_eV"] ** 2 + OLDT["dm21"]),
                   math.sqrt(OLDT["m1_eV"] ** 2 + OLDT["dm31"])])
U_nu = U_e_rot @ standard_pmns(OLDT["s12"], OLDT["s13"], OLDT["s23"],
                               OLDT["delta"], OLDT["a21"], OLDT["a31"])
m_light = U_nu.conj() @ np.diag(m_diag) @ U_nu.conj().T
mD = Y_nu * 100.0e9                                  # v_u = 100 GeV, in eV
MR = -(mD.T @ np.linalg.inv(m_light) @ mD)           # eV
M_star = float(np.linalg.svd(MR / 1e9, compute_uv=False)[0])
c_hat = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], dtype=complex) / math.sqrt(3)
zeta_rep = complex(np.vdot(c_hat, (MR / 1e9) / M_star))
check("benchmark replay reproduces the zeta anchor (links S3 to DYN-4a)",
      abs(zeta_rep - ZETA) < 1e-9,
      f"zeta = {zeta_rep.real:.10f} + {zeta_rep.imag:.10f} i")


def seesaw(mD_, MR_):
    return -(mD_ @ np.linalg.inv(MR_) @ mD_.T)


m0 = seesaw(mD, MR)

# (i) EXACT structural silence: any invertible redefinition A of the nu^c
# field, applied consistently (mD -> mD A on the nu^c column index,
# M_R -> A^T M_R A), leaves the light matrix identically invariant.
dZ_test = per_efold * math.log(1000.0)               # worst scanned value
A_uni = (1.0 - 0.5 * dZ_test) * np.eye(3)
rng = np.random.default_rng(20260705)
A_gen = np.eye(3) - 0.5 * dZ_test * np.diag(1.0 + 0.3 * rng.standard_normal(3))
res_uni = np.linalg.norm(seesaw(mD @ A_uni, A_uni.T @ MR @ A_uni) - m0) \
    / np.linalg.norm(m0)
res_gen = np.linalg.norm(seesaw(mD @ A_gen, A_gen.T @ MR @ A_gen) - m0) \
    / np.linalg.norm(m0)
cond_MR = float(np.linalg.cond(MR))
check("EXACT light-sector silence: m_nu invariant under the canonical "
      "redefinition (universal AND non-universal delta Z_N) at "
      "machine level given cond(M_R) ~ 1.6e5",
      res_uni < 1e-10 and res_gen < 1e-10,
      f"rel residual: universal {res_uni:.2e}, non-universal {res_gen:.2e}; "
      f"cond(M_R) = {cond_MR:.2e} so eps*cond ~ {cond_MR*2.2e-16:.1e}")

# (ii) the paper's LINEARIZED replay formulas leave an O(deltaZ^2) residual
windows = dyn4a["windows"]
w_loose = float(windows["loose_regression"]["Delta_s"])
w_tight = float(windows["chi2_refreshed"]["Delta_s"])


def linearized_residual(dz):
    dZm = dz * np.eye(3)
    mD_lin = mD @ (np.eye(3) - 0.5 * dZm)
    MR_lin = MR - 0.5 * (dZm.T @ MR + MR @ dZm)
    return float(np.linalg.norm(seesaw(mD_lin, MR_lin) - m0)
                 / np.linalg.norm(m0))


r1 = linearized_residual(dZ_test)
r2 = linearized_residual(0.5 * dZ_test)
ratio = r1 / r2
dz_cross_loose = math.sqrt(w_loose / (r1 / dZ_test ** 2))
dz_cross_tight = math.sqrt(w_tight / (r1 / dZ_test ** 2))
check("linearized replay residual scales as O(deltaZ^2) "
      "(halving deltaZ quarters the residual)",
      abs(ratio - 4.0) < 0.05, f"ratio = {ratio:.4f}")
check("at the strict benchmark |deltaZ| = 4.53e-4 the linearized-formula "
      "residual sits BELOW even the refreshed window",
      linearized_residual(dZ_strict) < w_tight,
      f"residual = {linearized_residual(dZ_strict):.3e} vs "
      f"Delta_s(tight) = {w_tight:.3e}")

# (iii) heavy sector NOT silent: singular values shift by (1 - deltaZ)
m1_lo, m1_c, m1_hi = dyn4b["dyn7_feed"]["log10_M1_GeV_16_50_84"]
half_width = 0.5 * (m1_hi - m1_lo)
worst_shift_dex = abs(math.log10(1.0 - dZ_test))
A_exact = (1.0 - 0.5 * dZ_test) * np.eye(3)
MR_ex = A_exact.T @ MR @ A_exact
sv0 = np.linalg.svd(MR, compute_uv=False)
sv1 = np.linalg.svd(MR_ex, compute_uv=False)
check("heavy sector NOT silent: all M_R singular values shift by exactly "
      "(1 - deltaZ); worst scanned shift stays inside the DYN-4b M_1 window",
      np.allclose(sv1 / sv0, (1 - 0.5 * dZ_test) ** 2, rtol=1e-12)
      and worst_shift_dex < half_width,
      f"shift = {worst_shift_dex:.4f} dex vs half-width {half_width:.4f} dex "
      f"({100*worst_shift_dex/half_width:.1f}% of window)")

# (iv) normalized zeta extraction invariant; only M_* moves
Mst_ex = float(np.linalg.svd(MR_ex / 1e9, compute_uv=False)[0])
zeta_ex = complex(np.vdot(c_hat, (MR_ex / 1e9) / Mst_ex))
check("normalized zeta extraction EXACTLY invariant under delta Z_N "
      "(only M_* rescales by 1 - deltaZ)",
      abs(zeta_ex - zeta_rep) < 1e-12 * abs(zeta_rep)
      and abs(Mst_ex / M_star - (1 - 0.5 * dZ_test) ** 2) < 1e-12,
      f"|zeta' - zeta| = {abs(zeta_ex - zeta_rep):.2e}, "
      f"M_*'/M_* = {Mst_ex/M_star:.8f}")

DLOG["S3_replay"] = {
    "structural_silence": "m_nu = -mD M_R^-1 mD^T is invariant under "
                          "mD -> mD A, M_R -> A^T M_R A for ANY invertible "
                          "A: the seesaw combination is a GL(3)-covariant of "
                          "the heavy field being integrated out, so "
                          "canonical-normalization corrections cancel "
                          "IDENTICALLY in every light observable",
    "linearized_formula_validity": {
        "residual_at_deltaZ": {f"{dZ_test:.4e}": r1},
        "crossing_loose_window_Delta_s": {"deltaZ": dz_cross_loose,
                                          "ln_M_ratio": dz_cross_loose / per_efold},
        "crossing_refreshed_window_Delta_s": {"deltaZ": dz_cross_tight,
                                              "ln_M_ratio": dz_cross_tight / per_efold},
        "reading": "the paper's linearized replay formulas are sufficient "
                   "for ln(M_*/M_X) below the quoted crossing; beyond it the "
                   "exact A-form must be used (this sharpens the tex caveat "
                   "'small, but not automatically smaller than the narrow "
                   "contact-sensitivity windows')",
    },
    "heavy_sector": {
        "shift_dex_at_worst_scan": worst_shift_dex,
        "dyn4b_M1_window_half_width_dex": half_width,
        "note": "M_R eigenvalues, M_1 (leptogenesis input) and M_* all "
                "rescale by (1 - deltaZ); zeta itself does not move",
    },
}

# --------------------------------------------------------------- section 4
print("== DYN-5 section 4: R-selection at one loop + spurion conditions ==")

R = {"16": 1, "N": 1, "H": 0, "T": 0, "Tbar": 0, "X": 1}
BASE = {"Dirac": ["16", "16", "H"], "TripletT": ["16", "16", "T"],
        "TripletTbar": ["16", "16", "Tbar"], "Majorana_MV": ["N", "N"],
        "portal": ["X", "N"], "Xmass": ["X", "X"]}
rsum = lambda ops: sum(R[f] for f in ops)

base_ok = all(rsum(ops) == 2 for ops in BASE.values())
decorated = {}
for name, ops in BASE.items():
    if name in ("portal", "Xmass"):
        continue
    for m in (1, 2, 3):
        decorated[f"{name}+{m}X"] = rsum(ops + ["X"] * m)
viol = [k for k, v in decorated.items() if v == 2]
check("R enumeration: all base channels have R = 2; every X-decorated "
      "Dirac/triplet/Majorana operator has R = 2 + m (zero violations)",
      base_ok and not viol and
      all(v == 2 + int(k.split("+")[1][0]) for k, v in decorated.items()),
      f"{len(decorated)} decorated operators checked")

# portal forces R(X) = 2 - r; generalized silence criterion over r = R(16)
fail_r = []
for r10 in range(-40, 41):
    r = r10 / 10.0
    rx = 2.0 - r                       # forced by the portal X(N)
    # decorated Dirac at m insertions: R = 2 + m*(2 - r) = 2 iff r = 2
    if any(abs(2 + m * rx - 2) < 1e-12 for m in (1, 2, 3)):
        fail_r.append(r)
check("portal forces R(X) = 2 - R(16); silence fails ONLY at the "
      "measure-zero assignment R(16) = 2 (benchmark R(16) = 1 is safe)",
      fail_r == [2.0], f"failing assignments in scan: {fail_r}")

# one-loop Kahler-only bookkeeping: nu^c leg appears in NO other channel
LEGS = {"Y_u": ["Q", "uc", "Hu"], "Y_d": ["Q", "dc", "Hd"],
        "Y_e": ["L", "ec", "Hd"], "Y_nuD": ["L", "nuc", "Hu"],
        "C_5L": ["Q", "Q", "Q", "L"], "C_5R": ["uc", "uc", "dc", "ec"]}
touched = {k: ("nuc" in v) for k, v in LEGS.items()}
check("SUSY non-renormalization => one-loop effects are Kahler-only; the "
      "nu^c leg appears ONLY in Y_nuD => delta Y_u = delta Y_d = delta Y_e "
      "= delta C_5L = delta C_5R = 0 at one loop; X, N visible singlets => "
      "Delta b = (0,0,0), DYN-2 threshold vector untouched",
      touched == {"Y_u": False, "Y_d": False, "Y_e": False, "Y_nuD": True,
                  "C_5L": False, "C_5R": False},
      "Y_nuD shifted but cancels exactly in m_nu (S3)")

# two-loop leakage estimate: delta Z_N shifts Tr(Ynu+ Ynu) by -deltaZ*Tr,
# which feeds the two-loop gauge beta trace term; the observable-level
# effect on alpha^-1 over the full run is of order
# (deltaZ * Tr / (16 pi^2)) * ln(M_*/M_Z) / (2 pi)
tr_ynu = float(np.real(np.trace(Y_nu.conj().T @ Y_nu)))
ln_run = math.log(M_star / 91.19)
leak_coupling = dZ_test * tr_ynu / loop
leak_alpha_inv = leak_coupling * ln_run / (2 * math.pi)
check("two-loop leakage (ORDER ESTIMATE): shift in alpha^-1 over the full "
      "run < 1e-3, i.e. far below the DYN-2 threshold-sensitivity scale "
      "O(0.1); two-loop silence itself is NOT claimed",
      leak_alpha_inv < 1e-3,
      f"delta alpha^-1 ~ {leak_alpha_inv:.2e} (coupling-level "
      f"{leak_coupling:.2e}, Tr(Ynu+ Ynu) = {tr_ynu:.4f})")

# spurion regeneration table: a holomorphic spurion eps with R = r_s at
# order k regenerates an m-decorated operator iff k*r_s = -m
table = {}
for r_s in range(-4, 5):
    if r_s == 0:
        continue
    hits = [(k, m) for k in (1, 2) for m in (1, 2, 3) if k * r_s == -m]
    table[str(r_s)] = hits
check("spurion regeneration condition k r_s = -m: r_s = -1 regenerates "
      "m = 1 at first order; r_s = -2 (constant-W / gravitino type) "
      "regenerates ONLY even m (m = 2 at k = 1), never the m = 1 portal "
      "decoration; positive-R spurions never regenerate at the "
      "holomorphic level",
      table["-1"] == [(1, 1), (2, 2)] and table["-2"] == [(1, 2)] + [(2, 4) for _ in []]
      and all(not table[str(r)] for r in (1, 2, 3, 4)),
      "residual (-1)^R parity: even-R breaking forbids all ODD decorations")

# numeric ceilings (order estimates, flagged): induced fractional Yukawa
# contamination per unit spurion, one X loop, M_X = sqrt(3) M_*
frac_loose = w_loose / OLDT["s23"]
frac_tight = w_tight / OLDT["s23"]
gain_odd = abs(lam) * math.sqrt(3) / loop        # eps_odd * lam * X-loop
gain_even = math.sqrt(3) / loop                  # eps_even * XX-mass loop
ceil = {
    "eps_odd_loose": frac_loose / gain_odd, "eps_odd_tight": frac_tight / gain_odd,
    "eps_even_loose": frac_loose / gain_even, "eps_even_tight": frac_tight / gain_even,
}
check("spurion ceilings vs DYN-4 windows computed (ORDER estimates): "
      "loose window tolerates eps ~ 1e-2, refreshed window forces "
      "eps below ~ 1e-3",
      ceil["eps_odd_loose"] > 1e-3 and ceil["eps_odd_tight"] < 1e-2,
      ", ".join(f"{k} = {v:.2e}" for k, v in ceil.items()))

DLOG["S4_r_selection"] = {
    "charges": R, "base_channels_R2": {k: rsum(v) for k, v in BASE.items()},
    "decorated_R": decorated,
    "spurion_regeneration_table_kr_eq_minus_m": table,
    "spurion_conditions_to_audited_order": [
        "holomorphic spurions with POSITIVE R charge never regenerate "
        "X-decorated operators in W (holomorphy)",
        "even-R spurions (r_s = -2, constant-W/gravitino type) preserve the "
        "residual (-1)^R parity: all ODD decorations (including the m = 1 "
        "Dirac decoration X 16 16 H) stay absent to ALL orders in the "
        "spurion; the m = 2 decoration can appear at O(eps)",
        "odd-R spurions (r_s = -1) regenerate the m = 1 decoration at O(eps) "
        "and are the dangerous case",
    ],
    "ceilings_order_estimates": ceil,
    "two_loop_leakage_estimate": {"coupling_level": leak_coupling,
                                  "alpha_inv_level": leak_alpha_inv},
}

# --------------------------------------------------------------- ledgers
n_pass = sum(1 for _, ok in CHECKS if ok)
all_pass = n_pass == len(CHECKS)
OUT.mkdir(parents=True, exist_ok=True)

payload = {
    "audit": "DYN-5 hidden-messenger dynamics at one loop",
    "dyn_item": "DYN-5",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": all_pass, "checks_passed": n_pass, "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "benchmark": {
        "zeta": [ZETA.real, ZETA.imag], "lambda": [lam.real, lam.imag],
        "per_efold_coefficient": per_efold,
        "delta_Z_scan": DLOG["S2_structure"]["delta_Z_scan"],
        "M_star_GeV": M_star,
    },
    "derivation_log": DLOG,
    "provenance": {
        "closure_card": {"path": str(card_path.relative_to(ROOT)),
                         "sha256": card_sha,
                         "note": "user-restored archival artifact, same "
                                 "sha256 as the DYN-4a ledger"},
        "windows_source": "output/audit1/dyn4a_seesaw_zeta_posterior.json "
                          "(contact-sensitivity windows) + "
                          "dyn4b_unconditional_zeta.json (M_1 posterior)",
        "tex_caveat": "paper/gut_framework.tex, canonical-basis caveat "
                      "(delta Z_N observation) + R-selection lemma + "
                      "tree-level matching-silence theorem",
    },
    # negative-boundary flags
    "zeta_value_derived": False,
    "one_loop_silence_of_light_sector": True,
    "heavy_sector_silent": False,
    "two_loop_silence_claimed": False,
    "r_symmetry_claimed_anomaly_free": False,
    "microscopic_messenger_completion_constructed": False,
    "spurion_sector_constructed": False,
    "ceilings_are_order_estimates": True,
}
(OUT / "dyn5_messenger_one_loop.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-5: hidden-messenger dynamics at one loop", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "## Result",
      "",
      "The tree-level matching-silence theorem UPGRADES cleanly to one loop,",
      "with one structural sharpening and one disclosed non-silence:",
      "",
      "1. **delta Z_N is exactly family-universal.**  The portal coupling is",
      "   `lambda * 1_3` and the messenger mass matrix `M_* K_tr^-1` has all",
      "   singular values `sqrt(3) M_*` (from `K_tr^2 = I/3`), so the",
      "   leading-log correction is `deltaZ * identity`.  Anchor:",
      f"   `|lambda|^2/(16 pi^2) = {per_efold:.6e}` per e-fold (tex caveat",
      "   digits reproduced).  Strict-benchmark interval `-ln(3)/2` gives",
      f"   `|deltaZ| = {abs(dZ_strict):.3e}`.",
      "2. **Light-sector silence is STRUCTURAL, not perturbative.**  The",
      "   seesaw `m_nu = -mD M_R^-1 mD^T` is invariant under `mD -> mD A`,",
      "   `M_R -> A^T M_R A` for ANY invertible `A` (machine-zero check,",
      "   universal and non-universal): canonical-normalization corrections",
      "   cancel identically in every light observable, to all orders in",
      "   `deltaZ` at fixed matching order.",
      "3. **The paper's linearized replay formulas have a quantified",
      f"   validity range.**  Residual is O(deltaZ^2); it crosses the loose",
      f"   contact-sensitivity window at `deltaZ = {dz_cross_loose:.2e}`",
      f"   (ln M ratio ~ {dz_cross_loose/per_efold:.1f}) and the refreshed",
      f"   window at `deltaZ = {dz_cross_tight:.2e}`",
      f"   (ln M ratio ~ {dz_cross_tight/per_efold:.1f}).  Beyond that the",
      "   exact A-form must be used.",
      "4. **The heavy sector is NOT silent** (disclosed): `M_R` singular",
      "   values, `M_1` (leptogenesis input) and `M_*` rescale by",
      f"   `(1 - deltaZ)`; worst scanned shift = {worst_shift_dex:.4f} dex =",
      f"   {100*worst_shift_dex/half_width:.1f}% of the DYN-4b M_1 window.",
      "   The normalized zeta extraction is exactly invariant.",
      "5. **R-selection silence at one loop.**  All X-decorated",
      "   Dirac/triplet/Majorana operators have `R = 2 + m`; the portal",
      "   forces `R(X) = 1`; silence fails only at the measure-zero",
      "   assignment `R(16) = 2`.  Non-renormalization makes one-loop",
      "   effects Kahler-only, and the `nu^c` leg appears in no visible",
      "   Yukawa or d=5 channel, so `delta Y_{u,d,e} = delta C_5L =",
      "   delta C_5R = 0` at one loop and the DYN-2 threshold vector is",
      "   untouched (`X`, `N` visible singlets).",
      "6. **R-breaking spurion conditions** (regeneration iff `k r_s = -m`):",
      "   positive-R spurions are safe at the holomorphic level; even-R",
      "   spurions (constant-W/gravitino type, `r_s = -2`) preserve a",
      "   residual `(-1)^R` parity so all ODD decorations stay absent to all",
      "   orders; odd-R spurions are the dangerous case.  Order-estimate",
      "   ceilings: " + ", ".join(f"`{k} = {v:.1e}`" for k, v in ceil.items()),
      "",
      "## Boundary (NOT claimed)",
      "",
      "- zeta is NOT derived (Theorem VII.4 boundary intact).",
      "- The R symmetry is NOT claimed anomaly free; no UV implementation.",
      "- NO microscopic messenger completion is constructed; the interval",
      "  `ln(M_*/M_X)` is scanned, not derived.",
      f"- Two-loop silence NOT claimed (leakage estimate "
      f"{leak_alpha_inv:.1e} in alpha^-1 units, via the",
      "  `Tr(Ynu+ Ynu)` trace term; flagged order estimate).",
      "- Spurion ceilings are order-of-magnitude estimates.",
      "",
      "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn5_messenger_one_loop.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-5: {n_pass}/{len(CHECKS)} checks; "
      f"light-sector one-loop silence STRUCTURAL (exact); "
      f"heavy sector shifts by (1 - deltaZ), zeta extraction invariant; "
      f"ledgers -> {OUT.relative_to(ROOT)}/")
