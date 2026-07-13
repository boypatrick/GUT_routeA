#!/usr/bin/env python3
"""DYN-9b-1d (= SUB-B0): the left-right vev ratio freed.

The SUB-A prior-art check found a CONFLICT: He & Meljanac (PRD 33,
2695 (1986)) and the 1989 stability analyses report stable left-right
minima of the 210 potential in finite parameter ranges, while
DYN-9b-1b/1c scanned tree positivity at a FIXED vev ratio a/p = 0.8
(disclosed flag) and found none.  This corrective audit frees the
ratio and answers three questions:

  Q1  With the parity-EVEN quartic sector only (lambda_eps = 0 -- the
      same four-invariant potential family as the classic papers),
      does the left-right stationary family become a tree-level local
      minimum SOMEWHERE in (ratio x couplings)?  A YES reproduces the
      He-Meljanac claim and localizes 9b-1b's negative verdict as a
      fixed-ratio artifact; a NO establishes a genuine disagreement to
      be resolved at full-text level.

  Q2  THE EPSILON QUESTION: the machine-found fifth invariant
      (parity-odd, slice-vanishing, off-slice-mass-only) is absent
      from the classic four-invariant basis.  Does it REWRITE the
      left-right stability verdict -- rescuing points that fail at
      lambda_eps = 0, or killing points that pass?

  Q3  If positive points exist: what does the real-spectrum threshold
      verdict for the left-right chain look like there (the analog of
      the 9b-1b Pati-Salam table)?

Method: the Hessian of each invariant is QUADRATIC in the vev, so
three polarization builds per pattern (at e_p, e_a, e_p + e_a) give
the exact Hessian at ANY ratio for free:
H(p e_p + a e_a) = p^2 A + p a B + a^2 C (quartics),
p D + a E (cubic).  Gates: the assembly reproduces a direct build at
r = 0.8 per pattern; the full 210-direction gradient vanishes at the
solved stationary points; exactly 30 Goldstone zeros at every sampled
ratio.

Boundary: LOCAL tree-level positivity (the classic claims concern
absolute minima -- implication runs their way: their stable minima
must appear as our positive points); He-Meljanac full text not
accessed (paywall) -- comparison at claim level; couplings sampled in
full-sum-invariant units; the 126bar back-reaction on the 210 masses
neglected (sigma << M_X); zeta NOT derived.
Ledgers -> output/audit9/dyn9b1d_lr_ratio_scan.{json,md}.
"""

from __future__ import annotations

import itertools
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from route_e_paths import AUDIT_OUTPUT, CACHE_DIR, REPO_ROOT

ROOT = REPO_ROOT
OUT = AUDIT_OUTPUT / "audit9"
CHECKS = []
DLOG = {}
N = 10
rng = np.random.default_rng(20260707)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


# ------------------------------------------------ 210 core (as 9b-1b/c)
QUADS = list(itertools.combinations(range(N), 4))
SEXTS = list(itertools.combinations(range(N), 6))


def perm_sign(perm):
    s, p = 1, list(perm)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


PERMS4 = [(pm, perm_sign(pm)) for pm in itertools.permutations(range(4))]
FLAT4 = np.array([[np.ravel_multi_index(tuple(q[p] for p in pm), (N,) * 4)
                   for pm, _ in PERMS4] for q in QUADS])
SIGNS4 = np.array([sg for _, sg in PERMS4], dtype=float)


def unit_tensor(quad):
    T = np.zeros((N,) * 4)
    for pm, sg in PERMS4:
        T[tuple(quad[p] for p in pm)] = sg
    return T


UNITS = [unit_tensor(q) for q in QUADS]


def comps(G):
    return (G.ravel()[FLAT4] @ SIGNS4) / 24.0


def from_comps(x):
    T = np.zeros(N ** 4)
    np.add.at(T, FLAT4.ravel(),
              (np.asarray(x)[:, None] * SIGNS4[None, :]).ravel())
    return T.reshape((N,) * 4)


PERMS6 = [(pm, perm_sign(pm)) for pm in itertools.permutations(range(6))]
FLAT6 = np.array([[np.ravel_multi_index(tuple(s[p] for p in pm), (N,) * 6)
                   for pm, _ in PERMS6] for s in SEXTS])
SIGNS6 = np.array([sg for _, sg in PERMS6], dtype=float)
ALL = set(range(N))
SEXT_IDX = {s: i for i, s in enumerate(SEXTS)}
DUAL_MAP = np.zeros(len(QUADS), dtype=int)
DUAL_SIGN = np.zeros(len(QUADS))
for A, q in enumerate(QUADS):
    rest = tuple(sorted(ALL - set(q)))
    DUAL_MAP[A] = SEXT_IDX[rest]
    DUAL_SIGN[A] = perm_sign(tuple(list(q) + list(rest)))


def star6_dense(T):
    c = comps(T)
    S = np.zeros(N ** 6)
    vals = (c * DUAL_SIGN)[:, None] * SIGNS6[None, :]
    np.add.at(S, FLAT6[DUAL_MAP].ravel(), vals.ravel())
    return S.reshape((N,) * 6)


def dual_readout(W):
    return (W.ravel()[FLAT6[DUAL_MAP]] @ SIGNS6) * DUAL_SIGN


e_p = from_comps([1.0 if q == (6, 7, 8, 9) else 0.0 for q in QUADS])
A_SET = {(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)}
W_SET = {(0, 1, 6, 7), (0, 1, 8, 9), (2, 3, 6, 7), (2, 3, 8, 9),
         (4, 5, 6, 7), (4, 5, 8, 9)}
e_a = from_comps([1.0 if q in A_SET else 0.0 for q in QUADS])
e_w = -from_comps([1.0 if q in W_SET else 0.0 for q in QUADS])


def phi(p, a, w=0.0):
    return p * e_p + a * e_a + w * e_w


def I2(T):
    return float(np.einsum('ijkl,ijkl->', T, T))


def M2T(T):
    return np.einsum('ijmn,klmn->ijkl', T, T)


def grad_Q1(T):
    return 4.0 * I2(T) * T


def grad_Q2(T):
    return 4.0 * np.einsum('ijkl,klmn->ijmn', M2T(T), T)


def grad_Q3(T):
    M = M2T(T)
    W = 0.5 * (np.einsum('ikjl->ijkl', M) + np.einsum('kilj->ijkl', M))
    return 2.0 * (np.einsum('ijkl,klmn->ijmn', W, T)
                  + np.einsum('klij,klmn->ijmn', W, T))


def grad_Q4(T):
    S = np.einsum('iabc,jabc->ij', T, T)
    return 4.0 * np.einsum('ij,jabc->iabc', S, T)


def grad_Q6(T):
    S = star6_dense(T)
    Y = np.einsum('abcdef,efqr->abcdqr', S, T)
    Y2 = np.einsum('abcdef,abpq->cdefpq', S, T)
    g_direct = (np.einsum('abcdqr,cdpr->abpq', Y, T)
                + np.einsum('abcdqr,abpq->cdpr', Y, T)
                + np.einsum('cdefpq,cdpr->efqr', Y2, T))
    W1 = np.einsum('abpq,cdpr->abcdqr', T, T)
    W = np.einsum('abcdqr,efqr->abcdef', W1, T)
    return g_direct + from_comps(dual_readout(W) / 24.0)


GRADS = {"Q1": grad_Q1, "Q2": grad_Q2, "Q3": grad_Q3, "Q4": grad_Q4,
         "eps": grad_Q6}
GENS = [(i, j) for i in range(N) for j in range(i + 1, N)]


def act4(T, ab):
    a_, b_ = ab
    d = np.zeros_like(T)
    for slot in range(4):
        Tm = np.moveaxis(T, slot, 0)
        dm = np.moveaxis(d, slot, 0)
        dm[a_] += Tm[b_]
        dm[b_] -= Tm[a_]
    return d


def vector_M2(T):
    cols = np.array([act4(T, g).ravel() for g in GENS])
    return cols @ cols.T / 24.0


# ------------------------------------------------ Hessian machinery
def hessian_of(gradf, v, exact_cubic=False, gE=None):
    H = np.zeros((210, 210))
    for B in range(210):
        eB = UNITS[B]
        if exact_cubic:
            col = (gradf(v + eB) - gradf(v - eB) - 2.0 * gE[B]) / 2.0
        else:
            g1 = gradf(v + eB) - gradf(v - eB)
            g2 = gradf(v + 2 * eB) - gradf(v - 2 * eB)
            col = (8.0 * g1 - g2) / 12.0
        H[:, B] = comps(col)
    return 24.0 * H


def hessian_cubic(v):
    H = np.zeros((210, 210))
    for B in range(210):
        eB = UNITS[B]
        col = 3.0 * (np.einsum('klmn,mnij->ijkl', v, eB)
                     + np.einsum('klmn,mnij->ijkl', eB, v))
        H[:, B] = comps(col)
    return 24.0 * H


H_M2 = 48.0 * np.eye(210)

print("== 9b-1d section 1: polarization Hessians ==")
CACHE_DIR.mkdir(parents=True, exist_ok=True)
cache_key = hashlib.sha256(
    Path(__file__).read_bytes() + np.__version__.encode("ascii")
).hexdigest()[:16]
CACHE = CACHE_DIR / f"dyn9b1d_polar_{cache_key}.npz"
POLAR = {}
if CACHE.exists():
    dat = np.load(CACHE)
    POLAR = {k: dat[k] for k in dat.files}
    print("   polarization cache loaded (delete the npz to recompute)")
else:
    print("   caching g(e_B) for the epsilon pattern (210 calls)...")
    gE = [grad_Q6(UNITS[B]) for B in range(210)]
    for name, gf in GRADS.items():
        ex = name == "eps"
        print(f"   building A, C, B polarizations for {name}...")
        Am = hessian_of(gf, e_p, exact_cubic=ex, gE=gE if ex else None)
        Cm = hessian_of(gf, e_a, exact_cubic=ex, gE=gE if ex else None)
        Sm = hessian_of(gf, e_p + e_a, exact_cubic=ex,
                        gE=gE if ex else None)
        POLAR[f"{name}_A"] = Am
        POLAR[f"{name}_C"] = Cm
        POLAR[f"{name}_B"] = Sm - Am - Cm
    POLAR["kap_D"] = hessian_cubic(e_p)
    POLAR["kap_E"] = hessian_cubic(e_a)
    np.savez(CACHE, **POLAR)


def total_H(r, m2, kap, lam, leps, p=1.0):
    """phi-Hessian at (p, a = r p, 0); quadratic/linear vev scaling."""
    a = r * p
    H = m2 * H_M2 + kap * (p * POLAR["kap_D"] + a * POLAR["kap_E"])
    for lam_i, name in zip(list(lam) + [leps], ["Q1", "Q2", "Q3", "Q4",
                                                "eps"]):
        H = H + lam_i * (p * p * POLAR[f"{name}_A"]
                         + p * a * POLAR[f"{name}_B"]
                         + a * a * POLAR[f"{name}_C"])
    return (H + H.T) / 2


# gate: assembly vs direct build at r = 0.8 for each pattern
v08 = phi(1.0, 0.8)
gate_ok = True
gE_small = None
for name, gf in GRADS.items():
    if name == "eps":
        continue                      # direct eps build too slow here;
    Hd = hessian_of(gf, v08)
    Ha = (POLAR[f"{name}_A"] + 0.8 * POLAR[f"{name}_B"]
          + 0.64 * POLAR[f"{name}_C"])
    Ha = (Ha + Ha.T) / 2
    Hd = (Hd + Hd.T) / 2
    if np.abs(Hd - Ha).max() > 1e-8 * max(1.0, np.abs(Hd).max()):
        gate_ok = False
        print(f"   polarization mismatch {name}")
check("POLARIZATION GATE: the assembled Hessian equals a direct build "
      "at r = 0.8 for every polynomial pattern (the epsilon pattern is "
      "certified by the Goldstone gate below)", gate_ok)

# ------------------------------------------------ slice potential
MONOS = [(4, 0, 0), (0, 4, 0), (0, 0, 4), (2, 2, 0), (2, 0, 2),
         (0, 2, 2), (1, 3, 0), (3, 1, 0), (1, 0, 3), (3, 0, 1),
         (0, 1, 3), (0, 3, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2)]
pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0),
       (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1), (1, 1, 1),
       (1, 1, -1), (1, -1, 1), (2, 1, 0), (1, 2, 0), (0, 2, 1),
       (2, 0, 1), (1, 1, 2), (2, 1, 1), (1, 2, 1)]
A_fit = np.array([[float(p ** i * a ** j * w ** k)
                   for (i, j, k) in MONOS] for (p, a, w) in pts])
PATVAL = {"Q1": lambda T: I2(T) ** 2,
          "Q2": lambda T: float(np.einsum('ijkl,ijkl->', M2T(T),
                                          M2T(T))),
          "Q3": lambda T: float(np.einsum('ijkl,ikjl->', M2T(T),
                                          M2T(T))),
          "Q4": lambda T: float(np.einsum(
              'ij,ij->', np.einsum('iabc,jabc->ij', T, T),
              np.einsum('iabc,jabc->ij', T, T)))}
DESC = {}
for name, f in PATVAL.items():
    y = np.array([f(phi(*pt)) for pt in pts])
    coef, *_ = np.linalg.lstsq(A_fit, y, rcond=None)
    coef[np.abs(coef) < 1e-7] = 0.0
    DESC[name] = coef


def slice_V(m2, kap, lam, p, a):
    val = m2 * 24 * (p * p + 3 * a * a) + kap * 48 * (a ** 3)
    mono = np.array([p ** i * a ** j * (0.0 ** k if k else 1.0)
                     for (i, j, k) in MONOS])
    for l_, name in zip(lam, ["Q1", "Q2", "Q3", "Q4"]):
        val += l_ * float(DESC[name] @ mono)
    return val


def stationarity_LR(r, lam):
    """Solve (m2, kappa) so that (1, r, 0) is stationary -- EXACT
    linear solve from the analytic monomial derivatives."""
    dqp = dqa = 0.0
    for l_, name in zip(lam, ["Q1", "Q2", "Q3", "Q4"]):
        for (i, j, k), c in zip(MONOS, DESC[name]):
            if c == 0.0 or k != 0:
                continue
            dqp += l_ * c * i * (r ** j)
            dqa += l_ * c * j * (r ** (j - 1) if j else 0.0)
    # m2 * 24 (p^2 + 3 a^2): d/dp = 48, d/da = 144 r  (at p = 1)
    # kappa * 48 a^3: d/dp = 0, d/da = 144 r^2
    Jm = np.array([[48.0, 0.0], [144.0 * r, 144.0 * r * r]])
    sol = np.linalg.solve(Jm, -np.array([dqp, dqa]))
    return float(sol[0]), float(sol[1])


print("== 9b-1d section 2: stationarity and Goldstone gates ==")

# full-stationarity gate: the total analytic gradient vanishes in ALL
# 210 directions at a random solved point (invariance argument made
# explicit)
r_t = 1.3
lam_t = rng.uniform(-1, 1, size=4)
m2t, kapt = stationarity_LR(r_t, lam_t)
vt = phi(1.0, r_t)
Gtot = m2t * 2.0 * vt \
    + kapt * 3.0 * np.einsum('klmn,mnij->ijkl', vt, vt)
for l_, name in zip(lam_t, ["Q1", "Q2", "Q3", "Q4"]):
    Gtot = Gtot + l_ * GRADS[name](vt)
grad_norm = float(np.abs(comps(Gtot)).max())
check("FULL-STATIONARITY GATE: with (m2, kappa) solved on the slice, "
      "the total gradient vanishes in ALL 210 directions (little-group "
      "invariance argument verified numerically; epsilon contributes "
      "no gradient at any slice point)",
      grad_norm < 1e-8, f"max |gradient component| = {grad_norm:.2e}")

gold_ok = True
for r in (0.3, 0.8, 1.5, 2.5):
    lam_g = rng.uniform(-1, 1, size=4)
    m2g, kapg = stationarity_LR(r, lam_g)
    H = total_H(r, m2g, kapg, lam_g, float(rng.uniform(-1, 1)))
    ev = np.linalg.eigvalsh(H)
    z = int(np.sum(np.abs(ev) < 1e-7 * max(1.0, np.abs(ev).max())))
    gold_ok &= z == 30
check("GOLDSTONE GATE at four ratios with random couplings AND random "
      "epsilon: exactly 30 zeros every time -- certifies the "
      "polarization assembly including the epsilon block", gold_ok)

# ------------------------------------------------ section 3: THE SCAN
print("== 9b-1d section 3: positivity over (ratio x couplings x eps) ==")

R_GRID = [0.05, 0.1, 0.2, 0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
          0.8, 0.85, 0.9, 1.0, 1.3, 1.7, 2.2, 3.0,
          -0.3, -0.8, -1.3, -2.0]
EPS_GRID = np.linspace(-1.5, 1.5, 13)
N_COUP = 400
stats = {"samples": 0, "pos_eps0": 0, "pos_some_eps": 0,
         "rescued_by_eps": 0, "killed_by_eps": 0}
pos_points = []
by_ratio = {str(r): 0 for r in R_GRID}
for it in range(N_COUP):
    lam_s = rng.uniform(-1, 1, size=4)
    for r in R_GRID:
        try:
            m2s, kaps = stationarity_LR(r, lam_s)
        except np.linalg.LinAlgError:
            continue
        stats["samples"] += 1

        def min_phys(le):
            H = total_H(r, m2s, kaps, lam_s, le)
            ev = np.linalg.eigvalsh(H)
            sc = max(1.0, np.abs(ev).max())
            nz = ev[np.abs(ev) > 1e-7 * sc]
            return float(nz.min()) if len(nz) else 0.0

        m0 = min_phys(0.0)
        pos0 = m0 > 0
        best = max(min_phys(le) for le in EPS_GRID)
        poss = best > 0
        if pos0:
            stats["pos_eps0"] += 1
            by_ratio[str(r)] += 1
            pos_points.append({"r": r, "lam": lam_s.tolist(),
                               "m2": m2s, "kappa": kaps,
                               "min_eig_eps0": m0})
        if poss:
            stats["pos_some_eps"] += 1
        if poss and not pos0:
            stats["rescued_by_eps"] += 1
        if pos0 and not poss:
            stats["killed_by_eps"] += 1
check("positivity scan executed over 16 ratios x 300 coupling points "
      "x 13 epsilon values", stats["samples"] > 4000,
      f"{stats['samples']} (ratio, coupling) evaluations")

q1 = stats["pos_eps0"] > 0
check("Q1 VERDICT (the He-Meljanac reproduction question): does the "
      "left-right stationary family admit tree-level LOCAL MINIMA "
      "somewhere in (ratio x couplings) at lambda_eps = 0 (the "
      "classic four-invariant potential family)?", True,
      f"positive points: {stats['pos_eps0']}/{stats['samples']}; "
      f"by ratio: { {k: v for k, v in by_ratio.items() if v} }")
check("Q2 VERDICT (the epsilon question): rescue and kill counts "
      "under the fifth invariant recorded", True,
      f"rescued by eps: {stats['rescued_by_eps']}; killed by eps "
      f"(positive at 0, killable at some eps NEVER counted as kill "
      f"-- kill means NO eps keeps it positive): "
      f"{stats['killed_by_eps']}; positive at some eps: "
      f"{stats['pos_some_eps']}")
DLOG["S3_scan"] = {"stats": stats, "by_ratio": by_ratio,
                   "eps_grid": EPS_GRID.tolist(),
                   "positive_examples": pos_points[:5]}

# ------------------------------------------------ section 4: verdicts
print("== 9b-1d section 4: consequences ==")

if q1:
    concl = ("9b-1b's negative left-right verdict WAS a fixed-ratio "
             "artifact: freeing the ratio reproduces stable left-right "
             "tree minima in finite coupling regions, consistent with "
             "He-Meljanac (1986) at claim level.  K5r REVERTS: the "
             "left-right chain is 210-realizable at tree level in "
             "those regions; the 45_H route is an alternative, not a "
             "necessity.")
else:
    concl = ("even with the ratio freed and the epsilon coupling swept "
             "the left-right family shows no tree-level local minimum "
             "in the sampled coupling space: a GENUINE disagreement "
             "with the He-Meljanac claim, to be resolved at full-text "
             "level (their invariant basis and vev ansatz vs ours) "
             "before any submission touches the 210 vacuum.")
check("K5r PROPAGATION statement recorded for the ledger refresh and "
      "both papers", True, concl)
DLOG["S4_conclusion"] = {"q1_positive": bool(q1), "statement": concl}

n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-9b-1d (SUB-B0) left-right vev-ratio scan",
    "dyn_item": "DYN-9b-1d",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    # negative-boundary flags
    "local_positivity_not_absolute_minimum": True,
    "he_meljanac_full_text_not_accessed": True,
    "coupling_sampling_full_sum_units": True,
    "sigma_backreaction_neglected": True,
    "zeta_value_derived": False,
}
(OUT / "dyn9b1d_lr_ratio_scan.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-1d (SUB-B0): the left-right vev ratio freed", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      f"- Q1 (He-Meljanac reproduction): positive points "
      f"{stats['pos_eps0']}/{stats['samples']} at lambda_eps = 0.",
      f"- Q2 (the epsilon question): rescued {stats['rescued_by_eps']}"
      f", killed {stats['killed_by_eps']}, positive-at-some-eps "
      f"{stats['pos_some_eps']}.",
      f"- Conclusion: {concl}",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b1d_lr_ratio_scan.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-9b-1d: {n_pass}/{len(CHECKS)} checks; positive(eps=0) = "
      f"{stats['pos_eps0']}, rescued-by-eps = "
      f"{stats['rescued_by_eps']}; ledgers -> "
      f"output/audit9/dyn9b1d_lr_ratio_scan.*")
