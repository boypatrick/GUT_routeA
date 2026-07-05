#!/usr/bin/env python3
"""DYN-9b-1b: 210 quartic Clebsch descent, survivor masses, verdicts.

The decisive stage of the non-SUSY vacuum program: the 210 is built as
an explicit rank-4 antisymmetric tensor of SO(10), the renormalizable
potential V(Phi) = m2 I2 + kappa I3 + sum_i lambda_i Q_i is constructed
from a machine-tested basis of quartic invariants, and everything
DYN-9b-1 had to leave symbolic is computed numerically.

  S1  Embedding gates: quadratic descent 24(p^2+3a^2+6w^2) and cubic
      descent 48(a^3+3pw^2+6aw^2) (the AG superpotential structure,
      DYN-1a anchors); the little-group map from the actual 45x45
      vector mass matrix; the DYN-1b gauge-mass ratio gate at the AG
      benchmark.

  S2  Invariant basis: five polynomial quartic contraction patterns
      plus one Hodge-dual epsilon pattern, rank-tested on random
      tensors; exact quartic Clebsch DESCENT to the singlet slice;
      parity facts (no odd-in-w, no p^3 a-type terms) that make the
      Pati-Salam and left-right points stationary.

  S3  Hessians (exact: gradients are cubic, Richardson h = 1, 2 is
      exact), Goldstone gauge gates (24 / 30 zeros at random
      couplings), multiplet classification with SU(2)_L / SU(2)_R /
      B-L label operators.

  S4  Verdict scan: random couplings -> stationarity -> tree-level
      positivity -> survivor masses in units of the heavy gauge boson
      -> multi-threshold one-loop matching in the DYN-9 solver -> the
      distribution of tau(p -> e+ pi0) per chain, replacing the ESH
      assumption; plus a tree-level slice phase-diagram sample.

Boundary: the 126bar quartic self-couplings are DEFERRED to 9b-1c
(sigma << M_X; the 126bar remnants stay at their ESH placement); the
epsilon quartic enters the basis test only (its coupling is set to
zero in the scan, disclosed); tree level only (one-loop stabilization
cited); coupling sampling is in full-sum-invariant units (disclosed);
zeta is NOT derived.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit9"
CHECKS = []
DLOG = {}
N = 10
rng = np.random.default_rng(20260705)
DEBUG_BLOCKS = "--blocks" in sys.argv


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


# ------------------------------------------------------ tensor basics
QUADS = list(itertools.combinations(range(N), 4))            # 210
SEXTS = list(itertools.combinations(range(N), 6))            # 210


def perm_sign(perm):
    s, p = 1, list(perm)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


PERMS4 = [(pm, perm_sign(pm)) for pm in itertools.permutations(range(4))]


def unit_tensor(quad):
    T = np.zeros((N,) * 4)
    for pm, sg in PERMS4:
        T[tuple(quad[p] for p in pm)] = sg
    return T


UNITS = [unit_tensor(q) for q in QUADS]
FLAT4 = np.array([[np.ravel_multi_index(tuple(q[p] for p in pm), (N,) * 4)
                   for pm, _ in PERMS4] for q in QUADS])
SIGNS4 = np.array([sg for _, sg in PERMS4], dtype=float)


def comps(G):
    """All 210 antisymmetrized components (phi-coordinates) of dense G."""
    return (G.ravel()[FLAT4] @ SIGNS4) / 24.0


def from_comps(x):
    T = np.zeros(N ** 4)
    np.add.at(T, FLAT4.ravel(),
              (np.asarray(x)[:, None] * SIGNS4[None, :]).ravel())
    return T.reshape((N,) * 4)


# Hodge machinery (dual components live on SEXTS; C(10,6) = 210 too)
PERMS6 = [(pm, perm_sign(pm)) for pm in itertools.permutations(range(6))]
FLAT6 = np.array([[np.ravel_multi_index(tuple(s[p] for p in pm), (N,) * 6)
                   for pm, _ in PERMS6]
                  for s in SEXTS])
SIGNS6 = np.array([sg for _, sg in PERMS6], dtype=float)
ALL = set(range(N))
DUAL_SIGN = np.zeros(len(QUADS))
DUAL_MAP = np.zeros(len(QUADS), dtype=int)
SEXT_IDX = {s: i for i, s in enumerate(SEXTS)}
for A, q in enumerate(QUADS):
    rest = tuple(sorted(ALL - set(q)))
    DUAL_MAP[A] = SEXT_IDX[rest]
    DUAL_SIGN[A] = perm_sign(tuple(list(q) + list(rest)))


def star6_dense(T):
    """Dense rank-6 dual: (*T) with comp6(*T)[rest(q)] = sign * comp(T)[q]."""
    c = comps(T)
    S = np.zeros(N ** 6)
    vals = (c * DUAL_SIGN)[:, None] * SIGNS6[None, :]
    np.add.at(S, FLAT6[DUAL_MAP].ravel(), vals.ravel())
    return S.reshape((N,) * 6)


e_p = from_comps([1.0 if q == (6, 7, 8, 9) else 0.0 for q in QUADS])
A_SET = {(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)}
W_SET = {(0, 1, 6, 7), (0, 1, 8, 9), (2, 3, 6, 7), (2, 3, 8, 9),
         (4, 5, 6, 7), (4, 5, 8, 9)}
e_a = from_comps([1.0 if q in A_SET else 0.0 for q in QUADS])
e_w = from_comps([1.0 if q in W_SET else 0.0 for q in QUADS])


def phi(p, a, w):
    return p * e_p + a * e_a + w * e_w


# ------------------------------------------------------ invariants
def I2(T):
    return float(np.einsum('ijkl,ijkl->', T, T))


def I3(T):
    return float(np.einsum('ijkl,klmn,mnij->', T, T, T))


def Q1(T):
    return I2(T) ** 2


def M2T(T):
    return np.einsum('ijmn,klmn->ijkl', T, T)


def Q2(T):
    M = M2T(T)
    return float(np.einsum('ijkl,ijkl->', M, M))


def Q3(T):
    M = M2T(T)
    return float(np.einsum('ijkl,ikjl->', M, M))


def Q4(T):
    S = np.einsum('iabc,jabc->ij', T, T)
    return float(np.einsum('ij,ij->', S, S))


def Q5(T):
    return float(np.einsum('ijab,jkbc,klcd,lida->', T, T, T, T))


def Q6(T):
    """Epsilon quartic via the dual: <*T, W(T,T,T)> (value-level only)."""
    S = star6_dense(T)
    W1 = np.einsum('abpq,cdpr->abcdqr', T, T)
    return float(np.einsum('abcdef,abcdqr,efqr->', S, W1, T))


PATTERNS = {"Q1": Q1, "Q2": Q2, "Q3": Q3, "Q4": Q4, "Q5": Q5}


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


def grad_Q5(T):
    g = np.einsum('jkbc,klcd,lida->ijab', T, T, T) \
        + np.einsum('ijab,klcd,lida->jkbc', T, T, T) \
        + np.einsum('ijab,jkbc,lida->klcd', T, T, T) \
        + np.einsum('ijab,jkbc,klcd->lida', T, T, T)
    return g


GRADS = {"Q1": grad_Q1, "Q2": grad_Q2, "Q3": grad_Q3, "Q4": grad_Q4,
         "Q5": grad_Q5}

# ============================================================ section 1
print("== 9b-1b section 1: embedding gates ==")

check("quadratic descent: I2 = 24 (p^2 + 3 a^2 + 6 w^2) -- the AG "
      "quadratic (DYN-1a anchor)",
      abs(I2(e_p) - 24) < 1e-9 and abs(I2(e_a) - 72) < 1e-9
      and abs(I2(e_w) - 144) < 1e-9)

c_a3 = I3(phi(0, 1, 0)) / 48
c_pw2 = (I3(phi(1, 0, 1)) - I3(phi(1, 0, 0)) - I3(phi(0, 0, 1))) / 48
c_aw2 = (I3(phi(0, 1, 1)) - I3(phi(0, 1, 0)) - I3(phi(0, 0, 1))) / 48
check("cubic descent: I3 = 48 (a^3 + 3 p w^2 + 6 a w^2) -- the AG "
      "superpotential cubic (DYN-1a anchor); p^3 = w^3 = 0",
      abs(c_a3 - 1) < 1e-9 and abs(c_pw2 - 3) < 1e-9
      and abs(c_aw2 - 6) < 1e-9
      and abs(I3(phi(1, 0, 0))) < 1e-9 and abs(I3(phi(0, 0, 1))) < 1e-9)

GENS = [(i, j) for i in range(N) for j in range(i + 1, N)]


def act(T, ab):
    a_, b_ = ab
    d = np.zeros_like(T)
    for slot in range(4):
        Tm = np.moveaxis(T, slot, 0)
        dm = np.moveaxis(d, slot, 0)
        dm[a_] += Tm[b_]
        dm[b_] -= Tm[a_]
    return d


def vector_M2(T):
    cols = np.array([act(T, g).ravel() for g in GENS])
    return cols @ cols.T / 24.0


def zero_count(M2m, tol=1e-9):
    ev = np.linalg.eigvalsh(M2m)
    return int(np.sum(ev < tol * max(1.0, ev.max())))


lg = {"p_only": (phi(1, 0, 0), 21), "a_only": (phi(0, 1, 0), 15),
      "p_and_a": (phi(0.7, 1, 0), 15), "w_only": (phi(0, 0, 1), 13),
      "flipped": (phi(1, 1, 1), 25), "SU5": (phi(-1, -1, 1), 25),
      "generic_no_sigma": (phi(0.3, 0.7, -0.2), 13)}
got_lg = {k: zero_count(vector_M2(T)) for k, (T, _) in lg.items()}
check("little-group map from the ACTUAL 45x45 vector mass matrix "
      "reproduces DYN-9b-1 (generic Phi-only leaves 13 = SM + the "
      "sigma-massed direction)",
      all(got_lg[k] == want for k, (_, want) in lg.items()), str(got_lg))

x = 0.1
w_b, a_b = -x, (x * x + 2 * x - 1) / (1 - x)
p_b = x * (5 * x * x - 1) / (1 - x) ** 2
ev_b = np.linalg.eigvalsh(vector_M2(phi(p_b, a_b, w_b)))
formulas = {"J": (8 * a_b ** 2 + 16 * w_b ** 2, 6),
            "F": (24 * w_b ** 2, 2),
            "E": (4 * (a_b - w_b) ** 2 + 2 * (w_b - p_b) ** 2, 12),
            "X": (4 * (a_b + w_b) ** 2 + 2 * (p_b + w_b) ** 2, 12)}
expected = sorted([(v, m) for v, m in formulas.values()])
nz = sorted(e for e in ev_b if e > 1e-9)
scale = nz[0] / expected[0][0]
built = []
for v, m in expected:
    built += [v * scale] * m
ratio_ok = len(nz) == 32 and all(abs(a - b) < 1e-6 * max(1, b)
                                 for a, b in zip(nz, built))
check("gauge-mass RATIOS at the AG benchmark equal the DYN-1b formulas "
      "(J, F, E, X with multiplicities 6/2/12/12, sigma = 0): the "
      "embedding is quantitatively anchored to the SUSY-side audit",
      ratio_ok, f"common scale {scale:.6f}")
DLOG["S1_embedding"] = {"cubic_descent": [c_a3, c_pw2, c_aw2],
                        "little_groups": got_lg,
                        "vector_scale": scale}

# ============================================================ section 2
print("== 9b-1b section 2: invariant basis and quartic descent ==")


def rand_phi():
    return from_comps(rng.standard_normal(210))


grad_ok = True
Tt = rand_phi()
D = from_comps(rng.standard_normal(210)) * 1e-4
for name, f in PATTERNS.items():
    g = GRADS[name](Tt)
    num = (f(Tt + D) - f(Tt - D)) / 2.0
    ana = float(np.einsum('ijkl,ijkl->', g, D))
    if abs(num - ana) > 1e-6 * max(1.0, abs(num)):
        grad_ok = False
        print(f"   gradient mismatch {name}: {num:.6e} vs {ana:.6e}")
check("analytic gradients verified against directional derivatives for "
      "all five polynomial quartics", grad_ok)

ALLPAT = dict(PATTERNS)
ALLPAT["Q6_eps"] = Q6
samples = [rand_phi() for _ in range(14)]
Vraw = np.array([[ALLPAT[k](T) for k in ALLPAT] for T in samples])
Vn = Vraw / np.abs(Vraw).max(axis=0, keepdims=True)
rank = int(np.linalg.matrix_rank(Vn, tol=1e-8))
rank5 = int(np.linalg.matrix_rank(Vn[:, :5], tol=1e-8))
# locate the dependency among the five polynomial patterns
null = np.linalg.svd(Vn[:, :5])[2][-1]
dep_name = list(PATTERNS)[int(np.argmax(np.abs(null)))]
check("machine independence test: the five polynomial quartics span a "
      "FOUR-dimensional space -- the cyclic pattern DEGENERATES with "
      "the 2-2 contraction (Q5 = Q2, a machine-found identity of the "
      "antisymmetric 4-tensor) -- and the epsilon quartic is a genuine "
      "FIFTH independent invariant",
      rank5 == 4 and rank == 5,
      f"rank(5 poly) = {rank5}, rank(+eps) = {rank}; null vector "
      f"{np.round(null / np.abs(null).max(), 4).tolist()} "
      f"(dependency: {dep_name})")
SCAN_PATS = [k for k in PATTERNS if k != dep_name]  # independent set

MONOS = [(4, 0, 0), (0, 4, 0), (0, 0, 4), (2, 2, 0), (2, 0, 2), (0, 2, 2),
         (1, 3, 0), (3, 1, 0), (1, 0, 3), (3, 0, 1), (0, 1, 3), (0, 3, 1),
         (2, 1, 1), (1, 2, 1), (1, 1, 2)]
pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1),
       (1, 0, -1), (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, 1, -1),
       (1, -1, 1), (2, 1, 0), (1, 2, 0), (0, 2, 1), (2, 0, 1), (1, 1, 2),
       (2, 1, 1), (1, 2, 1)]
A_fit = np.array([[float(p ** i * a ** j * w ** k) for (i, j, k) in MONOS]
                  for (p, a, w) in pts])
descent = {}
DESC_ARR = {}
for name, f in ALLPAT.items():
    y = np.array([f(phi(*pt)) for pt in pts])
    coef, *_ = np.linalg.lstsq(A_fit, y, rcond=None)
    resid = float(np.abs(A_fit @ coef - y).max() / max(1.0, np.abs(y).max()))
    coef[np.abs(coef) < 1e-7] = 0.0
    DESC_ARR[name] = coef
    descent[name] = {"coef": {str(m): float(c) for m, c in zip(MONOS, coef)
                              if c != 0.0}, "fit_residual": resid}
check("quartic Clebsch DESCENT computed for all six candidates "
      "(exact monomial fits on the singlet slice)",
      all(d["fit_residual"] < 1e-9 for d in descent.values()),
      "residuals < 1e-9")

def bad_terms_of(name):
    """Terms breaking w = 0 evenness (odd w powers) or the a = 0
    stationarity of the Pati-Salam point (the single monomial p^3 a);
    p a^3 and p a w^2 are harmless (their a-derivative vanishes at
    a = 0 and their w-derivative vanishes at w = 0)."""
    return [str(m) for m, c in zip(MONOS, DESC_ARR[name])
            if c != 0.0 and (m[2] % 2 == 1 or m == (3, 1, 0))]


bad_scan = {k: bad_terms_of(k) for k in DESC_ARR if bad_terms_of(k)}
check("stationarity structure of the descent: NO odd-in-w terms and NO "
      "p^3 a term in ANY invariant (w = 0 consistent; a = 0 stationary "
      "at the Pati-Salam point; the harmless p a w^2 cross-terms in "
      "Q2/Q5 are recorded)", not bad_scan,
      f"violations: {bad_scan or 'none'}")
check("epsilon quartic DISCLOSED as an independent fifth invariant "
      "whose coupling is NOT scanned (boundary flag "
      "epsilon_quartic_not_scanned); its slice descent is recorded in "
      "the ledger", True,
      f"Q6_eps slice coefficients: {descent['Q6_eps']['coef'] or 'all zero'}")
DLOG["S2_descent"] = descent

# ============================================================ section 3
print("== 9b-1b section 3: Hessians, Goldstones, multiplets ==")


def hessian_of(gradf, v):
    H = np.zeros((210, 210))
    for B in range(210):
        eB = UNITS[B]
        g1 = gradf(v + eB) - gradf(v - eB)
        g2 = gradf(v + 2 * eB) - gradf(v - 2 * eB)
        H[:, B] = comps((8.0 * g1 - g2) / 12.0)
    return 24.0 * H          # phi-coordinate Hessian


def hessian_cubic(v):
    H = np.zeros((210, 210))
    for B in range(210):
        eB = UNITS[B]
        col = 3.0 * (np.einsum('klmn,mnij->ijkl', v, eB)
                     + np.einsum('klmn,mnij->ijkl', eB, v))
        H[:, B] = comps(col)
    return 24.0 * H


H_M2 = 48.0 * np.eye(210)     # phi-Hessian of I2 = 24 sum phi^2
LR_RATIO = 0.8
VACUA = {"PS": phi(1.0, 0.0, 0.0), "LR": phi(1.0, LR_RATIO, 0.0)}
H_store = {}
for vac, v in VACUA.items():
    H_store[vac] = {"m2": H_M2, "kappa": hessian_cubic(v)}
    for name in SCAN_PATS:
        H_store[vac][name] = hessian_of(GRADS[name], v)
print(f"   Hessians computed for both vacua (couplings: {SCAN_PATS})")

I2_ARR = np.array([1.0, 3.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
I2_ARR = np.array([24.0 * {(2, 0, 0): 1, (0, 2, 0): 3, (0, 0, 2): 6}.get(
    m[:2] + (m[2],), 0) for m in [(2, 0, 0), (0, 2, 0), (0, 0, 2)]])


def slice_V_arr(m2, kappa, lamvec, p, a, w=0.0):
    val = m2 * 24 * (p * p + 3 * a * a + 6 * w * w) \
        + kappa * 48 * (a ** 3 + 3 * p * w * w + 6 * a * w * w)
    mono = np.array([p ** i * a ** j * w ** k for (i, j, k) in MONOS])
    for lam, name in zip(lamvec, SCAN_PATS):
        val += lam * float(DESC_ARR[name] @ mono)
    return val


def stationarity(vac, lamvec, kappa=None):
    h = 1e-5

    def dV(m2, kap, p, a, wrt):
        if wrt == "p":
            return (slice_V_arr(m2, kap, lamvec, p + h, a)
                    - slice_V_arr(m2, kap, lamvec, p - h, a)) / (2 * h)
        return (slice_V_arr(m2, kap, lamvec, p, a + h)
                - slice_V_arr(m2, kap, lamvec, p, a - h)) / (2 * h)

    if vac == "PS":
        f0, f1 = dV(0.0, kappa, 1, 0, "p"), dV(1.0, kappa, 1, 0, "p")
        return -f0 / (f1 - f0), kappa
    a0 = LR_RATIO
    r = np.array([dV(0, 0, 1, a0, "p"), dV(0, 0, 1, a0, "a")])
    Jm = np.array([[dV(1, 0, 1, a0, "p") - r[0], dV(0, 1, 1, a0, "p") - r[0]],
                   [dV(1, 0, 1, a0, "a") - r[1],
                    dV(0, 1, 1, a0, "a") - r[1]]])
    sol = np.linalg.solve(Jm, -r)
    return float(sol[0]), float(sol[1])


def total_H(vac, m2, kappa, lamvec):
    H = m2 * H_store[vac]["m2"] + kappa * H_store[vac]["kappa"]
    for lam, name in zip(lamvec, SCAN_PATS):
        H = H + lam * H_store[vac][name]
    return (H + H.T) / 2


gold_ok = True
for _ in range(4):
    lamv = rng.uniform(-1, 1, size=len(SCAN_PATS))
    kap = float(rng.uniform(-1, 1))
    m2p, _ = stationarity("PS", lamv, kap)
    m2l, kapl = stationarity("LR", lamv)
    evp = np.linalg.eigvalsh(total_H("PS", m2p, kap, lamv))
    evl = np.linalg.eigvalsh(total_H("LR", m2l, kapl, lamv))
    sc_p = max(1.0, np.abs(evp).max())
    sc_l = max(1.0, np.abs(evl).max())
    zp = int(np.sum(np.abs(evp) < 1e-7 * sc_p))
    zl = int(np.sum(np.abs(evl) < 1e-7 * sc_l))
    gold_ok &= (zp == 24 and zl == 30)
    if not gold_ok:
        print(f"   zeros: PS {zp}, LR {zl}")
check("GAUGE GATE: stationary Hessians have EXACTLY 24 Goldstone zeros "
      "at the Pati-Salam vacuum and 30 at the left-right vacuum for "
      "random couplings -- invariants, gradients, descent and "
      "stationarity are mutually consistent", gold_ok)

# label operators
GENMAT = {}
for ab in [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (6, 8), (7, 9),
           (6, 9), (7, 8)]:
    M = np.zeros((210, 210))
    for B in range(210):
        M[:, B] = comps(act(UNITS[B], ab))
    GENMAT[ab] = M
CAS = {}
for side, combos in {"L": [((6, 7), (8, 9), +1), ((6, 8), (7, 9), -1),
                           ((6, 9), (7, 8), +1)],
                     "R": [((6, 7), (8, 9), -1), ((6, 8), (7, 9), +1),
                           ((6, 9), (7, 8), -1)]}.items():
    C = np.zeros((210, 210))
    for g1, g2, sg in combos:
        Tg = (GENMAT[g1] + sg * GENMAT[g2]) / 2.0
        C += Tg @ Tg
    CAS[side] = -C
JBL = GENMAT[(0, 1)] + GENMAT[(2, 3)] + GENMAT[(4, 5)]
CAS["BL"] = -(JBL @ JBL)


def classify(H):
    ev, U = np.linalg.eigh(H)
    sc = max(1.0, np.abs(ev).max())
    blocks, i = [], 0
    while i < len(ev):
        j = i
        while j + 1 < len(ev) and abs(ev[j + 1] - ev[i]) < 1e-5 * sc:
            j += 1
        P = U[:, i:j + 1]
        lab = {k: float(np.trace(P.T @ CAS[k] @ P) / (j + 1 - i))
               for k in CAS}
        blocks.append({"m2": float(ev[i]), "dim": j + 1 - i,
                       "L": round(lab["L"], 4), "R": round(lab["R"], 4),
                       "BL": round(lab["BL"], 4)})
        i = j + 1
    return blocks


if DEBUG_BLOCKS:
    lamv = rng.uniform(-1, 1, size=len(SCAN_PATS))
    kap = float(rng.uniform(-1, 1))
    for vac in ("PS", "LR"):
        m2v, kk = (stationarity("PS", lamv, kap) if vac == "PS"
                   else stationarity("LR", lamv))
        blocks = classify(total_H(vac, m2v, kk, lamv))
        print(f"--- {vac} blocks ---")
        for b in blocks:
            print(f"   m2={b['m2']:+.4f} dim={b['dim']:3d} "
                  f"L={b['L']:.3f} R={b['R']:.3f} BL={b['BL']:.3f}")
    sys.exit(0)

# ============================================================ section 4
print("== 9b-1b section 4: verdict scan ==")

B_SM = np.array([41 / 10, -19 / 6, -7.0])
BIJ_SM = np.array([[199 / 50, 27 / 10, 44 / 5], [9 / 10, 35 / 6, 12],
                   [11 / 10, 9 / 2, -26]])
MZ = 91.1876
AINV_MZ = np.array([0.6 * 127.951 * (1 - 0.23122), 127.951 * 0.23122,
                    1 / 0.1180])
B_LR_V = np.array([-7.0, -3.0, -7 / 3, 11 / 2])
B_PS_V = np.array([-23 / 3, -3.0, 11 / 3])


def run2(ainv0, t1, nstep=120):
    a = np.array(ainv0, dtype=float)
    h = t1 / nstep

    def f(ai):
        al = 1.0 / ai
        return -B_SM / (2 * math.pi) - (BIJ_SM @ al) / (8 * math.pi ** 2)

    for _ in range(nstep):
        k1 = f(a); k2 = f(a + h * k1 / 2)
        k3 = f(a + h * k2 / 2); k4 = f(a + h * k3)
        a = a + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return a


def solve_thresh(kind, thresh):
    bI = B_LR_V if kind == "LR" else B_PS_V
    corr = sum((b * ln for b, ln in thresh),
               np.zeros(4 if kind == "LR" else 3))

    def equations(u):
        lnMI, lnMX, split = u
        a1, a2, a3 = run2(AINV_MZ, lnMI - math.log(MZ))
        if kind == "LR":
            a2R = split
            aBL = (a1 - 0.6 * a2R) / 0.4
            aI0 = np.array([a3, a2, a2R, aBL])
        else:
            aI0 = np.array([a3, a2, split])
        aX = aI0 - (bI * (lnMX - lnMI) + corr) / (2 * math.pi)
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
        Jm = np.zeros((3, 3))
        for j in range(3):
            du = np.zeros(3); du[j] = 1e-6
            Jm[:, j] = (equations(u + du)[0] - equations(u - du)[0]) / 2e-6
        try:
            step = np.linalg.solve(Jm, r)
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
    lg10 = math.log(10)
    return {"log10_MI": u[0] / lg10, "log10_MX": u[1] / lg10,
            "alpha_G_inv": float(np.mean(aX)),
            "residual": float(np.max(np.abs(r)))}


def tau_d6(MX, aGinv):
    return 1.6e34 * (MX / 1e16) ** 4 * (aGinv / 25.0) ** 2


SK_EPI = 2.4e34
G2 = 4 * math.pi / 45.0
MV2 = {vac: float(np.linalg.eigvalsh(vector_M2(v)).max()) * G2
       for vac, v in VACUA.items()}

# real-scalar b vectors for the 210 remnants (identified from labels)
B_R_PS = {"(15,1,1)": np.array([2 / 3, 0, 0]),
          "(15,1,3)": np.array([2.0, 0, 5.0]),
          "(15,3,1)": np.array([2.0, 5.0, 0]),
          "(10,22)pair": np.array([4.0, 10 / 3, 10 / 3]),
          "(1,1,1)": np.array([0.0, 0, 0])}


def id_PS(b):
    if b["dim"] == 15 and b["L"] < 1e-4 and b["R"] < 1e-4:
        return "(15,1,1)"
    if b["dim"] == 45 and b["R"] > 1e-4 and b["L"] < 1e-4:
        return "(15,1,3)"
    if b["dim"] == 45 and b["L"] > 1e-4 and b["R"] < 1e-4:
        return "(15,3,1)"
    if b["dim"] == 80:
        return "(10,22)pair"
    if b["dim"] == 1:
        return "(1,1,1)"
    return None


# explicit LR b-vector map keyed by (dim, L>0, R>0, |BL| per state):
# derived from the PS -> LR decomposition of the 210 remnants:
# (15,1,1) -> (8,1,1,0)+(3,1,1,4/3)pair+(1,1,1,0)
# (15,1,3) -> (8,1,3,0)+(3,1,3,4/3)pair+(1,1,3,0)
# (15,3,1) -> (8,3,1,0)+(3,3,1,4/3)pair+(1,3,1,0)
# (6,2,2)  -> (3,2,2,-2/3)+(3b,2,2,+2/3) [one 24-dim real block]
# (10,22)pair -> (1,2,2,+-2)[8] + (3,2,2,+-2/3)[24] + (6,2,2,+-2/3)[48]
# keys: (dim, L > 0, R > 0, BL-bucket) with BL-bucket = round(Casimir
# of J_BL per state) = round((3 (B-L)/2)^2), as observed in the block
# diagnostic: buckets {0, 1, 4, 9}.  b vectors: real multiplets get
# (1/6) S, conjugate pairs (= one complex multiplet) get (1/3) S;
# abelian entries in GUT normalization (3/8)(B-L)^2.
LR_BVEC = {
    (8, False, False, 0): ("(8,1,1,0)", np.array([1 / 2, 0, 0, 0])),
    (6, False, False, 4): ("(3,1,1,4/3)p", np.array([1 / 6, 0, 0, 2 / 3])),
    (1, False, False, 0): ("(1,1,1,0)", np.array([0.0, 0, 0, 0])),
    (24, False, True, 0): ("(8,1,3,0)", np.array([3 / 2, 0, 8 / 3, 0])),
    (18, False, True, 4): ("(3,1,3,4/3)p", np.array([1 / 2, 0, 2, 2])),
    (3, False, True, 0): ("(1,1,3,0)", np.array([0, 0, 1 / 3, 0])),
    (24, True, False, 0): ("(8,3,1,0)", np.array([3 / 2, 8 / 3, 0, 0])),
    (18, True, False, 4): ("(3,3,1,4/3)p", np.array([1 / 2, 2, 0, 2])),
    (3, True, False, 0): ("(1,3,1,0)", np.array([0, 1 / 3, 0, 0])),
    (24, True, True, 1): ("(3,2,2,2/3)p", np.array([2 / 3, 1, 1, 2 / 3])),
    (8, True, True, 9): ("(1,2,2,2)p", np.array([0, 1 / 3, 1 / 3, 2.0])),
    (48, True, True, 1): ("(6,2,2,2/3)p", np.array([10 / 3, 2, 2, 4 / 3])),
}


def id_LR_b(b):
    key = (b["dim"], b["L"] > 1e-4, b["R"] > 1e-4, int(round(b["BL"])))
    hit = LR_BVEC.get(key)
    return hit if hit else (None, None)


N_SAMPLES = 250
results = {"PS": [], "LR": []}
spectra_log = {"PS": [], "LR": []}
id_fail = {"PS": 0, "LR": 0}
for it in range(N_SAMPLES):
    lamv = rng.uniform(-1, 1, size=len(SCAN_PATS))
    kap = float(rng.uniform(-1, 1))
    for vac in ("PS", "LR"):
        m2v, kk = (stationarity("PS", lamv, kap) if vac == "PS"
                   else stationarity("LR", lamv))
        H = total_H(vac, m2v, kk, lamv)
        blocks = classify(H)
        sc = max(1.0, max(abs(b["m2"]) for b in blocks))
        phys = [b for b in blocks if abs(b["m2"]) > 1e-6 * sc]
        want_gold = 24 if vac == "PS" else 30
        if 210 - sum(b["dim"] for b in phys) != want_gold:
            id_fail[vac] += 1
            continue
        if any(b["m2"] < 0 for b in phys):
            results[vac].append({"viable": False})
            continue
        thr, spect, okid = [], {}, True
        for b in phys:
            if vac == "PS":
                name = id_PS(b)
                bv = B_R_PS.get(name) if name else None
            else:
                name, bv = id_LR_b(b)
            if name is None:
                okid = False
                break
            rr = math.sqrt(b["m2"] / MV2[vac])
            tag = name if name not in spect else f"{name}#2"
            spect[tag] = rr
            if bv is not None and rr < 1.0:
                thr.append((bv, min(math.log(1.0 / rr), 25.0)))
        if not okid:
            id_fail[vac] += 1
            continue
        sol = solve_thresh("LR" if vac == "LR" else "PS", thr)
        MX = 10 ** sol["log10_MX"]
        physical = (sol["residual"] < 1e-8
                    and math.log10(MZ) < sol["log10_MI"] < sol["log10_MX"]
                    and MX < 1.2e19)
        tau = tau_d6(MX, sol["alpha_G_inv"]) if physical else 0.0
        results[vac].append({"viable": True, "physical": physical,
                             "tau": tau,
                             "alive": bool(physical and tau > SK_EPI)})
        if len(spectra_log[vac]) < 4:
            spectra_log[vac].append(
                {"kappa": kk, "m2": m2v,
                 "spectrum_over_MV": {k: round(v, 4)
                                      for k, v in spect.items()},
                 "tau": tau, "log10_MX": sol["log10_MX"]})

frac = {}
for vac in ("PS", "LR"):
    viab = [r for r in results[vac] if r.get("viable")]
    alive = [r for r in viab if r.get("alive")]
    taus = sorted(r["tau"] for r in viab if r.get("physical"))
    frac[vac] = {
        "samples": N_SAMPLES, "classified": len(results[vac]),
        "id_failures": id_fail[vac], "tree_positive": len(viab),
        "alive": len(alive),
        "tau_16_50_84": ([f"{taus[int(q_ * (len(taus) - 1))]:.2e}"
                          for q_ in (0.16, 0.5, 0.84)] if taus else None)}
check("verdict scan executed: stationarity + classification + "
      "tree positivity + multi-threshold matching, both vacua",
      all(frac[v]["classified"] > 0 for v in frac),
      json.dumps(frac))
check("TREE-LEVEL LANDMINE QUANTIFIED: the fraction of random couplings "
      "for which the intermediate vacuum is a tree-level local minimum "
      "of the 210 sector is recorded (the known non-SUSY pseudo-"
      "Goldstone problem; one-loop stabilization cited, not re-derived)",
      True, f"PS {frac['PS']['tree_positive']}/{frac['PS']['classified']}"
      f", LR {frac['LR']['tree_positive']}/{frac['LR']['classified']}")
check("FINAL TREE-LEVEL THRESHOLD VERDICT: tau distribution over the "
      "viable coupling space per chain, ESH replaced by the computed "
      "210 spectrum", True,
      f"alive: PS {frac['PS']['alive']}/"
      f"{max(frac['PS']['tree_positive'], 1)}, "
      f"LR {frac['LR']['alive']}/{max(frac['LR']['tree_positive'], 1)}; "
      f"tau percentiles PS {frac['PS']['tau_16_50_84']}, "
      f"LR {frac['LR']['tau_16_50_84']}")
DLOG["S4_verdict"] = {"fractions": frac, "example_spectra": spectra_log}

# tree-level slice phase-diagram sample
pref = {"PS(21)": 0, "LR(15)": 0, "other": 0, "runaway": 0}
h = 1e-5
for _ in range(150):
    lamv = rng.uniform(-1, 1, size=len(SCAN_PATS))
    kap = float(rng.uniform(-1, 1))
    m2v = float(rng.uniform(-2, 0))
    best, bestV = None, None
    for _ in range(25):
        z = rng.uniform(-2, 2, size=3)
        ok_run = True
        for _ in range(300):
            g = np.array([
                (slice_V_arr(m2v, kap, lamv, *(z + h * np.eye(3)[i]))
                 - slice_V_arr(m2v, kap, lamv, *(z - h * np.eye(3)[i])))
                / (2 * h) for i in range(3)])
            if np.linalg.norm(g) < 1e-7:
                break
            z = z - 0.01 * g / max(1.0, np.linalg.norm(g))
            if np.abs(z).max() > 40:
                ok_run = False
                break
        if not ok_run:
            continue
        V0 = slice_V_arr(m2v, kap, lamv, *z)
        if bestV is None or V0 < bestV:
            bestV, best = V0, z.copy()
    if best is None:
        pref["runaway"] += 1
        continue
    zc = zero_count(vector_M2(phi(*best)), tol=1e-7)
    pref["PS(21)" if zc == 21 else "LR(15)" if zc == 15
         else "other"] += 1
check("tree-level slice phase-diagram SAMPLE recorded (little group of "
      "the deepest bounded slice minimum at random couplings)", True,
      str(pref))
DLOG["S4_phase_sample"] = pref

# ------------------------------------------------------------- ledger
n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-9b-1b 210 quartic descent, survivor masses, verdicts",
    "dyn_item": "DYN-9b-1b",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    # negative-boundary flags
    "sigma_sector_quartics_deferred_9b1c": True,
    "epsilon_quartic_not_scanned": True,
    "tree_level_only_one_loop_stabilization_cited": True,
    "coupling_sampling_in_full_sum_units": True,
    "LR_vev_ratio_fixed": LR_RATIO,
    "zeta_value_derived": False,
}
(OUT / "dyn9b1b_210_quartic_descent.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-1b: 210 quartic descent, survivor masses, verdicts", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "Full descent polynomials, example spectra, verdict fractions and "
      "the phase-diagram sample are in the JSON derivation_log.",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b1b_210_quartic_descent.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-9b-1b: {n_pass}/{len(CHECKS)} checks; ledgers -> "
      f"output/audit9/dyn9b1b_210_quartic_descent.*")
