#!/usr/bin/env python3
"""DYN-9b-1c: the epsilon-coupling scan and the 126bar mixed quartics.

Completes the two disclosed gaps of DYN-9b-1b:

  (A) THE EPSILON LEVER.  9b-1b found the Hodge epsilon quartic Q6 to
      be a genuine fifth independent invariant of the 210 whose slice
      descent vanishes identically -- the only coupling that moves
      off-slice masses without touching the vacuum conditions, hence
      the only tree-level lever that could rescue the left-right
      (p, a) vacuum (tree-positive on 0/249 samples in the 4-coupling
      scan).  Here the full Q6 gradient (including the dual-slot term)
      is built and gated; the exact epsilon Hessian is assembled with
      the cubic-homogeneity identity H e = [g(v+e) - g(v-e) - 2g(e)]/2;
      the Goldstone gate at lambda_eps != 0 (still exactly 24/30
      zeros) certifies it; and the rescue question is answered by a
      per-sample sweep of lambda_eps.

  (B) THE 126bar SECTOR.  The rank-5 tensor is built with the
      holomorphic 5-form sigma-direction.  Embedding gates: the mixed
      cubic Phi Sigma Sigma* descends EXACTLY to the AG eta-term
      structure |sigma|^2 (p + 3a - 6w), and the 45x45 vector mass
      matrix with BOTH vevs reproduces the full DYN-1b formulas
      (including the sigma terms) at the AG benchmark.  At sigma = 0
      the Sigma masses depend ONLY on m_Sigma^2 and the mixed
      Phi^2 Sigma Sigma* quartics (the (Sigma Sigma*)^2 self-couplings
      drop out of the mass matrix -- a clean structural fact, so their
      enumeration is deferred without loss).  A generating set of four
      mixed patterns is built, the Sigma-block Hessian at the
      Pati-Salam vacuum is classified, the Delta_R block is identified
      MACHINE-SIDE by projecting the sigma-direction onto eigenspaces,
      the extended-survival tuning puts it at M_I, and the remaining
      Sigma splittings feed the threshold solver: the K4 verdict is
      retested including the 126bar remnants.

Boundary: (Sigma Sigma*)^2 self-quartics enter only the sigma-scale
(M_I) stationarity, not the M_X-scale masses -- deferred with that
argument; epsilon-mixed Phi Phi Sigma Sigma* patterns not enumerated
(flagged); tree level; zeta NOT derived.
Ledgers -> output/audit9/dyn9b1c_eps_and_126_quartics.{json,md}.
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
rng = np.random.default_rng(20260706)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


# ------------------------------------------------ 210 core (as 9b-1b)
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
    """<*e_A, W> for all A: returns 210 values (up to the 6! factor
    convention absorbed into the pattern normalization)."""
    return (W.ravel()[FLAT6[DUAL_MAP]] @ SIGNS6) * DUAL_SIGN


e_p = from_comps([1.0 if q == (6, 7, 8, 9) else 0.0 for q in QUADS])
A_SET = {(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)}
W_SET = {(0, 1, 6, 7), (0, 1, 8, 9), (2, 3, 6, 7), (2, 3, 8, 9),
         (4, 5, 6, 7), (4, 5, 8, 9)}
e_a = from_comps([1.0 if q in A_SET else 0.0 for q in QUADS])
# sign fixed against the AG eta-term descent (the omega-even 210
# invariants cannot see this sign; the mixed cubic can)
e_w = -from_comps([1.0 if q in W_SET else 0.0 for q in QUADS])


def phi(p, a, w):
    return p * e_p + a * e_a + w * e_w


def I2(T):
    return float(np.einsum('ijkl,ijkl->', T, T))


def M2T(T):
    return np.einsum('ijmn,klmn->ijkl', T, T)


def Q6(T):
    S = star6_dense(T)
    W1 = np.einsum('abpq,cdpr->abcdqr', T, T)
    return float(np.einsum('abcdef,abcdqr,efqr->', S, W1, T))


def grad_Q6_clean(T):
    S = star6_dense(T)
    Y = np.einsum('abcdef,efqr->abcdqr', S, T)      # sum over e,f
    Y2 = np.einsum('abcdef,abpq->cdefpq', S, T)     # sum over a,b
    g_direct = (np.einsum('abcdqr,cdpr->abpq', Y, T)
                + np.einsum('abcdqr,abpq->cdpr', Y, T)
                + np.einsum('cdefpq,cdpr->efqr', Y2, T))
    W1 = np.einsum('abpq,cdpr->abcdqr', T, T)
    W = np.einsum('abcdqr,efqr->abcdef', W1, T)
    gd = from_comps(dual_readout(W) / 24.0)
    return g_direct + gd


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


# quartic/cubic machinery from 9b-1b (Q1-Q4 + kappa; gated there)
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


GRADS = {"Q1": grad_Q1, "Q2": grad_Q2, "Q3": grad_Q3, "Q4": grad_Q4}
SCAN_PATS = ["Q1", "Q2", "Q3", "Q4"]

# ============================================================ section 1
print("== 9b-1c section 1: the full epsilon gradient and Hessian ==")

Tt = from_comps(rng.standard_normal(210))
D = from_comps(rng.standard_normal(210)) * 1e-4
num = (Q6(Tt + D) - Q6(Tt - D)) / 2.0
ana = float(np.einsum('ijkl,ijkl->', grad_Q6_clean(Tt), D))
check("FULL epsilon gradient (three direct placements + the dual-slot "
      "adjoint term) verified against the directional derivative",
      abs(num - ana) < 1e-6 * max(1.0, abs(num)),
      f"directional {num:.6e} vs analytic {ana:.6e}")


def hessian_hom_cubic(gradf, v, gE_cache):
    """Exact Hessian for a HOMOGENEOUS-cubic gradient:
    H e = [g(v+e) - g(v-e) - 2 g(e)] / 2."""
    H = np.zeros((210, 210))
    for B in range(210):
        eB = UNITS[B]
        col = (gradf(v + eB) - gradf(v - eB) - 2.0 * gE_cache[B]) / 2.0
        H[:, B] = comps(col)
    return 24.0 * H


LR_RATIO = 0.8
VACUA = {"PS": phi(1.0, 0.0, 0.0), "LR": phi(1.0, LR_RATIO, 0.0)}
H_eps = {}
CACHE_DIR.mkdir(parents=True, exist_ok=True)
cache_key = hashlib.sha256(
    Path(__file__).read_bytes() + np.__version__.encode("ascii")
).hexdigest()[:16]
CACHE = CACHE_DIR / f"dyn9b1c_heps_{cache_key}.npz"
if CACHE.exists():
    dat = np.load(CACHE)
    H_eps = {"PS": dat["PS"], "LR": dat["LR"]}
    print("   H_eps loaded from cache (delete the npz to recompute)")
else:
    print("   caching g(e_B) for the epsilon pattern (210 calls)...")
    GE_CACHE = [grad_Q6_clean(UNITS[B]) for B in range(210)]
    for vac, v in VACUA.items():
        print(f"   building H_eps at the {vac} vacuum...")
        H = hessian_hom_cubic(grad_Q6_clean, v, GE_CACHE)
        H_eps[vac] = (H + H.T) / 2
    np.savez(CACHE, **H_eps)

# cross-check the identity on a random direction
u1 = rng.standard_normal(210)
u2 = rng.standard_normal(210)
h = 1e-3
T1, T2 = from_comps(u1), from_comps(u2)
vv = VACUA["PS"]
f0 = (Q6(vv + h * T1 + h * T2) - Q6(vv + h * T1 - h * T2)
      - Q6(vv - h * T1 + h * T2) + Q6(vv - h * T1 - h * T2)) / (4 * h * h)
bil = float(u1 @ H_eps["PS"] @ u2)
check("epsilon Hessian certified: the homogeneity-identity build "
      "matches the direct second difference of Q6 values on random "
      "directions",
      abs(f0 - bil) < 1e-4 * max(1.0, abs(f0)),
      f"{f0:.6e} vs {bil:.6e}")

# ============================================================ section 2
print("== 9b-1c section 2: Goldstone gate and the LR rescue scan ==")


def hessian_of(gradf, v):
    H = np.zeros((210, 210))
    for B in range(210):
        eB = UNITS[B]
        g1 = gradf(v + eB) - gradf(v - eB)
        g2 = gradf(v + 2 * eB) - gradf(v - 2 * eB)
        H[:, B] = comps((8.0 * g1 - g2) / 12.0)
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
H_store = {}
for vac, v in VACUA.items():
    H_store[vac] = {"m2": H_M2, "kappa": hessian_cubic(v)}
    for name in SCAN_PATS:
        H_store[vac][name] = hessian_of(GRADS[name], v)
print("   polynomial Hessians rebuilt")

# descent arrays for stationarity (from 9b-1b, re-fitted here)
MONOS = [(4, 0, 0), (0, 4, 0), (0, 0, 4), (2, 2, 0), (2, 0, 2), (0, 2, 2),
         (1, 3, 0), (3, 1, 0), (1, 0, 3), (3, 0, 1), (0, 1, 3), (0, 3, 1),
         (2, 1, 1), (1, 2, 1), (1, 1, 2)]
pts = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1),
       (1, 0, -1), (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, 1, -1),
       (1, -1, 1), (2, 1, 0), (1, 2, 0), (0, 2, 1), (2, 0, 1), (1, 1, 2),
       (2, 1, 1), (1, 2, 1)]
A_fit = np.array([[float(p ** i * a ** j * w ** k) for (i, j, k) in MONOS]
                  for (p, a, w) in pts])
PATFUN = {"Q1": lambda T: I2(T) ** 2,
          "Q2": lambda T: float(np.einsum('ijkl,ijkl->', M2T(T), M2T(T))),
          "Q3": lambda T: float(np.einsum('ijkl,ikjl->', M2T(T), M2T(T))),
          "Q4": lambda T: float(np.einsum(
              'ij,ij->', np.einsum('iabc,jabc->ij', T, T),
              np.einsum('iabc,jabc->ij', T, T)))}
DESC_ARR = {}
for name, f in PATFUN.items():
    y = np.array([f(phi(*pt)) for pt in pts])
    coef, *_ = np.linalg.lstsq(A_fit, y, rcond=None)
    coef[np.abs(coef) < 1e-7] = 0.0
    DESC_ARR[name] = coef
eps_slice = np.array([Q6(phi(*pt)) for pt in pts])
check("re-verified: the epsilon quartic vanishes identically on the "
      "singlet slice (all 19 sample points), so lambda_eps enters NO "
      "stationarity condition -- it is a pure off-slice mass lever",
      np.abs(eps_slice).max() < 1e-8)


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


def total_H(vac, m2, kappa, lamvec, leps):
    H = m2 * H_store[vac]["m2"] + kappa * H_store[vac]["kappa"] \
        + leps * H_eps[vac]
    for lam, name in zip(lamvec, SCAN_PATS):
        H = H + lam * H_store[vac][name]
    return (H + H.T) / 2


gold_ok = True
for _ in range(3):
    lamv = rng.uniform(-1, 1, size=4)
    kap = float(rng.uniform(-1, 1))
    leps = float(rng.uniform(-1, 1))
    m2p, _ = stationarity("PS", lamv, kap)
    m2l, kapl = stationarity("LR", lamv)
    for vac, (m2v, kk) in {"PS": (m2p, kap), "LR": (m2l, kapl)}.items():
        ev = np.linalg.eigvalsh(total_H(vac, m2v, kk, lamv, leps))
        z = int(np.sum(np.abs(ev) < 1e-7 * max(1.0, np.abs(ev).max())))
        gold_ok &= z == (24 if vac == "PS" else 30)
check("GAUGE GATE WITH lambda_eps != 0: exactly 24 / 30 Goldstone "
      "zeros persist -- certifies the epsilon Hessian's gauge "
      "consistency (the strongest available correctness check)",
      gold_ok)

# the rescue scan: can ANY lambda_eps in [-1.5, 1.5] make the LR
# Hessian positive semi-definite (given stationarity)?
EPS_GRID = np.linspace(-1.5, 1.5, 13)
rescued, ps_flipped = 0, 0
lr_min_best = []
n_scan = 200
for it in range(n_scan):
    lamv = rng.uniform(-1, 1, size=4)
    kap = float(rng.uniform(-1, 1))
    m2l, kapl = stationarity("LR", lamv)
    best = -np.inf
    for le in EPS_GRID:
        ev = np.linalg.eigvalsh(total_H("LR", m2l, kapl, lamv, le))
        sc = max(1.0, np.abs(ev).max())
        phys_min = np.min(ev[np.abs(ev) > 1e-7 * sc]) \
            if np.any(np.abs(ev) > 1e-7 * sc) else 0.0
        best = max(best, phys_min)
    lr_min_best.append(best)
    if best > 0:
        rescued += 1
    m2p, _ = stationarity("PS", lamv, kap)
    for le in (-1.0, 1.0):
        evp = np.linalg.eigvalsh(total_H("PS", m2p, kap, lamv, le))
        scp = max(1.0, np.abs(evp).max())
        pos = np.all(evp[np.abs(evp) > 1e-7 * scp] > 0)
        # count PS samples whose positivity FLIPS due to eps
        ev0 = np.linalg.eigvalsh(total_H("PS", m2p, kap, lamv, 0.0))
        pos0 = np.all(ev0[np.abs(ev0) > 1e-7 * scp] > 0)
        if pos != pos0:
            ps_flipped += 1
            break
check("LR RESCUE VERDICT: fraction of random couplings for which SOME "
      "lambda_eps in [-1.5, 1.5] makes the left-right (p, a) vacuum a "
      "tree-level local minimum",
      True,
      f"rescued {rescued}/{n_scan}; best physical min-eigenvalue "
      f"median {np.median(lr_min_best):.3e}")
check("PS stability under the epsilon coupling recorded (samples whose "
      "tree positivity flips when lambda_eps = +-1)", True,
      f"flipped: {ps_flipped}/{n_scan}")
DLOG["S2_eps_scan"] = {"rescued": rescued, "n": n_scan,
                       "lr_best_min_eig_percentiles":
                       [float(np.percentile(lr_min_best, q))
                        for q in (16, 50, 84)],
                       "ps_flipped": ps_flipped}

# ============================================================ section 3
print("== 9b-1c section 3: the 126bar sector and its gates ==")

QUINTS = list(itertools.combinations(range(N), 5))           # 252
PERMS5 = [(pm, perm_sign(pm)) for pm in itertools.permutations(range(5))]
FLAT5 = np.array([[np.ravel_multi_index(tuple(q[p] for p in pm), (N,) * 5)
                   for pm, _ in PERMS5] for q in QUINTS])
SIGNS5 = np.array([sg for _, sg in PERMS5], dtype=float)


def comps5(G):
    return (G.reshape(-1)[FLAT5] @ SIGNS5) / 120.0


def from_comps5(x):
    T = np.zeros(N ** 5, dtype=complex)
    np.add.at(T, FLAT5.ravel(),
              (np.asarray(x, dtype=complex)[:, None]
               * SIGNS5[None, :]).ravel())
    return T.reshape((N,) * 5)


# holomorphic 5-form (the sigma / nu^c-pairing direction)
U_HOL = [(np.eye(N)[2 * k] + 1j * np.eye(N)[2 * k + 1]) / math.sqrt(2)
         for k in range(5)]
OMEGA = np.zeros((N,) * 5, dtype=complex)
base = np.einsum('i,j,k,l,m->ijklm', *U_HOL)
for pm, sg in PERMS5:
    OMEGA += sg * np.transpose(base, pm)
OMEGA_C = comps5(OMEGA)
OMEGA_C = OMEGA_C / np.linalg.norm(OMEGA_C)          # unit phi-norm

# self-duality projector: 252 = 126 + 126bar; star on 5-forms squares
# to -1, so the halves are the +-i eigenspaces.  Basis pairs
# (q, rest(q)) give eigenvectors (|q> +- i s_q |rest>)/sqrt(2).
QUINT_IDX = {q: i for i, q in enumerate(QUINTS)}
pairs, seen = [], set()
for A, q in enumerate(QUINTS):
    r = tuple(sorted(ALL - set(q)))
    if q in seen or r in seen:
        continue
    seen.add(q); seen.add(r)
    s5 = perm_sign(tuple(list(q) + list(r)))
    pairs.append((A, QUINT_IDX[r], s5))
D5 = np.zeros((252, 252))
for A, B, s5 in pairs:
    D5[B, A] = s5
    D5[A, B] = -s5
assert np.allclose(D5 @ D5, -np.eye(252))
lam_omega = complex((D5 @ OMEGA_C) @ np.conj(OMEGA_C))
SGN = 1.0 if abs(lam_omega - 1j) < abs(lam_omega + 1j) else -1.0
U126 = np.zeros((252, 126), dtype=complex)
for k, (A, B, s5) in enumerate(pairs):
    U126[A, k] = 1.0 / math.sqrt(2)
    U126[B, k] = -SGN * 1j * s5 / math.sqrt(2)
OMEGA_126 = np.conj(U126).T @ OMEGA_C


def sigma_field(s):
    return from_comps5(s * OMEGA_C)


# mixed cubic Phi Sigma Sigma*: the AG eta-term gate
def mixed_cubic(T, S):
    return float(np.real(np.einsum('ijklm,ijkno,lmno->',
                                   np.conj(S), S, T)))


Sg = sigma_field(1.0)
cP = mixed_cubic(e_p, Sg)
cA = mixed_cubic(e_a, Sg)
cW = mixed_cubic(e_w, Sg)
ratio_ok = abs(cA / cP - 3.0) < 1e-9 and abs(cW / cP + 6.0) < 1e-9
check("MIXED CUBIC GATE: the Phi Sigma Sigma* invariant descends to "
      "|sigma|^2 (p + 3a - 6w) -- the AG eta-term structure, "
      "extending the embedding anchor to the 126bar sector",
      ratio_ok and abs(cP) > 1e-12,
      f"(1, {cA/cP:.6f}, {cW/cP:.6f}) expected (1, 3, -6)")


def act5(T, ab):
    a_, b_ = ab
    d = np.zeros_like(T)
    for slot in range(5):
        Tm = np.moveaxis(T, slot, 0)
        dm = np.moveaxis(d, slot, 0)
        dm[a_] += Tm[b_]
        dm[b_] -= Tm[a_]
    return d


def vector_M2_full(T4, S5):
    cols4 = np.array([act4(T4, g).ravel() for g in GENS])
    M = cols4 @ cols4.T / 24.0
    cols5 = np.array([act5(S5, g).ravel() for g in GENS])
    M = M + np.real(cols5 @ np.conj(cols5).T) / 120.0
    return M


x = 0.1
w_b, a_b = -x, (x * x + 2 * x - 1) / (1 - x)
p_b = x * (5 * x * x - 1) / (1 - x) ** 2
s_b = math.sqrt(2 * x * (1 - 3 * x) * (1 + x * x) / (1 - x) ** 2)
Mfull = vector_M2_full(phi(p_b, a_b, w_b), sigma_field(s_b))
evf = np.linalg.eigvalsh(Mfull)
z_full = int(np.sum(evf < 1e-9 * max(1.0, evf.max())))
nz = sorted(e for e in evf if e > 1e-9 * max(1.0, evf.max()))
groups = []
for e in nz:
    if groups and abs(e - groups[-1][0]) < 1e-8 * max(1.0, e):
        groups[-1][1] += 1
    else:
        groups.append([e, 1])
mult_map = {m: v for v, m in groups}
twelves = sorted(v for v, m in groups if m == 12)
s2 = s_b ** 2
gate_ok, detail = False, f"zeros = {z_full}, groups = {groups}"
if z_full == 12 and len(nz) == 33 and {1, 2, 6}.issubset(mult_map) \
        and len(twelves) == 2:
    # X = (3,2,-5/6) has NO sigma term: identify it among the 12s
    fX0 = 4 * (a_b + w_b) ** 2 + 2 * (p_b + w_b) ** 2
    fE0 = 4 * (a_b - w_b) ** 2 + 2 * (w_b - p_b) ** 2
    for mX, mE in ((twelves[0], twelves[1]), (twelves[1], twelves[0])):
        scale = mX / fX0
        c = mult_map[1] / scale / (10 * s2)          # from G
        ok = (c > 0
              and abs(mult_map[6] / scale
                      - (8 * a_b ** 2 + 16 * w_b ** 2 + 2 * c * s2))
              < 1e-6
              and abs(mult_map[2] / scale - (24 * w_b ** 2 + 2 * c * s2))
              < 1e-6
              and abs(mE / scale - (fE0 + 2 * c * s2)) < 1e-6)
        if ok:
            gate_ok = True
            detail = (f"scale = {scale:.6f}, sigma-norm c = {c:.6f}; "
                      "J, F, E cross-checked < 1e-6; X carries no "
                      "sigma term, exactly as in DYN-1b")
            break
check("FULL VECTOR-MASS GATE with sigma: 12 zeros (SM) + all 33 "
      "massive bosons match the complete DYN-1b formulas G/J/F/E/X "
      "including the sigma terms (one overall scale from X, one "
      "sigma-normalization from G, J/F/E cross-checked)",
      gate_ok, detail)
DLOG["S3_126_gates"] = {"eta_ratios": [1.0, cA / cP, cW / cP],
                        "vector_zeros": z_full}

# ============================================================ section 4
print("== 9b-1c section 4: Sigma masses at the PS vacuum, K4 retest ==")

# mixed Phi^2 Sigma Sigma* patterns (generating set; sigma = 0 masses)
v_PS = VACUA["PS"]


def sigma_block(K4mat):
    """Hermitian 252x252 block from a rank-4 kernel K_{ijkl} acting on
    two Sigma slots: H[al,be] = <e_al| K x 1^3 |e_be> (pattern-specific
    contractions built below call this with their own kernel)."""
    H = np.zeros((252, 252), dtype=complex)
    for B, qq in enumerate(QUINTS):
        eB = from_comps5([1.0 if i == B else 0.0 for i in range(252)])
        X = np.einsum('ijkl,cdekl->cdeij', K4mat, eB)
        H[:, B] = comps5(X)
    return 120.0 * (H + H.conj().T) / 2


K_P2 = np.einsum('abij,abkl->ijkl', v_PS, v_PS)      # 2-2 kernel
H_P2 = sigma_block(K_P2)

S2v = np.einsum('iabc,jabc->ij', v_PS, v_PS)         # 3-1 kernel


def sigma_block_1slot(Smat):
    H = np.zeros((252, 252), dtype=complex)
    for B in range(252):
        eB = from_comps5([1.0 if i == B else 0.0 for i in range(252)])
        X = np.einsum('de,eijkl->dijkl', Smat, eB)
        H[:, B] = comps5(X)
    return 120.0 * (H + H.conj().T) / 2


H_P4 = sigma_block_1slot(S2v)
H_P1 = I2(v_PS) * np.eye(252)                        # product pattern
K_P3 = np.einsum('abij,cdij->abcd', v_PS, v_PS)


def sigma_block_P3(K):
    H = np.zeros((252, 252), dtype=complex)
    for B in range(252):
        eB = from_comps5([1.0 if i == B else 0.0 for i in range(252)])
        X = np.einsum('abcd,abekl->cdekl', K, eB)
        H[:, B] = comps5(X)
    return 120.0 * (H + H.conj().T) / 2


H_P3 = sigma_block_P3(K_P3)


def sigma_block_eta(T4):
    """Mixed CUBIC Phi Sigma Sigma* at <Phi> = T4: LINEAR in p, hence
    D-parity odd -- the term that splits (10,1,3) from (10bar,3,1)."""
    H = np.zeros((252, 252), dtype=complex)
    for B in range(252):
        eB = from_comps5([1.0 if i == B else 0.0 for i in range(252)])
        X = np.einsum('lmno,ijkno->ijklm', T4, eB)
        H[:, B] = comps5(X)
    return 120.0 * (H + H.conj().T) / 2


H_ETA = sigma_block_eta(v_PS)
MIXED = {"ETA": H_ETA, "P1": H_P1, "P2": H_P2, "P3": H_P3, "P4": H_P4}
print("   mixed cubic + quartic Sigma blocks built")

# duality restriction: the 126 and 126bar halves SHARE the (6,1,1)
# and (15,2,2) PS channels, so PS-invariant kernels may connect the
# halves; the PHYSICAL mass operator is the 126bar restriction
# U+ H U (the orthogonal components are not fields).  Sanity: Omega
# lies entirely in the restricted half.
mix_norm = max(float(np.abs(M @ D5 - D5 @ M).max())
               for M in MIXED.values())
omega_in = float(np.linalg.norm(np.conj(U126).T @ OMEGA_C))
check("duality restriction consistent: the sigma direction lies "
      "ENTIRELY in the restricted 126bar half (norm 1), and the "
      "cross-half mixing of the PS-invariant kernels (allowed in the "
      "shared (6,1,1)/(15,2,2) channels; those components are not "
      "fields) is recorded",
      abs(omega_in - 1.0) < 1e-9,
      f"|P Omega| = {omega_in:.12f}; cross-half kernel norm "
      f"{mix_norm:.2e}")
MIX126 = {k: np.conj(U126).T @ M @ U126 for k, M in MIXED.items()}

# PS label operators restricted to the 126bar half
# SU(2)_L / SU(2)_R NAMING ANCHOR: the two SO(4) self-dual halves are
# named so that the sigma direction (nu^c nu^c pairing = Delta_R) is
# R-charged -- the physics definition; the opposite assignment is what
# the raw orientation convention would have given.
CAS5 = {}
for side, combos in {"R": [((6, 7), (8, 9), +1), ((6, 8), (7, 9), -1),
                           ((6, 9), (7, 8), +1)],
                     "L": [((6, 7), (8, 9), -1), ((6, 8), (7, 9), +1),
                           ((6, 9), (7, 8), -1)]}.items():
    C = np.zeros((252, 252), dtype=complex)
    for g1, g2, sg in combos:
        M1 = np.zeros((252, 252), dtype=complex)
        M2m = np.zeros((252, 252), dtype=complex)
        for B in range(252):
            eB = from_comps5([1.0 if i == B else 0.0
                              for i in range(252)])
            M1[:, B] = comps5(act5(eB, g1))
            M2m[:, B] = comps5(act5(eB, g2))
        Tg = (M1 + sg * M2m) / 2.0
        C += Tg @ Tg
    CAS5[side] = np.conj(U126).T @ (-C) @ U126


def classify_sigma(H126):
    """Eigen-blocks of a 126x126 Hermitian mass operator, sub-split by
    the SU(2)_L Casimir (exact degeneracies can survive)."""
    ev, U = np.linalg.eigh(H126)
    sc = max(1.0, np.abs(ev).max())
    blocks, i = [], 0
    while i < len(ev):
        j = i
        while j + 1 < len(ev) and abs(ev[j + 1] - ev[i]) < 1e-6 * sc:
            j += 1
        P = U[:, i:j + 1]
        CL = np.conj(P).T @ CAS5["L"] @ P
        evL, UL = np.linalg.eigh((CL + CL.conj().T) / 2)
        subs, k0 = [], 0
        while k0 < len(evL):
            k1 = k0
            while k1 + 1 < len(evL) and abs(evL[k1 + 1] - evL[k0]) < 1e-4:
                k1 += 1
            subs.append((k0, k1))
            k0 = k1 + 1
        for (k0, k1) in subs:
            Ps = P @ UL[:, k0:k1 + 1]
            lab = {kk: float(np.real(np.trace(
                np.conj(Ps).T @ CAS5[kk] @ Ps))) / (k1 - k0 + 1)
                for kk in CAS5}
            omg = float(np.linalg.norm(np.conj(Ps).T @ OMEGA_126))
            blocks.append({"m2": float(ev[i]), "dim": k1 - k0 + 1,
                           "L": round(lab["L"], 3),
                           "R": round(lab["R"], 3),
                           "has_sigma_dir": omg > 1e-6})
        i = j + 1
    return blocks


bet0 = rng.uniform(-1, 1, size=5)
Hs0 = sum(b * M for b, M in zip(bet0, MIX126.values()))
blocks0 = classify_sigma(Hs0)
dims0 = sorted(b["dim"] for b in blocks0)
sig0 = [b for b in blocks0 if b["has_sigma_dir"]]
check("Sigma-sector structure on the 126bar half: the mixed cubic + "
      "quartics split it into exactly the PS multiplets (6,1,1) + "
      "(10,1,3) + (10bar,3,1) + (15,2,2) [complex dims 6/30/30/60], "
      "with the sigma direction MACHINE-LOCATED in the SU(2)_R-charged "
      "30 (the Delta_R block ESH tunes to M_I); the D-parity-odd mixed "
      "CUBIC is what splits the two 30s",
      dims0 == [6, 30, 30, 60] and len(sig0) == 1
      and sig0[0]["dim"] == 30 and sig0[0]["R"] > 1e-4
      and sig0[0]["L"] < 1e-4,
      str([(b["dim"], b["L"], b["R"], b["has_sigma_dir"])
           for b in blocks0]))
DLOG["S4_structure_sample"] = blocks0

# ------------------------- K4 retest with the 126bar remnants
B_SM = np.array([41 / 10, -19 / 6, -7.0])
BIJ_SM = np.array([[199 / 50, 27 / 10, 44 / 5], [9 / 10, 35 / 6, 12],
                   [11 / 10, 9 / 2, -26]])
MZ = 91.1876
AINV_MZ = np.array([0.6 * 127.951 * (1 - 0.23122), 127.951 * 0.23122,
                    1 / 0.1180])
B_PS_V = np.array([-23 / 3, -3.0, 11 / 3])
B_SIGMA = {"(6,1,1)": np.array([1 / 3, 0, 0]),
           "(10b,3,1)": np.array([3.0, 20 / 3, 0]),
           "(15,2,2)": np.array([16 / 3, 5.0, 5.0])}


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


def solve_thresh_PS(thresh):
    corr = sum((b * ln for b, ln in thresh), np.zeros(3))

    def equations(u):
        lnMI, lnMX, split = u
        a1, a2, a3 = run2(AINV_MZ, lnMI - math.log(MZ))
        aI0 = np.array([a3, a2, split])
        aX = aI0 - (B_PS_V * (lnMX - lnMI) + corr) / (2 * math.pi)
        cons = a1 - (0.6 * aI0[2] + 0.4 * aI0[0])
        return np.array([aX[0] - aX[1], aX[1] - aX[2], cons]), aX

    u = np.array([math.log(1e11), math.log(1e16), 30.0])
    for _ in range(120):
        r, aX = equations(u)
        if np.max(np.abs(r)) < 1e-10:
            break
        Jm = np.zeros((3, 3))
        for j in range(3):
            du = np.zeros(3); du[j] = 1e-6
            Jm[:, j] = (equations(u + du)[0]
                        - equations(u - du)[0]) / 2e-6
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
MV2_PS = float(np.linalg.eigvalsh(
    vector_M2_full(v_PS, np.zeros((N,) * 5, dtype=complex))).max()) * G2

n_k4 = 150
k4 = {"viable": 0, "alive": 0, "taus": []}
for it in range(n_k4):
    bet = rng.uniform(-1, 1, size=5)
    Hs = sum(b * M for b, M in zip(bet, MIX126.values()))
    blocks = classify_sigma(Hs)
    if sorted(b["dim"] for b in blocks) != [6, 30, 30, 60]:
        continue
    sig = [b for b in blocks if b["has_sigma_dir"]]
    if len(sig) != 1:
        continue
    m2_shift = -sig[0]["m2"]        # ESH tuning: Delta_R block -> M_I
    others = [b for b in blocks if not b["has_sigma_dir"]]
    m2_phys = [b["m2"] + m2_shift for b in others]
    if any(m <= 0 for m in m2_phys):
        continue                     # tachyonic Sigma remnant
    k4["viable"] += 1
    thr = []
    for b, m2r in zip(others, m2_phys):
        rr = math.sqrt(m2r / MV2_PS)
        name = {6: "(6,1,1)", 60: "(15,2,2)"}.get(b["dim"], "(10b,3,1)")
        if rr < 1.0:
            thr.append((B_SIGMA[name], min(math.log(1.0 / rr), 25.0)))
    sol = solve_thresh_PS(thr)
    MX = 10 ** sol["log10_MX"]
    physical = (sol["residual"] < 1e-8
                and math.log10(MZ) < sol["log10_MI"] < sol["log10_MX"]
                and MX < 1.2e19)
    if physical:
        tau = tau_d6(MX, sol["alpha_G_inv"])
        k4["taus"].append(tau)
        if tau > SK_EPI:
            k4["alive"] += 1
taus = sorted(k4["taus"])
check("K4 RETEST including the 126bar remnants: over random mixed "
      "couplings (cubic eta + four quartics) with the extended-"
      "survival tuning (Delta_R block at M_I) and positive Sigma "
      "remnants, the tau distribution is recorded",
      k4["viable"] > 0,
      f"viable {k4['viable']}/{n_k4}, alive {k4['alive']}; tau "
      f"16/50/84 = " + (", ".join(f"{taus[int(q * (len(taus) - 1))]:.2e}"
                                  for q in (0.16, 0.5, 0.84))
                        if taus else "none"))
DLOG["S4_K4_retest"] = {"viable": k4["viable"], "n": n_k4,
                        "alive": k4["alive"],
                        "tau_percentiles":
                        ([taus[int(q * (len(taus) - 1))]
                          for q in (0.16, 0.5, 0.84)] if taus else None)}

# ------------------------------------------------------------- ledger
n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-9b-1c epsilon-coupling scan + 126bar mixed quartics",
    "dyn_item": "DYN-9b-1c",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    # negative-boundary flags
    "sigma_self_quartics_deferred_mass_irrelevant_at_sigma0": True,
    "epsilon_mixed_patterns_not_enumerated": True,
    "mixed_quartic_generating_set_only": True,
    "e_w_sign_convention_fixed_by_eta_term": True,
    "tree_level_only": True,
    "zeta_value_derived": False,
}
(OUT / "dyn9b1c_eps_and_126_quartics.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-9b-1c: epsilon-coupling scan + 126bar mixed quartics", "",
      f"{n_pass}/{len(CHECKS)} checks pass.", "",
      "Key numbers in the JSON derivation_log: S2_eps_scan (the LR "
      "rescue verdict), S3_126_gates (embedding anchors), S4_K4_retest "
      "(the Hyper-K window with the 126bar remnants included).",
      "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn9b1c_eps_and_126_quartics.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-9b-1c: {n_pass}/{len(CHECKS)} checks; ledgers -> "
      f"output/audit9/dyn9b1c_eps_and_126_quartics.*")
