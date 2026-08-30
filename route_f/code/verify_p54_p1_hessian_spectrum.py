#!/usr/bin/env python3
"""P1 stationary-point, Goldstone, full-Hessian, and spectrum verifier.

The field vector contains all 328 canonically normalized real coordinates:
54_R + 2*126_C + 2*10_C + 2*1_C.  The scalar potential is P54PQ-v2,
including chi7.  JAX differentiates and compiles the tensor action itself; no mass
formula is transcribed from the literature.

Install the small differentiator dependency with
    python3 -m pip install -r route_f/requirements-p54-p1.txt
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Callable

try:
    import jax
    import jax.numpy as anp
    jax.config.update("jax_enable_x64", True)
    grad = jax.grad
    hessian = jax.hessian
except ImportError as exc:  # pragma: no cover - explicit reproducibility gate
    raise SystemExit(
        "JAX is required; install route_f/requirements-p54-p1.txt"
    ) from exc

import numpy as np
from scipy.optimize import root

warnings.filterwarnings(
    "ignore", message="Casting complex values to real discards the imaginary part"
)


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "p54_p1_stationary_hessian_spectrum.tex"
SCRIPT = Path(__file__).resolve()

N_PHI = 54
N_SIGMA = 126
N_H = 10
N_REAL = 328
SL_PHI = slice(0, 54)
SL_SIGMA_RE = slice(54, 180)
SL_SIGMA_IM = slice(180, 306)
SL_H_RE = slice(306, 316)
SL_H_IM = slice(316, 326)
SL_S = slice(326, 328)

CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append({"group": group, "name": name, "pass": bool(condition), "detail": detail})


def permutation_sign(values: tuple[int, ...]) -> int:
    inversions = sum(values[i] > values[j] for i in range(len(values)) for j in range(i + 1, len(values)))
    return -1 if inversions % 2 else 1


def basis54() -> np.ndarray:
    basis = []
    for i in range(10):
        for j in range(i + 1, 10):
            matrix = np.zeros((10, 10))
            matrix[i, j] = matrix[j, i] = 1 / math.sqrt(2)
            basis.append(matrix)
    for k in range(1, 10):
        matrix = np.zeros((10, 10))
        matrix[:k, :k] += np.eye(k) / math.sqrt(k * (k + 1))
        matrix[k, k] = -k / math.sqrt(k * (k + 1))
        basis.append(matrix)
    result = np.asarray(basis)
    gram = np.einsum("aij,bij->ab", result, result)
    assert result.shape == (54, 10, 10)
    assert np.max(np.abs(gram - np.eye(54))) < 1e-13
    return result


def five_form_geometry() -> dict[str, np.ndarray | list[tuple[int, ...]]]:
    n = 10
    quints = list(itertools.combinations(range(n), 5))
    qindex = {q: i for i, q in enumerate(quints)}
    all_indices = set(range(n))
    pairs = []
    d5 = np.zeros((252, 252))
    seen: set[tuple[int, ...]] = set()
    for a, q in enumerate(quints):
        complement = tuple(sorted(all_indices - set(q)))
        if q in seen or complement in seen:
            continue
        seen.add(q)
        seen.add(complement)
        b = qindex[complement]
        sign = permutation_sign(tuple(list(q) + list(complement)))
        pairs.append((a, b, sign))
        d5[b, a] = sign
        d5[a, b] = -sign

    uhol = [(np.eye(n)[2 * k] + 1j * np.eye(n)[2 * k + 1]) / math.sqrt(2) for k in range(5)]
    omega_full = np.zeros((n,) * 5, dtype=complex)
    base = np.einsum("i,j,k,l,m->ijklm", *uhol)
    for pm in itertools.permutations(range(5)):
        omega_full += permutation_sign(pm) * np.transpose(base, pm)
    omega = np.asarray([omega_full[q] for q in quints])
    omega /= np.linalg.norm(omega)
    eigen = complex((d5 @ omega) @ np.conj(omega))
    duality = 1.0 if abs(eigen - 1j) < abs(eigen + 1j) else -1.0
    u126 = np.zeros((252, 126), dtype=complex)
    for k, (a, b, sign) in enumerate(pairs):
        u126[a, k] = 1 / math.sqrt(2)
        u126[b, k] = -duality * 1j * sign / math.sqrt(2)
    omega126 = np.conj(u126).T @ omega
    assert abs(np.linalg.norm(omega126) - 1) < 1e-12

    lookup = np.zeros(n**5, dtype=int)
    signs = np.zeros(n**5, dtype=float)
    for flat, indices in enumerate(itertools.product(range(n), repeat=5)):
        if len(set(indices)) < 5:
            continue
        ordered = tuple(sorted(indices))
        perm = tuple(ordered.index(value) for value in indices)
        lookup[flat] = qindex[ordered]
        signs[flat] = permutation_sign(perm)
    return {"quints": quints, "u126": u126, "omega": omega, "omega126": omega126, "lookup": lookup, "signs": signs}


B54 = basis54()
GEOM = five_form_geometry()
U126 = np.asarray(GEOM["u126"])
LOOKUP = np.asarray(GEOM["lookup"], dtype=int)
SIGNS = np.asarray(GEOM["signs"], dtype=float)
QUINTS = GEOM["quints"]


def unpack(x: anp.ndarray) -> tuple[anp.ndarray, ...]:
    phi54 = anp.einsum("a,aij->ij", x[SL_PHI], B54)
    z126 = (x[SL_SIGMA_RE] + 1j * x[SL_SIGMA_IM]) / anp.sqrt(2.0)
    z126bar = (x[SL_SIGMA_RE] - 1j * x[SL_SIGMA_IM]) / anp.sqrt(2.0)
    sigma_c = anp.dot(U126, z126)
    sigma_bar_c = anp.dot(anp.conj(U126), z126bar)
    sigma = anp.reshape(SIGNS * sigma_c[LOOKUP], (10, 10, 10, 10, 10))
    sigma_bar = anp.reshape(SIGNS * sigma_bar_c[LOOKUP], (10, 10, 10, 10, 10))
    h = (x[SL_H_RE] + 1j * x[SL_H_IM]) / anp.sqrt(2.0)
    hbar = (x[SL_H_RE] - 1j * x[SL_H_IM]) / anp.sqrt(2.0)
    singlet = (x[326] + 1j * x[327]) / anp.sqrt(2.0)
    singlet_bar = (x[326] - 1j * x[327]) / anp.sqrt(2.0)
    return phi54, sigma, sigma_bar, h, hbar, singlet, singlet_bar


def potential_factory(p: dict[str, float | complex]) -> Callable[[anp.ndarray], anp.ndarray]:
    e = lambda subs, *args: anp.einsum(subs, *args, optimize=True)
    f3, f4, f5 = math.factorial(3), math.factorial(4), math.factorial(5)

    def potential(x: anp.ndarray) -> anp.ndarray:
        a54, sig, sigb, hv, hvb, sv, svb = unpack(x)
        v = -p["mu2"] / 2 * e("ij,ij->", a54, a54)
        v += p["c"] / 3 * e("ij,jk,ki->", a54, a54, a54)
        aa = e("ij,ij->", a54, a54)
        v += p["a"] / 4 * aa * aa
        v += p["b"] / 2 * e("ij,jk,kl,li->", a54, a54, a54, a54)
        ss = e("ijklm,ijklm->", sig, sigb)
        v += -p["nu2"] / (2 * f5) * ss
        v += p["lambda0"] / (4 * f5**2) * ss * ss
        v += p["lambda2"] / f4**2 * e("ijklm,ijkln,opqrm,opqrn->", sig, sigb, sig, sigb)
        v += p["lambda4"] / (f3**2 * 2**2) * e("ijklm,ijkno,pqrlm,pqrno->", sig, sigb, sig, sigb)
        v += p["lambda4p"] / f3**2 * e("ijklm,ijkno,pqrln,pqrmo->", sig, sigb, sig, sigb)
        v += p["alpha"] / (2 * f5) * aa * ss
        v += p["beta"] / f3 * e("ij,kl,mnoik,mnojl->", a54, a54, sig, sigb)

        hh = e("i,i->", hv, hvb)
        v += -p["xi02"] * hh
        v += p["xi1"] * hh * hh
        v += p["xi2"] * e("i,i->", hv, hv) * e("j,j->", hvb, hvb)
        v += p["xi3"] * e("ij,i,j->", a54, hv, hvb)
        v += p["gamma1"] / f4 * e("ijklm,ijkln,m,n->", sig, sigb, hv, hvb)
        v += p["gamma2"] / f4 * e("ijklm,ijkln,n,m->", sig, sigb, hv, hvb)
        v += p["eta0"] / 2 * aa * hh
        eta1_expr = e("ijklm,ijkpq,lmpqn,n->", sig, sigb, sig, hv) / (f3**2 * 2**2)
        eta1_bar = e("ijklm,ijkpq,lmpqn,n->", sigb, sig, sigb, hvb) / (f3**2 * 2**2)
        v += p["eta1"] * eta1_expr + anp.conj(p["eta1"]) * eta1_bar
        v += p["eta2"] * e("ij,ik,j,k->", a54, a54, hv, hvb)
        eta3_expr = e("ijklm,ijkln,m,n->", sig, sig, hv, hv) / f4
        eta3_bar = e("ijklm,ijkln,m,n->", sigb, sigb, hvb, hvb) / f4
        v += p["eta3"] * eta3_expr + anp.conj(p["eta3"]) * eta3_bar

        v += -p["mus2"] * sv * svb + p["chi1"] * (sv * svb) ** 2
        v += p["chi2"] * ss * sv * svb + p["chi3"] * aa * sv * svb
        chi4_expr = e("ijklm,ijkln,mn->", sig, sig, a54) * sv / f4
        chi4_bar = e("ijklm,ijkln,mn->", sigb, sigb, a54) * svb / f4
        v += p["chi4"] * chi4_expr + anp.conj(p["chi4"]) * chi4_bar
        v += p["chi5"] * hh * sv * svb
        v += p["chi6"] * e("i,i->", hv, hv) * svb + anp.conj(p["chi6"]) * e("i,i->", hvb, hvb) * sv
        chi7_expr = e("ij,i,j->", a54, hv, hv) * svb
        chi7_bar = e("ij,i,j->", a54, hvb, hvb) * sv
        v += p["chi7"] * chi7_expr + anp.conj(p["chi7"]) * chi7_bar
        return anp.real(v)

    return potential


def chi7_only_factory(chi7: complex) -> Callable[[anp.ndarray], anp.ndarray]:
    """The new invariant alone, used for a cheap nonzero-Hessian audit."""
    def term(x: anp.ndarray) -> anp.ndarray:
        a54 = anp.einsum("a,aij->ij", x[SL_PHI], B54)
        hv = (x[SL_H_RE] + 1j * x[SL_H_IM]) / anp.sqrt(2.0)
        hvb = (x[SL_H_RE] - 1j * x[SL_H_IM]) / anp.sqrt(2.0)
        sv = (x[326] + 1j * x[327]) / anp.sqrt(2.0)
        svb = (x[326] - 1j * x[327]) / anp.sqrt(2.0)
        expr = anp.einsum("ij,i,j->", a54, hv, hv) * svb
        exprb = anp.einsum("ij,i,j->", a54, hvb, hvb) * sv
        return anp.real(chi7 * expr + anp.conj(chi7) * exprb)
    return term


def vacuum_vector(omega: float, sigma: float, vs: float) -> np.ndarray:
    x = np.zeros(N_REAL)
    matrix = omega * np.diag([-2 / 5] * 6 + [3 / 5] * 4)
    x[SL_PHI] = np.einsum("aij,ij->a", B54, matrix)
    coeff = sigma * np.asarray(GEOM["omega126"])
    x[SL_SIGMA_RE] = math.sqrt(2) * coeff.real
    x[SL_SIGMA_IM] = math.sqrt(2) * coeff.imag
    x[326] = vs
    return x


def benchmark() -> tuple[dict[str, float | complex], tuple[float, float, float]]:
    # Dimensionless diagnostic point.  It is deliberately not presented as a
    # phenomenological fit; the purpose is to exercise every Hessian block.
    omega, sigma, vs = 1.0, 0.35, 0.25
    p: dict[str, float | complex] = {
        "mu2": 0.0, "nu2": 0.0, "xi02": 0.14616372381894022, "mus2": 0.0,
        "c": -0.10, "xi3": -0.20,
        "a": 0.50, "b": 0.40,
        "lambda0": 0.90, "lambda2": 0.80, "lambda4": 0.60, "lambda4p": 0.40,
        "alpha": 0.25, "beta": 0.10,
        "xi1": 0.35, "xi2": 0.08,
        "gamma1": 0.12, "gamma2": -0.07,
        "eta0": 0.20, "eta1": 0.025, "eta2": 0.10, "eta3": 0.018,
        "chi1": 0.80, "chi2": 0.0010, "chi3": 0.05,
        "chi4": 0.030, "chi5": 0.20, "chi6": 0.025, "chi7": 0.030,
    }
    # These are derived from the tensor action in this card's normalization.
    radial_a = 36 / 25 * p["a"] + 42 / 125 * p["b"]
    d = p["alpha"] - p["beta"]
    p["mu2"] = p["c"] * omega / 5 + (5 / 3) * radial_a * omega**2 + d * sigma**2 + p["chi3"] * vs**2
    p["nu2"] = p["lambda0"] * sigma**2 + 12 / 5 * d * omega**2 + 120 * p["chi2"] * vs**2
    p["mus2"] = p["chi1"] * vs**2 + 120 * p["chi2"] * sigma**2 + 12 / 5 * p["chi3"] * omega**2
    return p, (omega, sigma, vs)


def radial_potential(r: np.ndarray, p: dict[str, float | complex]) -> float:
    w, s, v = r
    aq = 36 / 25 * p["a"] + 42 / 125 * p["b"]
    d = p["alpha"] - p["beta"]
    return float(
        -6 / 5 * p["mu2"] * w**2 + 4 / 25 * p["c"] * w**3 + aq * w**4
        -0.5 * p["nu2"] * s**2 + 0.25 * p["lambda0"] * s**4
        + 6 / 5 * d * w**2 * s**2
        -0.5 * p["mus2"] * v**2 + 0.25 * p["chi1"] * v**4
        + 60 * p["chi2"] * s**2 * v**2 + 6 / 5 * p["chi3"] * w**2 * v**2
    )


def radial_gradient(r: np.ndarray, p: dict[str, float | complex]) -> np.ndarray:
    w, s, v = r
    aq = 36 / 25 * p["a"] + 42 / 125 * p["b"]
    d = p["alpha"] - p["beta"]
    return np.array([
        -12 / 5 * p["mu2"] * w + 12 / 25 * p["c"] * w**2 + 4 * aq * w**3 + 12 / 5 * d * w * s**2 + 12 / 5 * p["chi3"] * w * v**2,
        -p["nu2"] * s + p["lambda0"] * s**3 + 12 / 5 * d * w**2 * s + 120 * p["chi2"] * s * v**2,
        -p["mus2"] * v + p["chi1"] * v**3 + 120 * p["chi2"] * s**2 * v + 12 / 5 * p["chi3"] * w**2 * v,
    ], dtype=float)


def stationary_points(p: dict[str, float | complex]) -> list[dict[str, Any]]:
    roots: list[np.ndarray] = []
    for start in itertools.product((-1.25, -0.4, 0.0, 0.4, 1.25), repeat=3):
        solution = root(lambda r: radial_gradient(r, p), np.asarray(start, dtype=float), method="lm")
        if np.linalg.norm(radial_gradient(solution.x, p)) > 2e-8:
            continue
        candidate = np.where(np.abs(solution.x) < 1e-8, 0.0, solution.x)
        if not any(np.linalg.norm(candidate - old) < 2e-5 for old in roots):
            roots.append(candidate)
    rows = []
    for r in roots:
        step = 2e-5
        h = np.zeros((3, 3))
        for i in range(3):
            ei = np.zeros(3); ei[i] = step
            for j in range(3):
                ej = np.zeros(3); ej[j] = step
                h[i, j] = (radial_potential(r + ei + ej, p) - radial_potential(r + ei - ej, p) - radial_potential(r - ei + ej, p) + radial_potential(r - ei - ej, p)) / (4 * step**2)
        rows.append({"omega": float(r[0]), "sigma": float(r[1]), "vs": float(r[2]), "V": radial_potential(r, p), "radial_hessian_eigenvalues": np.linalg.eigvalsh(h).tolist()})
    return sorted(rows, key=lambda row: row["V"])


def act5(tensor: np.ndarray, a: int, b: int) -> np.ndarray:
    delta = np.zeros_like(tensor)
    for slot in range(5):
        source = np.moveaxis(tensor, slot, 0)
        target = np.moveaxis(delta, slot, 0)
        target[a] += source[b]
        target[b] -= source[a]
    return delta


def act5_general(tensor: np.ndarray, generator: np.ndarray) -> np.ndarray:
    delta = np.zeros_like(tensor)
    for slot in range(5):
        source = np.moveaxis(tensor, slot, 0)
        transformed = np.tensordot(generator, source, axes=(1, 0))
        delta += np.moveaxis(transformed, 0, slot)
    return delta


def comp5(tensor: np.ndarray) -> np.ndarray:
    return np.asarray([tensor[q] for q in QUINTS])


def gauge_orbit(x0: np.ndarray) -> np.ndarray:
    a54, sigma, _, _, _, _, _ = unpack(anp.asarray(x0))
    a54 = np.asarray(a54, dtype=float)
    sigma = np.asarray(sigma, dtype=complex)
    columns = []
    for a in range(10):
        for b in range(a + 1, 10):
            generator = np.zeros((10, 10))
            generator[a, b] = 1 / math.sqrt(2)
            generator[b, a] = -1 / math.sqrt(2)
            d54 = generator @ a54 - a54 @ generator
            ds = comp5(act5(sigma, a, b)) / math.sqrt(2)
            dz = np.conj(U126).T @ ds
            column = np.zeros(N_REAL)
            column[SL_PHI] = np.einsum("aij,ij->a", B54, d54)
            column[SL_SIGMA_RE] = math.sqrt(2) * dz.real
            column[SL_SIGMA_IM] = math.sqrt(2) * dz.imag
            columns.append(column)
    return np.stack(columns, axis=1)


def pq_direction(x0: np.ndarray) -> np.ndarray:
    direction = np.zeros_like(x0)
    q_sigma = 2
    direction[SL_SIGMA_RE] = -q_sigma * x0[SL_SIGMA_IM]
    direction[SL_SIGMA_IM] = q_sigma * x0[SL_SIGMA_RE]
    q_h = -2
    direction[SL_H_RE] = -q_h * x0[SL_H_IM]
    direction[SL_H_IM] = q_h * x0[SL_H_RE]
    q_s = -4
    direction[326] = -q_s * x0[327]
    direction[327] = q_s * x0[326]
    return direction


def complex_generator_to_real(hermitian: np.ndarray) -> np.ndarray:
    """Real SO(10) generator for delta z=-i H z in five complex planes."""
    antihermitian = -1j * hermitian
    real = np.zeros((10, 10), dtype=float)
    for i in range(5):
        for j in range(5):
            a = antihermitian[i, j].real
            b = antihermitian[i, j].imag
            real[2 * i, 2 * j] = a
            real[2 * i, 2 * j + 1] = -b
            real[2 * i + 1, 2 * j] = b
            real[2 * i + 1, 2 * j + 1] = a
    assert np.max(np.abs(real + real.T)) < 1e-13
    return real


def su_generators(indices: tuple[int, ...]) -> list[np.ndarray]:
    """Hermitian fundamental generators with tr(T_a T_b)=delta_ab/2."""
    generators: list[np.ndarray] = []
    n = len(indices)
    for ii in range(n):
        for jj in range(ii + 1, n):
            i, j = indices[ii], indices[jj]
            symmetric = np.zeros((5, 5), dtype=complex)
            symmetric[i, j] = symmetric[j, i] = 0.5
            antisymmetric = np.zeros((5, 5), dtype=complex)
            antisymmetric[i, j] = -0.5j
            antisymmetric[j, i] = 0.5j
            generators.extend([symmetric, antisymmetric])
    for k in range(1, n):
        diagonal = np.zeros((5, 5), dtype=complex)
        normalization = math.sqrt(2 * k * (k + 1))
        for ii in range(k):
            diagonal[indices[ii], indices[ii]] = 1 / normalization
        diagonal[indices[k], indices[k]] = -k / normalization
        generators.append(diagonal)
    return [complex_generator_to_real(generator) for generator in generators]


def sm_generators() -> dict[str, list[np.ndarray]]:
    hypercharge = np.diag([-1 / 3, -1 / 3, -1 / 3, 1 / 2, 1 / 2]).astype(complex)
    return {
        "SU3": su_generators((0, 1, 2)),
        "SU2": su_generators((3, 4)),
        "Y": [complex_generator_to_real(hypercharge)],
    }


def representation_matrix(generator: np.ndarray) -> np.ndarray:
    """Generator on all 328 canonical real scalar coordinates."""
    result = np.zeros((N_REAL, N_REAL), dtype=float)
    for column, basis in enumerate(B54):
        variation = generator @ basis - basis @ generator
        result[SL_PHI, column] = np.einsum("aij,ij->a", B54, variation)

    qindex = {q: i for i, q in enumerate(QUINTS)}
    wedge = np.zeros((252, 252), dtype=float)
    for out_index, qout in enumerate(QUINTS):
        for slot in range(5):
            for source_index in range(10):
                coefficient = generator[qout[slot], source_index]
                if coefficient == 0:
                    continue
                ordered = list(qout)
                ordered[slot] = source_index
                if len(set(ordered)) < 5:
                    continue
                sorted_indices = tuple(sorted(ordered))
                perm = tuple(sorted_indices.index(value) for value in ordered)
                wedge[out_index, qindex[sorted_indices]] += coefficient * permutation_sign(perm)
    restricted = np.conj(U126).T @ wedge @ U126
    real, imag = restricted.real, restricted.imag
    result[SL_SIGMA_RE, SL_SIGMA_RE] = real
    result[SL_SIGMA_RE, SL_SIGMA_IM] = -imag
    result[SL_SIGMA_IM, SL_SIGMA_RE] = imag
    result[SL_SIGMA_IM, SL_SIGMA_IM] = real
    result[SL_H_RE, SL_H_RE] = generator
    result[SL_H_IM, SL_H_IM] = generator
    assert np.max(np.abs(result + result.T)) < 2e-12
    return result


def sm_representation_matrices() -> dict[str, list[np.ndarray]]:
    return {
        group: [representation_matrix(generator) for generator in generators]
        for group, generators in sm_generators().items()
    }


def su3_label(casimir: float) -> str:
    labels = {0.0: "1", 4 / 3: "3-pair", 3.0: "8", 10 / 3: "6-pair", 6.0: "10-pair", 16 / 3: "15-pair"}
    value, label = min(labels.items(), key=lambda row: abs(row[0] - casimir))
    return label if abs(value - casimir) < 2e-5 else f"C3={casimir:.6g}"


def classify_full_spectrum(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    absolute_tolerance: float,
    representations: dict[str, list[np.ndarray]],
) -> dict[str, Any]:
    rows = []
    maximum_leakage = 0.0
    start = 0
    while start < len(eigenvalues):
        stop = start + 1
        while stop < len(eigenvalues) and abs(eigenvalues[stop] - eigenvalues[start]) <= max(absolute_tolerance, 2e-7 * max(1.0, abs(eigenvalues[start]))):
            stop += 1
        block = eigenvectors[:, start:stop]
        restricted: dict[str, list[np.ndarray]] = {}
        casimirs: dict[str, np.ndarray] = {}
        for group, matrices in representations.items():
            restricted[group] = []
            for matrix in matrices:
                acted = matrix @ block
                small = block.T @ acted
                restricted[group].append(small)
                maximum_leakage = max(maximum_leakage, float(np.linalg.norm(acted - block @ small)))
            casimirs[group] = -sum(matrix @ matrix for matrix in restricted[group])
        combined = casimirs["SU3"] + math.sqrt(2) * casimirs["SU2"] + math.pi * casimirs["Y"]
        combined = (combined + combined.T) / 2
        values, vectors = np.linalg.eigh(combined)
        substart = 0
        while substart < len(values):
            substop = substart + 1
            while substop < len(values) and abs(values[substop] - values[substart]) < 2e-6:
                substop += 1
            selector = vectors[:, substart:substop]
            dimension = substop - substart
            c3 = float(np.trace(selector.T @ casimirs["SU3"] @ selector) / dimension)
            c2 = float(np.trace(selector.T @ casimirs["SU2"] @ selector) / dimension)
            y2 = max(0.0, float(np.trace(selector.T @ casimirs["Y"] @ selector) / dimension))
            j = max(0.0, (-1 + math.sqrt(max(1.0, 1 + 4 * c2))) / 2)
            su2_dimension = int(round(2 * j + 1))
            rows.append({
                "m2": float(np.mean(eigenvalues[start:stop])),
                "real_multiplicity": dimension,
                "SU3": su3_label(c3),
                "SU2_dimension": su2_dimension,
                "abs_hypercharge": math.sqrt(y2),
                "C3": c3,
                "C2": c2,
                "Y2": y2,
            })
            substart = substop
        start = stop
    return {
        "rows": rows,
        "real_dimension_sum": sum(row["real_multiplicity"] for row in rows),
        "maximum_invariance_leakage": maximum_leakage,
    }


def act_on_field_vector(vector: np.ndarray, generator: np.ndarray) -> np.ndarray:
    a54, sigma, _, hv, _, _, _ = unpack(anp.asarray(vector))
    a54 = np.asarray(a54, dtype=float)
    sigma = np.asarray(sigma, dtype=complex)
    hv = np.asarray(hv, dtype=complex)
    result = np.zeros(N_REAL)
    d54 = generator @ a54 - a54 @ generator
    result[SL_PHI] = np.einsum("aij,ij->a", B54, d54)
    ds = comp5(act5_general(sigma, generator))
    dz = np.conj(U126).T @ ds
    result[SL_SIGMA_RE] = math.sqrt(2) * dz.real
    result[SL_SIGMA_IM] = math.sqrt(2) * dz.imag
    dh = generator @ hv
    result[SL_H_RE] = math.sqrt(2) * dh.real
    result[SL_H_IM] = math.sqrt(2) * dh.imag
    return result


def classify_extra_zero_modes(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    zero_tolerance: float,
    symmetry_vectors: np.ndarray,
) -> dict[str, Any]:
    zero_vectors = eigenvectors[:, np.abs(eigenvalues) < zero_tolerance]
    symmetry_projector = symmetry_vectors @ np.linalg.pinv(symmetry_vectors, rcond=1e-11)
    residual = zero_vectors - symmetry_projector @ zero_vectors
    u, singular, _ = np.linalg.svd(residual, full_matrices=False)
    rank = int(np.sum(singular > 2e-8))
    extra = u[:, :rank]
    weights = {
        "Phi54": float(np.linalg.norm(extra[SL_PHI]) ** 2 / max(1, rank)),
        "Sigma126": float((np.linalg.norm(extra[SL_SIGMA_RE]) ** 2 + np.linalg.norm(extra[SL_SIGMA_IM]) ** 2) / max(1, rank)),
        "phi10": float((np.linalg.norm(extra[SL_H_RE]) ** 2 + np.linalg.norm(extra[SL_H_IM]) ** 2) / max(1, rank)),
        "S": float(np.linalg.norm(extra[SL_S]) ** 2 / max(1, rank)),
    }
    casimirs: dict[str, list[float]] = {}
    leakage: dict[str, float] = {}
    for group, generators in sm_generators().items():
        restricted = []
        maximum_leakage = 0.0
        for generator in generators:
            acted = np.column_stack([act_on_field_vector(extra[:, j], generator) for j in range(rank)])
            matrix = extra.T @ acted
            restricted.append(matrix)
            maximum_leakage = max(maximum_leakage, float(np.linalg.norm(acted - extra @ matrix)))
        casimir = -sum(matrix @ matrix for matrix in restricted)
        casimirs[group] = np.linalg.eigvalsh((casimir + casimir.T) / 2).tolist()
        leakage[group] = maximum_leakage
    return {"rank": rank, "field_weights": weights, "casimir_eigenvalues": casimirs, "invariance_leakage": leakage}


def group_eigenvalues(values: np.ndarray, abs_tol: float, rel_tol: float = 2e-7) -> list[dict[str, Any]]:
    values = np.sort(np.asarray(values, dtype=float))
    groups: list[list[float]] = []
    for value in values:
        if groups and abs(value - np.mean(groups[-1])) <= max(abs_tol, rel_tol * max(1.0, abs(value), abs(np.mean(groups[-1])))):
            groups[-1].append(float(value))
        else:
            groups.append([float(value)])
    return [{"m2": float(np.mean(group)), "multiplicity": len(group), "spread": float(max(group) - min(group))} for group in groups]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run() -> dict[str, Any]:
    started = time.time()
    p, vevs = benchmark()
    x0 = vacuum_vector(*vevs)
    potential = potential_factory(p)
    gradient = np.asarray(jax.jit(grad(potential))(anp.asarray(x0)), dtype=float)
    grad_scale = max(1.0, np.linalg.norm(x0), max(abs(float(v)) for v in p.values()))
    gradient_residual = float(np.linalg.norm(gradient) / grad_scale)
    check("stationarity", "the declared 328-field background is stationary", gradient_residual < 2e-8, f"scaled gradient={gradient_residual:.3e}")

    points = stationary_points(p)
    target = min(points, key=lambda row: (row["omega"] - vevs[0])**2 + (abs(row["sigma"]) - vevs[1])**2 + (abs(row["vs"]) - vevs[2])**2)
    check("stationarity", "radial solver recovers the declared nonzero branch", abs(target["omega"] - vevs[0]) < 2e-5 and abs(abs(target["sigma"]) - vevs[1]) < 2e-5 and abs(abs(target["vs"]) - vevs[2]) < 2e-5, f"root=({target['omega']:.6f},{target['sigma']:.6f},{target['vs']:.6f})")
    radial_global = bool(
        target["V"] <= points[0]["V"] + 2e-9
        and min(target["radial_hessian_eigenvalues"]) > 1e-7
    )
    check(
        "stationarity",
        "the declared orientation is the lowest enumerated radial branch",
        radial_global,
        f"V_target={target['V']:.12f}, V_min={points[0]['V']:.12f}",
    )

    h_start = time.time()
    hfull = np.asarray(jax.jit(hessian(potential))(anp.asarray(x0)), dtype=float)
    h_seconds = time.time() - h_start
    asymmetry = float(np.max(np.abs(hfull - hfull.T)))
    hfull = (hfull + hfull.T) / 2
    check("hessian", "the full 328x328 Hessian is symmetric", asymmetry < 2e-9, f"max asymmetry={asymmetry:.3e}")

    orbit = gauge_orbit(x0)
    rank_orbit = int(np.linalg.matrix_rank(orbit, tol=2e-10))
    check("goldstone", "the gauge orbit has exactly 33 broken directions", rank_orbit == 33, f"rank={rank_orbit}, unbroken={45-rank_orbit}")
    gauge_residual = float(np.linalg.norm(hfull @ orbit) / max(1.0, np.linalg.norm(hfull) * np.linalg.norm(orbit)))
    check("goldstone", "all broken-generator vectors are Hessian zero modes", gauge_residual < 3e-8, f"relative residual={gauge_residual:.3e}")

    qvec = pq_direction(x0)
    qphysical = qvec - orbit @ np.linalg.pinv(orbit, rcond=1e-11) @ qvec
    pq_residual = float(np.linalg.norm(hfull @ qphysical) / max(1.0, np.linalg.norm(hfull) * np.linalg.norm(qphysical)))
    goldstone_rank = int(np.linalg.matrix_rank(np.column_stack([orbit, qphysical]), tol=2e-10))
    check("goldstone", "one independent physical PQ Goldstone remains", goldstone_rank == 34 and np.linalg.norm(qphysical) > 1e-5, f"rank(gauge+PQ)={goldstone_rank}")
    check("goldstone", "the projected PQ direction is a Hessian zero mode", pq_residual < 3e-8, f"relative residual={pq_residual:.3e}")

    eig, eigvec = np.linalg.eigh(hfull)
    scale = max(1.0, float(np.max(np.abs(eig))))
    zero_tol = 3e-7 * scale
    nzero = int(np.sum(np.abs(eig) < zero_tol))
    nnegative = int(np.sum(eig < -zero_tol))
    check("spectrum", "the scalar eigenvalue census is 34 symmetry zeros plus one complex doublet", nzero == 38 and nnegative == 0, f"zero={nzero}, negative={nnegative}, tol={zero_tol:.3e}")
    spectrum_groups = group_eigenvalues(eig, abs_tol=zero_tol)
    symmetry_vectors = np.column_stack([orbit, qphysical])
    extra_zero = classify_extra_zero_modes(
        eig, eigvec, zero_tol, symmetry_vectors
    )
    doublet_casimirs = extra_zero["casimir_eigenvalues"]
    doublet_identified = bool(
        extra_zero["rank"] == 4
        and max(abs(value) for value in doublet_casimirs["SU3"]) < 3e-7
        and max(abs(value - 0.75) for value in doublet_casimirs["SU2"]) < 3e-7
        and max(abs(value - 0.25) for value in doublet_casimirs["Y"]) < 3e-7
        and max(extra_zero["invariance_leakage"].values()) < 3e-7
    )
    check(
        "doublet",
        "the four non-symmetry zero modes form one (1,2,+/-1/2) real multiplet",
        doublet_identified,
        f"rank={extra_zero['rank']}, Casimirs={doublet_casimirs}",
    )
    representations = sm_representation_matrices()
    irrep_spectrum = classify_full_spectrum(
        eig, eigvec, zero_tol, representations
    )
    irrep_labels_complete = bool(
        irrep_spectrum["real_dimension_sum"] == N_REAL
        and irrep_spectrum["maximum_invariance_leakage"] < 2e-7
        and all(not row["SU3"].startswith("C3=") for row in irrep_spectrum["rows"])
    )
    check(
        "spectrum",
        "the SM-Casimir irrep ledger covers all 328 real scalar coordinates",
        irrep_labels_complete,
        f"dimension={irrep_spectrum['real_dimension_sum']}, leakage={irrep_spectrum['maximum_invariance_leakage']:.3e}",
    )

    vector_m2 = np.linalg.eigvalsh(orbit.T @ orbit)
    vector_groups = group_eigenvalues(vector_m2, abs_tol=2e-9)
    check("vector", "vector spectrum has 12 massless and 33 massive generators", sum(row["multiplicity"] for row in vector_groups if abs(row["m2"]) < 2e-9) == 12 and sum(row["multiplicity"] for row in vector_groups if row["m2"] > 2e-9) == 33, str([(round(row['m2'], 8), row['multiplicity']) for row in vector_groups]))

    # Block-completeness and new-operator participation checks.
    field_counts = {"Phi54_real": 54, "Sigma126_real": 252, "phi10_real": 20, "S_real": 2}
    check("hessian", "field census is complete", sum(field_counts.values()) == N_REAL and hfull.shape == (N_REAL, N_REAL), json.dumps(field_counts, sort_keys=True))
    h_chi7 = np.asarray(
        jax.jit(hessian(chi7_only_factory(complex(p["chi7"]))))(anp.asarray(x0)),
        dtype=float,
    )
    chi7_delta = float(np.linalg.norm(h_chi7))
    check("hessian", "the newly required chi7 invariant enters the computed Hessian", chi7_delta > 1e-5, f"Frobenius delta={chi7_delta:.6e}")

    physical_positive = bool(nnegative == 0)
    p1_closed = bool(
        physical_positive and nzero == 38 and doublet_identified
        and irrep_labels_complete and radial_global
    )
    status = (
        "semidefinite_radial_global_vacuum_with_one_tuned_doublet"
        if p1_closed
        else ("candidate_local_minimum" if physical_positive else "stationary_saddle")
    )
    passed = sum(row["pass"] for row in CHECKS)
    return {
        "schema": "route-f-p54-p1-full-hessian-v1",
        "model_id": "P54PQ-v2",
        "date": "2026-08-30",
        "status": {
            "p1_stationary_background_built": gradient_residual < 2e-8,
            "full_328_real_hessian_built": True,
            "goldstone_alignment_closed": rank_orbit == 33 and goldstone_rank == 34,
            "benchmark_classification": status,
            "phenomenological_fit_claimed": False,
            "radial_global_within_enumerated_branches": radial_global,
            "one_light_doublet_tuned": doublet_identified,
            "p1_gate_fully_closed": p1_closed,
        },
        "field_census": field_counts,
        "vevs": {"omega": vevs[0], "sigma": vevs[1], "vs": vevs[2]},
        "parameters": {key: ([value.real, value.imag] if isinstance(value, complex) else value) for key, value in p.items()},
        "radial_stationary_points": points,
        "stationarity": {"full_gradient_norm": float(np.linalg.norm(gradient)), "scaled_gradient_residual": gradient_residual},
        "hessian": {"shape": list(hfull.shape), "symmetry_residual": asymmetry, "wall_seconds": h_seconds, "min_eigenvalue": float(eig[0]), "max_eigenvalue": float(eig[-1]), "zero_tolerance": zero_tol, "zero_count": nzero, "negative_count": nnegative, "eigenvalue_groups": spectrum_groups},
        "goldstone": {"gauge_orbit_rank": rank_orbit, "unbroken_generator_count": 45-rank_orbit, "gauge_hessian_residual": gauge_residual, "gauge_plus_physical_pq_rank": goldstone_rank, "physical_pq_norm": float(np.linalg.norm(qphysical)), "pq_hessian_residual": pq_residual},
        "light_doublet": extra_zero,
        "scalar_sm_irrep_spectrum": irrep_spectrum,
        "vector_spectrum_over_g2": vector_groups,
        "chi7_hessian_frobenius_delta": chi7_delta,
        "checks": CHECKS,
        "summary": {"passed": passed, "total": len(CHECKS), "all_pass": passed == len(CHECKS)},
        "runtime_seconds": time.time() - started,
        "sources": [{"path": str(TEX.relative_to(REPO)), "sha256": sha256(TEX)}, {"path": str(SCRIPT.relative_to(REPO)), "sha256": sha256(SCRIPT)}],
    }


def markdown(report: dict[str, Any]) -> str:
    h = report["hessian"]
    status = report["status"]
    lines = [
        "# P54PQ-v2 P1 stationary point, Hessian, and spectrum",
        "",
        f"Status: **{report['summary']['passed']}/{report['summary']['total']} checks passed**; benchmark is `{status['benchmark_classification']}`.",
        "",
        "The machine Hessian covers all `328 = 54 + 252 + 20 + 2` real scalar coordinates and includes the new `chi7` invariant.",
        "",
        f"- full-gradient residual: `{report['stationarity']['scaled_gradient_residual']:.3e}`",
        f"- gauge orbit: rank `{report['goldstone']['gauge_orbit_rank']}`; gauge plus physical PQ: rank `{report['goldstone']['gauge_plus_physical_pq_rank']}`",
        f"- scalar spectrum: `{h['zero_count']}` zero, `{h['negative_count']}` negative modes at tolerance `{h['zero_tolerance']:.3e}`",
        f"- Hessian wall time: `{h['wall_seconds']:.2f} s`",
        "",
        "## Scalar eigenvalue groups",
        "",
        "| m^2 | multiplicity | spread |",
        "|---:|---:|---:|",
    ]
    for row in h["eigenvalue_groups"]:
        lines.append(f"| {row['m2']:.9e} | {row['multiplicity']} | {row['spread']:.2e} |")
    lines.extend([
        "",
        "## Complete SM-Casimir spectrum",
        "",
        "`3-pair`, `6-pair`, and similar labels denote a representation plus its conjugate in the real Hessian.",
        "",
        "| m^2 | SU(3) | SU(2) dim | |Y| | real mult. |",
        "|---:|---|---:|---:|---:|",
    ])
    for row in report["scalar_sm_irrep_spectrum"]["rows"]:
        lines.append(f"| {row['m2']:.9e} | {row['SU3']} | {row['SU2_dimension']} | {row['abs_hypercharge']:.6g} | {row['real_multiplicity']} |")
    lines.extend(["", "## Checks", "", "| Group | Check | Result |", "|---|---|---|"])
    for row in report["checks"]:
        lines.append(f"| {row['group']} | {row['name']} | {'PASS' if row['pass'] else 'FAIL'} |")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report = run()
    (OUTPUT / "p54_p1_stationary_hessian_spectrum.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUTPUT / "p54_p1_stationary_hessian_spectrum.md").write_text(markdown(report), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    print(json.dumps(report["status"], sort_keys=True))
    if not report["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
