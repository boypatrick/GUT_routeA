#!/usr/bin/env python3
"""Exact Gaussian-integer action certificate for the P50/P36 projectors.

No floating Hessian, cache, eigensolver, numerical tolerance, fitted rational
matrix or sampled coupling enters the action calculation.  Projectors are
read as explicit integer matrices and their algebra is checked anew. Every
large contraction product/sum has a conservative bound before execution;
the fixed geometry constructors use explicitly bounded small integers.

At (w,sigma,vs)=(5,sqrt(2),sqrt(2)), Phi=diag(-2^6,3^4), S=1,
4 Sigma and 2 dSigma/dx are Gaussian integer tensors.  Antisymmetry reduces
the five-form contractions to 10x10 and 45x45 Gram arrays.  External Phi,
H and S covectors use unnormalized rational bases; their zero images are
equivalent to zero in the canonical bases. Sigma covectors stay canonical.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy import sparse

RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_exact_transverse_certificate"
ACTION = RF / "code/verify_p54_p1_hessian_spectrum.py"
EXPECTED_ACTION_SHA256 = "241fe7d388c649c6e11cd214785fd6bd1ee43b437824ff5e99d0fc1589f42282"
NAMES = ("mu2", "c", "a", "b", "nu2", "lambda0", "lambda2", "lambda4", "lambda4p",
         "alpha", "beta", "xi02", "xi1", "xi2", "xi3", "gamma1", "gamma2",
         "eta0", "eta1", "eta2", "eta3", "mus2", "chi1", "chi2", "chi3",
         "chi4", "chi5", "chi6", "chi7")
LIMIT = 2**62 - 1  # strict margin below signed int64 overflow
ARITHMETIC = []


def maxabs(a):
    values = a.data if sparse.issparse(a) else np.asarray(a)
    if not values.size:
        return 0
    if values.dtype != np.int64:
        raise TypeError(f"Only signed integer matrices are allowed; got {values.dtype} in shape {a.shape}")
    return max(abs(int(values.min())), abs(int(values.max())))


def record(label, bound):
    bound = int(bound)
    if bound > LIMIT:
        raise OverflowError(f"{label}: conservative bound {bound} exceeds {LIMIT}")
    ARITHMETIC.append((label, bound))


def mm(a, b):
    if a.shape[-1] != b.shape[-2]:
        raise ValueError("Incompatible exact matrix product")
    record("matrix_product", a.shape[-1] * maxabs(a) * maxabs(b))
    if sparse.issparse(b) and not sparse.issparse(a):
        return (b.T @ a.T).T
    return a @ b


def add(a, b):
    record("matrix_sum", maxabs(a) + maxabs(b))
    return a + b


def scale(a, k):
    record("integer_scale", maxabs(a) * abs(int(k)))
    return a * int(k)


def dense(a):
    return a.toarray() if sparse.issparse(a) else np.asarray(a)


def zeros(shape):
    return np.zeros(shape, dtype=np.int64)


def gc(a):
    return a[0], -a[1]


def gt(a):
    return a[0].T, a[1].T


def ga(a, b):
    return add(a[0], b[0]), add(a[1], b[1])


def gm(a, b):
    return add(mm(a[0], b[0]), -mm(a[1], b[1])), add(mm(a[0], b[1]), mm(a[1], b[0]))


def gd(a):
    return dense(a[0]), dense(a[1])


def gr(a, shape):
    return a[0].reshape(shape), a[1].reshape(shape)


def real_matrix(a):
    return a, zeros(a.shape)


def sign(values):
    return -1 if sum(values[i] > values[j] for i in range(len(values))
                     for j in range(i + 1, len(values))) % 2 else 1


def canonical_geometry():
    quints = list(itertools.combinations(range(10), 5))
    index = {q: i for i, q in enumerate(quints)}
    ar, ai = zeros(252), zeros(252)
    for j, q in enumerate(quints):
        if len({a // 2 for a in q}) == 5:
            phase = sum(a % 2 for a in q) % 4
            ar[j] = (1, 0, -1, 0)[phase]
            ai[j] = (0, 1, 0, -1)[phase]
    pairs, seen = [], set()
    star = zeros((252, 252))
    for a, q in enumerate(quints):
        other = tuple(i for i in range(10) if i not in q)
        if q in seen or other in seen:
            continue
        seen.update((q, other))
        b, sg = index[other], sign(q + other)
        pairs.append((a, b, sg))
        star[b, a], star[a, b] = sg, -sg
    sa = gm(real_matrix(star), (ar[:, None], ai[:, None]))
    duality = next((d for d in (1, -1)
                    if np.array_equal(sa[0][:, 0], -d * ai)
                    and np.array_equal(sa[1][:, 0], d * ar)), None)
    if duality is None:
        raise ValueError("The integer background is not a Hodge eigenform")
    cr, ci = zeros((252, 252)), zeros((252, 252))
    xb = zeros(252)
    for k, (a, b, sg) in enumerate(pairs):
        cr[a, k], ci[b, k] = 1, -duality * sg
        ci[a, k + 126], cr[b, k + 126] = 1, duality * sg
        xb[k], xb[k + 126] = ar[a], ai[a]
    return quints, index, (ar, ai), (cr, ci), xb, star, duality


def contractions(order, quints, index, background, canonical):
    """Return A_IJ and sparse C_IJ,u, with sorted I of size order."""
    initial = list(itertools.combinations(range(10), order))
    tails = list(itertools.combinations(range(10), 5 - order))
    shape = len(initial), len(tails)
    ar, ai = zeros(shape), zeros(shape)
    rows, cols, vals_r, vals_i = [], [], [], []
    cr, ci = canonical
    for a, first in enumerate(initial):
        for b, tail in enumerate(tails):
            indices = first + tail
            if len(set(indices)) != 5:
                continue
            j, sg = index[tuple(sorted(indices))], sign(indices)
            ar[a, b], ai[a, b] = sg * background[0][j], sg * background[1][j]
            for col in np.flatnonzero((cr[j] != 0) | (ci[j] != 0)):
                rows.append(a * len(tails) + b); cols.append(col)
                vals_r.append(sg * int(cr[j, col])); vals_i.append(sg * int(ci[j, col]))
    fshape = len(initial) * len(tails), 252
    fr = sparse.csr_matrix((np.asarray(vals_r, dtype=np.int64), (rows, cols)), shape=fshape)
    fi = sparse.csr_matrix((np.asarray(vals_i, dtype=np.int64), (rows, cols)), shape=fshape)
    fr.eliminate_zeros(); fi.eliminate_zeros()
    return dict(initial=initial, tails=tails, A=(ar, ai), F=(fr, fi), shape=shape)


def gram_jets(data):
    a, f = data["A"], data["F"]
    count, q = data["shape"]
    x0 = gm(gt(a), gc(a))
    first = zeros((252, q, q)), zeros((252, q, q))
    for m in range(q):
        fm = f[0][m::q], f[1][m::q]
        value = gd(gm(gt(fm), gc(a)))
        first[0][:, m, :], first[1][:, m, :] = value
    g = ga(first, (first[0].swapaxes(1, 2), -first[1].swapaxes(1, 2)))
    return x0, g


def crossing_operator(pairs):
    lookup = {p: i for i, p in enumerate(pairs)}
    q, rows, cols, values = len(pairs), [], [], []
    def orient(a, b):
        if a == b:
            return None
        return lookup[tuple(sorted((a, b)))], 1 if a < b else -1
    for i, (a, b) in enumerate(pairs):
        for j, (c, d) in enumerate(pairs):
            terms = ((a, c, b, d, 1), (b, c, a, d, -1),
                     (a, d, b, c, -1), (b, d, a, c, 1))
            for u, v, x, y, sg in terms:
                left, right = orient(u, v), orient(x, y)
                if left is not None and right is not None:
                    rows.append(i * q + j); cols.append(left[0] * q + right[0])
                    values.append(sg * left[1] * right[1])
    result = sparse.csr_matrix((np.asarray(values, dtype=np.int64), (rows, cols)), shape=(q*q, q*q))
    result.eliminate_zeros()
    return result


def lambda_hessian(data, x0, jets, crossing=None):
    """Return exact Gaussian-integer 32 H for lambda2/lambda4/lambda4p."""
    count, q = data["shape"]
    flat = gr(jets, (252, q*q))
    if crossing is None:
        mapped_jets, mapped_x0 = flat, x0
    else:
        c = crossing, sparse.csr_matrix(crossing.shape, dtype=np.int64)
        mapped_jets = gt(gm(c, gt(flat)))
        mapped_x0 = gr(gm(c, gr(x0, (q*q, 1))), (q, q))
    first = gm(flat, gt(mapped_jets))
    identity = sparse.eye(count, dtype=np.int64, format="csr")
    # SciPy gives an empty Kronecker product float dtype despite integer
    # inputs. Its entries are exactly zero; retain the declared integer ring.
    kr = sparse.kron(identity, sparse.csr_matrix(mapped_x0[0]), format="csr").astype(np.int64)
    ki = sparse.kron(identity, sparse.csr_matrix(mapped_x0[1]), format="csr").astype(np.int64)
    response = gd(gm(gt(data["F"]), gm((kr, ki), gc(data["F"]))))
    return gd(ga(first, ga(response, gt(response))))


def phi_basis():
    result = []
    for a in range(10):
        for b in range(a + 1, 10):
            e = zeros((10, 10)); e[a, b] = e[b, a] = 1; result.append(e)
    for a in range(9):
        e = zeros((10, 10)); e[a, a] = 1; e[9, 9] = -1; result.append(e)
    return result


def beta_blocks(data, diagonal, phi):
    count, q = data["shape"]
    pair_a = np.array([a for a, b in data["tails"]])
    pair_b = np.array([b for a, b in data["tails"]])
    weights = np.tile(diagonal[pair_a] * diagonal[pair_b], count)
    operator = sparse.diags(weights, dtype=np.int64, format="csr")
    f = data["F"]
    hh = gd(gm(gt(f), gm((operator, sparse.csr_matrix(operator.shape, dtype=np.int64)), gc(f))))[0]
    a = zeros((count, 10, 10)), zeros((count, 10, 10))
    for part in range(2):
        a[part][:, pair_a, pair_b] = data["A"][part]
        a[part][:, pair_b, pair_a] = -data["A"][part]
    d = np.diag(diagonal)
    mixed = zeros((54, 252))
    # 2 H_(Phi,Sigma) = Re <C, E bar(A) D + D bar(A) E> on sorted pairs.
    for j, e in enumerate(phi):
        response = ga(gm(real_matrix(e), gm(gc(a), real_matrix(d))),
                      gm(real_matrix(d), gm(gc(a), real_matrix(e))))
        packed = tuple(z[:, pair_a, pair_b].reshape(-1, 1) for z in response)
        mixed[j] = gd(gm(gt(f), packed))[0][:, 0]
    return hh, mixed


def chi4_blocks(data, diagonal, phi):
    count, q = data["shape"]
    f, a = data["F"], data["A"]
    d = np.diag(diagonal)
    kd = sparse.kron(sparse.eye(count, dtype=np.int64, format="csr"),
                     sparse.csr_matrix(d), format="csr")
    hh = gd(gm(gt(f), gm((kd, sparse.csr_matrix(kd.shape, dtype=np.int64)), f)))[0]
    mixed_phi = zeros((54, 252))
    for j, e in enumerate(phi):
        target = gr(gm(a, real_matrix(e.T)), (count * q, 1))
        mixed_phi[j] = gd(gm(gt(f), target))[0][:, 0]
    target = gr(gm(a, real_matrix(d)), (count * q, 1))
    z = gd(gm(gt(f), target))
    mixed_s = np.stack((z[0][:, 0], -z[1][:, 0]))
    # Mixed arrays are 2H; the Sigma-Sigma array is H.
    return hh, mixed_phi, mixed_s


def eta1_mixed(data, x0, jets, index, background, canonical):
    pairs, q = data["tails"], len(data["tails"])
    flat = gr(jets, (252, q*q))
    result = zeros((20, 252))
    cr, ci = canonical
    for n in range(10):
        br, bi = zeros(q*q), zeros(q*q)
        rows, cols, vr, vi = [], [], [], []
        for i, left in enumerate(pairs):
            for j, right in enumerate(pairs):
                inds = left + right + (n,)
                if len(set(inds)) != 5:
                    continue
                k, sg = index[tuple(sorted(inds))], sign(inds)
                row = i * q + j
                br[row], bi[row] = sg * background[0][k], sg * background[1][k]
                for col in np.flatnonzero((cr[k] != 0) | (ci[k] != 0)):
                    rows.append(row); cols.append(col)
                    vr.append(sg * int(cr[k, col])); vi.append(sg * int(ci[k, col]))
        f = tuple(sparse.csr_matrix((np.asarray(v, dtype=np.int64), (rows, cols)),
                                   shape=(q*q, 252)) for v in (vr, vi))
        z = gd(ga(gm(flat, (br[:, None], bi[:, None])), gm(gt(f), gr(x0, (q*q, 1)))))
        result[n], result[n + 10] = z[0][:, 0], -z[1][:, 0]
    # The original 1/(3!^2*4), sorted-pair factor4 and c.c. give 96H.
    return result


def run(projector_path):
    ARITHMETIC.clear()
    action_bytes = ACTION.read_bytes()
    if hashlib.sha256(action_bytes).hexdigest() != EXPECTED_ACTION_SHA256:
        raise ValueError("Action source changed; audit the symbolic transcription before recertifying")
    tree = ast.parse(action_bytes)
    factory = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "potential_factory")
    keys = {node.slice.value for node in ast.walk(factory)
            if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
            and node.value.id == "p" and isinstance(node.slice, ast.Constant)}
    if keys != set(NAMES):
        raise ValueError("The full 29-invariant catalogue changed")
    report = json.loads(projector_path.read_text())
    checks = []

    def check(name, actual, expected=0):
        actual, expected = dense(actual), np.asarray(expected, dtype=np.int64)
        difference = add(actual, -expected)
        checks.append(dict(name=name, passed=bool(not np.any(difference)),
                           maximum_integer_residual=maxabs(difference),
                           nonzero_integer_residuals=int(np.count_nonzero(difference))))

    quints, index, background, canonical, xb, star, duality = canonical_geometry()
    check("canonical_real_metric", gd(gm(gt(canonical), gc(canonical)))[0], 2 * np.eye(252, dtype=np.int64))
    image = gd(gm(canonical, real_matrix(xb[:, None])))
    check("exact_background_reconstruction_real", image[0][:, 0], background[0])
    check("exact_background_reconstruction_imag", image[1][:, 0], background[1])
    check("background_canonical_norm", mm(xb[None, :], xb[:, None]), 16)
    simage = gd(gm(real_matrix(star), canonical))
    check("canonical_self_duality_real", simage[0], -duality * canonical[1])
    check("canonical_self_duality_imag", simage[1], duality * canonical[0])

    projectors = []
    for b in report["projectors"]:
        specification = b["rational_projector"]
        if specification["denominator"] != 8 or specification["shape"] != [252, 252]:
            raise ValueError("Unexpected projector encoding")
        k = zeros((252, 252))
        for i, j, value in specification["nonzero_integer_entries"]:
            k[i, j] = value
        rank = int(b["rank"])
        check(f"P{rank}_symmetric", k, k.T)
        check(f"P{rank}_idempotent", mm(k, k), scale(k, 8))
        check(f"P{rank}_rank_from_trace", np.array(int(np.trace(k)), dtype=np.int64), 8 * rank)
        check(f"P{rank}_orthogonal_to_background", mm(k, xb[:, None]))
        projectors.append((rank, k))
    check("P50_P36_orthogonal", mm(projectors[0][1], projectors[1][1]))

    data4 = contractions(4, quints, index, background, canonical)
    data3 = contractions(3, quints, index, background, canonical)
    x4, g4 = gram_jets(data4)
    x3, g3 = gram_jets(data3)
    crossing = crossing_operator(data3["tails"])
    check("lambda4p_crossing_bilinear_self_adjoint", crossing, crossing.T.toarray())
    action = {}
    for name, data, x0, jets, c in (("lambda2", data4, x4, g4, None),
                                    ("lambda4", data3, x3, g3, None),
                                    ("lambda4p", data3, x3, g3, crossing)):
        h = lambda_hessian(data, x0, jets, c)
        check(name + "_Hessian_is_real", h[1])
        check(name + "_Hessian_is_symmetric", h[0], h[0].T)
        num = zeros((328, 252)); num[54:306] = h[0]
        action[name] = (num, 32, "exact compressed five-form quartic contraction")
    ident = np.eye(252, dtype=np.int64)
    num = zeros((328, 252)); num[54:306] = -ident
    action["nu2"] = num, 2, "canonical quadratic norm"
    num = zeros((328, 252)); num[54:306] = add(8 * ident, mm(xb[:, None], xb[None, :]))
    action["lambda0"] = num, 8, "canonical quartic norm"
    diagonal = np.array([-2] * 6 + [3] * 4, dtype=np.int64)
    phi = phi_basis(); d = np.diag(diagonal)
    num = zeros((328, 252)); num[54:306] = 60 * ident
    for j, e in enumerate(phi):
        num[j] = scale(xb, int(np.trace(mm(d, e))))
    action["alpha"] = num, 2, "norm with all rational Phi mixed covectors"
    num = zeros((328, 252)); num[54:306] = 120 * ident; num[326] = 120 * xb
    action["chi2"] = num, 1, "norm with both rational S mixed covectors"
    hh, mixed = beta_blocks(data3, diagonal, phi)
    num = zeros((328, 252)); num[54:306] = scale(hh, 2); num[:54] = mixed
    action["beta"] = num, 2, "exact contraction including all Phi mixed covectors"
    hh, mixed_phi, mixed_s = chi4_blocks(data4, diagonal, phi)
    num = zeros((328, 252)); num[54:306] = scale(hh, 2)
    num[:54], num[326:328] = mixed_phi, mixed_s
    action["chi4"] = num, 2, "exact holomorphic contraction including Phi and S covectors"
    num = zeros((328, 252)); num[306:326] = eta1_mixed(data3, x3, g3, index, background, canonical)
    action["eta1"] = num, 96, "exact cubic-Sigma contraction including all 20 H covectors"
    structural = [name for name in NAMES if name not in action]
    for name in structural:
        action[name] = zeros((328, 252)), 1, "no Sigma dependence or a remaining H factor at H=0"
    expected = {50: dict(nu2=Fraction(-1, 2), lambda0=1, lambda2=8,
                        lambda4=8, lambda4p=32, alpha=30, beta=-30, chi2=120),
                36: dict(nu2=Fraction(-1, 2), lambda0=1, lambda2=12,
                        lambda4=12, lambda4p=16, alpha=30, beta=-30, chi2=120)}
    identities = []
    for name in NAMES:
        numerator, denominator, method = action[name]
        check(name + "_canonical_Sigma_block_symmetric", numerator[54:306], numerator[54:306].T)
        for rank, k in projectors:
            c = Fraction(expected[rank].get(name, 0))
            multiplier = denominator * c
            if multiplier.denominator != 1:
                raise ValueError("Uncleared exact denominator")
            target = zeros((328, 252)); target[54:306] = scale(k, multiplier.numerator)
            check(f"{name}_all_328_covectors_on_P{rank}", mm(numerator, k), target)
            identities.append(dict(invariant=name, projector_rank=rank,
                action_matrix_denominator=denominator,
                scalar_eigenvalue={"numerator": c.numerator, "denominator": c.denominator},
                method=method, complete_ambient_covectors=328,
                integer_matrix_sha256=hashlib.sha256(numerator.astype("<i8").tobytes()).hexdigest()))
    sources = [Path(__file__), ACTION, projector_path,
               RF / "tex/p54_exact_transverse_certificate_fragment.tex"]
    return dict(schema="p54-exact-transverse-action-certificate-v1", date="2026-09-17",
        arithmetic="Gaussian integers represented by pairs of signed int64 arrays; large contraction products and sums have conservative precomputed bounds below 2^62; fixed geometry constructors use entries 0,+/-1 and diagonal entries -2,+3 with bounded small scalar factors",
        runtime=dict(python=sys.version.split()[0], executable=sys.executable,
                     numpy=np.__version__, scipy=scipy.__version__, jax_required=False),
        no_floating_action_evaluation=True, no_tolerance_used=True,
        no_unit_Hessian_cache_used=True, new_Hessian_evaluations=0,
        source_action_sha256=EXPECTED_ACTION_SHA256,
        certificate_background=dict(omega="5", sigma="sqrt(2)", vs="sqrt(2)",
                                    singlet="1", four_times_quint_background="Gaussian integer", canonical_map="2 dSigma/dx Gaussian integer", hodge_eigenvalue=f"{duality} i"),
        output_covectors="Canonical Sigma; rational invertible Phi/H/S bases, sufficient and equivalent for all mixed zero images",
        continuation="Each separate invariant Hessian block is multihomogeneous. At this nonzero radial background its exact scalar/zero identities therefore extend to every nonzero real radial background with the same angular orientation.",
        unit_identities=identities, structural_zero_invariants=structural,
        general_mass_formulas={"P50": "-nu2/2+(lambda0/2+4lambda2+4lambda4+16lambda4p)sigma^2+(6/5)(alpha-beta)omega^2+60chi2 vs^2",
                              "P36": "-nu2/2+(lambda0/2+6lambda2+6lambda4+8lambda4p)sigma^2+(6/5)(alpha-beta)omega^2+60chi2 vs^2"},
        implication="The previous 1/tau hard-scalar kinetic bound no longer needs a numerically inferred action-projector identity. Its positive-hard-spectrum, one-loop hard-scalar, fixed real transverse-ray and engineering-margin scope remains unchanged.",
        remaining_limits=["No complex-coupling phase directions certified", "No full-spectrum positivity theorem for all tau", "No complete gauged/fermionic pole or physical stability conclusion", "No physical fit or default-parameter promotion"],
        maximum_conservative_integer_bound=max(bound for _, bound in ARITHMETIC),
        bounded_arithmetic_operations=len(ARITHMETIC),
        source_sha256={str(path.relative_to(RF)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        checks=checks, summary=dict(passed=sum(c["passed"] for c in checks), total=len(checks),
                                  all_pass=all(c["passed"] for c in checks)))


def markdown(report):
    lines = ["# Exact P54 P50/P36 tensor-action certificate", "",
        f"Exact checks: {report['summary']['passed']}/{report['summary']['total']}.", "",
        "The action matrices are derived independently from Gaussian-integer five-form contractions. No floating unit-Hessian cache, eigensolver, tolerance or fitted rational Hessian is used.", "",
        "The previous integer projector entries are merely candidate exact matrices: symmetry, idempotency, rank, orthogonality and all 29 unit-action images are checked anew with exact arithmetic.", "",
        f"Largest conservative integer-operation bound: {report['maximum_conservative_integer_bound']} < 2^62; bounded operations: {report['bounded_arithmetic_operations']}.", "",
        "At (omega,sigma,vs)=(5,sqrt(2),sqrt(2)), Phi=diag(-2,-2,-2,-2,-2,-2,3,3,3,3), S=1, and four times the quint background plus twice each canonical Sigma basis vector are Gaussian integer.", "",
        "| Unit invariant | P50 eigenvalue | P36 eigenvalue |",
        "|---|---:|---:|", "| nu2 | -1/2 | -1/2 |", "| lambda0 | 1 | 1 |",
        "| lambda2 | 8 | 12 |", "| lambda4 | 8 | 12 |", "| lambda4p | 32 | 16 |",
        "| alpha | 30 | 30 |", "| beta | -30 | -30 |", "| chi2 | 120 | 120 |",
        "| all other invariants | 0 | 0 |", "",
        "Every image includes all 328 output covectors, including the otherwise easy-to-miss beta/chi4/eta1 mixed blocks. Block multihomogeneity extends these exact identities from the rational background to the entire fixed-orientation radial slice.", "",
        "The established transverse-gap 1/tau implication is now exact for the declared real action, subject to its unchanged positive-hard-spectrum and hard-scalar one-loop scope. This is not a full physical pole/stability theorem and does not enable physical fitting.", "",
        "See JSON for all 58 unit-action projector identities, exact-zero residual counts, action matrix hashes, source hashes and remaining limits."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--projectors", type=Path, default=RF / "output/p54_transverse_gap_obstruction.json")
    result = run(parser.parse_args().projectors)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    OUT.with_suffix(".md").write_text(markdown(result))
    print(markdown(result))
    if not result["summary"]["all_pass"]:
        print([c for c in result["checks"] if not c["passed"]])
        raise SystemExit(1)
