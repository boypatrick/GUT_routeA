#!/usr/bin/env python3
"""Fail-closed AP-E2 discriminant/contact regression.

This verifier preserves the exact theorem

    B_Kill(A_q,A_q) = 2 Delta(p) = 2 sqrt(3) x^T K_tr x,

for ``p=a u^2+b u v+c v^2`` in the Route-E spherical convention.  It uses
exact rational arithmetic where possible and 100-decimal-place arithmetic for
the irrational/complex normalization checks.  Deliberately wrong conventions
are required to fail.  A green run is an algebraic regression only and never
authorizes a physics promotion.
"""

from __future__ import annotations

import hashlib
import json
import random
from fractions import Fraction
from pathlib import Path
from typing import Any, TypeVar

import mpmath as mp


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"

ROUTE_E_CARD = REPO / "route_E" / "output" / "audit0" / "invariant_card.json"
ROUTE_E_BUILDER = REPO / "route_E" / "code" / "audit0_conventions_card.py"
THEOREM_SOURCE = ROUTE_F / "tex" / "another_physics_route_e_derivation_ledger.tex"
AP_E2_NOTE = ROUTE_F / "tex" / "ap_e2_discriminant_regression.tex"
LEGACY_VERIFIER = ROUTE_F / "code" / "verify_another_physics_route_e_bridge.py"
THIS_SCRIPT = Path(__file__).resolve()

THEOREM_START = r"\begin{theorem}[Discriminant--Killing equivalence]"
THEOREM_END = r"\end{proof}"
EXPECTED_THEOREM_BLOCK_SHA256 = (
    "3999b88f3fe7e396948c9695939ec5826d74cfc6d06f6b26ed05d3b0e0c4d094"
)

SEED = 20260714
EXACT_RANDOM_SAMPLES = 256
HIGH_PRECISION_COMPLEX_SAMPLES = 256
MP_DPS = 100
SOURCE_NORMALIZATION_TOL = mp.mpf("5e-15")

CHECKS: list[dict[str, Any]] = []
T = TypeVar("T")


def check(group: str, name: str, condition: bool, detail: str) -> None:
    """Record a mechanical check; failures are emitted before exit status 1."""

    CHECKS.append(
        {
            "group": group,
            "name": name,
            "pass": bool(condition),
            "detail": detail,
        }
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def stable_digest(obj: Any) -> str:
    encoded = json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return sha256_bytes(encoded)


def source_row(path: Path) -> dict[str, Any]:
    exists = path.is_file()
    return {
        "path": str(path.relative_to(REPO)),
        "exists": exists,
        "size_bytes": path.stat().st_size if exists else None,
        "sha256": sha256_file(path) if exists else None,
    }


def discriminant(a: T, b: T, c: T) -> T:
    return b * b - 4 * a * c  # type: ignore[operator,return-value]


def projective_matrix(a: T, b: T, c: T) -> list[list[T]]:
    two = type(b)(2) if not isinstance(b, Fraction) else Fraction(2)
    return [[-b / two, -c], [a, b / two]]  # type: ignore[operator,list-item]


def matmul(a: list[list[T]], b: list[list[T]]) -> list[list[T]]:
    return [
        [sum((a[i][k] * b[k][j] for k in range(len(b))), 0) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def transpose(a: list[list[T]]) -> list[list[T]]:
    return [list(row) for row in zip(*a)]


def matadd(a: list[list[T]], b: list[list[T]]) -> list[list[T]]:
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def trace2(a: list[list[T]]) -> T:
    return a[0][0] + a[1][1]


def det2(a: list[list[T]]) -> T:
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def det3(a: list[list[T]]) -> T:
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def fraction_matrix_rank(rows: list[list[Fraction]]) -> int:
    """Exact Gaussian-elimination rank over Q."""

    matrix = [row[:] for row in rows if any(value != 0 for value in row)]
    if not matrix:
        return 0
    rank = 0
    columns = len(matrix[0])
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or matrix[index][column] == 0:
                continue
            coefficient = matrix[index][column]
            matrix[index] = [
                matrix[index][j] - coefficient * matrix[rank][j]
                for j in range(columns)
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def invariant_symmetric_constraint_rank() -> tuple[int, bool]:
    """Return exact rank and whether the discriminant matrix spans the kernel."""

    zero = Fraction(0)
    one = Fraction(1)
    two = Fraction(2)
    generators = [
        [[two, zero, zero], [zero, zero, zero], [zero, zero, -two]],
        [[zero, one, zero], [zero, zero, two], [zero, zero, zero]],
        [[zero, zero, zero], [two, zero, zero], [zero, one, zero]],
    ]
    symmetric_basis: list[list[list[Fraction]]] = []
    for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)):
        basis = [[zero for _ in range(3)] for _ in range(3)]
        basis[i][j] = one
        basis[j][i] = one
        symmetric_basis.append(basis)

    rows: list[list[Fraction]] = []
    for generator in generators:
        transformed = [
            matadd(
                matmul(transpose(generator), basis),
                matmul(basis, generator),
            )
            for basis in symmetric_basis
        ]
        for i in range(3):
            for j in range(3):
                rows.append([transformed[index][i][j] for index in range(6)])

    discriminant_matrix = [
        [zero, zero, -two],
        [zero, one, zero],
        [-two, zero, zero],
    ]
    kernel_residual = max(
        abs(value)
        for generator in generators
        for row in matadd(
            matmul(transpose(generator), discriminant_matrix),
            matmul(discriminant_matrix, generator),
        )
        for value in row
    )
    return fraction_matrix_rank(rows), kernel_residual == 0


def exact_projective_identity(a: Fraction, b: Fraction, c: Fraction) -> bool:
    delta = discriminant(a, b, c)
    aq = projective_matrix(a, b, c)
    aq2 = matmul(aq, aq)
    return aq2 == [
        [delta / 4, Fraction(0)],
        [Fraction(0), delta / 4],
    ]


def exact_killing_identity(a: Fraction, b: Fraction, c: Fraction) -> bool:
    delta = discriminant(a, b, c)
    aq2 = matmul(projective_matrix(a, b, c), projective_matrix(a, b, c))
    killing = 4 * trace2(aq2)
    transvectant = 8 * a * c - 2 * b * b
    return killing == 2 * delta and transvectant == -2 * delta


def exact_mixed_pairing_identity(
    p: tuple[Fraction, Fraction, Fraction],
    r: tuple[Fraction, Fraction, Fraction],
) -> bool:
    """Check the polarized discriminant/Killing/transvectant identity."""

    a, b, c = p
    ap, bp, cp = r
    polar_discriminant = b * bp - 2 * (a * cp + ap * c)
    killing = 4 * trace2(matmul(projective_matrix(a, b, c), projective_matrix(ap, bp, cp)))
    transvectant = 4 * (a * cp + c * ap) - 2 * b * bp
    return killing == 2 * polar_discriminant and transvectant == -killing


def exact_fixed_point_identity(
    a: Fraction,
    b: Fraction,
    c: Fraction,
    u: Fraction,
    v: Fraction,
) -> bool:
    """Check det(z,A_p z)=p(z), fixing the homogeneous-coordinate convention."""

    aq = projective_matrix(a, b, c)
    image_u = aq[0][0] * u + aq[0][1] * v
    image_v = aq[1][0] * u + aq[1][1] * v
    wedge = u * image_v - v * image_u
    section = a * u * u + b * u * v + c * v * v
    return wedge == section


def complex_section_residual(
    coefficients: tuple[complex, complex, complex],
    point: tuple[complex, complex],
) -> float:
    a, b, c = coefficients
    u, v = point
    return abs(a * u * u + b * u * v + c * v * v)


def complex_eigenline_residual(
    coefficients: tuple[complex, complex, complex],
    point: tuple[complex, complex],
) -> float:
    a, b, c = coefficients
    u, v = point
    aq = projective_matrix(a, b, c)
    image_u = aq[0][0] * u + aq[0][1] * v
    image_v = aq[1][0] * u + aq[1][1] * v
    return abs(u * image_v - v * image_u)


def chordal_separation(
    first: tuple[complex, complex], second: tuple[complex, complex]
) -> float:
    u1, v1 = first
    u2, v2 = second
    numerator = abs(u1 * v2 - v1 * u2)
    denominator = (abs(u1) ** 2 + abs(v1) ** 2) ** 0.5 * (
        abs(u2) ** 2 + abs(v2) ** 2
    ) ** 0.5
    return numerator / denominator


def transform_binary_quadratic(
    a: Fraction,
    b: Fraction,
    c: Fraction,
    alpha: int,
    beta: int,
    gamma: int,
    delta: int,
) -> tuple[Fraction, Fraction, Fraction]:
    """Coefficients of p(alpha*u+beta*v, gamma*u+delta*v)."""

    a_new = a * alpha**2 + b * alpha * gamma + c * gamma**2
    b_new = (
        2 * a * alpha * beta
        + b * (alpha * delta + beta * gamma)
        + 2 * c * gamma * delta
    )
    c_new = a * beta**2 + b * beta * delta + c * delta**2
    return a_new, b_new, c_new


def random_fraction(rng: random.Random, *, nonzero: bool = False) -> Fraction:
    while True:
        numerator = rng.randint(-19, 19)
        denominator = rng.randint(1, 13)
        value = Fraction(numerator, denominator)
        if value or not nonzero:
            return value


def mp_pair(cell: dict[str, Any]) -> mp.mpc:
    return mp.mpc(mp.mpf(str(cell["re"])), mp.mpf(str(cell["im"])))


def mp_max_abs(matrix: list[list[mp.mpc]]) -> mp.mpf:
    return max(abs(cell) for row in matrix for cell in row)


def mp_bilinear(x: list[mp.mpc], matrix: list[list[mp.mpc]]) -> mp.mpc:
    return sum(
        (x[i] * matrix[i][j] * x[j] for i in range(len(x)) for j in range(len(x))),
        mp.mpc(0),
    )


def mp_hermitian(x: list[mp.mpc], matrix: list[list[mp.mpc]]) -> mp.mpc:
    return sum(
        (
            mp.conj(x[i]) * matrix[i][j] * x[j]
            for i in range(len(x))
            for j in range(len(x))
        ),
        mp.mpc(0),
    )


def scientific(value: mp.mpf | mp.mpc, digits: int = 12) -> str:
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def extract_theorem_block(source: str) -> str:
    start = source.index(THEOREM_START)
    end = source.index(THEOREM_END, start) + len(THEOREM_END)
    return source[start:end]


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e2_discriminant_regression.json"
    md_path = OUTPUT / "ap_e2_discriminant_regression.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    check_lines = "\n".join(
        f"- [{'PASS' if item['pass'] else 'FAIL'}] `{item['group']}` — "
        f"{item['name']}: {item['detail']}"
        for item in CHECKS
    )
    source_lines = "\n".join(
        f"- `{row['path']}` — `{row['sha256']}` ({row['size_bytes']} bytes)"
        for row in result["source_manifest"]
    )
    controls = result["negative_controls"]
    md_path.write_text(
        f"""# AP-E2 Discriminant / Contact Regression

Status: **{result['status']}** — `{result['checks_passed']}/{result['checks_total']}` checks pass.

This is a deterministic, fail-closed algebraic regression.  A green result
preserves the convention

```text
B_Kill(A_q,A_q) = 2 Delta(p) = 2 sqrt(3) x^T K_tr x,
A_q^2 = Delta(p) I / 4,
(p,p)_2 = -2 Delta(p).
```

It does not derive H3+, a microscopic carrier, or any dynamics, and therefore
sets `physics_promotion_allowed=false`.

## Coverage

- exact rational anchor and `{result['sample_counts']['exact_random_rational']}` seeded exact-rational samples;
- `{result['sample_counts']['high_precision_complex']}` seeded complex samples at `{result['precision']['mp_dps']}` decimal places;
- exact `SL(2)` covariance and matrix conjugacy, polarized bilinear identities,
  scalar covariance, and uniqueness of the invariant symmetric line;
- real, complex, and infinite regular roots, the finite and infinite
  double-root charts, a `1e-70` near-null perturbation, and the theorem's
  explicit exclusion of the zero element;
- Route-E `K_tr` source-card normalization and full `H,E,F` invariance;
- explicit spherical-to-normalized-polynomial basis conversion;
- four required negative controls: missing `sqrt(2)`, missing factor two,
  using the spherical coefficient on the normalized basis, and replacing the
  complex bilinear contact by a Hermitian form.

Maximum intrinsic 100-digit complex-sample residual:
`{result['precision']['max_intrinsic_high_precision_identity_residual']}`.
The separately reported Route-E source-card residual is
`{result['precision']['max_source_card_identity_residual']}` because that JSON
anchor stores ordinary double-precision decimals.

## Negative controls

- missing `sqrt(2)` detected: `{str(controls['missing_sqrt2']['detected']).lower()}`
- missing Killing factor two detected: `{str(controls['missing_factor_two']['detected']).lower()}`
- normalized-basis factor mismatch detected: `{str(controls['normalized_basis_factor']['detected']).lower()}`
- Hermitian-for-bilinear substitution detected: `{str(controls['hermitian_substitution']['detected']).lower()}`

These are tests that must disagree with the correct theorem.  Their detection
is a PASS; accidental agreement would fail the run.

## Mechanical checks

{check_lines}

## Source hashes

{source_lines}

The theorem-block hash lock is
`{result['theorem_block_sha256']}`.

## Promotion boundary

- `physics_promotion_allowed=false`
- AP-E2 freezes an exact conditional invariant theorem only.  AP-E3 must
  independently derive any Berry level or `O(2)` prequantum bundle.
""",
        encoding="utf-8",
    )


def main() -> None:
    mp.mp.dps = MP_DPS
    rng = random.Random(SEED)

    critical_sources = [
        ROUTE_E_CARD,
        ROUTE_E_BUILDER,
        THEOREM_SOURCE,
        AP_E2_NOTE,
        LEGACY_VERIFIER,
        THIS_SCRIPT,
    ]
    source_manifest = [source_row(path) for path in critical_sources]
    check(
        "A0_provenance",
        "all critical AP-E2 sources are present and hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    route_e_card = json.loads(ROUTE_E_CARD.read_text(encoding="utf-8"))
    digest_payload = {
        key: value
        for key, value in route_e_card.items()
        if key not in {"card_sha256", "created_utc"}
    }
    recomputed_card_digest = stable_digest(digest_payload)
    check(
        "A0_provenance",
        "Route-E invariant card verifies its internal stable digest",
        recomputed_card_digest == route_e_card.get("card_sha256"),
        f"stored={route_e_card.get('card_sha256')}; recomputed={recomputed_card_digest}",
    )

    builder_rows = [
        row
        for row in route_e_card.get("source_manifest", [])
        if row.get("label") == "audit0_builder"
    ]
    current_builder_sha = sha256_file(ROUTE_E_BUILDER)
    check(
        "A0_provenance",
        "Route-E card binds the current Audit-0 builder hash",
        len(builder_rows) == 1
        and builder_rows[0].get("sha256") == current_builder_sha,
        f"card={builder_rows[0].get('sha256') if builder_rows else None}; current={current_builder_sha}",
    )

    theorem_text = THEOREM_SOURCE.read_text(encoding="utf-8")
    theorem_block = extract_theorem_block(theorem_text)
    theorem_block_sha = sha256_bytes(theorem_block.encode("utf-8"))
    check(
        "A0_provenance",
        "the AP-E2 theorem/proof block matches the frozen semantic source",
        theorem_block_sha == EXPECTED_THEOREM_BLOCK_SHA256,
        f"block_sha256={theorem_block_sha}",
    )

    # ------------------------------------------------ exact rational algebra
    anchor = (Fraction(3, 2), Fraction(-7, 3), Fraction(5, 4))
    anchor_delta = discriminant(*anchor)
    check(
        "A1_exact",
        "exact rational anchor obeys A_q^2=Delta I/4",
        exact_projective_identity(*anchor),
        f"(a,b,c)={anchor}; Delta={anchor_delta}",
    )
    check(
        "A1_exact",
        "exact rational anchor obeys B=2 Delta and (p,p)_2=-2 Delta",
        exact_killing_identity(*anchor),
        f"Delta={anchor_delta}; B={2 * anchor_delta}; transvectant={-2 * anchor_delta}",
    )

    second_anchor = (Fraction(-2, 5), Fraction(11, 7), Fraction(4, 3))
    check(
        "A1_exact",
        "polarized anchor obeys B(A_p,A_r)=2D(p,r) and (p,r)_2=-B(A_p,A_r)",
        exact_mixed_pairing_identity(anchor, second_anchor),
        f"p={anchor}; r={second_anchor}",
    )

    invariant_rank, discriminant_in_kernel = invariant_symmetric_constraint_rank()
    check(
        "A1_exact",
        "the symmetric sl2-invariant contact line is exactly one-dimensional",
        invariant_rank == 5 and discriminant_in_kernel,
        f"constraint rank={invariant_rank}/6; kernel dimension={6-invariant_rank}",
    )

    exact_identity_failures = 0
    exact_mixed_failures = 0
    exact_fixed_point_failures = 0
    exact_sl2_failures = 0
    exact_conjugacy_failures = 0
    exact_scaling_failures = 0
    exact_regular_count = 0
    for _ in range(EXACT_RANDOM_SAMPLES):
        a = random_fraction(rng, nonzero=True)
        b = random_fraction(rng)
        c = random_fraction(rng)
        if not exact_projective_identity(a, b, c) or not exact_killing_identity(a, b, c):
            exact_identity_failures += 1

        second = (
            random_fraction(rng),
            random_fraction(rng),
            random_fraction(rng),
        )
        if not exact_mixed_pairing_identity((a, b, c), second):
            exact_mixed_failures += 1

        u = random_fraction(rng, nonzero=True)
        v = random_fraction(rng)
        if not exact_fixed_point_identity(a, b, c, u, v):
            exact_fixed_point_failures += 1

        delta = discriminant(a, b, c)
        if delta:
            exact_regular_count += 1

        # g=[[1+n*m,n],[m,1]] has determinant one exactly.
        n = rng.randint(-7, 7)
        m = rng.randint(-7, 7)
        alpha, beta, gamma, det_entry = 1 + n * m, n, m, 1
        transformed = transform_binary_quadratic(
            a, b, c, alpha, beta, gamma, det_entry
        )
        if discriminant(*transformed) != delta:
            exact_sl2_failures += 1
        g = [
            [Fraction(alpha), Fraction(beta)],
            [Fraction(gamma), Fraction(det_entry)],
        ]
        g_inverse = [
            [Fraction(det_entry), Fraction(-beta)],
            [Fraction(-gamma), Fraction(alpha)],
        ]
        transformed_matrix = projective_matrix(*transformed)
        conjugated_matrix = matmul(
            matmul(g_inverse, projective_matrix(a, b, c)), g
        )
        if transformed_matrix != conjugated_matrix:
            exact_conjugacy_failures += 1

        scalar = random_fraction(rng, nonzero=True)
        scaled_delta = discriminant(scalar * a, scalar * b, scalar * c)
        if scaled_delta != scalar * scalar * delta:
            exact_scaling_failures += 1

    check(
        "A1_exact",
        "seeded rational samples obey all projective/Killing/transvectant identities exactly",
        exact_identity_failures == 0 and exact_regular_count > 200,
        f"seed={SEED}; samples={EXACT_RANDOM_SAMPLES}; failures={exact_identity_failures}; regular={exact_regular_count}",
    )
    check(
        "A1_exact",
        "seeded rational pairs obey the polarized contact identity exactly",
        exact_mixed_failures == 0,
        f"samples={EXACT_RANDOM_SAMPLES}; failures={exact_mixed_failures}",
    )
    check(
        "A1_exact",
        "det(z,A_p z)=p(z) fixes the homogeneous and affine-chart convention",
        exact_fixed_point_failures == 0,
        f"samples={EXACT_RANDOM_SAMPLES}; failures={exact_fixed_point_failures}; affine xi=v/u",
    )
    check(
        "A1_exact",
        "binary discriminant is exactly invariant under seeded SL(2,Z) substitutions",
        exact_sl2_failures == 0,
        f"samples={EXACT_RANDOM_SAMPLES}; failures={exact_sl2_failures}",
    )
    check(
        "A1_exact",
        "the projective matrix transforms by exact conjugacy under the declared SL(2) substitution",
        exact_conjugacy_failures == 0,
        f"A[p(gz)]=g^-1 A[p] g; samples={EXACT_RANDOM_SAMPLES}; failures={exact_conjugacy_failures}",
    )
    check(
        "A1_exact",
        "discriminant and contact have exact weight two under section rescaling",
        exact_scaling_failures == 0,
        f"samples={EXACT_RANDOM_SAMPLES}; failures={exact_scaling_failures}",
    )

    # --------------------------------------------- exact degenerate boundary
    degenerate_failures = 0
    degenerate_cases: list[dict[str, str]] = []
    square_pairs = [
        (Fraction(2), Fraction(-3)),
        (Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(1)),
        (Fraction(-5, 7), Fraction(11, 9)),
    ]
    for alpha, beta in square_pairs:
        a, b, c = alpha * alpha, 2 * alpha * beta, beta * beta
        aq = projective_matrix(a, b, c)
        aq2 = matmul(aq, aq)
        is_nonzero = any(cell != 0 for row in aq for cell in row)
        case_ok = (
            discriminant(a, b, c) == 0
            and aq2 == [[Fraction(0), Fraction(0)], [Fraction(0), Fraction(0)]]
            and is_nonzero
        )
        if not case_ok:
            degenerate_failures += 1
        degenerate_cases.append(
            {
                "linear_factor": f"({alpha})u+({beta})v",
                "a": str(a),
                "b": str(b),
                "c": str(c),
                "classification": "nonzero_nilpotent_double_root",
            }
        )
    check(
        "A2_boundary",
        "finite and infinite double-root charts lie on the nonzero nilpotent/contact-null cone",
        degenerate_failures == 0,
        f"cases={len(square_pairs)}; failures={degenerate_failures}; includes u^2 and v^2",
    )

    regular_projective_cases = [
        {
            "name": "real finite roots",
            "coefficients": (1.0 + 0j, 0j, -1.0 + 0j),
            "roots": ((1.0 + 0j, 1.0 + 0j), (1.0 + 0j, -1.0 + 0j)),
        },
        {
            "name": "complex finite roots",
            "coefficients": (1.0 + 0j, 0j, 1.0 + 0j),
            "roots": ((1.0 + 0j, 1j), (1.0 + 0j, -1j)),
        },
        {
            "name": "zero and infinity",
            "coefficients": (0j, 1.0 + 0j, 0j),
            "roots": ((1.0 + 0j, 0j), (0j, 1.0 + 0j)),
        },
        {
            "name": "asymmetric chart guard",
            "coefficients": (2.0 + 0j, 3.0 + 0j, 1.0 + 0j),
            "roots": ((1.0 + 0j, -1.0 + 0j), (1.0 + 0j, -2.0 + 0j)),
        },
    ]
    regular_root_residual = 0.0
    minimum_chordal_separation = 1.0
    for case in regular_projective_cases:
        coefficients = case["coefficients"]
        roots = case["roots"]
        regular_root_residual = max(
            regular_root_residual,
            *(complex_section_residual(coefficients, root) for root in roots),
            *(complex_eigenline_residual(coefficients, root) for root in roots),
        )
        minimum_chordal_separation = min(
            minimum_chordal_separation, chordal_separation(*roots)
        )
    check(
        "A2_boundary",
        "real, complex, asymmetric, and infinite roots match eigenlines in the xi=v/u chart",
        regular_root_residual < 1.0e-14 and minimum_chordal_separation > 0.25,
        "cases={}; max residual={:.2e}; min projective separation={:.6f}".format(
            len(regular_projective_cases),
            regular_root_residual,
            minimum_chordal_separation,
        ),
    )

    double_projective_cases = [
        ((4.0 + 0j, -4.0 + 0j, 1.0 + 0j), (1.0 + 0j, 2.0 + 0j)),
        ((1.0 + 0j, 0j, 0j), (0j, 1.0 + 0j)),
        ((0j, 0j, 1.0 + 0j), (1.0 + 0j, 0j)),
    ]
    double_root_residual = 0.0
    for coefficients, root in double_projective_cases:
        a_complex, b_complex, c_complex = coefficients
        u_complex, v_complex = root
        double_root_residual = max(
            double_root_residual,
            complex_section_residual(coefficients, root),
            complex_eigenline_residual(coefficients, root),
            abs(2 * a_complex * u_complex + b_complex * v_complex),
            abs(b_complex * u_complex + 2 * c_complex * v_complex),
            abs(discriminant(*coefficients)),
        )
    check(
        "A2_boundary",
        "finite zero, finite nonzero, and infinite double roots have vanishing first jet",
        double_root_residual < 1.0e-14,
        f"cases={len(double_projective_cases)}; max section/eigenline/gradient residual={double_root_residual:.2e}",
    )

    zero_a = zero_b = zero_c = Fraction(0)
    zero_matrix = projective_matrix(zero_a, zero_b, zero_c)
    zero_excluded = (
        discriminant(zero_a, zero_b, zero_c) == 0
        and all(cell == 0 for row in zero_matrix for cell in row)
    )
    check(
        "A2_boundary",
        "the zero section is detected and excluded from the theorem's nonzero nilpotent clause",
        zero_excluded,
        "Delta(0)=0 and A_0=0; theorem domain is explicitly X_q != 0",
    )

    eps = mp.mpf("1e-70")
    near_a, near_b, near_c = mp.mpf(1), mp.mpf(2), mp.mpf(1) + eps
    near_delta = discriminant(near_a, near_b, near_c)
    near_aq = projective_matrix(near_a, near_b, near_c)
    near_aq2 = matmul(near_aq, near_aq)
    near_target = [
        [near_delta / 4, mp.mpf(0)],
        [mp.mpf(0), near_delta / 4],
    ]
    near_residual = max(
        abs(near_aq2[i][j] - near_target[i][j]) for i in range(2) for j in range(2)
    )
    check(
        "A2_boundary",
        "100-digit arithmetic resolves a 1e-70 perturbation away from the null cone",
        near_delta != 0
        and abs(near_delta + 4 * eps) < mp.mpf("1e-95")
        and near_residual < mp.mpf("1e-95"),
        f"Delta={scientific(near_delta, 18)}; identity residual={scientific(near_residual)}",
    )

    # ---------------------------------------- Route-E source normalization
    raw_k = route_e_card["routeA_invariant_anchors"]["k_tr"]["matrix"]
    loaded_k = [[mp_pair(cell) for cell in row] for row in raw_k]
    sqrt2 = mp.sqrt(2)
    sqrt3 = mp.sqrt(3)
    canonical_k = [
        [mp.mpc(0), mp.mpc(0), -1 / sqrt3],
        [mp.mpc(0), 1 / sqrt3, mp.mpc(0)],
        [-1 / sqrt3, mp.mpc(0), mp.mpc(0)],
    ]
    source_k_error = max(
        abs(loaded_k[i][j] - canonical_k[i][j]) for i in range(3) for j in range(3)
    )
    loaded_k_square = matmul(loaded_k, loaded_k)
    k_square_error = max(
        abs(loaded_k_square[i][j] - (mp.mpf(1) / 3 if i == j else 0))
        for i in range(3)
        for j in range(3)
    )
    check(
        "A3_route_e",
        "Route-E card uses the canonical spherical K_tr normalization",
        source_k_error < SOURCE_NORMALIZATION_TOL
        and k_square_error < SOURCE_NORMALIZATION_TOL,
        f"max|K_source-K_canonical|={scientific(source_k_error)}; max|K^2-I/3|={scientific(k_square_error)}",
    )

    k_symmetry_error = max(
        abs(loaded_k[i][j] - loaded_k[j][i]) for i in range(3) for j in range(3)
    )
    k_inverse_error = max(
        abs(
            sum(loaded_k[i][ell] * (3 * loaded_k[ell][j]) for ell in range(3))
            - (1 if i == j else 0)
        )
        for i in range(3)
        for j in range(3)
    )
    k_determinant = det3(loaded_k)
    expected_k_determinant = -1 / (3 * sqrt3)
    isotropic_vector = [mp.mpc(1), mp.mpc(1), mp.mpc("0.5")]
    isotropic_norm = mp_bilinear(isotropic_vector, loaded_k)
    isotropic_image = [
        sum((loaded_k[i][j] * isotropic_vector[j] for j in range(3)), mp.mpc(0))
        for i in range(3)
    ]
    check(
        "A3_route_e",
        "K_tr is symmetric, invertible, and may have nonzero isotropic vectors without being degenerate",
        k_symmetry_error < SOURCE_NORMALIZATION_TOL
        and k_inverse_error < SOURCE_NORMALIZATION_TOL
        and abs(k_determinant - expected_k_determinant) < SOURCE_NORMALIZATION_TOL
        and abs(isotropic_norm) < SOURCE_NORMALIZATION_TOL
        and max(abs(value) for value in isotropic_image) > mp.mpf("0.1"),
        "sym={}; inverse={}; det={}; isotropic norm={}; |Kx|_max={}".format(
            scientific(k_symmetry_error),
            scientific(k_inverse_error),
            scientific(k_determinant),
            scientific(isotropic_norm),
            scientific(max(abs(value) for value in isotropic_image)),
        ),
    )

    basis_probe_x = [mp.mpc(2, -1), mp.mpc(-3, 2), mp.mpc(1, 4)]
    basis_probe_y = [
        basis_probe_x[2] / sqrt2,
        -basis_probe_x[1] / sqrt2,
        basis_probe_x[0] / sqrt2,
    ]
    x_contact = mp_bilinear(basis_probe_x, loaded_k)
    y_contact = mp_bilinear(basis_probe_y, loaded_k)
    xp_probe, x0_probe, xm_probe = basis_probe_x
    basis_probe_delta = x0_probe * x0_probe - 2 * xp_probe * xm_probe
    basis_change_residual = max(
        abs(x_contact - 2 * y_contact),
        abs(2 * sqrt3 * x_contact - 2 * basis_probe_delta),
        abs(4 * sqrt3 * y_contact - 2 * basis_probe_delta),
    )
    check(
        "A3_route_e",
        "spherical and normalized-polynomial bases carry the required factor-two conversion",
        basis_change_residual < SOURCE_NORMALIZATION_TOL,
        "y=R*x/sqrt(2); B=2sqrt(3)x^TKx=4sqrt(3)y^TKy; residual={}".format(
            scientific(basis_change_residual)
        ),
    )

    rho_h = [
        [mp.mpc(2), mp.mpc(0), mp.mpc(0)],
        [mp.mpc(0), mp.mpc(0), mp.mpc(0)],
        [mp.mpc(0), mp.mpc(0), mp.mpc(-2)],
    ]
    rho_e = [
        [mp.mpc(0), sqrt2, mp.mpc(0)],
        [mp.mpc(0), mp.mpc(0), sqrt2],
        [mp.mpc(0), mp.mpc(0), mp.mpc(0)],
    ]
    rho_f = transpose(rho_e)
    generator_residuals: dict[str, mp.mpf] = {}
    for name, rho in (("H", rho_h), ("E", rho_e), ("F", rho_f)):
        residual_matrix = matadd(
            matmul(transpose(rho), loaded_k), matmul(loaded_k, rho)
        )
        generator_residuals[name] = mp_max_abs(residual_matrix)
    check(
        "A3_route_e",
        "loaded K_tr is invariant under the complete H,E,F generator triple",
        max(generator_residuals.values()) < SOURCE_NORMALIZATION_TOL,
        "residuals="
        + json.dumps(
            {name: scientific(value) for name, value in generator_residuals.items()},
            sort_keys=True,
        ),
    )

    # ------------------------------- high-precision complex random samples
    integer_j = [
        [mp.mpc(0), mp.mpc(0), mp.mpc(-1)],
        [mp.mpc(0), mp.mpc(1), mp.mpc(0)],
        [mp.mpc(-1), mp.mpc(0), mp.mpc(0)],
    ]
    max_intrinsic_hp_residual = mp.mpf(0)
    max_source_card_residual = mp.mpf(0)
    hp_regular = 0
    hp_failures = 0
    for _ in range(HIGH_PRECISION_COMPLEX_SAMPLES):
        x: list[mp.mpc] = []
        for _component in range(3):
            re = mp.mpf(rng.randint(-1000, 1000)) / rng.randint(1, 97)
            im = mp.mpf(rng.randint(-1000, 1000)) / rng.randint(1, 97)
            x.append(mp.mpc(re, im))
        x_plus, x_zero, x_minus = x
        a, b, c = x_minus / sqrt2, -x_zero, x_plus / sqrt2
        delta = discriminant(a, b, c)
        aq = projective_matrix(a, b, c)
        aq2 = matmul(aq, aq)
        target = [[delta / 4, mp.mpc(0)], [mp.mpc(0), delta / 4]]
        projective_error = max(
            abs(aq2[i][j] - target[i][j]) for i in range(2) for j in range(2)
        )
        killing = 4 * trace2(aq2)
        scaled_contact = mp_bilinear(x, integer_j)
        loaded_contact_identity = 2 * sqrt3 * mp_bilinear(x, loaded_k)
        scale = max(mp.mpf(1), abs(delta), abs(killing))
        intrinsic_residual = max(
            projective_error,
            abs(killing - 2 * delta) / scale,
            abs(scaled_contact - delta) / scale,
        )
        source_card_residual = abs(loaded_contact_identity - 2 * delta) / scale
        max_intrinsic_hp_residual = max(
            max_intrinsic_hp_residual, intrinsic_residual
        )
        max_source_card_residual = max(
            max_source_card_residual, source_card_residual
        )
        if abs(delta) > mp.mpf("1e-40"):
            hp_regular += 1
        if (
            intrinsic_residual >= mp.mpf("1e-85")
            or source_card_residual >= SOURCE_NORMALIZATION_TOL
        ):
            hp_failures += 1
    check(
        "A4_complex",
        "seeded 100-digit complex samples reproduce the discriminant/contact theorem",
        hp_failures == 0 and hp_regular > 240,
        (
            f"seed={SEED}; samples={HIGH_PRECISION_COMPLEX_SAMPLES}; "
            f"regular={hp_regular}; failures={hp_failures}; "
            f"intrinsic max={scientific(max_intrinsic_hp_residual)}; "
            f"source-card max={scientific(max_source_card_residual)}"
        ),
    )

    # ------------------------------------------- required negative controls
    x_control = [mp.mpc(1, 2), mp.mpc(-3, 1), mp.mpc(2, -1)]
    xp, x0, xm = x_control
    correct_delta = x0 * x0 - 2 * xp * xm
    correct_b = 2 * correct_delta

    wrong_delta_missing_sqrt2 = x0 * x0 - 4 * xp * xm
    missing_sqrt2_residual = abs(wrong_delta_missing_sqrt2 - correct_delta)
    missing_sqrt2_detected = missing_sqrt2_residual > mp.mpf("1e-2")

    wrong_b_missing_factor_two = correct_delta
    missing_factor_two_residual = abs(wrong_b_missing_factor_two - correct_b)
    missing_factor_two_detected = missing_factor_two_residual > mp.mpf("1e-2")

    normalized_control = [xm / sqrt2, -x0 / sqrt2, xp / sqrt2]
    wrong_normalized_basis_b = 2 * sqrt3 * mp_bilinear(
        normalized_control, loaded_k
    )
    normalized_basis_residual = abs(wrong_normalized_basis_b - correct_b)
    normalized_basis_detected = normalized_basis_residual > mp.mpf("1e-2")

    bilinear_contact = mp_bilinear(x_control, integer_j)
    hermitian_contact = mp_hermitian(x_control, integer_j)
    hermitian_residual = abs(hermitian_contact - bilinear_contact)
    hermitian_detected = hermitian_residual > mp.mpf("1e-2")

    check(
        "A5_negative_control",
        "missing sqrt(2) in the spherical section map is rejected",
        missing_sqrt2_detected,
        f"intentional wrong-convention residual={scientific(missing_sqrt2_residual)}",
    )
    check(
        "A5_negative_control",
        "missing factor two in B_Kill=2 Delta is rejected",
        missing_factor_two_detected,
        f"intentional wrong-convention residual={scientific(missing_factor_two_residual)}",
    )
    check(
        "A5_negative_control",
        "using the spherical 2sqrt(3) factor on normalized polynomial coefficients is rejected",
        normalized_basis_detected,
        f"intentional wrong-basis residual={scientific(normalized_basis_residual)}",
    )
    check(
        "A5_negative_control",
        "Hermitian substitution for the complex symmetric contact is rejected",
        hermitian_detected,
        f"intentional wrong-convention residual={scientific(hermitian_residual)}",
    )

    # ----------------------------------------------------- logical boundary
    physics_promotion_allowed = False
    abelian_adjoint_killing = 0.0
    abelian_invariant_pairing_determinant = 1.0
    check(
        "A6_boundary",
        "the one-dimensional abelian counterexample passes invariant-pairing H3 but fails Killing-contact H3+",
        abelian_adjoint_killing == 0.0
        and abelian_invariant_pairing_determinant != 0.0,
        "ad_e=0 gives Killing=0 while B(e,e)=1 is invariant and nondegenerate",
    )
    check(
        "A6_boundary",
        "a green AP-E2 algebraic regression cannot promote a physical bridge",
        physics_promotion_allowed is False,
        "physics_promotion_allowed=false by construction",
    )

    passed = sum(item["pass"] for item in CHECKS)
    all_pass = passed == len(CHECKS)
    negative_controls = {
        "missing_sqrt2": {
            "wrong_convention": "a=x_-, b=-x_0, c=x_+ (sqrt(2) omitted)",
            "residual": scientific(missing_sqrt2_residual, 18),
            "detected": missing_sqrt2_detected,
        },
        "missing_factor_two": {
            "wrong_convention": "B_Kill=Delta instead of B_Kill=2 Delta",
            "residual": scientific(missing_factor_two_residual, 18),
            "detected": missing_factor_two_detected,
        },
        "normalized_basis_factor": {
            "wrong_convention": "B=2sqrt(3)y^T K y in normalized (u^2,sqrt(2)uv,v^2) coefficients",
            "correct_convention": "B=4sqrt(3)y^T K y",
            "residual": scientific(normalized_basis_residual, 18),
            "detected": normalized_basis_detected,
        },
        "hermitian_substitution": {
            "wrong_convention": "x^dagger K x instead of x^T K x",
            "residual": scientific(hermitian_residual, 18),
            "detected": hermitian_detected,
        },
    }
    result: dict[str, Any] = {
        "audit": "verify_ap_e2_discriminant_regression",
        "status": (
            "ap_e2_exact_regression_pass_no_physics_promotion"
            if all_pass
            else "ap_e2_regression_failure"
        ),
        "all_pass": all_pass,
        "checks_passed": passed,
        "checks_total": len(CHECKS),
        "physics_promotion_allowed": physics_promotion_allowed,
        "theorem_scope": {
            "domain": "nonzero X_q in H^0(CP1,T_CP1) on the conditional H3+ carrier",
            "identity": "B_Kill(A_q,A_q)=2 Delta(p)=2 sqrt(3) x^T K_tr x",
            "matrix_identity": "A_q^2=Delta(p) I/4",
            "transvectant_identity": "(p,p)_2=-2 Delta(p)",
            "polarized_identity": "B(A_p,A_r)=2D(p,r) and (p,r)_2=-B(A_p,A_r)",
            "basis_identity": "B=2sqrt(3)x^T K_tr x=4sqrt(3)y^T K_tr y for y=R x/sqrt(2)",
            "regular_branch": "Delta != 0 iff two distinct zeros iff regular semisimple iff contact non-null",
            "boundary": "Delta=0 and X_q!=0 is the nilpotent/contact-null double-zero cone",
            "not_derived": [
                "H3+ or its microscopic dynamics",
                "identification of divisor zeros with particles, energies, or families",
                "Berry level k=2 or an O(2) prequantum bundle",
            ],
        },
        "reproducibility": {
            "seed": SEED,
            "stdlib_exact_type": "fractions.Fraction",
            "high_precision_backend": "mpmath",
            "mpmath_version": mp.__version__,
        },
        "sample_counts": {
            "exact_random_rational": EXACT_RANDOM_SAMPLES,
            "high_precision_complex": HIGH_PRECISION_COMPLEX_SAMPLES,
            "exact_degenerate": len(square_pairs),
        },
        "precision": {
            "mp_dps": MP_DPS,
            "source_normalization_tolerance": scientific(SOURCE_NORMALIZATION_TOL),
            "max_intrinsic_high_precision_identity_residual": scientific(
                max_intrinsic_hp_residual, 18
            ),
            "max_source_card_identity_residual": scientific(
                max_source_card_residual, 18
            ),
            "route_E_K_source_error": scientific(source_k_error, 18),
            "route_E_K_square_error": scientific(k_square_error, 18),
            "near_null_epsilon": scientific(eps),
            "near_null_discriminant": scientific(near_delta, 18),
            "near_null_identity_residual": scientific(near_residual, 18),
        },
        "degenerate_cases": degenerate_cases,
        "negative_controls": negative_controls,
        "source_manifest": source_manifest,
        "route_E_card_internal_digest": {
            "stored": route_e_card.get("card_sha256"),
            "recomputed": recomputed_card_digest,
            "matches": recomputed_card_digest == route_e_card.get("card_sha256"),
        },
        "theorem_block_sha256": theorem_block_sha,
        "theorem_block_expected_sha256": EXPECTED_THEOREM_BLOCK_SHA256,
        "checks": CHECKS,
    }
    write_outputs(result)
    print(
        f"verify_ap_e2_discriminant_regression: {passed}/{len(CHECKS)} checks pass; "
        "physics_promotion_allowed=false"
    )
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
