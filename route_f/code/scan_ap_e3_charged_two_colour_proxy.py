#!/usr/bin/env python3
"""Deterministic nonlinear-EFT proxy for the charged two-colour AP-E3 gate.

This program is deliberately *not* a lattice simulation of two-colour QCD.
It studies a declared low-energy proxy with a charge-two complex diquark and
an SU(2) mesonic Skyrme subsector.  It audits

* the homogeneous vacuum and its six-by-six Hessian;
* pion, radial-meson, and charged-diquark gaps;
* a nonlinear B=1 hedgehog boundary-value problem;
* Derrick's virial identity;
* the spherically symmetric meson Hessian;
* the most dangerous charged-diquark Hessian channel; and
* grid/volume convergence plus unstable negative controls.

A green run closes only a necessary classical-EFT gate.  It cannot establish
the strong phase of the microscopic gauge theory, the full three-dimensional
soliton Hessian, quantum determinant stability, or lattice QC2D closure.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy.integrate import simpson, solve_bvp
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh


REPO = Path(__file__).resolve().parents[2]
ROUTE_F = Path(__file__).resolve().parents[1]
OUTPUT = ROUTE_F / "output"
TEX = ROUTE_F / "tex" / "ap_e3_charged_two_colour_soliton_proxy.tex"
BIB = ROUTE_F / "tex" / "ap_e3_charged_two_colour_soliton_proxy.bib"
THIS_SCRIPT = Path(__file__).resolve()

MODEL_NAME = "charged-two-colour nonlinear EFT proxy; not lattice QC2D"
EPS = 1.0e-4
BVP_TOL = 3.0e-7


@dataclass(frozen=True)
class Model:
    """Dimensionless benchmark; f=e=1 fixes the displayed energy units."""

    f: float = 1.0
    e: float = 1.0
    lam: float = 6.0
    m_pi: float = 0.5
    mu_lift_sq: float = 0.56
    chi: float = -2.0
    q_diquark: int = 2

    @property
    def h(self) -> float:
        return self.f * self.m_pi**2

    @property
    def v_sq(self) -> float:
        return self.f**2 - self.m_pi**2 / self.lam

    @property
    def m_sigma_sq(self) -> float:
        return self.m_pi**2 + 2.0 * self.lam * self.f**2

    @property
    def m_diquark_sq(self) -> float:
        return self.m_pi**2 + self.mu_lift_sq

    @property
    def beta(self) -> float:
        return self.m_pi / (self.e * self.f)

    @property
    def alpha(self) -> float:
        return math.sqrt(self.m_diquark_sq) / (self.e * self.f)


CHECKS: list[dict[str, Any]] = []


def check(group: str, name: str, condition: bool, detail: str) -> None:
    CHECKS.append(
        {"group": group, "name": name, "pass": bool(condition), "detail": detail}
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_row(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPO)),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def potential(vector: np.ndarray, model: Model) -> float:
    """Linear-sigma proxy potential in canonical real coordinates.

    vector=(sigma, pi_1,pi_2,pi_3,d_1,d_2) and
    Delta=(d_1+i d_2)/sqrt(2), with U(1)_g charge two.
    """

    sigma = float(vector[0])
    radius_sq = float(np.dot(vector, vector))
    diquark_sq = float(np.dot(vector[4:], vector[4:]))
    return (
        0.25 * model.lam * (radius_sq - model.v_sq) ** 2
        - model.h * sigma
        + 0.5 * model.mu_lift_sq * diquark_sq
    )


def finite_difference_hessian(vector: np.ndarray, model: Model, step: float) -> np.ndarray:
    size = len(vector)
    hessian = np.zeros((size, size), dtype=float)
    centre = potential(vector, model)
    for first in range(size):
        unit_first = np.zeros(size)
        unit_first[first] = step
        hessian[first, first] = (
            potential(vector + unit_first, model)
            - 2.0 * centre
            + potential(vector - unit_first, model)
        ) / step**2
        for second in range(first + 1, size):
            unit_second = np.zeros(size)
            unit_second[second] = step
            value = (
                potential(vector + unit_first + unit_second, model)
                - potential(vector + unit_first - unit_second, model)
                - potential(vector - unit_first + unit_second, model)
                + potential(vector - unit_first - unit_second, model)
            ) / (4.0 * step**2)
            hessian[first, second] = value
            hessian[second, first] = value
    return hessian


def homogeneous_phase(mu_lift_sq: float, model: Model) -> dict[str, Any]:
    """Analytic homogeneous minimum on both sides of the diquark onset."""

    critical = -model.m_pi**2
    if mu_lift_sq >= critical:
        sigma = model.f
        diquark = 0.0
        phase = "mesonic_U1g_unbroken" if mu_lift_sq > critical else "critical"
    else:
        sigma = -model.h / mu_lift_sq
        radius_sq = model.v_sq - mu_lift_sq / model.lam
        diquark_sq = max(0.0, radius_sq - sigma**2)
        diquark = math.sqrt(diquark_sq)
        phase = "diquark_condensed_U1g_broken"
    probe = np.array([sigma, 0.0, 0.0, 0.0, diquark, 0.0])
    local_model = Model(
        f=model.f,
        e=model.e,
        lam=model.lam,
        m_pi=model.m_pi,
        mu_lift_sq=mu_lift_sq,
        chi=model.chi,
        q_diquark=model.q_diquark,
    )
    eigenvalues = np.linalg.eigvalsh(
        finite_difference_hessian(probe, local_model, step=5.0e-5)
    )
    return {
        "mu_lift_sq": mu_lift_sq,
        "critical_mu_lift_sq": critical,
        "phase": phase,
        "sigma": sigma,
        "diquark_magnitude_real_basis": diquark,
        "potential": potential(probe, local_model),
        "hessian_eigenvalues": eigenvalues.tolist(),
    }


def charged_coordinate_scan(model: Model) -> dict[str, Any]:
    """Scan explicit (e_g,v_g/Lambda_c) coordinates under a declared ansatz.

    The matching relation is not derived from the microscopic gauge theory:
      mu_lift^2/Lambda_c^2 = mu_PG,0^2 + c_g (e_g v_g/Lambda_c)^2.
    It is only a reproducible coordinate chart for testing the EFT phase gate.
    We set Lambda_c=f=1 in the benchmark units.
    """

    mu_pg_baseline_sq = -0.40
    matching_coefficient = 0.24
    gauge_couplings = [0.0, 0.25, 0.50, 0.75, 1.00]
    higgs_ratios = [0.50, 1.00, 2.00, 3.00]
    rows: list[dict[str, Any]] = []
    for gauge_coupling in gauge_couplings:
        for higgs_ratio in higgs_ratios:
            mu_lift_sq = mu_pg_baseline_sq + matching_coefficient * (
                gauge_coupling * higgs_ratio
            ) ** 2
            phase = homogeneous_phase(mu_lift_sq, model)
            rows.append(
                {
                    "e_g": gauge_coupling,
                    "v_g_over_Lambda_c": higgs_ratio,
                    "matched_mu_lift_sq_over_Lambda_c_sq": mu_lift_sq,
                    "phase": phase["phase"],
                    "meson_condensate_sigma_over_Lambda_c": phase["sigma"],
                    "diquark_condensate_real_basis_over_Lambda_c": phase[
                        "diquark_magnitude_real_basis"
                    ],
                    "minimum_vacuum_hessian_over_Lambda_c_sq": min(
                        phase["hessian_eigenvalues"]
                    ),
                    "vacuum_hessian_eigenvalues_over_Lambda_c_sq": phase[
                        "hessian_eigenvalues"
                    ],
                }
            )
    return {
        "matching_status": "declared_EFT_ansatz_not_UV_derived",
        "Lambda_c_over_f": 1.0,
        "mu_PG_baseline_sq_over_Lambda_c_sq": mu_pg_baseline_sq,
        "matching_coefficient_c_g": matching_coefficient,
        "critical_e_g_times_v_g_over_Lambda_c": math.sqrt(
            (-model.m_pi**2 - mu_pg_baseline_sq) / matching_coefficient
        ),
        "rows": rows,
    }


def origin_cubic(slope: float, beta: float) -> float:
    """Coefficient b in F=pi-a*x+b*x^3+O(x^5)."""

    numerator = slope * (2.0 * slope**4 + 4.0 * slope**2 + 3.0 * beta**2)
    return numerator / (30.0 * (1.0 + 2.0 * slope**2))


def profile_rhs(x: np.ndarray, y: np.ndarray, beta: float) -> np.ndarray:
    profile, derivative = y
    sine = np.sin(profile)
    cosine = np.cos(profile)
    leading = x**2 + 2.0 * sine**2
    second = (
        -2.0 * x * derivative
        - 2.0 * sine * cosine * (derivative**2 - 1.0)
        + 2.0 * sine**3 * cosine / x**2
        + beta**2 * x**2 * sine
    ) / leading
    return np.vstack((derivative, second))


def solve_hedgehog(beta: float, box: float) -> Any:
    mesh = np.linspace(EPS, box, max(500, int(28 * box)))
    radius_guess = 1.1
    guess = 2.0 * np.arctan(radius_guess / mesh) * (box - mesh) / (box - EPS)
    derivative_guess = np.gradient(guess, mesh)

    def ode(x: np.ndarray, y: np.ndarray, _parameter: np.ndarray) -> np.ndarray:
        return profile_rhs(x, y, beta)

    def boundary(
        left: np.ndarray, right: np.ndarray, parameter: np.ndarray
    ) -> np.ndarray:
        slope = float(parameter[0])
        cubic = origin_cubic(slope, beta)
        return np.array(
            [
                left[0] - (math.pi - slope * EPS + cubic * EPS**3),
                left[1] - (-slope + 3.0 * cubic * EPS**2),
                right[0],
            ]
        )

    solution = solve_bvp(
        ode,
        boundary,
        mesh,
        np.vstack((guess, derivative_guess)),
        p=np.array([2.0]),
        tol=BVP_TOL,
        max_nodes=50000,
        verbose=0,
    )
    if solution.status != 0:
        raise RuntimeError(f"hedgehog BVP failed: {solution.message}")
    return solution


def profile_observables(solution: Any, beta: float, box: float) -> dict[str, Any]:
    x = np.linspace(EPS, box, 12001)
    profile, derivative = solution.sol(x)
    sine = np.sin(profile)
    e_two_density = 0.5 * x**2 * derivative**2 + sine**2
    e_four_density = sine**2 * derivative**2 + 0.5 * sine**4 / x**2
    e_zero_density = beta**2 * x**2 * (1.0 - np.cos(profile))
    e_two = float(simpson(e_two_density, x=x))
    e_four = float(simpson(e_four_density, x=x))
    e_zero = float(simpson(e_zero_density, x=x))
    baryon = float(simpson(-2.0 * sine**2 * derivative / math.pi, x=x))
    virial = e_two - e_four + 3.0 * e_zero
    scale_curvature = 2.0 * e_four + 6.0 * e_zero
    return {
        "origin_slope": float(solution.p[0]),
        "bvp_nodes": int(len(solution.x)),
        "bvp_max_rms_residual": float(np.max(solution.rms_residuals)),
        "B_numeric": baryon,
        "B_boundary": 1.0,
        "B_residual": abs(baryon - 1.0),
        "energy_dimensionless_integral": e_two + e_four + e_zero,
        "energy_units_4pi_f_over_e": 4.0 * math.pi * (e_two + e_four + e_zero),
        "E2": e_two,
        "E4": e_four,
        "E0": e_zero,
        "derrick_virial": virial,
        "derrick_relative_residual": abs(virial) / (e_two + e_four + e_zero),
        "scale_second_derivative": scale_curvature,
        "no_skyrme_scale_derivative": e_two + 3.0 * e_zero,
    }


def radial_hessian_eigenvalues(
    solution: Any, beta: float, box: float, intervals: int, count: int = 4
) -> np.ndarray:
    """Generalized radial meson Hessian H eta=omega^2 A eta via linear FEM."""

    x = np.linspace(EPS, box, intervals + 1)
    profile, derivative = solution.sol(x)
    second = solution.sol(x, 1)[1]
    sine = np.sin(profile)
    cosine = np.cos(profile)
    cosine_two = np.cos(2.0 * profile)
    weight = x**2 + 2.0 * sine**2
    potential_second = (
        2.0 * cosine_two
        + 2.0 * (3.0 * sine**2 * cosine**2 - sine**4) / x**2
        + beta**2 * x**2 * cosine
    )
    effective = (
        potential_second
        - 2.0 * cosine_two * derivative**2
        - 4.0 * sine * cosine * second
    )

    spacing = float(x[1] - x[0])
    weight_mid = 0.5 * (weight[:-1] + weight[1:])
    potential_mid = 0.5 * (effective[:-1] + effective[1:])

    diagonal = np.zeros(intervals + 1)
    off_diagonal = -weight_mid / spacing + potential_mid * spacing / 6.0
    diagonal[:-1] += weight_mid / spacing + potential_mid * spacing / 3.0
    diagonal[1:] += weight_mid / spacing + potential_mid * spacing / 3.0

    mass_diagonal = np.zeros(intervals + 1)
    mass_off_diagonal = weight_mid * spacing / 6.0
    mass_diagonal[:-1] += weight_mid * spacing / 3.0
    mass_diagonal[1:] += weight_mid * spacing / 3.0

    stiffness = diags(
        [off_diagonal[1:-1], diagonal[1:-1], off_diagonal[1:-1]],
        offsets=[-1, 0, 1],
        format="csr",
    )
    mass = diags(
        [mass_off_diagonal[1:-1], mass_diagonal[1:-1], mass_off_diagonal[1:-1]],
        offsets=[-1, 0, 1],
        format="csr",
    )
    eigenvalues = eigsh(
        stiffness,
        k=count,
        M=mass,
        which="SA",
        tol=2.0e-10,
        maxiter=100000,
        return_eigenvectors=False,
    )
    return np.sort(eigenvalues)


def diquark_hessian_eigenvalues(
    solution: Any,
    alpha: float,
    chi: float,
    box: float,
    intervals: int,
    angular_momentum: int = 0,
    count: int = 4,
) -> np.ndarray:
    """Spectrum for u=x*d: -u''+[alpha^2+chi(1-cos F)+l(l+1)/x^2]u."""

    x = np.linspace(EPS, box, intervals + 1)
    profile = solution.sol(x)[0]
    interior = x[1:-1]
    spacing = float(x[1] - x[0])
    effective = (
        alpha**2
        + chi * (1.0 - np.cos(profile[1:-1]))
        + angular_momentum * (angular_momentum + 1.0) / interior**2
    )
    main = 2.0 / spacing**2 + effective
    off = -np.ones(len(main) - 1) / spacing**2
    operator = diags([off, main, off], offsets=[-1, 0, 1], format="csr")
    eigenvalues = eigsh(
        operator,
        k=count,
        which="SA",
        tol=2.0e-10,
        maxiter=100000,
        return_eigenvectors=False,
    )
    return np.sort(eigenvalues)


def convergence_case(
    model: Model, box: float, intervals: int, solution: Any | None = None
) -> dict[str, Any]:
    if solution is None:
        solution = solve_hedgehog(model.beta, box)
    observables = profile_observables(solution, model.beta, box)
    radial = radial_hessian_eigenvalues(
        solution, model.beta, box, intervals, count=4
    )
    diquark = diquark_hessian_eigenvalues(
        solution, model.alpha, model.chi, box, intervals, count=4
    )
    return {
        "box": box,
        "intervals": intervals,
        "spacing": (box - EPS) / intervals,
        "energy_dimensionless_integral": observables["energy_dimensionless_integral"],
        "B_numeric": observables["B_numeric"],
        "derrick_relative_residual": observables["derrick_relative_residual"],
        "radial_hessian_eigenvalues": radial.tolist(),
        "diquark_hessian_eigenvalues": diquark.tolist(),
        "radial_generalized_omega_sq_in_ef_sq_units": radial.tolist(),
        "diquark_omega_sq_in_ef_sq_units": diquark.tolist(),
    }


def write_outputs(result: dict[str, Any]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    json_path = OUTPUT / "ap_e3_charged_two_colour_proxy.json"
    md_path = OUTPUT / "ap_e3_charged_two_colour_proxy.md"
    json_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checks = "\n".join(
        f"- [{'PASS' if row['pass'] else 'FAIL'}] `{row['group']}` - "
        f"{row['name']}: {row['detail']}"
        for row in result["checks"]
    )
    phase_rows = "\n".join(
        f"| {row['mu_lift_sq']:.2f} | {row['phase']} | {row['sigma']:.6f} | "
        f"{row['diquark_magnitude_real_basis']:.6f} | {min(row['hessian_eigenvalues']):.6e} |"
        for row in result["homogeneous_phase_scan"]
    )
    convergence_rows = "\n".join(
        f"| {row['box']:.0f} | {row['intervals']} | "
        f"{row['energy_dimensionless_integral']:.9f} | {row['B_numeric']:.9f} | "
        f"{row['radial_hessian_eigenvalues'][0]:.9f} | "
        f"{row['diquark_hessian_eigenvalues'][0]:.9f} |"
        for row in result["volume_convergence"]
    )
    coordinate_rows = "\n".join(
        f"| {row['e_g']:.2f} | {row['v_g_over_Lambda_c']:.2f} | "
        f"{row['matched_mu_lift_sq_over_Lambda_c_sq']:.4f} | {row['phase']} | "
        f"{row['meson_condensate_sigma_over_Lambda_c']:.6f} | "
        f"{row['diquark_condensate_real_basis_over_Lambda_c']:.6f} | "
        f"{row['minimum_vacuum_hessian_over_Lambda_c_sq']:.3e} |"
        for row in result["charged_coordinate_scan"]["rows"]
    )
    md_path.write_text(
        f"""# AP-E3 charged two-colour nonlinear-EFT proxy

- Status: `{result['status']}`
- Model class: `{result['model_class']}`
- Checks: `{result['checks_passed']}/{result['checks_total']}`
- Nonlinear EFT proxy pass: `{str(result['nonlinear_eft_proxy_pass']).lower()}`
- Radial Hessian necessary gate: `{str(result['radial_hessian_necessary_gate_pass']).lower()}`
- Full three-dimensional Hessian closed: `{str(result['full_3d_hessian_closed']).lower()}`
- Lattice QC2D closed: `{str(result['lattice_qc2d_closed']).lower()}`
- Physics promotion allowed: `{str(result['physics_promotion_allowed']).lower()}`

This is a deterministic classical low-energy proxy, **not** a lattice or
continuum nonperturbative solution of the microscopic charged two-colour gauge
theory.  The word nonlinear refers to solving the full classical hedgehog BVP,
not to quantum nonperturbative closure.

## Homogeneous charged phase scan

| mu_lift^2 | phase | sigma | |d| | lowest vacuum Hessian |
|---:|---|---:|---:|---:|
{phase_rows}

The analytic transition is `mu_lift^2=-m_pi^2={result['vacuum']['critical_mu_lift_sq']:.6f}`.
At the benchmark, `(m_pi,m_sigma,m_Delta)={result['vacuum']['mass_gaps']}`.

## Declared charged-coordinate chart

The chart uses the explicitly non-derived matching ansatz
`mu_lift^2/Lambda_c^2 = mu_PG,0^2 + c_g (e_g v_g/Lambda_c)^2`.
It is an EFT sensitivity scan, not a microscopic prediction.

| e_g | v_g/Lambda_c | matched mu_lift^2 | phase | sigma | |d| | min Hessian |
|---:|---:|---:|---|---:|---:|---:|
{coordinate_rows}

## B=1 benchmark

- beta: `{result['model']['beta']}`
- alpha: `{result['model']['alpha']}`
- chi: `{result['model']['chi']}`
- origin slope: `{result['soliton']['origin_slope']}`
- numerical B: `{result['soliton']['B_numeric']}`
- B residual: `{result['soliton']['B_residual']}`
- dimensionless energy integral: `{result['soliton']['energy_dimensionless_integral']}`
- Derrick relative residual: `{result['soliton']['derrick_relative_residual']}`
- lowest radial squared frequency `omega^2/(ef)^2`: `{result['hessian']['radial_eigenvalues'][0]}`
- lowest charged-diquark squared frequency `omega^2/(ef)^2`: `{result['hessian']['diquark_l0_eigenvalues'][0]}`
- unstable negative-control squared frequency: `{result['hessian']['diquark_negative_control_eigenvalues'][0]}`

The radial eigenvalue is a finite-box necessary stability test.  It is not a
claim that every nonradial meson, vector, gauge, radial-sigma, or fermion-
determinant channel has been diagonalised.  Positivity in the charged scalar
`l=0` channel does imply positivity for `l>0` within this scalar block because
the centrifugal term is nonnegative.

## Volume convergence

| L | intervals | E | B | lambda_F,min | lambda_Delta,min |
|---:|---:|---:|---:|---:|---:|
{convergence_rows}

## Grid convergence certificate

- Observed radial order: `{result['grid_convergence_certification']['radial_observed_order']}`
- Observed diquark order: `{result['grid_convergence_certification']['diquark_observed_order']}`
- Radial Richardson extrapolate: `{result['grid_convergence_certification']['radial_richardson_extrapolate']}`
- Diquark Richardson extrapolate: `{result['grid_convergence_certification']['diquark_richardson_extrapolate']}`
- Radial fine-grid error estimate: `{result['grid_convergence_certification']['radial_richardson_fine_grid_error_estimate']}`
- Diquark fine-grid error estimate: `{result['grid_convergence_certification']['diquark_richardson_fine_grid_error_estimate']}`

## Mechanical checks

{checks}

## Remaining blockers

"""
        + "\n".join(f"- {item}" for item in result["remaining_blockers"])
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    model = Model()
    critical_sources = [THIS_SCRIPT, TEX, BIB]
    source_manifest = [source_row(path) for path in critical_sources]
    check(
        "P0_provenance",
        "proxy source, independent TeX derivation, and bibliography exist and are hashed",
        all(row["exists"] and row["sha256"] for row in source_manifest),
        f"hashed={sum(bool(row['sha256']) for row in source_manifest)}/{len(source_manifest)}",
    )

    # ------------------------------------------------ homogeneous vacuum and gap scan
    vacuum = np.array([model.f, 0.0, 0.0, 0.0, 0.0, 0.0])
    analytic_hessian = np.diag(
        [
            model.m_sigma_sq,
            model.m_pi**2,
            model.m_pi**2,
            model.m_pi**2,
            model.m_diquark_sq,
            model.m_diquark_sq,
        ]
    )
    numerical_hessian = finite_difference_hessian(vacuum, model, step=5.0e-5)
    hessian_residual = float(np.max(np.abs(numerical_hessian - analytic_hessian)))
    stationarity_residual = abs(
        model.lam * (model.f**2 - model.v_sq) * model.f - model.h
    )
    check(
        "P1_vacuum",
        "the declared mesonic vacuum solves the exact stationarity equation",
        stationarity_residual < 1.0e-14,
        f"residual={stationarity_residual:.3e}",
    )
    check(
        "P1_vacuum",
        "finite-difference and analytic six-field vacuum Hessians agree",
        hessian_residual < 2.0e-6,
        f"max residual={hessian_residual:.3e}",
    )
    check(
        "P1_vacuum",
        "benchmark pion, sigma, and doubly charged diquark gaps are positive",
        model.m_pi > 0.0
        and model.m_sigma_sq > 0.0
        and model.m_diquark_sq > 0.0,
        f"gaps=({model.m_pi:.6f},{math.sqrt(model.m_sigma_sq):.6f},{math.sqrt(model.m_diquark_sq):.6f})",
    )

    mu_scan = [-0.40, -0.30, -0.26, -0.25, -0.10, 0.00, 0.20, 0.56, 1.00]
    phase_scan = [homogeneous_phase(value, model) for value in mu_scan]
    coordinate_scan = charged_coordinate_scan(model)
    phase_labels_correct = all(
        (
            row["phase"] == "diquark_condensed_U1g_broken"
            if row["mu_lift_sq"] < -model.m_pi**2
            else row["phase"]
            in {"critical", "mesonic_U1g_unbroken"}
        )
        for row in phase_scan
    )
    critical_row = next(row for row in phase_scan if row["mu_lift_sq"] == -0.25)
    check(
        "P1_vacuum",
        "analytic scan finds diquark onset exactly at mu_lift^2=-m_pi^2",
        phase_labels_correct
        and abs(min(critical_row["hessian_eigenvalues"])) < 2.0e-6,
        f"critical={-model.m_pi**2:.6f}; critical Hessian min={min(critical_row['hessian_eigenvalues']):.3e}",
    )
    coordinate_phases = {row["phase"] for row in coordinate_scan["rows"]}
    check(
        "P1_vacuum",
        "the declared (e_g,v_g/Lambda_c) chart resolves both sides of the charged phase gate",
        "diquark_condensed_U1g_broken" in coordinate_phases
        and "mesonic_U1g_unbroken" in coordinate_phases,
        f"phases={sorted(coordinate_phases)}; matching={coordinate_scan['matching_status']}",
    )

    # ------------------------------------------------ nonlinear B=1 profile
    baseline_box = 18.0
    baseline_intervals = 1200
    baseline_solution = solve_hedgehog(model.beta, baseline_box)
    soliton = profile_observables(baseline_solution, model.beta, baseline_box)
    check(
        "P2_soliton",
        "nonlinear hedgehog BVP converges with the regular origin series",
        soliton["bvp_max_rms_residual"] < 1.05 * BVP_TOL,
        f"nodes={soliton['bvp_nodes']}; max rms={soliton['bvp_max_rms_residual']:.3e}",
    )
    check(
        "P2_soliton",
        "integrated baryon density agrees with the B=1 boundary formula",
        abs(soliton["B_numeric"] - 1.0) < 2.0e-7,
        f"B_numeric={soliton['B_numeric']:.12f}",
    )
    check(
        "P2_soliton",
        "Derrick virial identity E4=E2+3E0 is numerically satisfied",
        soliton["derrick_relative_residual"] < 2.0e-6,
        f"relative residual={soliton['derrick_relative_residual']:.3e}",
    )
    check(
        "P2_soliton",
        "scale Hessian is positive at the Skyrme saddle",
        soliton["scale_second_derivative"] > 0.0,
        f"d2E/dlambda2={soliton['scale_second_derivative']:.9f}",
    )
    check(
        "P2_soliton",
        "negative control: deleting E4 removes the finite-scale Derrick stationary point",
        soliton["no_skyrme_scale_derivative"] > 0.0,
        f"d(E2*lambda+E0*lambda^3)/dlambda|1={soliton['no_skyrme_scale_derivative']:.9f}",
    )

    # ------------------------------------------------ Hessian blocks
    radial = radial_hessian_eigenvalues(
        baseline_solution,
        model.beta,
        baseline_box,
        baseline_intervals,
        count=5,
    )
    diquark_l0 = diquark_hessian_eigenvalues(
        baseline_solution,
        model.alpha,
        model.chi,
        baseline_box,
        baseline_intervals,
        angular_momentum=0,
        count=5,
    )
    diquark_l1 = diquark_hessian_eigenvalues(
        baseline_solution,
        model.alpha,
        model.chi,
        baseline_box,
        baseline_intervals,
        angular_momentum=1,
        count=3,
    )
    negative_control = diquark_hessian_eigenvalues(
        baseline_solution,
        model.alpha,
        -5.0,
        baseline_box,
        baseline_intervals,
        angular_momentum=0,
        count=3,
    )
    free_box_threshold = model.alpha**2 + (math.pi / baseline_box) ** 2
    check(
        "P3_hessian",
        "all computed spherically symmetric meson Hessian eigenvalues are positive",
        float(radial[0]) > 0.0,
        f"lowest five={radial.tolist()}",
    )
    check(
        "P3_hessian",
        "the doubly charged l=0 soliton mode is positive and bound below its continuum threshold",
        0.0 < float(diquark_l0[0]) < model.alpha**2,
        f"lambda0={diquark_l0[0]:.9f}; continuum={model.alpha**2:.9f}; free-box={free_box_threshold:.9f}",
    )
    check(
        "P3_hessian",
        "centrifugal ordering makes l=0 the most dangerous charged scalar channel",
        float(diquark_l1[0]) > float(diquark_l0[0]),
        f"lambda_l0={diquark_l0[0]:.9f}; lambda_l1={diquark_l1[0]:.9f}",
    )
    check(
        "P3_hessian",
        "negative control: a stronger attractive core coupling produces a charged tachyon",
        float(negative_control[0]) < 0.0,
        f"chi=-5; lambda0={negative_control[0]:.9f}",
    )

    # ------------------------------------------------ discretisation and box convergence
    grid_sequence: list[dict[str, Any]] = []
    for intervals in [600, 1200, 2400]:
        grid_sequence.append(
            convergence_case(
                model,
                baseline_box,
                intervals,
                solution=baseline_solution,
            )
        )
    radial_grid_change = abs(
        grid_sequence[-1]["radial_hessian_eigenvalues"][0]
        - grid_sequence[-2]["radial_hessian_eigenvalues"][0]
    )
    diquark_grid_change = abs(
        grid_sequence[-1]["diquark_hessian_eigenvalues"][0]
        - grid_sequence[-2]["diquark_hessian_eigenvalues"][0]
    )
    radial_previous_change = abs(
        grid_sequence[-2]["radial_hessian_eigenvalues"][0]
        - grid_sequence[-3]["radial_hessian_eigenvalues"][0]
    )
    diquark_previous_change = abs(
        grid_sequence[-2]["diquark_hessian_eigenvalues"][0]
        - grid_sequence[-3]["diquark_hessian_eigenvalues"][0]
    )
    radial_observed_order = math.log2(radial_previous_change / radial_grid_change)
    diquark_observed_order = math.log2(
        diquark_previous_change / diquark_grid_change
    )
    radial_richardson = (
        grid_sequence[-1]["radial_hessian_eigenvalues"][0]
        + (
            grid_sequence[-1]["radial_hessian_eigenvalues"][0]
            - grid_sequence[-2]["radial_hessian_eigenvalues"][0]
        )
        / 3.0
    )
    diquark_richardson = (
        grid_sequence[-1]["diquark_hessian_eigenvalues"][0]
        + (
            grid_sequence[-1]["diquark_hessian_eigenvalues"][0]
            - grid_sequence[-2]["diquark_hessian_eigenvalues"][0]
        )
        / 3.0
    )
    check(
        "P4_convergence",
        "the N=600,1200,2400 sequence shows second-order Hessian convergence and a final change below 2e-5",
        radial_grid_change < 2.0e-5
        and diquark_grid_change < 2.0e-5
        and 1.8 < radial_observed_order < 2.2
        and 1.8 < diquark_observed_order < 2.2,
        f"delta_radial={radial_grid_change:.3e}, p={radial_observed_order:.4f}; "
        f"delta_diquark={diquark_grid_change:.3e}, p={diquark_observed_order:.4f}",
    )

    volume_sequence: list[dict[str, Any]] = []
    for box, intervals in [(14.0, 700), (18.0, 900), (22.0, 1100)]:
        solution = baseline_solution if box == baseline_box else None
        volume_sequence.append(convergence_case(model, box, intervals, solution))
    energy_volume_change = abs(
        volume_sequence[-1]["energy_dimensionless_integral"]
        - volume_sequence[-2]["energy_dimensionless_integral"]
    )
    diquark_volume_change = abs(
        volume_sequence[-1]["diquark_hessian_eigenvalues"][0]
        - volume_sequence[-2]["diquark_hessian_eigenvalues"][0]
    )
    radial_values = [row["radial_hessian_eigenvalues"][0] for row in volume_sequence]
    check(
        "P4_convergence",
        "the final box enlargement stabilises soliton energy and the localized diquark eigenvalue",
        energy_volume_change < 2.0e-6 and diquark_volume_change < 2.0e-5,
        f"delta_E={energy_volume_change:.3e}; delta_diquark={diquark_volume_change:.3e}",
    )
    check(
        "P4_convergence",
        "the lowest radial finite-box state decreases toward the pion continuum beta^2 without crossing zero",
        radial_values[0] > radial_values[1] > radial_values[2] > model.beta**2,
        f"sequence={radial_values}; beta^2={model.beta**2:.9f}",
    )

    all_pass = all(row["pass"] for row in CHECKS)
    nonlinear_eft_proxy_pass = all_pass
    radial_gate = bool(
        soliton["B_numeric"] > 0.999
        and radial[0] > 0.0
        and diquark_l0[0] > 0.0
        and negative_control[0] < 0.0
    )
    result: dict[str, Any] = {
        "schema_version": 1,
        "status": "pass" if all_pass else "fail",
        "all_pass": all_pass,
        "checks_passed": sum(row["pass"] for row in CHECKS),
        "checks_total": len(CHECKS),
        "model_class": MODEL_NAME,
        "nonlinear_eft_proxy_pass": nonlinear_eft_proxy_pass,
        "radial_hessian_necessary_gate_pass": radial_gate,
        "full_3d_hessian_closed": False,
        "quantum_nonperturbative_gauge_theory_closed": False,
        "lattice_qc2d_closed": False,
        "global_unwinding_excluded": False,
        "physics_promotion_allowed": False,
        "model": {
            **asdict(model),
            "h": model.h,
            "v_sq": model.v_sq,
            "m_sigma_sq": model.m_sigma_sq,
            "m_diquark_sq": model.m_diquark_sq,
            "beta": model.beta,
            "alpha": model.alpha,
        },
        "vacuum": {
            "critical_mu_lift_sq": -model.m_pi**2,
            "stationarity_residual": stationarity_residual,
            "analytic_hessian": analytic_hessian.tolist(),
            "finite_difference_hessian": numerical_hessian.tolist(),
            "hessian_max_residual": hessian_residual,
            "mass_gaps": [
                model.m_pi,
                math.sqrt(model.m_sigma_sq),
                math.sqrt(model.m_diquark_sq),
            ],
        },
        "homogeneous_phase_scan": phase_scan,
        "charged_coordinate_scan": coordinate_scan,
        "soliton": soliton,
        "hessian": {
            "operator_metadata": {
                "eigenvalue_interpretation": "dimensionless squared frequency omega^2/(e f)^2",
                "radial_problem": "generalized Sturm-Liouville H_F eta = omega^2 A eta with Dirichlet finite-box boundaries",
                "diquark_problem": "ordinary Sturm-Liouville H_Delta,l u = omega^2 u with Dirichlet finite-box boundaries",
                "legacy_eigenvalue_keys_retained_for_compatibility": True,
            },
            "radial_eigenvalues": radial.tolist(),
            "radial_generalized_omega_sq_in_ef_sq_units": radial.tolist(),
            "radial_finite_box_continuum_floor": model.beta**2,
            "diquark_l0_eigenvalues": diquark_l0.tolist(),
            "diquark_l1_eigenvalues": diquark_l1.tolist(),
            "diquark_l0_omega_sq_in_ef_sq_units": diquark_l0.tolist(),
            "diquark_l1_omega_sq_in_ef_sq_units": diquark_l1.tolist(),
            "diquark_continuum_threshold": model.alpha**2,
            "diquark_free_box_threshold": free_box_threshold,
            "diquark_negative_control_chi": -5.0,
            "diquark_negative_control_eigenvalues": negative_control.tolist(),
        },
        "grid_convergence": grid_sequence,
        "grid_convergence_certification": {
            "expected_order": 2.0,
            "radial_observed_order": radial_observed_order,
            "diquark_observed_order": diquark_observed_order,
            "radial_final_doubling_change": radial_grid_change,
            "diquark_final_doubling_change": diquark_grid_change,
            "radial_richardson_extrapolate": radial_richardson,
            "diquark_richardson_extrapolate": diquark_richardson,
            "radial_richardson_fine_grid_error_estimate": abs(
                radial_richardson
                - grid_sequence[-1]["radial_hessian_eigenvalues"][0]
            ),
            "diquark_richardson_fine_grid_error_estimate": abs(
                diquark_richardson
                - grid_sequence[-1]["diquark_hessian_eigenvalues"][0]
            ),
        },
        "volume_convergence": volume_sequence,
        "checks": CHECKS,
        "source_manifest": source_manifest,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "deterministic": True,
            "randomness_used": False,
        },
        "remaining_blockers": [
            "Replace the classical EFT proxy by a controlled lattice or other nonperturbative calculation of the charged SU(2)c gauge theory.",
            "Compute the full three-dimensional coupled meson, radial-sigma, vector, gauge, and ghost Hessian, including collective zero modes.",
            "Include the fermion determinant and quantum corrections; the present Hessian is tree-level classical.",
            "Prove that the microscopic charged theory dynamically selects the mesonic vacuum and match mu_lift_sq and chi from UV observables.",
            "Exclude global escape/unwinding paths in the faithful full target; positivity of local Hessian blocks is only a necessary condition.",
            "No Route-E family portal follows from this proxy scan.",
        ],
    }
    write_outputs(result)
    print(
        f"AP-E3 charged two-colour proxy: {result['status']} "
        f"({result['checks_passed']}/{result['checks_total']}); "
        f"radial gate={radial_gate}; promotion={result['physics_promotion_allowed']}"
    )


if __name__ == "__main__":
    main()
