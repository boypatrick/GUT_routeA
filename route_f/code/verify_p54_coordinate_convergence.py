#!/usr/bin/env python3
"""Frozen-action coordinate/preconditioning and spectral-series experiments.

No coupling scan, new Hessian, new physical pole, or loop-convergence claim.
The conformal-series bound is a theorem for the equal-mass one-loop master,
not for the unknown all-orders theory. Existing results are read, never edited.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import eigh

RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_coordinate_convergence"
mp.mp.dps = 60


def inverse_sqrt(a):
    values, vectors = np.linalg.eigh(a)
    if values.min() <= 0:
        raise ValueError("Whitening requires a positive matrix")
    return (vectors * values**(-.5)) @ vectors.T


def conformal_coordinate(z):
    z = mp.mpf(z)
    root = mp.sqrt(1+z/4)
    return z/(4*(root+1)**2)  # cancellation-free (root-1)/(root+1)


def bubble_exact(z):
    z = mp.mpf(z)
    if z == 0:
        return mp.mpf(0)
    return 2*(mp.sqrt((4+z)/z)*mp.asinh(mp.sqrt(z)/2)-1)


def bubble_taylor(z, n):
    z = mp.mpf(z)
    return mp.fsum((-1)**(k+1)*mp.factorial(k)**2*z**k /
                   (k*mp.factorial(2*k+1)) for k in range(1, n+1))


def bubble_conformal(z, n):
    w = conformal_coordinate(z)
    return mp.fsum(mp.mpf(8)*k*w**k/(4*k*k-1) for k in range(1, n+1))


def conformal_tail_bound(z, n):
    w = conformal_coordinate(z)
    k = n+1
    return mp.mpf(8)*k/(4*k*k-1)*w**k/(1-w)


def richardson(a, target, steps):
    values = np.linalg.eigvalsh(a)
    rate = 2/(values[0]+values[-1])
    q = (values[-1]-values[0])/(values[-1]+values[0])
    rhs = a@target
    x = np.zeros_like(target)
    energy0 = target@a@target
    rows = []
    for n in range(max(steps)+1):
        if n in steps:
            error = x-target
            rows.append(dict(iterations=n, relative_energy_error=float(
                np.sqrt(max(0., error@a@error)/energy0))))
        x += rate*(rhs-a@x)
    return q, rows


def run():
    gauge_path = RF/"output/p54_gauge_mixed_momentum.json"
    ren_path = RF/"output/p54_common_renormalization.json"
    gauge = json.loads(gauge_path.read_text())
    ren = json.loads(ren_path.read_text())
    if not gauge["summary"]["all_pass"]:
        raise ValueError("The gauge-kernel input has failing checks")
    for name, digest in gauge["source_sha256"].items():
        if hashlib.sha256((RF/name).read_bytes()).hexdigest() != digest:
            raise ValueError("Stale gauge-kernel input: "+name)
    g = np.asarray(gauge["radial_metric"])
    dz = np.asarray(gauge["hard_scalar_plus_gauge_kinetic"]["matrix"])
    kinetic = g+dz
    radii = np.asarray(gauge["vacuum"])
    base_eigs = eigh(dz, g, eigvals_only=True)
    checks = []

    def check(name, actual, expected=0., tol=2e-11):
        actual, expected = np.asarray(actual), np.asarray(expected)
        residual = float(np.linalg.norm(actual-expected)/max(1., np.linalg.norm(expected)))
        checks.append(dict(name=name, residual=residual, tolerance=tol, passed=bool(residual<tol)))

    transformations = [
        ("identity", np.eye(3)),
        ("u=x/pi; x=pi*u", np.pi*np.eye(3)),
        ("u=pi*x; x=u/pi", np.eye(3)/np.pi),
        ("relative-VEV or local log chart", np.diag(radii)),
        ("tree-metric whitening", inverse_sqrt(g)),
        ("known-local-kinetic whitening; not new loop counting", inverse_sqrt(kinetic)),
    ]
    cases = []
    for name, j in transformations:
        gp, zp, kp = j.T@g@j, j.T@dz@j, j.T@kinetic@j
        eigs = eigh(zp, gp, eigvals_only=True)
        check(name+": invariant relative-loop spectrum", eigs, base_eigs)
        eigen_rows = []
        for row in gauge["momenta"]:
            h = np.asarray(row["complete_radial_bosonic_kernel"])
            transformed = j.T@h@j
            got = eigh(transformed, gp, eigvals_only=True)
            check(name+f": invariant kernel spectrum at {row['pE2']}",
                  got, row["generalized_bosonic_eigenvalues"])
            # An actual source contraction is invariant when the source
            # is transformed too; this is stronger than comparing entries.
            source = np.array([1., -.3, .2])
            source_p = j.T@source
            response = source_p@np.linalg.solve(transformed, source_p)
            check(name+f": invariant source response at {row['pE2']}",
                  response, source@np.linalg.solve(h, source))
            eigen_rows.append(dict(pE2=row["pE2"], eigenvalues=got.tolist()))
        target = np.linalg.solve(j, np.array([1., -.3, .2]))
        q, iterations = richardson(kp, target, (0, 1, 2, 4, 8, 16, 32))
        for row in iterations:
            check(name+f": Richardson energy bound at {row['iterations']}",
                  max(0., row["relative_energy_error"]-q**row["iterations"]), tol=2e-12)
        cases.append(dict(name=name, field_jacobian=j.tolist(),
            transformed_metric=gp.tolist(), transformed_loop=zp.tolist(),
            raw_loop_operator_norm=float(np.linalg.norm(zp, 2)),
            metric_condition_number=float(np.linalg.cond(gp)),
            known_kinetic_condition_number=float(np.linalg.cond(kp)),
            generalized_loop_eigenvalues=eigs.tolist(),
            rho_loop=float(max(abs(eigs))), kernel_spectra=eigen_rows,
            Richardson_contraction_bound=float(q), iterations=iterations))
    check("tree whitening is identity", transformations[-2][1].T@g@transformations[-2][1], np.eye(3))
    check("known kinetic whitening is identity",
          transformations[-1][1].T@kinetic@transformations[-1][1], np.eye(3))

    # Nonlinear changes need the gradient term away from a stationary point.
    t = np.asarray(ren["broken_tadpoles"]["total"])
    loop_h = np.asarray(ren["radial_bosonic"]["scalar_loop_Hessian"])+np.asarray(
        ren["radial_bosonic"]["vector_loop_Hessian"])
    j = np.diag(radii)
    gradient_q = j.T@t
    ordinary_log_hessian = j.T@loop_h@j+np.diag(gradient_q)
    covariant_log_hessian = ordinary_log_hessian-np.diag(gradient_q)
    check("log-chart covariant Hessian congruence", covariant_log_hessian, j.T@loop_h@j)
    ct = np.asarray(ren["radial_bosonic"]["CT_Hessian"])
    check("common fixed-VEV tadpole cancellation", t+ct@radii)
    nonlinear = dict(chart="r_i=r0_i*exp(q_i); Gamma^i_ii=1",
        loop_gradient_q=gradient_q.tolist(), ordinary_loop_Hessian=ordinary_log_hessian.tolist(),
        covariant_loop_Hessian=covariant_log_hessian.tolist(),
        warning="Dropping the gradient/connection term off shell changes the calculation.")

    # Absorbing a known local one-loop insertion is not evidence about
    # the independent, uncomputed two-loop self-energy or vertices.
    normalized = inverse_sqrt(g)@dz@inverse_sqrt(g)
    exact_inverse = np.linalg.inv(np.eye(3)+normalized)
    neumann = []
    for n in (1, 2, 4, 8, 16):
        partial = sum(np.linalg.matrix_power(-normalized, k) for k in range(n+1))
        neumann.append(dict(order=n, matrix_inverse_error=float(np.linalg.norm(partial-exact_inverse, 2))))
    reorganization = dict(known_quadratic_inverse=exact_inverse.tolist(),
        original_Neumann_radius=float(max(abs(base_eigs))),
        relative_to_resummed_metric=(base_eigs/(1+base_eigs)).tolist(),
        original_Neumann_errors=neumann, independent_two_loop_remainder_known=False,
        physical_convergence_proven=False)
    check("original local inverse Neumann series does not converge",
          float(neumann[-1]["matrix_inverse_error"]<=neumann[0]["matrix_inverse_error"]))

    # Equal-mass bubbles are actual building blocks, not new parameter cards.
    small_vector = gauge["vector_group_masses2"][1]
    scalar50 = 12*radii[1]**2
    probes = [("lightest massive vector", small_vector, row["pE2"])
              for row in gauge["momenta"]]
    probes += [("certified P50 at unchanged tau=1", scalar50, .05),
               ("certified P50 at unchanged tau=1", scalar50, 1.)]
    spectral = []
    for name, mass2, s in probes:
        z = mp.mpf(s)/mp.mpf(mass2)
        exact = bubble_exact(z)
        quadrature = mp.quad(lambda x: mp.log1p(z*x*(1-x)), [0, .5, 1])
        check(f"equal-mass master vs independent integral {name} {s}",
              float(exact), float(quadrature), tol=1e-13)
        w = conformal_coordinate(z)
        check(f"conformal inverse map {name} {s}", float(16*w/(1-w)**2), float(z))
        rows = []
        for n in (2, 4, 8, 16, 32):
            tay, mapped = bubble_taylor(z, n), bubble_conformal(z, n)
            tail = exact-mapped
            bound = conformal_tail_bound(z, n)
            if not (-mp.mpf("1e-50") <= tail <= bound+mp.mpf("1e-50")):
                raise AssertionError("Analytic positive-tail bound violated")
            check(f"certified conformal tail {name} {s} N={n}",
                  float(max(mp.mpf(0), tail-bound)), tol=1e-13)
            rows.append(dict(terms=n, Taylor_value=float(tay), Taylor_abs_error=float(abs(tay-exact)),
                conformal_value=float(mapped), conformal_abs_error=float(abs(tail)),
                conformal_absolute_tail_bound=float(bound)))
        spectral.append(dict(name=name, mass2=float(mass2), pE2=s, z=float(z),
            pair_threshold4m2=4*float(mass2), derivative_expansion_ratio=float(z/4),
            inside_open_Taylor_disk=bool(z<4), conformal_w=float(w),
            exact_master=float(exact), approximations=rows))
    pilot = next(row for row in spectral if row["name"]=="lightest massive vector" and row["pE2"]==.05)
    check("pilot conformal N=8 is more accurate than Taylor N=8",
          float(pilot["approximations"][2]["conformal_abs_error"]>=pilot["approximations"][2]["Taylor_abs_error"]))
    check("pilot Taylor error grows outside its disk",
          float(pilot["approximations"][3]["Taylor_abs_error"]<=pilot["approximations"][2]["Taylor_abs_error"]))

    sources = [Path(__file__), gauge_path, ren_path]
    return dict(schema="p54-coordinate-and-spectral-convergence-v1", date="2026-09-17",
        scope="frozen action; coordinate invariants, linear solver and equal-mass one-loop series only",
        radial_metric=g.tolist(), local_hard_scalar_plus_gauge_insertion=dz.tolist(),
        baseline_generalized_eigenvalues=base_eigs.tolist(), coordinate_cases=cases,
        nonlinear_chart=nonlinear, quadratic_reorganization=reorganization,
        equal_mass_spectral_probes=spectral, pilot=pilot,
        no_new_action_or_coupling_card=True, new_Hessian_evaluations=0,
        default_parameters_changed=False, physical_fit_enabled=False,
        all_orders_loop_convergence_proven=False, full_kernel_conformal_replacement_enabled=False,
        scalar_soft_log_kept_nonlocal=True,
        source_sha256={str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        checks=checks, summary=dict(passed=sum(c["passed"] for c in checks), total=len(checks),
                                   all_pass=all(c["passed"] for c in checks)))


def markdown(r):
    lines = ["# Coordinate changes and convergence: frozen P54 experiment", "",
        f"Checks: {r['summary']['passed']}/{r['summary']['total']}. No physical convergence claim.", "",
        "| Chart | condition of known local kinetic matrix | invariant rho(loop/tree) |",
        "|---|---:|---:|"]
    for row in r["coordinate_cases"]:
        lines.append(f"| {row['name']} | {row['known_kinetic_condition_number']:.8g} | {row['rho_loop']:.8g} |")
    p = r["pilot"]
    lines += ["", "## Threshold-aware conformal pilot", "",
        f"Actual lightest-vector mass square={p['mass2']:.10g}; pE2={p['pE2']}; z={p['z']:.9g}; w={p['conformal_w']:.9g}.",
        "The Taylor disk ends at |z|=4. The mapped series has a proven absolute tail bound for this master.",
        "", "| Terms | Taylor absolute error | Conformal absolute error | Proven tail bound |",
        "|---:|---:|---:|---:|"]
    for row in p["approximations"]:
        lines.append(f"| {row['terms']} | {row['Taylor_abs_error']:.7g} | {row['conformal_abs_error']:.7g} | {row['conformal_absolute_tail_bound']:.7g} |")
    lines += ["", "## Interpretation", "",
        "- Pi unit changes cannot reduce the dimensionless loop/tree spectrum or change physical source responses.",
        "- Metric whitening improves linear algebra; even exact normalization of the known local kinetic matrix does not establish unknown higher-loop control.",
        "- The threshold-derived conformal variable genuinely enlarges the convergence domain of the equal-mass momentum series, without fitting coefficients.",
        "- Existing P-MOM1 already evaluates exact masters; this is a validated series representation, not a repair of its values.",
        "- Unequal-mass, mixed tensor and massless-log pieces need a coordinated kernel-level extension before using this as a replacement.",
        "- Retaining near-threshold modes/nonlocal kernels remains the structural option; Nielsen, full physical mixing, matching and finite seesaw gates stay open."]
    return "\n".join(lines)+"\n"


if __name__ == "__main__":
    result = run()
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    OUT.with_suffix(".md").write_text(markdown(result))
    print(markdown(result))
    failed = [c for c in result["checks"] if not c["passed"]]
    print("Failed checks:", failed)
    if failed:
        raise SystemExit(1)
