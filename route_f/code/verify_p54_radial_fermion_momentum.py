#!/usr/bin/env python3
"""Exact radial Majorana momentum subset and actual gauge vertex inventory.

Three unfitted Majorana Takagi masses remain inputs.  This is not a complete
gauged pole kernel, a fit, or a physical-stability decision.  No action or
benchmark parameters are modified.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh

import verify_p54_common_renormalization as common
from verify_p54_goldstone_ir import log_bubble

RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_radial_fermion_momentum"
LOOP = 16 * np.pi**2


def fixed_vev_majorana_kernel(masses, sigma, mu, s):
    """Sigma-coordinate Hessian at Euclidean p^2=s, including common CT.

    M_i=sigma*f_i.  Majorana multiplicity is two, not four.  The common
    fixed-VEV CT subtracts t_sigma/sigma at every momentum.
    """
    masses = np.asarray(masses, dtype=float)
    if np.any(masses < 0) or sigma <= 0 or mu <= 0 or s < 0:
        raise ValueError("Nonnegative masses/s and positive sigma/mu required")
    return float(-sum((m / sigma)**2 * (4*m*m+s)
                     * log_bubble(m*m, m*m, s, mu*mu)
                     for m in masses if m > 0) / LOOP)


def majorana_kinetic(masses, sigma, mu):
    """Finite MSbar derivative d Pi_F/d p_E^2 at zero; positive masses."""
    masses = np.asarray(masses, dtype=float)
    if np.any(masses < 0) or sigma <= 0 or mu <= 0:
        raise ValueError("Nonnegative masses and positive sigma/mu required")
    return float(-sum((m / sigma)**2 * (np.log(m*m/(mu*mu))+2/3)
                     for m in masses if m > 0) / LOOP)


def independent_spin_trace(masses, sigma, mu, s):
    """Independent Feynman-parameter numerator integration plus tadpole CT.

    1/2 Tr[S Y S Y] = 2 y^2 int (m^2-k.(k+p))/(A B).
    After shifting k, D=m^2+x(1-x)s.  The MS finite numerator integral is
    2 y^2 int D (1-3 log(D/mu^2))/(16 pi^2).
    """
    value = 0.
    for m in masses:
        if m == 0:
            continue
        y2 = (m/sigma)**2
        f = lambda x: (m*m+x*(1-x)*s) * (1-3*np.log((m*m+x*(1-x)*s)/(mu*mu)))
        raw = 2*y2*quad(f, 0, 1, epsabs=1e-12, epsrel=1e-12)[0]/LOOP
        tad_over_sigma = -2*y2*m*m*(np.log(m*m/(mu*mu))-1)/LOOP
        value += raw-tad_over_sigma
    return float(value)


def run():
    full_path = RF / "output/p54_full_doublet_cw.json"
    ren_path = RF / "output/p54_common_renormalization.json"
    full = json.loads(full_path.read_text())
    ren = json.loads(ren_path.read_text())
    r = np.array([full["vacuum"][k] for k in ("omega", "sigma", "vs")])
    mu = ren["mu"]
    p1 = common.yuk.module("radial_fermion_p1", RF / "code/verify_p54_p1_hessian_spectrum.py")
    cw = common.yuk.module("radial_fermion_cw", RF / "code/verify_p54_p2_bosonic_cw.py")
    checks = []

    def check(name, a, b=0., tol=1e-8):
        a, b = np.asarray(a), np.asarray(b)
        error = float(np.linalg.norm(a-b)/max(1., np.linalg.norm(a), np.linalg.norm(b)))
        checks.append(dict(name=name, residual=error, tolerance=tol, passed=error < tol))

    cap = 6*mu**4/(np.e*LOOP*r[1]**2)
    mass_cases = [np.array(x) for x in ([.02,.06,.11], [.25,.35,.45], [.6,.8,1.], [0.,.03,.2])]
    mass_cases.append(np.full(3, mu*np.exp(-.25)))
    rows = []
    for k, masses in enumerate(mass_cases):
        pi0 = fixed_vev_majorana_kernel(masses, r[1], mu, 0.)
        check(f"zero_momentum_matches_existing_CT_function_{k}", pi0,
              common.majorana_fixed_VEV_radial(masses, r[1], mu))
        check(f"existing_all_Majorana_bound_{k}", max(0., pi0-cap))
        kinetic = majorana_kinetic(masses, r[1], mu)
        step = min(masses[masses > 0]**2)*1e-4
        fd = (-3*pi0+4*fixed_vev_majorana_kernel(masses,r[1],mu,step)
              -fixed_vev_majorana_kernel(masses,r[1],mu,2*step))/(2*step)
        check(f"independent_kinetic_derivative_{k}", fd, kinetic, tol=2e-8)
        values = []
        for s in (0., .001, .01, .1, 1.):
            value = fixed_vev_majorana_kernel(masses,r[1],mu,s)
            check(f"independent_spin_trace_integral_{k}_{s}", value,
                  independent_spin_trace(masses,r[1],mu,s),tol=3e-10)
            values.append(dict(pE2=s, sigma_coordinate_kernel=value))
        rows.append(dict(synthetic_test_masses=masses.tolist(), finite_sigma_kinetic=kinetic,
                         generalized_radial_kinetic=kinetic/2, values=values))
    check("three_equal_Majorana_masses_saturate_zero_momentum_cap", rows[-1]["values"][0]["sigma_coordinate_kernel"],cap)
    check("all_zero_Yukawa_limit", fixed_vev_majorana_kernel([0.,0.,0.],r[1],mu,.1))

    # Actual canonical real representation tensors; no full scalar Hessian.
    radial = np.column_stack([p1.vacuum_vector(*e) for e in np.eye(3)])
    x = p1.vacuum_vector(*r)
    gens = cw.generators()
    v = cw.field_orbit(p1,x,gens)
    d = np.array([cw.field_orbit(p1,q,gens) for q in radial.T])
    g = full["scheme"]["mu_over_omega"] / r[0]
    m2 = g*g*v.T@v
    b = np.array([g*g*(dr.T@v+v.T@dr) for dr in d])
    c = np.array([[g*g*(di.T@dj+dj.T@di) for dj in d] for di in d])
    check("gauge_orbit_linear_reconstruction",np.einsum("r,ria->ia",r,d),v)
    check("gauge_mass_Euler_first_derivative",np.einsum("r,rab->ab",r,b),2*m2)
    check("gauge_mass_Euler_second_derivative",np.einsum("r,s,rsab->ab",r,r,c),2*m2)
    check("gauge_mass_matrix_symmetric",m2,m2.T)
    check("33_broken_gauge_generators",np.linalg.matrix_rank(m2,tol=1e-9),33)
    check("singlet_radial_gauge_vertex_vanishes",d[2])
    metric = radial.T@radial
    check("radial_metric",metric,np.diag([12/5,2,1]))
    for i in range(3):
        step = .031
        vp = cw.field_orbit(p1,x+step*radial[:,i],gens)
        vm = cw.field_orbit(p1,x-step*radial[:,i],gens)
        check(f"independent_gauge_first_derivative_{i}",g*g*(vp.T@vp-vm.T@vm)/(2*step),b[i])
        check(f"independent_gauge_second_derivative_{i}",g*g*(vp.T@vp+vm.T@vm-2*v.T@v)/(step*step),c[i,i])
    mixed_gram = g*g*np.einsum("ria,sia->rs",d,d)
    check("mixed_scalar_vector_vertex_Gram_positive",max(0.,-eigh(mixed_gram,metric,eigvals_only=True).min()))
    source_paths = [Path(__file__), Path(common.__file__), Path(p1.__file__), Path(cw.__file__),
                    RF/"code/verify_p54_goldstone_ir.py", full_path, ren_path]
    return dict(schema="p54-radial-fermion-momentum-v1", date="2026-09-13",
                scope="Exact one-loop radial Majorana subset, not complete gauged momentum kernel",
                vacuum=r.tolist(), radial_metric=metric.tolist(), mu=mu, g=g,
                zero_momentum_three_Majorana_cap=cap, synthetic_cases=rows,
                gauge_vertex_inventory=dict(vector_count=45, broken_count=33,
                    gauge_mass2=m2.tolist(), radial_first_mass2=b.tolist(), radial_second_mass2=c.tolist(),
                    derivative_scalar_vector_Gram=mixed_gram.tolist(),
                    derivative_scalar_vector_generalized_Gram=eigh(mixed_gram,metric,eigvals_only=True).tolist()),
                gauge_momentum_integrals_evaluated=False, ghost_xi_limit_numerically_verified=False,
                full_field_fermion_kernel_evaluated=False, full_stability_decided=False,
                default_parameters_changed=False, physical_fit_promoted=False,
                checks=checks, passed=sum(z["passed"] for z in checks), total=len(checks),
                source_hashes={str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths})


if __name__ == "__main__":
    report = run()
    OUT.with_suffix(".json").write_text(json.dumps(report,indent=2)+"\n")
    text = ["# P54 radial Majorana momentum subset", "",
            f"Checks: **{report['passed']}/{report['total']}**.", "",
            "The three Majorana masses are input parameters, not fitted values. All five numerical cards below are synthetic regression tests. No physical stability conclusion follows.", "",
            "Exact fixed-VEV radial kernel: `Pi_F(s)=-sum[(M_i/sigma)^2 (4 M_i^2+s) L(M_i^2,M_i^2;s)]/(16 pi^2)`.", "",
            "Only the sigma-sigma coordinate is nonzero. Its canonical generalized component has an extra factor 1/2 because G_sigma,sigma=2.", "",
            f"Zero-momentum universal three-Majorana cap: `{report['zero_momentum_three_Majorana_cap']:.12g}`.", "",
            "The actual gauge mass first/second vertices and mixed scalar-vector derivative Gram are stored in JSON. Their tensors are fixed by the existing canonical representation and coupling; they do not require a flavor fit.", "",
            "Still missing: finite vector/scalar-vector momentum integrals with dimensional rational terms, explicit ghost/gauge-fixing limit checks, common wave-function counterterms, Nielsen tests, complete field-space kernel and physical-pole analysis.", "",
            "```json",json.dumps(report["synthetic_cases"],indent=2),"```", ""]
    OUT.with_suffix(".md").write_text("\n".join(text))
    print(json.dumps({k:report[k] for k in ("passed","total","zero_momentum_three_Majorana_cap")}))
    if report["passed"] != report["total"]:
        raise SystemExit(1)
