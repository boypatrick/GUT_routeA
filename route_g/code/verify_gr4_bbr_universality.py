#!/usr/bin/env python3
"""Physical-source control: thermal E1 Stark shifts of Yb+ E2/E3 clocks.

Source-pinned coefficients, no fits, no new experimental observations.
Static leading-order prediction and explicitly conditional stress box.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from itertools import product
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT/"data"/"BBR_CLOCK_INPUTS.json"
OUT = ROOT/"output"/"gr4_bbr_universality"


def planck_moment(power):
    def integrand(x):
        if x == 0 or x > 700:
            return 0.
        return x**power/np.expm1(x)
    return quad(integrand, 0., np.inf, epsabs=2e-11, epsrel=2e-13)[0]


def ratio_offset(beta, u):
    """R(u)/R(0)-1, stable without subtracting nearly equal unit numbers."""
    b2, b3 = beta
    return (b3-b2)*u/(1+b2*u)


def ratio_modulation(beta, low, high):
    """R(high)/R(low)-1, exact for the retained static frequency model."""
    b2, b3 = beta
    return (b3-b2)*(high-low)/((1+b2*high)*(1+b3*low))


def run():
    cfg = json.loads(INPUT.read_text())
    constants = cfg["constants"]
    h, kb, c, eps0 = [constants[k] for k in ["h_J_s", "kB_J_K", "c_m_s", "epsilon0_F_m"]]
    nu = np.array([x["frequency_Hz"] for x in cfg["transitions"]])
    alpha = np.array([x["delta_alpha_SI"] for x in cfg["transitions"]])
    sigma = np.array([x["published_sigma_alpha_SI"] for x in cfg["transitions"]])
    old = sorted(p for p in (ROOT/"code").glob("verify_*.py") if p.name != Path(__file__).name)
    old += [ROOT/"output"/f"gr{i}_{name}.json" for i, name in
            [(1,"relational_clock"),(2,"two_event_records"),(3,"density_clocks")]]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in old}
    checks = []

    def check(name, ok, value=None):
        row = dict(name=name, passed=bool(ok))
        if value is not None:
            row["value"] = float(value)
        checks.append(row)

    def relative(name, actual, expected, tol=2e-11):
        actual, expected = np.asarray(actual), np.asarray(expected)
        error = float(np.max(abs(actual-expected))/max(float(np.max(abs(expected))), 1e-100))
        check(name, error < tol, error)

    moment3 = planck_moment(3)
    moment5 = planck_moment(5)
    relative("Planck:third_moment", moment3, np.pi**4/15)
    relative("Planck:fifth_moment", moment5, 8*np.pi**6/63)
    arad = 8*np.pi**5*kb**4/(15*h**3*c**3)
    integral_arad = 8*np.pi*kb**4*moment3/(h**3*c**3)
    relative("Planck:energy_density_coefficient", arad, integral_arad)
    u300 = arad*300**4
    rms300 = np.sqrt(u300/eps0)
    relative("Planck:published_rounded_rms", rms300, 831.9, tol=6e-5)
    relative("Planck:E_and_B_equal_energy", eps0*rms300**2/2, u300/2)
    relative("convention:E2_sign_reversal", alpha[0], -(-6.9e-40))
    check("source:both_measured_delta_alpha_positive", np.all(alpha > 0))
    check("source:both_positive_in_three_sigma_box", np.all(alpha-3*sigma > 0))
    # No copy of any claimed "latest" or combined covariance is introduced.
    beta = -alpha/(2*h*eps0*nu)
    beta_sigma = sigma/(2*h*eps0*nu)
    check("coupling:both_redshift", np.all(beta < 0))
    coefficient_ratio = beta[0]/beta[1]
    check("universality:fractional_coefficients_differ", coefficient_ratio > 7)
    relative("universality:ratio_vs_atomic_inputs", coefficient_ratio, alpha[0]*nu[1]/(alpha[1]*nu[0]))

    # Three-level diagonal lemma: delta H is f H0 + c I iff both
    # ground-to-excited gap shifts have the same fractional response.
    bare = np.array([0., 1., nu[1]/nu[0]])
    perturb = np.array([0., 1., alpha[1]/alpha[0]])
    check("operator:response_not_proportional",
          abs(perturb[2]-bare[2]) > .7)
    common = .07*bare+.43
    relative("operator:universal_gap_condition",
             (common[1:]-common[0])/(bare[1:]-bare[0]), [.07,.07])
    # Energy-origin shifts cannot change transition polarizabilities.
    absolute_alpha = np.array([3.,3.+alpha[0]/1e-40,3.+alpha[1]/1e-40])
    relative("operator:common_level_shift_cancels",
             np.diff(absolute_alpha[[0,1]])[0], alpha[0]/1e-40)

    rows = []
    for T in cfg["calculation"]["temperatures_K"]:
        u = arad*T**4
        shifts = -alpha*u/(2*h*eps0)
        fractions = shifts/nu
        relative(f"T={T}:energy_shift_route", fractions, beta*u)
        check(f"T={T}:positive_frequencies", np.all(1+fractions > 0))
        check(f"T={T}:ratio_positive", ratio_offset(beta,u) > 0)
        for i in range(2):
            # Integrate the thermal spectrum independently of the closed T^4 form.
            integral = 8*np.pi*(kb*T)**4/(h**3*c**3)*moment3
            numeric = -alpha[i]*integral/(2*h*eps0*nu[i])
            relative(f"T={T}:transition{i}:spectral_integral", numeric, fractions[i])
        rows.append(dict(temperature_K=T, thermal_density_J_m3=u, rms_field_V_m=np.sqrt(u/eps0),
                         shifts_Hz=shifts.tolist(), fractional_shifts=fractions.tolist(),
                         ratio_offset=ratio_offset(beta,u)))
    relative("static:fourth_power_scaling", np.array(rows[1]["shifts_Hz"])/rows[0]["shifts_Hz"],
             np.ones(2)*(310/300)**4)
    modulation = ratio_modulation(beta, rows[0]["thermal_density_J_m3"], rows[1]["thermal_density_J_m3"])
    # High-precision check is for stable evaluation, not fitting small residuals.
    with localcontext() as ctx:
        ctx.prec = 60
        b2,b3 = map(lambda x: Decimal(str(x)), beta)
        ul,uh = map(lambda x: Decimal(str(x)), [r["thermal_density_J_m3"] for r in rows])
        direct = ((1+b3*uh)/(1+b2*uh))/((1+b3*ul)/(1+b2*ul))-1
        relative("ratio:stable_formula", modulation, float(direct))
        for factor in [Decimal("0.9"),Decimal("1.1")]:
            original = (1+b3*uh)/(1+b2*uh)
            transformed = (factor*(1+b3*uh))/(factor*(1+b2*uh))
            check("ratio:common_factor_cancels:"+str(factor), original == transformed)
    # Frequency ratios cancel any shared accumulated physical time.
    for elapsed in [.2,1.,7.]:
        theta = 2*np.pi*nu*elapsed
        relative(f"relational:phase_ratio:{elapsed}", theta[1]/theta[0], nu[1]/nu[0])

    # Independent deterministic corners, not Gaussian sampling or a data fit.
    stress = cfg["calculation"]["stress"]
    ns, eta = stress["alpha_radius_in_published_sigma"], stress["independent_relative_dynamic_envelope"]
    worst = []
    for j, signs in enumerate(product([-1,1], repeat=4)):
        a = alpha+ns*sigma*np.array(signs[:2])
        dynamic = 1+eta*np.array(signs[2:])
        b = -a*dynamic/(2*h*eps0*nu)
        d = b[1]-b[0]
        check(f"stress:{j}:both_redshift", np.all(b < 0))
        check(f"stress:{j}:nonuniversal_positive_ratio", d > 0)
        worst.append(d)
    analytic_min = ((alpha[0]-ns*sigma[0])*(1-eta)/nu[0]
                    -(alpha[1]+ns*sigma[1])*(1+eta)/nu[1])/(2*h*eps0)
    relative("stress:analytic_corner_bound", min(worst), analytic_min)

    # Unknown cross-paper covariance is not silently assumed zero.
    sigma_diff_bounds = [abs(beta_sigma[0]-beta_sigma[1]), sum(beta_sigma)]
    check("uncertainty:ordered_covariance_bounds", 0 < sigma_diff_bounds[0] <= sigma_diff_bounds[1])
    for corr in [-1.,0.,1.]:
        variance = beta_sigma[0]**2+beta_sigma[1]**2-2*corr*np.prod(beta_sigma)
        check(f"uncertainty:correlation={corr}",
              sigma_diff_bounds[0]*(1-1e-14) <= np.sqrt(variance) <= sigma_diff_bounds[1]*(1+1e-14))
    eta3 = cfg["transitions"][1]["eta_300K"]
    corrected_e3 = rows[0]["shifts_Hz"][1]*(1+eta3)
    check("dynamic:E3_known_correction_small", abs(eta3) < .01)
    check("dynamic:E3_does_not_reverse_sign", corrected_e3 < 0)
    # The equality could only be repaired by an enormous response change,
    # not the independently declared +/-10% dynamic stress.
    required_dynamic_ratio = 1/coefficient_ratio
    check("dynamic:equality_outside_stress",
          required_dynamic_ratio < (1-eta)/(1+eta))

    # Elementary second-order sign check, not an ab initio Yb calculation.
    gap, dipole, field = 2., .3, .002
    matrix = np.array([[0.,-dipole*field],[-dipole*field,gap]])
    levels = np.linalg.eigvalsh(matrix)
    polar_ground = 2*dipole**2/gap
    predicted = -.5*polar_ground*field**2
    relative("perturbation:Stark_factor_and_sign", levels[0], predicted, tol=1e-6)
    check("perturbation:ground_shift_negative", levels[0] < 0)
    for rel,digest in hashes.items():
        check("preserved:"+rel, hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == digest)
    failed = [row["name"] for row in checks if not row["passed"]]
    result = dict(stage="G-R4", status="bounded-source-audit-done" if not failed else "failed",
                  input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
                  response_units="fractional clock shift / (J m^-3)",
                  beta=beta.tolist(), beta_sigma=beta_sigma.tolist(),
                  coefficient_magnitude_ratio=float(coefficient_ratio), predictions=rows,
                  modulation_300_to_310=modulation,
                  E3_dynamic_corrected_shift_300K_Hz=corrected_e3,
                  covariance_free_sigma_difference_bounds=sigma_diff_bounds,
                  conditional_stress=dict(min_beta_difference=min(worst),
                                          min_ratio_offset_at_300K_leading=min(worst)*u300,
                                          status="3 published sigma box plus assumed +/-10% dynamics; not confidence or certified dynamic bound"),
                  conclusions=["selected thermal electric-dipole response is not one universal fractional shift",
                               "both shift signs are fixed by independently measured polarizabilities",
                               "same-bath ratio is blind to an additional exactly common multiplicative factor",
                               "no new experimental observations, no universal-time exclusion, no hardware forecast",
                               "ordinary spatial radiation density is not silently identified with GR3 or Chen Ed"],
                  previous_file_hashes=hashes, checks=checks,
                  summary=dict(checks=len(checks),passed=len(checks)-len(failed),failed=failed))
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    report = f"""# G-R4: physical thermal source and cross-transition test

Status: {result['status']}. Checks: {len(checks)-len(failed)}/{len(checks)}.
No fitted parameters or new experimental data.

## Selected source and independently fixed sign

Same-ion 171Yb+ E2 and E3 clock transitions in one ideal isotropic
blackbody environment. Microscopic coupling is -d.E. Published static
polarizabilities, converted to alpha_excited-alpha_ground, are
6.9(1.4)e-40 and 0.888(0.016)e-40 J m^2/V^2.
Both are positive, so the thermal electric-dipole shifts are negative.
The 2005 E2 paper uses the opposite polarizability sign convention.
Source files, independent calibration methods and caveats are pinned
in data/BBR_CLOCK_INPUTS.json; these are not claimed to be the newest
or most precise possible constants.

## Static leading-order prediction at 300 K

Thermal electromagnetic density u={u300:.6g} J/m^3;
E_rms={rms300:.6g} V/m. This assumes conventional spatial QED and
does not identify u with Chen Ed or the previous energy-per-cell candidate.

| Transition | Shift / Hz | Fractional shift |
|---|---:|---:|
| E2 | {rows[0]['shifts_Hz'][0]:.6g} | {rows[0]['fractional_shifts'][0]:.6g} |
| E3 | {rows[0]['shifts_Hz'][1]:.6g} | {rows[0]['fractional_shifts'][1]:.6g} |

The fractional responses differ by a factor {coefficient_ratio:.4g}.
R=(nu_E3/nu_E2) gives R(300)/R(0)-1={rows[0]['ratio_offset']:.6g}.
An illustrative 300 -> 310 K modulation gives
R(310)/R(300)-1={modulation:.6g}, not an observed signal.

Static leading-order quantities are deliberately separated from full
finite-temperature predictions. The published E3 eta(300)=-0.0015
changes its 300 K shift to {corrected_e3:.6g} Hz. A complete E2 dynamic
response and response covariance are not supplied.
The unequal leading T^4 coefficients already prevent exact equality
over a low-temperature interval; T^6 terms cannot change that coefficient.

## Robustness, not a significance fit

Across the full 3-published-sigma polarizability box and independent
assumed +/-10% dynamic factors, the leading ratio offset at 300 K
stays above {min(worst)*u300:.6g}. This is a conditional deterministic
stress statement, NOT a confidence level or a derived bound on every
neglected atomic correction. Unknown cross-paper covariance is retained.

## Essential interpretation limit

The ordinary thermal Stark channel cannot be represented solely as
H -> f(u) H + c(u) I on these three clock levels.
It does NOT exclude an additional exactly common time factor:
nu_i(u)=F(u) nu_i^0[1+s_i(u)] makes F cancel from nu_E3/nu_E2.
Subtracting published BBR corrections from both outputs and then finding
a constant ratio is not an independent test of those corrections.
No raw temperature-modulated observations have been analyzed.

## Next useful work

Before an empirical test, specify actual common-bath radiometry and
interleaved raw E2/E3 readings, retain the deliberately modulated BBR
term, and control temperature-correlated non-BBR shifts.
A claimed common time factor needs an independently specified reference
outside the changed environment or another non-clock observable.
Alternatively test equal total radiation density with different spectra;
the full response is spectral, not a scalar-density law.

Complete formulas/proofs: tex/route_g_bbr_universality.tex.
Run: python3 route_g/code/verify_gr4_bbr_universality.py.
"""
    OUT.with_suffix(".md").write_text(report)
    print(json.dumps(result["summary"]))
    print("300 K:",rows[0])
    print("300 -> 310 K ratio modulation:",modulation)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
