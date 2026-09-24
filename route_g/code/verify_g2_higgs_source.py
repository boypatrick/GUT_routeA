#!/usr/bin/env python3
"""G2-H: an externally normalized on-shell Higgs source, with a yield ceiling.

The Run-2 luminosity is real and the Higgs production normalization is SM
theory, supported by measured production. Dark production is hypothetical.
No Gaussian, invented source luminosity, transport MC, or xenon relocation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

import verify_g2_xenon_feasibility as xe

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_higgs_source"
C = xe.CARD
SOURCE = dict(
    collision_energy_TeV=13.0, integrated_luminosity_fb_inverse=139.0,
    sigma_h_SM_pb=55.6, sigma_h_SM_error_pb=2.5,
    sigma_h_stress_pb=60.0,
    normalization_source="https://arxiv.org/html/2207.08615v2",
    normalization_scope="SM production prediction, not an independent measured cross section under changed invisible branching fractions.",
    invisible_source="https://arxiv.org/abs/2301.10731v2",
    invisible_allowance=.107, pump="off for this separate source branch",
    xenon_mass_g=2e6, xenon_radius_cm=47.9, xenon_height_cm=97.0,
    geometry_source="https://arxiv.org/html/1907.11485v2",
    geometry_scope="Rounded active-volume cylinder for a single-pass upper bound; XENON1T is not at the LHC.",
    column_stress_nuclei_cm2=2e24,
    scale_domain_GeV=[1.0, 30.0],
    source_is_original_G2_collision=False,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def column_card():
    r, h = SOURCE["xenon_radius_cm"], SOURCE["xenon_height_cm"]
    rho = SOURCE["xenon_mass_g"] / (np.pi*r*r*h)
    length = np.hypot(2*r, h)
    ncol = rho * length * C["Avogadro"] / C["molar_mass_g_mol"]
    return dict(density_g_cm3=float(rho), maximum_chord_cm=float(length),
                maximum_column_nuclei_cm2=float(ncol),
                face_area_cm2=float(np.pi*r*r))


def sigma_majorant(mass, coupling, fN=None):
    """Any p, any position-dependent efficiency <=1, retained elastic ROI.

    Integrating a mathematical constant majorant to Emax does NOT assume
    coherent nuclear physics or calibrated response outside the retained ROI.
    """
    fn = C["fN"] if fN is None else fN
    ma = xe.ISOTOPES * C["atomic_mass_GeV"]
    return float(coupling**2*fn**2*C["mN_GeV"]**2/(4*np.pi*C["mh_GeV"]**4)
                 * np.dot(xe.ABUNDANCE, xe.ISOTOPES**2*ma**2/(mass+ma)**2)
                 * C["GeV_minus2_to_cm2"])


def roi_sigma(p, mass, coupling, rows, acceptance, order=16):
    """Actual public energy-response fold per nucleus at a specified momentum.

    For validation and a unit-column kernel only: no claim that Higgs bosons
    at a hadron collider are all at rest, or that a narrow beam has this
    volume-averaged spatial response.
    """
    if p <= 0:
        return 0.0
    z, w = leggauss(order)
    lo = np.maximum(rows["energy_bin_start_kev"], .7)*1e-6
    hi0 = np.minimum(rows["energy_bin_end_kev"], 50)*1e-6
    value = 0.0
    for A, frac in zip(xe.ISOTOPES, xe.ABUNDANCE):
        emax = xe.endpoint(p, mass, float(A)*C["atomic_mass_GeV"])
        hi = np.minimum(hi0, emax)
        half = np.maximum(hi-lo, 0)/2
        points = lo[:, None] + half[:, None]*(1+z[None, :])
        integrals = half * (xe.nuclear_H(points, int(A)) @ w)
        value += frac*np.dot(acceptance, integrals)
    return float(value*coupling**2/p**2*C["GeV_minus2_to_cm2"])


def source_card(scale):
    tower = xe.tower_width_unit(scale)
    width_unit = tower["unit_width_GeV"]
    gamma_cap = SOURCE["invisible_allowance"]/(1-SOURCE["invisible_allowance"])*C["Gamma_SM_GeV"]
    lam = min(1, 50*scale*scale/C["v_GeV"]**2,
              np.sqrt(gamma_cap/width_unit) if width_unit else np.inf)
    gamma = lam*lam*width_unit
    total = C["Gamma_SM_GeV"]+gamma
    mode_unit = C["v_GeV"]**2/(16*np.pi*C["mh_GeV"])
    modes = [dict(j=r["j"], mass_GeV=r["mass_GeV"], beta_star=r["beta"],
                  p_star_GeV=C["mh_GeV"]*r["beta"]/2,
                  branching_fraction=lam*lam*mode_unit*r["beta"]/total)
             for r in tower["open_signed_modes"]]
    return lam, gamma/total, modes


def run():
    checks = []

    def check(name, ok, measured=None):
        result = dict(name=name, passed=bool(ok))
        if measured is not None:
            result["measured"] = float(measured)
        checks.append(result)

    prior_path = ROOT/"output"/"g2_xenon_feasibility.json"
    prior = json.loads(prior_path.read_text())
    dependencies = {"xenon_source": Path(xe.__file__), "xenon_result": prior_path}
    hashes = {k: sha(p) for k, p in dependencies.items()}
    check("prior:no_failed_checks", not prior["summary"]["failed"])
    check("prior:source_hash", sha(xe.__file__) == prior["verification"]["source_sha256"])
    for name, expected in xe.RAW_HASHES.items():
        check("data:"+name, sha(xe.DATA/name) == expected)
    rows = np.genfromtxt(xe.DATA/"s2_response_nr.csv", delimiter=",", names=True)
    bins = np.genfromtxt(xe.DATA/"s2_binning_info.csv", delimiter=",", names=True)
    indices = np.flatnonzero((bins["start_pe"] >= 150) & (bins["end_pe"] <= 3000))
    acceptance = sum(rows[f"s2_bin_{i:03d}"] for i in indices)
    check("response:bounded", np.min(acceptance) >= 0 and np.max(acceptance) <= 1)
    geom = column_card()
    col = geom["maximum_column_nuclei_cm2"]
    nh = SOURCE["integrated_luminosity_fb_inverse"]*SOURCE["sigma_h_SM_pb"]*1000
    nx_ceiling = 2*nh*SOURCE["invisible_allowance"]
    check("source:pb_fb_conversion", abs(nh/7728400-1) < 1e-14)
    check("source:two_particles_per_decay", abs(nx_ceiling/1653877.6-1) < 1e-14)
    check("geometry:stress_exceeds_nominal", SOURCE["column_stress_nuclei_cm2"] > col)
    threshold = C["mh_GeV"]/10
    pair_threshold = C["mh_GeV"]/(2*np.sqrt(26))
    lambda_global = 50*threshold**2/C["v_GeV"]**2
    # Continuous envelope: lambda <=50 Lambda^2/v^2, mass_j>=5 Lambda,
    # and Lambda^4/(5 Lambda+mA)^2 is strictly increasing. Its threshold
    # supremum is valid even though production vanishes exactly there.
    sigma_global = sigma_majorant(C["mh_GeV"]/2, lambda_global)
    n_global = nx_ceiling*col*sigma_global
    n_stress = (2*SOURCE["integrated_luminosity_fb_inverse"]*SOURCE["sigma_h_stress_pb"]*1000
                *SOURCE["invisible_allowance"]*SOURCE["column_stress_nuclei_cm2"]
                *sigma_majorant(C["mh_GeV"]/2, lambda_global, C["fN"]+C["fN_error"]))
    anchors = []
    for scale in [1., 2., 5., 10., 12., 12.26, 12.4, threshold*(1-1e-6), threshold, 15., 30.]:
        lam, branching, modes = source_card(scale)
        nx = 2*nh*branching
        weighted = sum(r["branching_fraction"]*sigma_majorant(r["mass_GeV"],lam) for r in modes)
        bound = 2*nh*col*weighted
        selected = [r for r in modes if abs(r["j"]) == 1]
        selected_b = sum(r["branching_fraction"] for r in selected)
        for r in modes:
            check(f"Lambda={scale:g}:energy:j={r['j']}",
                  abs(np.hypot(r["mass_GeV"],r["p_star_GeV"])/(C["mh_GeV"]/2)-1)<1e-13)
        check(f"Lambda={scale:g}:branch_sum", abs(sum(r["branching_fraction"] for r in modes)-branching)<1e-13)
        check(f"Lambda={scale:g}:branch_allowance", branching <= SOURCE["invisible_allowance"]*(1+1e-12))
        check(f"Lambda={scale:g}:continuous_yield_bound", bound <= n_global*(1+1e-12))
        check(f"Lambda={scale:g}:source_below_budget", nx <= nx_ceiling*(1+1e-12))
        check(f"Lambda={scale:g}:on_shell_gate", bool(modes) == (scale < threshold))
        check(f"Lambda={scale:g}:absj1_gate", bool(selected) == (scale < pair_threshold))
        by_j = {r["j"]: r for r in modes}
        check(f"Lambda={scale:g}:signed_symmetry",
              all(abs(r["branching_fraction"]-by_j[-r["j"]]["branching_fraction"])<1e-15 for r in modes))
        anchors.append(dict(Lambda_GeV=scale, lambda_X=lam, branching_fraction=branching,
                            total_X_particles=nx, absj1_particles=2*nh*selected_b,
                            open_signed_modes=len(modes), boost_geometry_independent_event_ceiling=bound,
                            purpose=("threshold regression only" if threshold*(1-1e-5)<scale<threshold
                                     else "coverage anchor"),
                            modes=modes))

    # Direct numerical checks of the majorant with the actual response,
    # including slow particles whose endpoint cuts through the calibration bins.
    folds = []
    for mass in [5., 25., 60., C["mh_GeV"]/2]:
        for p in [.003, .01, .03, .1, 1., 60., 1000.]:
            value = roi_sigma(p,mass,.01,rows,acceptance)
            coarse = roi_sigma(p,mass,.01,rows,acceptance,order=8)
            bound = sigma_majorant(mass,.01)
            check(f"kernel:m={mass:g}:p={p:g}:bound", value <= bound*(1+1e-12))
            check(f"kernel:m={mass:g}:p={p:g}:quadrature", abs(value-coarse) <= max(value,1e-100)*1e-11)
            folds.append(dict(mass_GeV=mass,p_GeV=p,coupling=.01,
                              response_fold_cm2=value,efficiency_one_majorant_cm2=bound))
    check("kernel:zero_momentum", roi_sigma(0,5,.01,rows,acceptance)==0)
    # Recover the old fixed-source result only as a regression normalization.
    D = sum(prior["S2_spectrum"]["nuclear_integral_per_bin"])
    direct = roi_sigma(1.1,np.sqrt(26),1.,rows,acceptance)
    check("kernel:old_factorized_normalization",
          abs(direct/(D/1.1**2*C["GeV_minus2_to_cm2"])-1)<1e-12)
    # Symmetric Higgs decay, any shared boost: full four-momentum/time
    # likelihood is identical for opposite field labels, not just S2.
    symmetry = []
    mass = np.sqrt(26)*5
    estar = C["mh_GeV"]/2
    pstar = np.sqrt(estar*estar-mass*mass)
    for rapidity in [0., .8, 2.]:
        for costheta in [-.8, .0, .7]:
            energy = np.cosh(rapidity)*estar + np.sinh(rapidity)*pstar*costheta
            pz = np.sinh(rapidity)*estar + np.cosh(rapidity)*pstar*costheta
            pt = pstar*np.sqrt(1-costheta*costheta)
            check(f"boost:y={rapidity}:c={costheta}:mass_shell",
                  abs((energy**2-pz**2-pt**2)/mass**2-1)<1e-12)
            beta = np.hypot(pz,pt)/energy
            # Both labels use the same mass, production density and vertex.
            eplus, eminus = energy, energy
            tplus = 1/beta-1
            tminus = 1/beta-1
            check(f"sign:y={rapidity}:c={costheta}:time_momentum_equal",
                  eplus == eminus and tplus == tminus)
            symmetry.append(dict(rapidity=rapidity,costheta=costheta,
                                 energy_GeV=energy,time_delay_in_L_over_c=tplus))
    check("global:decisive_yield_ceiling", n_global < 4e-6, n_global)
    check("global:stress_ceiling", n_stress < 6e-6, n_stress)
    check("global:source_not_old_preparation", not SOURCE["source_is_original_G2_collision"])
    check("dependencies:unchanged", all(sha(p)==hashes[k] for k,p in dependencies.items()))
    failed = [r["name"] for r in checks if not r["passed"]]
    result = dict(
        status="Bounded rejection of Run-2 on-shell Higgs production plus single-pass, retained-window xenon elastic readout, not of all Route G.",
        source_card=SOURCE, geometry=geom,
        source_budget=dict(higgs_bosons_SM=nh, maximum_X_particles=nx_ceiling,
                           maximum_absj1_scale_GeV=pair_threshold,
                           maximum_any_mode_scale_GeV=threshold,
                           uniform_face_fluence_ceiling_cm2_inverse=nx_ceiling/geom["face_area_cm2"],
                           fluence_scope="All produced particles artificially spread over the cylinder face; no geometric transmission loss. Not an LHC-to-XENON prediction or a time flux."),
        global_bound=dict(coupling_supremum=lambda_global,cross_section_majorant_cm2=sigma_global,
                          event_ceiling=n_global,stress_event_ceiling=n_stress,
                          required_source_multiplier_for_three_at_least=3/n_global,
                          event_ceiling_if_Binv_is_one=n_global/SOURCE["invisible_allowance"],
                          proof="sigma_ROI <= lambda^2 fN^2 mN^2/(4 pi mh^4) sum xi A^2 mA^2/(mX+mA)^2; use mX>=5Lambda and lambda<=50Lambda^2/v^2. Lambda^4/(5Lambda+mA)^2 increases up to mh/10. Multiply by 2 Nh Binv and maximum single-pass column.",
                          no_threshold_optimization="Independent envelope factors need not be simultaneously saturated. The source closes at the endpoint; no tuned threshold point is claimed."),
        anchors=anchors,kernel_validation=folds,sign_validation=symmetry,
        sign_result=dict(conditional_priors=[.5,.5],optimal_error=.5,
                         mutual_information_bits=0,
                         scope="Direct tree-level Higgs source and diagonal X-only probe, j=+-1. Any common four-dimensional time/direction/momentum response; not the original asymmetric Phi-X preparation or full loop-corrected theory."),
        exclusions_of_scope=["Off-shell h* production and other sources are not bounded by the on-shell source budget.",
                             "Only the elastic NR component in .7..50 keV is bounded; electronic, inelastic, and high-energy response are not calculated.",
                             "No repeated circulation/trapping, particle multiplication, or new interaction; stable X makes at most one straight traversal.",
                             "The global event bound uses response<=1, not the volume-average matrix at an invalid narrow-beam position.",
                             "Do not multiply the earlier kg-day exposure by this integrated source budget; actual timing/geometry correlations would only reduce the stated single-pass ceiling.",
                             "This physical source does not prepare a Phi beam, implement the original driven quartic, or demonstrate G2 mode conversion.",
                             "Rounded geometry, scalar current, tree matching and SM Higgs production premises remain conditional; the stress case is not a combined confidence bound."],
        summary=dict(checks=len(checks),passed=len(checks)-len(failed),failed=failed),
        checks=checks,verification=dict(source_sha256=sha(__file__),dependencies=hashes,
                                       response_sha256=xe.RAW_HASHES))
    OUT.with_suffix(".json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    lines = ["# G2-H: an actual accelerator source budget rejects the xenon rescattering chain", "",
             f"Verification: **{len(checks)-len(failed)}/{len(checks)} checks passed**.", "",
             "Use one 13 TeV LHC interaction point with the published Run-2 luminosity of 139 fb^-1 and SM Higgs production 55.6 pb. The luminosity is actual; the production cross section is a theory input supported by measurements; X production is hypothetical.", "",
             f"N_h = {nh:.7g}; N_X <= 2 N_h B_inv = **{nx_ceiling:.7g}** for B_inv<=.107. This is a particle budget across the full dataset, not a fabricated cm^-2 s^-1 flux.", "",
             f"A 2 tonne xenon cylinder with radius47.9cm and height97cm gives maximum chord {geom['maximum_chord_cm']:.3f}cm and column {col:.6g} nuclei/cm². Assigning EVERY produced particle that column and 100% readout efficiency gives the continuous-scale ceiling **N_NR < {n_global:.6g}**. The stress input (60pb, fN=.326, column2e24) gives {n_stress:.6g}.", "",
             "This excludes the chosen source/readout chain as a useful event-count experiment even before geometry, attenuation, timing overlap, backgrounds or calibration losses. It is not an exclusion of the particle theory, off-shell sources, or unmodeled recoil channels.", "",
             "## Same-coupling source and detection", "",
             "Gamma_j=lambda_X² v² beta_j/(16 pi mh); B_j=Gamma_j/(Gamma_SM+sum Gamma_j). Source and scattering strengths are not independent: event yield scales as lambda_X^4/(Gamma_SM+lambda_X² Gamma_unit) at fixed masses.", "",
             "| Lambda [GeV] | lambda_X ceiling | Produced X+anti-X | Produced abs(j)=1 | Arbitrary-boost event upper bound |",
             "|---:|---:|---:|---:|---:|"]
    for r in anchors:
        if r["purpose"] != "coverage anchor":
            continue
        lines.append(f"| {r['Lambda_GeV']:.8g} | {r['lambda_X']:.5g} | {r['total_X_particles']:.5g} | {r['absj1_particles']:.5g} | {r['boost_geometry_independent_event_ceiling']:.5g} |")
    lines += ["", f"The j=+-1 on-shell source closes at Lambda={pair_threshold:.6g}GeV, earlier than the j=0 threshold{threshold:g}GeV. Above the latter, all on-shell Higgs production is closed even where the old coupling/flux map looked less restrictive.", "",
              "## Mode information is source-dependent", "",
              "Direct scalar Higgs decay gives equal populations and identical four-dimensional distributions for j=+1 and -1. The diagonal probe is likewise identical. Thus timing, direction AND independent momentum measurements have optimal sign error50% for this source; they cannot lift an exact sign degeneracy. This differs from the old asymmetric collision source's14.057% S2 prior error.", "",
              "## Stop and redirect", "",
              "Do not run detailed beam transport or narrower-Gaussian scans for this rejected chain. A collider production/missing-momentum search avoids the second tiny scattering probability, but on-shell missing-Higgs production alone measures an inclusive width, not signed KK modes. A threshold/recoil-mass strategy or a genuinely mode-sensitive interaction is a separate next hypothesis, not an already feasible instrument.", "",
              "Full source, transport functional, bound, sign theorem and limitations: tex/route_g_higgs_source.tex. No physical Phi preparation or GeV pump has been built."]
    OUT.with_suffix(".md").write_text("\n".join(lines)+"\n")
    print(json.dumps(result["summary"]))
    print(f"N_h={nh:g}, N_X cap={nx_ceiling:g}, global count cap={n_global:.6g}, stress={n_stress:.6g}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
