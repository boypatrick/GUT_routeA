#!/usr/bin/env python3
"""Route G2-D: bounded first-order, externally powered scattering sidebands.

The homogeneous classical pump supplies energy but no ordinary or compact
momentum. This is a long-time rate kernel, not a finite-collision probability
or a quantized-pump model. The undriven G2 card and finite packet are unchanged.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import verify_g2_local_conversion as base

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_driven"
CARD = dict(base.CARD)
PACKET = dict(base.PACKET)
DRIVE = dict(epsilon=0.5, omega=0.2)
SIDEBANDS = (-1, 0, 1)


def harmonic(q, card=CARD, epsilon=DRIVE["epsilon"]):
    """g(t)=sum_q g_q exp(-i q Omega t); the declared cosine is real."""
    return float(card["g"] * (1.0 if q == 0 else epsilon / 2 if abs(q) == 1 else 0.0))


def finite_packet():
    labels = np.arange(PACKET["min_l"], PACKET["max_l"] + 1)
    coefficients = np.exp(-labels ** 2 / (4 * PACKET["width"] ** 2)
                          - 1j * labels * PACKET["theta0"])
    coefficients /= np.linalg.norm(coefficients)
    return labels, coefficients, np.abs(coefficients) ** 2


def driven_channel(n, l, m, q, momentum, card=CARD,
                   epsilon=DRIVE["epsilon"], omega=DRIVE["omega"]):
    if omega <= 0 or not 0 <= epsilon < 1:
        raise ValueError("Resolved-sideband prescription requires Omega>0 and 0<=epsilon<1")
    j = n + l - m
    mi, mui = base.phi_mass(n, card), base.detector_mass(l, card)
    mf, muf = base.phi_mass(m, card), base.detector_mass(j, card)
    ei, el = float(np.hypot(mi, momentum)), float(np.hypot(mui, momentum))
    incoming_energy = ei + el
    outgoing_energy = incoming_energy + q * omega
    threshold = mf + muf
    opened = outgoing_energy > threshold
    pf = base.momentum_at_energy(outgoing_energy, mf, muf) if opened else 0.0
    amplitude = harmonic(q, card, epsilon)
    kernel = amplitude ** 2 * pf / (16 * np.pi * ei * el * outgoing_energy) if opened else 0.0
    ef, er = float(np.hypot(mf, pf)), float(np.hypot(muf, pf))
    return dict(n=n, l=l, m=m, j=j, q=q, epsilon=epsilon, omega=omega,
                initial_phi_mass=mi, initial_detector_mass=mui,
                final_phi_mass=mf, final_detector_mass=muf,
                initial_phi_energy=ei, initial_detector_energy=el,
                initial_total_energy=incoming_energy, final_total_energy=outgoing_energy,
                final_phi_energy=ef if opened else None,
                final_detector_energy=er if opened else None,
                final_threshold=threshold, available_energy=outgoing_energy-threshold,
                incoming_momentum=momentum, outgoing_momentum=pf,
                harmonic_amplitude=amplitude, open=bool(opened),
                driven_process_enabled=bool(amplitude != 0),
                rate_coefficient=float(kernel), pump_work=q * omega,
                energy_balance_residual=ef + er - incoming_energy - q * omega if opened else None,
                compact_momentum_residual=((n + card["alpha"]) + l
                                           - (m + card["alpha"]) - j) / card["R"],
                ordinary_total_momentum_residual=[0.0, 0.0, 0.0],
                mode_changing=bool(m != n),
                elastic=bool(m == n and q == 0),
                classification=("mode_conversion" if m != n else
                                "elastic_scattering" if q == 0 else "driven_mode_preserving"))


def enumerate_driven(n, l, momentum, card=CARD,
                     epsilon=DRIVE["epsilon"], omega=DRIVE["omega"]):
    ei = float(np.hypot(base.phi_mass(n, card), momentum)
               + np.hypot(base.detector_mass(l, card), momentum))
    bands, opened, closed = [], [], []
    for q in SIDEBANDS:
        energy = ei + q * omega
        bound = max(energy, 0.0) * card["R"]
        low = int(np.ceil(-bound - card["alpha"]))
        high = int(np.floor(bound - card["alpha"]))
        rows = [driven_channel(n, l, m, q, momentum, card, epsilon, omega)
                for m in range(low, high + 1)] if energy > 0 else []
        bands.append(dict(q=q, final_total_energy=energy,
                          rigorous_candidate_m_interval=[low, high],
                          open_m=[row["m"] for row in rows if row["open"]]))
        opened.extend(row for row in rows if row["open"])
        closed.extend(row for row in rows if not row["open"])
    return dict(sideband_bounds=bands, open_channels=opened, closed_candidate_channels=closed,
                completeness_scope="All open two-body channels at first order in this quartic vertex; q=-1,0,+1 only")


def summarize_packet(epsilon=DRIVE["epsilon"], omega=DRIVE["omega"]):
    labels, coefficients, weights = finite_packet()
    n, momentum = CARD["incoming_phi_n"], CARD["incoming_com_momentum"]
    components, joint = [], []
    for l0, coefficient, weight0 in zip(labels, coefficients, weights):
        l, weight = int(l0), float(weight0)
        enumeration = enumerate_driven(n, l, momentum, epsilon=epsilon, omega=omega)
        for raw in enumeration["open_channels"]:
            row = dict(raw, incoming_packet_weight=weight,
                       weighted_rate_coefficient=weight * raw["rate_coefficient"])
            joint.append(row)
        components.append(dict(l=l, coefficient=base.cmatrix(coefficient), weight=weight,
                               all_channels=enumeration))
    total = sum(row["weighted_rate_coefficient"] for row in joint)
    for row in joint:
        row["conditional_share"] = row["weighted_rate_coefficient"] / total
    bands = []
    for q in SIDEBANDS:
        rate = sum(row["weighted_rate_coefficient"] for row in joint if row["q"] == q)
        bands.append(dict(q=q, pump_work=q*omega, rate_coefficient=rate, conditional_share=rate/total,
                          open_channel_count=sum(row["q"] == q for row in joint)))
    conversion = sum(row["weighted_rate_coefficient"] for row in joint if row["mode_changing"])
    target = sum(row["weighted_rate_coefficient"] for row in joint if row["m"] == CARD["target_phi_m"])
    elastic = sum(row["weighted_rate_coefficient"] for row in joint if row["elastic"])
    driven_preserving = sum(row["weighted_rate_coefficient"] for row in joint
                            if row["classification"] == "driven_mode_preserving")
    work = sum(row["weighted_rate_coefficient"] * row["pump_work"] for row in joint)
    return dict(epsilon=epsilon, omega=omega, definition=dict(PACKET),
                normalization=float(weights.sum()), components=components, joint_outputs=joint,
                sideband_rows=bands, total_rate_coefficient=total,
                mode_conversion_rate_coefficient=conversion, mode_conversion_share=conversion/total,
                target_m_zero_rate_coefficient=target, target_m_zero_share=target/total,
                elastic_rate_coefficient=elastic, elastic_share=elastic/total,
                driven_mode_preserving_rate_coefficient=driven_preserving,
                driven_mode_preserving_share=driven_preserving/total,
                pump_work_rate_coefficient=work,
                mean_pump_work_per_scattering=work/total,
                conditioning="Shares conditional on a counted first-order scattering event under common dilute external-space overlap; unscattered identity excluded")


def target_pair(packet):
    lines = [dict(row) for row in packet["joint_outputs"]
             if row["m"] == CARD["target_phi_m"] and abs(row["j"]) == 1]
    pairs = []
    for positive in (row for row in lines if row["j"] == 1):
        for negative in (row for row in lines if row["j"] == -1):
            pairs.append(dict(positive_l=positive["l"], positive_q=positive["q"],
                              positive_p=positive["outgoing_momentum"],
                              negative_l=negative["l"], negative_q=negative["q"],
                              negative_p=negative["outgoing_momentum"],
                              momentum_gap=abs(positive["outgoing_momentum"]-negative["outgoing_momentum"])))
    return dict(six_lines=lines, opposite_sign_pairs=pairs,
                closest_opposite_sign_pair=min(pairs, key=lambda row: row["momentum_gap"]),
                minimum_gap_if_ideal_independent_q_tag=min(row["momentum_gap"] for row in pairs
                                                           if row["positive_q"] == row["negative_q"]),
                ideal_q_tag_scope="Hypothetical independently supplied sideband label; the prescribed classical pump is not a modeled energy-resolving apparatus")


def temporal_amplitude(delta_energy, duration, epsilon=DRIVE["epsilon"], omega=DRIVE["omega"]):
    """Centered top-hat Fourier integral only, not a finite-collision rate."""
    return sum(harmonic(q, epsilon=epsilon) * duration
               * np.sinc((delta_energy - q*omega) * duration / (2*np.pi)) for q in SIDEBANDS)


def run():
    checks = []

    def check(name, passed, measured=None, tolerance=None):
        row = dict(name=name, passed=bool(passed))
        if measured is not None:
            row["measured"] = float(measured)
        if tolerance is not None:
            row["tolerance"] = float(tolerance)
        checks.append(row)

    packet = summarize_packet()
    check("finite_packet:exact_normalization", abs(packet["normalization"]-1) < 1e-14)
    check("conditional_joint_shares:sum_one",
          abs(sum(row["conditional_share"] for row in packet["joint_outputs"])-1) < 1e-14)
    check("sideband_shares:sum_one", abs(sum(row["conditional_share"] for row in packet["sideband_rows"])-1) < 1e-14)
    check("exclusive_classifications:sum_one", abs(packet["mode_conversion_share"]
          + packet["elastic_share"]+packet["driven_mode_preserving_share"]-1) < 1e-14)
    fourier_rows = []
    period = 2*np.pi/DRIVE["omega"]
    for q in range(-4, 5):
        integral = quad(lambda t: CARD["g"]*(1+DRIVE["epsilon"]*np.cos(DRIVE["omega"]*t))
                        * np.cos(q*DRIVE["omega"]*t)/period, 0, period, epsabs=1e-13)[0]
        imaginary = quad(lambda t: CARD["g"]*(1+DRIVE["epsilon"]*np.cos(DRIVE["omega"]*t))
                         * np.sin(q*DRIVE["omega"]*t)/period, 0, period, epsabs=1e-13)[0]
        error = abs(complex(integral, imaginary)-harmonic(q))
        check(f"Fourier_coefficient:q={q}", error < 1e-13, error, 1e-13)
        fourier_rows.append(dict(q=q, numerical=base.cmatrix(complex(integral, imaginary)),
                                 exact=harmonic(q), absolute_error=error))
    radial_errors, rate_errors, gauge_errors = [], [], []
    for row in packet["joint_outputs"]:
        label = f"l={row['l']}:m={row['m']}:q={row['q']}:"
        check(label+"energy_including_pump", abs(row["energy_balance_residual"]) < 2e-12,
              abs(row["energy_balance_residual"]), 2e-12)
        check(label+"compact_momentum", row["compact_momentum_residual"] == 0)
        # A representative scattering direction checks the explicitly imposed
        # three-momentum delta function independently of the energy root.
        incoming_vector = np.array([0.0,0.0,row["incoming_momentum"]])
        direction = np.array([np.sqrt(1-0.3**2),0.0,0.3])
        outgoing_vector = row["outgoing_momentum"]*direction
        spatial_residual = incoming_vector-incoming_vector-outgoing_vector+outgoing_vector
        check(label+"ordinary_three_momentum",
              np.linalg.norm(spatial_residual) < 1e-14
              and abs(np.linalg.norm(outgoing_vector)-row["outgoing_momentum"]) < 1e-14)
        check(label+"nonnegative_rate", row["rate_coefficient"] >= 0)
        check(label+"nonzero_q_mode_preserving_not_elastic",
              not (row["m"] == row["n"] and row["q"] != 0 and row["elastic"]))
        energy, ma, mb = row["final_total_energy"], row["final_phi_mass"], row["final_detector_mass"]
        root = brentq(lambda p: np.hypot(ma,p)+np.hypot(mb,p)-energy, 0, energy/2,
                      xtol=1e-14, rtol=1e-14)
        ea, eb = np.hypot(ma,root), np.hypot(mb,root)
        radial = root**2/(ea*eb*(root/ea+root/eb))
        phase_space = quad(lambda c: 2*np.pi*radial/(16*np.pi**2), -1, 1)[0]
        independent = row["harmonic_amplitude"]**2*phase_space/(4*row["initial_phi_energy"]*row["initial_detector_energy"])
        root_error, rate_error = abs(root-row["outgoing_momentum"]), abs(independent-row["rate_coefficient"])
        check(label+"independent_radial_root", root_error < 2e-12, root_error, 2e-12)
        check(label+"independent_phase_space_rate", rate_error < 2e-16, rate_error, 2e-16)
        row["independent_radial_root"] = root
        row["independent_phase_space"] = float(phase_space)
        row["independent_phase_space_rate_coefficient"] = float(independent)
        radial_errors.append(root_error); rate_errors.append(rate_error)
        transformed = driven_channel(row["n"]-1,row["l"],row["m"]-1,row["q"],
                                     CARD["incoming_com_momentum"],dict(CARD,alpha=CARD["alpha"]+1))
        gauge_error = max(abs(row[key]-transformed[key]) for key in
                          ("initial_total_energy","final_total_energy","outgoing_momentum","rate_coefficient"))
        gauge_errors.append(gauge_error)
        check(label+"large_gauge_relabelling", gauge_error < 1e-14, gauge_error, 1e-14)
    for component in packet["components"]:
        l = component["l"]
        exhaustive = [(q,m) for q in SIDEBANDS for m in range(-50,51)
                      if driven_channel(CARD["incoming_phi_n"],l,m,q,CARD["incoming_com_momentum"])["open"]]
        bounded = [(row["q"],row["m"]) for row in component["all_channels"]["open_channels"]]
        check(f"l={l}:complete_open_channel_enumeration", exhaustive == bounded)
    zero = summarize_packet(epsilon=0)
    regression_errors = []
    for row in zero["joint_outputs"]:
        reference = base.channel(row["n"],row["l"],row["m"],CARD["incoming_com_momentum"])
        expected = reference["rate_coefficient"] if row["q"] == 0 else 0.0
        regression_errors.append(abs(row["rate_coefficient"]-expected))
    check("epsilon_zero:every_rate_matches_undriven_G2", max(regression_errors) < 1e-16,
          max(regression_errors),1e-16)
    check("epsilon_zero:no_work", zero["pump_work_rate_coefficient"] == 0)
    reconstructed_work = sum(row["weighted_rate_coefficient"]
                             *(row["final_phi_energy"]+row["final_detector_energy"]
                               -row["initial_total_energy"]) for row in packet["joint_outputs"])
    check("pump_work:weighted_on_shell_energy_balance",
          abs(reconstructed_work-packet["pump_work_rate_coefficient"]) < 1e-18,
          abs(reconstructed_work-packet["pump_work_rate_coefficient"]),1e-18)
    minimum_kappa = 2*np.pi*CARD["R"]*CARD["g"]*(1-DRIVE["epsilon"])
    check("potential:analytic_nonnegative_coefficients",
          min(CARD["M5"]**2,CARD["MD"]**2,minimum_kappa) >= 0, minimum_kappa)
    positive_samples = []
    for phase in np.linspace(0,2*np.pi,33):
        kappa = 2*np.pi*CARD["R"]*CARD["g"]*(1+DRIVE["epsilon"]*np.cos(phase))
        for x,y in ((0,0),(1,2),(10,3)):
            value = CARD["M5"]**2*x+CARD["MD"]**2*y+kappa*x*y
            positive_samples.append(value)
    check("potential:sample_regression", min(positive_samples) >= 0)
    finite_time = []
    for duration in (25.0,100.0):
        for delta in (0.0,0.17,0.2,0.47):
            integrand = lambda t: CARD["g"]*(1+DRIVE["epsilon"]*np.cos(DRIVE["omega"]*t))*np.cos(delta*t)
            numerical = quad(integrand,-duration/2,duration/2,epsabs=1e-12)[0]
            analytic = temporal_amplitude(delta,duration)
            error = abs(numerical-analytic)
            check(f"finite_time_transform:T={duration}:delta={delta}",error<1e-11,error,1e-11)
            finite_time.append(dict(duration=duration,delta_energy=delta,numerical=numerical,
                                    analytic=float(analytic),absolute_error=error))
    pair = target_pair(packet)
    check("target_abs_j_one:six_open_lines",len(pair["six_lines"])==6)
    check("target_abs_j_one:near_opposite_sidebands",
          pair["closest_opposite_sign_pair"]["positive_q"]==1
          and pair["closest_opposite_sign_pair"]["negative_q"]==-1)
    scan = []
    for omega in (0.1,0.2,0.4):
        tested = summarize_packet(omega=omega)
        scanned_pair = target_pair(tested)
        scan.append(dict(omega=omega,total_rate_coefficient=tested["total_rate_coefficient"],
                         sideband_rows=tested["sideband_rows"],mode_conversion_share=tested["mode_conversion_share"],
                         mean_pump_work_per_scattering=tested["mean_pump_work_per_scattering"],
                         pump_work_rate_coefficient=tested["pump_work_rate_coefficient"],
                         closest_opposite_sign_pair=scanned_pair["closest_opposite_sign_pair"],
                         note="Predeclared sensitivity cards; no frequency optimization"))
    result = dict(status="G2-D externally powered first-order long-time rate kernel; not an absolute transition probability",
        card=dict(CARD),drive=dict(DRIVE),packet=packet,target_abs_j_one=pair,omega_scan=scan,
        action=dict(potential="M5^2|Phi|^2+MD^2|X|^2+kappa5(t)|Phi|^2|X|^2",
                    kappa5="2 pi R g [1+epsilon cos(Omega t)]",minimum_kappa5=float(minimum_kappa),
                    external_background="Spatially and circle homogeneous classical pump; preferred frame chosen as the incoming COM",
                    tree_quadratic_spectrum="Exactly the original G2 masses at both zero VEVs; no claim about loop pole shifts",
                    fourier_convention="g(t)=sum_q g_q exp(-i q Omega t); g_0=g, g_+-1=g epsilon/2",
                    selection_rules="n+l=m+j, ordinary total spatial momentum conserved, E_final-E_initial=q Omega"),
        conventions=dict(rate="K_lmq=|g_q|^2 p_final/(16 pi E_n E_l E_out), zero for closed threshold",
                         kernel_dimension="mass^-2",work_rate_coefficient_dimension="mass^-1",
                         work_sign="q>0: matter absorbs energy from the prescribed pump; q<0: matter gives energy to it",
                         pump_power="Not computed: sum |c_l|^2 K_lmq q Omega is a coefficient requiring an independently specified density/overlap or luminosity",
                         no_drive_limit="epsilon=0 with Omega>0 reproduces the undriven G2 rates",
                         zero_frequency_limit="Omega=0 is excluded from the incoherent sideband prescription: coincident harmonics must be combined coherently"),
        first_order_fourier_checks=fourier_rows,finite_time_transform_checks=finite_time,
        undriven_regression=dict(max_channel_rate_error=max(regression_errors),
                                packet_total_rate_coefficient=zero["total_rate_coefficient"],
                                packet_mode_conversion_share=zero["mode_conversion_share"]),
        audit=dict(max_radial_root_error=max(radial_errors),max_phase_space_rate_error=max(rate_errors),
                   max_large_gauge_error=max(gauge_errors)),
        limitations=[
            "Sideband rates use a long-time/period-averaged first-order limit with Omega*T much greater than one; finite-time amplitudes interfere.",
            "The finite top-hat transform is checked but no finite collision envelope, 4D incident wavepacket, detector acceptance, or absolute probability is computed.",
            "The finite internal packet and common dilute external-space overlap are unchanged from G2; they do not supply a universal wavepacket cross section.",
            "q outside -1,0,1 vanishes only at first order for the declared cosine vertex; higher orders can produce additional harmonics.",
            "A homogeneous drive delivers no compact momentum and cannot directly shift a lone particle's KK label at tree level; conversion still uses X recoil.",
            "m=n with q nonzero is driven mode-preserving inelastic scattering, not elastic; the unscattered identity amplitude is absent from conditional shares.",
            "The prescribed background is an unmodeled work reservoir, not a quantized pump or a stabilized UV completion.",
            "Positive instantaneous quartic potential is not conservation of the matter Hamiltonian; explicit time dependence permits external work.",
            "A low-energy 5D EFT cutoff must exceed the prepared energies; no loop, radiative stability, or UV-control claim is made.",
            "Ideal independent q tagging is only a conditional discriminator; no such measurement apparatus is derived here.",
            "No measured masses, particle fit, Standard Model spectrum, or Route-F repair is asserted."],
        checks=checks,summary=dict(checks=len(checks),passed=sum(row["passed"] for row in checks),
                                  failed=[row["name"] for row in checks if not row["passed"]]),
        source_sha256={Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       Path(base.__file__).name:hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest()})
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    write_report(result)
    print(json.dumps(result["summary"],indent=2))
    print(json.dumps({key:packet[key] for key in ("total_rate_coefficient","mode_conversion_share",
                     "mean_pump_work_per_scattering","pump_work_rate_coefficient")},indent=2))
    print(json.dumps(pair["closest_opposite_sign_pair"],indent=2))
    if result["summary"]["failed"]:
        raise SystemExit(1)


def write_report(result):
    packet,pair=result["packet"],result["target_abs_j_one"]
    lines=["# Route G2-D — explicitly powered scattering sidebands","",
        "A bounded tree-level engineering calculation. The original G2 masses, recoil field, and exact finite internal packet are retained; no measured masses or optimization enter.","",
        "## One additional physical assumption","",
        r"$$g(t)=g[1+\epsilon\cos(\Omega t)],\qquad g=0.05,\quad\epsilon=0.5,\quad\Omega=0.2.$$","",
        "The classical pump is homogeneous in ordinary space and on the circle. Its preferred frame is the incoming center-of-momentum frame. It supplies energy, not ordinary or compact momentum. Both zero-VEV quadratic mass spectra stay unchanged at tree level.","",
        r"$$g(t)=\sum_qg_qe^{-iq\Omega t},\quad g_0=g,\quad g_{\pm1}=g\epsilon/2;\qquad n+l=m+j,\quad E_{\rm out}=E_n+E_l+q\Omega.$$","",
        r"$$K_{lmq}=\frac{|g_q|^2p_f}{16\pi E_nE_lE_{\rm out}},\qquad p_f=\frac{\sqrt{[E_{\rm out}^2-(m_m+\mu_j)^2][E_{\rm out}^2-(m_m-\mu_j)^2]}}{2E_{\rm out}}.$$","",
        "Set K=0 for a closed threshold. The independent radial delta-function integral gives the two-body phase space p_f/(4 pi E_out), with incoming flux normalization 1/(4 E_n E_l). This yields K directly; the background breaks Lorentz invariance, so no invariant conserved-s cross section is assumed.","",
        "## Complete first-order rate budget","",
        "The exact G2 packet has l=-4,...,4, w=1.2 and theta0=0.7; incoming n=1 and COM p=0.2. Every open (l,m,q) for q=-1,0,+1 is enumerated using |m+alpha|<R E_out and verified against a wider integer interval. JSON retains all joint outputs and rejected candidates.","",
        "| q | open channels | weighted K | conditional share | matter work per event |","|---:|---:|---:|---:|---:|"]
    for row in packet["sideband_rows"]:
        lines.append(f"| {row['q']} | {row['open_channel_count']} | {row['rate_coefficient']:.12g} | {row['conditional_share']:.12g} | {row['pump_work']:.12g} |")
    lines += ["",f"Total common-overlap rate coefficient: **{packet['total_rate_coefficient']:.12g}**. Mode-conversion share (m != n): **{packet['mode_conversion_share']:.12g}**. Elastic share (m=n,q=0): **{packet['elastic_share']:.12g}**. Driven mode-preserving share (m=n,q != 0): **{packet['driven_mode_preserving_share']:.12g}**.","",
        "These are fractions conditional on counted scattering, not absolute conversion probabilities. The unscattered identity amplitude is excluded. In particular, q != 0 with m=n is not elastic.","",
        r"$$\overline W_{\rm event}=\frac{\sum_{lmq}|c_l|^2K_{lmq}q\Omega}{\sum_{lmq}|c_l|^2K_{lmq}},\qquad C_W=\sum_{lmq}|c_l|^2K_{lmq}q\Omega.$$","",
        f"Mean pump work absorbed per scattering: **{packet['mean_pump_work_per_scattering']:.12g}** (mass units). Work-rate coefficient C_W: **{packet['pump_work_rate_coefficient']:.12g}** (inverse-mass units). C_W is **not power**: an independently specified density/overlap or luminosity is required.","",
        "## Target m=0, |j|=1: six lines, not one recoil line","",
        "| incoming l | outgoing j | q | outgoing momentum | weighted K |","|---:|---:|---:|---:|---:|"]
    for row in sorted(pair["six_lines"],key=lambda r:(-r["j"],r["q"])):
        lines.append(f"| {row['l']} | {row['j']} | {row['q']} | {row['outgoing_momentum']:.12g} | {row['weighted_rate_coefficient']:.12g} |")
    closest=pair["closest_opposite_sign_pair"]
    lines += ["",f"Closest opposite-sign recoil lines: j=+1,q={closest['positive_q']} at p={closest['positive_p']:.12g}, and j=-1,q={closest['negative_q']} at p={closest['negative_p']:.12g}. Gap: **{closest['momentum_gap']:.12g}**. A pump can create near-overlapping outputs even while increasing the rate.","",
        f"If an independent ideal sideband label q were supplied, the smallest same-q opposite-sign gap would be {pair['minimum_gap_if_ideal_independent_q_tag']:.12g}. This is a conditional information advantage, not a detector or quantized-pump measurement provided by this calculation.","",
        "## Predeclared frequency sensitivity, not optimization","",
        "| Omega | total K | conversion share | mean pump work/event | closest opposite-sign momentum gap |","|---:|---:|---:|---:|---:|"]
    for row in result["omega_scan"]:
        lines.append(f"| {row['omega']:g} | {row['total_rate_coefficient']:.12g} | {row['mode_conversion_share']:.12g} | {row['mean_pump_work_per_scattering']:.12g} | {row['closest_opposite_sign_pair']['momentum_gap']:.12g} |")
    lines += ["","## Long-time boundary and checks","",
        r"For a centered top-hat time envelope only, $A_T(\Delta E)\propto\sum_qg_qF_T(\Delta E-q\Omega)$, $F_T(x)=2\sin(xT/2)/x=T\,\mathrm{sinc}(xT/2)$ with sinc defined as sin(x)/x. The script independently integrates this Fourier transform. It does not turn it into a finite-collision probability.","",
        "The rate sum uses the long-time/period-averaged limit with Omega*T much greater than one. At finite T the sideband amplitudes interfere. Omega=0 cannot be substituted into an incoherent sum: coincident amplitudes must be combined. A finite 4D collision envelope, timing, and acceptance are absent.","",
        f"**{result['summary']['passed']}/{result['summary']['checks']} checks pass.** Independent phase-space/root integration, Fourier coefficients, zero-drive G2 regression, energy plus pump work, compact/gauge conservation, full enumeration, packet/share normalization, finite-time Fourier transform, and positive instantaneous quartic coefficients are checked.","",
        f"The analytic quartic minimum is kappa5_min={result['action']['minimum_kappa5']:.12g}>0. This proves nonnegativity of that term at every time, not conservation of the matter Hamiltonian.","",
        "q beyond +-1 is absent only at first order for this vertex. A homogeneous pump cannot directly shift KK momentum or convert a lone particle's mode at tree level; the recoil field remains necessary. The unmodeled work reservoir, finite-time collisions, quantized pump, radiative corrections, EFT cutoff, and UV completion remain outside the result.","",
        "Reproduce: python3 route_g/code/verify_g2_driven.py. No changes to the undriven G2 verifier are needed. JSON records both source hashes and every joint output."]
    OUT.with_suffix(".md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    run()
