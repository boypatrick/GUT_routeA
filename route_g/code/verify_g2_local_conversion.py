#!/usr/bin/env python3
"""Route G G2: tree-level local, dynamical, recoil-carrying mode conversion.

Engineering example only: no observed masses, fit, background work reservoir,
static defect, loop claim, or identification with the Standard Model.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_local_conversion"
CARD = dict(R=1.0, M5=0.5, alpha=0.25, MD=5.0, g=0.05,
            incoming_phi_n=1, target_phi_m=0, incoming_com_momentum=0.2)
PACKET = dict(min_l=-4, max_l=4, width=1.2, theta0=0.7)


def cmatrix(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def phi_mass(n, card=CARD):
    return float(np.hypot(card["M5"], (n + card["alpha"]) / card["R"]))


def detector_mass(l, card=CARD):
    return float(np.hypot(card["MD"], l / card["R"]))


def momentum_at_energy(energy, mass_a, mass_b):
    if energy <= mass_a + mass_b:
        return 0.0
    # Factorized Kallen formula avoids subtraction of nearly equal squares.
    product = ((energy - mass_a - mass_b) * (energy + mass_a + mass_b)
               * (energy - mass_a + mass_b) * (energy + mass_a - mass_b))
    return float(np.sqrt(product) / (2 * energy))


def at_energy(energy, initial_a, initial_b, final_a, final_b, g):
    pi = momentum_at_energy(energy, initial_a, initial_b)
    pf = momentum_at_energy(energy, final_a, final_b)
    if pi <= 0:
        raise ValueError("This verifier requires a strictly above-threshold incoming state")
    sigma = g ** 2 * pf / (16 * np.pi * energy ** 2 * pi)
    vrel = pi / np.hypot(initial_a, pi) + pi / np.hypot(initial_b, pi)
    return dict(incoming_momentum=pi, outgoing_momentum=pf,
                sigma=float(sigma), relative_velocity=float(vrel),
                rate_coefficient=float(sigma * vrel))


def channel(n, l, m, momentum, card=CARD):
    j = n + l - m
    initial_a, initial_b = phi_mass(n, card), detector_mass(l, card)
    final_a, final_b = phi_mass(m, card), detector_mass(j, card)
    ei = float(np.hypot(initial_a, momentum))
    ed = float(np.hypot(initial_b, momentum))
    energy = ei + ed
    kinematics = at_energy(energy, initial_a, initial_b, final_a, final_b, card["g"])
    pf = kinematics["outgoing_momentum"]
    ef, er = float(np.hypot(final_a, pf)), float(np.hypot(final_b, pf))
    threshold = final_a + final_b
    opened = energy > threshold
    return dict(n=n, l=l, m=m, j=j, initial_phi_mass=initial_a,
                initial_detector_mass=initial_b, final_phi_mass=final_a,
                final_detector_mass=final_b, initial_phi_energy=ei,
                initial_detector_energy=ed, sqrt_s=energy, s=energy ** 2,
                final_threshold=threshold, available_energy=energy - threshold,
                open=opened, final_phi_energy=ef, final_detector_energy=er,
                energy_residual=ef + er - energy if opened else None,
                phi_energy_change=ef - ei if opened else None,
                detector_energy_change=er - ed if opened else None,
                compact_momentum_residual=((n + card["alpha"]) + l
                                           - (m + card["alpha"]) - j) / card["R"],
                **kinematics)


def enumerate_channels(n, l, momentum, card=CARD):
    energy = np.hypot(phi_mass(n, card), momentum) + np.hypot(detector_mass(l, card), momentum)
    # Any open final channel has m_phi < sqrt(s), hence |m+alpha| < R sqrt(s).
    # This deliberately loose finite interval is rigorous; exact thresholds filter it.
    bound = card["R"] * energy
    low = int(np.ceil(-bound - card["alpha"]))
    high = int(np.floor(bound - card["alpha"]))
    candidates = [channel(n, l, m, momentum, card) for m in range(low, high + 1)]
    return dict(rigorous_candidate_m_interval=[low, high],
                bound="Any open channel obeys |m+alpha| < R*sqrt(s).",
                open_channels=[row for row in candidates if row["open"]],
                closed_candidate_channels=[row for row in candidates if not row["open"]])


def run():
    checks = []

    def check(name, passed, measured=None, tolerance=None):
        row = dict(name=name, passed=bool(passed))
        if measured is not None:
            row["measured"] = float(measured)
        if tolerance is not None:
            row["tolerance"] = float(tolerance)
        checks.append(row)

    n, m, momentum = (CARD["incoming_phi_n"], CARD["target_phi_m"],
                      CARD["incoming_com_momentum"])
    kappa5 = 2 * np.pi * CARD["R"] * CARD["g"]
    overlaps = []
    for l in (-2, 0, 3):
        for target_m in (-2, 0, 1, 3):
            allowed_j = n + l - target_m
            for offset in (-2, -1, 0, 1, 2):
                j = allowed_j + offset
                exponent = n + l - target_m - j
                prefactor = kappa5 * CARD["R"] / (2 * np.pi * CARD["R"]) ** 2
                re = quad(lambda theta: prefactor * np.cos(exponent * theta),
                          0, 2 * np.pi, epsabs=1e-13)[0]
                im = quad(lambda theta: prefactor * np.sin(exponent * theta),
                          0, 2 * np.pi, epsabs=1e-13)[0]
                numerical = complex(re, im)
                exact = CARD["g"] if exponent == 0 else 0.0
                error = abs(numerical - exact)
                check(f"circle_overlap:l={l}:m={target_m}:j={j}", error < 1e-13, error, 1e-13)
                overlaps.append(dict(n=n, l=l, m=target_m, j=j,
                                     compact_integer_mismatch=exponent,
                                     overlap=cmatrix(numerical), exact=exact, error=error))

    labels = np.arange(PACKET["min_l"], PACKET["max_l"] + 1)
    coefficients = np.exp(-labels ** 2 / (4 * PACKET["width"] ** 2)
                          - 1j * labels * PACKET["theta0"])
    coefficients /= np.linalg.norm(coefficients)
    weights = np.abs(coefficients) ** 2
    check("finite_packet:exact_normalization", abs(weights.sum() - 1) < 1e-14,
          abs(weights.sum() - 1), 1e-14)
    components = []
    for index, l0 in enumerate(labels):
        l = int(l0)
        target = channel(n, l, m, momentum)
        enumeration = enumerate_channels(n, l, momentum)
        rows = enumeration["open_channels"]
        total = sum(row["rate_coefficient"] for row in rows)
        changed = sum(row["rate_coefficient"] for row in rows if row["m"] != n)
        components.append(dict(l=l, coefficient=cmatrix(coefficients[index]), weight=float(weights[index]),
                               target=target, all_channels=enumeration, total_rate_coefficient=total,
                               mode_changing_rate_coefficient=changed))
        check(f"l={l}:target_compact_conservation", target["compact_momentum_residual"] == 0,
              target["compact_momentum_residual"])
        for row in rows:
            prefix = f"l={l}:m={row['m']}:"
            check(prefix + "four_energy_conservation", abs(row["energy_residual"]) < 2e-12,
                  abs(row["energy_residual"]), 2e-12)
            check(prefix + "compact_momentum_conservation", row["compact_momentum_residual"] == 0)
            check(prefix + "threshold", row["available_energy"] > 0)
            energy, ma, mb = row["sqrt_s"], row["final_phi_mass"], row["final_detector_mass"]
            # Independent radial phase-space root and delta-function Jacobian.
            pf_root = brentq(lambda p: np.hypot(ma, p) + np.hypot(mb, p) - energy,
                             0, energy / 2, xtol=1e-14, rtol=1e-14)
            root_error = abs(pf_root - row["outgoing_momentum"])
            check(prefix + "radial_root_vs_Kallen", root_error < 2e-12, root_error, 2e-12)
            ea, eb = np.hypot(ma, pf_root), np.hypot(mb, pf_root)
            delta_jacobian = pf_root / ea + pf_root / eb
            radial = pf_root ** 2 / (ea * eb * delta_jacobian)
            phase_space = quad(lambda cos_theta: 2 * np.pi * radial / (16 * np.pi ** 2), -1, 1)[0]
            independent_sigma = CARD["g"] ** 2 * phase_space / (4 * row["incoming_momentum"] * energy)
            sigma_error = abs(independent_sigma - row["sigma"])
            check(prefix + "phase_space_sigma", sigma_error < 2e-15, sigma_error, 2e-15)
            inverse = at_energy(energy, ma, mb, row["initial_phi_mass"],
                                row["initial_detector_mass"], CARD["g"])
            balance = abs(row["incoming_momentum"] ** 2 * row["sigma"]
                          - row["outgoing_momentum"] ** 2 * inverse["sigma"])
            check(prefix + "inverse_detailed_balance", balance < 1e-16, balance, 1e-16)
            row["independent_phase_space"] = float(phase_space)
            row["radial_root_momentum"] = float(pf_root)
            row["inverse_sigma_at_same_s"] = inverse["sigma"]
            row["detailed_balance_residual"] = float(balance)
        # Independently scan a larger interval to confirm the analytic enumeration.
        wider = [candidate for candidate in range(-50, 51)
                 if channel(n, l, candidate, momentum)["open"]]
        enumerated = [row["m"] for row in rows]
        check(f"l={l}:complete_open_channel_enumeration", wider == enumerated)

    main_example = channel(n, 0, m, momentum)
    reverse_low_energy = channel(m, 1, n, momentum)
    check("example:target_open", main_example["open"])
    check("example:nonzero_conversion_rate", main_example["rate_coefficient"] > 0)
    check("reverse_at_same_small_momentum:closed", not reverse_low_energy["open"])
    check("reverse_at_same_small_momentum:zero_rate", reverse_low_energy["rate_coefficient"] == 0)
    zero_card = dict(CARD, g=0.0)
    check("zero_interaction:zero_rate", channel(n, 0, m, momentum, zero_card)["rate_coefficient"] == 0)

    rest_energy = phi_mass(n) + detector_mass(0)
    rest_final_p = momentum_at_energy(rest_energy, phi_mass(m), detector_mass(n - m))
    rate_limit = (CARD["g"] ** 2 * rest_final_p / (16 * np.pi * rest_energy ** 2)
                  * (1 / phi_mass(n) + 1 / detector_mass(0)))
    small_momentum_rows = []
    for incoming_p in (0.02, 0.002, 0.0002):
        row = channel(n, 0, m, incoming_p)
        relative_error = abs(row["rate_coefficient"] / rate_limit - 1)
        small_momentum_rows.append(dict(incoming_p=incoming_p, sigma=row["sigma"],
                                        rate_coefficient=row["rate_coefficient"],
                                        relative_error_to_finite_limit=relative_error))
    check("exothermic:finite_sigma_v_limit", small_momentum_rows[-1]["relative_error_to_finite_limit"] < 1e-7,
          small_momentum_rows[-1]["relative_error_to_finite_limit"], 1e-7)
    check("exothermic:decreasing_sigma_v_limit_error",
          all(a["relative_error_to_finite_limit"] > b["relative_error_to_finite_limit"]
              for a, b in zip(small_momentum_rows, small_momentum_rows[1:])))
    delta_mass = phi_mass(n) - phi_mass(m)
    critical_MD = ((n - m) ** 2 / CARD["R"] ** 2 - delta_mass ** 2) / (2 * delta_mass)
    threshold_rows = []
    for factor in (0.9, 1.0, 1.1):
        md = critical_MD * factor
        release = phi_mass(n) + md - phi_mass(m) - np.hypot(md, (n - m) / CARD["R"])
        expected = -1 if factor < 1 else (1 if factor > 1 else 0)
        passed = abs(release) < 1e-13 if expected == 0 else expected * release > 0
        check(f"heavy_recoil_threshold:factor={factor}", passed, release, 1e-13 if expected == 0 else None)
        threshold_rows.append(dict(factor_of_critical_mass=factor, MD=md,
                                   exact_rest_energy_release=float(release), expected_sign=expected))

    gauge_errors = []
    for component in components:
        original = component["target"]
        transformed = channel(n - 1, component["l"], m - 1, momentum,
                              dict(CARD, alpha=CARD["alpha"] + 1))
        error = max(abs(original[key] - transformed[key])
                    for key in ("sqrt_s", "final_threshold", "outgoing_momentum", "sigma", "rate_coefficient"))
        gauge_errors.append(error)
    check("large_gauge:alpha_plus_one_relabelled_modes", max(gauge_errors) < 1e-14,
          max(gauge_errors), 1e-14)

    target_rates = np.array([row["target"]["rate_coefficient"] for row in components])
    total_rates = np.array([row["total_rate_coefficient"] for row in components])
    changing_rates = np.array([row["mode_changing_rate_coefficient"] for row in components])
    packet_target = float(weights @ target_rates)
    packet_total = float(weights @ total_rates)
    packet_changing = float(weights @ changing_rates)
    rho = np.outer(coefficients, coefficients.conj())
    dephased = np.diag(weights)
    recoil_labels = labels + n - m
    # Rows carry distinct j; columns carry l. A is a declared rate map, not an
    # on-shell amplitude matrix at a common s: the packet components have distinct s_l.
    rate_map = np.diag(np.sqrt(target_rates))
    output = rate_map @ rho @ rate_map.conj().T
    output_dephased = rate_map @ dephased @ rate_map.conj().T
    inclusive = float(np.trace(output).real)
    inclusive_dephased = float(np.trace(output_dephased).real)
    check("orthogonal_recoil:coherent_inclusive_equals_weighted_sum", abs(inclusive - packet_target) < 1e-16,
          abs(inclusive - packet_target), 1e-16)
    check("orthogonal_recoil:dephasing_preserves_inclusive_target_rate", abs(inclusive - inclusive_dephased) < 1e-16,
          abs(inclusive - inclusive_dephased), 1e-16)
    # Off-diagonal entries of this bookkeeping map need not vanish. It does not
    # determine physical recoil coherences after tracing 3-momentum or timing.
    conditional_incoming_l = float((weights * target_rates) @ labels / packet_target)
    conditional_outgoing_j = float(np.diag(output).real @ recoil_labels / packet_target)
    shift_residual = abs(conditional_outgoing_j - conditional_incoming_l - (n - m))
    check("recoil:conditional_mean_shift", shift_residual < 1e-14, shift_residual, 1e-14)
    translated = coefficients * np.exp(-1j * labels * 1.3)
    translated_rho = np.outer(translated, translated.conj())
    translated_rate = float(np.trace(rate_map @ translated_rho @ rate_map.conj().T).real)
    check("incident_plane_mode:packet_translation_preserves_inclusive_rate",
          abs(translated_rate - packet_target) < 1e-16, abs(translated_rate - packet_target), 1e-16)

    angular_samples = 128
    theta = 2 * np.pi * np.arange(angular_samples) / angular_samples
    fourier = np.exp(1j * np.outer(theta, labels)) / np.sqrt(2 * np.pi * CARD["R"])
    energies = np.sqrt(CARD["MD"] ** 2 + (labels / CARD["R"]) ** 2 + momentum ** 2)
    mean_energy = float(weights @ energies)
    spreading = []
    for time in (0.0, 1.0, 5.0):
        evolved = coefficients * np.exp(-1j * energies * time)
        psi = fourier @ evolved
        hpsi = fourier @ (energies * evolved)
        norm = float(2 * np.pi * CARD["R"] * np.mean(abs(psi) ** 2))
        grid_energy = float((2 * np.pi * CARD["R"] * np.mean(psi.conj() * hpsi)).real)
        circle_moment_grid = 2 * np.pi * CARD["R"] * np.mean(abs(psi) ** 2 * np.exp(1j * theta))
        circle_moment_modes = np.sum(evolved[1:].conj() * evolved[:-1])
        moment_error = abs(circle_moment_grid - circle_moment_modes)
        check(f"packet:t={time}:norm_conserved", abs(norm - 1) < 1e-13, abs(norm - 1), 1e-13)
        check(f"packet:t={time}:mean_energy_conserved", abs(grid_energy - mean_energy) < 1e-12,
              abs(grid_energy - mean_energy), 1e-12)
        check(f"packet:t={time}:circle_moment_fourier_identity", moment_error < 1e-13, moment_error, 1e-13)
        spreading.append(dict(time=time, exact_finite_fourier_quadrature_samples=angular_samples,
                              norm=norm, mean_energy=grid_energy,
                              first_circular_moment=cmatrix(circle_moment_grid),
                              first_circular_moment_magnitude=float(abs(circle_moment_grid)),
                              circular_mean_angle=float(np.angle(circle_moment_grid))))
    check("packet:initial_localization", spreading[0]["first_circular_moment_magnitude"] > 0.8,
          spreading[0]["first_circular_moment_magnitude"])
    check("packet:initial_center", abs(spreading[0]["circular_mean_angle"] - PACKET["theta0"]) < 1e-13,
          abs(spreading[0]["circular_mean_angle"] - PACKET["theta0"]), 1e-13)
    check("packet:free_spreading_not_static_defect",
          spreading[2]["first_circular_moment_magnitude"] < spreading[0]["first_circular_moment_magnitude"])

    positive_potential_checks = []
    for phi_abs_squared, x_abs_squared in ((0, 0), (0, 2), (1, 0), (1, 2), (10, 3)):
        potential = (CARD["M5"] ** 2 * phi_abs_squared + CARD["MD"] ** 2 * x_abs_squared
                     + kappa5 * phi_abs_squared * x_abs_squared)
        positive_potential_checks.append(dict(phi_abs_squared=phi_abs_squared,
                                              x_abs_squared=x_abs_squared, potential=potential))
        check(f"potential:phi2={phi_abs_squared}:X2={x_abs_squared}", potential >= 0, potential)
    check("classical_potential:nonnegative_coefficients",
          min(CARD["M5"] ** 2, CARD["MD"] ** 2, kappa5) >= 0)

    result = dict(
        status="bounded tree-level G2 result; not a particle-physics fit or UV completion",
        action=dict(background="Minkowski(3,1) x circle of radius R; external flat U(1) holonomy alpha",
                    fields="Complex Phi of U(1) charge one; complex X neutral under that U(1); both VEVs zero",
                    potential="M5^2 |Phi|^2 + MD^2 |X|^2 + kappa5 |Phi|^2 |X|^2",
                    kinetic="Positive canonical 5D kinetic terms; Dtheta Phi=(partial_theta+i alpha)Phi; X periodic",
                    vertex="-i g delta_(n+l,m+j), g=kappa5/(2 pi R)",
                    validity="5D EFT at tree level; use below an independently chosen EFT cutoff; no loops or cutoff matching claimed"),
        conventions=dict(units="One arbitrary engineering mass unit; R in inverse mass; kappa5 has mass dimension -1; g dimension 0",
                         sigma="g^2 p_final / (16 pi s p_initial), distinct scalar particles, no identical-particle factor",
                         relative_velocity="COM Moller velocity p_initial/E_phi + p_initial/E_X",
                         packet_rate="sum_l |c_l|^2 sigma_l v_l under equal dilute external-space overlap; not a universal wavepacket cross section",
                         packet_position="Canonical angular one-particle wavefunction; not a covariant relativistic position density"),
        card=dict(**CARD, kappa5=float(kappa5)), main_example=main_example,
        closed_reverse_at_same_small_momentum=reverse_low_energy,
        exothermic_small_momentum=dict(analytic_finite_sigma_v_limit=float(rate_limit),
                                      rows=small_momentum_rows,
                                      note="sigma scales as 1/p_initial; sigma*v is finite. Neither is a per-particle probability."),
        heavy_recoil_threshold=dict(delta_phi_mass=delta_mass, analytic_critical_MD=float(critical_MD),
                                     formula="MD > ((n-m)^2/R^2 - (mn-mm)^2)/(2(mn-mm)), for l=0 and mn>mm",
                                     rows=threshold_rows),
        circle_overlap_checks=overlaps,
        packet=dict(**PACKET, definition="Finite EXACT normalized state, l=-4,...,4; c_l proportional to exp[-l^2/(4 w^2)-i l theta0]; no infinite-tail approximation",
                    normalization=float(weights.sum()), components=components,
                    mean_detector_energy=mean_energy, free_spreading=spreading,
                    target_rate_coefficient=packet_target, total_rate_coefficient=packet_total,
                    all_mode_changing_rate_coefficient=packet_changing,
                    target_fraction_of_all_scattering_rates=packet_target / packet_total,
                    target_fraction_of_mode_changing_rates=packet_target / packet_changing,
                    target_is_not_exclusive=True),
        recoil_rate_map=dict(row_j_labels=recoil_labels.tolist(), column_l_labels=labels.tolist(),
                             definition="A_(j,l)=sqrt(K_l) delta_(j,l+n-m), a common-overlap rate map, NOT a single-energy S-matrix",
                             A=cmatrix(rate_map), incoming_pure_density=cmatrix(rho),
                             outgoing_rate_matrix=cmatrix(output),
                             outgoing_dephased_rate_matrix=cmatrix(output_dephased),
                             inclusive_target_rate=inclusive, dephased_target_rate=inclusive_dephased,
                             translated_target_rate=translated_rate,
                             rate_conditioned_incoming_l_mean=conditional_incoming_l,
                             rate_conditioned_outgoing_j_mean=conditional_outgoing_j,
                             conditional_shift=conditional_outgoing_j - conditional_incoming_l,
                             note="Inclusive trace erases interference between orthogonal recoil labels. Off-diagonal entries in this bookkeeping map do not determine physical recoil coherence after tracing momentum or timing."),
        potential_samples=positive_potential_checks,
        limits=["Static homogeneous free backgrounds do not convert modes; this added local dynamical interaction does.",
                "No time-dependent background or external work reservoir is used. Recoil conserves energy and total compact momentum.",
                "Both masses and all channels follow the one shared radius and holonomy; masses were not individually adjusted.",
                "The initial packet is a prepared quantum state and spreads; it is not a permanent defect or derived localized ground state.",
                "Tree vacuum and quadratic G1 masses are unchanged at zero VEV; loop corrections, naturalness and stabilization remain open.",
                "This is not a microscopic derivation of G1's paired source probes, the cited paper's full claims, or a Standard Model spectrum."],
        checks=checks, summary=dict(checks=len(checks), passed=sum(row["passed"] for row in checks),
                                    failed=[row["name"] for row in checks if not row["passed"]]),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    write_report(result)
    print(json.dumps(result["summary"], indent=2))
    print(f"target l=0: sigma={main_example['sigma']:.12g}, K={main_example['rate_coefficient']:.12g}")
    print(f"packet: K_target={packet_target:.12g}, K_total={packet_total:.12g}, fraction={packet_target / packet_total:.12g}")
    if result["summary"]["failed"]:
        raise SystemExit(1)


def write_report(result):
    ex = result["main_example"]
    packet = result["packet"]
    summary = result["summary"]
    lines = [
        "# Route G G2 — local dynamical mode conversion", "",
        "A bounded tree-level engineering calculation, not an empirical particle model. No measured masses or fitted parameters enter.", "",
        "## Shared action and interaction", "",
        "Add one positive-kinetic complex bulk detector field X, neutral under Phi's external U(1), and the nonnegative potential +kappa5 |Phi|^2 |X|^2. Both VEVs are zero, so the G1 quadratic masses remain unchanged at tree level. The local 5D interaction is not a pinned defect or an arbitrarily selected projector.", "",
        r"$$m_n^2=M_5^2+(n+\alpha)^2/R^2,\quad \mu_l^2=M_D^2+l^2/R^2,\quad g=\kappa_5/(2\pi R).$$", "",
        r"$$\Phi_n+X_l\longrightarrow\Phi_m+X_j,\qquad \mathcal M=-g\,\delta_{n+l,m+j}.$$", "",
        "The action conserves both complex fields' net U(1) charges, four-momentum, and total compact momentum. This channel has one particle of each species before and after; total particle-plus-antiparticle number is not an exact symmetry of the full field theory. Changing Phi's internal mode is balanced by detector recoil.", "",
        "One declared card: R=1, M5=0.5, alpha=0.25, MD=5, g=0.05, initial COM momentum p=0.2; n=1 and target m=0. All quantities use a common arbitrary mass unit. kappa5 has dimension -1, g is dimensionless, sigma and K have dimension -2.", "",
        r"$$p_f=\frac{\sqrt{[s-(m_m+\mu_j)^2][s-(m_m-\mu_j)^2]}}{2\sqrt{s}},\quad \sigma=\frac{g^2}{16\pi s}\frac{p_f}{p_i},\quad K=\sigma\left(\frac{p_i}{E_\Phi}+\frac{p_i}{E_X}\right).$$", "",
        "Set sigma=K=0 when the final threshold is closed. These expressions are for distinct scalars, with no identical-particle factor.", "",
        "## Explicit nonzero process", "",
        "| Quantity | Engineering value |", "|---|---:|",
    ]
    for label, key in (("initial Phi mass m1", "initial_phi_mass"), ("final Phi mass m0", "final_phi_mass"),
                       ("initial X mass mu0", "initial_detector_mass"), ("final X mass mu1", "final_detector_mass"),
                       ("sqrt(s)", "sqrt_s"), ("final momentum", "outgoing_momentum"),
                       ("sigma", "sigma"), ("K = sigma v_rel", "rate_coefficient"),
                       ("Phi energy change", "phi_energy_change"), ("detector energy change", "detector_energy_change")):
        lines.append(f"| {label} | {ex[key]:.12g} |")
    reverse = result["closed_reverse_at_same_small_momentum"]
    lines += [
        "", "This is genuine scattering from one Phi tree-level mass eigenstate to another; it is not a passive basis change. The same tree spectrum supplies both states. Interacting loop pole shifts are not computed. Phi's lost energy is detector recoil/energy, not missing energy.", "",
        f"The inverse channel Phi0 + X1 -> Phi1 + X0 at the *same small p=0.2* has sqrt(s)={reverse['sqrt_s']:.12g}, below threshold {reverse['final_threshold']:.12g}, so it is closed. At the same s as the forward reaction, inverse detailed balance is satisfied. These are different energy preparations, not a violation of reversibility.", "",
        f"For l=0 the analytic rest-threshold condition is MD > [((n-m)/R)^2-(mn-mm)^2]/[2(mn-mm)] = {result['heavy_recoil_threshold']['analytic_critical_MD']:.12g}; the chosen MD=5 is well above it. Direct threshold checks at 0.9, 1.0 and 1.1 times this critical mass agree. Heavy recoil costs sqrt(MD^2+((n-m)/R)^2)-MD, which approaches (n-m)^2/(2 MD R^2), not the massless emitted-mode cost |n-m|/R.", "",
        f"As p_initial approaches zero, the exothermic cross section scales as 1/p_initial, but sigma*v_rel has the finite limit {result['exothermic_small_momentum']['analytic_finite_sigma_v_limit']:.12g}. The smallest tested momentum 0.0002 agrees to relative error {result['exothermic_small_momentum']['rows'][-1]['relative_error_to_finite_limit']:.6g}. A cross section is not a probability.", "",
        "## A finite localized dynamical detector state", "",
        r"$$c_l=Z^{-1/2}\exp[-l^2/(4w^2)-il\theta_0],\quad l=-4,\ldots,4,\quad w=1.2,\quad\theta_0=0.7.$$", "",
        "This is an exactly normalized finite prepared state; no approximation to an infinite Gaussian is claimed. Each component has the same 3-momentum magnitude but different s_l. The free packet spreads. Its canonical angular one-particle wavefunction is not a covariant relativistic position density.", "",
        "| time | norm | mean detector energy | magnitude of first circular moment | mean angle |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in packet["free_spreading"]:
        lines.append(f"| {row['time']:g} | {row['norm']:.12g} | {row['mean_energy']:.12g} | {row['first_circular_moment_magnitude']:.12g} | {row['circular_mean_angle']:.12g} |")
    lines += [
        "", "For the fixed target m=0, recoil j=l+1 distinguishes incoming components. After tracing recoil, the inclusive rate kernel is diagonal in l. A pure localized packet and its momentum-dephased mixture have equal inclusive target rates. The illustrative rate matrices retain different off-diagonal entries, but do not determine physical recoil coherence after unobserved momenta or timing are traced. Translating this packet cannot change that rate against a delocalized incident Phi_n mode.", "",
        r"$$K_{\rm packet}=\sum_l|c_l|^2K_l,\qquad A_{jl}=\sqrt{K_l}\,\delta_{j,l+1},\qquad \mathrm{tr}(A\rho A^\dagger)=\sum_lK_l\rho_{ll}.$$", "",
        "This is a **common-overlap dilute rate coefficient**, not a universal wavepacket cross section, scattering probability, or a G1 source-visibility weight. The rate map A is not a single-energy S-matrix: the components have different s_l. Conditional incoming and outgoing compact momenta must both be weighted by the target event rate.", "",
        "| initial l | probability weight | sqrt(s_l) | target K_l (m=0) | all open final m | all-channel K_l |",
        "|---:|---:|---:|---:|---|---:|",
    ]
    for row in packet["components"]:
        opened = ", ".join(str(r["m"]) for r in row["all_channels"]["open_channels"])
        lines.append(f"| {row['l']} | {row['weight']:.8g} | {row['target']['sqrt_s']:.10g} | {row['target']['rate_coefficient']:.10g} | {opened} | {row['total_rate_coefficient']:.10g} |")
    recoil = result["recoil_rate_map"]
    lines += [
        "", f"Weighted target coefficient: **{packet['target_rate_coefficient']:.12g}**. All-channel coefficient: **{packet['total_rate_coefficient']:.12g}**. Target fraction of all scattering rates: **{packet['target_fraction_of_all_scattering_rates']:.12g}**. These fractions depend on this prepared state and common-overlap prescription; the chosen m=0 channel is not exclusive. The total includes elastic scattering m=n, not the probability of remaining in that KK sector: the latter also includes the unscattered identity amplitude.", "",
        f"Conditional mean incoming detector l={recoil['rate_conditioned_incoming_l_mean']:.12g}; conditional mean outgoing j={recoil['rate_conditioned_outgoing_j_mean']:.12g}; shift={recoil['conditional_shift']:.12g}=n-m. Comparing outgoing mean with the *unconditioned* incident mean would incorrectly mix recoil with event-selection bias.", "",
        "Every open channel is enumerated by the rigorous finite candidate bound |m+alpha| < R sqrt(s_l), followed by the exact two-body threshold. JSON also records rejected candidates and every open channel's energies, cross section, independent radial phase-space check, and inverse rate.", "",
        "## Verification and scope", "",
        f"**{summary['passed']}/{summary['checks']} numerical/code checks passed.** This certifies this implementation, not the truth of the model as a particle theory.", "",
        "Checks cover normalized circle overlap and forbidden channels; independent radial-root/phase-space versus Kallen kinematics; complete open-channel enumeration; energy and compact-momentum conservation; inverse detailed balance at common s; zero coupling; large-gauge relabelling; exact finite-packet normalization, free spreading and energy; inclusive recoil trace, dephasing, translation and conditional momentum shift.", "",
        "The positivity statement is analytic: all three potential coefficients are nonnegative. Sampled potential values are regression checks, not a numerical proof for all field amplitudes.", "",
        "Limits: 5D EFT and tree amplitudes only; the EFT cutoff must exceed the prepared scattering energies, and no loop matching/naturalness assertion is made. Radius/holonomy stabilization, localized-state preparation, spin/chirality, three families, the full CERN-paper claims, and a Route-F or Standard Model bridge remain open. This does not derive the G1 paired probe operator from the new detector field.", "",
        "Reproduce: python3 route_g/code/verify_g2_local_conversion.py. No plot, lattice scan, per-particle fitting, or external work source is used. JSON retains the source SHA-256 and every check.", "",
        "Normalization references: [PDG Kinematics, equations 49.27--49.33](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-kinematics.pdf); [Tong QFT, sections 3.4 and 3.6](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S3.html). Full action and packet derivation: ../tex/route_g_local_conversion.tex.", "",
    ]
    OUT.with_suffix(".md").write_text("\n".join(lines))


if __name__ == "__main__":
    run()
