#!/usr/bin/env python3
"""Route G G1: bounded free-circle spectrum and fixed-source-response checks.

No measured particle masses, Route-F benchmark, fit, or per-mode projector.
The Gaussian toy model has no dynamical gauge field or physical S-matrix.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g1_circle_spectrum"
CARDS = (
    dict(id="circle_trivial", R=1.0, M5=0.5, alpha=0.0),
    dict(id="circle_quarter", R=1.0, M5=0.5, alpha=0.25),
    dict(id="circle_half", R=1.0, M5=0.5, alpha=0.5),
    dict(id="circle_larger", R=2.0, M5=0.5, alpha=0.25),
)
LAMBDA = 1.0


def cmatrix(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def modes(card, integers):
    q = np.asarray(integers, dtype=int) + card["alpha"]
    mass_squared = card["M5"] ** 2 + q ** 2 / card["R"] ** 2
    phase = np.exp(1j * np.pi * q)
    coupling = LAMBDA / np.sqrt(4 * np.pi * card["R"]) * np.stack((1 + phase, 1 - phase))
    return q, mass_squared, coupling


def response_truncated(card, p_squared, cutoff):
    _, mass_squared, coupling = modes(card, np.arange(-cutoff, cutoff + 1))
    return (coupling / (mass_squared + p_squared)) @ coupling.conj().T


def response_exact(card, p_squared):
    """Independent infinite-circle image sum; beta must be positive.

    S0 = sum_n 1/((n+alpha)^2+beta^2).
    Spi = sum_n exp(i*pi*(n+alpha))/((n+alpha)^2+beta^2).
    With a=exp(-pi*beta), t=2pi*alpha and D=1+a^4-2a^2*cos(t):
    S0=(pi/beta)*(1-a^4)/D;
    Spi=(pi/beta)*a*(1-a^2)*(1+exp(i*t))/D.
    """
    beta = card["R"] * np.sqrt(card["M5"] ** 2 + p_squared)
    if beta <= 0:
        raise ValueError("Image-sum implementation requires beta > 0")
    a = np.exp(-np.pi * beta)
    angle = 2 * np.pi * card["alpha"]
    denominator = 1 + a ** 4 - 2 * a ** 2 * np.cos(angle)
    s0 = np.pi / beta * (1 - a ** 4) / denominator
    spi = np.pi / beta * a * (1 - a ** 2) * (1 + np.exp(1j * angle)) / denominator
    return LAMBDA ** 2 * card["R"] / (2 * np.pi) * np.array([
        [s0 + spi.real, 1j * spi.imag],
        [-1j * spi.imag, s0 - spi.real],
    ])


def tail_bound(card, cutoff):
    """Rigorous absolute bound for every response entry, for p_E^2>=0.

    |g_a,n g_b,n*| <= lambda^2/(pi R);
    m_n^2+p_E^2 >= (|n|-|alpha|)^2/R^2.
    Bound the two omitted positive tails by integrals from N to infinity.
    """
    if cutoff <= abs(card["alpha"]):
        raise ValueError("Tail bound requires N > |alpha|")
    return 2 * LAMBDA ** 2 * card["R"] / (np.pi * (cutoff - abs(card["alpha"])))


def lattice(card, sites):
    if sites % 2:
        raise ValueError("An antipodal probe requires an even site count")
    spacing = 2 * np.pi * card["R"] / sites
    transport = np.zeros((sites, sites), dtype=complex)
    for j in range(sites):
        transport[j, (j + 1) % sites] = np.exp(2j * np.pi * card["alpha"] / sites)
    mass = ((2 * np.eye(sites) - transport - transport.conj().T) / spacing ** 2
            + card["M5"] ** 2 * np.eye(sites))
    return mass, transport, spacing


def lattice_probes(transport, spacing):
    sites = transport.shape[0]
    wilson = np.prod([transport[j, j + 1] for j in range(sites // 2)])
    probes = np.zeros((2, sites), dtype=complex)
    probes[:, 0] = LAMBDA / np.sqrt(2 * spacing)
    probes[:, sites // 2] = np.array([1, -1]) * LAMBDA * wilson / np.sqrt(2 * spacing)
    return probes, wilson


def run():
    checks = []

    def check(name, passed, value=None, tolerance=None):
        row = dict(name=name, passed=bool(passed))
        if value is not None:
            row["measured"] = float(value)
        if tolerance is not None:
            row["tolerance"] = float(tolerance)
        checks.append(row)

    cards, finite_difference, response_convergence = [], [], []
    for card in CARDS:
        integers = np.arange(-3, 4)
        q, mass_squared, coupling = modes(card, integers)
        rows = []
        for k, n in enumerate(integers):
            residue = np.outer(coupling[:, k], coupling[:, k].conj())
            weights = LAMBDA ** 2 / (2 * np.pi * card["R"]) * (
                1 + np.array([1, -1]) * np.cos(np.pi * q[k]))
            mineig = np.linalg.eigvalsh(residue).min()
            error = np.max(np.abs(residue.diagonal().real - weights))
            total_error = abs(np.trace(residue).real - LAMBDA ** 2 / (np.pi * card["R"]))
            check(f"{card['id']}:n={n}:PSD_residue", mineig >= -1e-14, mineig, 1e-14)
            check(f"{card['id']}:n={n}:weight_formula", error < 1e-14, error, 1e-14)
            check(f"{card['id']}:n={n}:weight_sum_rule", total_error < 1e-14, total_error, 1e-14)
            rows.append(dict(n=int(n), q=float(q[k]), mass_squared=float(mass_squared[k]),
                             mass=float(np.sqrt(mass_squared[k])), coupling=cmatrix(coupling[:, k]),
                             residue=cmatrix(residue), plus_weight=float(weights[0]),
                             minus_weight=float(weights[1]),
                             normalized_relative_probe_weights=(weights * np.pi * card["R"] / LAMBDA ** 2).tolist()))
        second_difference_error = np.max(np.abs(np.diff(mass_squared, n=2) - 2 / card["R"] ** 2))
        check(f"{card['id']}:signed_branch_second_difference",
              second_difference_error < 1e-12, second_difference_error, 1e-12)
        _, wide_mass, _ = modes(card, np.arange(-32, 33))
        cards.append(dict(**card, lambda_common=LAMBDA, modes=rows,
                          lowest_eight_mass_squared=np.sort(wide_mass)[:8].tolist(),
                          signed_branch_second_difference=2 / card["R"] ** 2))
        previous_error = None
        for sites in (32, 64, 128):
            mass, transport, spacing = lattice(card, sites)
            numerical = np.linalg.eigvalsh(mass)[:8]
            integers = np.arange(-sites // 2, sites // 2)
            expected = np.sort(modes(card, integers)[1])[:8]
            error = float(np.max(np.abs(numerical - expected)))
            discrete = np.sort(card["M5"] ** 2 + 4 / spacing ** 2 * np.sin(
                np.pi * (integers + card["alpha"]) / sites) ** 2)[:8]
            independent_error = float(np.max(np.abs(numerical - discrete)))
            check(f"{card['id']}:N={sites}:independent_lattice_diagonalization",
                  independent_error < 5e-11, independent_error, 5e-11)
            ratio = None if previous_error is None else previous_error / error
            if ratio is not None:
                check(f"{card['id']}:N={sites}:second_order_convergence", 3.75 < ratio < 4.1, ratio)
            probe, wilson = lattice_probes(transport, spacing)
            selected = np.arange(-3, 4)
            fourier = np.exp(2j * np.pi * np.outer(np.arange(sites), selected) / sites) / np.sqrt(sites)
            probe_error = float(np.max(np.abs(probe @ fourier - modes(card, selected)[2])))
            check(f"{card['id']}:N={sites}:canonical_probe_normalization",
                  probe_error < 1e-13, probe_error, 1e-13)
            finite_difference.append(dict(
                card_id=card["id"], sites=sites, physical_spacing=spacing,
                lowest_eight_eigenvalues=numerical.tolist(),
                analytic_continuum_eigenvalues=expected.tolist(), max_absolute_error=error,
                previous_error_over_current=ratio, discrete_symbol_error=independent_error,
                canonical_probe_error=probe_error, path_wilson_line=cmatrix(wilson)))
            previous_error = error
        for p_squared in (0.0, 0.7, 4.0):
            exact = response_exact(card, p_squared)
            previous, errors = None, []
            for cutoff in (16, 32, 64, 128, 256):
                response = response_truncated(card, p_squared, cutoff)
                bound = tail_bound(card, cutoff)
                error = float(np.max(np.abs(response - exact)))
                smallest = float(np.linalg.eigvalsh(response).min())
                remainder_smallest = float(np.linalg.eigvalsh(exact - response).min())
                prefix = f"{card['id']}:p2={p_squared}:cutoff={cutoff}:"
                check(prefix + "rigorous_tail_bound", error <= bound + 1e-13, error / bound, 1)
                check(prefix + "PSD_response", smallest >= -1e-13, smallest, 1e-13)
                check(prefix + "PSD_omitted_response", remainder_smallest >= -1e-13, remainder_smallest, 1e-13)
                response_convergence.append(dict(
                    card_id=card["id"], p_squared=p_squared, cutoff_abs_n=cutoff,
                    response=cmatrix(response), image_sum_exact_response=cmatrix(exact),
                    max_entry_error_against_image_sum=error,
                    change_from_previous_cutoff=None if previous is None else float(np.max(np.abs(response - previous))),
                    rigorous_per_entry_tail_bound=bound,
                    smallest_response_eigenvalue=smallest,
                    smallest_omitted_response_eigenvalue=remainder_smallest))
                errors.append(error)
                previous = response
            check(f"{card['id']}:p2={p_squared}:decreasing_cutoff_error",
                  all(x > y for x, y in zip(errors, errors[1:])))

    integers = np.arange(-8, 9)
    _, _, coupling = modes(CARDS[0], integers)
    forbidden = max(np.max(np.abs(coupling[0, integers % 2 != 0])),
                    np.max(np.abs(coupling[1, integers % 2 == 0])))
    check("alpha_zero:even_odd_selection_rule", forbidden < 1e-14, forbidden, 1e-14)
    rng = np.random.default_rng(20260922)
    card = CARDS[1]
    q0, mass0, coupling0 = modes(card, integers)
    shifted = dict(card, alpha=card["alpha"] + 1)
    q1, mass1, coupling1 = modes(shifted, integers - 1)
    large_gauge = dict(
        alpha_original=card["alpha"], alpha_transformed=shifted["alpha"],
        original_n_window=[-8, 8], mapped_n_window=[-9, 7],
        q_error=float(np.max(np.abs(q1 - q0))),
        mass_squared_error=float(np.max(np.abs(mass1 - mass0))),
        probe_coupling_error=float(np.max(np.abs(coupling1 - coupling0))),
        note="Relabel n->n-1 with alpha->alpha+1; never compare identical hard-cutoff windows.")
    check("large_gauge:matched_window_spectrum_and_probes",
          max(large_gauge[k] for k in ("q_error", "mass_squared_error", "probe_coupling_error")) < 1e-13)

    sites = 24
    mass, transport, spacing = lattice(card, sites)
    chi = rng.normal(size=sites)
    gauge = np.diag(np.exp(1j * chi))
    changed_transport = gauge @ transport @ gauge.conj().T
    changed_mass = ((2 * np.eye(sites) - changed_transport - changed_transport.conj().T)
                    / spacing ** 2 + card["M5"] ** 2 * np.eye(sites))
    probes, wilson = lattice_probes(transport, spacing)
    changed_probes, changed_wilson = lattice_probes(changed_transport, spacing)
    original_response = probes @ np.linalg.solve(mass + 0.7 * np.eye(sites), probes.conj().T)
    changed_response = changed_probes @ np.linalg.solve(changed_mass + 0.7 * np.eye(sites), changed_probes.conj().T)
    local_gauge = dict(
        sites=sites, random_seed=20260922,
        operator_covariance_error=float(np.max(np.abs(changed_mass - gauge @ mass @ gauge.conj().T))),
        probe_covariance_error=float(np.max(np.abs(changed_probes - np.exp(1j * chi[0]) * probes @ gauge.conj().T))),
        wilson_covariance_error=float(abs(changed_wilson - np.exp(1j * (chi[0] - chi[sites // 2])) * wilson)),
        spectrum_error=float(np.max(np.abs(np.linalg.eigvalsh(mass) - np.linalg.eigvalsh(changed_mass)))),
        response_error=float(np.max(np.abs(original_response - changed_response))))
    for key, value in local_gauge.items():
        if key.endswith("_error"):
            check("local_U1:" + key, value < 1e-11, value, 1e-11)

    _, mass_squared, probes = modes(card, np.arange(-4, 5))
    kinetic = np.diag(mass_squared + 0.7)
    original = probes @ np.linalg.solve(kinetic, probes.conj().T)
    basis_changes = []
    for label in ("real_orthogonal", "complex_unitary"):
        matrix = rng.normal(size=(9, 9))
        if label == "complex_unitary":
            matrix = matrix + 1j * rng.normal(size=(9, 9))
        basis, _ = np.linalg.qr(matrix)
        changed_kinetic = basis.conj().T @ kinetic @ basis
        changed_probes = probes @ basis
        changed = changed_probes @ np.linalg.solve(changed_kinetic, changed_probes.conj().T)
        spectral_error = float(np.max(np.abs(np.linalg.eigvalsh(changed_kinetic) - np.linalg.eigvalsh(kinetic))))
        response_error = float(np.max(np.abs(changed - original)))
        check(label + ":spectrum", spectral_error < 1e-11, spectral_error, 1e-11)
        check(label + ":response", response_error < 1e-11, response_error, 1e-11)
        basis_changes.append(dict(type=label, spectrum_error=spectral_error, response_error=response_error))

    # Degenerate eigenvectors have no unique individual residue: sum the projector.
    degenerate_tests = []
    for alpha, selected in ((0.0, [-1, 1]), (0.5, [-1, 0])):
        degenerate_card = dict(card, alpha=alpha)
        _, masses, coupling = modes(degenerate_card, selected)
        unitary, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        total_residue = coupling @ coupling.conj().T
        changed_coupling = coupling @ unitary
        changed_residue = changed_coupling @ changed_coupling.conj().T
        error = float(np.max(np.abs(total_residue - changed_residue)))
        singular_values = np.linalg.svd(coupling, compute_uv=False)
        rank = int(np.sum(singular_values > 1e-12))
        expected_rank = 1 if alpha == 0 else 2
        check(f"degenerate_alpha={alpha}:summed_projector_residue", error < 1e-13, error, 1e-13)
        check(f"degenerate_alpha={alpha}:observable_coupling_rank", rank == expected_rank, rank)
        if alpha == 0:
            dark_norm = float(np.linalg.norm(coupling @ np.array([1, -1]) / np.sqrt(2)))
            check("alpha_zero:dark_sine_standing_wave", dark_norm < 1e-13, dark_norm, 1e-13)
        degenerate_tests.append(dict(alpha=alpha, mode_pair=selected, common_mass_squared=float(masses[0]),
                                     summed_residue=cmatrix(total_residue), summed_residue_basis_error=error,
                                     coupling_singular_values=singular_values.tolist(),
                                     observable_coupling_rank=rank, dark_subspace_dimension=2-rank,
                                     note="Positive branch residues do not forbid dark combinations inside a degenerate eigenspace."))
    changes = []
    for alternative, label in ((CARDS[0], "physical_holonomy_change"),
                               (CARDS[2], "physical_holonomy_change"),
                               (CARDS[3], "physical_radius_change")):
        ref = np.sort(modes(card, np.arange(-16, 17))[1])[:8]
        alt = np.sort(modes(alternative, np.arange(-16, 17))[1])[:8]
        change = float(np.max(np.abs(ref - alt)))
        check(alternative["id"] + ":physical_background_changes_spectrum", change > 1e-4, change)
        changes.append(dict(type=label, from_card=card["id"], to_card=alternative["id"],
                            max_low_spectrum_change=change,
                            max_response_change_at_p_squared_0_7=float(np.max(np.abs(
                                response_exact(card, 0.7) - response_exact(alternative, 0.7))))))
    result = dict(
        project="Route G", stage="G1 free-circle fixed-source-response prototype",
        action=dict(
            spacetime="Fixed M4 x S1; theta in [0,2pi); eta=(+---); ds^2=eta dx dx-R^2 dtheta^2",
            field="Periodic charged complex scalar, M5^2>=0, no chiral fermions",
            covariant_derivative="Dtheta=partial_theta+i alpha",
            lagrangian="R*(|partial_mu Phi|^2-R^(-2)|Dtheta Phi|^2-M5^2|Phi|^2)",
            holonomy="exp(2 pi i alpha), flat nondynamical U1 background",
            mode_normalization="Phi=sum_n phi_n exp(i n theta)/sqrt(2 pi R)",
            fixed_probes="O_plus/minus=[Phi(0)+/-exp(i pi alpha)Phi(pi)]/sqrt(2); interaction=lambda*J_a^dagger*O_a+h.c.",
            probe_path="Same oriented half-circle 0 to pi for both channels and every mode",
            source_charge="O_a and J_a transform with the same charge at theta=0; J_a^dagger O_a is invariant",
            units="One arbitrary reference mass unit; [R]=-1, [Phi]=3/2, [phi_n]=1, [J_a]=5/2; lambda=1 is dimensionless",
            radius_notation="R is a compactification radius, not the cited paper's R dimension label",
            phase_hypothesis="alpha is a declared holonomy; no proved mapping to Ei, Ep, En exists"),
        cards=cards, finite_difference=finite_difference,
        invariance_checks=dict(large_gauge=large_gauge, local_gauge=local_gauge,
                               passive_basis_changes=basis_changes, degenerate_subspaces=degenerate_tests),
        physical_background_changes=changes, response_convergence=response_convergence,
        limitations=[
            "No measured masses, flavor, unification, seesaw, or proton data are inputs.",
            "Negative mode number is not negative energy; coordinate changes do not change mass.",
            "Relative probe weights are not probabilities or decay branching fractions.",
            "Each Fourier branch has nonzero combined weight, but alpha=0 degenerate sine combinations are dark to both antipodal probes.",
            "The two nonlocal external probes are prescribed, not uniquely derived local detectors.",
            "Fixed-source response is not a physical S-matrix with a dynamical gauge sector.",
            "R and alpha are external backgrounds; no stabilization or real-time transition has been derived.",
            "The numerical KK cutoff is not an interacting high-dimensional UV regulator.",
            "This scalar tower has no SM content, three chiral families, or Route-F Yukawa derivation.",
            "Positive Gaussian residues do not validate an unspecified enlarged theory."],
        references=[
            "https://arxiv.org/abs/0908.0333, section 8 (circle compactification)",
            "https://doi.org/10.1093/ptep/ptab085, section 2 (Wilson-line mass shifts)"],
        source_sha256={"code/verify_g1_circle_spectrum.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        checks=checks, summary=dict(checks=len(checks), passed=sum(c["passed"] for c in checks),
                                    all_pass=all(c["passed"] for c in checks),
                                    physical_model_validated=False, particle_spectrum_fitted=False))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    write_report(result)
    print(json.dumps(result["summary"], indent=2))
    if not result["summary"]["all_pass"]:
        raise SystemExit("Failed checks: " + ", ".join(c["name"] for c in checks if not c["passed"]))


def write_report(result):
    lines = [
        "# Route G G1: common circle spectrum and fixed probe channels", "",
        "This is a bounded free Gaussian-model verification, not a particle fit. No measured masses or Route-F benchmark values were used.", "",
        "## Declared model", "",
        r"For a periodic complex scalar on fixed $M_4\times S^1$, $D_\theta=\partial_\theta+i\alpha$ and $\Phi=\sum_n\phi_n e^{in\theta}/\sqrt{2\pi R}$ give $m_n^2=M_5^2+(n+\alpha)^2/R^2$.", "",
        r"The common source operators are $O_\pm=[\Phi(0)\pm W(0,\pi)\Phi(\pi)]/\sqrt2$, $W=e^{i\pi\alpha}$. The interaction $\sum_a\lambda J_a^\dagger O_a+\mathrm{h.c.}$ gives $g_{\pm,n}=\lambda[1\pm e^{i\pi(n+\alpha)}]/\sqrt{4\pi R}$.", "",
        r"Every residue $\rho_n=g_ng_n^\dagger$ is positive semidefinite. Its diagonal obeys $|g_+|^2+|g_-|^2=\lambda^2/(\pi R)$ for each chosen Fourier momentum branch. Normalized relative probe weights are $(1\pm\cos\pi(n+\alpha))/2$, not probabilities or branching fractions.", "",
        r"At $\alpha=0$, these same probes select even/odd mode numbers. The signed-branch second difference is $m_{n+1}^2-2m_n^2+m_{n-1}^2=2/R^2$; this is not a statement about consecutively sorted degenerate levels.", "",
        r"The source $J_a$ and charged, gauge-covariant operator $O_a$ transform with the same charge at theta=0, so $J_a^\dagger O_a$ is invariant. With dimensionless lambda, $[J_a]=5/2$. This is a fixed-background source response, not a dynamical gauge-theory S-matrix. The paired nonlocal probes have not been derived from a local detector theory.", "",
        "Here R is compactification radius, not the cited paper's dimension label R; alpha is a modeling hypothesis, not a proved identification with Ei, Ep, or En.", "",
        "## Four predeclared engineering cards", "",
        "One arbitrary common mass unit is used; R is inverse mass. M5=0.5 and common source normalization lambda=1 are fixed.", "",
        "| Card | R | alpha | Lowest four mass-squared values (multiplicity retained) |",
        "|---|---:|---:|---|",
    ]
    for card in result["cards"]:
        entries = ", ".join(f"{x:.8g}" for x in card["lowest_eight_mass_squared"][:4])
        lines.append(f"| {card['id']} | {card['R']:g} | {card['alpha']:g} | {entries} |")
    lines += [
        "", "These are background cards, not measured particles or a dynamical transition. Negative mode number is not negative energy.", "",
        "## Independent numerical checks", "",
        "Link-covariant finite-difference Laplacians are built and diagonalized at 32, 64, and 128 sites. Arbitrary site phases transform both links and the same physical probe Wilson line. Passive real and complex mode-basis changes transform both the kernel and the probe map. Under alpha -> alpha+1, the mode and cutoff window are relabelled n -> n-1.", "",
        "An arbitrary eigenbasis within a degenerate eigenspace is not fixed by mass alone: unitary basis changes preserve its summed residue/projector. Particular states can still be distinguished by other operators or preparation. The residue sums are explicitly tested at alpha=0 and alpha=1/2.", "",
        r"**Important visibility qualification:** at $\alpha=0$, the columns for $n=\pm k$ are equal, so their sine standing-wave combination is dark to both antipodal probes. The coupling map on that two-dimensional eigenspace has rank one, although its pole is visible through the other combination. At $\alpha=1/2$, the tested degenerate pair has rank two. Therefore the branch-weight sum rule is not a proof that every degenerate state is visible. The invariant diagnostic is $\dim\ker G_E=\dim E-\operatorname{rank}G_E$.", "",
        "| Card | Low-spectrum error: 32 sites | 64 sites | 128 sites | Last error ratio |",
        "|---|---:|---:|---:|---:|",
    ]
    for card in result["cards"]:
        rows = [r for r in result["finite_difference"] if r["card_id"] == card["id"]]
        lines.append(f"| {card['id']} | {rows[0]['max_absolute_error']:.6g} | {rows[1]['max_absolute_error']:.6g} | {rows[2]['max_absolute_error']:.6g} | {rows[2]['previous_error_over_current']:.6g} |")
    lines += [
        "", "An error ratio approaching four is second-order discretization convergence, not a particle-physics uncertainty.", "",
        "## Infinite image sum and rigorous truncation bound", "",
        r"For $p_E^2\ge0$, $G_{ab}=\sum_n g_{a,n}g_{b,n}^*/(p_E^2+m_n^2)$. Set $\beta=R\sqrt{M_5^2+p_E^2}>0$, $a=e^{-\pi\beta}$, $t=2\pi\alpha$, $D=1+a^4-2a^2\cos t$. Independent image sums give", "",
        r"$$S_0=\frac\pi\beta\frac{1-a^4}{D},\qquad S_\pi=\frac\pi\beta\frac{a(1-a^2)(1+e^{it})}{D},$$",
        r"$$G=\frac{\lambda^2R}{2\pi}\begin{pmatrix}S_0+\Re S_\pi&i\Im S_\pi\\-i\Im S_\pi&S_0-\Re S_\pi\end{pmatrix}.$$", "",
        r"Here $S_0=\sum_n[(n+\alpha)^2+\beta^2]^{-1}$ and $S_\pi=\sum_n e^{i\pi(n+\alpha)}[(n+\alpha)^2+\beta^2]^{-1}$.", "",
        "The script compares this expression with finite mode sums at cutoffs 16,32,64,128,256 and p_E^2=0,0.7,4 for every card.", "",
        r"Since $|g_{a,n}g_{b,n}^*|\le\lambda^2/(\pi R)$ and $p_E^2+m_n^2\ge(|n|-|\alpha|)^2/R^2$, every entry for $N>|\alpha|$ obeys", "",
        r"$$|G_{ab}-G^{(N)}_{ab}|\le\frac{2\lambda^2R}{\pi}\sum_{n=N+1}^\infty\frac1{(n-|\alpha|)^2}\le\frac{2\lambda^2R}{\pi(N-|\alpha|)}.$$", "",
        "This is a rigorous bound for the declared free response, not a statistical interval or interacting UV regulator.", "",
        "| Card | p_E^2 | Largest entry error at cutoff 256 | Rigorous bound |",
        "|---|---:|---:|---:|",
    ]
    for row in result["response_convergence"]:
        if row["cutoff_abs_n"] == 256:
            lines.append(f"| {row['card_id']} | {row['p_squared']:g} | {row['max_entry_error_against_image_sum']:.6g} | {row['rigorous_per_entry_tail_bound']:.6g} |")
    summary = result["summary"]
    lines += [
        "", "## Result and boundary of the claim", "",
        f"**{summary['passed']}/{summary['checks']} numerical/code checks passed.** This does not validate a physical theory or establish an empirical mass spectrum.", "",
        "Passive coordinate/field-basis changes leave poles and complete probe response invariant. Physically inequivalent holonomy or radius cards change poles. Fixed physical probes can weight modes differently without moving a pole.", "",
        "Still absent: dynamical stabilization, local detector realization, transition amplitudes, spin/chirality/three families, a Standard Model spectrum, and a derivation of Route-F Yukawa matrices. None is implicitly promoted.", "",
        "Reproduce: python3 route_g/code/verify_g1_circle_spectrum.py. JSON retains complex residues, responses, every check, seed, and source hash.", "",
        "References: [Tong, section 8](https://arxiv.org/abs/0908.0333); [Yamada, section 2](https://doi.org/10.1093/ptep/ptab085).", "",
    ]
    OUT.with_suffix(".md").write_text("\n".join(lines))


if __name__ == "__main__":
    run()
