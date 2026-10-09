#!/usr/bin/env python3
"""All-channel contact readout, exact stationary dark clocks and preparation limit.

No fits. Phase effects, full spectral projectors and conditional states retain
all nine modes. A uniform unread-update prior is an explicit assumption.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from verify_gr9_event_network import F, laplacian, spectrum

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr10_relational_readout"
E9 = np.eye(9, dtype=complex)
N = np.diag(np.arange(3))
DPLUS = (E9[:, 2]-2*E9[:, 4]+E9[:, 6])/np.sqrt(6)
DMINUS = (E9[:, 2]-E9[:, 6])/np.sqrt(2)


def proj(v):
    return np.outer(v, v.conj())


def word(phi):
    return np.exp(1j*np.arange(3)*phi)


def conditional(rho, phi, m):
    w = word(phi)
    r = np.einsum("a,abcd,c->bd", w.conj(), rho.reshape(3, 3, 3, 3), w)/m
    p = float(np.trace(r).real)
    if p <= 0:
        raise ValueError("Undefined conditional state at zero probability")
    return r/p, p


def joint(rho, phases):
    m = len(phases)
    return np.array([[np.vdot(np.kron(word(a), word(b)),
                             rho@np.kron(word(a), word(b))).real/m**2
                      for b in phases] for a in phases])


def marginal_b(rho):
    return np.einsum("abad->bd", rho.reshape(3, 3, 3, 3))


def trace_distance(a, b):
    return float(np.linalg.svd(a-b, compute_uv=False).sum()/2)


def readout_metrics(rho, phases):
    reference, _ = conditional(rho, 0, len(phases))
    residuals, likelihoods = [], []
    for phi in phases:
        r, p = conditional(rho, phi, len(phases))
        u = np.diag(word(phi))
        residuals.append(trace_distance(r, u@reference@u.conj().T))
        likelihoods.append(len(phases)*p)
    p = joint(rho, phases)
    delta = phases[None, :]-phases[:, None]
    return {"seven_reading_rigidity_residual": max(residuals),
            "scaled_A_likelihood_range": max(likelihoods)-min(likelihoods),
            "relative_first_harmonic_real": float(np.sum(p*np.cos(delta))),
            "relative_second_harmonic_real": float(np.sum(p*np.cos(2*delta)))}


def energy_projectors(h):
    values, vectors = np.linalg.eigh(h)
    groups = []
    for index, value in enumerate(values):
        if not groups or abs(value-values[groups[-1][0]]) > 1e-10:
            groups.append([index])
        else:
            groups[-1].append(index)
    return values, [vectors[:, ix]@vectors[:, ix].conj().T for ix in groups]


def dephase_energy(a, projectors):
    return sum(p@a@p for p in projectors)


def run():
    card = json.loads((ROOT/"data/RELATIONAL_READOUT_CARD.json").read_text())
    parent = json.loads((ROOT/card["parent_card"]).read_text())
    m = card["phase_outcomes"]
    phases = 2*np.pi*np.arange(m)/m
    theta = card["theta_over_pi"]*np.pi
    alpha = parent["alpha_per_assumed_update"]
    lam = spectrum(theta)
    h0 = np.diag((lam[:, None]+lam[None, :]).ravel())
    contact = np.ones((9, 9))/9
    pdark = proj(DPLUS)+proj(DMINUS)
    rho0 = proj(E9[:, 4])
    effects = [proj(word(phi))/m for phi in phases]
    joint_effects = [np.kron(a, b) for a in effects for b in effects]
    checks, rows = [], []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, atol=2e-11):
        check(name, np.allclose(a, b, atol=atol, rtol=2e-11))

    for rel, sha in card["preserved_sha256"].items():
        check("preserved:"+rel, hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == sha)
    close("seven_outcome_POVM_completeness", sum(effects), np.eye(3))
    check("phase_effect_positivity", min(np.linalg.eigvalsh(x).min() for x in effects) > -1e-13)
    close("joint_POVM_completeness", sum(joint_effects), E9)
    basis = [np.diag(np.eye(3)[k]).astype(complex) for k in range(3)]
    for k in range(3):
        for l in range(k+1, 3):
            e = np.outer(np.eye(3)[k], np.eye(3)[l])
            basis += [e+e.T, -1j*(e-e.T)]
    check("full_local_observation_algebra_rank9", np.linalg.matrix_rank(np.array(basis).reshape(9, 9)) == 9)
    check("full_joint_observation_algebra_rank81",
          np.linalg.matrix_rank(np.array([np.kron(a, b) for a in basis for b in basis]).reshape(81, 81)) == 81)
    check("phase_POVM_not_tomographically_complete",
          np.linalg.matrix_rank(np.array(effects).reshape(m, 9)) == 5)
    close("dark_projector_idempotent", pdark@pdark, pdark)
    close("dark_projector_rank2", np.trace(pdark), 2)
    close("dark_contact_annihilation", contact@pdark, np.zeros((9, 9)))
    close("initial_dark_weight_two_thirds", np.trace(pdark@rho0), 2/3)
    close("hypothetical_dark_filter_normalized_output", pdark@rho0@pdark/(2/3), proj(DPLUS))

    # Fixed total-index states have an exact conditional orbit, without n.
    dark_report = {}
    for name, psi in [("symmetric", DPLUS), ("antisymmetric", DMINUS)]:
        rho = proj(psi)
        p = joint(rho, phases)
        delta = phases[:, None]-phases[None, :]
        analytic = ((2/3)*(1-np.cos(delta))**2 if name == "symmetric"
                    else 2*np.sin(delta)**2)/m**2
        close(name+":joint_analytic", p, analytic)
        close(name+":both_phase_marginals_uniform", [p.sum(axis=0), p.sum(axis=1)], np.ones((2, m))/m)
        r0, _ = conditional(rho, 0, m)
        for r, phi in enumerate(phases):
            rcond, prob = conditional(rho, phi, m)
            u = np.diag(word(phi))
            close(name+f":conditional_unitary:{r}", rcond, u@r0@u.conj().T)
        unread = sum(conditional(rho, phi, m)[0]*conditional(rho, phi, m)[1] for phi in phases)
        close(name+":unread_A_is_partial_trace", unread, marginal_b(rho))
        diagonal = np.diag(np.diag(rho))
        close(name+":same_populations_dephased_joint_is_flat", joint(diagonal, phases), np.ones((m, m))/m**2)
        close(name+":dephasing_preserves_local_states", marginal_b(diagonal), marginal_b(rho))
        dark_report[name] = readout_metrics(rho, phases)
    close("symmetric:two_harmonics", [dark_report["symmetric"]["relative_first_harmonic_real"],
          dark_report["symmetric"]["relative_second_harmonic_real"]], [-2/3, 1/6])
    close("antisymmetric:twofold_ambiguity", [dark_report["antisymmetric"]["relative_first_harmonic_real"],
          dark_report["antisymmetric"]["relative_second_harmonic_real"]], [0, -.5])
    # General dark coherences and mixed states, not only two real pure examples.
    dk = np.column_stack([DPLUS, DMINUS])
    for index, a in enumerate([np.array([[.5, .2j], [-.2j, .5]]),
                               proj(np.array([1, 1j])/np.sqrt(2))]):
        rho = dk@a@dk.conj().T
        close(f"mixed_dark:{index}:rigidity", readout_metrics(rho, phases)["seven_reading_rigidity_residual"], 0)

    secular_rows = []
    for x in [0., .7, np.pi, 2*np.pi]:
        z = np.exp(-1j*x)
        a, b = (z-1)/3, (z+2)/3
        psi = a*E9[:, 2]+b*E9[:, 4]+a*E9[:, 6]
        delta = phases[:, None]-phases[None, :]
        close(f"secular:{x}:all_three_channels_joint", joint(proj(psi), phases), abs(b+2*a*np.cos(delta))**2/m**2)
        metric = readout_metrics(proj(psi), phases)
        close(f"secular:{x}:harmonics", [metric["relative_first_harmonic_real"],
              metric["relative_second_harmonic_real"]], [-abs(a)**2, abs(a)**2])
        secular_rows.append({"exchange_phase_x": x, "each_outlet_probability": abs(a)**2, **metric})

    bright = (E9[:, 2]+E9[:, 4]+E9[:, 6])/np.sqrt(3)
    weak_avg = (2/3)*proj(DPLUS)+(1/3)*proj(bright)
    close("weak_average:harmonics", [readout_metrics(weak_avg, phases)["relative_first_harmonic_real"],
          readout_metrics(weak_avg, phases)["relative_second_harmonic_real"]], [-2/9, 2/9])
    _, zero_projectors = energy_projectors(h0)
    close("zero_coupling_before_averaging_leaves_input", dephase_energy(rho0, zero_projectors), rho0)
    check("weak_limit_and_infinite_average_do_not_commute", trace_distance(weak_avg, rho0) > .4)
    weak_distances = []
    for gamma in parent["local_contact_gamma_fixtures"]:
        tag = f"gamma={gamma}"
        h = h0+gamma*contact
        # Independent site construction pins the unchanged physical contact.
        pv = np.diag([1., 0., 0.])
        l = laplacian(theta)
        site = np.kron(l, np.eye(3))+np.kron(np.eye(3), l)+gamma*np.kron(pv, pv)
        ff = np.kron(F, F)
        close(tag+":site_contact_equivalence", ff.conj().T@site@ff, h)
        close(tag+":dark_full_eigenstates", h@dk, 4*dk)
        close(tag+":dark_weight_is_conserved", h@pdark-pdark@h, np.zeros((9, 9)))
        values, projectors = energy_projectors(h)
        check(tag+":no_integer_update_spectral_alias", alpha*(values.max()-values.min()) < 2*np.pi)
        n = round(np.pi/(gamma/3*alpha))
        psi = expm(-1j*alpha*n*h)@E9[:, 4]
        actual = proj(psi)
        averaged = dephase_energy(rho0, projectors)
        weak_distances.append(trace_distance(averaged, weak_avg))
        close(tag+":unread_update_stationarity", h@averaged-averaged@h, np.zeros((9, 9)))
        close(tag+":averaged_trace", np.trace(averaged), 1)
        check(tag+":averaged_positivity", np.linalg.eigvalsh(averaged).min() > -1e-12)
        close(tag+":initial_dark_overlap_cannot_purify", abs(np.vdot(DPLUS, psi))**2, 2/3)
        close(tag+":antisymmetric_state_inaccessible", abs(np.vdot(DMINUS, psi))**2, 0)
        invariant_effects = [dephase_energy(e, projectors) for e in joint_effects]
        close(tag+":invariant_POVM_completeness", sum(invariant_effects), E9)
        close(tag+":invariant_effect_commutators", [h@e-e@h for e in invariant_effects], np.zeros((m*m, 9, 9)))
        check(tag+":invariant_effects_positive", min(np.linalg.eigvalsh(e).min() for e in invariant_effects) > -1e-12)
        for state_name, rho in [("exchange", actual), ("unread_updates", averaged), ("dark", proj(DPLUS))]:
            p = joint(rho, phases)
            close(tag+":"+state_name+":probability_sum", p.sum(), 1)
            check(tag+":"+state_name+":probability_positive", p.min() > -1e-12)
            conditional_probs = np.array([conditional(rho, phi, m)[1] for phi in phases])
            close(tag+":"+state_name+":marginal_likelihood", p.sum(axis=1), conditional_probs)
            unread = sum(conditional(rho, phi, m)[0]*conditional(rho, phi, m)[1] for phi in phases)
            close(tag+":"+state_name+":all_outcomes_retained", unread, marginal_b(rho))
            metric = readout_metrics(rho, phases)
            if state_name != "dark":
                check(tag+":"+state_name+":not_exact_free_rigidity", metric["seven_reading_rigidity_residual"] > 1e-4)
            if state_name != "exchange":
                p_inv = np.array([np.trace(rho@e).real for e in invariant_effects]).reshape(m, m)
                close(tag+":"+state_name+":invariant_readout_same_statistics", p_inv, p)
            rows.append({"gamma": gamma, "state": state_name, "exchange_update_n": n,
                         "stationarity_commutator_norm": float(np.linalg.norm(h@rho-rho@h)),
                         "joint_phase_probabilities": p.tolist(), **metric})
        # Exact projected stationary equation, retaining the contact's nonlocal term.
        vals, vecs = np.linalg.eigh(h)
        near = np.argmin(abs(vals-(4+gamma/3)))
        v = vecs[:, near].reshape(3, 3)
        pv_modes = np.ones((3, 3))/3
        chi0 = word(0).conj()@v/np.sqrt(3)
        for phi in [0., .37, 1.4]:
            w = word(phi)
            chi = w.conj()@v/np.sqrt(3)
            derivative = (-1j*np.arange(3)*w.conj())@v/np.sqrt(3)
            f = np.sum(w.conj())/3
            rhs = (vals[near]-2*lam[0])*chi-np.sqrt(3)*N@chi-gamma*f*pv_modes@chi0
            close(tag+f":conditional_stationary_equation:{phi}", 1j*np.sqrt(3)*derivative, rhs)

    check("full_average_approaches_weak_limit_on_fixed_fixtures",
          weak_distances[0] > weak_distances[1] > weak_distances[2])

    # The symmetry-protected antisymmetric pair survives flux detuning;
    # the symmetric three-channel dark eigenstate requires equal spacing.
    for theta in [np.pi/3, 2*np.pi/3]:
        lam_alt = spectrum(theta)
        hh = np.diag((lam_alt[:, None]+lam_alt[None, :]).ravel())+.09*contact
        close(f"detuned:{theta}:antisymmetric_exact", hh@DMINUS, (lam_alt[0]+lam_alt[2])*DMINUS)
        mean = np.vdot(DPLUS, hh@DPLUS).real
        check(f"detuned:{theta}:symmetric_not_stationary", np.linalg.norm(hh@DPLUS-mean*DPLUS) > .1)

    failed = [c["name"] for c in checks if not c["passed"]]
    summary = {"checks": len(checks), "passed": len(checks)-len(failed), "failed": failed}
    report = {"stage": card["stage"], "status": card["status"], "summary": summary,
              "dark_states": dark_report, "secular_fixtures": secular_rows,
              "weak_limit_trace_distances_in_gamma_card_order": weak_distances,
              "full_contact_readout": rows, "checks": checks, "boundary": card["not_claimed"]}
    OUT.with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n")
    lines = ["# G-R10 all-channel relational readout", "",
             f"{summary['passed']}/{summary['checks']} checks pass. No fitted input or physical clock-unit calibration.", "",
             "## New result: the extra channel supports a stationary conditional clock", "",
             "At theta=pi/2 the SAME full local contact has two exact E=4 dark states.",
             "Dplus=(|0,2>-2|1,1>+|2,0>)/sqrt(6) is annihilated by the contact.",
             "For the fixed seven-outcome phase POVM all A readings have probability 1/7,",
             "while rho(B|phi)=exp(i N_B phi) rho(B|0) exp(-i N_B phi).",
             "This relation contains no update number n. It is finite conditional quantum",
             "evolution, not a new Page-Wootters mechanism, physical seconds or an arrow.", "",
             "The full local observation algebra has dimension 9 and its joint product 81.",
             "The phase POVM alone spans only 5 local dimensions and is not tomography.",
             "No exit is discarded: summing ALL A results gives the ordinary B partial trace.",
             "For Dplus the relative phase density is (1-cos(delta))^2/(3*pi); it is broad,",
             "not an exact timestamp. Local dephasing leaves all populations/marginals",
             "unchanged but makes the joint phase table flat. This is a coherence control,",
             "not a unique entanglement witness. Dminus has a pi-periodic twofold ambiguity.", "",
             "## Same contact, three physically different states", "",
             "The exchange input |1,1>, its uniformly unread-update average, and a prepared",
             "pure dark state must not be identified. The average uses an assumed Cesaro prior.",
             "It is stationary but need not follow the bare-gap rigid conditional orbit.",
             "The residual below is a maximum over the seven FIXED readings, not a continuous",
             "optimization or a refutation of general interacting-clock descriptions.", "",
             "| gamma | state | seven-reading trace-distance residual | scaled A likelihood range |",
             "|---|---|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['gamma']} | {r['state']} | {r['seven_reading_rigidity_residual']:.7g} | {r['scaled_A_likelihood_range']:.7g} |")
    lines += ["", "An exact projected stationary equation includes a contact term linking",
              "the conditional amplitude at phi to its value at phi=0. It vanishes for",
              "the dark state, not for generic dressed states. No interacting term is discarded.",
              "Energy-averaged joint effects are positive, complete and stationary and",
              "give the same statistics on stationary states; their autonomous local",
              "apparatus implementation and record energetics are NOT established.", "",
              "## Preparation limit and next bounded decision", "",
              "The original input has conserved dark weight 2/3, not 1. Fixed unitary contact",
              "cannot purify it to Dplus. A hypothetical dark-energy filter succeeds with",
              "probability 2/3, but is an additional instrument, not dynamics already derived.",
              "Dminus survives equal-copy flux detuning but is inaccessible from the symmetric",
              "input. Dplus supplies a three-level clock but requires the equal-gap resonance.",
              "The gamma->0 positive-coupling averaged state is (2/3)|Dplus><Dplus|",
              "+(1/3)|bright><bright|; taking gamma=0 BEFORE averaging instead leaves |1,1>.",
              "Thus weak coupling and infinite unread-update averaging cannot be interchanged.", "",
              "Next choose one bounded physical test: specify a preparation plus one-record",
              "instrument, keeping its energy/reference/outcome cost, or retain the unfiltered",
              "mixed stationary state and test its genuinely interacting conditional equation.",
              "Do not add a loss bath, discard an outcome, or change the interaction without",
              "explicitly declaring the added mechanism. No Route-F gate is promoted.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
