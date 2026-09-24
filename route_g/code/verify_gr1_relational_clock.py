#!/usr/bin/env python3
"""G-R1: a finite Page-Wootters clock, not a derivation of Chen's Ed.

All numbers are dimensionless algebra tests. No physical time, spatial
volume, detector calibration, particle mass or gravity is inferred.
Prior KK/TOF files are hashed but not modified.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "gr1_relational_clock"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def projector(v):
    return np.outer(v, v.conj())


def clock_ket(d, phi):
    return np.exp(1j*np.arange(d)*phi)/np.sqrt(d)


def model(amplitudes, gap=1.0):
    a = np.asarray(amplitudes, dtype=complex)
    a = a/np.linalg.norm(a)
    d = len(a)
    ec = gap*(d-1-np.arange(d))
    es = gap*np.arange(d)
    total = (ec[:, None]+es[None, :]).ravel()
    constraint = np.diag(total-gap*(d-1))
    state = np.diag(a).ravel()
    rho = projector(state)
    mixture = np.diag(np.diag(rho))
    return dict(d=d, a=a, ec=ec, es=es, total=total,
                K=constraint, psi=state, rho=rho, mix=mixture, gap=gap)


def reductions(rho, d):
    tensor = rho.reshape(d, d, d, d)
    return np.einsum("abcb->ac", tensor), np.einsum("abad->bd", tensor)


def conditioned(rho, d, phi):
    """Clock effect |phi><phi|; density prefactor cancels in the conditional."""
    v = clock_ket(d, phi)
    unnormalized = np.einsum("a,abcd,c->bd", v.conj(),
                            rho.reshape(d, d, d, d), v)
    p = float(np.trace(unnormalized).real)
    if p <= 0:
        raise ValueError("Zero-probability clock conditioning")
    return unnormalized/p, p


def twirl(effect, generators, order):
    """Exact Fourier average for this integer-spectrum finite model."""
    phases = 2*np.pi*np.arange(order)/order
    delta = generators[:, None]-generators[None, :]
    kernel = np.exp(-1j*phases[:, None, None]*delta[None, :, :]).mean(axis=0)
    return effect*kernel


def run():
    checks = []

    def check(name, ok, value=None):
        row = dict(name=name, passed=bool(ok))
        if value is not None:
            row["value"] = float(value)
        checks.append(row)

    old_files = sorted((ROOT/"code").glob("verify_g[12]_*.py"))
    old_files += [ROOT/"output"/"g2_higgs_source.json", ROOT/"output"/"g2_tof_mode_audit.json"]
    old_hashes = {str(p.relative_to(ROOT)): sha(p) for p in old_files}
    cards = [
        ("three_level_uniform", np.ones(3), 1.0),
        ("four_level_nonuniform", np.sqrt([.4, .3, .2, .1])
         * np.exp(1j*np.array([0, .3, -.5, .8])), 1.7),
    ]
    maxima = dict(constraint_norm=0., conditional_state_error=0.,
                  conditional_equation_error=0., invariant_effect_error=0.)
    phase_rows = []
    for name, a, gap in cards:
        m = model(a, gap)
        d, K, psi = m["d"], m["K"], m["psi"]
        c_res = np.linalg.norm(K@psi)
        maxima["constraint_norm"] = max(maxima["constraint_norm"], float(c_res))
        check(name+":positive_local_hamiltonians", min(m["ec"].min(), m["es"].min()) >= 0)
        check(name+":normalized", abs(np.vdot(psi, psi)-1) < 2e-14)
        check(name+":constraint", c_res < 2e-14, c_res)
        check(name+":mixture_constraint", np.linalg.norm(K@m["mix"]) < 2e-14)
        rc, rs = reductions(m["rho"], d)
        mc, ms = reductions(m["mix"], d)
        check(name+":identical_marginals", np.allclose(rc, mc) and np.allclose(rs, ms))
        check(name+":identical_joint_energy_statistics",
              np.array_equal(np.diag(m["rho"]), np.diag(m["mix"])))
        grid = 2*np.pi*np.arange(d)/d
        clock_basis = np.column_stack([clock_ket(d, p) for p in grid])
        check(name+":orthonormal_clock_grid", np.linalg.norm(clock_basis.conj().T@clock_basis-np.eye(d)) < 2e-14)
        # Uniform quadrature is exact because all clock Fourier modes have |n-m|<d.
        phases = 2*np.pi*np.arange(2*d+1)/(2*d+1)
        povm_integral = sum(d*projector(clock_ket(d, p))/(2*d+1) for p in phases)
        check(name+":continuous_clock_POVM", np.linalg.norm(povm_integral-np.eye(d)) < 2e-14)
        alpha = .713
        v0 = clock_ket(d, .27)
        rotated = np.exp(-1j*m["ec"]*alpha/gap)*v0
        covariant = np.exp(-1j*(d-1)*alpha)*clock_ket(d, .27+alpha)
        check(name+":clock_covariance", np.linalg.norm(rotated-covariant) < 2e-14)
        gauge = np.exp(-1j*np.diag(K)*alpha/gap)
        check(name+":global_stationary", np.linalg.norm(gauge*psi-psi) < 2e-14)
        dial = clock_basis@np.diag(grid)@clock_basis.conj().T
        comm = dial@np.diag(m["ec"])-np.diag(m["ec"])@dial
        check(name+":finite_commutator_trace_zero", abs(np.trace(comm)) < 1e-14)
        check(name+":not_canonical_time_operator", np.linalg.norm(comm-1j*np.eye(d)) > 1)
        history = sum(np.kron(clock_ket(d, p), m["a"]*np.exp(-1j*np.arange(d)*p))
                      for p in grid)/np.sqrt(d)
        check(name+":history_identity", np.linalg.norm(history-psi) < 2e-14)
        # A naive product of clock dial and system superposition is NOT physical here.
        product = np.kron(clock_ket(d, 0), m["a"])
        check(name+":arbitrary_product_fails_constraint", np.linalg.norm(K@product) > .1)
        plus = np.ones(d)/np.sqrt(d)
        measurement = projector(plus)
        for phi in [0., np.pi/3, 2*np.pi/3, np.pi, .713]:
            rho_phi, p = conditioned(m["rho"], d, phi)
            sigma_phi, p_mix = conditioned(m["mix"], d, phi)
            expected = m["a"]*np.exp(-1j*np.arange(d)*phi)
            err = np.linalg.norm(rho_phi-projector(expected))
            maxima["conditional_state_error"] = max(maxima["conditional_state_error"], float(err))
            tag = f"{name}:phi={phi:.6g}"
            check(tag+":conditional_state", err < 2e-14, err)
            check(tag+":clock_probability", abs(p-1/d) < 2e-14 and abs(p_mix-p) < 2e-14)
            check(tag+":dephased_state_constant", np.linalg.norm(sigma_phi-np.diag(abs(m["a"])**2)) < 2e-14)
            # Independently differentiate the directly projected global ket.
            step = 2e-6
            def project_ket(x):
                return np.sqrt(d)*np.einsum("a,ab->b", clock_ket(d, x).conj(), psi.reshape(d, d))
            derivative = (project_ket(phi+step)-project_ket(phi-step))/(2*step)
            eq_err = np.linalg.norm(1j*gap*derivative-m["es"]*expected)
            maxima["conditional_equation_error"] = max(maxima["conditional_equation_error"], float(eq_err))
            check(tag+":conditional_Schrodinger", eq_err < 2e-9, eq_err)
            # Twirled effects commute with the constraint and keep physical probabilities.
            effect = np.kron(projector(clock_ket(d, phi)), measurement)
            invariant = twirl(effect, np.diag(K)/gap, 4*d+1)
            inv_err = np.linalg.norm(K@invariant-invariant@K)
            maxima["invariant_effect_error"] = max(maxima["invariant_effect_error"], float(inv_err))
            check(tag+":invariant_joint_effect", inv_err < 2e-14 and np.linalg.eigvalsh(invariant).min() > -2e-14)
            check(tag+":twirl_preserves_probability",
                  abs(np.trace(m["rho"]@(effect-invariant))) < 2e-14)
            if name == "three_level_uniform":
                measured = float(np.trace(rho_phi@measurement).real)
                mixed = float(np.trace(sigma_phi@measurement).real)
                exact = float((1+2*np.cos(phi))**2/9)
                check(tag+":analytic_fringe", abs(measured-exact) < 2e-14 and abs(mixed-1/3) < 2e-14)
                phase_rows.append(dict(phi_over_pi=float(phi/np.pi),
                                       coherent_probability=measured,
                                       dephased_probability=mixed,
                                       exact_coherent_probability=exact))
        # Global constraint scaling and common energy rescaling don't change phase correlations.
        scaled = model(a, 2*gap)
        check(name+":common_energy_scaling_same_state", np.array_equal(scaled["rho"], m["rho"]))
        check(name+":constraint_scaling_same_kernel",
              np.allclose(scaled["K"], 2*K) and np.linalg.norm(scaled["K"]@psi) < 4e-14)
        shift = .73
        shifted_K = np.diag((m["ec"][:, None]+shift+m["es"][None, :]).ravel()
                            -(gap*(d-1)+shift))
        check(name+":energy_zero_shift_same_constraint", np.linalg.norm(shifted_K-K) < 2e-14)

    m = model(np.ones(3))
    plus = np.ones(3)/np.sqrt(3)
    coherent, _ = conditioned(m["rho"], 3, 0)
    dephased, _ = conditioned(m["mix"], 3, 0)
    mean_energy = float(np.dot(np.diag(m["rho"]).real, m["total"]))
    variance = float(np.dot(np.diag(m["rho"]).real, (m["total"]-mean_energy)**2))
    check("counterexample:coherent_purity", abs(np.trace(coherent@coherent)-1) < 2e-14)
    check("counterexample:dephased_purity", abs(np.trace(dephased@dephased)-1/3) < 2e-14)
    check("counterexample:fixed_total_energy", abs(mean_energy-2) < 2e-14 and variance < 2e-14)
    check("counterexample:nonzero_probability_difference",
          abs(np.trace((coherent-dephased)@projector(plus))-2/3) < 2e-14)
    # Prescribed finite dephasing examples, not fits or a detector-noise scan.
    for eta in [0., .5, 1.]:
        rho = eta*m["rho"]+(1-eta)*m["mix"]
        r, _ = conditioned(rho, 3, np.pi/3)
        actual = float(np.trace(r@projector(plus)).real)
        check(f"coherence_fraction:{eta}", abs(actual-(eta*4/9+(1-eta)/3)) < 2e-14)
    # Algebraic check of the interaction projection identity. V is arbitrary,
    # and this state is NOT asserted to solve the interacting constraint.
    V = np.kron(np.array([[0, 1, 0], [1, 0, .2], [0, .2, 0]]),
                np.array([[1, .3, 0], [.3, 0, .1], [0, .1, -1]]))
    phi = .41
    lhs = np.sqrt(3)*clock_ket(3, phi).conj()@(V@m["psi"]).reshape(3, 3)
    rhs = np.zeros(3, dtype=complex)
    for p in 2*np.pi*np.arange(7)/7:
        kernel = np.einsum("a,abcd,c->bd", clock_ket(3, phi).conj(),
                           V.reshape(3, 3, 3, 3), clock_ket(3, p))
        rhs += 3/7*kernel@(m["a"]*np.exp(-1j*np.arange(3)*p))
    check("interaction:projection_identity_only", np.linalg.norm(lhs-rhs) < 2e-14)
    check("interaction:not_a_stationary_completion", np.linalg.norm((m["K"]+V)@m["psi"]) > .1)
    check("prior:all_files_unchanged", all(sha(ROOT/p) == h for p, h in old_hashes.items()))
    failed = [r["name"] for r in checks if not r["passed"]]
    result = dict(
        status="G-R1 bounded-done: finite, noninteracting relational-clock theorem and density-only counterexample.",
        assumptions=["Finite clock/system factorization is given, not derived.",
                     "Local Hamiltonians have complementary equally spaced nonnegative energies.",
                     "The matched joint state, Born rule and covariant clock measurement are specified.",
                     "Clock phase is a readout label, not a physical extra dimension or external elapsed time.",
                     "Seconds require calibration tau=hbar*phi/epsilon; no numerical physical epsilon is chosen.",
                     "No memory register, repeated-measurement experiment, time arrow, spacetime metric or Ed dynamics is supplied."],
        example=dict(d=3, gap_model_units=1, energy_clock=[2, 1, 0],
                     energy_system=[0, 1, 2], total_energy=mean_energy,
                     total_energy_variance=variance, rows=phase_rows,
                     coherent_conditional_purity=float(np.trace(coherent@coherent).real),
                     dephased_conditional_purity=float(np.trace(dephased@dephased).real)),
        counterexample_scope="Same local reduced density matrices and complete joint energy statistics, different clock-conditioned phase-sensitive probabilities. Refutes identifying an energy density alone with the full relational history, not every possible definition of Chen's unspecified Ei/Ed.",
        not_proved=["Ed increases imply clock slowdown", "mass spectrum or negative mass",
                    "emergent space or general relativity", "unbounded event count or time arrow",
                    "an actual sensor/source", "sequential measurement statistics"],
        error_maxima=maxima, summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed),
        checks=checks, verification=dict(source_sha256=sha(__file__), unchanged_comparison_files=old_hashes))
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# G-R1: a finite relational clock and an energy-density-only counterexample", "",
             f"Verification: **{len(checks)-len(failed)}/{len(checks)} checks passed**.", "",
             "The relational-time line is now Route G's priority. Older KK/TOF files are retained unchanged as comparison models. No dimensional energy, measured particle mass, external elapsed time or spatial geometry is fitted.", "",
             "## Exact finite construction", "",
             "For n=0,...,d-1, HC|n>=epsilon(d-1-n)|n>, HS|n>=epsilon n|n>. Both energies are nonnegative. The state |Psi>=sum a_n |n,n> satisfies (HC+HS-(d-1)epsilon)|Psi>=0.",
             "Clock |phi>=sum exp(i n phi)|n>/sqrt(d); its d grid readings are orthonormal, while the continuous phase dial is a normalized POVM, not uncountably many orthogonal states.",
             "Conditioning gives |psi(phi)>=sum a_n exp(-i n phi)|n>, so i epsilon d_phi |psi>=HS|psi>. Only after calibration tau=hbar phi/epsilon does this become the usual Schrödinger equation. No geometric rotation is introduced.", "",
             "## Three-level example and exact counterexample", "",
             "Uniform amplitudes, epsilon=1 in model units. Compare the coherent |Psi><Psi| with sum |n,n><n,n|/3. Both have identical local marginals, identical joint energy statistics, total energy2 and variance0. Measure |+>=(|0>+|1>+|2>)/sqrt(3) conditionally:", "",
             "| clock phase / pi | coherent P(+|phi) | dephased P(+|phi) |",
             "|---:|---:|---:|"]
    for row in phase_rows:
        lines.append(f"| {row['phi_over_pi']:.7g} | {row['coherent_probability']:.9g} | {row['dephased_probability']:.9g} |")
    lines += ["", "Exact laws: P_coherent=(1+2 cos(phi))²/9; P_dephased=1/3. At phi=0 the probabilities differ by2/3 despite identical energy data. The missing information is joint coherence, not a decimal correction to the energy spectrum.",
              "Even after supplying a fixed volume V0, mean energy/V0 would be identical. This rules out that density-only identification, not a richer Ei/Ed that explicitly includes relational coherences and has its own dynamics.", "",
              "## Guardrails that passed", "",
              "- Constraint, positivity, clock-grid completeness, continuous POVM and history-state identities.",
              "- Direct projection versus unitary conditional state, finite-difference conditional equation, arbitrary normalized nonuniform amplitudes.",
              "- Constraint-invariant joint effects obtained by exact Fourier twirling preserve the same conditional probabilities.",
              "- Finite clock commutator has zero trace; no impossible exact canonical time operator is assumed.",
              "- Common energy scaling changes calibration but not phase-conditioned correlations; a local energy-zero shift with matching total shift leaves the constraint unchanged.",
              "- A naive clock/system product generally fails the constraint. An arbitrary clock interaction does not preserve the old stationary state.",
              "- All comparison code and selected earlier result hashes are unchanged.", "",
              "## Not closed / next step", "",
              "The model assumes quantum mechanics, the tensor split, Hamiltonians, a state and a clock POVM. It is an exact conditional construction, not an unconditional derivation of time from nothing. The clock is cyclic and single-event conditioning is not a sequential measurement history.",
              "G-R2 should add a minimal explicit record/measurement model before interpreting the original T as accumulated event count. G-R3 can then define an Ei/Ed candidate and one shared clock-matter constraint, and compare two physical clocks. Do not insert a fitted slowdown function or relabel the phase as an extra spatial coordinate.",
              "Full derivation and dictionary: tex/route_g_relational_clock.tex and RELATIONAL_TIME_CARD.md. The earlier Higgs/xenon rejection is neither undone nor a bound on this unspecified new realization."]
    OUT.with_suffix(".md").write_text("\n".join(lines)+"\n")
    print(json.dumps(result["summary"]))
    print(json.dumps(maxima))
    print("Exact counterexample: same energies, P(+|phi=0)=1 versus 1/3.")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
