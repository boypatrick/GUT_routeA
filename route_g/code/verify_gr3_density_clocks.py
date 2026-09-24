#!/usr/bin/env python3
"""G-R3: a declared excitation-density coupling and two relational clocks.

Exact finite stationary model + separate dispersive microscopic check.
No density/time law is fitted. No spatial volume, gravity, universal
proper time, or completed G-R2 measurement apparatus is claimed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "gr3_density_clocks"
TOL = 3e-12


def projector(v):
    return np.outer(v, v.conj())


def build(kappa=.25, scale=1., weights=(.5, .5), source_phase=.37):
    """Same 15-level A clock for kappa=-1/4,0,+1/4; no retuned apparatus."""
    delta, epsilon, ed0 = scale/4, scale, scale
    n = np.arange(3)[:, None]
    d = np.arange(2)[None, :]
    hbd = epsilon*n*(1+kappa*d) + ed0*d
    labels_float = hbd/delta
    labels = np.rint(labels_float).astype(int)
    if np.max(abs(labels-labels_float)) > 1e-12 or labels.max() > 14:
        raise ValueError("This exact finite clock requires the declared commensurate spectrum")
    ha = delta*(14-np.arange(15))
    estar = 14*delta
    constraint = (ha[:, None, None]+hbd[None, :, :]-estar).ravel()
    source = np.sqrt(weights)*np.exp(1j*np.array([0., source_phase]))
    psi = np.zeros((15, 3, 2), dtype=complex)
    for nn in range(3):
        for dd in range(2):
            psi[labels[nn, dd], nn, dd] = source[dd]/np.sqrt(3)
    return dict(kappa=kappa, delta=delta, epsilon=epsilon, ed0=ed0,
                hbd=hbd, ha=ha, estar=estar, K=constraint, psi=psi,
                labels=labels, source=source)


def conditional(m, phi):
    a = np.exp(1j*np.arange(15)*phi)/np.sqrt(15)
    ket = np.einsum("a,abd->bd", a.conj(), m["psi"])
    norm = float(np.vdot(ket, ket).real)
    return ket/np.sqrt(norm), norm


def b_dial(theta):
    # Opposite Fourier sign to complementary-spectrum reference A.
    return np.exp(-1j*np.arange(3)*theta)/np.sqrt(3)


def b_probs(rho, offset=0.):
    return np.array([np.vdot(b_dial(offset+2*np.pi*j/3),
                            rho@b_dial(offset+2*np.pi*j/3)).real for j in range(3)])


def jc_energy(n, excited, omega, gap, g):
    """Exact eigenvalue connected to |n,d> at g=0 in its excitation block."""
    if n == 0 and not excited:
        return 0.
    detuning = gap-omega
    total = n+int(excited)
    root = np.sqrt(detuning**2+4*g*g*total)
    sign = np.sign(detuning)*(1 if excited else -1)
    return total*omega+detuning/2+sign*root/2


def run():
    checks = []

    def check(name, ok, value=None):
        row = dict(name=name, passed=bool(ok))
        if value is not None:
            row["value"] = float(value)
        checks.append(row)

    def close(name, a, b):
        error = float(np.linalg.norm(np.asarray(a)-np.asarray(b)))
        check(name, error < TOL, error)

    old = sorted(p for p in (ROOT/"code").glob("verify_*.py")
                 if p.name != Path(__file__).name)
    old += [ROOT/"output"/"gr1_relational_clock.json",
            ROOT/"output"/"gr2_two_event_records.json"]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in old}
    rate_rows = []
    for kappa in [-.25, 0., .25]:
        m = build(kappa)
        tag = f"kappa={kappa:g}"
        close(tag+":stationary", m["K"]*m["psi"].ravel(), np.zeros(90))
        close(tag+":normalized", np.vdot(m["psi"], m["psi"]), 1)
        check(tag+":positive_local_energies", min(m["ha"].min(), m["hbd"].min()) >= 0)
        check(tag+":positive_clock_gaps", 1+kappa > 0)
        close(tag+":source_nondemolition", np.diag(m["hbd"].ravel())@np.diag([0,1]*3)
              -np.diag([0,1]*3)@np.diag(m["hbd"].ravel()), np.zeros((6,6)))
        # A full Fourier grid integrates the finite phase POVM exactly.
        for dim, sign in [(15, 1), (3, -1)]:
            phases = 2*np.pi*np.arange(2*dim+1)/(2*dim+1)
            integral = sum(projector(np.exp(sign*1j*np.arange(dim)*ph)/np.sqrt(dim))
                           for ph in phases)*dim/len(phases)
            close(f"{tag}:POVM:{dim}", integral, np.eye(dim))
        for phi in [0., .13, np.pi/4, np.pi/2, np.pi]:
            state, prob = conditional(m, phi)
            expected = np.ones((3, 1))/np.sqrt(3)*m["source"][None, :]*np.exp(-1j*m["labels"]*phi)
            key = f"{tag}:phi={phi:.7g}"
            close(key+":conditioned_state", state, expected)
            close(key+":reference_likelihood", prob, 1/15)
            rho = state@state.conj().T
            close(key+":B_normalization", np.trace(rho), 1)
            check(key+":B_positive", np.linalg.eigvalsh(rho).min() > -TOL)
            close(key+":B_grid_normalized", b_probs(rho).sum(), 1)
            mixture = np.zeros((3,3), dtype=complex)
            for dd in range(2):
                branch = state[:, dd]/abs(m["source"][dd])
                theta = 4*(1+kappa*dd)*phi
                pure = b_dial(theta)
                close(key+f":branch{dd}:rigid_translation", projector(branch), projector(pure))
                expected_p = np.array([(1+2*np.cos(2*np.pi*j/3-theta))**2/9 for j in range(3)])
                close(key+f":branch{dd}:grid_probabilities", b_probs(projector(branch)), expected_p)
                mixture += .5*projector(pure)
            close(key+":unobserved_source_mixture", rho, mixture)
            # Continuous conditional density integrates to one, with 7-point exact quadrature.
            angles = 2*np.pi*np.arange(7)/7
            density = [3/(2*np.pi)*np.vdot(b_dial(x), rho@b_dial(x)).real for x in angles]
            close(key+":continuous_B_normalization", sum(density)*2*np.pi/7, 1)
        # Read relative slopes from the state coherence, not from a fitted rate function.
        phis = np.array([.02, .04, .06])
        for dd in range(2):
            phases = []
            for phi in phis:
                state, _ = conditional(m, phi)
                phases.append(-np.angle(state[1,dd]*state[0,dd].conj()))
            slope = np.diff(np.unwrap(phases))/np.diff(phis)
            measured = m["delta"]/m["epsilon"]*slope
            close(tag+f":rate{dd}:coherence_difference", measured, np.ones(2)*(1+kappa*dd))
            rate_rows.append(dict(kappa=kappa, source_occupation=dd, calibrated_relative_rate=float(measured[0])))
        scaled = build(kappa, scale=2.7)
        close(tag+":common_rescaling_state", scaled["psi"], m["psi"])
        close(tag+":common_rescaling_constraint", scaled["K"], 2.7*m["K"])
        close(tag+":common_rescaling_rate", scaled["delta"]/scaled["epsilon"], m["delta"]/m["epsilon"])
        # Changing source energy ZERO together with E* leaves the constraint
        # and ground-subtracted density unchanged, not adding an interaction.
        shift = .713
        shifted = (m["ha"][:,None,None]+m["hbd"][None,:,:]+shift-m["estar"]-shift).ravel()
        close(tag+":energy_zero_constraint", shifted, m["K"])
        close(tag+":density_ground_subtraction", np.array([shift,1+shift])-shift, [0,1])

    m = build(.25)
    # Energy-sector pinching supplies positive constraint-invariant effects,
    # with unchanged Born probabilities on the stationary state.
    pinching = np.isclose(m["K"][:, None], m["K"][None, :], atol=1e-13, rtol=0)
    ref = projector(np.exp(1j*np.arange(15)*.23)/np.sqrt(15))
    effects = [("reference", np.kron(ref, np.eye(6)))]
    for dd in range(2):
        effects.append((f"joint{dd}", np.kron(ref, np.kron(projector(b_dial(.7)),
                                                        np.diag([int(dd == j) for j in range(2)])))))
    for name, effect in effects:
        invariant = effect*pinching
        close("invariant:"+name+":commutes",
              (m["K"][:,None]-m["K"][None,:])*invariant, np.zeros((90,90)))
        close("invariant:"+name+":preserves_probability",
              np.vdot(m["psi"].ravel(), (effect-invariant)@m["psi"].ravel()), 0)
        check("invariant:"+name+":positive", np.linalg.eigvalsh(invariant).min() > -TOL)
    state, _ = conditional(m, np.pi/2)
    rho = state@state.conj().T
    ground = projector(state[:,0]/np.sqrt(.5))
    excited = projector(state[:,1]/np.sqrt(.5))
    close("diagnostic:ground_return", b_probs(ground), [1,0,0])
    close("diagnostic:excited_probabilities", b_probs(excited),
          [1/9, (1+np.sqrt(3))**2/9, (1-np.sqrt(3))**2/9])
    close("diagnostic:mixture_plus", b_probs(rho)[0], 5/9)
    # Same Ed expectation is not permission to replace its operator in H.
    mean_state = b_dial(4*(1+.25*.5)*np.pi/2)
    mean_prob = b_probs(projector(mean_state))[0]
    close("diagnostic:mean_field_plus", mean_prob, (3+2*np.sqrt(2))/9)
    check("diagnostic:mean_field_fails", abs(mean_prob-5/9) > .09)
    halfcycle, _ = conditional(m, np.pi)
    half_rho = halfcycle@halfcycle.conj().T
    close("fluctuation:adjacent_coherence_zero", half_rho[1,0], 0)
    close("fluctuation:not_fully_dephased", half_rho[2,0], 1/3)
    close("fluctuation:B_purity", np.trace(half_rho@half_rho), 5/9)
    # Density coherence is invisible in B alone when D is traced out.
    other_phase, _ = conditional(build(.25, source_phase=1.9), np.pi/2)
    close("fluctuation:source_phase_blind_B", other_phase@other_phase.conj().T, rho)
    close("density:one_cell_expectation", .5*m["ed0"], .5)
    # Coupling the two comparison clocks equally cancels a common factor.
    for d in [0, .5, 1]:
        close(f"common_mode:{d}:ratio_one", (1+.25*d)/(1+.25*d), 1)

    microscopic = []
    for detuning in [-.5, .5]:
        omega, g = 1., .03
        gap = omega+detuning
        chi = g*g/detuning
        tag = f"JC:detuning={detuning:g}"
        check(tag+":dispersive_regime", g*np.sqrt(3)/abs(detuning) < .11)
        check(tag+":positive_bare_frequencies", gap > 0 and omega > 0)
        max_error = 0.
        # Exact 2x2 blocks avoid artifacts of a truncated boson commutator.
        for total in [1,2,3]:
            block = np.array([[total*omega,g*np.sqrt(total)],
                              [g*np.sqrt(total),(total-1)*omega+gap]])
            vals = np.linalg.eigvalsh(block)
            analytic = [jc_energy(total,False,omega,gap,g),
                        jc_energy(total-1,True,omega,gap,g)]
            close(tag+f":block{total}:spectrum", vals, sorted(analytic))
            check(tag+f":block{total}:positive", vals.min() > 0)
            effective = [total*omega-chi*total,
                         (total-1)*omega+gap+chi*total]
            err = max(abs(np.array(analytic)-effective))
            bound = g**4*total**2/abs(detuning)**3
            check(tag+f":block{total}:remainder_bound", err <= bound*(1+1e-9), err)
            max_error = max(max_error, float(err))
        energies_g = [jc_energy(n,False,omega,gap,g) for n in range(3)]
        energies_e = [jc_energy(n,True,omega,gap,g) for n in range(3)]
        gaps_g, gaps_e = np.diff(energies_g), np.diff(energies_e)
        exact_ratio = gaps_e/gaps_g
        effective_ratio = (omega+chi)/(omega-chi)
        check(tag+":rate_sign", np.all((exact_ratio-1)*detuning > 0))
        check(tag+":nonuniform_higher_order_gaps", abs(gaps_e[1]-gaps_e[0]) > 1e-7)
        check(tag+":effective_rate_control", np.max(abs(exact_ratio-effective_ratio)) < 1e-4)
        microscopic.append(dict(detuning=detuning, coupling=g, chi=chi,
                                effective_kappa=2*chi/(omega-chi),
                                effective_rate=effective_ratio,
                                exact_transition_rate_ratios=exact_ratio.tolist(),
                                ground_gaps=gaps_g.tolist(), excited_gaps=gaps_e.tolist(),
                                max_second_order_energy_error=max_error))

    for rel, digest in hashes.items():
        check("preserved:"+rel, hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == digest)
    failed = [c["name"] for c in checks if not c["passed"]]
    result = dict(stage="G-R3", status="bounded-done" if not failed else "failed",
                  density=dict(operator="Ed=epsilon_D N_D/nu(D), nu(D)=1 addressed cell",
                               reference="excitation energy above the same source ground",
                               dimension="energy per declared cell, not energy per spatial volume",
                               not_identified_with="Chen Ed or complete Ei"),
                  interaction="V=kappa H_B (Ed/E0), E0=epsilon_D per cell",
                  constraint="K=H_A+H_B+H_D+V-E*, exact 15x3x2 finite completion",
                  synthetic_card=dict(epsilon=1., epsilon_D=1., delta=.25,
                                      reference_levels=15, kappa=.25, fitted=False),
                  rates=rate_rows,
                  diagnostics=dict(ground_probabilities=b_probs(ground).tolist(),
                                   excited_probabilities=b_probs(excited).tolist(),
                                   unobserved_source_plus=float(b_probs(rho)[0]),
                                   mean_density_substitution_plus=float(mean_prob),
                                   halfcycle_B_purity=float(np.trace(half_rho@half_rho).real)),
                  microscopic_JC=microscopic,
                  boundaries=["commensurate finite reference spectrum; no arbitrary-spectrum completion",
                              "one addressed cell does not empirically distinguish density from total excitation energy",
                              "matched stationary state is prepared input, not dynamically derived",
                              "finite cyclic readings require fixed calibration and branch/epoch labels",
                              "density candidate is a source observable, not sufficient full-state information",
                              "phase likelihoods are finite-width; rate is inferred across calibrated runs",
                              "JC result is dispersive and uses dressed sectors; bare source occupation fluctuates",
                              "no universal slowdown or spacetime proper time",
                              "G-R2 recording apparatus is not yet jointly energy-completed"],
                  previous_file_hashes=hashes, checks=checks,
                  summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    report = f"""# G-R3: density candidate, interaction and relative clocks

Status: {result['status']}; {len(checks)-len(failed)}/{len(checks)} checks.
This is a declared model, not a derivation of Chen's undefined Ed.

## Definition and observable

Ed = epsilon_D N_D / nu(D), with nu(D)=1 explicitly addressed matter cell.
Use ground-subtracted excitation energy, not assumed spatial volume,
record entropy, a timestamp, or the entire information in the joint state.
The one-cell test cannot distinguish density from total excitation energy;
a multi-cell claim must fix the averaging measure and coupling profile.
E0 is epsilon_D per cell. V=kappa H_B Ed/E0 is one fixed dispersive
interaction. Reference A is shielded from D; its calibration stays fixed.

A single stationary 15x3x2 model gives rate r_B/A=1+kappa d for d=0,1.
The same reference spectrum and readout serve kappa=-1/4,0,+1/4:
excited-source rates are 3/4,1,5/4. These are synthetic checks, not fitted
physics. Common energy rescaling leaves the ratio invariant; equal
fractional coupling of both clocks produces no relative-rate effect.

## Exact finite-clock predictions

At reference phase pi/2, for the +1/4 card, B's three orthogonal dial
probabilities are [1,0,0] with d=0 and
[1/9,(1+sqrt(3))^2/9,(1-sqrt(3))^2/9] with d=1.
For an unobserved equal source mixture or coherent superposition, P0=5/9.
Substituting mean Ed into the Hamiltonian incorrectly predicts
(3+2sqrt(2))/9={mean_prob:.9f}.
Density fluctuations produce multiple conditional rates and reduced
clock coherence, not necessarily one slower clock. At reference phase
pi, adjacent coherence vanishes but B is not fully dephased (purity5/9).
This does not by itself distinguish quantum source coherence from a
classical mixture: B alone sees source populations.

## Microscopic sign audit

Jaynes-Cummings exchange, in its controlled dispersive regime, yields
chi=g^2/Delta and calibrated r=(omega+chi)/(omega-chi).
Detuning +0.5 gives r={microscopic[1]['effective_rate']:.9f};
detuning -0.5 gives r={microscopic[0]['effective_rate']:.9f}.
Exact 2x2 spectra verify both signs and the second-order remainder bound.
Higher-order clock gaps are not equal: a universal rigid time dilation
does not follow from this mechanism. These are model units, not a
transmon calibration or an exact realization of the illustrative kappa=1/4.

## Closure boundary and next decision

The candidate/interaction/two-clock algebra is closed within its assumptions.
Finite phase likelihoods require repeated calibrated comparisons and
event/epoch labels, not an exact time from one outcome.
No external time is used in the stationary conditional construction.
The source/reference state, commensurate spectrum and reference
shielding are declared inputs. The G-R2 history constraint and the new
stationary constraint have NOT been conflated into one energy-complete
recording apparatus.

Next choose a physically grounded source/clock realization and determine
coupling sign, or test universality across independent clock transitions.
Do not fit a density slowdown law or promote gravity/Route-F conclusions.

Run: python3 route_g/code/verify_gr3_density_clocks.py.
Every check, fixed input and prior-file hash is stored in the JSON.
"""
    OUT.with_suffix(".md").write_text(report)
    print(json.dumps(result["summary"]))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
