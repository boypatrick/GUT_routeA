#!/usr/bin/env python3
"""One fixed C1 relay pulse, with full 21-state energy and controller work.

Same peak coupling and leading preparation area as the G-R12 square pulse.
Only solver/quad refinements are varied; no physical pulse search is performed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.linalg import expm

from verify_gr12_correlated_preparation import local_matrices, relays
from verify_gr10_relational_readout import DPLUS, DMINUS, proj
from verify_gr11_preparation_record import make_instrument

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr13_smooth_switching"


def envelope(t, duration):
    if t <= 0 or t >= duration:
        return 0., 0.
    return np.sin(np.pi*t/duration)**2, np.pi/duration*np.sin(2*np.pi*t/duration)


def cf4_propagate(hamiltonian, initial, duration, steps):
    """Independent unitary fourth-order commutator-free Magnus propagation."""
    a1, a2 = (3-2*np.sqrt(3))/12, (3+2*np.sqrt(3))/12
    c1, c2 = .5-np.sqrt(3)/6, .5+np.sqrt(3)/6
    dt = duration/steps
    state = initial.copy()
    for j in range(steps):
        h1, h2 = hamiltonian((j+c1)*dt), hamiltonian((j+c2)*dt)
        # Earlier-weighted exponential acts first (on the right).
        for h in (a2*h1+a1*h2, a1*h1+a2*h2):
            e, v = np.linalg.eigh(h)
            state = v@(np.exp(-1j*dt*e)*(v.conj().T@state))
    return state


def expectation(state, operator):
    return float(np.vdot(state, operator@state).real)


def run():
    card = json.loads((ROOT/"data/SMOOTH_SWITCHING_CARD.json").read_text())
    parent = json.loads((ROOT/card["parent_card"]).read_text())
    previous = json.loads((ROOT/card["parent_result"]).read_text())
    gamma, delta, g = card["gamma"], card["delta"], card["g_peak"]
    tau = previous["duration"]
    duration = 8*tau/3
    ff, _, density, current, h0, vertex, _ = local_matrices()
    _, b = relays(current, g, delta)
    h = h0+gamma*vertex
    z = np.zeros((9, 12))
    hoff = np.block([[h, z], [z.T, (4+delta)*np.eye(12)]])
    relay = np.block([[np.zeros((9, 9)), b.conj().T], [b, np.zeros((12, 12))]])
    compensation = np.zeros((21, 21), complex)
    compensation[:9, :9] = np.sqrt(2)*g*density
    initial = np.r_[ff[:, 4], np.zeros(12)]
    target, odd = ff@DPLUS, ff@DMINUS
    shifted = hoff-4*np.eye(21)  # Common phase removed, not an energy omission.
    checks = []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, atol=2e-9):
        check(name, np.allclose(a, b, atol=atol, rtol=2e-9))

    for path, digest in card["preserved_sha256"].items():
        check("preserved:"+path, hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest)
    close("fixed_gamma", gamma, parent["mediator_gamma"])
    close("fixed_g", g, previous["g"])
    close("fixed_Delta", delta, 80)
    check("three_numerical_CF4_refinements", len(card["solver_checks"]["cf4_steps"]) == 3)
    close("same_full_mediator_identity", -b.conj().T@b/delta,
          g*current-np.sqrt(2)*g*(np.eye(9)+density))
    close("compensated_leading_operator", -b.conj().T@b/delta+compensation[:9, :9],
          g*current-np.sqrt(2)*g*np.eye(9))
    close("relay_Hermitian", relay, relay.conj().T)
    close("compensation_Hermitian", compensation, compensation.conj().T)
    for t in [-1., 0., duration, duration+1]:
        close(f"off_and_zero_slope:{t}", envelope(t, duration), [0., 0.])
    close("same_peak", envelope(duration/2, duration)[0], 1)
    close("same_effective_area_analytic", g*duration*3/8, g*tau)
    close("duration_tradeoff", duration/tau, 8/3)
    check("not_rounded_to_an_update", abs(duration/.25-round(duration/.25)) > .1)

    def hamiltonian(t):
        f, _ = envelope(t, duration)
        return shifted+f*relay+f*f*compensation

    def rhs(t, y):
        psi = y[:21]
        f, df = envelope(t, duration)
        pb = df*expectation(psi, relay)
        pc = 2*f*df*expectation(psi, compensation)
        return np.r_[-1j*hamiltonian(t)@psi, pb, pc]

    settings = card["solver_checks"]
    initial_aug = np.r_[initial, 0j, 0j]
    solutions = []
    for name, rtol, atol, step in [
        ("primary", settings["dop853_rtol"], settings["dop853_atol"], .02),
        ("refined", settings["refined_rtol"], settings["refined_atol"], .01),
    ]:
        sol = solve_ivp(rhs, (0, duration), initial_aug, method="DOP853",
                        rtol=rtol, atol=atol, max_step=step, dense_output=True)
        check(name+":solver_completed", sol.success)
        close(name+":unit_norm", np.vdot(sol.y[:21, -1], sol.y[:21, -1]), 1)
        solutions.append(sol)
    sol = solutions[-1]
    out = sol.y[:21, -1]
    work_components = sol.y[21:, -1].real
    ode_error = float(np.linalg.norm(solutions[0].y[:21, -1]-out))
    check("tolerance_refinement_state", ode_error < 1e-10)
    close("tolerance_refinement_work", solutions[0].y[21:, -1], sol.y[21:, -1], 1e-10)
    f_target = float(abs(np.vdot(target, out[:9]))**2)
    paux = float(np.linalg.norm(out[9:])**2)
    eaux = (4+delta)*paux
    ein = expectation(initial, hoff)
    eout = expectation(out, hoff)
    close("net_work_equals_endpoint_energy", work_components.sum(), eout-ein, 1e-10)
    close("all_probabilities_retained", np.linalg.norm(out[:9])**2+paux, 1)
    close("odd_clock_not_populated", np.vdot(odd, out[:9]), 0)

    cf_rows = []
    for steps in settings["cf4_steps"]:
        other = cf4_propagate(hamiltonian, initial, duration, steps)
        err = float(np.linalg.norm(other-out))
        pa = float(np.linalg.norm(other[9:])**2)
        close(f"CF4:{steps}:norm", np.vdot(other, other), 1, 1e-10)
        check(f"CF4:{steps}:residual_resolved", abs(pa-paux)/paux < .02)
        cf_rows.append({"steps": steps, "state_error_against_refined_ODE": err,
                        "auxiliary_probability": pa,
                        "fidelity": float(abs(np.vdot(target, other[:9]))**2)})
    check("CF4_fourth_order_convergence", all(cf_rows[i]["state_error_against_refined_ODE"]/
          cf_rows[i+1]["state_error_against_refined_ODE"] > 12 for i in range(2)))
    check("CF4_final_state_agreement", cf_rows[-1]["state_error_against_refined_ODE"] < 5e-11)

    quadratures = []
    for count in settings["quadrature_points"]:
        ts = np.linspace(0, duration, count)
        ys = sol.sol(ts)
        states = ys[:21]
        fs = np.sin(np.pi*ts/duration)**2
        dfs = np.pi/duration*np.sin(2*np.pi*ts/duration)
        eb = np.einsum("it,ij,jt->t", states.conj(), relay, states).real
        ec = np.einsum("it,ij,jt->t", states.conj(), compensation, states).real
        pb, pc = dfs*eb, 2*fs*dfs*ec
        pt = pb+pc
        result = {"points": count,
                  "work_relay": float(simpson(pb, x=ts)),
                  "work_compensation": float(simpson(pc, x=ts)),
                  "positive_common_bus": float(simpson(np.maximum(pt, 0), x=ts)),
                  "returned_common_bus": float(simpson(np.maximum(-pt, 0), x=ts)),
                  "positive_separate_ports": float(simpson(np.maximum(pb, 0)+np.maximum(pc, 0), x=ts)),
                  "returned_separate_ports": float(simpson(np.maximum(-pb, 0)+np.maximum(-pc, 0), x=ts))}
        quadratures.append(result)
        close(f"quad:{count}:relay_work", result["work_relay"], work_components[0], 1e-10)
        close(f"quad:{count}:compensation_work", result["work_compensation"], work_components[1], 1e-10)
        close(f"quad:{count}:same_area", simpson(g*fs**2, x=ts), g*tau)
        instantaneous_energy = (np.einsum("it,ij,jt->t", states.conj(), hoff, states).real
                                +fs*eb+fs**2*ec)
        check(f"quad:{count}:first_law_along_whole_path",
              np.max(abs(instantaneous_energy-ein-ys[21:].real.sum(axis=0))) < 1e-9)
    for key in ["positive_common_bus", "returned_common_bus", "positive_separate_ports", "returned_separate_ports"]:
        close("quadrature_refinement:"+key, quadratures[0][key], quadratures[1][key], 2e-8)
    quad = quadratures[-1]
    net = float(work_components.sum())
    close("positive_returned_net", quad["positive_separate_ports"]-quad["returned_separate_ports"], net)
    close("common_bus_net", quad["positive_common_bus"]-quad["returned_common_bus"], net)
    transient_aux = np.sum(abs(states[9:])**2, axis=0)
    check("transient_mediator_active_not_removed", transient_aux.max() > .001)

    # Reconstruct the one exact square comparator, without changing any parameters.
    square = expm(-1j*tau*(shifted+relay+compensation))@initial
    square_ports = np.array([[expectation(initial, relay), expectation(initial, compensation)],
                             [-expectation(square, relay), -expectation(square, compensation)]])
    sq_positive = float(np.maximum(square_ports, 0).sum())
    sq_returned = float(np.maximum(-square_ports, 0).sum())
    square_row = next(r for r in previous["mediator"] if r["compensated"] and r["delta"] == delta)
    close("square_reference_fidelity", abs(np.vdot(target, square[:9]))**2, square_row["fidelity_unconditional"])
    close("square_reference_auxiliary", np.linalg.norm(square[9:])**2, square_row["auxiliary_probability"])
    close("square_reference_work", square_ports.sum(), square_row["net_work_on_full_model"])
    close("square_complete_work", square_ports.sum(), expectation(square, hoff)-ein)
    check("fidelity_maintained_and_improved", f_target > square_row["fidelity_unconditional"])
    check("endpoint_residue_lower", paux < square_row["auxiliary_probability"])
    check("auxiliary_energy_lower", eaux < square_row["auxiliary_energy"])
    check("net_work_lower", net < square_row["net_work_on_full_model"])
    check("positive_input_separate_ports_lower", quad["positive_separate_ports"] < sq_positive)
    # Important countercheck: not EVERY resource measure is lower.
    check("absolute_work_not_claimed_lower", quad["positive_separate_ports"]+quad["returned_separate_ports"] > sq_positive+sq_returned)
    check("duration_not_claimed_equal", duration > tau)

    _, _, _, ks_mode, *_ = make_instrument()
    ks = [ff@a@ff.conj().T for a in ks_mode]
    record = [float(np.linalg.norm(ks[0]@out[:9])**2+paux)]
    record += [float(np.linalg.norm(a@out[:9])**2) for a in ks[1:]]
    close("GR11_record_all_outcomes", sum(record), 1)
    close("GR11_success_is_unconditional_fidelity", sum(record[1:]), f_target)
    success_state = sum(proj(a@out[:9]) for a in ks[1:])/f_target
    close("GR11_backreaction_not_removed", success_state, (57/82)*proj(target)+(25/82)*proj(odd))
    for wait in [1.7, 13.]:
        future = expm(-1j*wait*shifted)@out
        close(f"post_pulse:{wait}:fidelity_preserved", abs(np.vdot(target, future[:9]))**2, f_target)
        close(f"post_pulse:{wait}:aux_population_preserved", np.linalg.norm(future[9:])**2, paux)

    # Conditional overhead threshold, not a measured controller model or a scan.
    extra_duration = duration-tau
    threshold_zero = (sq_positive-quad["positive_separate_ports"])/extra_duration
    threshold_recovery_coefficient = (quad["returned_separate_ports"]-sq_returned)/extra_duration
    check("nonrecovering_ideal_input_improves", threshold_zero > 0)
    snapshots = []
    for fraction in [0., .25, .5, .75, 1.]:
        t = fraction*duration
        y = sol.sol(t)
        snapshots.append({"fraction": fraction, "f": envelope(t, duration)[0],
                          "fidelity": float(abs(np.vdot(target, y[:9]))**2),
                          "auxiliary_probability": float(np.linalg.norm(y[9:21])**2),
                          "work_relay": float(y[21].real), "work_compensation": float(y[22].real)})
    failed = [x["name"] for x in checks if not x["passed"]]
    summary = {"checks": len(checks), "passed": len(checks)-len(failed), "failed": failed}
    result = {"stage": "G-R13", "summary": summary, "fixed_parameters": {"gamma": gamma, "delta": delta, "g": g},
              "square_duration": tau, "smooth_duration": duration, "effective_area": g*tau,
              "square": {**square_row, "positive_work_separate_ports": sq_positive,
                         "returned_work_separate_ports": sq_returned, "component_quenches": square_ports.tolist()},
              "smooth": {"fidelity_unconditional": f_target, "auxiliary_probability": paux,
                         "auxiliary_energy": eaux, "auxiliary_excess_energy": delta*paux,
                         "work_relay": float(work_components[0]), "work_compensation": float(work_components[1]),
                         "net_work": net, "energy_change": eout-ein, **quad,
                         "max_sampled_transient_aux_probability": float(transient_aux.max())},
              "overhead_threshold_model_power": {"eta_zero_intercept": threshold_zero,
                                                 "eta_coefficient": threshold_recovery_coefficient},
              "solver_validation": {"ODE_state_difference": ode_error, "CF4": cf_rows, "quadratures": quadratures},
              "snapshots": snapshots, "checks": checks, "limits": card["not_claimed"]}
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# G-R13 one smooth pulse on the fixed 21-state mediator", "",
             f"{summary['passed']}/{summary['checks']} checks pass. One physical pulse, no optimization.", "",
             "Fix gamma=.03, Delta=80 and g_peak=.1. The relay envelope is f=sin²(pi*s/T);",
             "the compensation is f²*Vc0, not f*Vc0. Peak coupling and effective preparation",
             f"area are unchanged. Duration increases from {tau:g} to {duration:.10g}, a factor 8/3.",
             "The continuous control parameter is assumed, not derived time or a rounded update count.", "",
             "| quantity | compensated square | one smooth pulse |", "|---|---:|---:|"]
    for name, left, right in [
        ("unconditional preparation fidelity", square_row["fidelity_unconditional"], f_target),
        ("endpoint auxiliary probability", square_row["auxiliary_probability"], paux),
        ("endpoint auxiliary energy", square_row["auxiliary_energy"], eaux),
        ("net controller work on full model", square_row["net_work_on_full_model"], net),
        ("positive input, separate control ports", sq_positive, quad["positive_separate_ports"]),
        ("returned work, separate control ports", sq_returned, quad["returned_separate_ports"]),
        ("absolute work, separate control ports", sq_positive+sq_returned, quad["positive_separate_ports"]+quad["returned_separate_ports"]),
        ("duration", tau, duration),
    ]:
        lines.append(f"| {name} | {left:.10g} | {right:.10g} |")
    lines += ["", "## Complete energy ledger", "",
              f"Integrated relay work = {work_components[0]:.12g}; compensation work = {work_components[1]:.12g}.",
              "Their sum matches the endpoint energy change, with the first law also checked",
              "throughout the pulse. Compensation is not free or silently cancelled from the ledger.",
              "Negative net work means energy is returned from this prepared input; it does not",
              "mean zero input work, free initial-state preparation, zero dissipation or a free controller.",
              f"The sampled transient auxiliary probability reaches {transient_aux.max():.6g},",
              "so the mediator remains active even though the endpoint residue is suppressed.", "",
              "## What improved and what did not", "",
              "For this single fixed case, fidelity, endpoint residue and ideal positive input work",
              "all improve. Duration and absolute exchanged work increase. Actual hardware energy",
              "consumption has not been supplied. If both protocols share return efficiency eta and",
              "constant overhead power p0, define cost=W_plus-eta*W_minus+p0*duration using separate ports.",
              f"The smooth pulse costs less only if p0 < {threshold_zero:.9g} + {threshold_recovery_coefficient:.9g}*eta",
              "in model power units. No eta, p0 or physical scale is fitted or inferred.", "",
              "## Verification and stopping boundary", "",
              "Adaptive DOP853 tolerance refinement, independent unitary CF4 step refinement,",
              "and power quadrature refinement test the same waveform. The tiny endpoint residue",
              "is nonzero and numerically resolved, not a universal scaling law or a noise tolerance.",
              "The unchanged G-R11 reader keeps every outcome and its 57/82 versus 25/82 disturbance.",
              "This one-pulse question is bounded-done. Do not begin pulse optimization. A real",
              "controller's losses/overhead or an explicit autonomous trigger are separate physical",
              "choices; no time-arrow, Ed or Route-F claim is promoted.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2))
    print(json.dumps(result["smooth"], indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
