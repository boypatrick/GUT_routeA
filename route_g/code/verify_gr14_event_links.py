#!/usr/bin/env python3
"""One event rule: clock drift, loop phase, visibility and memory dependence.

Fixed fixtures, not a fit. Exact reduced channels are checked against full
system-register unitary matrices. Every register branch is retained.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr14_event_links"
PI1 = np.diag([0., 1.])


def characteristic(p, x):
    return 1-p+p*np.exp(-1j*x)


def fractional_rate(p, x):
    """Continuous phase branch at x=0, used here only for |x|<pi.

    Zeros have no readable phase. Not a general phase-unwrapping routine.
    """
    if x == 0:
        return float(p)
    if not abs(x) < np.pi:
        raise ValueError("Use the connected phase branch |x|<pi")
    z = characteristic(p, x)
    return None if abs(z) < 1e-12 else float(-np.angle(z)/x)


def kron_all(items):
    result = np.array([1.], complex)
    for item in items:
        result = np.kron(result, item)
    return result


def full_protocol(charges, p, a, n, reuse=False):
    """Independent dense unitary construction, system first, registers last."""
    q = np.diag(charges)
    d = len(charges)
    nr = 1 if reuse else n
    register = np.array([np.sqrt(1-p), np.sqrt(p)], complex)
    initial = np.kron(np.ones(d)/np.sqrt(d), kron_all([register]*nr))
    unitary = np.eye(d*2**nr, dtype=complex)
    for step in range(n):
        factors = [np.eye(2)]*nr
        factors[0 if reuse else step] = PI1
        gate = expm(-1j*a*np.kron(q, kron_all(factors)))
        unitary = gate@unitary
    psi = unitary@initial
    amplitudes = psi.reshape(d, 2**nr)
    reduced = amplitudes@amplitudes.conj().T
    return reduced, unitary, initial, psi


def channel(charges, p, a, rho):
    v = expm(-1j*a*np.diag(charges))
    return (1-p)*rho+p*v@rho@v.conj().T


def pack_complex(value):
    return {"real": float(np.real(value)), "imag": float(np.imag(value))}


def run():
    card = json.loads((ROOT/"data/EVENT_LINK_CARD.json").read_text())
    a = np.pi*card["a_over_pi"]
    qa, qb = card["clock_charges"]
    n = card["events_in_joint_test"]
    checks, rows, errors = [], [], []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, actual, expected, tol=3e-12):
        error = float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))
        if tol == 3e-12:
            errors.append(error)
        check(name, np.isfinite(error) and error < tol)

    for path, digest in card["preserved_sha256"].items():
        check("preserved:"+path, hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest)
    close("primitive_three_cycle", np.exp(-3j*a), 1)

    for p in card["p_fixtures"]:
        tag = f"p={p:g}"
        z = characteristic(p, a)
        r, w = abs(z), z**3
        states = np.array([[np.sqrt(1-p), np.sqrt(p)*np.exp(-1j*j*a)] for j in range(3)])
        gram = states.conj()@states.T
        close(tag+":Gram_diagonal", np.diag(gram), np.ones(3))
        close(tag+":three_equal_oriented_links", [gram[0, 1], gram[1, 2], gram[2, 0]], [z]*3)
        check(tag+":Gram_positive", np.linalg.eigvalsh(gram).min() > -3e-12)
        close(tag+":rank_two_Gram_determinant", np.linalg.det(gram), 0)
        close(tag+":one_parameter_visibility", r*r, 1-3*p*(1-p))
        close(tag+":loop_Gram_constraint", 1-3*r*r+2*np.real(w), 0)
        phases = np.exp(1j*np.array([.37, -.81, 1.19]))
        gauged = phases.conj()[:, None]*gram*phases[None, :]
        close(tag+":closed_loop_rephasing", gauged[0, 1]*gauged[1, 2]*gauged[2, 0], w)
        close(tag+":orientation_reversal", gram[1, 0]*gram[2, 1]*gram[0, 2], w.conjugate())

        # CPTP check on two clocks; the SAME event register couples to both.
        charges = np.array([0., qb, qa, qa+qb])
        unit = expm(-1j*a*np.diag(charges))
        kraus = [np.sqrt(1-p)*np.eye(4), np.sqrt(p)*unit]
        close(tag+":Kraus_trace_preserving", sum(k.conj().T@k for k in kraus), np.eye(4))
        choi = sum(np.outer(k.ravel(), k.ravel().conj()) for k in kraus)
        check(tag+":Choi_positive", np.linalg.eigvalsh(choi).min() > -3e-12)
        rho, u, initial, psi = full_protocol(charges, p, a, n)
        close(tag+":full_unitarity", u.conj().T@u, np.eye(len(u)))
        close(tag+":joint_normalized", np.vdot(psi, psi), 1)
        occupancy = np.array([bin(j).count("1") for j in range(2**n)])
        probabilities = np.sum(abs(psi.reshape(4, 2**n))**2, axis=0)
        close(tag+":mean_active_events", probabilities@occupancy, n*p)
        close(tag+":active_event_variance", probabilities@((occupancy-n*p)**2), n*p*(1-p))
        charge_full = np.kron(np.diag(charges), np.eye(2**n))
        close(tag+":bare_charge_conservation", u.conj().T@charge_full@u, charge_full)
        close(tag+":all_branches_trace", np.trace(rho), 1)
        predicted = characteristic(p, a*(charges[:, None]-charges[None, :]))**n/4
        close(tag+":full_vs_characteristic_matrix", rho, predicted)
        iterated = np.ones((4, 4), complex)/4
        for _ in range(n):
            iterated = channel(charges, p, a, iterated)
        close(tag+":full_vs_iterated_Kraus", rho, iterated)
        tensor = rho.reshape(2, 2, 2, 2)
        ra, rb = np.trace(tensor, axis1=1, axis2=3), np.trace(tensor, axis1=0, axis2=2)
        close(tag+":clock_A_marginal", 2*ra[1, 0], characteristic(p, a*qa)**n)
        close(tag+":clock_B_marginal", 2*rb[1, 0], characteristic(p, a*qb)**n)
        close(tag+":relative_clock_coherence", 4*rho[2, 1], characteristic(p, a*(qa-qb))**n)
        close(tag+":global_reversibility", u.conj().T@psi, initial)
        # Commuting random product unitaries cannot create A-B entanglement
        # from this product input. This also checks the physical claim boundary.
        partial_transpose = tensor.transpose(0, 3, 2, 1).reshape(4, 4)
        check(tag+":two_clock_PPT", np.linalg.eigvalsh(partial_transpose).min() > -3e-12)

        # Actual loop interferometer: balanced reference/loop paths, no filters.
        loop, _, _, _ = full_protocol([0., 1.], p, a, 3)
        close(tag+":loop_full_unitary_coherence", 2*loop[1, 0], w)
        close(tag+":loop_equal_populations_no_postselection", np.diag(loop), [.5, .5])
        for angle in [0., np.pi/2, np.pi, 3*np.pi/2]:
            plus = np.array([1, np.exp(1j*angle)])/np.sqrt(2)
            prob = np.vdot(plus, loop@plus).real
            close(tag+f":fringe:{angle:.3f}", prob, (1+np.real(np.exp(-1j*angle)*w))/2)
            check(tag+f":both_outputs_physical:{angle:.3f}", 0 <= prob <= 1+1e-12)
        reused, _, _, _ = full_protocol([0., 1.], p, a, 3, reuse=True)
        close(tag+":reuse_same_register_loop_is_identity", 2*reused[1, 0], 1)
        # Joint clock/path response is NOT the product of the two marginals.
        joint, _, _, _ = full_protocol([0., 1., qa, qa+1], p, a, 3)
        cjoint = characteristic(p, a*(qa+1))**3
        cmargproduct = characteristic(p, a*qa)**3*w
        close(tag+":joint_clock_loop_prediction", 4*joint[3, 0], cjoint)
        pauli_x = np.array([[0., 1.], [1., 0.]])
        pauli_y = np.array([[0., -1j], [1j, 0.]])
        close(tag+":joint_real_quadrature", np.trace(joint@(np.kron(pauli_x, pauli_x)
              -np.kron(pauli_y, pauli_y))), cjoint.real)
        close(tag+":joint_imaginary_quadrature", np.trace(joint@(np.kron(pauli_y, pauli_x)
              +np.kron(pauli_x, pauli_y))), cjoint.imag)
        if 0 < p < 1:
            check(tag+":joint_is_not_independent_noise", abs(cjoint-cmargproduct) > .01)
            check(tag+":fresh_and_reused_loop_differ", abs(w-1) > .1)
        equal, _, _, _ = full_protocol([0., qa, qa, 2*qa], p, a, 3)
        close(tag+":equal_clock_relative_channel_protected", 4*equal[2, 1], 1)
        clock_rows = [{"q": q, "mean_event_fraction": p,
                       "fractional_phase_rate": fractional_rate(p, a*q),
                       "visibility_after_three": float(abs(characteristic(p, a*q))**3)}
                      for q in [qa, qb]]
        rows.append({"p": p, "r_link": float(r), "loop_W": pack_complex(w),
                     "loop_phase_principal": float(np.angle(w)), "loop_visibility": float(abs(w)),
                     "clocks": clock_rows, "joint_clock_loop": pack_complex(cjoint),
                     "incorrect_product_of_marginals": pack_complex(cmargproduct),
                     "reuse_loop_W": pack_complex(2*reused[1, 0])})

    # Exact balanced rule gives common fractional rate on the connected band.
    for x in [.1, .5, 1., 2., 3.]:
        close(f"balanced:factorization:{x}", characteristic(.5, x), np.exp(-.5j*x)*np.cos(x/2))
        close(f"balanced:exact_common_rate:{x}", fractional_rate(.5, x), .5)
    close("balanced:unreadable_phase_zero", characteristic(.5, np.pi), 0)
    close("balanced:loop_minus_one_eighth", characteristic(.5, a)**3, -.125)
    check("generic_rate_not_universal", abs(fractional_rate(.25, a*qa)-fractional_rate(.25, a*qb)) > .01)

    # Cumulant expansion predicts the error without a fit. Halve the numerical
    # argument, not a physical parameter search. Error orders: rate x^2, next x^4.
    p = .25
    c3 = p*(1-p)*(1-2*p)/6
    small = []
    for x in [.08, .04, .02]:
        actual = fractional_rate(p, x)
        variance = -np.log(abs(characteristic(p, x)))
        close(f"small_x:quadratic_rate_coefficient:{x}", (p-actual)/x**2, c3, 1e-5)
        close(f"small_x:dephasing_coefficient:{x}", variance/x**2, p*(1-p)/2, 1e-5)
        small.append({"x": x, "rate": actual, "quadratic_residual": abs(actual-(p-c3*x*x))})
    check("small_x:fourth_order_rate_remainder", 14 < small[0]["quadratic_residual"]/small[1]["quadratic_residual"] < 18)

    failed = [c["name"] for c in checks if not c["passed"]]
    result = {"id": "G-R14", "summary": {"total": len(checks), "passed": len(checks)-len(failed), "failed": failed},
              "max_exact_matrix_or_identity_error": max(errors) if errors else 0,
              "fixtures": rows, "small_argument_checks": small, "checks": checks,
              "scope": card["status"], "not_claimed": card["not_claimed"]}
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# G-R14 Event links from one binary register rule", "",
             f"{result['summary']['passed']}/{len(checks)} checks pass. Fixed fixtures; no physical fit.", "",
             "U(Q)=exp(-i a Q tensor |1><1|), a=2 pi/3; fresh registers have active probability p.",
             "Every response comes from F_p(x)=1-p+p exp(-i x); shared events use the SUM of probe charges.", "",
             "| p | mean progress N | phase rate A / bare A | phase rate B / bare B | loop phase (rad) | loop visibility |",
             "|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        ca, cb = row["clocks"]
        lines.append(f"| {row['p']:.2f} | {row['p']:.2f} | {ca['fractional_phase_rate']:.9f} | {cb['fractional_phase_rate']:.9f} | {row['loop_phase_principal']:.9f} | {row['loop_visibility']:.9f} |")
    lines += ["", "## Exact constraints and the universality boundary", "",
              "With z=F_p(2 pi/3), r=|z| and W=z^3=r^3 exp(i Phi),",
              "r^2=1-3p(1-p) and 1-3r^2+2r^3 cos(Phi)=0.",
              "Clock coherence is F_p(a q)^n, not an independently chosen damping function.",
              "A common mean progress p does not generally give a common fractional phase rate.",
              "For all sufficiently small q it is exact only at p=0, 1/2 or 1; at p=1/2 it holds",
              "on |a q|<pi, with visibility |cos(a q/2)|^n. Phase is undefined at zeros.", "",
              "## The memory counterprotocol and joint prediction", "",
              "Reusing one record on three links gives F_p(3a)=1, instead of F_p(a)^3.",
              "At p=1/2, fresh records give phase pi and visibility 1/8; the reused record gives phase 0 and visibility 1.",
              "The marginals of one event are identical. The inter-event memory rule changes, not a fitted phase.",
              "A simultaneous clock/loop coherence is F_p(a(q+1))^3, not F_p(a q)^3 F_p(a)^3.",
              "Equal-frequency clocks sharing the same records retain unit relative coherence even when their individual dials fade.", "",
              "## Interpretation and stopping boundary", "",
              "This is a bounded quantum collision model inspired by event comparison, not a proof of the source theory.",
              "The three links are not spatial dimensions. The reduced loop phase is operational interference, not a Maxwell field.",
              "Update ordering, register supply, universal coupling, probe charges, quantum probabilities and readout remain inputs.",
              "Tracing records is reversible loss of access; no irreversible decoherence or autonomous event generation is proved.",
              "For the commuting observables tested, classical correlated binary noise reproduces the reduced predictions.",
              "Next constrain memory correlations or require cross-probe agreement without independent phase/noise fits.",
              "G-R13 controller calibration remains an implementation side branch, not a prerequisite.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(result["summary"], indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
