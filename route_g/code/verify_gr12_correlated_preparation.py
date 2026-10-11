#!/usr/bin/env python3
"""Local conditional hopping, full-contact preparation and finite relay ledger.

No fitting. The new hopping, relay phases and pulse control are explicit inputs.
Auxiliary exits remain counted; small population is not negligible energy.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from verify_gr10_relational_readout import DPLUS, DMINUS, proj
from verify_gr11_preparation_record import make_instrument

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr12_correlated_preparation"


def local_matrices():
    eye = np.eye(3)
    shift = np.roll(eye, 1, axis=0)
    f = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
    ff = np.kron(f, f)
    lzero = 2*eye-shift-shift.T
    k = np.kron(lzero, eye)+np.kron(eye, lzero)
    density = np.diag([float(j == l) for j in range(3) for l in range(3)])
    current = 1j*(k@density-density@k)/np.sqrt(2)
    l = 2*eye-np.exp(1j*np.pi/6)*shift-np.exp(-1j*np.pi/6)*shift.T
    h0 = np.kron(l, eye)+np.kron(eye, l)
    contact = np.diag(np.eye(9)[0])
    return ff, k, density, current, h0, contact, np.kron(shift, shift)


def relays(current, g, delta):
    edges = [(x, y) for x in range(9) for y in range(x+1, 9)
             if abs(current[x, y]) > 1e-12]
    omega = np.sqrt(g*delta/np.sqrt(2))
    b = np.zeros((len(edges), 9), complex)
    for e, (x, y) in enumerate(edges):
        b[e, x] = omega
        b[e, y] = -g*delta*current[x, y]/omega
    return edges, b


def evolve(h, psi, s):
    return expm(-1j*s*h)@psi


def run():
    card = json.loads((ROOT/"data/CORRELATED_PREPARATION_CARD.json").read_text())
    parent = json.loads((ROOT/card["parent_card"]).read_text())
    g, alpha = parent["exchange_g"], parent["alpha_per_assumed_update"]
    beta = float(np.arcsin(1/np.sqrt(3)))
    n = round(beta/(g*alpha))
    duration = n*alpha
    ff, k, d, c, h0, vertex, translation = local_matrices()
    eye = np.eye(9)
    m, symmetric = eye[:, 4], (eye[:, 2]+eye[:, 6])/np.sqrt(2)
    initial, target, odd = ff@m, ff@DPLUS, ff@DMINUS
    triple = ff[:, [2, 4, 6]]
    p = triple@triple.conj().T
    swap = eye.reshape(3, 3, 3, 3).transpose(1, 0, 2, 3).reshape(9, 9)
    _, pd_mode, _, ks_mode, *_ = make_instrument()
    pd = ff@pd_mode@ff.conj().T
    ks = [ff@a@ff.conj().T for a in ks_mode]
    checks, direct_rows, mediator_rows = [], [], []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, atol=5e-10):
        check(name, np.allclose(a, b, atol=atol, rtol=5e-10))

    for path, sha in card["preserved_sha256"].items():
        check("preserved:"+path, hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha)
    close("density_is_projector", d@d, d)
    close("current_Hermitian", c, c.conj().T)
    close("joint_translation", translation@c, c@translation)
    close("swap_symmetry", swap@c, c@swap)
    close("triple_closed_for_added_current", (eye-p)@c@p, np.zeros((9, 9)))
    cmode = ff.conj().T@c@ff
    pmode = ff.conj().T@p@ff
    desired = 1j*(np.outer(m, symmetric)-np.outer(symmetric, m))
    close("exact_preparation_block", pmode@cmode@pmode, desired)
    close("old_resonant_kinetic_commutator_cannot_replace_K",
          p@(h0@d-d@h0)@p, np.zeros((9, 9)))
    check("not_density_only", np.linalg.norm(c-np.diag(np.diag(c))) > 1)
    cycle = np.roll(np.eye(3), 1, axis=0)
    block = triple.conj().T@c@triple
    check("translation_does_not_imply_circulant", np.linalg.norm(cycle@block-block@cycle) > 1)
    close("odd_dark_annihilated", c@odd, np.zeros(9))
    for x in range(9):
        for y in range(x+1, 9):
            if abs(c[x, y]) > 1e-12:
                a, b = divmod(x, 3); ap, bp = divmod(y, 3)
                check(f"local_support:{x}:{y}", (a == ap) != (b == bp))
                close(f"occupation_condition:{x}:{y}", abs(d[x, x]-d[y, y]), 1)
    for chi in [.37, np.pi/2]:
        u = expm(1j*chi*d)
        plus, minus = u@k@u.conj().T, u.conj().T@k@u
        close(f"finite_Peierls_identity:{chi}", -(plus-minus)/(2*np.sqrt(2)*np.sin(chi)), c)
    exact = evolve(h0+g*c, initial, beta/g)
    close("ideal_exact_state", proj(exact), proj(target))
    for area in [0, .3, beta, beta+.1, beta+np.pi/2]:
        state = evolve(h0+g*c, initial, area/g)
        close(f"ideal_area_response:{area}", abs(np.vdot(target, state))**2, np.cos(area-beta)**2)
    wrong = evolve(h0-g*c, initial, beta/g)
    close("wrong_phase_is_not_a_coordinate_change", abs(np.vdot(target, wrong))**2, 1/9)
    close("integer_update_rounding", g*duration-beta, .1*25*.25-beta)

    for gamma in [0.]+parent["local_contact_gamma_fixtures"]:
        h = h0+gamma*vertex
        hon = h+g*c
        out = evolve(hon, initial, duration)
        fidelity = abs(np.vdot(target, out))**2
        leakage = float(np.linalg.norm((eye-p)@out)**2)
        win = float(g*np.vdot(initial, c@initial).real)
        wout = float(-g*np.vdot(out, c@out).real)
        ein, ef = float(np.vdot(initial, h@initial).real), float(np.vdot(out, h@out).real)
        tag = f"direct:gamma={gamma}"
        close(tag+":norm", np.vdot(out, out), 1)
        check(tag+":improved_unconditional_preparation", fidelity > .99)
        close(tag+":hold_energy_conserved", np.vdot(out, hon@out), np.vdot(initial, hon@initial))
        close(tag+":on_work_zero", win, 0)
        close(tag+":complete_work_ledger", ef-ein, win+wout)
        close(tag+":odd_dark_remains_zero", np.vdot(odd, out), 0)
        close(tag+":dark_weight_is_fidelity", np.vdot(out, pd@out), fidelity)
        close(tag+":GR11_success_is_fidelity", sum(np.linalg.norm(a@out)**2 for a in ks[1:]), fidelity)
        success = sum(proj(a@out) for a in ks[1:])/fidelity
        close(tag+":GR11_backreaction_retained", success,
              (57/82)*proj(target)+(25/82)*proj(odd))
        for wait in [1.7, 13.]:
            future = evolve(h, out, wait)
            close(tag+f":fidelity_preserved_after_switch_off:{wait}", abs(np.vdot(target, future))**2, fidelity)
        continued = evolve(hon, initial, duration+np.pi/(2*g))
        check(tag+":continued_control_not_storage", abs(np.vdot(target, continued))**2 < .02)
        if gamma:
            check(tag+":not_full_bare_energy_conserving", np.linalg.norm(h@c-c@h) > .001)
            close(tag+":old_invariant_unitary_bound", np.vdot(initial, pd@initial), 2/3)
        direct_rows.append({"gamma": gamma, "fidelity_unconditional": float(fidelity),
                            "outside_triple_probability": leakage,
                            "initial_energy": ein, "final_energy": ef,
                            "work_on": win, "work_off": wout, "net_work_on_system": win+wout})

    # An explicit finite auxiliary completion. No auxiliary outcome is discarded.
    gamma = card["mediator_gamma"]
    h = h0+gamma*vertex
    low_identity = np.eye(9)
    for compensated in [False, True]:
        compensation = np.sqrt(2)*g*d if compensated else np.zeros((9, 9))
        heff = h+g*c-np.sqrt(2)*g*(low_identity+d)+compensation
        effective_f = float(abs(np.vdot(target, evolve(heff, initial, duration)))**2)
        if compensated:
            close("compensated_effective_matches_direct", effective_f,
                  next(r["fidelity_unconditional"] for r in direct_rows if r["gamma"] == gamma))
        for delta in card["mediator_detunings"]:
            tag = f"mediator:comp={compensated}:Delta={delta}"
            edges, b = relays(c, g, delta)
            auxdim = len(edges)
            check(tag+":twelve_auxiliary_edges", auxdim == 12)
            close(tag+":SW_hopping_AND_Stark", -b.conj().T@b/delta,
                  g*c-np.sqrt(2)*g*(low_identity+d))
            counts = np.count_nonzero(abs(b) > 1e-12, axis=0)
            close(tag+":degree_controls_Stark", counts, 2+2*np.diag(d))
            z = np.zeros((9, auxdim))
            hoff = np.block([[h, z], [z.T, (4+delta)*np.eye(auxdim)]])
            hon = np.block([[h+compensation, b.conj().T], [b, (4+delta)*np.eye(auxdim)]])
            initial_full, target_full = np.r_[initial, np.zeros(auxdim)], np.r_[target, np.zeros(auxdim)]
            out = evolve(hon, initial_full, duration)
            e, v = np.linalg.eigh(hon)
            independent = (v*np.exp(-1j*duration*e))@v.conj().T@initial_full
            close(tag+":exponential_crosscheck", out, independent)
            close(tag+":unit_norm", np.vdot(out, out), 1)
            close(tag+":hold_energy", np.vdot(out, hon@out), np.vdot(initial_full, hon@initial_full))
            f = float(abs(np.vdot(target_full, out))**2)
            auxp = float(np.vdot(out[9:], out[9:]).real)
            systemp = float(np.vdot(out[:9], out[:9]).real)
            close(tag+":all_probabilities", systemp+auxp, 1)
            close(tag+":odd_sector_zero", np.vdot(odd, out[:9]), 0)
            check(tag+":unconditional_not_postselected", f <= systemp+1e-10)
            pulse = hon-hoff
            won = float(np.vdot(initial_full, pulse@initial_full).real)
            woff = float(-np.vdot(out, pulse@out).real)
            en_initial = float(np.vdot(initial_full, hoff@initial_full).real)
            en_system = float(np.vdot(out[:9], h@out[:9]).real)
            en_aux = (4+delta)*auxp
            close(tag+":work_ledger_including_aux", en_system+en_aux-en_initial, won+woff)
            close(tag+":on_work_analytic", won, np.sqrt(2)*g/3 if compensated else 0)
            # Extend the existing reader: auxiliary population is a failure,
            # never a successful-system renormalization.
            record = [np.linalg.norm(ks[0]@out[:9])**2+auxp]
            record += [np.linalg.norm(a@out[:9])**2 for a in ks[1:]]
            close(tag+":record_complete_with_aux_failure", sum(record), 1)
            close(tag+":readout_success_unconditional", sum(record[1:]), f)
            check(tag+":aux_energy_not_neglected", en_aux > .03)
            mediator_rows.append({"compensated": compensated, "delta": delta,
                                  "effective_reference_fidelity": effective_f,
                                  "fidelity_unconditional": f, "auxiliary_probability": auxp,
                                  "system_probability": systemp, "auxiliary_energy": en_aux,
                                  "system_energy_unnormalized": en_system,
                                  "work_on": won, "work_off": woff,
                                  "net_work_on_full_model": won+woff})
    failed = [x["name"] for x in checks if not x["passed"]]
    summary = {"checks": len(checks), "passed": len(checks)-len(failed), "failed": failed}
    data = {"stage": "G-R12", "summary": summary, "g": g, "beta": beta,
            "updates": n, "duration": duration, "ideal_duration": beta/g,
            "direct": direct_rows, "mediator": mediator_rows, "checks": checks,
            "limits": card["not_claimed"]}
    OUT.with_suffix(".json").write_text(json.dumps(data, indent=2)+"\n")
    lines = ["# G-R12 local correlated hopping and preparation energy ledger", "",
             f"{summary['passed']}/{summary['checks']} checks pass. No fitted pulse or physical units.", "",
             "The current C=i[K,D]/sqrt(2) is local on configuration-graph edges,",
             "joint-translation and swap invariant, but not a density-only interaction.",
             "It supplies the exact desired two-state preparation block without inserting",
             "a Fourier-mode projector. Its coefficient/phase and control remain new inputs.", "",
             f"All fixtures use inherited g={g}, alpha={alpha}, and the same n={n} updates",
             f"(duration {duration}; ideal continuous duration {beta/g:.9g}). No optimization.", "",
             "## Original contact retained throughout the direct pulse", "",
             "| gamma | unconditional Dplus fidelity | outside triple | net controller work |",
             "|---|---:|---:|---:|"]
    for row in direct_rows:
        lines.append(f"| {row['gamma']} | {row['fidelity_unconditional']:.9f} | {row['outside_triple_probability']:.6g} | {row['net_work_on_system']:.8g} |")
    lines += ["", "The pulse is not energy invariant under the old full H. Controller work",
              "equals the full system energy change. Fidelity is preserved after switch-off,",
              "but continued control removes the prepared state; autonomous isolation is absent.",
              "The unchanged G-R11 reader succeeds with this fidelity, and its unread-success",
              "57/82 versus 25/82 backreaction remains. No failed branch is renormalized away.", "",
              "## Explicit 21-state relay completion", "",
              "Twelve off-resonant auxiliaries give -B†B/Delta = g*C - sqrt(2)*g*(I+D).",
              "The density-dependent Stark term is leading order, not a small correction.",
              "Compensation +sqrt(2)*g*D is separately declared; it is not silently dropped.", "",
              "| compensated | Delta | effective F | full unconditional F | auxiliary probability | auxiliary energy | net controller work |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for r in mediator_rows:
        lines.append(f"| {r['compensated']} | {r['delta']:g} | {r['effective_reference_fidelity']:.8f} | {r['fidelity_unconditional']:.8f} | {r['auxiliary_probability']:.6g} | {r['auxiliary_energy']:.6g} | {r['net_work_on_full_model']:.6g} |")
    lines += ["", "Small auxiliary population does not make the energy ledger negligible:",
              "its energy is (4+Delta)*P_aux. The displayed positive full-model switching",
              "work cannot be replaced by the direct effective pulse's negative work.",
              "The sudden-switch protocol counts all high-energy exits, not just a normalized",
              "low-energy conditional state. No measured hardware or autonomous drive is claimed.", "",
              "## Bounded decision", "",
              "The local-operator existence and declared relay audit are complete. The old",
              "2/3 obstruction is specific to a contact family, not joint translation alone.",
              "A physical controller/mediator choice is still required for an implementation.",
              "If pursuing it, test one energy-accounted smooth switching protocol or a passing",
              "trigger; otherwise close this as a controlled preparation model. Do not add",
              "another parameter fit, claim a time arrow, or promote any Route-F gate.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
