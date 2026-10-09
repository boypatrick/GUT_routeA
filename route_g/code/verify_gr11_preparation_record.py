#!/usr/bin/env python3
"""One energy-compatible preparation and one relative-phase record.

An explicit collective apparatus coupling is added; its switching, physical
implementation and pure blank are not derived. All failed preparations count.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from verify_gr10_relational_readout import (
    DPLUS, DMINUS, E9, proj, word, energy_projectors, dephase_energy,
    readout_metrics,
)
from verify_gr9_event_network import spectrum

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr11_preparation_record"


def entropy_bits(p):
    p = np.asarray(p)
    p = p[p > 0]
    return float(-np.sum(p*np.log2(p)))


def make_instrument(m=7):
    d_basis = np.column_stack([DPLUS, DMINUS])
    pd = d_basis@d_basis.conj().T
    effects, kraus = [np.zeros((9, 9), complex)], [E9-pd]
    for d in range(1, m):
        delta = 2*np.pi*d/m
        g = np.array([np.sqrt(2/3)*(np.cos(delta)-1),
                      1j*np.sqrt(2)*np.sin(delta)])
        c = proj(g)/m
        effects.append(d_basis@c@d_basis.conj().T)
        # Rank-one positive square root, not a matrix elementwise square root.
        kraus.append(d_basis@(c/np.sqrt(np.trace(c).real))@d_basis.conj().T)
    memory = np.eye(m+1)
    blank_embedding = np.kron(E9, memory[:, 0, None])
    w = sum(np.kron(k, memory[:, l+1, None]) for l, k in enumerate(kraus))
    a = w@blank_embedding.conj().T-blank_embedding@w.conj().T
    generator = 1j*a
    u = (np.eye(9*(m+1))-blank_embedding@blank_embedding.conj().T
         -w@w.conj().T+a)
    return d_basis, pd, effects, kraus, blank_embedding, w, generator, u


def run():
    card = json.loads((ROOT/"data/PREPARATION_RECORD_CARD.json").read_text())
    readout = json.loads((ROOT/card["parent_card"]).read_text())
    network = json.loads((ROOT/readout["parent_card"]).read_text())
    m = readout["phase_outcomes"]
    dim_m = card["memory_dimension"]
    phases = 2*np.pi*np.arange(m)/m
    db, pd, effects, ks, j, w, g, u = make_instrument(m)
    rho0, rhod = proj(E9[:, 4]), proj(DPLUS)
    checks, gamma_rows = [], []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, tol=4e-11):
        check(name, np.allclose(a, b, atol=tol, rtol=4e-11))

    for rel, sha in card["preserved_sha256"].items():
        check("preserved:"+rel, hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == sha)
    check("one_distinct_blank_plus_seven_records", dim_m == m+1 == 8)
    close("relative_effects_sum_dark_projector", sum(effects), pd)
    check("relative_effects_positive", min(np.linalg.eigvalsh(c).min() for c in effects) > -1e-12)
    close("all_outcomes_instrument_complete", sum(k.conj().T@k for k in ks), E9)
    close("success_effect_square_roots", [k@k for k in ks[1:]], effects[1:])
    close("isometry", w.conj().T@w, E9)
    close("blank_orthogonal_to_every_record", j.conj().T@w, np.zeros((9, 9)))
    close("explicit_unitary", u.conj().T@u, np.eye(9*dim_m))
    close("unitary_implements_isometry", u@j, w)
    close("Hermitian_apparatus_generator", g.conj().T, g)
    close("pulse_exponential", expm(-1j*np.pi*g/2), u)

    probs = np.array([np.trace(k@rho0@k.conj().T).real for k in ks])
    conditional_q = np.array([(2/(3*m))*(1-np.cos(phi))**2 for phi in phases[1:]])
    close("failure_probability_not_discarded", probs[0], 1/3)
    close("success_probability", probs[1:].sum(), 2/3)
    close("record_probabilities_exact", probs, np.r_[1/3, (2/3)*conditional_q])
    close("record_probability_normalization", probs.sum(), 1)
    bright = (E9[:, 2]+E9[:, 4]+E9[:, 6])/np.sqrt(3)
    close("failure_conditional_state", ks[0]@rho0@ks[0]/probs[0], proj(bright))
    success_unread = sum(k@rhod@k.conj().T for k in ks[1:])
    # Exact seven-root sum: T7'(-2)/(1-T7(-2)) = 497/123.
    cosines = np.cos(phases)
    close("Chebyshev_reciprocal_sum", np.sum(1/(2+cosines)), 497/123)
    close("unread_success_state_backaction", success_unread,
          (57/82)*proj(DPLUS)+(25/82)*proj(DMINUS))
    close("unread_success_fidelity", np.vdot(DPLUS, success_unread@DPLUS), 57/82)
    close("unread_success_retains_conditional_orbit",
          readout_metrics(success_unread, phases)["seven_reading_rigidity_residual"], 0)
    for d, k in enumerate(ks[1:], 1):
        post = k@rho0@k.conj().T/probs[d]
        close(f"outcome:{d}:stays_dark", pd@post@pd, post)
        close(f"outcome:{d}:pure_branch", np.trace(post@post), 1)
        close(f"outcome:{d}:conditional_orbit",
              readout_metrics(post, phases)["seven_reading_rigidity_residual"], 0)
    # The apparatus measures a state-dependent observable, not a preloaded
    # classical outcome table for Dplus alone.
    minus_prob = np.array([np.trace(k@proj(DMINUS)@k.conj().T).real for k in ks])
    check("dark_input_dependence_is_nontrivial", np.linalg.norm(minus_prob-probs) > .2)

    final = w@E9[:, 4]
    joint_rho = proj(final)
    pointer = np.einsum("al,am->lm", final.reshape(9, dim_m), final.reshape(9, dim_m).conj())
    close("record_diagonals", np.diag(pointer), np.r_[0., probs])
    check("premeasurement_not_assumed_classical", np.linalg.norm(pointer-np.diag(np.diag(pointer))) > .01)
    # The new record can be erased by continuing the SAME coupling.
    for area in [0., .31, np.pi/2, np.pi]:
        state = expm(-1j*area*g)@j@E9[:, 4]
        close(f"unswitched:blank_probability:{area}", np.linalg.norm(j.conj().T@state)**2,
              np.cos(area)**2)
    close("same_pulse_twice_returns_blank_ray", u@u@j, -j)

    for gamma in network["local_contact_gamma_fixtures"]:
        tag = f"gamma={gamma}"
        lam = spectrum(np.pi/2)
        h = np.diag((lam[:, None]+lam[None, :]).ravel())+gamma*np.ones((9, 9))/9
        htot = np.kron(h, np.eye(dim_m))
        close(tag+":all_Kraus_commute_energy", [h@k-k@h for k in ks], np.zeros((m, 9, 9)))
        close(tag+":energy_intertwiner", htot@w, w@h)
        close(tag+":unitary_energy_conserving", u@htot-htot@u, np.zeros_like(u))
        close(tag+":generator_energy_conserving", g@htot-htot@g, np.zeros_like(g))
        # Compare to the complete G-R10 stationary joint effects before
        # compression. This is an exact coarse-graining, not postselection.
        _, projectors = energy_projectors(h)
        local = [proj(word(phi))/m for phi in phases]
        compressed = []
        same_d = []
        for d in range(m):
            block = [pd@dephase_energy(np.kron(local[r], local[(r+d) % m]), projectors)@pd
                     for r in range(m)]
            compressed.append(sum(block))
            same_d.extend([m*b-effects[d] for b in block])
        close(tag+":all_49_effects_compress_exactly", compressed, effects)
        close(tag+":common_origin_label_carries_no_dark_state_information", same_d, np.zeros((m*m, 9, 9)))
        for pindex, p in enumerate(projectors):
            before = np.trace(rho0@p)
            after = np.vdot(final, np.kron(p, np.eye(dim_m))@final)
            close(tag+f":full_energy_distribution:{pindex}", before, after)
        energy_before = float(np.trace(rho0@h).real)
        energy_after = float(np.vdot(final, htot@final).real)
        failure_energy = float(np.trace(ks[0]@rho0@ks[0]@h).real/probs[0])
        success_energy = float(sum(np.trace(k@rho0@k.conj().T@h).real for k in ks[1:])/(2/3))
        close(tag+":energy_initial", energy_before, 4+gamma/9)
        close(tag+":energy_success", success_energy, 4)
        close(tag+":energy_failure", failure_energy, 4+gamma/3)
        close(tag+":energy_all_branches", energy_after, (2/3)*success_energy+(1/3)*failure_energy)
        for wait in [.73, 5.2]:
            free = expm(-1j*wait*h)
            close(tag+f":waiting_has_no_phase_alignment_requirement:{wait}",
                  [np.trace(k@free@rho0@free.conj().T@k.conj().T) for k in ks], probs)
            future = np.kron(free, np.eye(dim_m))@final
            future_probs = np.sum(abs(future.reshape(9, dim_m))**2, axis=0)
            close(tag+f":record_preserved_after_decoupling:{wait}", future_probs, np.r_[0., probs])
        noncomm = np.linalg.norm(h@np.kron(local[0], np.eye(3))-np.kron(local[0], np.eye(3))@h)
        check(tag+":separate_local_phase_is_not_invariant", noncomm > 1)
        gamma_rows.append({"gamma": gamma, "initial_energy": energy_before,
                           "success_energy": success_energy, "failure_energy": failure_energy,
                           "all_branch_energy": energy_after, "local_phase_commutator_norm": float(noncomm)})

    shannon = entropy_bits(probs)
    close("classical_record_entropy_chain_rule", shannon,
          entropy_bits([1/3, 2/3])+(2/3)*entropy_bits(conditional_q))
    full_pair_probs = np.r_[1/3, np.repeat((2/3)*conditional_q/m, m)]
    saved_bits = entropy_bits(full_pair_probs)-shannon
    close("redundant_common_origin_entropy", saved_bits, (2/3)*np.log2(m))
    failed = [c["name"] for c in checks if not c["passed"]]
    summary = {"checks": len(checks), "passed": len(checks)-len(failed), "failed": failed}
    results = {"stage": card["stage"], "status": card["status"], "summary": summary,
               "record_probabilities_failure_then_d1_to_d6": probs.tolist(),
               "conditional_success_phase_probabilities": conditional_q.tolist(),
               "unread_success_fidelity_exact": "57/82", "unread_success_fidelity": 57/82,
               "classical_record_entropy_bits": shannon, "redundant_origin_bits_saved_per_attempt": saved_bits,
               "conditional_reset_bound_in_kBT": shannon*np.log(2),
               "energy_ledger": gamma_rows, "checks": checks, "boundary": card["not_claimed"]}
    OUT.with_suffix(".json").write_text(json.dumps(results, indent=2)+"\n")
    lines = ["# G-R11 one preparation and one relative-phase record", "",
             f"{summary['passed']}/{summary['checks']} checks pass. No physical clock scale or hardware fit.", "",
             "## A complete abstract instrument, not an autonomous apparatus", "",
             "Keep the G-R10 system contact. Add an explicitly new COLLECTIVE apparatus coupling.",
             "Its 72-dimensional unitary commutes with the full system Hamiltonian and writes",
             "one eight-level degenerate register: a distinct blank, failure, and d=1,...,6.",
             "This is a relative-phase record, NOT separate local A/B phase measurements.",
             "The common phase-origin label is uniformly redundant in the dark sector:",
             "compressing all 49 invariant joint effects retains the full relative instrument.",
             "d=0 has identically zero dark-sector effect. No possible failure or exit is discarded.", "",
             "## Success, failure and backreaction", "",
             "The declared |1,1> input succeeds with probability 2/3 and fails with 1/3.",
             "Success includes the relative-phase readout, so it does not leave pure Dplus",
             "undisturbed. If the relative record is unread, the success state is exactly",
             "(57/82)|Dplus><Dplus| + (25/82)|Dminus><Dminus|.",
             "Every successful branch remains in the stationary dark subspace and retains",
             "the conditional phase orbit, with changed clock state and correlations.", "",
             "| label | unconditional probability |", "|---|---:|"]
    for label, p in zip(["failure"]+[f"d={d}" for d in range(1, m)], probs):
        lines.append(f"| {label} | {p:.9g} |")
    lines += ["", "## Energy and information ledger", "",
              "Bare mean energy is 4+gamma/9 initially, 4 after successful recording, and",
              "4+gamma/3 on failure. Keeping both branches conserves not only the mean but",
              "the full total-energy distribution. The ideal register Hamiltonian is zero.",
              "This does not make blank memory, apparatus synthesis, reset or switching free.", "",
              f"After classicalizing the register its Shannon entropy is {shannon:.9g} bits.",
              f"Dropping only the redundant origin label saves {saved_bits:.9g} bits per attempt.",
              "For a thermal reset with degenerate memory, no accessible side information",
              "and cyclic controls, the ideal lower bound is k_B*T*ln(2) times that entropy.",
              "No temperature or energy unit is supplied. This is not a universal write-energy",
              "cost; coherent registers and accessible side information require separate accounting.", "",
              "## Phase reference and the remaining physical boundary", "",
              "Energy-conserving processing with a stationary ancilla and invariant pointer",
              "can induce only effects commuting with H. The separate A phase effects do not",
              "commute with H; this construction does not implement them without another resource.",
              "The relative effects do commute with H. No external energy-phase reference is",
              "needed for their abstract implementation, but the internal calibration and",
              "collective coupling are specified, not derived from the original vertex contact.",
              "U describes the interaction-picture instrument. The physical propagator also",
              "has free system evolution, which leaves record and energy probabilities unchanged",
              "but need not leave the failed system state unchanged. The pulse parameter is assumed.",
              "A pulse of area pi/2 writes the record. Continuing the SAME interaction to",
              "area pi returns the whole input to blank (up to a phase). Records persist only",
              "after the declared decoupling and under isolated-memory future dynamics.",
              "The premeasurement can retain coherence between outcome labels; irreversibility",
              "and a measured hardware lifetime are not claimed. Memory can be reread, but",
              "remeasuring the system is not required to give the same relative outcome.", "",
              "## Decision", "",
              "The finite energy-compatible instrument-existence task is bounded-done.",
              "Do not add another abstract prerequisite proof. If autonomous time generation",
              "remains the goal, next choose a physical trigger/self-decoupling mechanism and",
              "account for its energy and records. Otherwise keep this as an effective",
              "relational-clock measurement model, with its switching assumption explicit.",
              "G-R9b is an independent retained side branch; no oriented contact is imported.",
              "No Route-F gate or emergent mass/spacetime claim is promoted.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
