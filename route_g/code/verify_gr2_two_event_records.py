#!/usr/bin/env python3
"""G-R2: two stored event records in a bounded, finite history constraint.

No Ed, fitted clock rate, spacetime metric, apparatus energy completion,
or irreversible memory is inferred. The old noninteracting G-R1
Hamiltonian constraint is NOT asserted to generate these recording gates.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from verify_gr1_relational_clock import clock_ket, projector

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "gr2_two_event_records"
TOL = 2e-12


def kron(*args):
    result = np.array([1.0])
    for arg in args:
        result = np.kron(result, arg)
    return result


def swap_blank(outcome):
    result = np.eye(3, dtype=complex)
    result[[0, outcome + 1]] = result[[outcome + 1, 0]]
    return result


def card():
    eye = np.eye(3, dtype=complex)
    plus = np.ones(3, dtype=complex) / np.sqrt(3)
    p = [projector(eye[:, 0]), eye - projector(eye[:, 0])]
    q = [projector(plus), eye - projector(plus)]
    u = np.diag(np.exp(-2j * np.pi * np.arange(3) / 3))
    wi = sum(kron(p[a], swap_blank(a), eye) for a in range(2))
    wf = sum(kron(q[b], eye, swap_blank(b)) for b in range(2))
    uw = kron(u, eye, eye)
    # This is a second instrument with the SAME initial POVM effects.
    flip = np.diag([1., 1., -1.])
    feedback = sum(kron(flip if a == 2 else eye,
                        projector(eye[:, a]), eye) for a in range(3))
    return dict(eye=eye, plus=plus, p=p, q=q, u=u, uw=uw,
                wi=wi, wf=wf, feedback=feedback, flip=flip)


def history(initial, links):
    slices = [initial]
    for v in links:
        slices.append(v @ slices[-1])
    d, steps = len(initial), len(slices)
    prop = np.zeros((d * steps, d * steps), dtype=complex)
    for q, v in enumerate(links):
        left = slice(q * d, (q + 1) * d)
        right = slice((q + 1) * d, (q + 2) * d)
        prop[left, left] += np.eye(d) / 2
        prop[right, right] += np.eye(d) / 2
        prop[right, left] -= v / 2
        prop[left, right] -= v.conj().T / 2
    constraint = prop.copy()
    constraint[:d, :d] += np.eye(d) - projector(initial)
    state = np.concatenate(slices) / np.sqrt(steps)
    return slices, state, prop, constraint


def pointer_probabilities(state):
    return (abs(state.reshape(3, 3, 3)) ** 2).sum(axis=0)


def reduced_initial_pointer(state):
    tensor = state.reshape(3, 3, 3)
    return np.einsum("saf,sbf->ab", tensor, tensor.conj())


def sequential(rho, initial_kraus, final_kraus, u):
    return np.array([[np.trace(n @ u @ m @ rho @ m.conj().T
                              @ u.conj().T @ n.conj().T).real
                      for n in final_kraus] for m in initial_kraus])


def run():
    checks = []

    def check(name, ok, value=None):
        row = dict(name=name, passed=bool(ok))
        if value is not None:
            row["value"] = float(value)
        checks.append(row)

    def close(name, actual, expected):
        error = float(np.linalg.norm(np.asarray(actual) - np.asarray(expected)))
        check(name, error < TOL, error)

    old_paths = sorted((ROOT / "code").glob("verify_g[12]_*.py"))
    old_paths += [ROOT / "code" / "verify_gr1_relational_clock.py",
                  ROOT / "output" / "gr1_relational_clock.json"]
    old_hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in old_paths}
    m = card()
    e, p, q, u = m["eye"], m["p"], m["q"], m["u"]
    wi, wf, uw = m["wi"], m["wf"], m["uw"]
    blank = e[:, 0]
    initial = kron(m["plus"], blank, blank)
    links = [wi, uw, wf @ uw, uw, uw]
    slices, omega, prop, constraint = history(initial, links)
    for name in ["wi", "wf", "uw", "feedback"]:
        v = m[name]
        close(name + ":unitary", v.conj().T @ v, np.eye(27))
    for name, projectors in [("initial", p), ("final", q)]:
        close(name + ":complete", sum(projectors), e)
        close(name + ":orthogonal", projectors[0] @ projectors[1], np.zeros((3, 3)))
        for j, effect in enumerate(projectors):
            close(f"{name}:{j}:idempotent", effect @ effect, effect)
    close("initial:rank_two_second_effect", np.trace(p[1]), 2)
    close("free:three_tick_recurrence", np.linalg.matrix_power(u, 3), e)
    basis = np.column_stack([clock_ket(3, 2 * np.pi * k / 3) for k in range(3)])
    close("clock:GR1_phase_grid", basis.conj().T @ basis, e)

    for j, state in enumerate(slices):
        close(f"slice:{j}:normalization", np.vdot(state, state), 1)
        close(f"slice:{j}:history_conditioning", omega.reshape(6, 27)[j] * np.sqrt(6), state)
        close(f"slice:{j}:clock_probability", np.linalg.norm(omega.reshape(6, 27)[j]) ** 2, 1/6)
        if j < 5:
            close(f"link:{j}:propagation", links[j] @ state, slices[j+1])
    close("constraint:hermitian", constraint, constraint.conj().T)
    close("constraint:zero_residual", constraint @ omega, np.zeros(162))
    eigenvalues, eigenvectors = np.linalg.eigh(constraint)
    check("constraint:positive", eigenvalues[0] > -TOL, eigenvalues[0])
    check("constraint:unique_kernel", np.count_nonzero(abs(eigenvalues) < TOL) == 1)
    close("constraint:ground_state", abs(np.vdot(eigenvectors[:, 0], omega)), 1)
    # Independent sum-of-squares construction of K_prop.
    direct = np.zeros_like(prop)
    for j, v in enumerate(links):
        a = np.zeros((27, 162), dtype=complex)
        a[:, j*27:(j+1)*27] = -v
        a[:, (j+1)*27:(j+2)*27] = np.eye(27)
        direct += a.conj().T @ a / 2
    close("constraint:sum_of_squares", direct, prop)

    rho = projector(m["plus"])
    delta = u @ u
    joint = sequential(rho, p, q, delta)
    exact = np.array([[1., 2.], [1., 5.]]) / 9
    close("joint:analytic_table", joint, exact)
    close("joint:pointer_matches_Kraus", pointer_probabilities(slices[3])[1:, 1:], joint)
    close("joint:normalization", joint.sum(), 1)
    close("joint:initial_marginal", joint.sum(axis=1), [1/3, 2/3])
    close("joint:conditional_plus", joint[:, 0] / joint.sum(axis=1), [1/3, 1/6])
    no_initial = float(np.trace(q[0] @ delta @ rho @ delta.conj().T).real)
    close("backreaction:no_initial_record", no_initial, 0)
    close("backreaction:recorded_plus", joint[:, 0].sum(), 2/9)
    rho_unread = sum(a @ rho @ a for a in p)
    close("backreaction:unread_initial_instrument",
          np.trace(q[0] @ delta @ rho_unread @ delta.conj().T), 2/9)
    alternative_kraus = [p[0], m["flip"] @ p[1]]
    for a in range(2):
        close(f"alternative:{a}:same_effect",
              alternative_kraus[a].conj().T @ alternative_kraus[a], p[a])
    alternative = sequential(rho, alternative_kraus, q, delta)
    close("alternative:analytic_table", alternative, np.array([[1., 2.], [3., 3.]]) / 9)
    alternate_state = wf @ uw @ uw @ m["feedback"] @ wi @ initial
    close("alternative:explicit_pointer", pointer_probabilities(alternate_state)[1:, 1:], alternative)
    close("alternative:plus", alternative[:, 0].sum(), 4/9)

    # Pointer sectors include a distinct blank state; records have fixed
    # protocol stamps i=(epoch0,k1), f=(epoch1,k0), not arbitrary timestamps.
    filled = np.diag([0., 1., 1.])
    count = kron(e, filled, e) + kron(e, e, filled)
    counts = []
    for j, expected in enumerate([0, 1, 1, 2, 2, 2]):
        state = slices[j]
        measured = float(np.vdot(state, count @ state).real)
        counts.append(measured)
        close(f"count:{j}:mean", measured, expected)
        close(f"count:{j}:sharp", count @ state, expected * state)
        if j >= 1:
            close(f"memory:{j}:initial_reduced_state",
                  reduced_initial_pointer(state), reduced_initial_pointer(slices[1]))
        if j >= 3:
            close(f"memory:{j}:joint_preserved", pointer_probabilities(state),
                  pointer_probabilities(slices[3]))
    for j in range(1, 5):
        for a in range(3):
            pointer_a = kron(e, projector(e[:, a]), e)
            close(f"memory:link{j}:commutes_A{a}",
                  links[j] @ pointer_a - pointer_a @ links[j], np.zeros((27, 27)))
    for j in [3, 4]:
        for b in range(3):
            pointer_f = kron(e, e, projector(e[:, b]))
            close(f"memory:link{j}:commutes_F{b}",
                  links[j] @ pointer_f - pointer_f @ links[j], np.zeros((27, 27)))
    # Nonselective pointer readout changes possible coherences, not statistics.
    final_rho = projector(slices[3])
    dephased = np.zeros_like(final_rho)
    for a in range(3):
        for b in range(3):
            sector = kron(e, projector(e[:, a]), projector(e[:, b]))
            dephased += sector @ final_rho @ sector
            close(f"pointer_read:{a}:{b}:repeatability", sector @ sector, sector)
    close("pointer_read:statistics_unchanged", np.diag(dephased), np.diag(final_rho))
    close("pointer_read:trace", np.trace(dephased), 1)

    # Forget epoch, condition only on the recurrent k=0 clock effect.
    phase_zero_rho = (projector(slices[0]) + projector(slices[3])) / 2
    final_valid = kron(e, e, filled)
    alias_probability = np.trace(phase_zero_rho @ final_valid).real
    close("alias:terminal_record_given_phase_zero", alias_probability, .5)
    close("alias:terminal_record_given_phase_zero_epoch_one",
          np.vdot(slices[3], final_valid @ slices[3]), 1)
    # Full periodic closure is incompatible with a fresh, durable record.
    cycle = wf @ uw @ uw @ wi
    closure = float(np.linalg.norm(cycle @ initial - initial))
    close("cycle:blank_to_filled_orthogonal", np.vdot(initial, cycle @ initial), 0)
    close("cycle:periodic_closure_fails", closure, np.sqrt(2))
    close("cycle:reversal_restores_blank", cycle.conj().T @ cycle @ initial, initial)
    # A bare system energy-conserving measurement is NOT supplied.
    bare_hs = kron(np.diag([0., 1., 2.]), e, e)
    close("energy:initial_gate_commutes_Hs", wi @ bare_hs - bare_hs @ wi, np.zeros((27, 27)))
    commutator = float(np.linalg.norm(wf @ bare_hs - bare_hs @ wf))
    check("energy:final_gate_requires_apparatus_completion", commutator > .1, commutator)

    # Generality checks at fixed dimension, not a fitted parameter scan.
    v = np.array([1., 2j, -.7], dtype=complex)
    v /= np.linalg.norm(v)
    for name, input_rho in [("nonuniform_pure", projector(v)),
                            ("mixed", .37 * projector(v) + .63 * np.diag([.2, .3, .5]))]:
        work_rho = kron(input_rho, projector(blank), projector(blank))
        evolved = cycle @ work_rho @ cycle.conj().T
        expected = sequential(input_rho, p, q, delta)
        close(name + ":normalized", expected.sum(), 1)
        check(name + ":nonnegative", expected.min() > -TOL)
        for a in range(2):
            for b in range(2):
                sector = kron(e, projector(e[:, a+1]), projector(e[:, b+1]))
                close(f"{name}:pointer:{a}:{b}", np.trace(evolved @ sector), expected[a, b])

    for rel, digest in old_hashes.items():
        check("preserved:" + rel, hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest)
    failed = [row["name"] for row in checks if not row["passed"]]
    result = dict(stage="G-R2", status="bounded-done" if not failed else "failed",
                  scope="Finite two-event protocol constraint; not an Ed/rate/energy-complete apparatus",
                  assumptions=["ordinary quantum mechanics and Born rule", "three-reading clock plus one epoch bit",
                               "two blank-or-binary memories with fixed protocol stamps",
                               "declared ordered gates and pure input boundary",
                               "isolated pointers after writing, no noise or reset",
                               "K_hist replaces, rather than preserves, the G-R1 noninteracting constraint"],
                  dimensions=dict(clock=3, epoch=2, system=3, initial_memory=3,
                                  final_memory=3, total=162),
                  event_stamps=dict(initial=dict(q=1, epoch=0, phase_index=1),
                                    final=dict(q=3, epoch=1, phase_index=0)),
                  clock_labels=[0, 1, 2, 0, 1, 2], event_counts=counts,
                  joint_probabilities=joint.tolist(), exact_joint_fractions=[["1/9", "2/9"], ["1/9", "5/9"]],
                  no_initial_record_plus=no_initial,
                  alternate_instrument_probabilities=alternative.tolist(),
                  phase_only_final_valid_probability=float(alias_probability),
                  periodic_closure_residual=closure,
                  constraint_residual=float(np.linalg.norm(constraint @ omega)),
                  constraint_first_nonzero_eigenvalue=float(eigenvalues[1]),
                  final_gate_bare_energy_commutator=commutator,
                  previous_file_hashes=old_hashes, checks=checks,
                  summary=dict(total=len(checks), passed=len(checks)-len(failed), failed=failed))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    report = f"""# G-R2: one initial record and one final readout

Status: {result['status']}. Checks: {len(checks)-len(failed)}/{len(checks)}.
All quantities are dimensionless; no clock slowdown is fitted.

## Exact two-event test

Rows: initial outcome 0,1. Columns: terminal outcome +,-.

| Initial / terminal | + | - | Initial marginal |
|---|---:|---:|---:|
| 0 | 1/9 | 2/9 | 1/3 |
| 1 | 1/9 | 5/9 | 2/3 |

P(+) = 2/9 with the initial record; P(+) = 0 without it.
A different initial instrument with the SAME POVM effects gives P(+) = 4/9.
The initial effects alone cannot determine a sequential experiment.

## What is stored

Initial memory labels include the fixed stamp (epoch0,k1); terminal
labels include (epoch1,k0). These are predeclared scheduled-event stamps,
not an autonomous recorder of an arbitrary unknown clock reading.
Both outcome records are preserved after q=3 while their pointers are isolated.
Clock labels are 0,1,2,0,1,2; stored-event counts are 0,1,1,2,2,2.
Count refers only to the two declared events, not all phase ticks.

Conditioning on phase k=0 alone mixes the blank and completed histories:
P(terminal valid | k=0) = 1/2. Including epoch1 makes it 1.
A completely periodic blank-to-filled work state would require
W_cycle chi0 = chi0, but their overlap is zero and closure residual is sqrt(2).
This is not a no-go theorem for cyclic clocks with resetting or longer memory.

## Constraint and scope

The 162-dimensional positive open-history constraint has a unique zero
mode. Residual: {result['constraint_residual']:.3e}.
Its first nonzero eigenvalue is {result['constraint_first_nonzero_eigenvalue']:.6g}
in arbitrary constraint units, NOT a physical clock energy gap.
K_hist includes ordered recording gates and an initial boundary;
it replaces the noninteracting G-R1 constraint, rather than proving that
the old Hamiltonian automatically records events.

The final gate does not commute with the bare system Hamiltonian
(commutator norm {commutator:.6g}). No energy-conserving apparatus,
microscopic gate duration, thermodynamic irreversibility or physical
memory lifetime is certified. Reversing the protocol erases the records.
The pointer probabilities are classical-readable; complete decoherence
of all pointer density-matrix coherences is not assumed or proved.

## Next, not executed

G-R3: explicitly choose an Ei/Ed object and measure, then one interaction
and a two-clock observable. No density law, relative clock rate, metric,
mass formula, detector scan, or Route-F promotion is included here.

Run: python3 route_g/code/verify_gr2_two_event_records.py.
The JSON stores every check and hashes of unchanged prior code/results.
The TeX/PDF give the complete finite construction and its assumptions.
"""
    OUT.with_suffix(".md").write_text(report)
    print(json.dumps(result["summary"]))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
