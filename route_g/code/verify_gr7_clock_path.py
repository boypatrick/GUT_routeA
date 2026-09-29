#!/usr/bin/env python3
"""Frequency, phase and visibility from one conservative clock/path action.

Full joint density matrices check the analytic overlap and distinguish
entanglement from mixed-state visibility loss. No noise fit or hardware model.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr7_clock_path"


def joint_state(rho_clock, ua, ub):
    d = len(rho_clock)
    control = np.zeros((2*d, 2*d), dtype=complex)
    control[:d, :d], control[d:, d:] = ua, ub
    initial = np.kron(np.ones((2, 2))/2, rho_clock)
    return control @ initial @ control.conj().T


def path_trace(rho, d):
    return np.trace(rho.reshape(2, d, 2, d), axis1=1, axis2=3)


def negativity(rho, d):
    pt = rho.reshape(2, d, 2, d).transpose(2, 1, 0, 3).reshape(2*d, 2*d)
    return float(-np.minimum(np.linalg.eigvalsh(pt), 0).sum())


def dephase_clock(rho, d):
    out = np.zeros_like(rho)
    for j in range(d):
        proj = np.kron(np.eye(2), np.diag(np.arange(d) == j))
        out += proj @ rho @ proj
    return out


def run():
    card = json.loads((ROOT/"data/CLOCK_PATH_CARD.json").read_text())
    cfg = json.loads((ROOT/"data/BBR_CLOCK_INPUTS.json").read_text())
    cst = cfg["constants"]
    h, kb, c, eps = [cst[k] for k in
                     ("h_J_s", "kB_J_K", "c_m_s", "epsilon0_F_m")]
    nu = np.array([v["frequency_Hz"] for v in cfg["transitions"]])
    alpha = np.array([v["delta_alpha_SI"] for v in cfg["transitions"]])
    beta = -alpha/(2*h*eps*nu)
    arad = 8*np.pi**5*kb**4/(15*h**3*c**3)
    ex = card["illustration"]
    hold, cold, hot = ex["hold_s"], ex["cold_K"], ex["hot_K"]
    gamma = ex["g_scale_m_s2"]*ex["height_separation_m"]/c**2
    b0, b1 = beta*arad*cold**4, beta*arad*hot**4
    weights = np.array([beta[1], -beta[0]])/(beta[1]-beta[0])
    checks = []

    def check(name, ok):
        checks.append(dict(name=name, passed=bool(ok)))

    def close(name, a, b, scale=1):
        check(name, np.allclose(np.asarray(a)/scale, np.asarray(b)/scale,
                               rtol=3e-11, atol=3e-12))

    for rel, sha in card["preserved_sha256"].items():
        check("preserved:"+rel,
              hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == sha)
    check("scope:no_measurement_or_new_law",
          not card["empirical_data_ready"] and not card["new_time_law"]
          and not card["hardware_claim"])
    plus = np.ones(2)/np.sqrt(2)
    pure = np.outer(plus, plus)
    X = np.array([[0, 1], [1, 0]])
    Z = np.diag([1, -1])
    H = (X+Z)/np.sqrt(2)
    for delta in (0., .37, np.pi, 2*np.pi):
        ua, ub = np.diag([1., np.exp(-1j*delta)]), np.eye(2)
        rho = joint_state(pure, ua, ub)
        k = np.trace(pure @ ub.conj().T @ ua)
        analytic = np.exp(-1j*delta/2)*np.cos(delta/2)
        close(f"overlap:{delta}:analytic", k, analytic)
        close(f"partial_trace:{delta}", path_trace(rho, 2)[0, 1], k/2)
        close(f"entanglement:{delta}:negativity", negativity(rho, 2),
              abs(np.sin(delta/2))/2)
        close(f"probability:{delta}:trace", np.trace(rho), 1)
        # Direct projective probability checks the fringe phase sign.
        phi = .43
        v = np.array([1., np.exp(1j*phi)])/np.sqrt(2)
        p = np.vdot(v, path_trace(rho, 2) @ v).real
        close(f"fringe:{delta}:projection", p, (1+np.real(np.exp(1j*phi)*k))/2)
    ua = np.diag([1., np.exp(-.73j)])
    ub = np.eye(2)
    k = np.trace(pure @ ua)
    close("derivation:common_final_clock_unitary_preserves_overlap",
          np.trace(pure @ (H@ub).conj().T @ (H@ua)), k)
    close("derivation:controlled_inverse_restores_overlap",
          np.trace(pure @ (ub.conj().T@ub).conj().T @ (ua.conj().T@ua)), 1)
    # General unequal populations, including an energy eigenstate.
    for p in (.0, .2, 1.):
        psi = np.array([np.sqrt(1-p), np.sqrt(p)])
        rho0 = np.outer(psi, psi)
        kk = np.trace(rho0 @ ua)
        close(f"population:{p}:visibility_squared", abs(kk)**2,
              1-4*p*(1-p)*np.sin(.73/2)**2)
    # Shared ground: a three-level system is not a product of two clocks.
    state3 = np.array([1., np.exp(.2j), np.exp(-.3j)])/np.sqrt(3)
    p3 = np.outer(state3, state3.conj())
    u3 = np.diag(np.exp(-1j*np.array([0., .4, 1.2])))
    r3 = joint_state(p3, u3, np.eye(3))
    k3 = np.trace(p3 @ u3)
    close("qutrit:population_overlap", k3, np.sum(np.diag(u3))/3)
    close("qutrit:partial_trace", 2*path_trace(r3, 3)[0, 1], k3)
    close("qutrit:Schmidt_negativity", negativity(r3, 3),
          np.sqrt(1-abs(k3)**2)/2)
    kk_bad_product = abs(np.cos(np.pi/2))**2
    kk_shared_ground = abs((1+2*np.exp(-1j*np.pi))/3)
    check("qutrit:not_product_of_E2_E3_qubits",
          abs(kk_shared_ground-kk_bad_product) > .3)
    # Exact counterexample: same path signal, different joint entanglement.
    # Ub=Z gives the standard two-qubit cluster state from |+>|+>.
    rpure = joint_state(pure, np.eye(2), Z)
    rmixed = joint_state(np.eye(2)/2, np.eye(2), Z)
    close("counterexample:same_path_signal", path_trace(rpure, 2),
          path_trace(rmixed, 2))
    close("counterexample:pure_entangled", negativity(rpure, 2), .5)
    close("counterexample:mixed_separable", negativity(rmixed, 2), 0)
    witness_operator = (np.eye(4)-np.kron(X, Z)-np.kron(Z, X))/2
    wpure = np.trace(witness_operator @ rpure).real
    wmixed = np.trace(witness_operator @ rmixed).real
    close("witness:pure_negative", wpure, -.5)
    close("witness:mixed_not_negative", wmixed, 0)
    # The separable bound follows analytically from Cauchy-Schwarz, not this fixture.
    for a, b in ((plus, plus), (np.array([1., 0.]), plus),
                 (np.array([1., 1j])/np.sqrt(2), plus)):
        prod = np.kron(a, b)
        check("witness:product_fixture:"+str(len(checks)),
              np.vdot(prod, witness_operator@prod).real >= -1e-12)
    dephased = dephase_clock(rpure, 2)
    close("no_signalling:local_clock_dephasing_keeps_path",
          path_trace(dephased, 2), path_trace(rpure, 2))
    close("no_signalling:local_dephasing_removes_entanglement",
          negativity(dephased, 2), 0)
    # Scalar energy-zero redistribution: Ua -> exp(-ia)Ua,
    # Ub -> exp(-ib)Ub requires zeta -> zeta+(a-b), not a physical change.
    a, b, zeta = .2, -.7, .11
    knew = np.trace(pure @ (np.exp(-1j*b)*ub).conj().T @ (np.exp(-1j*a)*ua))
    close("action:scalar_phase_redistribution",
          np.exp(1j*(zeta+a-b))*knew, np.exp(1j*zeta)*k)
    # A large-amplitude fixture tests the derivative/integral identity
    # without subtracting optical phases of order 1e13 radians.
    rates = np.array([.12, -.03, .05])
    dt = np.array([.2, .5, .3])
    close("action:piecewise_phase_integral",
          2*np.pi*(rates@dt), sum(2*np.pi*r*t for r, t in zip(rates, dt)))
    small = .01
    loss = 2*np.sin(small/4)**2
    # Analytic Taylor remainder, not an exact equality at finite phase:
    # 0 <= delta^2/8 - (1-cos(delta/2)) <= delta^4/384.
    check("small_phase:quadratic_with_bounded_remainder",
          0 <= small**2/8-loss <= small**4/384*(1+1e-6))
    # Four hold kernels. Stable expansion avoids subtracting two optical
    # frequencies whose tiny difference would otherwise round to zero.
    physical = []
    for name, g, ba in (
        ("identical", 0., b0),
        ("gravity_only", gamma, b0),
        ("thermal_only", 0., b1),
        ("gravity_and_thermal", gamma, b1),
    ):
        df = nu*(g+(ba-b0)+g*ba)
        delta = 2*np.pi*df*hold
        normalized_phase = delta/(2*np.pi*nu)
        proper_difference = g*hold
        close("phase:frequency_integral:"+name,
              delta/(2*np.pi*hold), df)
        close("phase:thermal_rejected_proper_time:"+name,
              weights@normalized_phase, proper_difference, 1e-18)
        for i in range(2):
            kap = (1+np.exp(-1j*delta[i]))/2
            close(f"physical:{name}:{i}:visibility",
                  abs(kap), abs(np.cos(delta[i]/2)))
        physical.append(dict(
            case=name, differential_frequency_Hz=df.tolist(),
            relative_clock_phase_rad=delta.tolist(),
            ideal_two_level_visibility=np.abs(np.cos(delta/2)).tolist(),
            ideal_two_level_visibility_loss=(2*np.sin(delta/4)**2).tolist(),
            ideal_pure_state_negativity=(np.abs(np.sin(delta/2))/2).tolist(),
            reconstructed_proper_time_difference_s=float(weights@normalized_phase),
            ideal_equal_qutrit_visibility=float(abs((1+np.exp(-1j*delta).sum())/3))))
    close("regression:gravity_differential_frequency",
          physical[1]["differential_frequency_Hz"], nu*gamma, 1)
    close("regression:GR6_thermal_frequency",
          np.array(physical[2]["differential_frequency_Hz"])/nu,
          json.loads((ROOT/"output/gr6_equivalence_baseline.json").read_text())
          ["log_thermal_contrast_E2_E3"], 1e-16)
    first_zero = 1/(2*nu*gamma*(1+b0))
    check("feasibility:first_zero_seconds_not_ms", np.all(first_zero > 6))
    check("scope:real_radiative_visibility_not_assigned",
          "static BBR" in card["required_next_inputs"][-1])
    failed = [v["name"] for v in checks if not v["passed"]]
    result = dict(
        stage="G-R7", status="conservative-kernel-bounded-done; full-interferometer-open",
        measurement_claim=False, physical_total_fringe_phase_available=False,
        physical_open_system_visibility_available=False,
        card_sha256=hashlib.sha256((ROOT/"data/CLOCK_PATH_CARD.json").read_bytes()).hexdigest(),
        hold_s=hold, illustrative_gamma=gamma,
        physical_hold_cases=physical, gravity_first_visibility_zero_hold_s=first_zero.tolist(),
        pure_mixed_counterexample=dict(path_visibility=0, pure_negativity=.5,
                                      mixed_negativity=0, pure_witness=wpure, mixed_witness=wmixed),
        summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed),
        checks=checks)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    rows = []
    for rec in physical[1:]:
        for i, label in enumerate(("E2", "E3")):
            rows.append(f"| {rec['case']} | {label} | "
                        f"{rec['differential_frequency_Hz'][i]:+.5g} | "
                        f"{rec['relative_clock_phase_rad'][i]:+.5g} | "
                        f"{rec['ideal_two_level_visibility_loss'][i]:.4g} |")
    table = "\n".join(rows)
    OUT.with_suffix(".md").write_text(f"""# G-R7: frequency, phase and clock-path correlations

Status: conservative kernel bounded-done; full apparatus and empirical gates OPEN.
{len(checks)-len(failed)}/{len(checks)} checks pass. No fit or new time law.

## One evolution, three linked observables

Same G-R4/G-R6 gaps: theta_ir=2*pi*nu_i*integral[(1+b_i(T_r))*d_tau_r].
The phase derivative gives frequency. The overlap
kappa=Tr(rho_clock U_B^dagger U_A) fixes the internal contribution
to fringe phase and ideal visibility. It is not an independently fitted
cosine. The full port phase also contains scalar motional/ground/guide
and laser contributions, which are not supplied by differential polarizabilities.

For (g+e_i)/sqrt(2), kappa=exp(-i*delta_i/2)*cos(delta_i/2).
The ideal visibility is |cos(delta_i/2)|; the phase is undefined at zeros.
Pure joint-state negativity is |sin(delta_i/2)|/2.
For one shared-ground qutrit use (1+exp(-i*delta_2)+exp(-i*delta_3))/3,
not the product of two independent qubit overlaps.

## Conditional hold calculation, not an apparatus forecast

Hold {hold} s, illustrative 1 m separation and 300/310 K baths.
Other trajectory segments and control phases are not invented.
Numbers below are for the ideal separate two-level preparations.

| Hold case | Clock | Delta frequency (Hz) | Delta clock phase (rad) | 1 - ideal visibility |
|---|---|---:|---:|---:|
{table}

The static-response phase combination
sum_i w_i*delta_i/(2*pi*nu_i) = Delta tau,
with the unchanged G-R5 weights, rejects the same leading thermal
integral even with path-dependent proper-time weighting. This requires
matched trajectories/readout conventions and calibrated relative phases.
It is not a linear combination of visibilities, nor evidence for a
free-standing time law. Internal scalar energy-zero conventions cannot
change the total predicted fringe when the scalar action is transformed too.

## Visibility does not certify entanglement

For an algebra fixture Ub=Z, Ua=I, the pure input |+> gives a
maximally entangled joint state with negativity 1/2. The energy-dephased
input I/2 gives a separable joint state with negativity zero. Both have
the same completely incoherent reduced path state and zero path visibility.
The joint witness (I-X_path*Z_clock-Z_path*X_clock)/2 is -1/2 for the
pure fixture and zero for the mixed fixture; the separable bound is zero.
This fixture is not a realized pi-phase Yb+ experiment.

A common final clock unitary, or any unconditioned trace-preserving local
clock channel, cannot change the reduced path state. Clock decoherence
is therefore not automatically an extra multiplicative path-visibility
factor. A path-controlled inverse can undo the ideal correlation, but
requires specified controls; postselection must retain its success rates.

## Physical limits and next useful action

At 1 m the ideal gravity-only first visibility zeros need about
{first_zero[0]:.3g} s (E2) and {first_zero[1]:.3g} s (E3).
Published E2 excited-state lifetimes are on the tens-of-ms scale.
Do not extrapolate that unperturbed E2 hold to seconds. The lifetime
references in CLOCK_PATH_CARD.json are cautions, not a pooled decay rate
or complete environmental model. E3 longevity alone does not certify
charged-ion path coherence, guide closure, or laser stability.

The conservative real BBR shift does not determine dissipative field
correlations, photon which-path information or spontaneous decay.
Only under additional factorization assumptions may a separate motional/
environment overlap multiply kappa. No physical open-system visibility is
assigned. A conventional freely falling interferometer can also cancel the
putative uniform-field proper-time contrast; derive complete paths and
pulse phases before substituting g*height*time.

Next prioritize an explicit closed-path/guide/readout sequence and
state-resolved differential phase (linear in delta), with coherent vs
dephased controls and joint correlations when claiming entanglement.
Visibility loss is quadratic for small delta and is a consistency
observable, not the sole or first discovery channel.
No quantum-gravity, new Ed law, experiment, or Route-F promotion is claimed.
Derivation: tex/route_g_clock_path.tex; sources: SOURCES.md.
""")
    print(json.dumps(result["summary"]))
    print("gravity hold:", physical[1])
    print("first visibility zeros (s):", first_zero)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
