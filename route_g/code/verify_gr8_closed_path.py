#!/usr/bin/env python3
"""Closed guided protocol: action, pulses, mass closure and joint witnesses.

Synthetic design inputs; no environmental fit or claimed hardware. The finite
mass harmonic diagnostic is not an all-orders relativistic completion.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp

from verify_gr7_clock_path import joint_state, negativity, path_trace

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr8_closed_path"
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.])
plus = np.ones(2)/np.sqrt(2)


def beam(alpha):
    return (I-1j*(np.cos(alpha)*X+np.sin(alpha)*Y))/np.sqrt(2)


def rotation(pauli, angle):
    return np.cos(angle/2)*I-1j*np.sin(angle/2)*pauli


def run():
    card = json.loads((ROOT/"data/CLOSED_PATH_CARD.json").read_text())
    bbr = json.loads((ROOT/"data/BBR_CLOCK_INPUTS.json").read_text())
    cfg, constants = card["geometry"], bbr["constants"]
    d, ramp, hold, freq, g, temp, mass = [cfg[k] for k in (
        "height_separation_m", "ramp_s", "hold_s", "trap_frequency_Hz",
        "g_m_s2", "temperature_K", "mass_kg")]
    h, c, kb, eps = [constants[k] for k in
                      ("h_J_s", "c_m_s", "kB_J_K", "epsilon0_F_m")]
    hb, omega, total = h/(2*np.pi), 2*np.pi*freq, 2*ramp+hold
    cuts = [0, ramp, ramp+hold, total]
    checks = []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, atol=2e-11, rtol=2e-10):
        check(name, np.allclose(a, b, atol=atol, rtol=rtol))

    def shape(u):
        return 10*u**3-15*u**4+6*u**5

    def path(t):
        if t <= ramp:
            u, sign, offset = t/ramp, 1, 0
        elif t < ramp+hold:
            return d, 0., 0.
        else:
            u, sign, offset = (t-ramp-hold)/ramp, -1, 1
        return (d*(offset+sign*shape(u)),
                sign*d*(30*u**2-60*u**3+30*u**4)/ramp,
                sign*d*(60*u-180*u**2+120*u**3)/ramp**2)

    def integrate(f):
        return sum(quad(f, a, b, epsabs=1e-13, epsrel=1e-12)[0]
                   for a, b in zip(cuts[:-1], cuts[1:]))

    for rel, sha in card["preserved_sha256"].items():
        check("preserved:"+rel, hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == sha)
    check("scope:no_empirical_hardware_or_new_law", not any(card[k] for k in
          ("empirical_data_ready", "hardware_claim", "new_time_law")))
    for t, target in zip(cuts, ((0, 0, 0), (d, 0, 0), (d, 0, 0), (0, 0, 0))):
        close("trajectory:endpoint:"+str(t), path(t), target)
    area = d*(ramp+hold)
    close("area:includes_both_ramps", integrate(lambda t: path(t)[0]), area)
    v2 = integrate(lambda t: (path(t)[1]/2)**2)
    close("kinetic:exact_integral", v2, 5*d*d/(7*ramp))
    close("symmetry:kinetic_difference_zero", integrate(
        lambda t: (path(t)[1]/2)**2-(-path(t)[1]/2)**2), 0.)
    # L at the packet center includes the actual linear guide potential.
    action = []
    for sign in (1, -1):
        act = integrate(lambda t: .5*(sign*path(t)[1]/2)**2
                        +(sign*path(t)[2]/2)*(sign*path(t)[0]/2))
        action.append(mass*act)
    close("ground:action_difference", (action[0]-action[1])/hb, 0.)
    close("ground:each_action", np.array(action)/hb, -.5*mass*v2/hb)
    # Explicit change of scalar potential, with no force change, is physical.
    centered_change = -mass*g*area/hb
    check("guide:force_does_not_determine_phase", abs(centered_change) > 1.)
    close("commensurate:ramp", np.exp(1j*omega*ramp), 1)
    close("commensurate:hold", np.exp(1j*omega*hold), 1)
    for trig in (np.sin, np.cos):
        close("closure:first_order_"+trig.__name__, integrate(
            lambda t: path(t)[2]*trig(omega*(total-t))), 0., atol=2e-12)

    # Independent dimensionless ODE checks both trajectories and their action.
    eta_stress = .002
    tau, w = total/ramp, omega*ramp

    def dimensionless(u):
        q, v, a = path(u*ramp)
        return q/d, v*ramp/d, a*ramp*ramp/d

    def rhs(u, state):
        q, _, acc = dimensionless(u)
        z, vel, zc, vc, integral, action_diff = state
        # Independent L_A-L_B integral, in units m*d^2/ramp.
        lagdiff = ((1+eta_stress)*vc*vel-w*w*zc*(z-q)
                   +acc*zc-eta_stress*g*ramp*ramp/d*z)
        return [vel, (w*w*(q-z)+acc)/(1+eta_stress), vc,
                (-w*w*zc-eta_stress*g*ramp*ramp/d)/(1+eta_stress), z, lagdiff]

    state = np.zeros(6)
    for a, b in zip(np.array(cuts[:-1])/ramp, np.array(cuts[1:])/ramp):
        sol = solve_ivp(rhs, (a, b), state, method="DOP853", rtol=2e-12, atol=2e-13)
        check("ODE:success:"+str(a), sol.success)
        state = sol.y[:, -1]
    Om = omega/np.sqrt(1+eta_stress)
    ic = integrate(lambda t: path(t)[2]*np.cos(Om*(total-t)))
    iss = integrate(lambda t: path(t)[2]*np.sin(Om*(total-t)))
    zrel = -eta_stress/(1+eta_stress)*iss/Om
    vrel = -eta_stress/(1+eta_stress)*ic
    zc = -eta_stress*g/omega**2*(1-np.cos(Om*total))
    vc = -eta_stress*g*Om/omega**2*np.sin(Om*total)
    effective_area = area+eta_stress/omega**2*ic
    close("ODE:finite_mass_separation", state[0], zrel/d, atol=2e-11)
    close("ODE:finite_mass_velocity", state[1], vrel*ramp/d, atol=2e-10)
    close("ODE:common_center", state[2:4], [zc/d, vc*ramp/d], atol=2e-12)
    close("ODE:area", state[4], effective_area/(d*ramp), atol=2e-11)
    phase_direct = mass*((1+eta_stress)*(zc*vrel-vc*zrel)-eta_stress*g*effective_area)/hb
    fourier_real = integrate(lambda t: path(t)[2]*np.cos(Om*t))
    phase_exact = -eta_stress*mass*g*area/hb-mass*eta_stress**2*g*fourier_real/(hb*omega**2)
    close("phase:Gaussian_boundary_term_identity", phase_direct, phase_exact)
    phase_ode = mass*d*d/(ramp*hb)*(state[5]-(1+eta_stress)*state[3]*state[0])
    close("phase:independent_action_ODE", phase_ode, phase_exact, atol=1e-7)

    # Matrix pulse sequence fixes sign and all impulse-limit control phases.
    rho0 = np.outer(plus, plus)
    W = (np.eye(4)-np.kron(X, X)+np.kron(Z, Y)+np.kron(Y, Z))/4
    noncomm = np.kron(Z, Y)@np.kron(Z, Z)-np.kron(Z, Z)@np.kron(Z, Y)
    close("readout:noncommuting_joint_control", noncomm, 2j*np.kron(I, X))
    for pauli, pulse in ((X, rotation(Y, -np.pi/2)),
                         (Y, rotation(X, np.pi/2)), (Z, I)):
        close("readout:rotate_then_Z:"+str(len(checks)), pulse.conj().T@Z@pulse, pauli)
    witness_fixtures = []
    for delta in (0., .21, 1.3, np.pi):
        for alpha in (0., .4, 1.7):
            a_s = .17
            clock_prep = rotation(Y, np.pi/2)
            initial = np.kron(np.array([1., 0.]), np.array([1., 0.]))
            D = np.diag([1., np.exp(-1j*delta), 1., 1.])
            psi = np.kron(beam(alpha).conj().T, I)@D@np.kron(beam(a_s), clock_prep)@initial
            joint_port = abs(psi.reshape(2, 2))**2
            pred = np.array([[(1+np.cos(alpha-a_s))/4,
                              (1+np.cos(alpha-a_s-delta))/4],
                             [(1-np.cos(alpha-a_s))/4,
                              (1-np.cos(alpha-a_s-delta))/4]])
            close(f"pulses:energy_resolved:{delta}:{alpha}", joint_port, pred)
        U = np.diag(np.exp(1j*delta*np.diag(np.kron(Z, Z))/4))
        pure = U@np.kron(rho0, rho0)@U.conj().T
        mixed = U@np.kron(rho0, I/2)@U.conj().T
        common_phase = .71
        propagator = np.diag([1., np.exp(-1j*(common_phase+delta/2)),
                              1., np.exp(-1j*(common_phase-delta/2))])
        path_correction = np.diag([1., 1j])@np.diag(np.exp(1j*delta*np.diag(Z)/4))
        clock_correction = np.diag([1., np.exp(1j*common_phase)])
        prepared_path = beam(0)@np.array([1., 0.])
        local = np.kron(path_correction, clock_correction)
        for initial_clock, target, label in ((rho0, pure, "pure"), (I/2, mixed, "mixed")):
            raw = propagator@np.kron(np.outer(prepared_path, prepared_path.conj()), initial_clock)@propagator.conj().T
            close(f"canonical:full_pulse_{label}:{delta}", local@raw@local.conj().T, target)
        close("pure_mixed:path_identity:"+str(delta), path_trace(pure, 2), path_trace(mixed, 2))
        close("pure:negativity:"+str(delta), negativity(pure, 2), abs(np.sin(delta/2))/2)
        close("mixed:separable:"+str(delta), negativity(mixed, 2), 0)
        wp, wm = [float(np.trace(W@r).real) for r in (pure, mixed)]
        close("witness:pure:"+str(delta), wp, -np.sin(delta/2)/2)
        close("witness:mixed:"+str(delta), wm, (1-np.sin(delta/2))/4)
        witness_fixtures.append(dict(delta_rad=delta, pure_W=wp, dephased_W=wm))
    # Positivity on product states follows from W=P_negative^T_C; test identity.
    ket = (np.kron(plus, np.array([1., -1.])/np.sqrt(2))
           -1j*np.kron(np.array([1., -1.])/np.sqrt(2), plus))/np.sqrt(2)
    projector = np.outer(ket, ket.conj())
    Wpt = projector.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    close("witness:partial_transpose_separable_proof", W, Wpt)

    # High precision is only for evaluating an analytically O(eta^2) closure
    # defect without floating-point cancellation; no fitted input or scan.
    mp.mp.dps = 65
    mm, hh, cc, gg, dd, tr, th = map(lambda v: mp.mpf(str(v)),
                                     (mass, h, c, g, d, ramp, hold))
    ww, TT, hbar = 2*mp.pi*mp.mpf(str(freq)), 2*tr+th, hh/(2*mp.pi)
    AA = dd*(tr+th)
    rows = []
    arad = 8*np.pi**5*kb**4/(15*h**3*c**3)
    for item in bbr["transitions"]:
        name, nu0 = item["name"], item["frequency_Hz"]
        b = -item["delta_alpha_SI"]*arad*temp**4/(2*h*eps*nu0)
        nuT = mp.mpf(str(nu0))*(1+mp.mpf(str(b)))
        eta = hh*nuT/(mm*cc**2)
        OO = ww/mp.sqrt(1+eta)
        # Exact outbound/return Fourier factorization, including transport.
        outbound = dd/tr*mp.quad(lambda u: (60*u-180*u*u+120*u**3)*mp.exp(1j*OO*tr*u), [0, .5, 1])
        FF = outbound*(1-mp.exp(1j*OO*(tr+th)))
        convolution = mp.exp(1j*OO*TT)*mp.conj(FF)
        dz = -eta/(1+eta)*mp.im(convolution)/OO
        dv = -eta/(1+eta)*mp.re(convolution)
        dp = mm*(1+eta)*dv
        co, si, mj = mp.cos(OO*TT), mp.sin(OO*TT), mm*(1+eta)
        vx, vp = hbar/(2*mm*ww), hbar*mm*ww/2
        Vxx = co**2*vx+(si/(mj*OO))**2*vp
        Vpp = (mj*OO*si)**2*vx+co**2*vp
        Vxp = co*si*(vp/(mj*OO)-mj*OO*vx)
        loss_exponent = (dp*dp*Vxx+dz*dz*Vpp-2*dz*dp*Vxp)/(2*hbar*hbar)
        delta = 2*mp.pi*nuT*gg*AA/cc**2
        correction = mm*eta**2*gg*mp.re(FF)/(hbar*ww**2)
        total_delta = delta+correction
        r = mp.exp(-loss_exponent)
        kappa = (1+r*mp.exp(-1j*total_delta))/2
        check(name+":finite_mass_closure_not_overclaimed", loss_exponent > 0)
        check(name+":closure_negligible_in_declared_model", loss_exponent < mp.mpf("1e-24"))
        check(name+":phase_mass_correction_bounded", abs(correction/delta) < mp.mpf("1e-15"))
        close(name+":frequency_phase_area", float(delta), 2*np.pi*float(nuT)*g*area/c**2)
        ideal_loss = 2*mp.sin(delta/4)**2
        ideal_N = abs(mp.sin(delta/2))/2
        # Conservative iid shot-noise illustration for the explicit witness,
        # ideal three-setting variances; no flux or apparatus sensitivity claim.
        n3 = mp.ceil(9*mp.cos(delta/2)**2/(8*ideal_N**2))
        rows.append({
            "clock": name, "thermal_fraction": b, "eta": float(eta),
            "hold_dnu_Hz": float(nuT)*g*d/c**2,
            "total_clock_delta_rad": float(delta),
            "port_phase_excited_minus_ground_rad": -float(total_delta),
            "internal_fringe_phase_rad": float(mp.arg(kappa)),
            "ideal_visibility_loss": float(ideal_loss),
            "ideal_negativity": float(ideal_N),
            "ideal_witness_pure": -float(ideal_N),
            "ideal_witness_dephased": float((1-mp.sin(delta/2))/4),
            "finite_mass_residual_position_m": float(dz),
            "finite_mass_residual_velocity_m_s": float(dv),
            "motional_overlap_loss": float(-mp.expm1(-loss_exponent)),
            "finite_mass_phase_correction_rad": float(correction),
            "witness_three_sigma_shots_per_setting_ideal": int(n3)
        })
    summary = {"checks": len(checks), "passed": sum(v["passed"] for v in checks),
               "failed": [v["name"] for v in checks if not v["passed"]]}
    results = {"stage": "G-R8", "status": card["status"], "summary": summary,
               "protocol": {"total_s": total, "area_m_s": area,
                            "delta_tau_s": g*area/c**2,
                            "ground_phase_rad": 0.0,
                            "each_ground_action_over_hbar": action[0]/hb,
                            "different_centered_guide_ground_phase_rad": centered_change,
                            "max_arm_speed_m_s": 15*d/(16*ramp),
                            "max_arm_acceleration_m_s2": 5*np.sqrt(3)*d/(3*ramp*ramp),
                            "ground_wavepacket_sigma_m": np.sqrt(hb/(2*mass*omega))},
               "clocks": rows, "witness_fixtures": witness_fixtures, "checks": checks,
               "boundary": "Impulse-limit conservative effective model; no guide field calibration, finite pulse, environment, lifetime, count-rate or apparatus claim."}
    OUT.with_suffix(".json").write_text(json.dumps(results, indent=2)+"\n")
    lines = ["# G-R8: specified closed guided protocol", "",
             f"{summary['passed']}/{summary['checks']} checks pass. Synthetic control card, not observations.", "",
             "## One complete sequence and phase convention", "",
             "Prepare the clock at the common height; split two orthogonal guide modes;",
             "5 ms smooth outward motion, 10 ms hold, 5 ms return; apply the inverse",
             "mode coupler and energy-resolved readout at the common height.",
             "The symmetric 1 mm guide separation encloses A=d*(hold+ramp)=1.5e-5 m s.",
             "Both ramps contribute. The guide is defined with a global-z linear",
             "support/acceleration potential and zero branch scalar offsets.",
             "Its ground-state phase difference vanishes by the complete action,",
             "not by deleting the ground rest energy. Total port phases are",
             "Phi_g=alpha_f-alpha_s and Phi_e=Phi_g-delta_i in this model.", "",
             "| Clock | Total delta (rad) | Ideal 1-V | Ideal negativity |", "|---|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['clock']} | {row['total_clock_delta_rad']:.8g} | {row['ideal_visibility_loss']:.4g} | {row['ideal_negativity']:.4g} |")
    lines += ["", "## Closure is checked, not assumed", "",
              "The compensated harmonic guide closes ground wave packets exactly.",
              "The same guide with the excited inertial mass is not exactly closed",
              "to all orders. Integer oscillator periods cancel its first-order",
              "residual; finite-mass residuals and Gaussian overlaps are in the JSON.",
              "An independent ODE at an explicitly artificial mass difference checks",
              "the convolution and phase boundary term. This mass-only diagnostic",
              "does not establish all-orders general relativity or real trap accuracy.", "",
              "## Pure/dephased controls and joint readout", "",
              "After calibrated local phase rotations, the pure state is",
              "exp(i*delta*Zp*Zc/4)|++>. Its energy-dephased control has the same",
              "path visibility but is separable. The witness",
              "W=(I-Xp*Xc+Zp*Yc+Yp*Zc)/4 is a partial transpose of a positive",
              "projector and is nonnegative on every separable state.",
              "For positive small delta its ideal pure value is -sin(delta/2)/2;",
              "the dephased value is (1-sin(delta/2))/4. XX, ZY and YZ use",
              "incompatible local axes but commute globally. Add ZZ as an actually",
              "noncommuting joint control: [ZY,ZZ]=2i IX. Analysis pulses followed",
              "by local population detection implement the specified observables.", "",
              "## Physical priority, not another derivation gate", "",
              "The full guide potential, not just the trajectories or forces, fixes",
              "the common phase. Replacing z by z-q in its linear term leaves forces",
              "unchanged but changes the ground fringe. The declared support field",
              "does not cancel the clock-energy gravitational coupling.",
              "Prioritize calibrated state-resolved phase. Ideal witness shot counts",
              "are already enormous at this geometry, before finite lifetime, controller",
              "records, pulse errors, radiation noise and count-rate limitations.",
              "Do not infer operational feasibility or observed entanglement from the",
              "small mass-closure residual. E2's finite lifetime still matters.",
              "Next select/calibrate one actual guide/coupler or reject this protocol",
              "on a measured noise budget. No new time law or Route-F promotion.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(summary))
    if summary["failed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
