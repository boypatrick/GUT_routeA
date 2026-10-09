#!/usr/bin/env python3
"""G-R9: fixed spectral network and ordered-update candidate, not emergent time.

No fit. Independent site-basis matrix evolution checks the exact closed
four-level sector. n is an assumed update count, not an observed time.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr9_event_network"
I3 = np.eye(3, dtype=complex)
S = np.roll(I3, 1, axis=0)
F = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)


def spectrum(theta):
    return 2-2*np.cos(theta/3-2*np.pi*np.arange(3)/3)


def laplacian(theta):
    return 2*I3-np.exp(1j*theta/3)*S-np.exp(-1j*theta/3)*S.conj().T


def sector(theta, g):
    """Fixed logical |00>,|01>,|10>,|11> = |0,1>,|0,2>,|1,1>,|1,2>."""
    lam = spectrum(theta)
    da, db = np.diff(lam)
    h = np.diag([0., db, da, da+db]).astype(complex)
    h[1, 2] = h[2, 1] = g
    return h, da, db


def site_generator(theta, g):
    l = laplacian(theta)
    left = np.kron(F[:, 1], F[:, 1])
    right = np.kron(F[:, 0], F[:, 2])
    v = g*(np.outer(left, right.conj())+np.outer(right, left.conj()))
    return np.kron(l, I3)+np.kron(I3, l)+v


def coherence(psi):
    a = psi.reshape(2, 2)
    ra, rb = a@a.conj().T, a.T@a.conj()
    return 2*ra[1, 0], 2*rb[1, 0]


def run():
    card = json.loads((ROOT/"data/EVENT_NETWORK_CARD.json").read_text())
    alpha, g = card["alpha_per_assumed_update"], card["exchange_g"]
    checks, rows = [], []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, atol=3e-11):
        check(name, np.allclose(a, b, atol=atol, rtol=3e-11))

    for rel, sha in card["preserved_sha256"].items():
        check("preserved:"+rel, hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == sha)
    close("Fourier_unitary", F.conj().T@F, I3)
    close("zero_flux_spectrum", np.linalg.eigvalsh(laplacian(0)), [0, 3, 3])
    close("quarter_flux_spectrum", np.linalg.eigvalsh(laplacian(np.pi/2)),
          [2-np.sqrt(3), 2, 2+np.sqrt(3)])
    rng = np.random.default_rng(20261009)
    embedding = np.column_stack([np.kron(F[:, a], F[:, b]) for a, b in
                                 [(0, 1), (0, 2), (1, 1), (1, 2)]])
    plus = np.ones(4, dtype=complex)/2
    exchange_in = np.array([0, 0, 1, 0], dtype=complex)
    # Fixed fixtures / property checks; not a parameter fit or physical scan.
    for theta in np.pi*np.linspace(*card["theta_domain_over_pi"], 9):
        tag = f"theta={theta:.6f}"
        l, lam = laplacian(theta), spectrum(theta)
        b = I3-np.exp(1j*theta/3)*S
        close(tag+":positive_edge_factorization", b.conj().T@b, l)
        close(tag+":analytic_spectrum", np.linalg.eigvalsh(l), lam)
        close(tag+":fixed_eigenvectors", l@F, F*lam)
        close(tag+":trace", np.trace(l), 6)
        close(tag+":local_diagonals", np.diag(l), [2, 2, 2])
        check(tag+":ordered_simple_modes", np.all(np.diff(lam) > 0))
        phases = np.exp(1j*rng.normal(size=3))
        d = np.diag(phases)
        transformed = d@l@d.conj().T
        close(tag+":gauge_spectrum", np.linalg.eigvalsh(transformed), lam)
        # The ordered product of oriented hopping amplitudes is invariant.
        flux = (-l[1, 0])*(-l[2, 1])*(-l[0, 2])
        close(tag+":loop_phase", flux, np.exp(1j*theta))
        close(tag+":loop_gauge", (-transformed[1, 0])*(-transformed[2, 1])
              *(-transformed[0, 2]), flux)
        # Passive permutation is no more a spectrum change than rephasing.
        close(tag+":permutation", np.linalg.eigvalsh(S@l@S.conj().T), lam)
        u = expm(-1j*alpha*l)
        close(tag+":unitary_update", u.conj().T@u, I3)
        close(tag+":persistent_modes", u@F, F*np.exp(-1j*alpha*lam))
        # Spectral mode stability under a genuinely changed Hermitian rule.
        e = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
        e = (e+e.conj().T)/2
        eps = .05*(lam[1]-lam[0])
        e *= eps/np.linalg.norm(e, 2)
        vals, vecs = np.linalg.eigh(l+e)
        check(tag+":Weyl_gap", vals[1]-vals[0] >= lam[1]-lam[0]-2*eps-1e-12)
        off = np.linalg.norm((I3-np.outer(F[:, 0], F[:, 0].conj()))@vecs[:, 0])
        check(tag+":mode_angle_bound", off <= eps/(lam[1]-lam[0]-eps)+1e-12)
        h, da, db = sector(theta, g)
        full = site_generator(theta, g)
        origin = lam[0]+lam[1]
        close(tag+":exact_sector_not_RWA", full@embedding, embedding@(h+origin*np.eye(4)))
        close(tag+":coupled_spectrum", np.linalg.eigvalsh(h),
              [0, (da+db)/2-np.hypot((da-db)/2, g),
               (da+db)/2+np.hypot((da-db)/2, g), da+db])
        check(tag+":positive_sector_in_domain", da*db > g*g)
        close(tag+":selected_total_excitation_conserved",
              h@np.diag([0, 1, 1, 2])-np.diag([0, 1, 1, 2])@h, np.zeros((4, 4)))
        fullu = expm(-1j*alpha*full)
        for n in [0, 1, 7, 31, 63]:
            s = alpha*n  # Interpolation coordinate, not independently measured time.
            predicted = expm(-1j*s*h)@exchange_in
            actual = np.linalg.matrix_power(fullu, n)@embedding@exchange_in
            close(tag+f":full_vs_reduced:n={n}", actual,
                  np.exp(-1j*s*origin)*embedding@predicted)
            omega = np.hypot((da-db)/2, g)
            p = g*g/omega**2*np.sin(omega*s)**2
            close(tag+f":exchange_probability:n={n}", abs(predicted[1])**2, p)
            free = expm(-1j*s*sector(theta, 0)[0])@plus
            close(tag+f":free_clock_readout:n={n}", coherence(free),
                  [np.exp(-1j*da*s), np.exp(-1j*db*s)])

    for frac in card["fixtures_theta_over_pi"]:
        theta = np.pi*frac
        h, da, db = sector(theta, g)
        detuning = da-db
        omega = np.hypot(detuning/2, g)
        n_peak = round(np.pi/(2*alpha*omega))
        pmax = g*g/omega**2  # Continuous envelope; integer updates need not attain it.
        sampled = abs((np.linalg.matrix_power(expm(-1j*alpha*h), n_peak)@exchange_in)[1])**2
        check(f"fixture:{frac}:sample_below_envelope", sampled <= pmax+1e-12)
        close(f"fixture:{frac}:detuning_identity", detuning, -6*np.cos((theta-2*np.pi)/3))
        rows.append({"theta_over_pi": frac, "eigenvalues": spectrum(theta).tolist(),
                     "gap_A": da, "gap_B": db, "free_ratio_B_over_A": db/da,
                     "detuning": detuning, "continuous_exchange_envelope": pmax,
                     "near_first_peak_update_n": n_peak, "sampled_exchange_probability": sampled})

    # Actual interacting clock signals, not a free-gap ratio asserted for all g.
    h, da, db = sector(np.pi/2, g)
    for n in [0, 1, 11, 31, 63, 89]:
        s = alpha*n
        expected = np.exp(-1j*np.sqrt(3)*s)*np.cos(g*s)
        close(f"resonance:interacting_two_dials:{n}", coherence(expm(-1j*s*h)@plus),
              [expected, expected])
    close("resonance:visibility_can_vanish", coherence(expm(-1j*np.pi/(2*g)*h)@plus), [0, 0])
    # This zero occurs in the interpolation; no claim that an integer sample hits it.
    s = 8.
    detuned, da, db = sector(np.pi/3, g)
    interacting = np.array(coherence(expm(-1j*s*detuned)@plus))
    free = np.exp(-1j*s*np.array([da, db]))
    check("detuned:free_clock_formula_is_not_exact", np.linalg.norm(interacting-free) > .01)
    # Same graph does not choose a generator: f(L)=L+eta L^2 is another fixed rule.
    lam = spectrum(np.pi/2)
    eta = .2
    f = lam+eta*lam**2
    alternate_ratio = np.diff(f)[1]/np.diff(f)[0]
    check("rule_nonuniqueness:same_eigenvectors_different_ratio", abs(alternate_ratio-1) > .1)
    close("rule_nonuniqueness:functional_calculus", F.conj().T@(laplacian(np.pi/2)
          +eta*laplacian(np.pi/2)@laplacian(np.pi/2))@F, np.diag(f))
    # Probe changes weights, never the underlying eigenvalues.
    close("site_probe:equal_weights", abs(F[0, :])**2, np.ones(3)/3)
    close("mode_probe:one_visible_mode", abs(F.conj().T@F[:, 0])**2, [1, 0, 0])
    check("opposite_fractional_response_not_universal_time",
          rows[0]["gap_A"] > rows[1]["gap_A"] > rows[2]["gap_A"]
          and rows[0]["gap_B"] < rows[1]["gap_B"] < rows[2]["gap_B"])
    # Fixed vertex contact exposes the missing degenerate channel. This is
    # an alternative interaction, not quietly the same selected-mode model.
    contact_rows = []
    h0 = np.diag((lam[:, None]+lam[None, :]).ravel())
    site0 = np.diag([1., 0., 0.])
    ff = np.kron(F, F)
    contact_matrix = ff.conj().T@np.kron(site0, site0)@ff
    close("contact:fixed_vertex_all_channels", contact_matrix, np.ones((9, 9))/9)
    swap = np.eye(9).reshape(3, 3, 9).transpose(1, 0, 2).reshape(9, 9)
    resonant = np.eye(9)[:, [2, 4, 6]]
    close("contact:three_degenerate_states", resonant.conj().T@h0@resonant, 4*np.eye(3))
    close("contact:rank_one_secular_block", resonant.conj().T@contact_matrix@resonant,
          np.ones((3, 3))/9)
    target_sector = np.eye(9)[:, [1, 2, 4, 5]]
    omitted = np.eye(9)-target_sector@target_sector.T
    close("contact:nonclosed_two_clock_sector_norm",
          np.linalg.norm(omitted@contact_matrix@target_sector, 2), np.sqrt(20)/9)
    for gamma in card["local_contact_gamma_fixtures"]:
        tag = f"contact:gamma={gamma}"
        total = h0+gamma*contact_matrix
        close(tag+":exchange_symmetry", total@swap, swap@total)
        v = gamma/9
        n = round(np.pi/(3*v*alpha))
        state = expm(-1j*n*alpha*total)@np.eye(9)[:, 4]
        p = abs(state)**2
        secular = 4/9*np.sin(3*v*n*alpha/2)**2
        close(tag+":equal_outputs_exact", p[2], p[6])
        check(tag+":single_target_exact_half_bound", max(p[2], p[6]) <= .5+1e-12)
        check(tag+":weak_contact_secular_prediction", abs(p[2]-secular) < gamma**2)
        # Secular expression is separately checked against its exact 3x3 block.
        small = expm(-1j*n*alpha*v*np.ones((3, 3)))@np.array([0, 1, 0])
        close(tag+":secular_analytic", abs(small)**2, [secular, 1-2*secular, secular])
        leakage = max(0., float(1-p[[2, 4, 6]].sum()))
        contact_rows.append({"gamma": gamma, "update_n": n,
                             "P_0_2": float(p[2]), "P_2_0": float(p[6]),
                             "P_1_1": float(p[4]), "off_resonant_leakage": leakage,
                             "secular_single_target_probability": float(secular)})
    check("contact:fixture_error_decreases_in_weak_limit",
          all(abs(a["P_0_2"]-a["secular_single_target_probability"]) >
              abs(b["P_0_2"]-b["secular_single_target_probability"])
              for a, b in zip(contact_rows[:-1], contact_rows[1:])))
    failed = [v["name"] for v in checks if not v["passed"]]
    summary = {"checks": len(checks), "passed": len(checks)-len(failed), "failed": failed}
    result = {"stage": "G-R9", "status": card["status"], "summary": summary,
              "fixtures": rows, "local_contact_fixtures": contact_rows,
              "alternative_rule_ratio_at_pi_over_2": alternate_ratio,
              "detuned_readout_deviation_at_s8": float(np.linalg.norm(interacting-free)),
              "checks": checks, "boundary": card["not_claimed"]}
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# G-R9 fixed network and two event-dial candidate", "",
             f"{summary['passed']}/{summary['checks']} checks pass. All inputs are dimensionless toy assumptions.", "",
             "## Result", "",
             "One fixed graph rule produces persistent spectral rays, two independent-copy",
             "transition dials and an exact exchange sector. A loop phase, not spatial",
             "rotation, changes the spectrum at unchanged edge magnitudes and trace.",
             "Ordered updates and mode-selective interaction are supplied: physical time,",
             "saved event records and that interaction's microscopic origin are NOT derived.", "",
             "| theta/pi | gap A | gap B | free B/A | exchange envelope | integer n | sampled P |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['theta_over_pi']:.6g} | {r['gap_A']:.8g} | {r['gap_B']:.8g} | "
                     f"{r['free_ratio_B_over_A']:.8g} | {r['continuous_exchange_envelope']:.8g} | "
                     f"{r['near_first_peak_update_n']} | {r['sampled_exchange_probability']:.8g} |")
    lines += ["", "The envelope is a bound over a continuous interpolation, not a guaranteed",
              "integer-update maximum or a fitted experimental probability. g=0.1 and alpha=0.25",
              "are declared examples. Full 9x9 site-basis evolution checks the invariant 4x4 sector.", "",
              "At theta=pi/2 both free gaps are sqrt(3). In the same interacting model,",
              "the two plus-state dial coherences are exp(-i sqrt(3) s) cos(g s).",
              "Equal dial phases thus coexist with reduced contrast; at zero contrast",
              "a phase clock is unreadable. This is coherent exchange, not dissipative synchronization.", "",
              "## What changed physically", "",
              "Increasing the loop phase across the fixed fixtures slows A but speeds B",
              "per assumed update. Equal trace and local magnitudes are not sufficient",
              "clock data. Trace/3 is NOT defined as Ed. No universal slowing law follows.",
              "One-for-one resonance is isolated at theta=pi/2; no stable mechanism selects",
              "that value. General rational common-multiple exchange would need additional",
              "higher-order channels and is not derived by this model.", "",
              f"An equally admissible alternative generator L+0.2 L^2 gives B/A={alternate_ratio:.8g}",
              "at the same graph theta=pi/2. The graph alone does not determine the update law.", "",
              "## Natural local-contact audit: the ideal two-channel sector fails", "",
              "A fixed contact gamma*P_vertex0 tensor P_vertex0 has every Fourier-basis",
              "matrix entry gamma/9. At resonance |0,2>, |1,1>, |2,0> are ALL degenerate.",
              "The contact does not preserve the selected four-level two-dial sector.",
              "Swap symmetry gives equal probabilities in the two outgoing targets, so",
              "either single target has probability <=1/2 EXACTLY at all updates.",
              "The leading weak-contact resonant block gives each target an envelope 4/9",
              "(combined 8/9), not unit transfer. The full 9x9 calculation retains leakage.", "",
              "| gamma | update n | P(0,2) | P(2,0) | off-resonant leakage |",
              "|---|---:|---:|---:|---:|"]
    for r in contact_rows:
        lines.append(f"| {r['gamma']} | {r['update_n']} | {r['P_0_2']:.9g} | "
                     f"{r['P_2_0']:.9g} | {r['off_resonant_leakage']:.6g} |")
    lines += ["", "These are counterchecks, not a fitted repair. The selected-mode exchange",
              "cannot be identified with this local contact. No claim of uniqueness of",
              "this contact or no-go theorem for all possible microscopic rules is made.", "",
              "## Next bounded decision", "",
              "Prefer the honest three-channel local-contact model for further event/readout",
              "work. Derive its full two-copy observation algebra and ask whether one",
              "relational clock description survives channel leakage and unread outcomes.",
              "Alternatively supply an independent symmetry/physical source that forbids",
              "the extra channel before reusing the ideal selected-mode exchange.",
              "Do not launch another apparatus scan or assert emergent seconds/mass.",
              "G-R1/G-R2 records and stationary constraints are NOT silently combined with",
              "this ordered-update model. G-R4-G-R8 remain controls, not mainline gates.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
