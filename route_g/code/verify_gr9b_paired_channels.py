#!/usr/bin/env python3
"""G-R9b: paired exchange channels of three-node rings, and two ways to one outlet.

Side branch of G-R9, written alongside (not replacing) G-R10. No fit.
numpy only: propagators use the Hermitian eigendecomposition and are cross-checked by an independent Taylor
scaling-and-squaring exponential. n is an assumed update count, not an
observed time. Ring labels a,b are energy-sorted modes of each copy.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/gr9b_paired_channels"
OMEGA = np.exp(2j*np.pi/3)


def shift(n):
    return np.roll(np.eye(n, dtype=complex), 1, axis=0)


def fourier(n):
    return np.exp(2j*np.pi*np.outer(np.arange(n), np.arange(n))/n)/np.sqrt(n)


def spectrum(n, theta):
    """Eigenvalues in Fourier-label order k=0..n-1."""
    return 2-2*np.cos((theta-2*np.pi*np.arange(n))/n)


def laplacian(n, theta):
    s = shift(n)
    return 2*np.eye(n)-np.exp(1j*theta/n)*s-np.exp(-1j*theta/n)*s.conj().T


def sorted_modes(n, theta):
    lam = spectrum(n, theta)
    order = np.argsort(lam)
    return lam[order], fourier(n)[:, order], order


def gaps3(theta):
    lam = np.sort(spectrum(3, theta))
    return lam[1]-lam[0], lam[2]-lam[1]


def bisect(f, lo, hi, iters=200):
    flo = f(lo)
    for _ in range(iters):
        mid = (lo+hi)/2
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo+hi)/2


def propagate(h, s):
    e, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*s*e))@v.conj().T


def taylor_expm(a):
    """Independent exponential: scaling and squaring of a Taylor series."""
    norm = np.linalg.norm(a, 1)
    k = max(0, int(np.ceil(np.log2(max(norm, 1e-300)))) + 4)
    b = a/2**k
    out = np.eye(len(a), dtype=complex)
    term = out.copy()
    for j in range(1, 30):
        term = term@b/j
        out = out+term
    for _ in range(k):
        out = out@out
    return out


def amplitudes(h, psi0, targets, s):
    """<t|exp(-i s H)|psi0> for each target and each s, via eigendecomposition."""
    e, v = np.linalg.eigh(h)
    c0 = v.conj().T@psi0
    phase = np.exp(-1j*np.outer(e, np.atleast_1d(s)))
    return [(t.conj()@v)@(c0[:, None]*phase) for t in targets]


def density_contact(c):
    return np.diag(np.asarray(c, float).ravel()).astype(complex)


def oriented_contact(h, gamma):
    c = np.array([[h[(j-l) % 3] for l in range(3)] for j in range(3)], float)
    return gamma*density_contact(c)


def rank_one_spectrum(eps, weights, gamma):
    """Bright eigenvalues of H0+gamma|v><v| from 1=gamma sum_k W_k/(E-eps_k).

    eps: distinct H0 shell energies (sorted); weights: |v|^2 per shell.
    Returns bright roots and the squared overlap of each root with |v>.
    """
    f = lambda e: 1-gamma*np.sum(weights/(e-eps))
    roots = []
    edges = list(eps)+[eps[-1]+gamma*weights.sum()+1.0]
    for lo, hi in zip(edges[:-1], edges[1:]):
        span = hi-lo
        roots.append(bisect(f, lo+1e-14*max(1, span), hi-1e-14*max(1, span)))
    roots = np.array(roots)
    overlap = np.array([1/(gamma**2*np.sum(weights/(r-eps)**2)) for r in roots])
    return roots, overlap


def run():
    card = json.loads((ROOT/"data/PAIRED_CHANNEL_CARD.json").read_text())
    alpha = card["alpha_per_assumed_update"]
    rng = np.random.default_rng(card["seed"])
    checks = []

    def check(name, ok):
        checks.append({"name": name, "passed": bool(ok)})

    def close(name, a, b, atol=3e-11):
        check(name, np.allclose(a, b, atol=atol, rtol=3e-11))

    for rel, sha in card["preserved_sha256"].items():
        check("preserved:"+rel, hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() == sha)

    # ---- 1. Rings, link sign flip and the spectral mirror -------------------
    grid = np.pi*np.linspace(*card["pairing_grid_over_pi"][:2], card["pairing_grid_over_pi"][2])
    for n in (3, 4):
        f = fourier(n)
        close(f"N={n}:Fourier_unitary", f.conj().T@f, np.eye(n))
        d = np.diag(np.exp(2j*np.pi*np.arange(n)/n))
        for theta in grid[::8]:
            tag = f"N={n}:theta={theta:.6f}"
            l = laplacian(n, theta)
            close(tag+":Fourier_eigenvectors", l@f, f*spectrum(n, theta))
            # Flipping the sign of every link adds pi per link, n*pi per loop.
            close(tag+":link_sign_flip", 4*np.eye(n)-l, laplacian(n, theta+n*np.pi))
            close(tag+":flux_period_is_gauge", laplacian(n, theta+2*np.pi), d@l@d.conj().T)
            close(tag+":conjugation", l.conj(), laplacian(n, -theta))
        if n == 4:
            z = np.diag((-1.)**np.arange(4))
            close("N=4:bipartite_self_mirror", z@laplacian(4, 1.1)@z, 4*np.eye(4)-laplacian(4, 1.1))
    for theta in grid:
        close(f"N=3:theta={theta:.6f}:spectral_mirror",
              np.sort(spectrum(3, np.pi-theta)), np.sort(4-spectrum(3, theta)))

    # ---- 2. Channel pairing for any two flux 3-rings ------------------------
    g01 = np.array([gaps3(t)[0] for t in grid])
    g12 = np.array([gaps3(t)[1] for t in grid])
    close("pairing:lower_gap_closed_form", g01, 2*np.sqrt(3)*np.sin((np.pi-grid)/3))
    close("pairing:upper_gap_closed_form", g12, 2*np.sqrt(3)*np.sin(grid/3))
    check("pairing:lower_gap_strictly_decreasing", np.all(np.diff(g01) < 0))
    check("pairing:upper_gap_strictly_increasing", np.all(np.diff(g12) > 0))
    close("pairing:gap_exchange_under_mirror", [gaps3(np.pi-t)[1] for t in grid], g01)
    pairing_rows = []
    for theta_a, lower in zip(grid, g01):
        tag = f"pairing:theta_A={theta_a:.6f}"
        theta_b = bisect(lambda t: gaps3(t)[1]-lower, 1e-12, np.pi-1e-12)
        close(tag+":resonant_partner_is_mirror", theta_b, np.pi-theta_a, atol=1e-11)
        la, lb = np.sort(spectrum(3, theta_a)), np.sort(spectrum(3, theta_b))
        e11, e02, e20 = la[1]+lb[1], la[0]+lb[2], la[2]+lb[0]
        close(tag+":target_resonant", e02, e11, atol=1e-10)
        close(tag+":reverse_channel_degenerate", e20, e11, atol=1e-10)
        pairing_rows.append({"theta_A": float(theta_a), "theta_B": float(theta_b),
                             "E02_minus_E11": float(e02-e11), "E20_minus_E11": float(e20-e11)})

    # ---- 3. Any real density-density contact sees the pair equally ----------
    f3 = fourier(3)
    idx = lambda a, b, nb=3: nb*a+b
    for sample in range(card["random_contact_samples"]):
        tag = f"equal_coupling:sample={sample}"
        c = rng.normal(size=(3, 3))
        ff = np.kron(f3, f3)  # sorted = Fourier order for 3-rings with 0<theta<pi
        vm = ff.conj().T@density_contact(c)@ff
        m12, m23, m31 = vm[idx(0, 2), idx(1, 1)], vm[idx(1, 1), idx(2, 0)], vm[idx(2, 0), idx(0, 2)]
        z = np.sum(c*OMEGA**np.subtract.outer(np.arange(3), np.arange(3)))/9
        close(tag+":target_coupling_formula", m12, z)
        close(tag+":equal_magnitudes", abs(vm[idx(2, 0), idx(1, 1)]), abs(m12))
        close(tag+":resonant_block_circulant", [m23, m31], [z, z])
        close(tag+":equal_diagonal_shift", np.diag(vm)[[idx(0, 2), idx(1, 1), idx(2, 0)]],
              np.full(3, c.sum()/9))
        sym = (c+c.T)/2
        vs = ff.conj().T@density_contact(sym)@ff
        close(tag+":symmetric_contact_real_loop", vs[idx(0, 2), idx(1, 1)].imag, 0)
        # Identical rings + symmetric contact: swap symmetry, hence G-R9's exact 1/2.
        swap9 = np.eye(9).reshape(3, 3, 9).transpose(1, 0, 2).reshape(9, 9)
        l2 = laplacian(3, np.pi/2)
        hs = np.kron(l2, np.eye(3))+np.kron(np.eye(3), l2)+density_contact(sym)
        close(tag+":identical_rings_symmetric_contact_swap", hs@swap9, swap9@hs)

    # ---- 4. Rank-one (vertex) contact: exact dark-state bound ---------------
    def vertex_model(theta_a, nb, theta_b, gamma):
        la, fa, _ = sorted_modes(3, theta_a)
        lb, fb, _ = sorted_modes(nb, theta_b)
        h = (np.kron(laplacian(3, theta_a), np.eye(nb))+np.kron(np.eye(3), laplacian(nb, theta_b))
             + gamma*np.kron(np.diag(np.eye(3)[0]), np.diag(np.eye(nb)[0])))
        ket = lambda a, b: np.kron(fa[:, a], fb[:, b])
        energies = (la[:, None]+lb[None, :]).ravel()
        return h, ket, energies

    def secular_check(tag, h, energies, gamma, nsites, e0, bright_weight):
        shells = []
        for e in np.sort(energies):
            if not shells or abs(e-shells[-1][0]) > 1e-9:
                shells.append([e, 0])
            shells[-1][1] += 1
        eps = np.array([s[0] for s in shells])
        w = np.array([s[1] for s in shells])/nsites  # |v_ab|^2 = 1/(N_A N_B)
        roots, ov = rank_one_spectrum(eps, w, gamma)
        dark = np.concatenate([[s[0]]*(s[1]-1) for s in shells if s[1] > 1] or [[]])
        close(tag+":secular_equation_spectrum", np.sort(np.concatenate([roots, dark])),
              np.linalg.eigvalsh(h), atol=1e-9)
        # |<n|b>|^2 for the normalized bright vector of the E0 shell.
        b_ov = bright_weight*ov/(roots-e0)**2/(1/gamma)**2
        close(tag+":bright_weights_sum_to_one", b_ov.sum(), 1, atol=1e-9)
        return roots, b_ov

    vertex_rows = []
    gammas = card["vertex_gamma_fixtures"]+card["vertex_gamma_strong_fixtures"]
    pair_thetas = [np.pi/2]+[np.pi*x for x in card["mirror_pair_fixtures_theta_A_over_pi"]]
    for theta_a in pair_thetas:
        for gamma in gammas:
            tag = f"vertex33:theta_A={theta_a:.6f}:gamma={gamma}"
            h, ket, energies = vertex_model(theta_a, 3, np.pi-theta_a, gamma)
            la = np.sort(spectrum(3, theta_a))
            e0 = la[1]+np.sort(spectrum(3, np.pi-theta_a))[1]
            i, t, r = ket(1, 1), ket(0, 2), ket(2, 0)
            bright = (t+i+r)/np.sqrt(3)
            for dvec in [(i-t)/np.sqrt(2), (i+t-2*r)/np.sqrt(6)]:
                close(tag+":dark_state_exact", h@dvec, e0*dvec, atol=1e-10)
            roots, b_ov = secular_check(tag, h, energies, gamma, 9, e0, 3/9)
            s_max = max(12*np.pi/gamma, 60.)
            s = np.linspace(0, s_max, 6001)
            a_t, a_r, g_bb = amplitudes(h, i, [t, r, bright], s)
            g_bb = g_bb*np.sqrt(3)  # <b|U|b> from <b|U|i>, since <b|U|d>=0 and <b|i>=1/sqrt3
            g_analytic = (b_ov[:, None]*np.exp(-1j*np.outer(roots, s))).sum(0)
            close(tag+":return_amplitude_from_secular_roots", g_bb, g_analytic, atol=1e-8)
            close(tag+":dark_state_formula", a_t, (g_bb-np.exp(-1j*e0*s))/3, atol=1e-10)
            close(tag+":outlets_identical_amplitudes", a_r, a_t, atol=1e-10)
            p_t, p_r = abs(a_t)**2, abs(a_r)**2
            check(tag+":exact_bound_4_over_9", p_t.max() <= 4/9+1e-12 and p_r.max() <= 4/9+1e-12)
            if gamma < .1:
                check(tag+":bound_nearly_reached", p_t.max() > 4/9-2e-3)
            vertex_rows.append({"theta_A": float(theta_a), "gamma": gamma,
                                "max_P_0_2_sampled": float(p_t.max()), "max_P_2_0_sampled": float(p_r.max())})

    # ---- 5. Exit A: distinguishable spectra (3-ring + 4-ring) ---------------
    theta_a = theta_a_fix = np.pi*card["theta_A_over_pi"]
    la, fa, _ = sorted_modes(3, theta_a)
    target_gap = la[1]-la[0]
    middle_gap4 = lambda t: np.diff(np.sort(spectrum(4, t)))[1]
    theta_b = bisect(lambda t: middle_gap4(t)-target_gap, 1e-12, np.pi-1e-12)
    close("exitA:resonant_flux_closed_form", theta_b, 4*np.arcsin(np.sqrt(3)/4), atol=1e-11)
    lb = np.sort(spectrum(4, theta_b))
    close("exitA:B_spectrum_closed_form", lb,
          [2-np.sqrt(13)/2, 2-np.sqrt(3)/2, 2+np.sqrt(3)/2, 2+np.sqrt(13)/2])
    energies = la[:, None]+lb[None, :]
    e0 = energies[1, 1]
    shell = np.argwhere(abs(energies-e0) < 1e-9)
    check("exitA:initial_shell_is_exactly_initial_and_target",
          sorted(map(tuple, shell.tolist())) == [(0, 2), (1, 1)])
    other = abs(energies-e0)
    other[1, 1] = other[0, 2] = np.inf
    close("exitA:reverse_detuning_closed_form", energies[2, 0]-e0, (3*np.sqrt(3)-np.sqrt(13))/2)
    close("exitA:min_other_detuning_is_reverse", other.min(), abs(energies[2, 0]-e0))
    # Whole resonant domain: theta_A from the 4-ring limit (lower gap 2 sqrt2) to pi.
    theta_c = bisect(lambda t: gaps3(t)[0]-2*np.sqrt(2), 1e-12, np.pi-1e-12)
    domain_rows, reverse_min = [], np.inf
    for ta in np.linspace(theta_c, np.pi, 2001)[1:-1]:
        lo = np.sort(spectrum(3, ta))
        tb = 4*np.arcsin((lo[1]-lo[0])/4)
        hi = np.sort(spectrum(4, tb))
        en = lo[:, None]+hi[None, :]
        cnt = int(np.sum(abs(en-en[1, 1]) < 1e-9))
        reverse_min = min(reverse_min, abs(en[2, 0]-en[1, 1]))
        domain_rows.append(cnt)
    check("exitA:domain_initial_shell_always_two", all(c == 2 for c in domain_rows))
    check("exitA:domain_reverse_detuning_above_0.3", reverse_min > .3)

    exit_a_rows = []
    for gamma in card["vertex_gamma_fixtures"]:
        tag = f"exitA:gamma={gamma}"
        h, ket, en = vertex_model(theta_a, 4, theta_b, gamma)
        i, t, r = ket(1, 1), ket(0, 2), ket(2, 0)
        dark = (i-t)/np.sqrt(2)
        close(tag+":dark_state_exact", h@dark, e0*dark, atol=1e-10)
        roots, b_ov = secular_check(tag, h, en, gamma, 12, e0, 2/12)
        # Bright root near E0: E0 + (gamma/6)(1 + kappa*gamma + O(gamma^2)).
        off = en[abs(en-e0) > 1e-9]
        kappa = np.sum(1/(e0-off))/12
        near = roots[np.argmin(abs(roots-e0-gamma/6))]
        check(tag+":bright_root_second_order", abs((near-e0)/(gamma/6)-1-kappa*gamma) < 2*gamma**2)
        s_weak = 6*np.pi/gamma  # bright shift gamma/6 -> P=sin^2(gamma s/12)
        ns = np.arange(int(.9*s_weak/alpha), int(1.12*s_weak/alpha))
        a_t, a_r, a_b = amplitudes(h, i, [t, r, (i+t)/np.sqrt(2)], alpha*ns)
        g_bb = a_b*np.sqrt(2)
        close(tag+":dark_state_formula", a_t, (g_bb-np.exp(-1j*e0*alpha*ns))/2, atol=1e-10)
        p_t, p_r = abs(a_t)**2, abs(a_r)**2
        k = int(p_t.argmax())
        n_best = int(ns[k])
        # Independent propagator at the selected update count.
        u1 = taylor_expm(-1j*alpha*h)
        un = np.linalg.matrix_power(u1, n_best)
        close(tag+":independent_propagator", abs(t.conj()@un@i)**2, p_t[k], atol=1e-9)
        check(tag+":near_unit_transfer", 1-p_t[k] < .1*gamma**2)
        check(tag+":reverse_channel_suppressed", p_r[k] < .05*gamma**2)
        check(tag+":beats_identical_ring_bound", p_t[k] > 4/9+.5)
        exit_a_rows.append({"gamma": gamma, "weak_peak_update": s_weak/alpha, "best_update_n": n_best,
                            "bright_root_kappa": float(kappa),
                            "P_0_2": float(p_t[k]), "P_2_0": float(p_r[k]),
                            "one_minus_P_over_gamma2": float((1-p_t[k])/gamma**2)})
    check("exitA:deficit_decreases_in_weak_limit",
          all(a["P_0_2"] < b["P_0_2"] for a, b in zip(exit_a_rows[:-1], exit_a_rows[1:])))

    # ---- 6. Exit B: oriented contact on identical rings ---------------------
    ff = np.kron(f3, f3)
    t33 = np.kron(shift(3), shift(3))
    swap = np.eye(9).reshape(3, 3, 9).transpose(1, 0, 2).reshape(9, 9)
    triple = [idx(0, 2), idx(1, 1), idx(2, 0)]
    h_or, h_sym = card["oriented_profile_h"], card["symmetric_profile_h"]
    z_unit = (h_or[0]+h_or[1]*OMEGA+h_or[2]*OMEGA**2)/3
    close("exitB:loop_phase_minus_pi_over_2", 3*np.angle(z_unit), -np.pi/2)
    close("exitB:coupling_magnitude", abs(z_unit), 1/np.sqrt(3))
    # Secular chiral triangle: exp(-i M t*) is a cyclic permutation.
    m = np.zeros((3, 3), complex)
    m[0, 1] = m[1, 2] = m[2, 0] = z_unit
    m = m+m.conj().T
    t_star_unit = 2*np.pi/(3*np.sqrt(3)*abs(z_unit))
    close("exitB:secular_perfect_transfer", abs(propagate(m, t_star_unit)[0, 1]), 1)
    exit_b_rows = []
    for gamma in card["oriented_gamma_fixtures"]:
        for theta_a in [np.pi/2, np.pi/3]:
            tag = f"exitB:gamma={gamma}:theta_A={theta_a:.6f}"
            v = oriented_contact(h_or, gamma)
            h = (np.kron(laplacian(3, theta_a), np.eye(3))+np.kron(np.eye(3), laplacian(3, np.pi-theta_a))+v)
            vm = ff.conj().T@v@ff
            close(tag+":joint_translation_symmetry", v@t33, t33@v)
            sector = [a*3+b for a in range(3) for b in range(3) if (a+b) % 3 == 2]
            check(tag+":momentum_sector_is_resonant_triple", sorted(sector) == sorted(triple))
            rest = [k for k in range(9) if k not in triple]
            leak = vm[np.ix_(rest, triple)]
            close(tag+":triple_closed_exactly", leak, np.zeros_like(leak))
            check(tag+":not_swap_symmetric", np.linalg.norm(h@swap-swap@h) > .1*gamma)
            i, t, r = ff[:, idx(1, 1)], ff[:, idx(0, 2)], ff[:, idx(2, 0)]
            t_star = t_star_unit/gamma
            a_t, a_r = amplitudes(h, i, [t, r], np.array([t_star]))
            close(tag+":unit_transfer_at_t_star", abs(a_t[0])**2, 1, atol=1e-10)
            n_star = int(round(t_star/alpha))
            un = np.linalg.matrix_power(taylor_expm(-1j*alpha*h), n_star)
            p_n = abs(t.conj()@un@i)**2
            close(tag+":independent_propagator", p_n, abs(amplitudes(h, i, [t], np.array([alpha*n_star]))[0][0])**2, atol=1e-9)
            check(tag+":integer_update_near_unit", p_n > 1-1e-5)
            if theta_a == np.pi/2:
                exit_b_rows.append({"gamma": gamma, "t_star": t_star, "update_n": n_star,
                                    "P_0_2_at_n": float(p_n), "P_2_0_at_n": float(abs(r.conj()@un@i)**2)})
    # Reversed orientation sends the quantum to the other outlet.
    v_rev = oriented_contact([h_or[0], h_or[2], h_or[1]], .03)
    h_rev = np.kron(laplacian(3, np.pi/2), np.eye(3))+np.kron(np.eye(3), laplacian(3, np.pi/2))+v_rev
    a_r = amplitudes(h_rev, ff[:, idx(1, 1)], [ff[:, idx(2, 0)]], np.array([t_star_unit/.03]))[0]
    close("exitB:reversed_orientation_reverse_outlet", abs(a_r[0])**2, 1, atol=1e-10)
    # G-R10's dark clock D+ is orthogonal to T; any translation-covariant contact
    # conserves the T-weight 1/3 of |1,1>, so no such contact prepares D+ beyond 2/3.
    d_plus = (ff[:, idx(0, 2)]-2*ff[:, idx(1, 1)]+ff[:, idx(2, 0)])/np.sqrt(6)
    t_vec = (ff[:, idx(0, 2)]+ff[:, idx(1, 1)]+ff[:, idx(2, 0)])/np.sqrt(3)
    profiles = [h_or, h_sym, [h_or[0], h_or[2], h_or[1]]]+[list(rng.normal(size=3)) for _ in range(4)]
    for k, prof in enumerate(profiles):
        tag = f"exitB:Dplus_obstruction:profile={k}"
        hp = (np.kron(laplacian(3, np.pi/2), np.eye(3))+np.kron(np.eye(3), laplacian(3, np.pi/2))
              + oriented_contact(prof, .05))
        s = np.linspace(0, 600, 3001)
        a_d, a_tt = amplitudes(hp, ff[:, idx(1, 1)], [d_plus, t_vec], s)
        close(tag+":T_weight_conserved", abs(a_tt)**2, np.full(len(s), 1/3), atol=1e-10)
        check(tag+":Dplus_fidelity_at_most_2_over_3", (abs(a_d)**2).max() <= 2/3+1e-10)
    # Symmetric (unoriented) translation-covariant contact: closed triple, real loop, 4/9.
    v_s = oriented_contact(h_sym, .03)
    hs = np.kron(laplacian(3, np.pi/2), np.eye(3))+np.kron(np.eye(3), laplacian(3, np.pi/2))+v_s
    close("exitB:symmetric_contact_swap_symmetric", hs@swap, swap@hs)
    s = np.linspace(0, 4000, 8001)
    a_t, = amplitudes(hs, ff[:, idx(1, 1)], [ff[:, idx(0, 2)]], s)
    check("exitB:symmetric_contact_bound_4_over_9", (abs(a_t)**2).max() <= 4/9+1e-12)
    check("exitB:symmetric_contact_bound_reached", (abs(a_t)**2).max() > 4/9-1e-3)

    failed = [c["name"] for c in checks if not c["passed"]]
    summary = {"checks": len(checks), "passed": len(checks)-len(failed), "failed": failed}
    result = {"stage": "G-R9b", "status": card["status"], "summary": summary,
              "exitA": {"theta_A": theta_a_fix, "theta_B": theta_b, "B_spectrum": lb.tolist(),
                        "reverse_detuning": float(energies[2, 0]-e0),
                        "domain_theta_A": [theta_c, np.pi], "domain_min_reverse_detuning": float(reverse_min),
                        "fixtures": exit_a_rows},
              "exitB": {"profile_h": h_or, "loop_phase_over_pi": float(3*np.angle(z_unit)/np.pi),
                        "fixtures": exit_b_rows},
              "vertex_rank_one_fixtures": vertex_rows,
              "pairing_rows": pairing_rows,
              "checks": checks, "boundary": card["not_claimed"]}
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# G-R9b paired channels and two ways to one outlet", "",
             f"{summary['passed']}/{summary['checks']} checks pass. All inputs are dimensionless toy assumptions; no fit.", "",
             "## Result", "",
             "G-R9's extra outlet is not an accident of theta=pi/2. Flipping the sign of every",
             "link maps L to 4-L; on a three-node ring this shifts the loop phase by pi, so",
             "spec L(pi-theta) = 4 - spec L(theta). The lower gap of one 3-ring equals the upper",
             "gap of its mirror partner and vice versa. Hence for ANY two flux 3-rings, resonance",
             "of |1,1> with |0,2> forces theta_B = pi - theta_A and makes |2,0> degenerate too.",
             f"Checked on {len(pairing_rows)} theta_A values in (0,pi); partner found by bisection.", "",
             "Every real density-density contact couples |1,1> to the two outlets with equal",
             "magnitude, and the resonant 3x3 block is circulant with loop phase Phi=3 arg z.",
             "Symmetric contacts give Phi=0.", "",
             "## Exact bound for the G-R9 vertex contact", "",
             "The vertex contact is rank one. Two combinations of the degenerate triple are exact",
             "dark eigenstates, so P(0,2) = |G_bb(s) - exp(-i E0 s)|^2 / 9 with |G_bb| <= 1.",
             "Therefore P(0,2) <= 4/9 EXACTLY, for every gamma and every update, for identical",
             "rings and for every mirror pair. This sharpens G-R9's exact 1/2 and promotes its",
             "weak-contact 4/9 envelope to an exact bound. General rank-one form: P <= 4 w_i w_t.",
             "Both outlets have IDENTICAL amplitudes at every update, also for mirror pairs,",
             "which have no swap symmetry: the vertex contact cannot prefer either outlet.", "",
             "| theta_A | gamma | max sampled P(0,2) | max sampled P(2,0) |", "|---:|---:|---:|---:|"]
    for row in vertex_rows:
        lines.append(f"| {row['theta_A']:.6f} | {row['gamma']} | {row['max_P_0_2_sampled']:.9f} | "
                     f"{row['max_P_2_0_sampled']:.9f} |")
    lines += ["", "## Exit A: distinguishable spectra (3-ring + 4-ring)", "",
              "A: 3-ring at theta_A=pi/2 (gaps sqrt3, sqrt3). B: 4-ring with",
              f"theta_B = 4 arcsin(sqrt3/4) = {theta_b:.12f}; its spectrum is",
              "{2 -+ sqrt13/2, 2 -+ sqrt3/2}, so B's middle gap is sqrt3.",
              "Only |1,1> and |0,2> share the initial energy; the reverse outlet is detuned by",
              f"(3 sqrt3 - sqrt13)/2 = {energies[2, 0]-e0:.9f}. Over the whole resonant domain",
              f"theta_A in ({theta_c:.6f}, pi) the reverse detuning stays >= {reverse_min:.4f}.",
              "With the SAME G-R9 vertex contact, P(0,2) = |G_bb - exp(-i E0 s)|^2 / 4 can reach one.", "",
              "| gamma | weak-peak update | best n | P(0,2) | P(2,0) | (1-P)/gamma^2 |", "|---:|---:|---:|---:|---:|---:|"]
    for row in exit_a_rows:
        lines.append(f"| {row['gamma']} | {row['weak_peak_update']:.1f} | {row['best_update_n']} | "
                     f"{row['P_0_2']:.9f} | {row['P_2_0']:.3g} | {row['one_minus_P_over_gamma2']:.4f} |")
    lines += ["", "## Exit B: oriented contact on identical rings", "",
              "Contact gamma sum h((j-l) mod 3) P_j x P_l with h=(2,0,1) keeps both copies",
              "identical but orients the link. It commutes with the joint translation, so the",
              "resonant triple is an exactly closed sector (no leakage at any gamma). Its loop",
              "phase is -pi/2: a chiral triangle, which transfers |1,1> -> |0,2> with unit",
              "probability at t* = 2 pi/(3 gamma). Reversing the orientation selects |2,0>.",
              "The symmetric profile h=(1,0,0) has loop phase 0 and stays at 4/9.", "",
              "| gamma | t* | update n | P(0,2) at n | P(2,0) at n |", "|---:|---:|---:|---:|---:|"]
    for row in exit_b_rows:
        lines.append(f"| {row['gamma']} | {row['t_star']:.4f} | {row['update_n']} | "
                     f"{row['P_0_2_at_n']:.10f} | {row['P_2_0_at_n']:.3g} |")
    lines += ["", "Integer-update values differ from one only through the sampling offset of t*.", "",
              "## What changed physically", "",
              "One-for-one exchange needs the two copies to be distinguishable: by their spectra",
              "(Exit A) or by an orientation of the contact (Exit B). Identical 3-rings, or a",
              "mirror pair, joined by an unoriented local contact split the quantum between two",
              "outlets: exactly for identical copies (swap, <=1/2) and for the vertex contact",
              "(<=4/9), and in the weak-contact block otherwise.",
              "Both exits still contain one supplied choice (a resonant flux; an oriented",
              "contact ratio). Neither is derived from a dynamical principle.", "",
              "## Next bounded decision", "",
              "G-R10 built the all-channel readout for the vertex contact and found a conserved",
              "2/3 dark weight of |1,1>. Lemma 3 is the same structure: |1,1> has bright weight 1/3.",
              "Exit B does not lift G-R10's preparation obstruction: every translation-covariant",
              "contact conserves the |T> weight 1/3 of |1,1>, and D+ is orthogonal to |T>, so",
              "the fidelity to D+ never exceeds 2/3 (checked for oriented, symmetric and random",
              "profiles). Exits A and B are declared alternatives for one-outlet EXCHANGE only.",
              "G-R10's own next choice (a declared preparation/record instrument, or the",
              "unfiltered stationary mixture) is unchanged by this side branch.",
              "No physical seconds, records, metric or particle spectrum are claimed.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(summary, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    run()
