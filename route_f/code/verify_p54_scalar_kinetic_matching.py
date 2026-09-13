#!/usr/bin/env python3
"""Actual scalar-cubic hard contribution to the four-doublet kinetic matrix.

Frozen P54 tree propagators and cubic vertices, at the broken background.
This is a one-step hard-minus-soft scalar bubble, NOT a completed two-site
matching, gauge/Nielsen completion, or pole-mass calculation. No new Hessian
evaluation is permitted: all vertices use the previously verified cache.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"
CACHE = ROOT / "tmp/p54_full_doublet_cw"
OUT = RF / "output/p54_scalar_kinetic_matching"
LOOP2 = 32 * math.pi**2


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def decode(x):
    return np.asarray(x["real"]) + 1j * np.asarray(x["imag"])


def cj(x):
    x = np.asarray(x)
    return {"real": x.real.tolist(), "imag": x.imag.tolist()}


def bubble_slope(x, y):
    """Integral t(1-t)/[(1-t)x+ty]; reject the purely massless kernel."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    if np.any(x < 0) or np.any(y < 0) or np.any(x+y == 0):
        raise ValueError("nonnegative squared masses with one nonzero mass required")
    s, delta = x+y, (x-y)/(x+y)
    out = np.zeros_like(s)
    near = abs(delta) < .2
    # Symmetric expansion also covers exactly degenerate masses.
    for n in range(12):
        out[near] += 2/s[near] * delta[near]**(2*n) / (2*(2*n+1)*(2*n+3))
    edge = ((x == 0) | (y == 0)) & ~near
    out[edge] = 1/(2*s[edge])
    far = ~(near | edge)
    xx, yy = x[far], y[far]
    out[far] = (xx**2-yy**2-2*xx*yy*np.log(xx/yy))/(2*(xx-yy)**3)
    return out


def gram(jets, weight):
    return np.einsum("aij,bij,ij->ab", jets, jets, weight, optimize=True)/LOOP2


def complexify(k):
    n = len(k)//2
    return (k[:n, :n]+k[n:, n:])/2 + .5j*(k[n:, :n]-k[:n, n:])


def realify(k):
    return np.block([[k.real, -k.imag], [k.imag, k.real]])


def run():
    p1path = RF / "code/verify_p54_p1_hessian_spectrum.py"
    fullpath = RF / "output/p54_full_doublet_cw.json"
    phasepath = RF / "output/p54_common_yukawa_phase.json"
    localpath = RF / "output/p54_self_consistent_light.json"
    interfacepath = RF / "code/verify_p54_finite_yukawa_interface.py"
    p1 = module("p54_kinetic_p1", p1path)
    interface = module("p54_kinetic_interface", interfacepath)
    full, phase, local = [json.loads(p.read_text()) for p in (fullpath, phasepath, localpath)]
    pars, v = full["tree_parameters"], full["vacuum"]
    x0 = p1.vacuum_vector(v["omega"], v["sigma"], v["vs"])
    keybase = hashlib.sha256(p1path.read_bytes()+json.dumps(pars, sort_keys=True).encode()).digest()
    cache_keys = []

    def hess(x):
        key = hashlib.sha256(keybase+np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        path = CACHE / (key+".npz")
        if not path.exists():
            raise RuntimeError("Expected existing Hessian cache; refusing new evaluation: "+key)
        cache_keys.append(key)
        with np.load(path) as data:
            return data["h"]

    h0 = hess(x0)
    b = decode(full["basis"]["vectors_in_328_real_complexification"])
    qr, qi = math.sqrt(2)*b.real, -math.sqrt(2)*b.imag
    reps = p1.sm_representation_matrices()
    ry = reps["Y"][0]
    rphase, rcharge = expm(-math.pi*ry), expm(math.pi*reps["SU2"][0])
    tr = np.array([(hess(x0+.02*qr[:,a])-hess(x0-.02*qr[:,a]))/.04 for a in range(4)])
    ti = np.array([rphase @ t @ rphase.T for t in tr])
    independent_i = (hess(x0+.01*qi[:,0])-hess(x0-.01*qi[:,0]))/.02
    independent_r = (hess(x0+.01*qr[:,0])-hess(x0-.01*qr[:,0]))/.02
    neutral = np.concatenate([tr, ti])
    charged = np.array([rcharge @ t @ rcharge.T for t in neutral])
    jets = np.concatenate([neutral, charged])
    qneutral = np.column_stack([qr, qi])
    qfull = np.column_stack([qneutral, rcharge @ qneutral])

    lam, u = np.linalg.eigh(h0)
    hard = lam > 1e-8
    masses2 = np.where(hard, lam, 0.)
    xx, yy = np.broadcast_arrays(masses2[:, None], masses2[None, :])
    include = hard[:, None] | hard[None, :]
    weight = np.zeros_like(xx)
    weight[include] = bubble_slope(xx[include], yy[include])
    jet_eig = np.array([u.T @ t @ u for t in jets])
    k16 = gram(jet_eig, weight)
    k8 = k16[:8, :8]
    k4 = complexify(k8)
    hh = hard[:, None] & hard[None, :]
    mixed = include & ~hh
    k4_hh = complexify(gram(jet_eig[:8], weight*hh))
    k4_hl = complexify(gram(jet_eig[:8], weight*mixed))
    c = decode(full["retuned_bosonic_eigenpair"]["light_coefficients"])
    delta_z = float(np.vdot(c, k4 @ c).real)
    checks = []

    def check(name, residual, tolerance=1e-9):
        checks.append({"name": name, "residual": float(residual),
                       "tolerance": tolerance, "pass": bool(residual < tolerance)})

    check("exactly 290 positive and 38 soft scalar directions", abs(int(hard.sum())-290))
    check("all soft masses vanish within the tree tolerance", float(max(abs(lam[~hard]))), 1e-8)
    check("sixteen canonical real components of four weak doublets", np.linalg.norm(qfull.T@qfull-np.eye(16)))
    check("hypercharge transports real to imaginary cubic jets", np.linalg.norm(independent_i-ti[0]), 1e-7)
    check("independent cubic derivative steps agree", np.linalg.norm(independent_r-tr[0]), 1e-7)
    check("charged rotation fixes the SM-invariant vacuum", np.linalg.norm(rcharge@x0-x0))
    check("complex kinetic matrix is Hermitian", np.linalg.norm(k4-k4.conj().T))
    check("eight-real kinetic matrix has the full complex block structure", np.linalg.norm(k8-realify(k4)))
    check("charged and neutral kinetic tensors coincide", np.linalg.norm(k16[8:,8:]-k8))
    check("no neutral-charged kinetic mixing", np.linalg.norm(k16[:8,8:]))
    check("hard-heavy and mixed hard-soft pieces exhaust the result", np.linalg.norm(k4-k4_hh-k4_hl))
    check("scalar-cubic contribution is a positive weighted Gram matrix", max(0., -np.linalg.eigvalsh(k16).min()))
    ward = max(np.linalg.norm((qfull.T@r@qfull).T@k16+k16@(qfull.T@r@qfull))
               for r in reps["SU2"]+[ry])
    check("all four electroweak Ward identities hold on the actual real doublet plane", ward)
    quadratures = []
    for x, y in ((.3,.3),(.3,.3000000001),(.3,.45),(0.,.3),(1e-5,8.),(7.,1e-7)):
        exact = float(bubble_slope(x,y))
        numeric = quad(lambda t:t*(1-t)/((1-t)*x+t*y),0,1,epsabs=1e-12,epsrel=1e-11)[0]
        quadratures.append({"x":x,"y":y,"analytic":exact,"quadrature":numeric})
    check("independent Feynman-parameter integration including degenerate and massless limits",
          max(abs(r["analytic"]-r["quadrature"])/max(1,abs(r["analytic"])) for r in quadratures), 1e-10)
    check("single real-heavy toy has g^2/(192 pi^2 M^2)",
          abs(float(bubble_slope(.7,.7))/LOOP2 - 1/(192*math.pi**2*.7)))
    try:
        bubble_slope(0,0)
        rejected = False
    except ValueError:
        rejected = True
    check("purely soft IR bubble is rejected by the hard-kernel API", 0 if rejected else 1)

    # Independent nonzero-momentum integral, with no subtraction of two UV
    # divergent numbers. Ordered i,j count and the real-scalar 1/2 factor
    # follow directly from 1/2 Tr log, not from an empirical normalisation.
    qlight = np.r_[c.real, c.imag]
    tlight = np.einsum("a,aij->ij", qlight, jet_eig[:8])
    wa = (tlight*tlight)[include]
    x, y = xx[include], yy[include]
    nodes, weights = np.polynomial.legendre.leggauss(96)
    nodes, weights = (nodes+1)/2, weights/2
    momentum_rows = []
    for p2 in (1e-5, 5e-6):
        finite = sum(w*np.dot(wa, np.log1p(t*(1-t)*p2/((1-t)*x+t*y)))
                     for t,w in zip(nodes,weights))/LOOP2
        momentum_rows.append({"p_E2_over_omega2":p2,"Pi_difference_over_p2":float(finite/p2),
                              "relative_to_analytic":float(abs(finite/p2-delta_z)/delta_z)})
    check("nonzero Euclidean momentum difference tends to the computed slope",
          momentum_rows[-1]["relative_to_analytic"], 2e-4)
    check("halving momentum squared improves the derivative approximation",
          0 if momentum_rows[-1]["relative_to_analytic"]<momentum_rows[0]["relative_to_analytic"] else 1)

    # Use the previously tested complex canonical-matching interface. These
    # are actual scalar legs on old synthetic h/f inputs, not a full matching.
    local0 = local["cases"][0]
    h, f = decode(local0["h_raw"]), decode(local0["f_raw"])
    z = np.eye(4)+k4
    zhalf = interface.hermitian_power(z,.5)
    ccan = zhalf@c/math.sqrt(1+delta_z)
    d = decode(full["retuned_bosonic_eigenpair"]["matrix_over_omega2"])
    yrows = {}
    for species, coeff in phase["raw_copy_coefficients"].items():
        yy = decode(coeff["h_raw"])[:,None,None]*h+decode(coeff["f_raw"])[:,None,None]*f
        conjugate = species in ("d", "e", "down", "charged_lepton")
        # The phase dictionary's down/e coefficients contract with c*, not c.
        conjugate = species in ("Yd", "Ye") or conjugate
        scalar_k = k4.conj() if conjugate else k4
        cc = ccan.conj() if conjugate else ccan
        oldc = c.conj() if conjugate else c
        transformed = interface.finite_match(yy,np.zeros((3,3)),np.zeros((3,3)),scalar_k)
        before = np.einsum("a,aij->ij",oldc,yy)
        after = np.einsum("a,aij->ij",cc,transformed["vertices"])
        err = np.linalg.norm(after-before/math.sqrt(1+delta_z))
        check(species+": scalar leg and light-vector normalisation cancel copy rotations",err)
        yrows[species] = {"before":cj(before),"scalar_leg_only_after":cj(after),"identity_residual":float(err)}
    dc = interface.hermitian_power(z,-.5)@d@interface.hermitian_power(z,-.5)
    check("canonical congruence preserves the tuned zero mode",np.linalg.norm(dc@ccan))
    check("canonically normalised light vector has unit norm",abs(np.linalg.norm(ccan)-1))

    passed = sum(row["pass"] for row in checks)
    report = {"date":"2026-09-06","scope":"actual one-step scalar-cubic hard kinetic subset at frozen broken background",
              "convention":"Euclidean inverse propagator p^2(I+K)+D; K positive for real scalar cubic bubbles",
              "source_sha256":{str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (p1path,fullpath,phasepath,localpath,interfacepath,Path(__file__).resolve())},
              "cache":{"hits":len(cache_keys),"evaluated":0,"keys":cache_keys},
              "scalar_masses":{"hard_count":int(hard.sum()),"soft_count":int((~hard).sum()),
                               "hard_min_squared":float(lam[hard].min()),"hard_max_squared":float(lam.max())},
              "K4_complex":cj(k4),"K4_hard_hard":cj(k4_hh),"K4_hard_soft":cj(k4_hl),
              "K4_eigenvalues":np.linalg.eigvalsh(k4).tolist(),"K16_real":k16.tolist(),
              "bosonic_light_delta_Z":delta_z,"bosonic_light_scalar_leg_factor":1/math.sqrt(1+delta_z),
              "canonical_light_vector":cj(ccan),"canonical_curvature":cj(dc),
              "momentum_regression":momentum_rows,"kernel_quadrature":quadratures,"Yukawa_scalar_leg_only":yrows,
              "not_included":["gauge/vector/ghost and fermion wave functions","finite Yukawa vertex subsets other than scalar leg",
                              "two-site upper-Wilson subtraction","field-dependent higher-derivative operator basis",
                              "momentum-dependent pole and Nielsen identity","a physical Higgs or flavor fit"],
              "checks":checks,"summary":{"passed":passed,"total":len(checks),"all_pass":passed==len(checks)}}
    OUT.with_suffix(".json").write_text(json.dumps(report,indent=2)+"\n")
    lines = ["# P54 actual scalar-cubic kinetic matching subset","","Date: 2026-09-06.","",
             "The calculation reuses the frozen tree action, 290 massive scalar propagators and all 38 soft directions. It includes mixed hard-soft loops and subtracts purely soft loops. It is a one-step broken-background scalar contribution, not the completed staged lower threshold. Massless Goldstone propagators inherit background-field Landau xi=0; this scalar subset is not gauge independent. With exactly zero soft masses, the mixed bubble's soft-region momentum expansion is scaleless in dimensional regularization. Finite soft masses or a staged EFT require a new subtraction, not this shortcut.","",
             "## Derivation","",
             r"For $V\supset\tfrac12\eta_i(M^2_{ij}+T_{aij}q_a)\eta_j$, expand $\Gamma_1=\tfrac12\mathrm{Tr}\log(-D^2+M^2+T_aq_a+\cdots)$. The bubble is $-\tfrac14\mathrm{Tr}(GTqGTq)$; the inverse two-point function therefore has the following finite Euclidean momentum difference:","",
             r"$$\Pi_{ab}(p_E^2)-\Pi_{ab}(0)=\frac1{32\pi^2}\sum_{ij}^{\rm hard}T_{aij}T_{bji}\int_0^1\log\!\left[1+\frac{t(1-t)p_E^2}{(1-t)m_i^2+tm_j^2}\right]dt.$$","",
             r"Thus $K_{ab}=\partial_{p_E^2}\Pi_{ab}(0)$ is a positive weighted Gram matrix. The ordered pair sum contains at least one hard index. Quartic scalar tadpoles have no momentum derivative. A one-loop mass counterterm inside this already one-loop bubble would be a selected two-loop contribution and is not inserted.","",
             r"$$I(x,y)=\int_0^1\frac{t(1-t)dt}{(1-t)x+ty}=\frac{x^2-y^2-2xy\log(x/y)}{2(x-y)^3},\quad I(x,x)=\frac1{6x},\quad I(x,0)=\frac1{2x}.$$","",
             r"The stable near-degenerate series is $I=\frac1m\sum_{n\ge0}d^{2n}/[2(2n+1)(2n+3)]$, with $m=(x+y)/2$ and $d=(x-y)/(x+y)$. A single real-heavy field with interaction $ghS^2/2$ gives $K_h=g^2/(192\pi^2M^2)$.","",
             "## Actual tensors and covariant completion","",
             r"Cubic jets are exact central differences of the quartic action's Hessian. Hypercharge transports the real jets to imaginary jets; an independent imaginary cache check verifies this. SU(2) generates the charged partners, and all 16 real external directions are assembled. The resulting full matrix obeys the four electroweak generator Ward identities.","",
             r"At quadratic order in the doublets, the covariant kinetic operator is $(D_\mu H_a)^\dagger K_{ab}(D^\mu H_b)$. Its three- and four-point gauge couplings are fixed by this operator. This completion does not replace an independent calculation of gauge-loop diagrams or the Nielsen identity, and higher-field derivative operators remain separate.","",
             "```text",np.array2string(k4.real,precision=10),"```","",
             f"The imaginary matrix norm is {np.linalg.norm(k4.imag):.6g}; the four eigenvalues are {np.linalg.eigvalsh(k4).tolist()}.","",
             f"On the actual bosonic-improved light vector, delta Z = {delta_z:.12g}; the scalar-leg-only normalisation factor is {1/math.sqrt(1+delta_z):.12g}.","",
             r"Canonical matching uses $Z=1+K$, $D_c=Z^{-1/2}DZ^{-1/2}$ and $c_c=Z^{1/2}c/\sqrt{c^\dagger Zc}$. It follows that $Y_c(c_c)=Y(c)/\sqrt{c^\dagger Zc}$, with conjugation for down/e. The exact square roots test the algebra; their higher powers are not a two-loop prediction.","",
             "## Validation and scope","",f"{passed}/{len(checks)} checks pass; {len(cache_keys)} old Hessian cache hits, zero new Hessians.","",
             "| Check | Residual | Pass |","|---|---:|:---:|"]
    lines += [f"| {r['name']} | {r['residual']:.3e} | {'yes' if r['pass'] else 'NO'} |" for r in checks]
    lines += ["","The nonzero-momentum Feynman-parameter integral is independently evaluated at two momenta. This closes a scalar momentum-dependence subset, not Higgs pole matching. In particular the heavy-doublet gap used in the original curvature analysis is not reinterpreted as a pole gap.",""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps({"summary":report["summary"],"cache_hits":len(cache_keys),"delta_Z":delta_z,
                      "failed":[r for r in checks if not r["pass"]]},indent=2))
    return report


if __name__ == "__main__":
    result = run()
    raise SystemExit(0 if result["summary"]["all_pass"] else 1)
