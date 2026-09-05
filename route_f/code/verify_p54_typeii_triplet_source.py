#!/usr/bin/env python3
"""Same-action scalar triplet source along the bosonic-improved Higgs direction.

K and V''' are tree-action tensors at the historical P2 background.  Only
the neutral Higgs direction comes from the fixed-VEV bosonic CW matrix.
This mixed-order scalar diagnostic does not determine a fermionic Clebsch,
a physical type-II neutrino matrix, or a loop-accurate Weinberg coefficient.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "route_f"
CACHE = ROOT / "tmp/p54_full_doublet_cw"
OUTPUT = RF / "output/p54_typeii_triplet_source.json"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def herm(x):
    return (x + x.conj().T) / 2


def from_cjson(x):
    return np.asarray(x["real"]) + 1j * np.asarray(x["imag"])


def cjson(x):
    return {"real": np.asarray(x).real.tolist(), "imag": np.asarray(x).imag.tolist()}


def run():
    started = time.time()
    p1path = RF / "code/verify_p54_p1_hessian_spectrum.py"
    inputpath = RF / "output/p54_full_doublet_cw.json"
    p1 = module("p54_typeii_p1", p1path)
    data = json.loads(inputpath.read_text())
    parameters = data["tree_parameters"]
    vevs = data["vacuum"]
    x0 = p1.vacuum_vector(vevs["omega"], vevs["sigma"], vevs["vs"])
    potential = p1.potential_factory(parameters)
    hessian_fn = p1.jax.jit(p1.hessian(potential))
    keybase = hashlib.sha256(
        p1path.read_bytes() + json.dumps(parameters, sort_keys=True).encode()
    ).digest()
    CACHE.mkdir(parents=True, exist_ok=True)
    counts = {"cache_hits": 0, "evaluated": 0}

    def hessian(x, label):
        tick = time.time()
        key = hashlib.sha256(keybase + np.asarray(x, dtype="<f8").tobytes()).hexdigest()
        path = CACHE / (key + ".npz")
        if path.exists():
            with np.load(path) as cached:
                result = cached["h"]
            counts["cache_hits"] += 1
        else:
            result = np.asarray(hessian_fn(p1.anp.asarray(x)), dtype=float)
            result = (result + result.T) / 2
            np.savez_compressed(path, h=result)
            counts["evaluated"] += 1
        print(f"{label}: {time.time() - tick:.2f}s", flush=True)
        return result

    h0 = hessian(x0, "tree background")
    reps = p1.sm_representation_matrices()
    eye = np.eye(p1.N_REAL)
    c3 = -sum(r @ r for r in reps["SU3"])
    c2 = -sum(r @ r for r in reps["SU2"])
    ry = reps["Y"][0]
    r3 = reps["SU2"][-1]
    # R=-iY: iR_Y=+1 and iR_3=-1 specify the neutral component.
    filt = c3 @ c3 + (c2 - 2 * eye) @ (c2 - 2 * eye)
    filt = filt + (1j * ry - eye) @ (1j * ry - eye)
    filt = filt + (1j * r3 + eye) @ (1j * r3 + eye)
    values, vectors = np.linalg.eigh(herm(filt))
    neutral = vectors[:, abs(values) < 1e-8]
    if neutral.shape[1] != 2:
        raise ValueError(f"Expected two complex neutral triplet copies, got {neutral.shape}")
    p54 = np.zeros_like(h0)
    p54[p1.SL_PHI, p1.SL_PHI] = np.eye(54)
    p126 = np.zeros_like(h0)
    p126[p1.SL_SIGMA_RE, p1.SL_SIGMA_RE] = np.eye(126)
    p126[p1.SL_SIGMA_IM, p1.SL_SIGMA_IM] = np.eye(126)
    sectors = []
    sector_ranks = []
    for projector in (p54, p126):
        weights, coeff = np.linalg.eigh(herm(neutral.conj().T @ projector @ neutral))
        sector_ranks.append(int(np.sum(weights > 1e-8)))
        if abs(weights[-1] - 1) > 1e-9 or np.linalg.norm(weights[:-1]) > 1e-9:
            raise ValueError("54 and 126 must each supply one neutral triplet copy")
        vector = neutral @ coeff[:, -1]
        vector *= np.exp(-1j * np.angle(vector[np.argmax(abs(vector))]))
        sectors.append(vector)
    bt = np.column_stack(sectors)
    tr = np.column_stack([math.sqrt(2) * bt.real, -math.sqrt(2) * bt.imag])
    # Load the exact exported complex doublet basis, not its magnitudes.
    basis_record = data["basis"]
    if "complex_embedding_328x4" not in basis_record:
        candidates = [v for v in basis_record.values() if isinstance(v, dict)
                      and "real" in v and np.asarray(v["real"]).shape == (328, 4)]
        if len(candidates) != 1:
            raise ValueError("Cannot identify the unique exported 328x4 complex doublet basis")
        bd = from_cjson(candidates[0])
    else:
        bd = from_cjson(basis_record["complex_embedding_328x4"])
    light = from_cjson(data["retuned_bosonic_eigenpair"]["light_coefficients"])
    q = math.sqrt(2) * (bd @ light).real
    k = (tr.T @ h0 @ tr)
    k = (k + k.T) / 2
    kc = herm(bt.conj().T @ h0 @ bt)
    if np.linalg.eigvalsh(k).min() <= 0:
        raise ValueError("The complete tree triplet block is not positive")

    def source(step):
        hp = hessian(x0 + step * q, f"Higgs +{step}")
        hm = hessian(x0 - step * q, f"Higgs -{step}")
        full = ((hp - hm) / (2 * step)) @ q
        return full, tr.T @ full

    full_j, j = source(.02)
    full_j_check, j_check = source(.01)
    response = -.5 * np.linalg.solve(k, j)
    full_response = tr @ response
    source54 = j.copy()
    source54[[1, 3]] = 0
    source126 = j - source54
    response54 = -.5 * np.linalg.solve(k, source54)
    response126 = -.5 * np.linalg.solve(k, source126)
    response_no_mixing = np.zeros(4)
    for indices in ([0, 2], [1, 3]):
        response_no_mixing[indices] = -.5 * np.linalg.solve(k[np.ix_(indices, indices)], j[indices])
    step_error = float(np.linalg.norm(j - j_check) / max(np.linalg.norm(j), 1e-14))
    full_step_error = float(np.linalg.norm(full_j - full_j_check) / max(np.linalg.norm(full_j), 1e-14))
    canonical_error = float(np.linalg.norm(tr.T @ tr - np.eye(4)))
    real_k = np.block([[kc.real, -kc.imag], [kc.imag, kc.real]])
    leakage = float(np.linalg.norm(h0 @ tr - tr @ k) / max(np.linalg.norm(h0 @ tr), 1e-14))
    equation_error = float(np.linalg.norm(k @ response + j / 2))
    projected_response_error = float(np.linalg.norm(tr.T @ (h0 @ full_response + full_j / 2)))
    electric_error = float(np.linalg.norm((ry + r3) @ q))
    checks = [
        {"name": "two complex neutral Y=1 T3=-1 triplet copies", "pass": neutral.shape[1] == 2},
        {"name": "one triplet copy each from 54 and 126", "pass": sector_ranks == [1, 1]},
        {"name": "complete canonical four-real neutral triplet plane", "pass": canonical_error < 1e-10},
        {"name": "corrected complex Higgs gives unit neutral real direction", "pass": abs(q @ q - 1) < 1e-10 and electric_error < 1e-10},
        {"name": "tree Hessian preserves the complete triplet plane", "pass": leakage < 1e-10},
        {"name": "real and complex triplet mass representations agree", "pass": np.linalg.norm(k - real_k) < 1e-10},
        {"name": "complete tree triplet mass matrix is positive", "pass": np.linalg.eigvalsh(k).min() > 1e-8},
        {"name": "cubic source agrees at two finite-difference steps", "pass": step_error < 1e-8 and full_step_error < 1e-8},
        {"name": "induced triplet response solves the projected field equation", "pass": max(equation_error, projected_response_error) < 1e-11},
        {"name": "separate source responses add to full mixed response", "pass": np.linalg.norm(response54 + response126 - response) < 1e-12},
    ]
    for row in checks:
        row["pass"] = bool(row["pass"])
    sources = (p1path, inputpath, Path(__file__).resolve())
    return {
        "schema": "p54-same-action-typeii-triplet-source-v1", "date": "2026-09-05",
        "scope": {
            "light_direction": "fixed-VEV bosonic-improved eigenvector with all exported complex phases",
            "heavy_K_and_cubic_J": "tree action at the historical P2 background",
            "perturbative_order": "mixed-order scalar diagnostic, not a loop-accurate Wilson coefficient",
            "fermionic_Clebsch_computed": False, "typeII_neutrino_matrix_computed": False,
            "global_flavor_fit_performed": False,
        },
        "conventions": {
            "generator": "R=-iY; neutral triplet has Y=+1,T3=-1",
            "real_basis_order": ["54_R", "126_R", "54_I", "126_I"],
            "source_definition": "J_A=V'''[e_A,q,q]; q=sqrt(2) Re(B_doublet c_bosonic)",
            "response_definition": "z_A/v^2=-(K^-1 J)_A/2 for H0=v/sqrt(2)",
            "units": "omega=1: K in omega^2, J in omega, response in 1/omega",
            "basis_phase": "each pure-representation neutral complex vector has its largest component real positive",
        },
        "vacuum": vevs,
        "complex_triplet_basis_328x2": cjson(bt),
        "complex_triplet_K_over_omega2": cjson(kc),
        "real_triplet_K_over_omega2": k.tolist(),
        "triplet_mass_eigenvalues_over_omega2": np.linalg.eigvalsh(kc).tolist(),
        "complex_54_126_mass_mixing_over_omega2": cjson(kc[0, 1]),
        "canonical_cubic_source_over_omega": j.tolist(),
        "canonical_induced_z_over_v2_times_omega": response.tolist(),
        "complex_induced_amplitude_over_v2_times_omega": cjson((response[:2] + 1j * response[2:]) / math.sqrt(2)),
        "sector_norms": {
            "source_54": float(np.linalg.norm(j[[0, 2]])), "source_126": float(np.linalg.norm(j[[1, 3]])),
            "response_54": float(np.linalg.norm(response[[0, 2]])), "response_126": float(np.linalg.norm(response[[1, 3]])),
            "response_126_from_54_source": float(np.linalg.norm(response54[[1, 3]])),
            "response_126_from_126_source": float(np.linalg.norm(response126[[1, 3]])),
            "response_126_ignoring_mass_mixing": float(np.linalg.norm(response_no_mixing[[1, 3]])),
        },
        "numerics": {"source_step_relative_error": step_error, "full_source_step_relative_error": full_step_error,
                     "canonical_basis_error": canonical_error, "Higgs_electric_charge_error": electric_error,
                     "triplet_Hessian_leakage": leakage, "projected_field_equation_error": projected_response_error},
        "checks": checks,
        "summary": {"passed": sum(row["pass"] for row in checks), "total": len(checks), "all_pass": all(row["pass"] for row in checks)},
        "cache": counts, "runtime_seconds": time.time() - started,
        "sources": [{"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources],
    }


def main():
    result = run()
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"P54 type-II scalar source: {result['summary']['passed']}/{result['summary']['total']}")
    print(json.dumps(result["sector_norms"], sort_keys=True))
    if not result["summary"]["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
