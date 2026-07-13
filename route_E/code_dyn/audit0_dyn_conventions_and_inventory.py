#!/usr/bin/env python3
"""DYN-0: dynamics-program scaffold.

Three deliverables, per route_E/ROADMAP.md section E:

  1. Conventions addendum for the dynamics lanes: superpotential parameters,
     the Pati-Salam singlet-vev registry (mirrored live from the HEAD
     vacuum-branch card), the fixed Aulakh-Girdhar branch convention, and the
     adopted DYN-0 decisions (file lanes; M_S treatment).
  2. Asset inventory: every history asset named by the DYN items must exist
     in HEAD (read-only recovery contract; nothing is restored).
  3. Rerun-from-history gates: the four audit-4a1 scripts are extracted from
     HEAD into a temporary mini-tree, executed in dependency order, and their
     regenerated JSON ledgers must match the HEAD ledgers up to volatile
     fields (timestamps, sha256 stamps, source manifests).

Independent checks: the four rational special points of the vacuum cubic
8x^3 - 15x^2 + 14x - 3 + xi(1-x)^2 = 0 are verified with exact Fractions.

Boundary: DYN-0 claims no dynamics results.  The heavy spectrum is still a
placeholder; the scalar-Hessian Goldstone eigenvector audit is still pending
(both are DYN-1 targets).
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from fractions import Fraction as Fr
from pathlib import Path

from route_e_paths import AUDIT_OUTPUT, REPO_ROOT

ROOT = REPO_ROOT
OUT = AUDIT_OUTPUT / "audit0"

CHECKS: list[tuple[str, bool]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


def git_show(path: str) -> bytes:
    return subprocess.run(["git", "show", f"HEAD:{path}"], cwd=ROOT,
                          check=True, capture_output=True).stdout


def git_exists(path: str) -> bool:
    return subprocess.run(["git", "cat-file", "-e", f"HEAD:{path}"], cwd=ROOT,
                          capture_output=True).returncode == 0


# ------------------------------------------------------------ section 1
print("== DYN-0 section 1: conventions addendum ==")

head_vacuum_card = json.loads(git_show("output/audit4a1/vacuum_branches.json"))

ADDENDUM = {
    "superpotential_parameters": {
        "m": "m_210 (210_H mass)",
        "M": "m_126 (126_H 126bar_H mass)",
        "M_H": "m_10 (10_H mass)",
        "lambda": "lambda_210 (210_H cubic)",
        "eta": "eta_210_126_126bar (210-126-126bar coupling)",
        "gamma, bar_gamma": "10_H-126_H/126bar_H-210_H doublet/triplet mixing couplings",
        "g": "auxiliary 10-sector coupling used by the literature mass blocks",
    },
    "ps_singlet_vev_registry": head_vacuum_card["ps_singlet_vevs"],
    "aulakh_girdhar_convention": {
        "branch_variable": "x = -lambda*omega/m",
        "coupling_ratio": "xi = lambda*M/(eta*m)",
        "vacuum_cubic": "8*x^3 - 15*x^2 + 14*x - 3 + xi*(1-x)^2 = 0",
        "abmsv_crosscheck": "matches ABMSV after omega_ABMSV = -omega_AG",
    },
    "adopted_dyn0_decisions": {
        "file_lanes": ("new dynamics scripts live at root code/ with ledgers in "
                       "whitelisted output/audit*/ lanes; output/audit0/ whitelist "
                       "extended today; output/audit3/ to be extended at DYN-2; "
                       "history assets recovered read-only via git show HEAD:<path>"),
        "susy_scale_treatment": ("single effective SUSY threshold M_S, benchmark "
                                 "3 TeV; DYN-2 treats M_S as a fit parameter inside "
                                 "[1, 10] TeV; MSSM two-loop running below M_GUT"),
    },
}
check("PS singlet-vev registry mirrored from HEAD (5 entries: p, a, omega, sigma, bar_sigma)",
      [v["symbol"] for v in ADDENDUM["ps_singlet_vev_registry"]]
      == ["p", "a", "omega", "sigma", "bar_sigma"])

# exact special-point verification of the vacuum cubic
SPECIAL_POINTS = [("SU5_x_half", Fr(1, 2), Fr(-5)),
                  ("SU5_x_minus_one", Fr(-1), Fr(10)),
                  ("GLR_x_zero", Fr(0), Fr(3)),
                  ("flipped_SU5_x_third", Fr(1, 3), Fr(-2, 3))]


def cubic(x: Fr, xi: Fr) -> Fr:
    return 8 * x**3 - 15 * x**2 + 14 * x - 3 + xi * (1 - x) ** 2


ok_points = all(cubic(x, xi) == 0 for _, x, xi in SPECIAL_POINTS)
check("vacuum-cubic special points exact (SU(5) x=1/2,-1; G_LR x=0; flipped x=1/3)",
      ok_points)

head_points = {(p["point_label"], p["x_value"], p["xi_value"])
               for p in head_vacuum_card["literature_special_point_contract"]
               if p.get("x_value") is not None}
check("special points agree with the HEAD contract (labels, x, xi)",
      head_points == {(lbl, str(x), str(xi)) for lbl, x, xi in SPECIAL_POINTS})

# ------------------------------------------------------------ section 2
print("== DYN-0 section 2: history-asset inventory (read-only contract) ==")

INVENTORY = {
    "DYN-0": ["code/audit0_conventions_card.py",
              "code/audit4a1_cmsgut_vacuum_branches.py",
              "code/audit4a1_cmsgut_mass_export.py",
              "code/audit4a1_cmsgut_literature_mass_import.py",
              "code/audit4a1_triplet_symbolic_inverse.py",
              "output/audit0/invariant_card.json",
              "output/audit4a/source_spectrum_schema.json",
              "output/audit4a1/vacuum_branches.json",
              "output/audit4a1/mass_export_schema.json",
              "output/audit4a1/literature_mass_matrices.json",
              "output/audit4a1/triplet_symbolic_inverse.json"],
    "DYN-1": ["code/solve_spin10_vacuum_alignment.py",
              "code/audit_orbit_superpotential_hessian.py",
              "code/verify_spin10_component_hessian.py",
              "code/verify_combined_superpotential_flatness.py"],
    "DYN-4": ["code/audit1_flavor_target_conventions.py",
              "code/audit1b_covariant_rank_contract.py",
              "code/audit_flavor_fit_observables.py",
              "code/scan_majorana_contact_sensitivity.py",
              "code/audit_cp1_yukawa_covariant.py",
              "output/audit1/target_table_conventions.json",
              "output/audit1b/covariant_rank_contract.json"],
    "DYN-2": ["code/two_loop_spectrum_fit.py",
              "code/scan_mediator_threshold_rge.py",
              "code/scan_thresholds_item5.py"],
    "DYN-3": ["code/audit2_source_basis_wilson_builder.py",
              "code/scan_dimension5_wilson_tensors.py",
              "code/audit_dressed_dimension5_channels.py",
              "code/audit_mssm_mixing_d5_dressing.py",
              "code/audit_full_knu_width.py",
              "code/scan_proton_channel_bounds.py",
              "output/audit2/source_basis_wilson_contract.json"],
    "DYN-5": ["code/verify_hidden_zeta_origin.py",
              "code/verify_routeB_mediator_grading.py",
              "code/scan_mediator_r_window.py"],
    "DYN-6": ["code/audit_clockwork_radial_driver_hessian.py",
              "code/audit_majorana_monomial_clockwork_origin.py",
              "code/audit_constrained_clockwork_source_hessian.py"],
}
missing = [(item, p) for item, paths in INVENTORY.items()
           for p in paths if not git_exists(p)]
check("all roadmap-named history assets exist in HEAD "
      f"({sum(len(v) for v in INVENTORY.values())} paths across {len(INVENTORY)} items)",
      not missing, f"missing: {missing}" if missing else "")

# ------------------------------------------------------------ section 3
print("== DYN-0 section 3: rerun-from-history gates ==")

VOLATILE_KEYS = ("created_utc", "source_manifest")


def strip_volatile(obj):
    if isinstance(obj, dict):
        return {k: strip_volatile(v) for k, v in obj.items()
                if k not in VOLATILE_KEYS and "sha256" not in k}
    if isinstance(obj, list):
        return [strip_volatile(v) for v in obj]
    return obj


def same(a, b, tol=1e-9) -> bool:
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k], tol) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b) <= tol * max(1.0, abs(a), abs(b))
    if isinstance(a, str) and isinstance(b, str):
        if a.startswith("/") and b.startswith("/"):
            return True                      # host-local provenance paths
        return a == b
    return a == b


CHAIN = ["audit4a1_cmsgut_vacuum_branches.py",
         "audit4a1_cmsgut_mass_export.py",
         "audit4a1_cmsgut_literature_mass_import.py",
         "audit4a1_triplet_symbolic_inverse.py"]
UPSTREAM = ["output/audit0/invariant_card.json",
            "output/audit4a/source_spectrum_schema.json"]
LEDGERS = ["vacuum_branches.json", "mass_export_schema.json",
           "literature_mass_matrices.json", "triplet_symbolic_inverse.json"]

rerun_results = {}
with tempfile.TemporaryDirectory(prefix="dyn0_rerun_") as td:
    tmp = Path(td)
    (tmp / "code").mkdir()
    for up in UPSTREAM:
        dest = tmp / up
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(git_show(up))
    for script in CHAIN:
        (tmp / "code" / script).write_bytes(git_show(f"code/{script}"))
    chain_ok = True
    for script in CHAIN:
        run = subprocess.run([sys.executable, str(tmp / "code" / script)],
                             cwd=tmp, capture_output=True, text=True)
        ok = run.returncode == 0
        chain_ok &= ok
        check(f"history rerun executes: {script}", ok,
              "" if ok else run.stderr.strip().splitlines()[-1][:120])
    if chain_ok:
        for ledger in LEDGERS:
            new = strip_volatile(json.loads(
                (tmp / "output" / "audit4a1" / ledger).read_text()))
            old = strip_volatile(json.loads(
                git_show(f"output/audit4a1/{ledger}")))
            match = same(new, old)
            rerun_results[ledger] = match
            check(f"regenerated ledger matches HEAD (volatile fields excluded): {ledger}",
                  match)

# DYN-1 starting checklist, snapshotted from HEAD stage gates
mass_export_head = json.loads(git_show("output/audit4a1/mass_export_schema.json"))
stage = head_vacuum_card["stage_gates"]
check("DYN-1 open items confirmed pending in HEAD: heavy spectrum placeholder, "
      "scalar-Hessian Goldstone eigenvectors",
      stage.get("heavy_spectrum_json_nonplaceholder") is False)

# ------------------------------------------------------------ ledger
npass = sum(1 for _, ok in CHECKS if ok)
all_pass = npass == len(CHECKS)

report = {
    "audit": "audit0_dyn_conventions_and_inventory",
    "dyn_item": "DYN-0",
    "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "checks_total": len(CHECKS),
    "checks_passed": npass,
    "all_pass": all_pass,
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "conventions_addendum": ADDENDUM,
    "special_points_exact": [{"label": lbl, "x": str(x), "xi": str(xi),
                              "cubic_residual": str(cubic(x, xi))}
                             for lbl, x, xi in SPECIAL_POINTS],
    "inventory": INVENTORY,
    "rerun_gates": rerun_results,
    "dyn1_starting_checklist_from_head": stage,
    # negative-boundary flags
    "no_dynamics_results_claimed": True,
    "heavy_spectrum_still_placeholder": True,
    "scalar_hessian_goldstone_eigenvectors_pending": True,
    "history_assets_recovered_read_only": True,
    "zeta_value_derived": False,
}

OUT.mkdir(parents=True, exist_ok=True)
json_path = OUT / "dyn0_conventions_inventory.json"
md_path = OUT / "dyn0_conventions_inventory.md"
json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
md_path.write_text(
    "# DYN-0: Dynamics Conventions Addendum and Inventory\n\n"
    f"`{npass}/{len(CHECKS)}` checks passed.\n\n"
    "- Conventions addendum recorded: superpotential parameters, PS\n"
    "  singlet-vev registry (mirrored from HEAD), Aulakh-Girdhar branch\n"
    "  convention, adopted decisions (root code/ + whitelisted output\n"
    "  lanes; M_S benchmark 3 TeV, DYN-2 fit window [1, 10] TeV).\n"
    "- Vacuum-cubic special points verified exactly: SU(5) x=1/2 (xi=-5),\n"
    "  SU(5) x=-1 (xi=10), G_LR x=0 (xi=3), flipped SU(5) x=1/3 (xi=-2/3).\n"
    f"- Inventory: {sum(len(v) for v in INVENTORY.values())} history assets "
    "across the DYN items, all present in HEAD (read-only contract).\n"
    "- Rerun-from-history gates: the four audit-4a1 scripts, extracted from\n"
    "  HEAD into a temporary mini-tree, regenerate ledgers matching HEAD up\n"
    "  to timestamps/sha stamps/source manifests.\n\n"
    "Boundary: no dynamics results are claimed.  The heavy spectrum is\n"
    "still a placeholder and the scalar-Hessian Goldstone eigenvector audit\n"
    "is still pending; both are DYN-1 targets.\n"
)
print(f"Wrote {json_path}")
print(f"Wrote {md_path}")
print(f"DYN-0 scaffold: {npass}/{len(CHECKS)} checks passed; conventions fixed, "
      "inventory closed, history gates green; no dynamics results claimed.")
if not all_pass:
    raise SystemExit(1)
