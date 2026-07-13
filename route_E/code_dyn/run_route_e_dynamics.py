#!/usr/bin/env python3
"""Fail-closed DAG runner for the recovered Route-E dynamics audits.

The 19 scientific scripts and their path resolver are copied byte-for-byte
into an isolated ``workspace/route_E/code_dyn`` mirror before execution.
Their ``AUDIT_OUTPUT`` therefore resolves to ``workspace/output`` and cannot
overwrite repository ledgers.  Inputs are copied as a minimal snapshot, while
Git-history reads are pinned to the recorded commit by
``GIT_DIR``/``GIT_WORK_TREE``.

Profiles:

* ``quick``: dependency smoke path ending at DYN-1b, DYN-4c, and DYN-9b-1;
* ``lane``: one or more explicit ``--lane`` targets plus transitive inputs;
* ``full``: all 19 recovered scripts plus two fail-closed validity guards,
  including the DYN-8 collector (21 nodes total).

The added guards are labelled ``guard`` rather than ``recovered`` in the
manifest.  Selecting DYN-5 targets DYN-5V automatically; selecting DYN-7 or
DYN-9b-3 likewise includes DYN-7F.

RE-SC3/4/5 are optional theory lanes but conditional hard dependencies of the
scripts which read them.  In ``--string-cards optional`` mode a missing card
blocks its consumer and descendants without crashing.  In ``required`` mode
the same absence fails preflight and no scientific script starts.

Exit codes are 0 for a mechanically complete run (or complete dry-run plan),
1 for a script/ledger failure, and 2 for preflight failure or an explicitly
incomplete optional run.  A ledger's ``all_pass`` value is recorded only as
mechanical evidence.  The top-level outcome is therefore
``mechanical_success``, never an unqualified physics success.  Physics
interpretation is loaded independently from
``dyn_claim_registry.json`` and is never promoted by this runner.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1
RUNNER_DIR = Path(__file__).resolve().parent
DEFAULT_REPO_ROOT = RUNNER_DIR.parents[1]
DEFAULT_REGISTRY = RUNNER_DIR / "dyn_claim_registry.json"
QUICK_TARGETS = ("dyn-1b", "dyn-4c", "dyn-9b-1")


@dataclass(frozen=True)
class Node:
    node_id: str
    lane_id: str
    script: str
    deps: tuple[str, ...]
    ledger: str
    outputs: tuple[str, ...]
    inputs: tuple[str, ...] = ()
    external_deps: tuple[str, ...] = ()
    node_kind: str = "recovered"


def _outputs(prefix: str) -> tuple[str, str]:
    return (prefix + ".json", prefix + ".md")


NODES: tuple[Node, ...] = (
    Node("dyn-0", "DYN-0", "audit0_dyn_conventions_and_inventory.py", (),
         "output/audit0/dyn0_conventions_inventory.json",
         _outputs("output/audit0/dyn0_conventions_inventory")),
    Node("dyn-1a", "DYN-1a", "audit4a_dyn1a_vacuum_goldstone_audit.py",
         ("dyn-0",), "output/audit4a/dyn1a_vacuum_goldstone.json",
         _outputs("output/audit4a/dyn1a_vacuum_goldstone")
         + ("output/audit4a/heavy_spectrum.json",)),
    Node("dyn-1b", "DYN-1b", "audit4a_dyn1b_full_spectrum.py",
         ("dyn-1a",), "output/audit4a/dyn1b_full_spectrum.json",
         _outputs("output/audit4a/dyn1b_full_spectrum")
         + ("output/audit4a/heavy_spectrum.json",)),
    Node("dyn-2", "DYN-2", "audit3_dyn2_thresholds_unification.py",
         ("dyn-1b",), "output/audit3/dyn2_thresholds_unification.json",
         _outputs("output/audit3/dyn2_thresholds_unification")),
    Node("dyn-3", "DYN-3", "audit2_dyn3_proton_d5_kill_criterion.py",
         ("dyn-2",), "output/audit2/dyn3_proton_d5_kill_criterion.json",
         _outputs("output/audit2/dyn3_proton_d5_kill_criterion"),
         ("full_knu", "knu_target", "closure_card")),
    Node("dyn-2b", "DYN-2b", "audit3_dyn2b_rescue_scan.py",
         ("dyn-2", "dyn-3"), "output/audit3/dyn2b_rescue_scan.json",
         _outputs("output/audit3/dyn2b_rescue_scan")),
    Node("dyn-4a", "DYN-4a", "audit1_dyn4a_seesaw_replay_zeta_posterior.py",
         ("dyn-0",), "output/audit1/dyn4a_seesaw_zeta_posterior.json",
         _outputs("output/audit1/dyn4a_seesaw_zeta_posterior"),
         ("closure_card", "majorana_rank")),
    Node("dyn-4b", "DYN-4b", "audit1_dyn4b_refreshed_card_unconditional_zeta.py",
         ("dyn-4a",), "output/audit1/dyn4b_unconditional_zeta.json",
         _outputs("output/audit1/dyn4b_unconditional_zeta"),
         ("closure_card",)),
    Node("dyn-4c", "DYN-4c", "audit1_dyn4c_kernel_dirac_refit.py",
         ("dyn-4a", "dyn-4b"), "output/audit1/dyn4c_kernel_dirac_refit.json",
         _outputs("output/audit1/dyn4c_kernel_dirac_refit"),
         ("closure_card",)),
    Node("dyn-5", "DYN-5", "audit5_dyn5_messenger_one_loop.py",
         ("dyn-4a", "dyn-4b"), "output/audit5/dyn5_messenger_one_loop.json",
         _outputs("output/audit5/dyn5_messenger_one_loop"),
         ("closure_card",)),
    Node("dyn-5v", "DYN-5V", "audit5_dyn5_model_validity.py",
         ("dyn-5",), "output/audit5/dyn5_model_validity.json",
         _outputs("output/audit5/dyn5_model_validity"),
         node_kind="guard"),
    Node("dyn-7", "DYN-7", "audit7_dyn7_leptogenesis_argzeta.py",
         ("dyn-4b",), "output/audit7/dyn7_leptogenesis_argzeta.json",
         _outputs("output/audit7/dyn7_leptogenesis_argzeta"),
         ("closure_card",)),
    Node("dyn-7f", "DYN-7F", "audit7_dyn7_flavor_regime_gate.py",
         ("dyn-7",), "output/audit7/dyn7_flavor_regime_gate.json",
         _outputs("output/audit7/dyn7_flavor_regime_gate"),
         node_kind="guard"),
    Node("dyn-9", "DYN-9", "audit9_dyn9_nonsusy_intermediate.py",
         ("dyn-2b",), "output/audit9/dyn9_nonsusy_intermediate.json",
         _outputs("output/audit9/dyn9_nonsusy_intermediate")),
    Node("dyn-9b-1", "DYN-9b-1", "audit9_dyn9b1_nonsusy_vacuum_thresholds.py",
         ("dyn-9",), "output/audit9/dyn9b1_nonsusy_vacuum_thresholds.json",
         _outputs("output/audit9/dyn9b1_nonsusy_vacuum_thresholds")),
    Node("dyn-9b-1b", "DYN-9b-1b", "audit9_dyn9b1b_210_quartic_descent.py",
         ("dyn-9b-1",), "output/audit9/dyn9b1b_210_quartic_descent.json",
         _outputs("output/audit9/dyn9b1b_210_quartic_descent")),
    Node("dyn-9b-1c", "DYN-9b-1c", "audit9_dyn9b1c_eps_and_126_quartics.py",
         ("dyn-9b-1b",), "output/audit9/dyn9b1c_eps_and_126_quartics.json",
         _outputs("output/audit9/dyn9b1c_eps_and_126_quartics")),
    Node("dyn-9b-1d", "DYN-9b-1d", "audit9_dyn9b1d_lr_ratio_scan.py",
         ("dyn-9b-1c",), "output/audit9/dyn9b1d_lr_ratio_scan.json",
         _outputs("output/audit9/dyn9b1d_lr_ratio_scan")),
    Node("dyn-9b-2", "DYN-9b-2", "audit9_dyn9b2_nonsusy_flavor_refit.py",
         ("dyn-4a", "dyn-4b", "dyn-9b-1d"),
         "output/audit9/dyn9b2_nonsusy_flavor_refit.json",
         _outputs("output/audit9/dyn9b2_nonsusy_flavor_refit"),
         ("closure_card",), ("RE-SC3",)),
    Node("dyn-9b-3", "DYN-9b-3", "audit9_dyn9b3_nonsusy_leptogenesis.py",
         ("dyn-7f", "dyn-9b-2"),
         "output/audit9/dyn9b3_nonsusy_leptogenesis.json",
         _outputs("output/audit9/dyn9b3_nonsusy_leptogenesis"),
         ("closure_card",)),
    Node("dyn-8", "DYN-8", "audit8_dyn8_falsifiability_collection.py",
         ("dyn-2", "dyn-2b", "dyn-3", "dyn-4b", "dyn-4c", "dyn-5v",
          "dyn-7f", "dyn-9", "dyn-9b-1", "dyn-9b-1b", "dyn-9b-1c",
          "dyn-9b-1d", "dyn-9b-2", "dyn-9b-3"),
         "output/audit8/dyn8_falsifiability_collection.json",
         _outputs("output/audit8/dyn8_falsifiability_collection"),
         ("route_e_core", "route_e_closure", "route_e_tex", "paper_tex"),
         ("RE-SC3", "RE-SC4", "RE-SC5")),
)

NODE_BY_ID = {node.node_id: node for node in NODES}

# Only the immutable/read-only inputs actually opened outside generated output.
ARCHIVAL_INPUTS: dict[str, tuple[str, str]] = {
    "closure_card": (
        "route_E/output/publication_closure_card/publication_closure_card.json",
        "route_E/output/publication_closure_card/publication_closure_card.json"),
    "majorana_rank": (
        "route_E/output/source_majorana_texture_rank/summary.json",
        "route_E/output/source_majorana_texture_rank/summary.json"),
    "full_knu": ("route_E/output/full_knu_width/summary.json",
                 "route_E/output/full_knu_width/summary.json"),
    "knu_target": ("route_E/output/knu_target_map/summary.json",
                   "route_E/output/knu_target_map/summary.json"),
    "route_e_core": ("route_E/output/route_e_first_principles.json",
                     "route_E/output/route_e_first_principles.json"),
    "route_e_closure": ("route_E/output/route_e_dependency_closure.json",
                        "route_E/output/route_e_dependency_closure.json"),
    "route_e_tex": ("route_E/tex/route_e_first_principles.tex",
                    "route_E/tex/route_e_first_principles.tex"),
    "paper_tex": ("paper/gut_framework.tex", "paper/gut_framework.tex"),
}

EXTERNAL_ARTIFACTS: dict[str, tuple[str, str]] = {
    "RE-SC3": ("route_d/output/d3_instanton_majorana_pricing.json",
               "route_d/output/d3_instanton_majorana_pricing.json"),
    "RE-SC4": ("route_d/output/d4_stueckelberg_protection.json",
               "route_d/output/d4_stueckelberg_protection.json"),
    "RE-SC5": ("route_d/output/d5_susy_breaking_bridge.json",
               "route_d/output/d5_susy_breaking_bridge.json"),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    data = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    with temp.open("wb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def run_git(repo_root: Path, *args: str) -> tuple[int, str]:
    proc = subprocess.run(["git", *args], cwd=repo_root, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          check=False)
    return proc.returncode, proc.stdout.strip()


def git_provenance(repo_root: Path) -> dict[str, Any]:
    rc_head, head = run_git(repo_root, "rev-parse", "HEAD")
    rc_status, status = run_git(repo_root, "status", "--porcelain=v1",
                                "--untracked-files=normal")
    return {
        "commit": head if rc_head == 0 else None,
        "head_query_exit_code": rc_head,
        "dirty": bool(status) if rc_status == 0 else None,
        "status_porcelain": status.splitlines() if status else [],
        "status_sha256": sha256_bytes((status + "\n").encode()),
        "status_query_exit_code": rc_status,
    }


def runtime_provenance() -> dict[str, Any]:
    result: dict[str, Any] = {
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
    }
    try:
        import numpy as np
        result["numpy_version"] = np.__version__
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            np.show_config()
        config = buffer.getvalue()
        result["numpy_config"] = config
        result["numpy_config_sha256"] = sha256_bytes(config.encode())
    except Exception as exc:  # pragma: no cover - reported, scripts then fail
        result["numpy_import_error"] = repr(exc)
    try:
        import scipy
        result["scipy_version"] = scipy.__version__
    except Exception as exc:  # pragma: no cover - reported, relevant lane fails
        result["scipy_import_error"] = repr(exc)
    return result


def load_claim_registry(path: Path) -> tuple[dict[str, Any], str]:
    raw = path.read_bytes()
    registry = json.loads(raw)
    claims = registry.get("claims", {})
    missing = sorted(set(NODE_BY_ID) - set(claims))
    extra = sorted(set(claims) - set(NODE_BY_ID))
    allowed = set(registry.get("allowed_statuses", []))
    invalid = sorted(node_id for node_id, claim in claims.items()
                     if claim.get("status") not in allowed)
    if missing or extra or invalid:
        raise ValueError("invalid claim registry: "
                         f"missing={missing}, extra={extra}, "
                         f"invalid_status={invalid}")
    return registry, sha256_bytes(raw)


def dependency_closure(targets: Iterable[str]) -> set[str]:
    selected: set[str] = set()

    def visit(node_id: str) -> None:
        if node_id in selected:
            return
        if node_id not in NODE_BY_ID:
            raise KeyError(f"unknown lane {node_id!r}")
        for dep in NODE_BY_ID[node_id].deps:
            visit(dep)
        selected.add(node_id)

    for target in targets:
        visit(target)
    return selected


def normalized_lane(raw: str) -> str:
    lane = raw.strip().lower().replace("_", "-")
    lane = re.sub(r"^dyn(?=[0-9])", "dyn-", lane)
    if lane not in NODE_BY_ID:
        raise argparse.ArgumentTypeError(
            f"unknown lane {raw!r}; use --list-lanes")
    return lane


def verify_static_dag() -> None:
    seen: set[str] = set()
    for node in NODES:
        missing = set(node.deps) - seen
        if missing:
            raise RuntimeError(
                f"DAG order invalid at {node.node_id}: {sorted(missing)}")
        seen.add(node.node_id)
    recovered = sum(node.node_kind == "recovered" for node in NODES)
    guards = sum(node.node_kind == "guard" for node in NODES)
    if (recovered, guards) != (19, 2):
        raise RuntimeError(
            f"node inventory changed: recovered={recovered}, guards={guards}")


def selected_nodes(profile: str, lanes: list[str]) -> list[Node]:
    if profile == "lane":
        if not lanes:
            raise ValueError("--profile lane requires at least one --lane")
        guard_targets = {"dyn-5": "dyn-5v", "dyn-7": "dyn-7f"}
        targets = [guard_targets.get(lane, lane) for lane in lanes]
    elif lanes:
        raise ValueError("--lane is only valid with --profile lane")
    elif profile == "quick":
        targets = list(QUICK_TARGETS)
    else:
        targets = list(NODE_BY_ID)
    selected = dependency_closure(targets)
    return [node for node in NODES if node.node_id in selected]


def copy_snapshot(source: Path, destination: Path) -> dict[str, Any]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    destination.chmod(0o444)
    digest = sha256_file(destination)
    if digest != sha256_file(source):
        raise OSError(f"snapshot digest mismatch: {source}")
    return {
        "source": str(source),
        "staged": str(destination),
        "bytes": destination.stat().st_size,
        "sha256": digest,
    }


def rel_to(path: Path, base: Path) -> str:
    try:
        return path.relative_to(base).as_posix()
    except ValueError:
        return str(path)


def validate_outputs(workspace: Path, node: Node) -> tuple[list[str],
                                                               dict[str, Any],
                                                               bool | None]:
    errors: list[str] = []
    artifacts: dict[str, Any] = {}
    for relative in node.outputs:
        path = workspace / relative
        if not path.is_file():
            errors.append(f"missing expected output: {relative}")
            continue
        if path.stat().st_size == 0:
            errors.append(f"empty expected output: {relative}")
        artifacts[relative] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }

    ledger_all_pass: bool | None = None
    ledger_path = workspace / node.ledger
    if ledger_path.is_file():
        try:
            ledger = json.loads(ledger_path.read_text())
            value = ledger.get("all_pass")
            ledger_all_pass = value if isinstance(value, bool) else None
            if ledger_all_pass is not True:
                errors.append("ledger all_pass is not literal true")
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"ledger JSON invalid: {exc}")
    return errors, artifacts, ledger_all_pass


def update_manifest(path: Path, manifest: dict[str, Any]) -> None:
    manifest["manifest_updated_at"] = utc_now()
    atomic_write_json(path, manifest)


def summarize_physics(records: dict[str, dict[str, Any]]) -> dict[str, Any]:
    counts: dict[str, int] = {}
    blocking_claims: dict[str, list[str]] = {}
    conditionally_promotable = {"proved", "conditional", "no-go"}
    for node_id, record in records.items():
        status = record["physics_status"]
        counts[status] = counts.get(status, 0) + 1
        blockers = [item["id"] for item in
                    record["physics_assessment"].get("blockers", [])
                    if item.get("severity") == "blocking"]
        if blockers:
            blocking_claims[node_id] = blockers
    non_promotable = sorted(
        node_id for node_id, record in records.items()
        if record["physics_status"] not in conditionally_promotable)
    return {
        "status_counts": dict(sorted(counts.items())),
        "blocking_claims": blocking_claims,
        "non_promotable_nodes": non_promotable,
        "physics_promotion_allowed": not non_promotable,
        "note": "Mechanical all_pass never changes this registry-derived summary."
    }


def execute_node(node: Node, record: dict[str, Any], workspace: Path,
                 run_dir: Path, env: dict[str, str], timeout: float,
                 python: str) -> bool:
    stdout_path = run_dir / "logs" / f"{node.node_id}.stdout.log"
    stderr_path = run_dir / "logs" / f"{node.node_id}.stderr.log"
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    command = [python, str(workspace / "route_e" / "code_dyn" / node.script)]
    record.update({
        "mechanical_status": "running",
        "command": command,
        "started_at": utc_now(),
    })
    start = time.monotonic()
    timed_out = False
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        try:
            proc = subprocess.run(command, cwd=workspace, env=env,
                                  stdout=stdout, stderr=stderr,
                                  timeout=timeout, check=False)
            returncode: int | None = proc.returncode
        except subprocess.TimeoutExpired:
            returncode = 124
            timed_out = True
    elapsed = time.monotonic() - start

    errors, artifacts, ledger_all_pass = validate_outputs(workspace, node)
    if returncode != 0:
        errors.insert(0, f"nonzero exit code: {returncode}")
    record.update({
        "finished_at": utc_now(),
        "elapsed_seconds": round(elapsed, 6),
        "exit_code": returncode,
        "timed_out": timed_out,
        "stdout": {
            "path": rel_to(stdout_path, run_dir),
            "bytes": stdout_path.stat().st_size,
            "sha256": sha256_file(stdout_path),
        },
        "stderr": {
            "path": rel_to(stderr_path, run_dir),
            "bytes": stderr_path.stat().st_size,
            "sha256": sha256_file(stderr_path),
        },
        "ledger_all_pass": ledger_all_pass,
        "output_artifacts": artifacts,
        "validation_errors": errors,
        "mechanical_status": "passed" if not errors else (
            "failed_timeout" if timed_out else "failed"),
    })
    return not errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("quick", "lane", "full"),
                        default="quick")
    parser.add_argument("--lane", action="append", type=normalized_lane,
                        default=[], help="Exact target such as DYN-4c; repeatable")
    parser.add_argument("--list-lanes", action="store_true")
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT)
    parser.add_argument("--output-root", type=Path,
                        default=Path(tempfile.gettempdir()) / "route_e_dyn_runs",
                        help="Parent for a newly created run-id directory")
    parser.add_argument("--run-id", help="Safe unique directory name")
    parser.add_argument("--claim-registry", type=Path,
                        default=DEFAULT_REGISTRY)
    parser.add_argument("--string-cards", choices=("optional", "required"),
                        default="optional")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--timeout-seconds", type=float, default=7200.0)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--python-hash-seed", default="0")
    parser.add_argument(
        "--allow-legacy-global-cache", action="store_true",
        help="Permit scripts still containing hard-coded Path('/tmp') caches")
    return parser


def main(argv: list[str] | None = None) -> int:
    verify_static_dag()
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.list_lanes:
        for node in NODES:
            print(f"{node.lane_id:10s} [{node.node_kind:9s}] {node.script}")
        return 0
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")

    try:
        nodes = selected_nodes(args.profile, args.lane)
    except ValueError as exc:
        parser.error(str(exc))

    repo_root = args.repo_root.expanduser().resolve()
    source_dir = repo_root / "route_e" / "code_dyn"
    registry_path = args.claim_registry.expanduser().resolve()
    try:
        registry, registry_sha = load_claim_registry(registry_path)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        print(f"claim registry error: {exc}", file=sys.stderr)
        return 2

    run_id = args.run_id or (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + f"-{os.getpid()}")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}", run_id):
        parser.error("--run-id must match [A-Za-z0-9][A-Za-z0-9._-]{0,127}")
    output_root = args.output_root.expanduser().resolve()
    run_dir = output_root / run_id
    try:
        run_dir.relative_to(source_dir)
    except ValueError:
        pass
    else:
        parser.error("--output-root/run-id must not be inside route_E/code_dyn")
    output_root.mkdir(parents=True, exist_ok=True)
    try:
        run_dir.mkdir()
    except FileExistsError:
        print(f"refusing to reuse existing run directory: {run_dir}",
              file=sys.stderr)
        return 2

    workspace = run_dir / "workspace"
    manifest_path = run_dir / "run_manifest.json"
    claims = registry["claims"]
    node_records: dict[str, dict[str, Any]] = {}
    for node in nodes:
        claim = claims[node.node_id]
        node_records[node.node_id] = {
            "lane_id": node.lane_id,
            "node_kind": node.node_kind,
            "script": f"route_E/code_dyn/{node.script}",
            "dependencies": list(node.deps),
            "external_dependencies": list(node.external_deps),
            "expected_outputs": list(node.outputs),
            "mechanical_status": "pending",
            "ledger_all_pass": None,
            "physics_status": claim["status"],
            "physics_assessment": claim,
        }

    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "runner": "route_E/code_dyn/run_route_e_dynamics.py",
        "runner_sha256": sha256_file(Path(__file__).resolve()),
        "claim_registry": {
            "path": rel_to(registry_path, repo_root),
            "sha256": registry_sha,
            "reviewed_at": registry.get("reviewed_at"),
        },
        "run_id": run_id,
        "run_directory": str(run_dir),
        "workspace": str(workspace),
        "started_at": utc_now(),
        "finished_at": None,
        "profile": args.profile,
        "requested_lanes": args.lane,
        "selected_order": [node.node_id for node in nodes],
        "node_inventory": {
            "total": len(NODES),
            "recovered_scripts": sum(
                node.node_kind == "recovered" for node in NODES),
            "added_validity_guards": sum(
                node.node_kind == "guard" for node in NODES),
            "selected_recovered": sum(
                node.node_kind == "recovered" for node in nodes),
            "selected_guards": sum(
                node.node_kind == "guard" for node in nodes),
        },
        "string_card_policy": args.string_cards,
        "dry_run": args.dry_run,
        "timeout_seconds": args.timeout_seconds,
        "legacy_global_cache_allowed": args.allow_legacy_global_cache,
        "repo_root": str(repo_root),
        "git": git_provenance(repo_root),
        "runtime": runtime_provenance(),
        "environment_controls": {
            "PYTHONHASHSEED": args.python_hash_seed,
            "TMPDIR": "workspace/tmp",
            "ROUTE_E_REPO_ROOT": str(repo_root),
            "ROUTE_E_OUTPUT_ROOT": "workspace/output",
            "ROUTE_E_CACHE_DIR": "workspace/cache",
            "note": "Scripts may not consume all controls; staged route_e_paths.py still resolves every declared output inside the workspace."
        },
        "preflight_errors": [],
        "preflight_warnings": [],
        "source_files": {},
        "input_snapshot": {},
        "external_artifacts": {},
        "nodes": node_records,
        "physics_summary": summarize_physics(node_records),
        "outcome": "preflight",
    }
    update_manifest(manifest_path, manifest)

    # ------------------------------- fail-closed preflight before any script
    required_input_ids = sorted({item for node in nodes for item in node.inputs})
    external_ids = sorted({item for node in nodes
                           for item in node.external_deps})
    resolver_source = source_dir / "route_e_paths.py"
    if not resolver_source.is_file():
        manifest["preflight_errors"].append(
            f"missing shared path resolver: {resolver_source}")
    for node in nodes:
        source = source_dir / node.script
        if not source.is_file():
            manifest["preflight_errors"].append(
                f"missing source for {node.node_id}: {source}")
            continue
        source_bytes = source.read_bytes()
        manifest["source_files"][node.node_id] = {
            "path": rel_to(source, repo_root),
            "bytes": len(source_bytes),
            "sha256": sha256_bytes(source_bytes),
        }
        if (b'Path("/tmp")' in source_bytes
                or b"Path('/tmp')" in source_bytes):
            manifest["source_files"][node.node_id][
                "hard_coded_global_tmp_cache"] = True
            if not args.allow_legacy_global_cache:
                node_records[node.node_id]["mechanical_status"] = (
                    "blocked_unsafe_global_cache")
                node_records[node.node_id]["blocking_reason"] = (
                    "hard-coded /tmp cache; rerun only with explicit "
                    "--allow-legacy-global-cache or repair the script")

    for input_id in required_input_ids:
        source_rel, _ = ARCHIVAL_INPUTS[input_id]
        source = repo_root / source_rel
        if not source.is_file():
            manifest["preflight_errors"].append(
                f"missing required archival input {input_id}: {source_rel}")

    missing_optional: set[str] = set()
    for external_id in external_ids:
        source_rel, _ = EXTERNAL_ARTIFACTS[external_id]
        source = repo_root / source_rel
        present = source.is_file()
        manifest["external_artifacts"][external_id] = {
            "source": source_rel,
            "present": present,
            "policy": args.string_cards,
            "consumers": [node.node_id for node in nodes
                          if external_id in node.external_deps],
            "sha256": sha256_file(source) if present else None,
        }
        if not present and args.string_cards == "required":
            manifest["preflight_errors"].append(
                f"missing required external artifact {external_id}: {source_rel}")
        elif not present:
            missing_optional.add(external_id)
            manifest["preflight_warnings"].append(
                f"optional {external_id} absent; its consumers are blocked")

    for node in nodes:
        missing = sorted(set(node.external_deps) & missing_optional)
        if missing:
            node_records[node.node_id]["mechanical_status"] = (
                "blocked_missing_optional_external")
            node_records[node.node_id]["blocking_external_dependencies"] = missing

    if manifest["preflight_errors"]:
        for record in node_records.values():
            if record["mechanical_status"] == "pending":
                record["mechanical_status"] = "not_run_preflight_failure"
        manifest["outcome"] = "preflight_failed"
        manifest["finished_at"] = utc_now()
        update_manifest(manifest_path, manifest)
        print(f"PRE-FLIGHT FAILED; manifest: {manifest_path}", file=sys.stderr)
        return 2

    # ------------------------------------------------------- isolated snapshot
    staged_code_dir = workspace / "route_e" / "code_dyn"
    staged_code_dir.mkdir(parents=True)
    (workspace / "output").mkdir()
    (workspace / "tmp").mkdir()
    (workspace / "cache").mkdir()
    for node in nodes:
        source = source_dir / node.script
        staged = staged_code_dir / node.script
        copied = copy_snapshot(source, staged)
        copied["source"] = rel_to(source, repo_root)
        copied["staged"] = rel_to(staged, run_dir)
        manifest["source_files"][node.node_id]["staged"] = copied["staged"]
        manifest["source_files"][node.node_id]["staged_sha256"] = copied["sha256"]

    resolver_staged = staged_code_dir / "route_e_paths.py"
    resolver_copy = copy_snapshot(resolver_source, resolver_staged)
    manifest["path_resolver"] = {
        "source": rel_to(resolver_source, repo_root),
        "staged": rel_to(resolver_staged, run_dir),
        "bytes": resolver_copy["bytes"],
        "sha256": resolver_copy["sha256"],
    }

    for input_id in required_input_ids:
        source_rel, stage_rel = ARCHIVAL_INPUTS[input_id]
        copied = copy_snapshot(repo_root / source_rel, workspace / stage_rel)
        copied["source"] = source_rel
        copied["staged"] = rel_to(workspace / stage_rel, run_dir)
        manifest["input_snapshot"][input_id] = copied
    for external_id in external_ids:
        source_rel, stage_rel = EXTERNAL_ARTIFACTS[external_id]
        source = repo_root / source_rel
        if source.is_file():
            copied = copy_snapshot(source, workspace / stage_rel)
            copied["source"] = source_rel
            copied["staged"] = rel_to(workspace / stage_rel, run_dir)
            manifest["external_artifacts"][external_id]["staged"] = copied

    update_manifest(manifest_path, manifest)

    if args.dry_run:
        for node in nodes:
            record = node_records[node.node_id]
            if record["mechanical_status"].startswith("blocked_"):
                continue
            blocked_deps = [dep for dep in node.deps
                            if dep in node_records and
                            node_records[dep]["mechanical_status"].startswith(
                                "blocked_")]
            if blocked_deps:
                record["mechanical_status"] = "blocked_upstream"
                record["blocking_dependencies"] = blocked_deps
            else:
                record["mechanical_status"] = "planned"
        incomplete = any(record["mechanical_status"].startswith("blocked_")
                         for record in node_records.values())
        manifest["outcome"] = (
            "dry_run_incomplete" if incomplete else "dry_run_ready")
        manifest["finished_at"] = utc_now()
        update_manifest(manifest_path, manifest)
        print(f"{manifest['outcome']}; manifest: {manifest_path}")
        return 2 if incomplete else 0

    # ---------------------------------------------------------- sequential DAG
    env = os.environ.copy()
    env.update({
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONHASHSEED": args.python_hash_seed,
        "TMPDIR": str(workspace / "tmp"),
        "ROUTE_E_REPO_ROOT": str(repo_root),
        "ROUTE_E_OUTPUT_ROOT": str(workspace / "output"),
        "ROUTE_E_CACHE_DIR": str(workspace / "cache"),
        "GIT_DIR": str(repo_root / ".git"),
        "GIT_WORK_TREE": str(repo_root),
    })
    aborted = False
    for node in nodes:
        record = node_records[node.node_id]
        if aborted:
            if record["mechanical_status"] == "pending":
                record["mechanical_status"] = "not_run_after_failure"
            continue
        if record["mechanical_status"].startswith("blocked_"):
            continue
        blocked_deps = [dep for dep in node.deps
                        if dep in node_records and
                        node_records[dep]["mechanical_status"] != "passed"]
        if blocked_deps:
            record["mechanical_status"] = "blocked_upstream"
            record["blocking_dependencies"] = blocked_deps
            update_manifest(manifest_path, manifest)
            continue

        print(f"[{node.lane_id}] running {node.script}", flush=True)
        update_manifest(manifest_path, manifest)
        ok = execute_node(node, record, workspace, run_dir, env,
                          args.timeout_seconds, args.python)
        update_manifest(manifest_path, manifest)
        if not ok:
            # Global abort is stricter than merely skipping descendants: no
            # independent lane runs after a failed numerical assertion.
            aborted = True
            print(f"[{node.lane_id}] FAILED; aborting remaining DAG",
                  file=sys.stderr, flush=True)
        else:
            print(f"[{node.lane_id}] passed mechanically; physics status: "
                  f"{record['physics_status']}", flush=True)

    statuses = [record["mechanical_status"] for record in node_records.values()]
    if any(status.startswith("failed") for status in statuses):
        manifest["outcome"] = "failed"
        exit_code = 1
    elif all(status == "passed" for status in statuses):
        manifest["outcome"] = "mechanical_success"
        exit_code = 0
    else:
        manifest["outcome"] = "incomplete"
        exit_code = 2
    manifest["finished_at"] = utc_now()
    update_manifest(manifest_path, manifest)
    print(f"{manifest['outcome']}; manifest: {manifest_path}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
