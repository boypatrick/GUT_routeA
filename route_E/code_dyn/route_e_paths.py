"""Canonical, case-sensitive paths shared by the Route-E DYN audits.

The DYN sources were moved from the repository-level ``code/`` directory to
``route_E/code_dyn/``.  Deriving every path from ``parents[1]`` therefore no
longer identifies the repository root.  Keep the layout knowledge here so the
individual numerical audits remain independent of their current working
directory.
"""

import os
from pathlib import Path


_CODE_DYN_ROOT = Path(__file__).resolve().parent
ROUTE_E_ROOT = _CODE_DYN_ROOT.parent
REPO_ROOT = ROUTE_E_ROOT.parent

# Route-E's recovered, read-only publication inputs live in this tree.
ARCHIVE_OUTPUT = ROUTE_E_ROOT / "output"

# Generated DYN ledgers retain the repository-level output/audit* convention
# by default.  The environment override is used by the isolated DAG runner
# and makes direct lane-level smoke tests non-destructive when requested.
AUDIT_OUTPUT = Path(
    os.environ.get("ROUTE_E_OUTPUT_ROOT", str(REPO_ROOT / "output"))
).expanduser().resolve()

# Expensive deterministic tensor caches are never global /tmp files.  A
# content-addressing suffix is owned by each consumer; this module only fixes
# a run-local/cache-root contract.
CACHE_DIR = Path(
    os.environ.get("ROUTE_E_CACHE_DIR", str(REPO_ROOT / "tmp" / "route_e_dyn_cache"))
).expanduser().resolve()


if ROUTE_E_ROOT.name != "route_e":
    raise RuntimeError(
        "Route-E path resolution requires the exact lowercase directory "
        f"name 'route_e'; resolved {ROUTE_E_ROOT!s}"
    )
