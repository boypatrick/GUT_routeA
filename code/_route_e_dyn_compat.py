"""Compatibility delegation for the historical root DYN entry points."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path


def delegate(current_file: str) -> int:
    """Execute the single canonical implementation under route_E/code_dyn."""

    current = Path(current_file).resolve()
    target = current.parents[1] / "route_E" / "code_dyn" / current.name
    if not target.is_file():
        raise FileNotFoundError(
            f"canonical Route-E dynamics implementation is missing: {target}"
        )
    if target.samefile(current):
        raise RuntimeError(f"refusing recursive dynamics delegation: {target}")

    target_dir = str(target.parent)
    sys.path.insert(0, target_dir)
    try:
        runpy.run_path(str(target), run_name="__main__")
    finally:
        if sys.path and sys.path[0] == target_dir:
            sys.path.pop(0)
        else:
            try:
                sys.path.remove(target_dir)
            except ValueError:
                pass
    return 0
