#!/usr/bin/env python3
"""Compatibility entry point for the canonical Route-E dynamics implementation."""

from _route_e_dyn_compat import delegate

raise SystemExit(delegate(__file__))
