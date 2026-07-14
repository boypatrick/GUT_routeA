# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

A physics research repository (not a software product): a lean "Route-A" paper for a face-projection Spin(10) GUT effective framework, plus verification sandboxes. The single source of truth for all claims is `paper/gut_framework.tex` (RevTeX 4.2, PRD preprint style). Everything else exists to verify, or to bound, what the paper is allowed to claim.

## Commands

Build the paper (run inside `paper/`):

```
pdflatex gut_framework.tex
bibtex gut_framework
pdflatex gut_framework.tex
pdflatex gut_framework.tex
```

Bibliography is `refs.bib` with `apsrev4-2` style. LaTeX build products (`*.aux`, `*.log`, `*.out`, `*.blg`, `*Notes.bib`) are gitignored.

Run the Route-D / Route-E verification audits (plain python3; numpy is available, sympy is NOT installed):

```
python3 route_d/code/verify_d1_f_theory_local_placement.py
python3 route_d/code/verify_d2_e3_instanton_zeta.py
python3 route_d/code/verify_d3_instanton_majorana_pricing.py
python3 route_d/code/verify_d4_stueckelberg_protection.py
python3 route_d/code/verify_d5_susy_breaking_bridge.py
python3 route_E/code/verify_route_e_first_principles.py
python3 route_E/code/verify_route_e_dependency_closure.py
python3 code/audit0_dyn_conventions_and_inventory.py
python3 code/audit4a_dyn1a_vacuum_goldstone_audit.py
python3 code/audit4a_dyn1b_full_spectrum.py
python3 code/audit3_dyn2_thresholds_unification.py
python3 code/audit1_dyn4a_seesaw_replay_zeta_posterior.py
python3 code/audit1_dyn4b_refreshed_card_unconditional_zeta.py
python3 code/audit7_dyn7_leptogenesis_argzeta.py
python3 code/audit2_dyn3_proton_d5_kill_criterion.py
python3 code/audit3_dyn2b_rescue_scan.py
python3 code/audit9_dyn9_nonsusy_intermediate.py
python3 code/audit5_dyn5_messenger_one_loop.py
python3 code/audit8_dyn8_falsifiability_collection.py
python3 code/audit9_dyn9b1_nonsusy_vacuum_thresholds.py
python3 code/audit9_dyn9b1b_210_quartic_descent.py
python3 code/audit9_dyn9b1c_eps_and_126_quartics.py
python3 code/audit9_dyn9b2_nonsusy_flavor_refit.py
python3 code/audit9_dyn9b3_nonsusy_leptogenesis.py
python3 code/audit9_dyn9b1d_lr_ratio_scan.py
python3 code/audit1_dyn4c_kernel_dirac_refit.py
```

The DYN-4a and DYN-5 scripts require the archival output set restored under `route_E/output/` (user-managed; sha256 recorded in their ledgers).

The `audit0_dyn_conventions_and_inventory.py` script is DYN-0 of the dynamics program (`route_E/ROADMAP.md` section E).  Canonical implementations live in case-sensitive `route_E/code_dyn/`; root `code/audit*_dyn*.py` paths are compatibility delegates and must not acquire independent physics logic.  Ledgers remain in whitelisted `output/audit*/` lanes.  Prefer `python3 route_E/code_dyn/run_route_e_dynamics.py`; its claim registry, not `all_pass`, controls promotion.

The Route-E paper (`route_E/tex/route_e_first_principles.tex`) and the submission letter draft (`route_E/tex/three_families_killing_contact.tex`) build with the same pdflatex/bibtex cycle run inside `route_E/tex/`; their bibliography resolves to `../../paper/refs.bib`. The submission program (units, prior-art verdicts, checklists) is tracked in `route_E/SUBMISSION.md` and ROADMAP section G.

There is no test framework. Each verifier is a standalone script that writes a JSON + Markdown ledger pair into the sibling `output/` directory and prints a one-line summary. The JSON is the machine-readable claim record and includes explicit negative-boundary flags (e.g. `global_e3_divisor_constructed: false`). New checks should follow this pattern: script in `code/`, paired `.json`/`.md` ledger in `output/`, boundary flags stating what is NOT claimed.

## Architecture: the route system

The project is organized as parallel "routes" with strict claim boundaries:

- **Route A (theorem core, main paper):** three representation-geometric results only — (1) the Spin(10) half-spinor 16 restricted to the SM face gives exactly Q, L, u^c, d^c, ν^c, e^c ("six faces"); (2) three families protected by dim H⁰(CP¹, O(2)) = 3; (3) the Majorana contact direction is the unique SL(2)-invariant second transvectant K_tr. Everything beyond these is a conditional input, an optional extension, or a deferred audit.
- **Route B (optional extension, in the paper):** a hidden sterile messenger X generates the contact via a Schur complement, λ² = ζ. Explicitly not part of the theorem core.
- **Route C (amplitude bootstrap):** exists only in git history; currently inactive.
- **Route D (`route_d/`, sandbox):** string-constrained interpretations (D1: local F-theory D5→E6 placement; D2: E3-instanton reading of ζ K_tr). Promotion into the paper requires: a precise statement, an assumptions-vs-derived list, a reproducible script, and a non-overclaiming impact note (see `route_d/README.md`). Until promoted, Route-D material enters the paper only as the optional appendix and must never be cited as evidence.
- **Route E (`route_E/`, sandbox):** a first-principles reconstruction of the Route-A skeleton from P0 consistency, P1 face projection, P2 self-carried index, and P3 canonical contact.  The enumerated Spin(10)-16 result, forced ν^c, genus ladder, and bound `N_fam <= 3` are theorem-level.  The stronger `g=0/N=3`, two-center, and Killing-contact claims are conditional on an added semisimplicity/non-degenerate-Killing axiom: a one-dimensional abelian algebra is a counterexample to the former Cartan inference because it admits non-degenerate invariant forms up to scale.  The two-model boundary theorem still proves ζ underivable.  Dynamics promotion is separately fail-closed by `route_E/code_dyn/dyn_claim_registry.json`; RE-SC3/4/5 are present but unpromoted, DYN-5 is invalid pending rederivation, DYN-9b-2 is preliminary, and DYN-7/9b-3 require flavored thermal inputs.
- **Companion audits (deferred):** full flavor fit, d=5 proton-decay Wilson tensors, thresholds, source-sector UV origin. The dependency graph is fixed by Audit 0 as `0 → (1 ∥ 4a) → 3 → 2`, with 4b as an optional Route-D/global-geometry track.

## Claim discipline (the load-bearing convention)

The paper's central result is deliberately bounded: "face-projection motivation + conditional Spin(10) EFT ≠ unconditional PSLT-only GUT proof" (boxed corollary). Every statement carries a status tag — proved/constructed, conditional input, deferred audit, optional UV interpretation, or not claimed (see the Lean-Note Status Ledger section of the tex). When editing the paper or adding results, never move a conditional input into the theorem core without a verified derivation, and keep the "Not claimed" list intact.

## Fixed conventions

- Benchmark contact coefficient: ζ = 0.1076472949 + 0.0736514853i. These digits are reproducibility anchors shared bit-identically between the tex tables and the verifier scripts; only the first few significant figures are physically meaningful.
- K_tr normalization: K_tr² = I/3, K_tr⁻¹ = 3 K_tr, in the normalized (ψ₀, ψ₁, ψ₂) family basis; K_tr = (1/√3)[[0,0,−1],[0,1,0],[−1,0,0]].
- Hypercharge: Y = T_{3R} + (B−L)/2, GUT normalization k_Y = 5/3.
- CMSGUT vacuum convention (Audit 4a.1): Aulakh–Girdhar x = −λω/m, ξ = λM/(ηm), cubic 8x³ − 15x² + 14x − 3 = −ξ(1−x)², cross-checked against the ABMSV convention after the ω-sign change.

## Git state caveats

- The working tree intentionally shows many deleted-but-uncommitted files (`code/audit_*.py`, `code/verify_*.py`, `REPRODUCIBILITY.md`, `roadmap.md`, `output/`): heavy audit history is preserved in git history while the tree is kept lean. Do not restore, commit, or "clean up" these deletions without explicit direction.
- The paper's "Minimal Reproducibility Manifest" still references `code/verify_d5_half_spin_hypercube.py`, `code/verify_hidden_zeta_origin.py`, `output/...`, and `roadmap.md`; those paths currently exist only in git history (`git show HEAD:<path>`).
- `.gitignore` whitelists specific `output/` subdirectories (audit05, audit1, audit1b, audit4a, audit4a1, d5_half_spin, hidden_zeta_origin). Keep new outputs inside whitelisted directories or extend the whitelist deliberately.
