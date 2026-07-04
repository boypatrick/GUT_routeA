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
python3 route_E/code/verify_route_e_first_principles.py
python3 route_E/code/verify_route_e_dependency_closure.py
python3 code/audit0_dyn_conventions_and_inventory.py
python3 code/audit4a_dyn1a_vacuum_goldstone_audit.py
python3 code/audit4a_dyn1b_full_spectrum.py
python3 code/audit3_dyn2_thresholds_unification.py
python3 code/audit1_dyn4a_seesaw_replay_zeta_posterior.py
python3 code/audit1_dyn4b_refreshed_card_unconditional_zeta.py
python3 code/audit7_dyn7_leptogenesis_argzeta.py
```

The DYN-4a script requires the archival output set restored under `route_E/output/` (user-managed; sha256 recorded in its ledger).

The last script is DYN-0 of the dynamics program (`route_E/ROADMAP.md` section E): new dynamics-lane scripts live at root `code/` with ledgers in whitelisted `output/audit*/` lanes; history assets are recovered read-only via `git show HEAD:<path>`, never restored into the tree.

The Route-E paper (`route_E/tex/route_e_first_principles.tex`) builds with the same pdflatex/bibtex cycle run inside `route_E/tex/`; its bibliography resolves to `../../paper/refs.bib`.

There is no test framework. Each verifier is a standalone script that writes a JSON + Markdown ledger pair into the sibling `output/` directory and prints a one-line summary. The JSON is the machine-readable claim record and includes explicit negative-boundary flags (e.g. `global_e3_divisor_constructed: false`). New checks should follow this pattern: script in `code/`, paired `.json`/`.md` ledger in `output/`, boundary flags stating what is NOT claimed.

## Architecture: the route system

The project is organized as parallel "routes" with strict claim boundaries:

- **Route A (theorem core, main paper):** three representation-geometric results only — (1) the Spin(10) half-spinor 16 restricted to the SM face gives exactly Q, L, u^c, d^c, ν^c, e^c ("six faces"); (2) three families protected by dim H⁰(CP¹, O(2)) = 3; (3) the Majorana contact direction is the unique SL(2)-invariant second transvectant K_tr. Everything beyond these is a conditional input, an optional extension, or a deferred audit.
- **Route B (optional extension, in the paper):** a hidden sterile messenger X generates the contact via a Schur complement, λ² = ζ. Explicitly not part of the theorem core.
- **Route C (amplitude bootstrap):** exists only in git history; currently inactive.
- **Route D (`route_d/`, sandbox):** string-constrained interpretations (D1: local F-theory D5→E6 placement; D2: E3-instanton reading of ζ K_tr). Promotion into the paper requires: a precise statement, an assumptions-vs-derived list, a reproducible script, and a non-overclaiming impact note (see `route_d/README.md`). Until promoted, Route-D material enters the paper only as the optional appendix and must never be cited as evidence.
- **Route E (`route_E/`, sandbox):** a first-principles reconstruction of the Route-A skeleton from a four-principle set (P0 consistency, P1 face projection, P2 self-carried index, P3 canonical contact). It upgrades to theorems: enumerated Spin(10)-16 uniqueness with ν^c forced, N_fam ≤ 3 with g=0/N=3 selected by the Cartan criterion, the two-center divisor as Poincaré–Hopf, K_tr = Killing form (B = 2√3 K_tr), Route-B form-uniqueness, and a two-model boundary theorem proving ζ underivable (witness family is the internal Route-B EFT at variable coupling — no Route-D dependency). A dependency-closure audit machine-verifies the former classical inputs (w0 involutions rank 4–15 + exceptional, exact cubic anomaly traces incl. the vanishing tensor of the Fock-built 16, exceptional minima); `route_E/ROADMAP.md` tracks all dependency items (R1–R15) and the paper's provenance ledger tags every ingredient as machine-verified / textbook-cited / retained input (RI-1..6) / benchmark corroboration. Same promotion discipline as Route D: the main paper's ledger is relocated only if Route E is accepted; until then it must not be cited as evidence.
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
