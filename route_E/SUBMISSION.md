# Submission program notes (SUB lane)

Working document for the submission program (ROADMAP section G).
Created 2026-07-06.

## SUB-A: prior-art due diligence, round 2 (2026-07-06)

Adversarial search verdicts.  Every candidate-novelty claim gets:
**CLEAR** (no prior hit found), **ADJACENT** (related mechanisms exist,
ours differs -- cite and distinguish), **CONFLICT** (prior literature
overlaps or contradicts -- corrective action required before any
submission).

### (i) The three-point assembly (SUB-1 content)

| claim | verdict | notes |
|---|---|---|
| Family symmetry GENERATED as the automorphism algebra of its own carrier (not imposed) | **ADJACENT** | Octonion/sedenion/triality routes to "why three" exist (S3 from Aut(S), Spin(8) triality, exceptional Jordan algebra) but use a DIFFERENT mechanism (discrete outer automorphisms of a division-algebra structure, not the automorphism Lie algebra of the family carrier with the Cartan criterion).  SU(2)-flavor models exist but the symmetry is imposed with flavons.  SUB-1 must cite and distinguish in one paragraph. |
| N_fam <= 3 a priori with equality by the Cartan criterion (g = 0 via dim H^0(CP^1, O(2)) = 3) | **ADJACENT** | The paradigm "generation number = compactification index" is classical (CHSW, heterotic N_gen = index/c_3 -- already cited).  The specific self-carried H^0(CP^1, O(2)) construction with the a-priori bound was not found. |
| Majorana contact direction = Killing form of the family algebra | **CLEAR** | No hit for Killing-form Majorana textures; texture literature (zeros, cyclic, etc.) is unrelated. |
| Second-transvectant reading of the contact | **CLEAR** | "transvectant" does not appear in the flavor-physics literature at all. |

Assembly verdict: **the combination remains unclaimed in the searched
literature**; SUB-1 is viable with an ADJACENT-work paragraph
(octonionic/triality family symmetries; index-theorem generation
counting; SU(2) flavon models).

### (ii) The non-SUSY 210 vacuum (SUB-2 content) -- CONFLICT FOUND

- **He & Meljanac, PRD 33, 2695 (1986)** ("Symmetry breaking of SO(10)
  by the 210-dimensional Higgs boson and Michel's conjecture"):
  minimization of the 210 potential; finds finite parameter ranges
  whose minima break to SU(2)xSU(2)xSU(4) **or
  SU(2)xSU(2)xSU(3)xU(1)** (the left-right group).
- **Meljanac & O'Raifeartaigh-type stability analyses, PRD 39, 3110
  and PRD 40, 2098 (1989)**; also Nucl. Phys. B (1986): stability of
  the 210 breaking; the left-right scenario "can satisfy stability
  requirements" in a class of potentials.
- Modern context: Jarkovska-Malinsky-Susic
  ([2304.14227](https://arxiv.org/abs/2304.14227)) establish tree-level
  tachyon trouble and the one-loop rescue for the 45+126 minimal model
  (consistent with our landmine framing); Ohlsson et al.
  ([2212.11315](https://arxiv.org/abs/2212.11315)) treat 210/45 + 126
  chains at two loops (running-level, not potential-level).

CONSEQUENCES:
1. The DYN-9b-1b claim "the 210-only left-right vacuum is NEVER a
   tree-level minimum" **conflicts at abstract level** with the
   1986-89 results.  The suspect is exactly our disclosed hole: the
   left-right vev RATIO was FIXED at a/p = 0.8
   (`LR_vev_ratio_fixed` flag), while the classic analyses minimize
   over the full ratio.  The PS-side conclusions are NOT affected
   (single-vev pattern, no ratio freedom): K4 and its triple-test
   stand.
2. **Corrective audit SUB-B0 is mandatory before any submission**:
   free the left-right vev ratio in the 9b-1b/1c machinery, re-scan
   tree positivity over (ratio x couplings x epsilon), and either
   reproduce the He-Meljanac stable left-right regions or localize the
   disagreement.  Until then, K5r's "adjoint route required" narrative
   is DOWNGRADED to "ratio-conditional, under correction" in the next
   ledger refresh.
3. SUB-2 in its planned form ("tree-level verdicts" as headline) is
   NOT submittable.  It survives only if SUB-B0 turns it into a
   reproduce-and-extend paper.  One potentially NEW technical point
   emerged: the quartic-invariant literature counts **4** independent
   210 quartics (e.g. the ||(phi phi)_45||^2, ||..._210||^2,
   ||..._1050||^2, ||phi||^4 basis), which matches our machine count of
   the parity-EVEN invariants -- but our Hodge-epsilon contraction is a
   machine-independent FIFTH invariant (SO(10)-only, parity-odd,
   vanishing identically on the singlet slice, contributing only to
   off-slice masses, i.e. to STABILITY).  If the 1986-89 "most
   general" potentials omitted it, their stability conclusions were
   derived without one invariant -- to be verified against the full
   texts (paywalled; needs library access or inspire scan) before any
   claim is worded.

### (iii) Boundary-theorem genre

Not exhaustively searched (low priority-risk: methodological claim,
not a priority-sensitive result).  No hit encountered incidentally.
SUB-3's framing should present it as unusual practice without claiming
absolute novelty.

### Search log

Queries run 2026-07-06 (WebSearch): automorphism family symmetry +
Killing form; Killing-form Majorana texture; H^0(CP^1,O(2)) three
generations; He-Meljanac 210; non-SUSY 210 tree vacuum; 210 quartic
invariant count; transvectant in flavor physics; Jarkovska-Malinsky
210/126.  Key sources:
- https://link.aps.org/doi/10.1103/PhysRevD.33.2695
- https://journals.aps.org/prd/abstract/10.1103/PhysRevD.39.3110
- https://journals.aps.org/prd/abstract/10.1103/PhysRevD.40.2098
- https://arxiv.org/abs/2304.14227
- https://arxiv.org/abs/2212.11315
- https://arxiv.org/pdf/gr-qc/9512033 (the 4-invariant 210 potential)

## SUB-B0 outcome (2026-07-06): conflict RESOLVED as our artifact

`code/audit9_dyn9b1d_lr_ratio_scan.py` (7/7): freeing the left-right
vev ratio reproduces stable tree-level left-right minima at 8/8800
sampled points, all at ratios 0.6-0.7 -- consistent with He-Meljanac
(1986) at claim level; the 9b-1b negative at fixed a/p = 0.8 was a
ratio artifact.  The epsilon question is answered NO: the fifth
invariant neither rescues nor kills any sampled point.  K5r reverted
in the ledger and both papers.

Consequences for the submission units:
- SUB-2 is now honestly a **reproduce-and-extend** paper: reproduces
  the classic left-right/Pati-Salam tree structure with modern
  machine-verified machinery, extends it with (i) the machine
  invariant-basis analysis (Q5 = Q2 identity; the parity-odd fifth
  invariant, shown stability-IRRELEVANT in the sampled space -- an
  honest negative), (ii) exact Goldstone-gated Hessians, (iii) the
  126bar mixed-sector spectrum, (iv) the threshold/proton pipeline
  coupling the vacuum to Hyper-K phenomenology.  Viable but weaker
  than hoped; full-text check of He-Meljanac remains advisable for
  wording.
- SUB-3 is UNBLOCKED (the K5r correction has propagated).

## SUB-B outcome (2026-07-06): the theta_23 attack is defused

`code/audit1_dyn4c_kernel_dirac_refit.py` (7/7): the 6.9 sigma pull is
a frozen-anchor property -- the Majorana sector refits the NuFIT 6.0
centrals IDENTICALLY at unchanged Dirac kernels (chi^2_osc = 0 by
construction), a 5% charged-lepton kernel perturbation absorbs the
shift alone, and contact essentiality survives 1000 kernel-perturbed
refits up to 30% (cf never below 0.01).  Referee-facing sentence: the
benchmark card is a reproducibility anchor; the framework fits current
oscillation data exactly, and its one falsifiable contact statement is
robust to both data refreshes and kernel deformations.

## SUB-1 status (2026-07-06): draft complete

`route_E/tex/three_families_killing_contact.tex` -- "Three fermion
families from the automorphism algebra of their own carrier, and the
Killing-form Majorana contact".  Builds clean (0 errors, 14 references
resolved; the bibtex control-entry warning is the benign shared-refs
duplication).  8 pp preprint ~ 5 pp two-column.

Checklist before submission (author-side):
- [ ] read-through and voice pass;
- [ ] venue: PRD (letter-style) vs PLB (needs elsarticle reformat);
      recommendation: arXiv (hep-ph, cross-list hep-th/math-ph) first,
      then PRD;
- [ ] cover letter: lead with the assembly novelty (SUB-A verdicts:
      Killing contact CLEAR, transvectant CLEAR, self-carried family
      ADJACENT-with-distinction), the machine-verification
      reproducibility, and the explicit not-claimed discipline;
- [ ] decide whether to cite the companion reconstruction note as an
      arXiv preprint (requires SUB-3 posted first) or as a repository
      reference (current draft: repository only).

## Revised submission order

SUB-A (done) -> SUB-B0 (done; conflict resolved) -> SUB-B (done;
theta_23 defused) -> SUB-1 (draft done; author review pending) ->
SUB-3 (arXiv flagship -- posting it FIRST also gives SUB-1 its
companion citation) -> SUB-2 (reproduce-and-extend form, optional).
