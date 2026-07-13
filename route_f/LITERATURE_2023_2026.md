# 2023--2026 Primary-Literature Matrix

Checked on 2026-07-13.  These papers are inputs for Route-F design, not proof
that the present model inherits their results.

| Source | Main result relevant to Route F | Required response |
|---|---|---|
| Jarkovska, Malinsky, Susic, [The trouble with the minimal renormalizable SO(10) GUT](https://arxiv.org/abs/2304.14227) (2023) | The minimal non-SUSY `45 + 126 + 10_C` model has a low-scale/proton problem and difficulty supporting a perturbative SM-like Higgs doublet. | Do not select this field content by minimality alone; require a full light-doublet and perturbativity gate. |
| Haba, Shimizu, Yamada, [Neutrino Mass in Non-Supersymmetric SO(10) GUT](https://arxiv.org/abs/2304.06263) (2023) | A `54 + 126 + 10_C` two-step branch performs gauge-unification and fermion/neutrino fits, including both mass orderings. | Use as the independent F-54 comparison branch and reproduce its scale/matching assumptions. |
| Djouadi, Fonseca, Ouyang, Raidal, [Non-supersymmetric SO(10) models with Gauge and Yukawa coupling unification](https://arxiv.org/abs/2212.11315), revised/published 2023 | Two-loop gauge/Yukawa running for Pati--Salam and left-right intermediate groups shows strong proton and threshold sensitivity. | Route F must use two-loop running and explicit matching; one-loop central values are triage only. |
| Fu et al., [Testing Realistic SO(10) SUSY GUTs with Proton Decay and Gravitational Waves](https://arxiv.org/abs/2308.05799) (2023) | Demonstrates the appropriate joint scale of analysis: two-loop RGEs, thresholds, flavor, leptogenesis, dark matter, proton decay, and gravitational waves. | Use as a benchmark for scope and uncertainty discipline; do not transfer its SUSY conclusions to the non-SUSY branch. |
| Kaladharan and Saad, [Fermion mass, Axion dark matter, and Leptogenesis in SO(10) GUT](https://arxiv.org/abs/2308.04497) (2023) | Shows how vectorlike fermion/scalar extensions can relieve minimal Yukawa-fit tension and simultaneously affect axion/leptogenesis predictions. | If the Route-F minimal flavor fit fails, pre-register this as an extension rather than adding ad hoc texture parameters. |
| Babu, Di Bari, Fong, Saad, [Leptogenesis in SO(10) with Minimal Yukawa sector](https://arxiv.org/abs/2409.03840) (2024) | A `10 + 120 + overline{126}` Yukawa sector fits masses/mixings, gives `N_2`-dominated leptogenesis, and checks unification/proton limits. | Use as the primary comparison for the F-210 flavor sector; reproduce full RG and density-matrix/flavored assumptions before claiming agreement. |
| Bertucci et al., [Positivity Bounds on Massive Vectors](https://arxiv.org/abs/2402.13327) (2024) | Massive-vector positivity requires nontrivial crossing functions, sum rules, and null constraints. | Replace Route-C constructed PSD residues with amplitude-extracted residues and crossing/dispersion tests. |
| Sato, [Majorana Fermion Zero Modes and Anomalies in String and M Theories](https://arxiv.org/abs/2404.18138) (2024 thesis) | Gives a conceptual nonperturbative zero-mode-parity/anomaly criterion for pointlike solitons and brane configurations. | Use only as an anomaly cross-check; it is not E3/M5 charged-zero-mode cohomology or a Pfaffian calculation. |
| Super-Kamiokande, [Search for proton decay via p to e+ eta and mu+ eta](https://arxiv.org/abs/2409.19633) (2024) | Updated channel-specific null search and lifetime limits illustrate that a theory needs a channel table, not one generic lifetime. | Refresh all experimental inputs and publish branching fractions and correlated theory errors by channel. |
| LEGEND Collaboration, [First Results on the Search for Lepton Number Violating Neutrinoless Double Beta Decay with LEGEND-200](https://arxiv.org/abs/2505.10440) (2025) | No signal; reports a `76Ge` half-life limit and nuclear-matrix-element-dependent effective-mass range. | Treat absence of a signal as a likelihood, not proof of Dirac neutrinos; propagate nuclear matrix element uncertainty. |

## Literature-driven design decisions

1. Keep both F-210 and F-54 until the same two-loop/threshold/proton pipeline
   compares them.
2. Use `10 + 120 + overline{126}` as a motivated flavor candidate, but count
   parameters and demand out-of-fit predictions.
3. Add massive-vector crossing/positivity to Route C only after full amplitudes
   are available.
4. Keep D1/D2 conditional: recent zero-mode consistency results strengthen the
   list of gates but do not construct the missing global geometry.
5. Refresh experimental likelihoods at execution time.  Numerical limits are
   versioned inputs, not permanent theorem constants.

The scoped recent search did not identify a 2023--2026 primary paper that
already supplies the exact global `Spin(10)` F-theory geometry, quantized flux,
massless hypercharge, exotic lifting, and E3/M5 Pfaffian needed here.  Route F
should therefore retain the relevant foundational literature rather than use a
recent but only adjacent zero-mode paper as a substitute.
