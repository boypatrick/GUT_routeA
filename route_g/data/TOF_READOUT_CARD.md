# G2-T measurement input card

Date: 2026-09-23. Requirements and algebraic examples, **not calibration data**.
No action, coupling, source strength or raw detector response was changed.

| Input | Current status | Permitted use |
|---|---|---|
| X mass tower and neutrality | Retained conditional Route-G action | Mass-shell and identifiability derivation |
| Production timestamp t0 and covariance with momentum | Not supplied | Keep in general joint source measure; do not assume zero jitter |
| Flight path and incoming direction | Not measured for X | Symbols in likelihood; illustrative L=10 m only |
| Independent single-X momentum or energy | No specified physical sensor | Do not import charged-track resolution or infer p from assumed m |
| Current S2-only response | Actual energy-to-S2 matrix, no recoil vector or tagged X TOF | Retain G2-X bounded use, not a new timing/direction calibration |
| Directional recoil energy and vector | Not supplied | Same-action elastic inversion is only a mathematical alternative |
| Target isotope per event | Unobserved in natural Xe | A=132 is an algebraic example; data would need an isotope-mixture likelihood |
| Phi source/incoming and outgoing mode labels | No physical preparation/tag | Revalidate old conditional G2-R inference only |
| Sign-sensitive extra coupling/field | Not activated | Interference criterion and proposal gate only |

Numerical illustration:

- Lambda=5 GeV, E=mh/2=62.54 GeV, Higgs at rest, L=10 m.
  No assertion that LHC Higgs bosons are all at rest.
- Error box: |delta p|/p<=0.1%, |delta T|<=20 ps, |delta L|<=1 mm.
  These are declared bounded requirements, not measured standard deviations.
- Directional case: target M=132u, u=.93149410242 GeV, ER=10 keV.
  Target initially at rest and elastic; angle-only error budget holds beta
  and energy exact. Motion, binding, isotope mixture and joint response
  are not certified at the resulting sub-arcsecond requirement.
- The four-atom time/momentum source in the verifier is synthetic and
  explicitly correlated. It demonstrates selection reweighting, not flux.

Primary references checked on 2026-09-23:

- [PDG Kinematics, 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf),
  sections49.1 and49.4: mass shell, beta=p/E, proper time and decay kinematics.
- [PDG Particle Detectors at Accelerators, 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-particle-detectors-accel.pdf),
  section35.12, Eq.35.66: curvature momentum measurement requires charge.
  No published timing or charged-track performance is assigned to X.
- [XENON1T S2-only primary paper](https://arxiv.org/html/1907.11485v2),
  detector response and drift-time discussion: the retained response is
  projected onto S2, not a directional nuclear recoil or tagged X-flight
  measurement. Raw data provenance remains in xenon1t_s2only/PROVENANCE.md.

The source/readout yield rejection stays in force: any extra time,
momentum or direction gate can only reduce the G2-H3.9596863e-6 ceiling.
An alternative actual source, target or interaction requires its own
linked production/detection budget; this card does not invent one.
