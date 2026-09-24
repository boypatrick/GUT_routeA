# Route G sources and interpretation ledger

Reviewed: 2026-09-24.
Current mainline: relational quantum clocks; old KK/scattering sources retained below.

## G-R4 physical-source and cross-transition audit

Selected source: one isotropic thermal electromagnetic bath, coupled
through -d.E to the E2 and E3 transitions of the same 171Yb+ ion.
This is a conventional-spacetime control, not a new definition of Chen Ed.
All numerical inputs and source PDF hashes are pinned in
[BBR_CLOCK_INPUTS.json](data/BBR_CLOCK_INPUTS.json).

- T. Schneider, E. Peik and Chr. Tamm, *Sub-Hertz Optical Frequency
  Comparisons between Two Trapped 171Yb+ Ions*, Phys.Rev.Lett.94,230801(2005),
  [DOI](https://doi.org/10.1103/PhysRevLett.94.230801),
  [institution-hosted primary PDF](https://www.ptb.de/cms/fileadmin/internet/fachabteilungen/abteilung_4/4.4_zeit_und_frequenz/pdf/PRL_94_230801.pdf).
  Independently applied-field E2 scalar polarizability difference:
  alpha_ground-alpha_excited=-6.9(1.4)e-40 J m^2/V^2.
  **Sign converted** to excited-minus-ground +6.9(1.4)e-40 here.
  Original formula/convention inspected, including the rendered page.
- N. Huntemann et al., *Single-Ion Atomic Clock with 3e-18 Systematic
  Uncertainty*, Phys.Rev.Lett.116,063001(2016),
  [DOI](https://doi.org/10.1103/PhysRevLett.116.063001),
  [institution-hosted primary PDF](https://www.ptb.de/cms/fileadmin/internet/fachabteilungen/abteilung_4/4.4_zeit_und_frequenz/pdf/2016_Huntemann_PhysRevLett.116.063001.pdf).
  Independent near-IR Stark-response measurements constrain the E3
  static excited-minus-ground scalar difference +0.888(0.016)e-40.
  The published eta(300 K)=-0.0015(7) is kept separately from our
  static leading table; no full E2 dynamic response is invented.
- A. Tofful et al., *171Yb+ optical clock with 2.2e-18 systematic
  uncertainty and absolute frequency measurements*, Metrologia61,045001(2024),
  [DOI](https://doi.org/10.1088/1681-7575/ad53cd),
  [institutional primary PDF](https://eprintspublications.npl.co.uk/10218/1/eid10218.pdf).
  Section3.6 independently checks the BBR sign/formula and
  E_rms=831.9(T/300 K)^2 V/m convention; relevant rendered page checked.
  Its apparatus radiometry and uncertainty are NOT imported as our
  hardware performance. Corrected frequency outputs are NOT treated
  as raw temperature-cycle observations.
- M. Filzinger et al., *Improved limits on the coupling of ultralight
  bosonic dark matter to photons from optical atomic clock comparisons*
  (2023), [primary abstract](https://arxiv.org/abs/2301.03433).
  Supports the feasibility of interleaved E2/E3 interrogation of the
  same ion as a protocol, not a new blackbody modulation dataset.
- Rounded frequencies only: E2 688.359 THz is consistent with
  [Tamm et al., Phys.Rev.A89,023820(2014)](https://www.ptb.de/cms/fileadmin/internet/fachabteilungen/abteilung_4/4.4_zeit_und_frequenz/pdf/2014_Tamm_PhysRevA.89.023820.Cs-Based_optical_frequency_measurement_.pdf);
  E3 642.121 THz with
  [Lange et al., primary text](https://arxiv.org/html/2010.06620v2).
  We do not claim new recommended absolute frequencies or optimize their digits.

Input-selection caution: the current abstract of C. F. A. Baynham et al.,
*Measurement of differential polarizabilities at a mid-infrared wavelength
in 171Yb+*, [arXiv:1801.10134](https://arxiv.org/abs/1801.10134), explicitly
warns about its uncertainty characterization/applicability. Do not pool
its tighter-looking errors without resolving that warning. Our pinned
older independent inputs are a transparent structure test, not a claim
to the newest or most precise polarizability combination.

The source values fix signs; our Planck integral, diagonal-level lemma,
ratio formulas and stress box are separately derived. The +/-10%
dynamic envelope is an explicit conditional assumption, not sourced
atomic theory. No temperature-cycle data, fitted slowdown, autonomous
recording apparatus, emergent space or empirical universal-time effect
is claimed. In particular an additional exactly common multiplicative
factor is unidentifiable from this same-bath ratio alone.

## G-R3 density candidate and relative-clock audit

- A. R. H. Smith and M. Ahmadi, *Quantizing time: Interacting clocks
  and systems*, Quantum3,160(2019),
  [primary text](https://arxiv.org/html/1712.00081v3).
  Generic reference-clock interactions can produce a nonlocal
  conditional evolution kernel. G-R3 avoids claiming a general theorem:
  its reference A is decoupled from the explicitly interacting BD
  sector, whose finite stationary solution is derived directly.
- A. Blais, A. L. Grimsmo, S. M. Girvin and A. Wallraff,
  *Circuit quantum electrodynamics*, Rev.Mod.Phys.93,025005(2021),
  [publisher paper](https://doi.org/10.1103/RevModPhys.93.025005),
  [full text](https://harvest.aps.org/v2/journals/articles/10.1103/RevModPhys.93.025005/fulltext).
  SectionsIII.B and AppendixB.2.b establish the Jaynes-Cummings
  spectrum and dispersive shift chi=g²/Delta. The review also warns
  that an ideal two-level model is not generally a quantitative transmon
  model. Our note independently derives the commutator, exact blocks
  and error bound as a sign/structure audit, not a device forecast.

Our counting-measure Ed definition, rational finite reference spectrum,
shielding, matched state and model-unit inputs are declared choices.
They are not inferred from either source or established by Chen.
The exact mean-density counterexample and fixed-calibration comparison
are calculations in this model; the original universal slowdown claim
remains unproved. No physical density, clock frequency or coupling was
fitted. See [DENSITY_CLOCK_CARD.md](DENSITY_CLOCK_CARD.md).

## G-R2 two-event records: sources and added assumptions

- V. Giovannetti, S. Lloyd and L. Maccone, *Quantum Time* (2015),
  [primary text](https://arxiv.org/html/1504.04215v3), Measurements,
  Eqs.31–46. Explicit memories give sequential instrument statistics;
  later readout need not overwrite the first record. Our derivation
  uses this established measurement structure, not a newly discovered
  Born rule. We do not import an infinite canonical clock.
- D. Aharonov et al., *Adiabatic Quantum Computation is Equivalent to
  Standard Quantum Computation*, SIAM J.Comput.37,166–194(2007),
  [primary paper](https://arxiv.org/abs/quant-ph/0405098).
  This is the history-state propagation-constraint construction used
  as mathematical background. Our small open-path sum of squares and
  its one-dimensional kernel are derived explicitly in the note.
  No adiabatic implementation or locality theorem is claimed for our
  clock-and-apparatus model.

The six-configuration history, epoch bit, blank input, two fixed event
stamps and scheduled recording gates are our declared protocol choices.
They are not extracted from Chen, measured in a clock experiment, or
claimed as dimension-minimal. The exact three-level probability tables
and blank-to-filled cycle obstruction are direct algebraic diagnostics,
not empirical data. G-R2 replaces the noninteracting G-R1 constraint;
it does not supply apparatus energetics or a density-dependent rate.
See [the G-R2 derivation](tex/route_g_two_event_records.tex).

## G-R1 relational-time sources and interpretation boundary

- Chen's final five-page document below was reread in full, with all
  pages visually inspected during the time-concept audit. Sections1.2,
  1.8 and1.9 / Figures2,9,10 motivate examining information, event
  quantities and clocks. They do not specify the Hilbert tensor split,
  a constraint operator, clock POVM or the spectra used in G-R1.
- E. Moreva et al., *Time from quantum entanglement: an experimental
  illustration*, Phys.Rev.A89,052122(2014),
  [primary text](https://arxiv.org/html/1310.4691v1).
  Read the conditional-state construction, photon implementation and
  distinction between single-clock conditioning and two-time measurements.
  The laboratory experiment illustrates the mechanism; it does not prove
  Chen's Ed interpretation, derive spacetime or supply our model's data.
- V. Giovannetti, S. Lloyd and L. Maccone, *Quantum Time*,
  Phys.Rev.D92,045033(2015),
  [primary text](https://arxiv.org/html/1504.04215v3), sectionsI–II.
  Constraint projection, energy-origin freedom and the need to specify
  measurement records inform the audit. Our finite cyclic construction
  does not borrow the infinite ideal clock's exact canonical commutator.
  Those sequential-measurement results were not implemented in G-R1;
  G-R2 above now gives a separate bounded two-record construction.
- P. A. Höhn, A. R. H. Smith and M. P. E. Lock,
  *The Trinity of Relational Quantum Dynamics*, Phys.Rev.D104,066001(2021),
  [primary paper and abstract](https://arxiv.org/abs/1912.00033v3).
  Used for the distinction between relational conditional descriptions,
  covariant clock POVMs and invariant observables. G-R1's finite
  Fourier-twirled effects are derived directly; no general equivalence
  theorem for arbitrary interacting clocks is asserted.

The finite positive-energy complementary spectra, chosen state, clock
measurement, and coherent/dephased comparison are explicitly stated in
[the G-R1 TeX](tex/route_g_relational_clock.tex). They are mathematical
examples, not extracted experimental results or a novel invention of
Page-Wootters time. Ei is only tentatively compared with joint state
information plus an event algebra. G-R1 assigned no energy-density,
entropy-density or gravitational meaning to Ed. The separate G-R3
candidate above now supplies a counting-measure excitation density,
without equating it with the original claim.

The exact same-energy/different-history example excludes a density-only
identification in this model; it does not exclude every possible richer
definition of Ei/Ed. Common phase labels and calibrated readings are not
automatically spacetime time coordinates. See
[RELATIONAL_TIME_CARD.md](RELATIONAL_TIME_CARD.md) for the full dictionary.

## Original motivation, not established dynamics

BoYu Chen, *Energy spectrum for elementary particles*, final version,
five pages, available through the 2022 conference contribution:

[Exact final PDF](https://indico.cern.ch/event/1109513/contributions/4969967/attachments/2480973/4410560/Energy%20spectrum%20for%20elementary%20particles_final.pdf)

The downloaded file was initially read in full as text with Figures1,8,9
inspected; the later time-concept audit reread and visually inspected
all five pages, especially Figures2,9,10.
File length: 506899 bytes.
SHA-256: d86c71cee1d61c1af2ec33b3084ff4ed626e8d05551bd5685327f9e821f47e7e.
The PDF is not duplicated into the repository.

Relevant motivation:

- Pages1–2, section1.2 and Figure2: time associated with Ei accumulation
  and Ed, not an assumed external time coordinate.
- Page 1, section 1.1: internal phase/state dependence of mass appearance.
- Page 4, section 1.8 and Figure 8: elementary particles described through
  energy-information levels, of which only some are visible.
- Page 4, Figure 9 and accompanying text: internal deformation and mass.

G1 supplies a possible conventional field-theory realization of a limited
part of that intuition. The source does not derive the G1 action,
compactification radius, holonomy or probe couplings. Ei, Ep and En are
NOT silently identified with fields, positive/negative KK labels or
positive/negative Hamiltonian energies. Its time, gravity, entanglement,
negative-mass and superluminal claims are not established by G1.

## Established mechanisms used in the prototype

1. David Tong, *Lectures on String Theory*, section 8:
   [author's lecture text](https://arxiv.org/html/0908.0333v3#S8).
   Used for the compactification mechanism and the distinction between
   physical compact momenta and a mere change of description.
   G1 has no string oscillators or winding sectors and claims no string UV
   completion. Its mass formula is also derived directly from its action.

2. Atsuyuki Yamada, *The UV sensitivity of the Higgs potential in
   Gauge–Higgs Unification*, PTEP 2021, 093B01, section 2, especially
   equation (10):
   [original research article](https://doi.org/10.1093/ptep/ptab085).
   Used for the Wilson-line shift of KK momenta. G1 keeps the connection
   external and does not import that paper's gauge dynamics or Higgs
   potential. No all-orders finiteness or stabilization theorem is assumed.

The two antipodal source operators, their normalization and the resulting
two-channel residue sum rule are the explicitly declared G1 apparatus.
They are not quoted as a prediction of either source above.

## Relation to prior project work

Route-F P-SPEC1 treated common four-dimensional poles and interaction-defined
channels without physical extra dimensions. G1 instead postulates one
physical compact dimension and derives its tower. Similar linear-algebra
tests may be reused as principles, but G1 does not inherit a physical
Spin(10) embedding, fitted masses or a solved Route-F flavor problem.

## G2 local conversion: normalization sources and new assumptions

- Particle Data Group, *Kinematics*, Review of Particle Physics (2024),
  sections 49.5--49.5.1, equations 49.27--49.33:
  [official review](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-kinematics.pdf).
  Used to check invariant flux and two-body cross-section normalization.
- David Tong, *Quantum Field Theory*, sections 3.4 and 3.6:
  [author's lecture notes](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S3.html).
  Used to check contact-vertex and S-matrix conventions.

The extra neutral complex detector X, positive local quartic interaction,
its engineering mass/coupling, and finite detector packet are declared
G2 assumptions. The compact selection rule, heavy-recoil threshold and
inclusive-recoil calculation are derived from that action in the G2 TeX.
These sources do not establish that nature contains this detector or
that the original energy-spectrum document uniquely implies this model.
No novelty claim is made for KK scattering or standard two-body
kinematics. The purpose is a concrete physical mechanism with an explicit
energy/momentum ledger.

## G2-R joint readout (2026-09-23)

The phase-space normalization is inherited from G2 and checked against
the same PDG review. Joint rate normalization, mass-sign coarse graining
and kinematic sign reconstruction are derived for the unchanged G2
action and preparation. No source is cited as experimental evidence
for the predicted lines or as a design for measuring compact momentum.
Reconstruction uses the assumed mass law and conservation, so it cannot
serve as an independent verification of those assumptions.

## G2-S/D readout and activated drive (2026-09-23)

- M. Moskalets and M. Büttiker, *Floquet scattering theory of quantum
  pumps*, Physical Review B 66, 205320 (2002), section II, equations
  (1) and (6):
  [primary manuscript](https://arxiv.org/abs/cond-mat/0208356),
  [journal record](https://doi.org/10.1103/PhysRevB.66.205320).
  Used for the established energy-sideband and incoherent incoming-flow
  framework. Its mesoscopic model is not imported as our relativistic
  action, nor as evidence for an extra dimension.
- The same PDG kinematics review above supplies state/phase-space
  normalization. The driven two-body kernel is derived directly with
  initial energies retained in the flux and pump-shifted final energy.

The user explicitly activated the time-dependent option on 2026-09-23.
The selected local modulation g(t)=g[1+epsilon cos(Omega t)], epsilon=.5,
Omega=.2, is a new declared engineering assumption; it is not uniquely
implied by the energy-spectrum paper. The truncated-normal incoherent
beam, Gaussian readout, perfect mass labels, equal acceptance and 5%
Bayes-error decision threshold are likewise stated choices, not
experimentally calibrated facts. The finite-width and Floquet rates are
conditional, long-time/tree-level calculations. No quantum-pump model,
finite-collision absolute probability or all-orders UV claim is inherited
from a citation. No novelty claim is made for standard Floquet sidebands
or Bayesian classification.

## G2-M readout contract (2026-09-23)

Particle Data Group, *Particle Detectors at Accelerators* (2024),
section 35.13, equations (35.59)--(35.63):
[official review](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-particle-detectors-accel.pdf).
Used only to distinguish a real momentum-measurement design, with
charge, field, geometry and tracking/material errors, from an abstract
Gaussian response. We import no quoted instrument resolution or
calibration into this arbitrary-unit scalar model. Its external U(1)
is not silently identified with electromagnetism.

The affine calibration box, frozen-classifier analysis and symmetric
ternary tag response are declared engineering choices. The continuous
CDF bound, affine oracle identity, channel-degradation proof and
same-data no-information identity are derived in the new TeX.
No source establishes the existence of an instrument for X or an
event-resolved work counter for the prescribed classical pump.

## G2-P physical anchors and candidate visible probe (2026-09-23)

- J. M. Cline, K. Kainulainen, P. Scott and C. Weniger,
  *Update on scalar singlet dark matter*, Phys. Rev. D 88, 055025 (2013),
  [primary manuscript](https://arxiv.org/abs/1306.4710),
  [full text](https://arxiv.org/html/1306.4710v3), equations (3), (23),
  (26). Used to cross-check Higgs-density coupling and nucleon-scattering
  normalization. Our complex-scalar decay has distinguishable final
  particles, unlike the identical real-scalar convention in their (3).
  No old exclusion region or dark-matter abundance is imported.
- ATLAS Collaboration, *Combination of searches for invisible decays of
  the Higgs boson using 139 fb^-1*, Phys. Lett. B 842, 137963 (2023),
  [primary experimental paper](https://arxiv.org/abs/2301.10731).
  An example of a physical constraint on a proposed Higgs coupling, not
  Route-G calibration. We do not claim it is the newest limit or apply
  its single-species interpretation directly to a full KK tower.

The scale-identifiability proof, signed-mode anchor reconstruction,
common-shift closure, brane/uniform comparison and tower-width/readout
compatibility formulas are derived in our TeX for the explicitly stated
candidate. The Higgs extension is only a 4D EFT option, not established
by the motivating paper or a completed 5D embedding. All new numerical
cards are synthetic. No real energy scale, nucleon form-factor estimate,
measured response or apparatus performance is assigned by these sources.

## G2-X physical feasibility with published response (2026-09-23)

These inputs now supersede the deliberately unfilled G2-P intake for
one conditional candidate; they do not measure the Route-G radius.

- Particle Data Group, *Electroweak Model and Constraints on New
  Physics*, §10.4.2, [currently served review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf).
  The served file carries the April 2026 footer/revised-November-2025
  text despite its 2025 URL. We use mh=125.08 GeV, v≈246 GeV and the
  **SM-predicted**, not measured, Higgs width 4.10 MeV. Constants are
  checked against the [PDG physical constants review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-phys-constants.pdf).
- ATLAS Collaboration, [2023 invisible-Higgs combination](https://arxiv.org/abs/2301.10731),
  cited above. Observed B_inv<0.107 at 95% CL is applied to the full
  open signed complex-X tower, conditional on SM production/visible
  width and invisible escaping final states in the drive-off vacuum.
  This is not a newest-world-limit or complete-collider-exclusion claim.
- XENON Collaboration, *Light Dark Matter Search with Ionization
  Signals in XENON1T*, [2019 primary paper](https://arxiv.org/abs/1907.11485v2),
  and its [official S2-only response release](https://github.com/XENON1T/s2only_data_release).
  We directly fold the NR matrix at pinned commit
  `5a364bc8709f2561e5a013ddea6993a5a7c8e313`. The full-volume selection
  response and 356770 kg-day exposure are used once; no halo exclusion,
  observed-event fit, background-free confidence limit or invented
  efficiency is imported. The actual local files, hashes, license
  metadata and bounded use are in [data provenance](data/xenon1t_s2only/PROVENANCE.md).
- M. Hoferichter, P. Klos, J. Menéndez and A. Schwenk, *Improved limits
  for Higgs-portal dark matter from LHC searches*, PRL119,181803 (2017),
  [primary manuscript](https://arxiv.org/abs/1708.02245): fN=.308±.018.
  Its small effective two-nucleon contribution is not counted twice.
  This nuisance is not a full finite-q nuclear/detector uncertainty.
- J. D. Lewin and P. F. Smith, Astropart. Phys.6,87 (1996),
  [primary paper](https://hepwww.pp.rl.ac.uk/groups/ukdmc/pub/papers/journal/app6-87.pdf),
  equations(4.7),(4.10–11): explicit Helm conventions. We retain only
  low-transfer elastic recoils within the response window, not an
  uncontrolled coherent extrapolation to MeV/GeV recoil endpoints.
- CIAAW, [xenon isotope composition and atomic weight](https://ciaaw.org/xenon.htm):
  natural-isotope number fractions and molar mass131.293 g/mol.
  Individual nuclear masses are explicitly approximated by A×u.

The selected Lambda=1–30 GeV interval spans the known Higgs threshold;
its limits and lambda_X scan cap are declared coverage choices, not
external measurements of the new sector. XENON's calibration-informed
response is real, but it does not supply the hypothetical source flux,
mode selection, pump, transmission or momentum measurement. The exact
factorization/no-sign-information result is derived in the G2-X TeX
for the stated scalar-current model and retained true-energy component;
it is not attributed to an experimental collaboration.

## G2-H physical source budget (2026-09-23)

- ATLAS Collaboration, JHEP05(2023)028,
  [arXiv:2207.08615v2](https://arxiv.org/html/2207.08615v2), sections1,2,5.
  The139 fb^-1 luminosity is actual;55.6±2.5 pb is its quoted SM
  Higgs production prediction. The agreeing visible-channel measurement
  assumes SM branching fractions and is not independently imported as
  a dark-source normalization after altering them.
- The [ATLAS invisible combination](https://arxiv.org/abs/2301.10731v2)
  supplies the same conditional allowance, not an observed X population.
- The [XENON1T paper](https://arxiv.org/html/1907.11485v2) supplies rounded
  active-cylinder geometry for a deliberately maximal single-pass column.
  It does not propose an LHC installation or a dark-particle beam.

The [external source card](data/HIGGS_SOURCE_CARD.md) separates inputs,
theory normalization, stress choices and missing transport. New bounds,
the linked-coupling yield and the source-specific sign-information theorem
are our derivations within the stated tree model. No source claims the
predicted dark particle exists, validates the original GeV drive, or
establishes sensitivity to a KK tower. No newest-limit claim is made.

## G2-T timing, momentum and signed-mode audit (2026-09-23)

- Particle Data Group, [Kinematics, 2025 update](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf),
  sections49.1 and49.4: the mass shell, beta=p/E, proper time and
  two-body relations. TOF inversion, finite-error boxes, the general
  source/response gate and the directional-recoil inversion are derived
  in our note for the stated assumptions, not claimed as new relativity.
- Particle Data Group, [Particle Detectors at Accelerators, 2025 update](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-particle-detectors-accel.pdf),
  section35.12, Eq.35.66: charged-particle magnetic-curvature momentum
  measurement. Its quoted performance is NOT a momentum sensor for
  neutral X. No real timing, momentum or direction calibration is used.
- [XENON1T S2-only primary paper](https://arxiv.org/html/1907.11485v2),
  detector response and drift-time discussion. The selected response is
  projected onto S2; it is not a tagged source-to-X flight time or a
  nuclear recoil direction. The pinned raw response is unchanged.

The [G2-T input card](data/TOF_READOUT_CARD.md) separates source facts
from requirements and synthetic examples. The source-sign inversion
is explicitly a revalidation of G2-R, not a new source implementation.
The even/odd interference criterion is elementary amplitude algebra,
not evidence for an activated new interaction. No hardware feasibility,
novel fundamental law or empirical KK-mode identification is claimed.
