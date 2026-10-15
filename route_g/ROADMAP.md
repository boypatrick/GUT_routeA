# Route G Roadmap — Relational-time mainline and comparison models

Created: 2026-09-22
G1 verified and PDF visually reviewed: 2026-09-23
Current priority (2026-10-13): calibrate the actual controller cost while
keeping the G-R13 model and single pulse fixed. The model calculation is
complete and pushed as `5e1d8f9`; empirical calibration is OPEN pending a
specified hardware implementation and measurements. Its lower model work
does not establish lower electrical consumption. No pulse scan or autonomous
trigger study starts by default; there is no new mathematical prerequisite.
G-R4–G-R8 are retained physical comparison models. Their missing data and
apparatus calibration are not prerequisites for this event-network line.
No physical time, mass spectrum or universal density law is claimed.
Earlier G0–G2-T: preserved comparison/regression models, not the current mainline.
Status labels: open, in-progress, bounded-done, rejected-for-this-model.

## Authority and research discipline

The user explicitly chose Scheme A and requested a new Route G.
This is not the historical Route A and does not replace Route F.
Keep all input assumptions and claimed interpretations visible.
Do not import arbitrary particle-specific radii, projectors or mass shifts.
Do not turn every unfinished UV detail into another prerequisite paper.
No existing P54 parameters, fits or regression results are changed.

## Governing update — event-first candidate and a missing channel (2026-10-09)

The user approved returning to the original paper's events, relational
structure and compatible exchange, without initially assuming physical
space, seconds or a fitted density/slowdown function. G-R6–G-R8 assumed
a metric: their successful algebra does not derive that metric or the
source's Ei/Ed concept. Preserve, but demote, the apparatus mainline.

### Controller cost calibration requires a physical implementation

Decision date: 2026-10-13. **Selected next task: actual controller cost.
Status: open, awaiting hardware identification and measurement inputs.**
The G-R12/G-R13 cards specify no physical energy/time scale or controller.
The existing cold-atom citation supports the operator class, not our
21-state realization or its measured power. The stored xenon and optical-clock
datasets belong to different comparison models and cannot fill this gap.

To start a defensible calibration, supply one existing platform/controller
and any available raw records. Missing fields can then be requested for
that apparatus; do not invent them or select a new platform silently.

1. **Implementation and units:** identify the system and controller,
   the physical energy/frequency and time mapping, the number of controlled
   systems per attempt, and how the relay and compensation controls implement
   B=f B0 and Vc=f² Vc0. A pulse-generator model alone is not this mapping.
2. **Energy boundary:** provide synchronized, calibrated signed input-power
   records (or voltage/current records with measurement configuration and
   bandwidth) for both control channels and shared equipment. Include idle
   baselines, the complete attempt duration, and controller/storage state
   before and after. Count shared supply inputs once, not again as outputs.
3. **Return and reset:** identify whether returned energy reaches a reusable
   supply/storage device or is dissipated. Match initial/final stored energy
   or account for its change. Report pulse-only incremental energy separately
   from whole-cycle consumption including initialization and reset; do not
   identify G-R13 negative work with recoverable wall-supply energy.
4. **Matched performance:** retain the fixed G-R13 pulse and the G-R12
   reference, all attempts, measured preparation fidelity/residue and
   repeated-run uncertainty. The square reference is an ideal limit; any
   actual finite-rise waveform must be recorded and its implementation error
   reported, not silently treated as identical to that limit.

The supplied model criterion p0 < .008864075 + .016753645 eta remains
conditional, in model power units. Neither p0 nor eta is a measured fit,
and this two-parameter diagnostic is not a substitute for the physical
energy boundary above. No additional TeX derivation, numerical pulse scan
or hardware-success claim is needed before receiving the missing input.
If no apparatus exists, platform selection is a separate user decision;
a literature feasibility estimate would not be actual calibration.
The autonomous trigger with energy/reflection/return channels remains a
conditional alternative only if autonomous events become the chosen aim.

### G-R13 — One smooth pulse, compensated full-model energy ledger (2026-10-11)

Status: **bounded-done for the requested single-pulse comparison**.
G-R12 was committed and pushed first as `7d6d425` on
`codex/route-f-invariant-control`. G-R13 was subsequently committed/pushed
as `5e1d8f9` on 2026-10-13. No physical implementation is invented.

1. **One pulse, fixed model.** Keep gamma=.03, Delta=80 and g_peak=.1,
   all 21 states, old contact, relay edges and phases. Choose only
   f(s)=sin²(pi*s/T), B(s)=f B0, compensation=f² Vc0. Its compact extension
   is C1, with f=f'=0 at both boundaries; not C-infinity. Set T=50/3 from
   g integral(f²)=g tau=.625, tau=6.25. Same peak and effective area,
   but 8/3 longer. No pulse, amplitude, detuning or timing optimization.
   The continuous evolution parameter is supplied, not rounded to updates
   or promoted to emergent seconds.
2. **Why this can work.** The exact leading identity
   -f² B0†B0/Delta+f² Vc0=g f² C−sqrt(2)g f² I fixes the compensation.
   Two integrations by parts remove terminal boundary terms for the
   exact auxiliary amplitude. This explains suppression without claiming
   a uniform detuning theorem; actual dynamics remain the full 21 states.
3. **Simultaneous improvement — computed, not postselected.**
   F: .996311117 → .998120722; terminal P_aux: .003335852 → 5.02e-13;
   E_aux: .280211579 → 4.22e-11. The mediator is active during the pulse:
   sampled peak P_aux=.00184626. G-R11 still retains all failures and the
   exact 57/82 versus 25/82 unread-success backreaction.
4. **Full work, including compensation.** Integrate
   P_B=f' <V_B> and P_c=2 f f' <Vc0>. W_B=−.0381709643,
   W_c=+.0348565906; their sum −.00331437377 equals endpoint energy change.
   The square reference is +.263536868. Separate-port positive input
   decreases .263585530 → .171251419, so this is not just netting against
   returned energy. First-law equality also holds throughout the pulse.
5. **Tradeoff, not free control.** Returned work rises to .174565792;
   absolute separate-port exchange rises .263634192 → .345817211.
   Duration rises 6.25 → 16.6667. For the explicitly conditional cost
   C=Wplus−eta Wminus+p0 T, smooth wins only if
   p0 < .008864075 + .016753645 eta (model power, common 0≤eta≤1).
   Actual overhead, waveform losses, initial-state preparation and reset
   have not been calibrated. Negative net work is not zero apparatus cost.

70/70 new checks pass; **6,447/6,447 checks across 26 verifiers** pass.
G-R3/G-R4 reports gain only the two new verifier-preservation checks.
The combined nine-page PDF compiled successfully and was visually
reviewed on every page. Independent DOP853 tolerance refinement and a
unitary two-exponential method resolve the tiny terminal residue; the
1024/2048/4096-step refinements are NUMERICAL, not three physical pulses.
The coarse ODE tolerance was tightened after its initial error exceeded
1e-10; no physical parameter was changed. Parent G-R12 code/card/result
hashes stay unchanged. Derivation is appended as Sections 7–9 to the SAME
TeX, preserving the original six sections instead of opening another paper.

**Historical next choices; controller calibration selected on 2026-10-13:** calibrate a specific controller's supplied,
returned and overhead energy against the criterion above; or, if autonomous
events are the objective, specify a passing trigger with energy, reflection
and return channels. Either choice needs a physical input. Do not launch
more Gaussian/pulse scans, claim irreversible records, identify Ed or
promote any Route-F gate. The requested one-pulse test is finished.

Artifacts: [fixed pulse card](data/SMOOTH_SWITCHING_CARD.json),
[results and work ledger](output/gr13_smooth_switching.md),
[code](code/verify_gr13_smooth_switching.py),
[same TeX, G-R13 addendum](tex/route_g_correlated_preparation.tex),
[combined PDF](output/pdf/route_g_correlated_preparation.pdf).

### G-R12 — Local correlated preparation and finite mediator energy

Status: **bounded-done for local-operator existence, full-contact pulse
tests and finite-mediator accounting**. Autonomous control and a physical
hardware/energy calibration remain **open**. Preceding G-R9–G-R11 and
the retained G-R9b side branch were committed/pushed as `bcdbfb1` before
this stage began. G-R12 itself was subsequently committed/pushed as
`7d6d425` on 2026-10-11, before G-R13.

1. **Correct the symmetry scope — proved.** Joint translation acts as a
   scalar on the resonant triple. The density-contact family h(j−l) is
   circulant, but translation symmetry alone does not imply this form.
   The original rank-one 4/9 theorem is unchanged. G-R9b files are kept
   intact; this clarification supersedes any broader reading of its prose.
2. **A local added generator — constructed.** D=sum P_j tensor P_j,
   K=L3(0) tensor I+I tensor L3(0), C=i[K,D]/sqrt(2). It is Hermitian,
   swap and joint-translation invariant, and supported on 12 undirected
   nearest-hop configuration edges. On the resonant triple it is the
   required i(|11><A|−|A><11|), A=(|02>+|20>)/sqrt(2).
   K is a NEW hopping quadrature: the old resonant H0 commutator gives
   zero within this triple. Its phase/sign is not derived from symmetry.
3. **Original contact retained — tested.** Use inherited g=0.1 and
   alpha=0.25; fix n=25 before evaluating any gamma. Full nine-state
   fidelities at gamma=0.09/0.03/0.01 are 0.997141/0.999601/0.999875.
   All leakage and sudden on/off work are retained. This pulse does not
   commute with the old full H; there is no free violation of the 2/3
   energy-invariant preparation bound. Continued control loses the target.
4. **Mediator and Stark term — explicit.** One auxiliary configuration
   per allowed edge gives a 21-state model with E_aux=4+Delta. Eliminating
   it gives -B†B/Delta=gC−sqrt(2)g(I+D), with a bounded low-energy Schur
   remainder. Compensation +sqrt(2)gD is a separately declared term.
   The three detunings 20/40/80 are fixed diagnostics, not pulse fits.
5. **Small probability is not small energy — verified.** With compensation,
   gamma=0.03 and Delta=80, unconditional F=0.996311, P_aux=0.00333585,
   E_aux=0.280212 and net controller work=+0.263537. Without compensation
   the effective F is only 0.951831. A suddenly switched effective model's
   negative work cannot replace the full relay's positive work ledger.
   No auxiliary exit is discarded or normalized away.
6. **Existing readout reused — verified.** Swap symmetry leaves Dminus
   unoccupied before readout. The G-R11 instrument therefore succeeds with
   the unconditional preparation fidelity; auxiliary population counts
   as failure. Its exact unread-success 57/82 versus 25/82 mixture remains.
   Preparation, single-outlet transport and saved records stay distinct.

189 new checks pass; **6,375/6,375 checks across 25 verifiers** pass.
G-R9b and G-R11 code/card/result hashes remain unchanged. The six-page
standalone TeX compiled in the built-in editor; exported PDF was visually
reviewed on all pages. The cited cold-atom experiment supports the operator
class, not these particular relays, coefficients or autonomous controls.

**Historical next choice, now tested in G-R13 for the fixed mediator:** a single smooth on/off protocol for a
chosen mediator/control implementation, retaining compensation and high-energy
work; or a passing trigger with explicit energy, reflection and return.
Do not optimize hundreds of pulse shapes or add another prerequisite proof.
The local-generator existence task is closed. No density-only time law,
emergent seconds, irreversible arrow or Route-F gate is promoted.

Artifacts: [card](data/CORRELATED_PREPARATION_CARD.json),
[report](output/gr12_correlated_preparation.md),
[code](code/verify_gr12_correlated_preparation.py),
[TeX](tex/route_g_correlated_preparation.tex),
[PDF](output/pdf/route_g_correlated_preparation.pdf).

### G-R11 — One preparation and one relative-phase record

Status: **bounded-done for the finite energy-compatible instrument**;
**autonomous apparatus, isolation and arrow remain open**.

1. **Readout compression — proved.** In the unchanged two-dimensional
   dark subspace P, the 49 energy-invariant joint effects satisfy
   P Fbar(r,r+d) P=C_d/7. The common-origin r is uniformly redundant for
   every dark state, not merely Dplus. Summing r retains the reduced
   joint Lüders instrument. C_0=0; the six other C_d sum to P. This is
   a relative joint observable, not separate local clock measurements.
2. **Complete instrument and energy ledger — done.** K_failure=I−P and
   K_d=sqrt(C_d) form a complete instrument commuting with H. An explicit
   72×72 unitary writes a degenerate eight-level memory (distinct blank,
   failure, d=1,...,6). It conserves all energy spectral probabilities.
   From |1,1>, success is 2/3 and failure 1/3; their energies are 4 and
   4+gamma/3, averaging to 4+gamma/9. The initial input and pure blank
   are supplied. This collective coupling is NOT derived from the old
   vertex contact or G-R9b's oriented interaction.
3. **Backreaction and resources — calculated.** Unread successful output
   is (57/82)|Dplus><Dplus|+(25/82)|Dminus><Dminus|, analytically derived
   from a seven-root Chebyshev sum. The conditional orbit survives with
   a changed seed. A classicalized record has entropy 2.27003421247 bits;
   a reset bound is stated only for a thermal, degenerate, no-accessible-
   side-information, cyclic-controller setting. No write-energy or
   actual heat cost is inferred from this entropy alone.
4. **Reference boundary — proved.** Stationary ancilla + invariant pointer
   + conserving unitary imply invariant effects. Bare local phase effects
   fail this condition; the relative effects pass and have an explicit
   implementation. No external energy-phase asymmetry resource is needed
   for this abstract instrument. Internal calibration and pulse control
   are still supplied. Interaction-picture formulas do not remove the
   physical free-evolution factor or derive the pulse parameter as time.
5. **Continuous-gate countercheck — decisive.** A pulse area a yields
   P_blank=cos²(a). Area pi/2 writes; area pi restores blank. Only after
   declared decoupling does isolated-memory evolution preserve records.
   Pointer coherence, reversible erasure and nonrepeatability of the
   system POVM remain explicit. No actual lifetime or arrow is proved.

120 new checks pass; all **6,184 checks across 24 verifiers pass**, including
the independently retained G-R9b side branch and predecessor hash checks.
G-R10's code, card and JSON hashes are unchanged. The five-page TeX compiled
in the built-in editor; exported PDF was visually checked on every page.

**Historical next decision; G-R12 now supplies a controlled local preparation audit:**

- Effective-model option: stop here with a controlled relational-clock
  measurement, and clearly retain its apparatus assumptions.
- Autonomous-event option: specify a realizable recorder interaction and
  a trigger that leaves its contact region after writing. A single-pass
  trigger is a concrete candidate, not yet selected or validated. Include
  its energy, reflected/returning branches and finite-memory recurrence;
  do not silently impose a one-way trajectory or a free external clock.
  A successful finite control model would still not prove fundamental
  time emergence or the source's Ed identification.

Do not automatically add a bath, replace the contact, fit a time law, or
promote any Route-F gate. The invariant-instrument existence task is closed.

Artifacts: [card](data/PREPARATION_RECORD_CARD.json),
[report](output/gr11_preparation_record.md),
[code](code/verify_gr11_preparation_record.py),
[TeX](tex/route_g_preparation_record.tex),
[PDF](output/pdf/route_g_preparation_record.pdf).

### G-R10 — Full observation algebra and stationary conditional clock

Status: **bounded-done for the all-channel readout and its limits**.
No new interaction, bath, selected outlet or empirical fit is introduced.
The G-R9 contact, its three gamma fixtures and all predecessor files are
preserved. The retained graph resonance is an input, not dynamically selected.

- **Full measurement accounting — done.** Number projectors and pair
  quadratures span the 9-dimensional local / 81-dimensional joint
  observation algebra. A fixed seven-outcome phase POVM is positive
  and complete but not tomographically complete. Its joint table and
  remote conditional states retain all nine modes and all outcomes.
  Ignoring A exactly gives the B partial trace, not a selected branch.
- **Exact stationary dark clock — done for this candidate.** The triple
  |0,2>, |1,1>, |2,0> has a two-dimensional contact-null subspace at
  theta=pi/2. Dplus=(|0,2>−2|1,1>+|2,0>)/sqrt(6) has H Dplus=4 Dplus
  for every gamma. It gives p_A=1/7 and the exact conditional orbit
  rho(B|phi)=exp(i N_B phi) rho(B|0) exp(−i N_B phi), without n.
  Its continuous relative-phase density is (1−cos delta)^2/(3 pi),
  a broad cyclic reading rather than a sharp/unbounded timestamp.
  This realizes established conditional-time reasoning in the unchanged
  contact model; it is not a new general Page–Wootters mechanism.
- **Different states kept separate — done.** The actual exchange state,
  the uniform Cesaro unread-update mixture and a prepared pure dark
  state are not identified. At gamma=0.09 the first two give seven-reading
  bare-gap rigidity residuals 0.0277141 and 0.0108673; the dark state
  gives exactly zero. These do not refute general interacting clocks.
  An exact conditional eigenstate equation retains the contact term
  coupling chi(phi) to chi(0), together with the reading likelihood.
  Averaged joint effects commute with H and preserve stationary-state
  statistics, but their autonomous/local apparatus implementation is open.
- **Autonomous pure-state preparation — impossible under this fixed
  unitary contact from |1,1>.** Its conserved dark weight is 2/3.
  A hypothetical energy filter gives Dplus with probability 2/3 but adds
  a measurement instrument and a rejected outcome. Dminus is robust
  against equal-copy flux detuning but inaccessible from this symmetric
  input, and its continuous phase has a twofold ambiguity. Dplus itself
  needs the equal-gap resonance. No dissipation or state selector is hidden.

151 new checks pass, including full site/mode matrices, conditional
orbits, all-outcome sums, stationary effects, the exact contact kernel,
dephased controls and the order-of-limits distinction for weak coupling
versus infinite uniform averaging. The latter prior remains an assumption.
All 5,474 checks across 22 verifiers pass, including two added predecessor
hash-preservation checks. G-R9's code, card and result hashes are unchanged.
The five-page PDF is compiled and visually reviewed on every page.

**Historical next test, now addressed by G-R11:** one preparation plus one-record instrument,
with explicit energy/reference costs and all outcomes. Do not assume a
fresh external clock, filtered success or bath is free. Alternatively
retain the unfiltered stationary mixed state and test its interacting
conditional equation. If an operational completion requires an external
time reference, retain this result as a relational-clock representation,
not a claim that time itself has been generated. No metric/mass/Route-F
gate is promoted; do not start new parameter scans to bypass this decision.

Artifacts: [card](data/RELATIONAL_READOUT_CARD.json),
[report](output/gr10_relational_readout.md),
[code](code/verify_gr10_relational_readout.py),
[TeX](tex/route_g_relational_readout.tex),
[PDF](output/pdf/route_g_relational_readout.pdf).

### G-R9b — Paired channels and two ways to one outlet (side branch, 2026-10-09)

Status: **bounded-done for the declared finite constructions**. Written in
parallel with G-R10 by a separate session. Its density-contact scope is
clarified by G-R12 without changing the original files or restricted results.

1. **The extra outlet is structural — proved.** Flipping every link sign
   maps L to 4−L; on a 3-ring this adds pi to the loop phase, so
   spec L(pi−theta)=4−spec L(theta). The sorted gaps are
   2 sqrt3 sin((pi−theta)/3) and 2 sqrt3 sin(theta/3). For ANY two flux
   3-rings, |1,1>↔|0,2> resonance forces theta_B=pi−theta_A and makes
   |2,0> degenerate as well. Every real density-density contact couples
   the two outlets with equal magnitude; the resonant block is circulant.
2. **Exact 4/9 for the G-R9 vertex contact — proved.** The contact is rank
   one. Its two dark states give P(0,2)=P(2,0)=|G−exp(−i E0 s)|²/9 ≤ 4/9
   at every gamma and every update, also for mirror pairs without swap
   symmetry. This promotes G-R9's weak-contact envelope to an exact bound.
   It is the same structure as G-R10's conserved 2/3 dark weight.
3. **Exit A, distinguishable spectra — done.** A 3-ring at theta=pi/2 and
   a 4-ring at theta_B=4 arcsin(sqrt3/4), with the unchanged vertex
   contact: only the target shares the input energy (reverse outlet
   detuned by (3 sqrt3−sqrt13)/2). Full 12×12 evolution at integer
   updates gives P(0,2)=0.999994927 at gamma=0.01.
4. **Exit B, oriented contact — done.** h=(2,0,1) on identical rings gives
   a closed chiral triangle with loop phase −pi/2 and unit transfer at
   t*=2pi/(3 gamma). Any translation-covariant contact conserves the |T>
   weight 1/3 of |1,1>, so it cannot prepare G-R10's D+ beyond fidelity 2/3.

586 checks pass (numpy only). A resonant flux (A) or an orientation and
weight ratio (B) is supplied, not derived. No time, record or mass claim.
Artifacts: [input card](data/PAIRED_CHANNEL_CARD.json),
[report](output/gr9b_paired_channels.md), [code](code/verify_gr9b_paired_channels.py),
[TeX](tex/route_g_paired_channels.tex), [PDF](output/pdf/route_g_paired_channels.pdf).

### G-R9 — Fixed spectrum, two dials and local-contact countercheck

Status: **bounded-done for the declared finite candidate and audit**.
Emergent time and microscopic selection remain **open**; this is not a
new claim that all preceding conceptual gaps have closed.

1. **Graph and persistence — done.** For S|j>=|j+1 mod 3>,
   L=2I−exp(i theta/3)S−exp(−i theta/3)S† is positive with
   lambda_k=2−2 cos((theta−2 pi k)/3). Local magnitudes and trace are
   fixed while the gauge-invariant loop phase changes the spectrum.
   Mode/gap perturbation bounds are proved. Persistence is not attraction.
2. **Two-copy dials and selected exchange — done under explicit inputs.**
   A uses modes (0,1), B uses (1,2) on independent tensor factors.
   U=exp(−i alpha H), alpha=1/4, presupposes ordered updates n, not seconds.
   Fixed Fourier-mode exchange g=0.1 has an exact invariant four-state
   sector and P<=4g²/(delta²+4g²). One-for-one resonance is theta=pi/2,
   not selected dynamically. At n=63 P=0.99998233 in this idealization.
   For plus-state input the interacting dial coherences both equal
   exp(−i sqrt(3) s) cos(g s): equal phase is not permanent visibility.
3. **Universal slowdown and uniqueness — not established.** With fixed
   trace, one free gap increases while the other decreases. Trace is
   not identified as Ed. Replacing the update generator by the equally
   positive fixed rule L+0.2 L² changes B/A from 1 to 1.4766271 at the
   same graph. The graph alone does not determine event dynamics.
4. **Vertex-contact completion of selected exchange — rejected for this
   contact.** gamma P_vertex0 tensor P_vertex0 has all mode-basis entries
   gamma/9. At resonance, |0,2>, |1,1>, |2,0> are degenerate; the selected
   four-state sector is not closed. Exact swap symmetry gives equal
   outgoing probabilities, each <=1/2. The weak-contact resonant block
   gives 4/9 per target (8/9 combined); this sharper bound is not claimed
   for the full nine-state generator. At gamma=0.09, n=419 the full model
   gives 0.44430738 to each target with 0.000205484 outside the triple.

342 new checks pass; all 5,321 checks across 21 verifiers pass.
The earlier G-R3/G-R4 preservation audits each add one check for the
new verifier; their physical predictions and source code are unchanged.
The six-page PDF is compiled and visually reviewed page by page.
Full site-basis propagation, analytic sector
formulae and contact leakage are checked independently. The six-page
derivation and fixed input card avoid arbitrary particle masses,
dimensional identifications or nonlinear time-law fits.

**Historical next choice, now addressed by G-R10:** retain the honest three-channel contact and
derive its full observation algebra/conditional readout, including
unread channels. An alternative is an independently motivated symmetry
or interaction that forbids the second outlet before reinstating the
selected-mode idealization. Do not launch a precision/parameter scan
to hide this structural difference. Neither choice yet derives event
ordering, permanent records, a spatial metric or a particle spectrum.

Artifacts: [input card](data/EVENT_NETWORK_CARD.json),
[report](output/gr9_event_network.md), [code](code/verify_gr9_event_network.py),
[TeX](tex/route_g_event_network.tex), [PDF](output/pdf/route_g_event_network.pdf).
Older milestone instructions below describe their historical scope;
they do not override this current priority. All Route-F gates stay unchanged.

## New governing decision — G-R relational time (2026-09-23)

The user selected the approach closer to the source document: one
subsystem's events label another subsystem's changes, without assuming
a common physical external time. The old KK model explicitly assumed
Minkowski spacetime; G2-T calculated travel within that spacetime. Those
are not derivations of the source's time concept.

The original sections1.2 and1.8 use Ei accumulation, Ed and a time-event
quantity T. These are not silently identified with an energy eigenvalue,
ordinary energy density, flight time or a KK circle. The source's
negative-energy information and dimensional projection are likewise not
assigned to the clock Hamiltonians or Hilbert dimension.

This is a new conditional quantum construction, **not an edit of the
old KK action**. Preserve all earlier outputs and conditional exclusions.
No new source, detector scan, extra spatial dimension or Route-F bridge
is a prerequisite. Use [the interpretation card](RELATIONAL_TIME_CARD.md).

### G-R1 — Finite Page-Wootters clock and density-only counterexample

Status: bounded-done. Assumed inputs are ordinary quantum mechanics,
a clock/system tensor factorization, two complementary finite energy
spectra, a fixed-total-energy constraint, a matched state and a clock POVM.
They are not derived from Chen's text or from gravity.

For n=0,...,d−1, HC=epsilon(d−1−n), HS=epsilon n and E*=(d−1)epsilon,
the normalized state sum a_n |n,n> obeys K|Psi>=0 for
K=HC+HS−E*. Both local Hamiltonians are nonnegative.
Condition on |phi>=sum exp(i n phi)|n>/sqrt(d):
|psi(phi)>=sum a_n exp(−i n phi)|n> and
i epsilon d_phi |psi>=HS|psi>.
Only after calibration tau=hbar phi/epsilon does this become the
usual time-unit Schrödinger equation. No absolute physical gap is fitted.

The d clock-grid events are orthogonal; the continuous dial is a POVM,
not uncountably many orthogonal readings. The finite clock is cyclic
and cannot obey an exact canonical time/energy commutator. An averaged
joint effect commutes with K and reproduces the conditional probabilities.
One-event statistics are certified, not sequential measurement records.

New diagnostic result: a uniform three-level coherent state and its
energy-dephased mixture have identical local reduced density matrices,
complete joint energy statistics and total energy2epsilon with zero
variance. Yet for a fixed equal-superposition system projector,
P_coh(+|phi)=(1+2 cos(phi))²/9 and P_diag(+|phi)=1/3.
At phi=0 these are1 versus1/3. Thus an ordinary energy density alone
does not specify the relational history. This does not refute a richer
Ei/Ed containing additional state/coherence data, nor prove density
causes slowdown. Candidate Ei must retain relational information if
it is to encode this kind of conditional dynamics.

Common energy rescaling changes the conversion from phase to calibrated
time but not phase correlations; it is not relative clock dilation
without a reference. An energy-zero shift with a matching E* shift
also preserves K in this nongravitational model.

A general clock interaction yields an integral kernel in clock labels.
It does not automatically give a local Ed-dependent clock speed; the
joint state must be solved again and normalized conditioning rechecked.
The interaction projection identity is checked, but no interacting
stationary completion is claimed.

Artifacts: [interpretation card](RELATIONAL_TIME_CARD.md),
[report](output/gr1_relational_clock.md),
[code](code/verify_gr1_relational_clock.py),
[TeX](tex/route_g_relational_clock.tex),
[PDF](output/pdf/route_g_relational_clock.pdf).
New checks109/109; all4115 previous checks remain passing, total4224/4224.
The five-page PDF was compiled twice and every page visually reviewed;
no final layout or reference warnings remain. These checks certify
finite algebra and implementation, not physical realization.

### G-R2 — One explicit record and one final readout

Status: **bounded-done (2026-09-24)**. No G-R3 interaction is activated.

One initial unitary record and one terminal unitary readout act on a
three-level system and two blank-or-binary pointers. The same finite
model gives p(a,b)=tr[N_b U M_a rho M_a† U† N_b†]; it does not multiply
unperturbed one-reading probabilities. A fixed three-level example gives
the joint table [[1,2],[1,5]]/9 and terminal P(+)=2/9, versus zero when
the initial interaction is omitted. A different initial instrument
with exactly the same POVM effects gives P(+)=4/9. Thus both the event
effects AND their recording backreaction must be specified.

**Explicit change of assumptions:** an open history-state propagation
constraint K_hist replaces the noninteracting G-R1 constraint. It has
a unique zero mode on 162 dimensions, with ordered gates, an input
boundary and one epoch bit. G-R1's free U is retained, but its old
fixed-total-local-energy proof is not asserted for this new interaction.
The final readout does not commute with bare HS; apparatus/work
energetics and microscopic gate durations are not supplied.

The cyclic labels k=0,1,2,0,1,2 coexist with saved-event counts
N_evt=0,1,1,2,2,2. The two pointers carry fixed scheduled stamps
(epoch0,k1) and (epoch1,k0), not an arbitrary unknown timestamp.
Once written, their outcomes persist under the declared future links.
Conditioning only on recurrent k=0 mixes blank and completed histories;
including the epoch or completed-record marker resolves this ambiguity.
Full-state recurrence with the same blank input cannot also retain a
fresh nonblank record: the two states are orthogonal (closure residual
sqrt(2)). This does not rule out cyclic clocks with reset, longer cycles
or additional memory. No irreversible arrow, universal event counter
or empirical memory lifetime is derived.

Artifacts: [report](output/gr2_two_event_records.md),
[code](code/verify_gr2_two_event_records.py),
[TeX](tex/route_g_two_event_records.tex),
[PDF](output/pdf/route_g_two_event_records.pdf).
New checks139/139; complete Route-G regression4363/4363.
The five-page PDF is compiled twice and visually reviewed page by page.
These checks certify the declared finite protocol, not Chen's Ed law
or a physically energy-complete autonomous clock.

This stage is closed within its stated scope. Do not expand it into
an unlimited measurement-theory prerequisite. The subsequently authorized
G-R3 candidate/interaction/two-clock test is recorded below.

### G-R3 — An Ei/Ed candidate and an observable clock comparison

Status: **bounded-done (2026-09-24), for one declared candidate**.
No universal density-time law or original-paper identification is proved.

Define Ed=epsilon_D N_D/nu(D), with nu(D)=1 addressed matter cell and
ground-subtracted excitation energy. Units are energy per specified
cell, not assumed spatial volume. Full information still includes
the joint state, event algebra, instruments and reference measure;
Ed is only one observable. The one-cell example cannot experimentally
distinguish density from total excitation energy.

Choose V=kappa H_B Ed/E0, with E0=epsilon_D per cell and a shielded
reference A. The common stationary constraint K=H_A+H_B+H_D+V-E*
has an exact 15x3x2 matched state, positive energies, normalized phase
POVMs and constraint-invariant joint effects. The reference spectrum
is fixed for the synthetic kappa=-1/4,0,+1/4 audit cards. Commensurate
energies, preparation and phase calibration are assumptions, not fitted
data or a solution for arbitrary interacting clocks.

The conditional likelihood peak gives r_B/A=1+kappa d under fixed
calibration. At excited source occupation d=1 the cards yield 3/4,1,5/4.
Common energy rescaling leaves the ratio unchanged; equal fractional
response of two clocks cancels. A single finite-clock outcome does not
give an exact rate or an unambiguous cycle count.

**New physical insight:** density does not fix the direction of the
rate change, and density fluctuations do not generally define one
clock at the mean density. For the +1/4 card an unread equal source
mixture gives a specified outcome P0=5/9; replacing Ed by its mean
incorrectly gives (3+2sqrt(2))/9. Reduced visibility and multiple rate
branches must be retained. B-only data do not distinguish source
coherence from a classical mixture with the same populations.

A separate Jaynes-Cummings exchange calculation derives chi=g²/Delta
and kappa_eff=2chi/(omega-chi). Exact excitation-block spectra verify
that either detuning sign is possible and control the dispersive
energy remainder. This does not implement the rational kappa=1/4 card.
Higher-order gaps are transition dependent; dressed density/readout
operators matter. The result is an environment-dependent clock shift,
not universal time dilation, gravity or a fitted slowdown law.

Artifacts: [density/convention card](DENSITY_CLOCK_CARD.md),
[report](output/gr3_density_clocks.md),
[code](code/verify_gr3_density_clocks.py),
[TeX](tex/route_g_density_clocks.tex),
[PDF](output/pdf/route_g_density_clocks.pdf).
New checks273/273; complete Route-G regression4636/4636.
The six-page PDF is compiled twice and visually checked on every page.

The G-R2 history constraint and this stationary energy constraint are
not silently merged: an autonomous energy-complete recording apparatus
remains absent. That is a limitation, not a reason to reopen every
measurement-theory detail before reporting the bounded result.

### G-R4 — Independently fixed physical source and universality criterion

Status: **bounded-done (2026-09-24), source/response audit only**.
No new clock experiment, physical apparatus or universal-time effect.

Choose one isotropic blackbody radiation bath and two transitions of
the same 171Yb+ ion: E2 and E3 from the common S ground level.
H_int=-d.E fixes the mechanism. Independently measured differential
scalar polarizabilities, converted to alpha_excited-alpha_ground, are
6.9(1.4)e-40 and 0.888(0.016)e-40 J m^2/V^2. Both positive means both
thermal electric-dipole shifts are negative. The E2 paper originally
uses the opposite sign convention. These are pinned historical inputs,
not a claimed newest calibration or a fit to a temperature-cycle signal.

**Explicit boundary:** conventional radiation density u_gamma=a_rad T^4
has units J/m^3 and assumes physical space. It is NOT silently equated
to G-R3 energy per operational cell or Chen Ed. This is a physical
control on the interpretation, not a derived realization of emergent time.

From second-order perturbation theory derive
delta_nu_i=-(2h epsilon0)^(-1) integral Delta_alpha_i(omega) u_omega d_omega.
The leading fractional coefficient is
beta_i=-Delta_alpha_i(0)/(2h epsilon0 nu_i). At 300 K, the static scalar
electric-dipole prediction gives fractions -5.235e-16 and -7.223e-17:
their magnitudes differ by about7.25, not by numerical rounding.
The illustrative ratio modulation R(310 K)/R(300 K)-1 is +6.325e-17,
where R=nu_E3/nu_E2. This is a prediction, not an observed signal.

For a connected graph of levels, equal fractional gap shifts require
delta_H_secular=f H0+C I. The independently fixed thermal coefficients
fail that necessary condition. Thus this selected thermal response
cannot be explained solely as one universal multiplicative time factor.
An additional exactly common factor F(u) remains **unidentifiable**:
nu_i=F(u) nu_i^0(1+s_i) cancels F from the ratio exactly. Neither
nonconstant nor constant same-bath ratios alone settle universal time.

The leading T^4 coefficient mismatch is separate from full dynamic
corrections. The published E3 correction is recorded; complete E2
spectral response, covariance and apparatus controls remain absent.
A +/-3 published-sigma coefficient box plus assumed +/-10% independent
dynamic factors preserves the mismatch at all16 corners. This is a
conditional stress bound, not statistical significance or a derived
bound on omitted physics. Do not reuse BBR-corrected outputs as raw
temperature-response observations.

Artifacts: [pinned inputs](data/BBR_CLOCK_INPUTS.json),
[report](output/gr4_bbr_universality.md),
[code](code/verify_gr4_bbr_universality.py),
[TeX](tex/route_g_bbr_universality.tex),
[PDF](output/pdf/route_g_bbr_universality.pdf).
New checks91/91; complete Route-G regression4728/4728.
The old G-R3 script now counts274 checks because its automatic file-
preservation inventory also includes G-R4; its original273-check stage
snapshot above is historical. No earlier verifier code was changed.
The five-page PDF is compiled twice and visually checked page by page.
These validate the calculation, not an empirical discovery.

### G-R5 — Acquired data, independent reference, and empirical eligibility

Status (2026-09-24): **data acquisition/audit and reference model bounded-done;
requested two-temperature empirical validation OPEN**.
Checkpoint before this work: commit9b66092 pushed to
origin/codex/route-f-invariant-control, preserving milestones through G-R4.

Downloaded and checksum-verified two primary releases. The PTB2021
archive has11 published E3/E2 points. The ROCIT2025 release supplies
34 selected same-package E2/E3 files:2,324,201 rows, all flagged valid,
zero nonfinite values and duplicate timestamps. All selected days and
release metadata were audited. The small PTB archive and exact first
ROCIT day (51,990 rows) plus YAML are preserved; every selected full-
release member hash is recorded. No physical frequency slope was fitted.

**Missing:** paired bath-state labels, independent ion-site radiation
calibration and a per-run BBR correction ledger. The released processed
comparator is not certified as an uncorrected servo record. Its direction
is E2 versus E3 and its scaled-transfer-beat convention must be respected.
The older PTB ratio offset is absolute, not fractional. Validity does not
imply sample independence. This is a scope mismatch for our experiment,
not bad-data discovery or proof that suitable data do not exist anywhere.
The September2025 deposit supersedes the paper's earlier on-request
availability statement. Motional temperature is not radiation temperature.

**Reference selected:** a second independent Yb+ E3 clock outside the
target enclosure, in its own monitored enclosure, with a comb and
phase-stabilized optical paths. Both target transitions are compared
to it. Its BBR response, other clock shifts, possible common-factor
response, drift, gravity and optical-transfer terms remain explicit.
This specifies a measurement model, not a constructed apparatus or
physical implementation of G-R1--G-R3. The observable is a change in
log(F_target/F_reference), never an unreferenced absolute time factor.

Derived differential and common-response channels. At static leading
order, E2/E3 weights(-0.16004,1.16004) cancel the target T^4 term while
retaining the common target-reference response. They do NOT remove
reference drift or a common path/gravity bias. An unconstrained common
bias is exactly degenerate with the desired common signal. A derived
third ratio adds no independent equation. Radiometry must remain
independent of the frequency law under test; two channels fitted to
two unknowns cannot validate that law by themselves.

For an illustrative300->310 K change, the pinned polarizability errors
alone give a2.17e-18 to2.59e-18 composite standard-error range over
unknown correlation, not a confidence interval or achieved total error.
Thus the retained inputs do not support a sub-1e-18 common-effect claim.
An ABBA sequence cancels offset/linear drift but not quadratic drift
or heater-correlated bias; matched effective windows and covariance
are required. No new stochastic scan or free time-rate law is introduced.

Artifacts: [source provenance](data/clock_comparison_audit/PROVENANCE.md),
[full archive audit](output/gr5_clock_data_full_audit.json),
[reference/data contract](data/EXTERNAL_REFERENCE_CARD.json),
[request checklist](data/clock_comparison_audit/DATA_REQUEST.md),
[report](output/gr5_external_reference.md),
[TeX](tex/route_g_external_reference.tex),
[PDF](output/pdf/route_g_external_reference.pdf).
New model checks46/46. Current full Route-G regression4776/4776;
G-R3/G-R4 each automatically add one preservation check for the new script.
The four-page PDF is compiled twice and visually checked page by page.

### G-R6 — Independent thermal and gravitational controls (2026-09-29)

Status: **theory baseline bounded-done; empirical validation OPEN**.
User-approved scope: retain E2/E3 thermal response and add a known
gravitational-potential comparison, independently controlled.
G-R4 inputs and G-R4/G-R5 verifier files are pinned by hash and unchanged.
The independent E3 reference remains fixed and fully modeled.

In a static metric, d tau=N dt and received clock ratios contain N_S/N_R.
Taking logarithms gives q_i=G+B_i(T_S)-B_3(T_R)+ordinary/link terms.
G is specified by the metric, not a freely fitted slowdown factor.
Two local temperatures crossed with two potential settings give
temperature, height and mixed contrasts. After controlled corrections:
thermal_i=Delta B_i, height_i=Delta ln N, mixed_i=0.
The E3/E2 internal ratio rejects common gravity but retains the
differential thermal response. A raw fractional mixed product is
ordinary multiplication, not a log-interaction anomaly.

The illustrative 300/310 K, 1 m example uses conventional g=9.80665
m/s^2 only as a scale: gravity +1.091e-16 in both channels; heating
-7.337e-17 in E2 and -1.012e-17 in E3. No site geodesy or new clock
observations exist for this card. One cannot derive DeltaPhi from the
same clocks and then call the gravitational prediction independently tested.

As a distinct real external scale, Grotti et al. (2024) report an
independent geodetic potential magnitude of3915.88(0.30) m^2/s^2.
Its universal redshift prediction is4.3570e-14 with3.34e-18 propagated
potential uncertainty. Their clocks were Sr, not the proposed E2/E3
experiment; no published performance or clock-inferred potential is
used as our calibration. This anchor does not fill the four-setting gate.

**New useful insight:** matching the local radiation conditions between
heights cancels even uncertain fixed BBR shifts. The older polarizability
error is not automatically a floor on this height test. Upper-site +0.1 K
mismatch at both setpoints would instead leak about -7.34e-19 (E2) and
-1.01e-19 (E3). Those are conditional sensitivities, not achieved errors.
Nonthermal transport, geometry, reference and link systematics remain.

A common height-locked nuisance has the same design column as a common
redshift anomaly. Only independently bounded nuisances permit a physical
interpretation. Nonzero mixed contrast rejects controlled separability,
not automatically EEP. Four means with a fully fitted interaction model
are saturated; repeat blocks and retain covariance before significance.
Matched sampling windows, reversals and return-to-start controls matter.
An ideal symmetric eight-plateau sequence rejects linear drift but not
quadratic or state-locked drift. No hardware timing precision is assigned.

Independent thermostats do not imply one global passive thermal bath.
The latter obeys Tolman-Ehrenfest T_local*N=constant; it must not be
imposed simultaneously with independent local setpoints. This remains
a classical-background redshift/LPI baseline, not a test of every EEP
component or a quantum mass-operator/off-diagonal experiment.

Artifacts: [input/acquisition card](data/EEP_BASELINE_CARD.json),
[numerical report](output/gr6_equivalence_baseline.md),
[verifier](code/verify_gr6_equivalence_baseline.py),
[TeX](tex/route_g_equivalence_baseline.tex),
[PDF](output/pdf/route_g_equivalence_baseline.pdf).
Completion audit (2026-09-29): 49/49 new implementation checks and
4827/4827 checks across all 18 Route-G verifiers pass. G-R3/G-R4 each
automatically add one preservation check for the new verifier; no old
physics formula or pinned thermal input changed. The four-page PDF was
compiled locally and visually reviewed page by page without layout or
reference warnings. The built-in compiler did not return a confirmed
result, so PDF validation used the local compiler. No empirical gate
is promoted by software tests.

### G-R7 — Same-action phase and clock-path correlations (2026-09-29)

Status: **conservative quantum kernel bounded-done; complete physical
interferometer and empirical validation OPEN**.
This is the user's approved quantum extension, not another fitted
time-rate formula. The G-R4 input card and G-R6 card/verifier are
hash-pinned and preserved. The proper-time conservative action gives
theta_ir=2*pi*nu_i*integral[(1+b_i)*d_tau_r]; its derivative recovers
the G-R6 frequency. The joint controlled unitary gives
kappa=Tr(rho_clock U_B^dagger U_A), the internal fringe phase and
ideal visibility |kappa|. A full scalar path/ground/guide/laser phase
is not fixed by differential polarizabilities and remains explicit.

For pure (g+e_i)/sqrt(2), kappa=exp(-i delta_i/2)*cos(delta_i/2),
and joint-state negativity is |sin(delta_i/2)|/2. A common-ground
(g+e2+e3)/sqrt(3) instead has kappa=(1+exp(-i delta2)+exp(-i delta3))/3;
it is not a product of independently invented qubits.
Frequency, phase and visibility are linked consistency predictions,
not three freely adjustable fits.

**Exact diagnostic counterexample:** with Ua=I and Ub=Z, pure and
energy-dephased clock preparations have identical zero path visibility.
The pure joint state has negativity1/2; the dephased one is explicitly
separable. The joint witness (I-X_path Z_clock-Z_path X_clock)/2
distinguishes this calibrated fixture. State-resolved energy fringes
alone still do not certify entanglement; joint noncommuting measurements
or tomography and preparation characterization are required.

A common final unitary, or unconditioned trace-preserving channel acting
only on the clock, leaves the path marginal unchanged. It can change
entanglement without changing path visibility. Thus an assumed Ramsey
decay factor cannot simply be multiplied into the path signal.
Reversing the ideal controlled unitary can undo its correlations,
but needs actual controls; postselection is not unconditional recovery.
The conservative real BBR response supplies no physical radiative
correlation function, emission which-path overlap or full dissipator.

**Same-weight phase check:** x_i=delta_i/(2*pi*nu_i)=Delta_tau+beta_i*J,
where J=integral(lambda_A*u_A-lambda_B*u_B)dt. The unchanged G-R5
weights yield sum w_i*x_i=Delta_tau within the static response.
This requires calibrated unwrapped phases, matched histories and
independent radiometry; it is not a weighted sum of visibilities or
an empirical determination of proper time from missing data.

Illustrative 1 m/10 ms hold: gravity-only Delta_nu=(0.07511,0.07006) Hz,
delta=(0.004719,0.004402) rad, ideal visibility loss=(2.784,2.422)e-6.
Thermal-only and combined controls are also computed, without scans.
Complete launch/transport/recombination phases and wave-packet closure
are not supplied by this hold. Standard freely falling geometries may
cancel uniform-field sensitivity; paths and initialization matter.
At 1 m the ideal E2/E3 first contrast zeros need6.66/7.14 s.
E2's tens-of-ms upper-state lifetime makes such unperturbed E2
extrapolation inappropriate. Isotope lifetime values are not pooled
into a made-up environmental rate. E3 longevity alone does not establish
charged-ion path coherence or hardware feasibility.

Artifacts: [control card](data/CLOCK_PATH_CARD.json),
[numerical report](output/gr7_clock_path.md),
[verifier](code/verify_gr7_clock_path.py),
[TeX](tex/route_g_clock_path.tex),
[PDF](output/pdf/route_g_clock_path.pdf).
All 66 new checks and all 4895 checks across the 19 Route-G verifiers
pass. The total includes two automatically added preservation checks
in G-R3/G-R4. The five-page PDF was compiled locally and every rendered
page was visually reviewed; the final compiler log has no layout warnings.
These are algebra/numerical/preservation checks, not experimental evidence.

Next quantum step: specify one closed-path/guide/pulse sequence and
state-resolved differential phase readout, with coherent/dephased
preparation controls. The phase is first order in small delta while
visibility loss is second order. Joint measurements follow if
entanglement is claimed. Do not prioritize forcing full visibility
collapse, treating all noise as time dilation, or fitting another Ed law.
This does not require relabeling the G-R6 ordinary two-site comparison
as a coherent spatial-superposition experiment.

### G-R8 — Closed guided protocol, total phase and joint controls (2026-09-29)

Status: **specified conservative effective sequence bounded-done;
real controls, environmental dynamics and empirical gates OPEN**.
Prior G-R5--G-R7 work was committed and pushed as `2aab1f1` before this step.
The G-R4 input card and G-R7 card/verifier remain hash-pinned.

Choose two orthogonal guide modes at common initial/final height.
A quintic 5 ms outward ramp reaches symmetric +/-0.5 mm, followed by
10 ms hold and 5 ms return. The actual guide potential is
V_r=m*omega^2*(z-q_r)^2/2-m*(qddot_r+g)*z, with explicitly zero scalar
arm offsets. The 1 kHz compensated harmonic dynamics closes the ground
wave packets. Preparation and both path couplers occur at common height;
their ideal pulse matrices and phase conventions are included.

The full scalar action cancels the ground arm phase. The enclosed area
is d*(hold+ramp)=1.5e-5 m s, not d*hold. The E2/E3 phases are
(7.07889,6.60339)e-6 rad, and total port phases are
Phi_g=alpha_f-alpha_s, Phi_e=Phi_g-delta. Static matched BBR factors
are retained. No temperature mismatch or new clock law is introduced.

**New physical insight:** specifying the same paths and forces does
not determine the scalar fringe. Replacing z by z-q_r in the linear
guide potential leaves forces unchanged but adds a different arm
offset and about -3.96e5 rad of ground phase. This is not an arbitrary
phase convention unless coupler phases are transformed consistently.
The chosen cancellation is conditional on the declared actual potential.

The same guide is audited with m_e=m*(1+eta). Commensurate ramp/hold
oscillator periods cancel the first-order mass-closure residual.
Exact finite-mass harmonic Gaussian residuals, endpoint phase and
overlap are calculated and checked against an independent ODE.
They are tiny for the card, but nonzero and not an all-orders GR result.
This addresses the ideal closure question without claiming physical
trap precision at the many decimal places of the diagnostic.

Pure and energy-dephased preparations again have identical path signals.
After explicit local phase calibration, W=(I-XX+ZY+YZ)/4 is a partial
transpose of a positive projector, hence nonnegative on separable states.
Its ideal pure value is -sin(delta/2)/2; the dephased value is
(1-sin(delta/2))/4. XX, ZY, YZ commute globally but need incompatible
local settings. Added ZZ genuinely fails to commute with ZY. Population
readout following declared local rotations implements each setting.

**Feasibility decision, not another mathematical gate:** for this card
the ideal witness requires about (3.59,4.13)e11 independent shots per
setting for a three-standard-error separation, before systematic errors,
lifetime, loss and rate constraints. This is not a confidence bound or
experiment. The phase itself is also small. An unknown state-dependent
guide term is degenerate with the desired shift; averaged differential
frequency bias must be below roughly (5.6,5.3)e-5 Hz just to be smaller
than the signal. Calibrate or reject a specific guide/coupler candidate.
Do not force an entanglement claim or cure poor feasibility by another
time law. G-R5/G-R6 empirical gaps and Route-F promotion gates are unchanged.

Artifacts: [control card](data/CLOSED_PATH_CARD.json),
[numerical report](output/gr8_closed_path.md),
[verifier](code/verify_gr8_closed_path.py),
[TeX](tex/route_g_closed_path.tex),
[PDF](output/pdf/route_g_closed_path.pdf).
All 80 new checks and all 4977 checks across the 20 Route-G verifiers
pass (including two new automatic G-R3/G-R4 preservation checks).
The six-page TeX compiled successfully with the desktop editor; a local
PDF export was also compiled and all six rendered pages visually reviewed.
No final compiler layout warnings. Numerical agreement is not empirical
verification. The G-R8 checkpoint includes the protocol, derivation,
numerical results and synchronized roadmaps.

### Separate empirical action — not another mathematical prerequisite

Obtain the specific missing thermal/clock/correction records using the
saved G-R5 checklist, supplemented by the G-R6 independent potential
and four-setting acquisition fields, or conduct a genuinely controlled
acquisition with independent radiometry and geodesy. No lab contact
has been sent. Existing public
frequency comparisons must not be relabeled as the requested experiment.
Without those records, empirical validation stays open; no slope,
significance, common-time bound or physical-fit promotion is authorized
by the present data. If external coordination is desired, obtain user
direction before contacting researchers.

With qualifying records, first test the differential thermal prediction
and matched-temperature gravitational prediction with matching windows
and known corrections, then test mixed residuals and common residuals
against reference, link, transport and heater controls.
Keep source-only versus common-factor claims separate. No new detector
scan, emergent-space proof, or arbitrary Ed slowdown fit is a substitute.

### Beyond the bounded clock program

Space, Lorentzian causal structure, universal proper time, gravity,
mass generation, original negative-mass/time-jump claims, and a physical
implementation remain unproved. They are not all made mandatory gates
for G-R2, and no existing Route-F result is promoted.

## Preserved comparison program — G0–G2-T

The historical priorities and next steps below apply only within their
declared conventional-spacetime branch. They do not supersede G-R or
silently become prerequisites for relational time. G3/G4 at the end
likewise remain optional continuations of the old KK prototype.

## G0 — Action and interpretation card

Status: bounded-done.

- One fixed M4 × S1 background, one complex scalar, positive kinetic term.
- Shared R, M5 and alpha; alpha modulo integers is a flat holonomy.
- One pair of source couplings at 0 and pi, with a specified path and
  parallel transport; no per-mode detector tuning.
- Four numerical engineering cards, chosen without experimental masses.
- Scope: free-field spectrum and linear source response only. No dynamical
  gravity, gauge bosons, radion potential or stabilization is assumed.
- The Chen final document's common-energy-level idea motivates the
  candidate; its other assertions are not treated as established axioms.

## G1 — Spectrum, visibility and covariance

Status: bounded-done; implementation independently rerun, 331/331 checks.

Derive and verify:

1. The KK mass tower and positive-energy Hamiltonian.
2. The fixed probes' full 2 × 2 residue matrix, including cross-response.
3. Passive basis covariance, small-gauge covariance and integer
   large-gauge relabeling.
4. A finite-difference spectrum independent of the analytic formula,
   and a bounded KK response truncation.
5. The two physically meaningful restrictions:
   - signed-branch squared masses have second difference 2/R²;
   - the two probe residues have a nonzero fixed sum for every signed
     Fourier branch, not every superposition in a degenerate subspace.

Do not misread adjacent SORTED energy levels as adjacent signed KK labels.
At degeneracy, compare spectral projectors/residue sums, not arbitrary
eigenvectors. Relative probe weights are not branching fractions.

Completion requires one executable report with spectra, response and
failed-claim flags. Numerical checks certify the implementation, not the
existence of an extra dimension in nature.

Delivered:

- One concise action-to-response derivation, TeX and typeset PDF.
- Four predeclared cards, full complex two-channel residues and responses.
- Independent 32/64/128-site diagonalization. The last spacing-halving
  error ratios are 3.9846--3.9882, as required by the O(h²) symbol.
- Gauge-transformed response error <=1.4e-15 in the recorded run;
  passive-basis response error <=2.8e-16.
- Infinite image-sum reference and rigorous per-entry truncation bound,
  tested at 16/32/64/128/256 modes on each side and three momenta.
- A physically important correction: at alpha=0 the n=±1 coupling map
  has rank 1, so a sine standing wave is dark to BOTH probes. At alpha=1/2
  the n=-1,0 map has rank 2. Full pole residues, not an arbitrary choice
  of degenerate eigenvectors, are basis-invariant.

This milestone is closed within its declared free-field scope. Do not
reopen it merely to gain more decimal precision or fit particle masses.

## G2 — One genuine transition, not a change of coordinates

Status: bounded-done (tree-level EFT); 262/262 independently rerun
implementation checks. G1 regression: 331/331.

Chosen action: add a neutral complex dynamical field X with positive
kinetic term, mass MD and potential kappa5 |Phi|² |X|², kappa5>=0.
This is local in five dimensions. Localization belongs to the detector
quantum state, not an externally pinned coupling or infinite-mass sink.
The classical vacuum remains Phi=X=0; G1 is the tree quadratic baseline.

The exact circle overlap gives g=kappa5/(2 pi R) and

\[
\phi_n+X_l\longrightarrow\phi_m+X_j,\quad j=l+n-m,\quad
\mathcal M_{\rm tree}=-g.
\]

Every component obeys four-momentum and compact-momentum conservation.
The open-channel cross section is g² p_f/(16 pi s p_i); no identical
particle factor is needed. Inverse detailed balance uses the SAME s.

The earlier G1 candidate remains rejected kinematically: for the
bulk vertex chi Phi-dagger Phi with a massless neutral five-dimensional
chi, KK momentum conservation requires q=n-m and m_chi=|n-m|/R.
The triangle inequality gives

\[
m_n < m_m+|n-m|/R \quad(M_5>0,\ n\ne m).
\]

Thus the proposed one-to-two decay is closed. At M5=0 it is forbidden
or threshold-only, not an open two-body decay. This excludes that process,
not all interactions or scattering. No stable interacting completion
is inferred from a cubic vertex alone.

The G2 repair is physical: for l=0 the detector's internal recoil costs
sqrt(MD²+(n-m)²/R²)-MD rather than |n-m|/R. With
Delta=mn-mm>0 the rest-limit channel is exothermic iff

\[
M_D>\frac{(n-m)^2/R^2-\Delta^2}{2\Delta}.
\]

The engineering card R=1,M5=.5,alpha=.25,MD=5,g=.05,p_i=.2 gives
Phi_1+X_0 -> Phi_0+X_1 with p_f=1.02175458 and
sigma=6.27164081e-6. Its rest-limit critical MD is only 0.24146563:
this opening is not a near-threshold decimal adjustment.
The released Phi rest energy .78727421 pays .09901951 internal detector
recoil plus .68825469 additional ordinary kinetic energy.

Preparation: a finite normalized detector packet with l=-4,...,4,
c_l proportional to exp[-l²/(4*1.2²)] exp[-i l*.7]. It freely spreads.
The exact first-order packet amplitude and traced-recoil probability
are derived in the new TeX. Different l give orthogonal outgoing j
for fixed n,m. The inclusive rate cannot measure the initial packet
phase or center against a single incident KK plane wave.

Numerical scope: component cross sections, all open final modes,
conservation laws, inverse detailed balance, localization diagnostics
and a common-overlap averaged sigma*v. Do NOT label the latter an
absolute finite-encounter probability; four-dimensional envelopes and
encounter geometry are still needed. The five-dimensional EFT has no
certified loop/cutoff uncertainty or radiative naturalness result.
No interacting all-orders G1 identity or physical particle fit is promoted.

Results: the common-overlap packet rate is 1.11933296e-6 for target m=0
and 1.85661445e-6 summed over all scattering channels, including elastic.
The 0.60288929 ratio is a prepared-state rate fraction, not an absolute
encounter probability. Conditional incoming/outgoing detector mean
momenta differ by exactly +1; weighting both by the selected event rate
avoids confusing preparation bias with recoil.
Free packet localization decreases from .91640536 to .48041872 in
the first circular moment between t=0 and 5, with energy and norm
unchanged. This is not a pinned defect or a derived localized ground state.

Artifacts: [derivation](tex/route_g_local_conversion.tex),
[PDF](output/pdf/route_g_local_conversion.pdf),
[numerical report](output/g2_local_conversion.md),
[executable](code/verify_g2_local_conversion.py).

### G2-R — Recoil-resolved continuation

Status: bounded-done at tree level (2026-09-23). The reproducible verifier
passes 357/357 checks, independently rerun; the four-page derivation PDF
has been compiled and visually reviewed. This closes the joint prediction
under the declared preparation, not an experimental detector model.

Keep the G2 action, g, all masses, alpha, common incoming momentum and
finite detector packet unchanged. For every open channel, compute
W(m,j)=w_l K(l,m), l=m+j-n. Normalize separately over all scattering
events and over mode-changing events m!=n. Elastic scattering is not
the identity/no-encounter probability.

Deliver full signed joint matrices, mass-only folded matrices, marginals,
Bayes conditionals, event-conditioned incoming l, ordinary recoil
energies/momenta and the isotropic tree angular law. All 25 open channels
are retained. The all-scattering output weights for m=(-2,-1,0,1) are
(.00785985,.27119405,.60288929,.11805680); m=1 is elastic.
Mode-changing events comprise .88194320 of scattering events under
the common-overlap prescription, not of arbitrary incident packets.

New physical insight: the X mass alone determines |j| and loses its
sign. At this nondegenerate Phi card, a known incoming p_i and the
outgoing ordinary p_f can reconstruct the sign for m!=n. With d=m-n,

\[
l^2=R^2[(E_{\rm out}-E_n)^2-M_D^2-p_i^2],
\qquad j=(l^2-d^2-|j|^2)/(2d).
\]

This is MODEL-DEPENDENT reconstruction. Do not infer j from this formula
and then claim the assumed conservation law has been independently
verified. Independent observables are the predicted joint line positions,
fractions and angular dependence. The m=0, |j|=1 pair has ordinary
momenta 1.02175458 and 1.36191317 (gap .34015860). A deterministic
momentum error below half the gap is a sufficient ideal resolution
requirement, not an assigned experimental uncertainty.
Elastic events cannot resolve the sign this way. At a degenerate Phi
holonomy, even the initial output-mode assignment needs a new audit.

No packet phase or localization is witnessed by this diagonal joint
distribution: the dephased preparation gives the same prediction.
No actual detector response, efficiency, loop-error model or absolute
4D encounter probability has been calculated. No parameters are fitted.

Artifacts: [joint report](output/g2_recoil_joint.md),
[TeX](tex/route_g_recoil_joint.tex),
[PDF](output/pdf/route_g_recoil_joint.pdf),
[executable](code/verify_g2_recoil_joint.py).

The finite-readout continuation and activated drive are now recorded
below. They do not overwrite this undriven joint prediction.

### G2-S — Finite-resolution recoil discrimination

Status: bounded-done (2026-09-23); 272/272 implementation checks.

Use a classical incoherent incoming COM momentum ensemble with
parent-Gaussian location p0=.2, widths sigma_i=0,.02,.1,.2,.4, truncated
to p>=0. sigma_i is not the truncated distribution's actual deviation.
Use a Gaussian real-valued reconstructed momentum estimator with
sigma_d=.02,.05,.1,.2,.4; a negative fitted estimator is allowed.
Mass labels m=0 and |j|=1 are assumed exact, with equal acceptance,
common COM frame and the same external-overlap prescription.

Integrate the unnormalized rates w_l f(p) K_lmq(p) R(x|p_f) first.
Then normalize the selected pair, derive the posterior and compute the
optimal equal-cost sign error integral min(h_plus,h_minus)/Z.
The engineering decision criterion is <=5%, not a discovery threshold.
Compare to the event-selected majority-only baseline (~22.5%).
The 50-card scan covers both undriven and driven cases, with ideal
independent q-tagging as a mathematical information bound only.

| sigma_i | sigma_d | Undriven error | Driven q-unknown error | Ideal q-tag error |
|---:|---:|---:|---:|---:|
| .1 | .1 | 3.77384% | 5.94485% | 3.77589% |
| .2 | .1 | 4.84993% | 6.83823% | 4.85210% |
| .4 | .1 | 11.0294% | 12.2710% | 11.0331% |

At (.1,.1), the undriven card passes the declared 5% criterion, but
the driven card without q information fails. At (.1,.05) the driven
untagged error drops to 2.20225%. Extra detector Gaussian noise cannot
improve optimal classification (L1 contraction); beam broadening also
changes selected priors, so do not assert the same general monotonicity.

Numerical controls: exact sharp-Gaussian Bayes comparison, G2-R hash/rate
regressions, truncated-normal moments, nonlinear vs narrow-width
derivatives, quadrature-order/boundary-grid comparison, and analytic
incident/response tail bounds. Independent adaptive integration confirms
the main card. Root searches are convergence-tested, not certified
exhaustive for arbitrary mixture parameters. Gaussian tails are bounded
under the formal tree kernel, not evidence for infinite UV validity.

Artifacts: [report](output/g2_resolution.md),
[code](code/verify_g2_resolution.py), [full data](output/g2_resolution.json).

### G2-D — Activated time-dependent energy reservoir

Status: bounded-done at first Born/long-time order (2026-09-23);
565/565 checks. User activation supersedes the former backup-only status.

Choose the minimal explicit local modulation
g(t)=g[1+epsilon cos(Omega t)], epsilon=.5, Omega=.2.
It is spatially/circle homogeneous, with pump rest frame equal to the
incoming COM. The instantaneous quartic remains positive; free masses
are unchanged. The new external reservoir is a hypothesis, not a
closed UV completion or dynamically selected geometry.

At one vertex g0=g, g(+/-1)=g epsilon/2, Eout=Ein+q Omega, and
K_lmq=|gq|^2 p_f/(16 pi E_n E_l Eout). Spatial and compact momenta
remain conserved. Every event records pump work q Omega; matter alone
is not energetically closed. The epsilon=0 rate regresses to G2.

All 66 open rows are retained (14/25/27 in q=-1/0/+1).
Total K=2.10563491267e-6; mode-conversion share=.87381701.
Mean absorbed pump work/scattering=.00900097215 and work coefficient
1.89527612162e-8 are not a luminosity-normalized power.
Mode-preserving m=n,q!=0 events are inelastic; only m=n,q=0 is elastic.

The selected two passive lines become six components. Opposite-sign
sidebands at p=1.20257369 and 1.18924247 have gap .01333122.
No exact degeneracy occurs at the chosen card; perfect sharp-beam
readout still distinguishes the full catalogue. The passive inversion
must include unknown q, and finite resolution can worsen classification.
Omega=.1,.2,.4 are sensitivity cards, not a frequency fit.
The analytic exact-crossing frequency .19243951642 is a diagnostic,
not used to tune the background.

The code checks first-order Fourier coefficients, independent radial
phase space, full enumeration, gauge/compact/spatial conservation,
on-shell pump work, positive quartic and finite top-hat Fourier factors.
The latter are NOT finite-collision probabilities. Fixed Omega>0 and
long-time rates are assumed; Omega=0 requires coherent harmonic addition.
Coherent finite beams/pulses need amplitude-level interference and cannot
be replaced by this classical momentum histogram.

Artifacts: [drive report](output/g2_driven.md),
[code](code/verify_g2_driven.py), [full data](output/g2_driven.json),
[shared derivation](tex/route_g_resolution_drive.tex),
[PDF](output/pdf/route_g_resolution_drive.pdf).

The readout choice is made in G2-M below. Finite 4D packets remain
unnecessary for these conditional long-time rates; they are needed
for absolute encounters, finite pulses or coherent interference.

Checkpoint 6bbd6fa saved and pushed G2-R before this continuation,
as requested. G2-S/D does not promote a Route-F fit or identify SM fields.

### G2-M — Momentum-only readout selected; noisy tag comparison

Status: bounded-done as an engineering contract (2026-09-23).
Frozen-readout checks: 154/154. Noisy-tag checks: 125/125.
Both verifiers were independently rerun; the four-page TeX/PDF was
compiled and visually reviewed.
Physical apparatus realization remains open, not claimed by a Gaussian
response or the name "detector" attached to the recoil field X.

Main contract: unchanged driven action, sigma_i=.1, nominal sigma_d=.05,
same ideal m=0, |j|=1 selection. Freeze the nominal rule
j_hat=+1 iff x<1.25892910. Nominal total error is 2.20225%, with
class-conditional errors (1.0315%,6.2279%) for true signs (+1,-1).
Aggregate <=5% does NOT establish >=95% recall for each sign.
Do not silently change the loss function or priors to hide this.

Calibration response: x=(1+a)p_f+b+N(0,sigma^2).
The continuous stress box is |a|<=.01, |b|<=.02, .04<=sigma<=.06,
in model units. The classifier is NOT retrained at each unknown
calibration. The worst of 27 sampled cards is 3.8648%, not a proof
of a global supremum. The unsplit CDF bound is loose (5.3902%).
One predetermined 2x2x2 shared-nuisance cover reduces the upper envelope
to about 4.46%, sufficient for the overall 5% criterion.
Kink-aware integration is convergence-tested with a separate analytic
tail bound, not interval-certified. No repeated refinement or fit.

Oracle retraining is displayed separately. With true a,b known,
an affine transformation makes oracle error depend only on sigma/(1+a).
This is not deployed robustness. At a=b=0 the nominal 5% noise limits
are about .0907 with retraining and .0891 with the frozen rule;
they are not guarantees for arbitrary calibration.

Alternative: T_eta(z|q)=(1-eta)delta_zq+eta/3, true mistag=2eta/3.
At sigma_d=.1 the 5% crossing lies in eta=[.1758,.1768].
This requires the FULL specified matrix, not merely 88.2% tag accuracy:
constant z=0 already scores 88.98% but adds no information and leaves
5.94485% recoil-sign error. Tests preserve rates/priors and verify
stochastic degradation and ideal/no-tag limits. These tag comparisons
assume a known calibrated matrix and its optimized classifier,
not robustness against unknown tag drift.
At primary sigma_d=.05 no extra tag is needed for the overall criterion.

A tag made only from the same measured x cannot add information;
a drive clock records phase, not event work. Extra sensor data need
one joint likelihood over shared latent momentum, not a product of
separately marginalized likelihoods. No real sensor coupling, physical
energy unit, mass-confusion/acceptance calibration or quantum pump
counter has been derived.

Artifacts: [frozen-readout report](output/g2_readout_contract.md),
[code](code/verify_g2_readout_contract.py),
[noisy-tag report](output/g2_noisy_tag.md),
[tag code](code/verify_g2_noisy_tag.py),
[derivation](tex/route_g_readout_contract.tex),
[PDF](output/pdf/route_g_readout_contract.pdf).

Next physical handoff: choose a concrete energy scale/realization and
sensor interaction, or provide an actual apparatus's beam/response
calibration. Decide explicitly whether the target is overall error or
per-sign recall. Do not introduce another mathematical prerequisite
or finer scan in place of that missing physical choice.
Only an absolute-yield/finite-pulse/coherence question activates a
four-dimensional packet and pump-envelope specification.

### G2-P — Physical energy anchor and visible-probe candidate

Status: bounded-done as a scale/coupling audit (2026-09-23);
physical realization and measured calibration remain open.
Scale checks: 877/877; candidate probe/matching checks: 113/113.
Both were independently rerun. Existing G1/G2/G2-R/G2-S/G2-D/G2-M/tag
regressions also pass (2066/2066 combined). The four-page TeX/PDF was
compiled and every page visually reviewed. These are implementation
checks, not empirical validation or calibrated hardware qualification.

The current dimensionless forecast has a scale degeneracy. With
Lambda=1/R, energies and resolutions scale as Lambda, K as Lambda^-2,
mean work/event as Lambda, and normalized classification error is
unchanged. The 877/877 synthetic checks rescale every one of the 66
driven open channels and recover the same 2.20225% error. No eV/GeV
choice or measured calibration is inferred from this exercise.

Three externally calibrated consecutive signed-mode mass² values
determine A=Lambda², B=2 alpha Lambda², C=M_eff²+alpha² Lambda².
The fourth must satisfy Delta³(m_n²)=0 within specified covariance and
theory error. One mass only anchors Lambda if all dimensionless ratios
and its mode identification are already assumed. Sorted mass order is
not signed-mode order. A common positive mass² shift preserves this
closure test and cannot repair a wrong spacing pattern.

Candidate, not yet adopted: an SM-singlet tower with a uniform
`-H†H(lambda_Phi sum|phi_n|² + lambda_X sum|x_j|²)` interaction.
It is a KK-diagonal, charge-preserving 4D EFT. Writing one 4D Higgs
uniformly around the circle is not a completed local 5D SM embedding.
A brane-localized alternative is off diagonal, mixes masses and lets
the apparatus absorb compact momentum; old full-system selection rules
cannot be retained without accounting for that change.

The candidate audit derives common EWSB shifts, the finite-t free-nucleon
kernel, its NR limit, and the complete kinematically open complex-tower
Higgs width. Their shared coupling yields a conditional upper bound on
scattering from an applicable invisible-width budget. No cosmological
abundance, experimental limit value or instrument performance is imported.
At fixed bare masses the bound is implicit because thresholds change.

Two corrections prevent false promotion:

1. Keeping physical MD²=25 Lambda² with the retained action (no X
   self-quartic) requires lambda_X<=50 Lambda²/v². A negative underlying
   X mass runs away along H=Phi=0; adding a stabilizer changes the model.
2. If both portals are nonzero, Higgs exchange adds an elastic Phi-X
   amplitude that interferes with -g. Old inclusive fractions then need
   recomputation. An X-only coupling at a stated matching scale is a
   minimal candidate; a zero Phi portal is not radiatively protected.

The collision and readout verifiers are unchanged. These normalization
checks are synthetic, not evidence of physical realizability. A nucleon
cross section is not an instrument response: one recoil deposit retains
an unknown angle; nuclear response, efficiencies, geometry, flux and
mass-label confusion remain physical inputs. The drive is still a
prescribed reservoir, not an event-resolved work sensor.

Deliverables: [physical input card](PHYSICAL_INPUT_CARD.md),
[scale code/report](output/g2_scale_anchor.md),
[probe code/report](output/g2_probe_coupling.md),
[TeX](tex/route_g_physical_anchor.tex),
[PDF](output/pdf/route_g_physical_anchor.pdf).

Next: choose a neutral-sector physical scale and target to compare a
required interaction strength with stability and experimental constraints;
alternatively supply actual calibration data, or explicitly select an
analogue-mode platform with a more limited interpretation. Do not replace
this physical choice by a finer resolution scan or an unrequested UV
construction. If the required coupling is incompatible with the declared
constraints, reject that candidate. No Route-F promotion follows.

### G2-X — External scale/target choice and actual response fold

Status: bounded-done (2026-09-23); a conditional feasibility calculation,
not a completed production/detection apparatus.

Authorized choice: X-only circle-uniform Higgs coupling at a stated tree
matching prescription, physical intercept MD=5 Lambda, Lambda=1–30 GeV,
natural liquid xenon. External motivation is the Higgs mass threshold
Lambda_h=mh/10=12.508 GeV, not a fit of a separate radius to each particle.
Input mh125.08 GeV, v246 GeV, GammaSM4.10 MeV comes from the pinned PDG
review provenance; use ATLAS2023 observed B_inv<.107 at95%CL with SM
production/escape/drive-off assumptions. Sum all open complex X modes.

Actual response is prioritized: use the official XENON1T S2-only2019/2020
release, commit5a364bc8709f2561e5a013ddea6993a5a7c8e313. Keep raw CSVs,
hashes, license metadata and usage instructions. Restrict true NR to
.7–50 keV and full S2 bins29–198 (150.027–3000 PE). The exact raw search
exposure is356770 kg day, and the matrix already includes selection
losses. No extra fiducial fraction, fake calibration, halo density,
halo-WIMP exclusion or new Gaussian-resolution scan is introduced.

The coherent scalar-current/Helm kernel keeps finite incident kinematics
and natural-isotope weights. It is only used at admitted low momentum
transfer; unmodeled low/high recoil contributions are not assigned zero
physical response or included through uncontrolled total-cross-section
fractions. Source weights remain the rate-conditioned G2 driven mixture.

Map three distinct items: conditional Higgs exclusion; failure of the
retained action when lambda_X>50 Lambda²/v²; and the at-detector selected
particle flux required for three expected accepted-window recoils.
No actual flux is known. The third item is not discovery sensitivity,
95% exclusion, or an achievable beam. Response, form-factor, matching
and source uncertainties are not replaced by numerical quadrature error.

New readout result: p_X,min=.8091444644 Lambda gives a recoil endpoint
above the complete response window throughout the chosen range. Thus
d_sigma/d_ER=(lambda_X²/p_X²)H_A(ER); every common energy-only response
preserves identical normalized S2 shapes for j=+1 and j=-1. Accepted
priors are reweighted by1/p_X², but S2 provides no further sign information
within this kernel/window. Accepted priors are 85.943%/14.057%; the best
S2-only error is **14.057%**, merely the minority prior. The minimum
endpoint is 9.543 MeV. Monotonicity in Lambda bounds the whole interval,
not just numerical sample points. Reject transferring the old2.20% momentum
classifier claim to this sensor; it is not repaired by finer Gaussian fits.

The folded cross section is 7.99989e-39 cm² times
lambda_X²/(Lambda/GeV)². At Lambda=1, 5, 30 GeV the three-event
requirements at the conditional coupling ceiling are 5.07e9, 2.15e10,
4.32e6 cm⁻²s⁻¹. The minimum here optimizes only the declared coupling
domain, not production or full experimental sensitivity. The new
response/kinematics/tower verifier passes 573/573; all prior 3056
checks pass unchanged, total **3629/3629**. The figure uses analytic
formula rendering, not additional Gaussian scans. Data, formula,
implementation and plot hashes are preserved for reproduction.
The six-page derivation PDF and standalone figure were rendered and
visually reviewed; the final TeX build has no layout/reference warnings.

Physical source issue exposed by setting units: the activated drive
Omega=.2–6 GeV corresponds to frequencies about4.84e22–1.45e24 Hz.
A prescribed pump is not a realized supply or an ordinary oscillator.
The drive remains in the mathematical branch, but no hardware promotion
is made. X-only is not loop-protected; this is a tree feasibility test,
not a claim of completed quantum matching or certified UV control.

Next decision is now between two physical aims, not another resolution
scan: (1) event-count search, requiring production/at-detector flux and
background inference; (2) sign identification, requiring an independent
timing/momentum channel, direction or endpoint-sensitive target/readout.
If a real source cannot meet the required flux, reject that realization;
if sign separation is mandatory, this retained S2-only channel fails it.

Artifacts: [figure](output/figures/g2_xenon_feasibility.png),
[report](output/g2_xenon_feasibility.md),
[TeX](tex/route_g_xenon_feasibility.tex),
[PDF](output/pdf/route_g_xenon_feasibility.pdf),
[raw-response provenance](data/xenon1t_s2only/PROVENANCE.md).
Route-F physics and fits remain independent and unpromoted.

### G2-H — Physical source budget and production/readout compatibility

Status: bounded-done (2026-09-23). The on-shell-Higgs/single-pass-xenon
retained elastic NR realization is **rejected-for-this-model** at the
published Run-2 source budget. This is not a theory-wide no-go.

Source selection: ordinary pp→h+anything, h→X_j anti-X_j, no classical
drive and no assumption of the old selected Phi-X collision ensemble.
Use one actual13 TeV interaction-point dataset,139 fb^-1, and the
SM production input55.6 pb from ATLAS2023. The measured visible-channel
extraction assumes SM branching ratios; it is not treated as a separate
model-independent normalization after changing invisible decays.
The retained conditional B_inv<=.107 allowance gives
N_h=7.7284e6 and N_X+anti-X<=1.6538776e6. No actual dark flux is observed.

Source weights now follow B_j=Gamma_j/(GammaSM+sum Gamma), not the
old rate-conditioned collision priors. The source-to-fluence expression
includes the actual Higgs boost distribution, geometry/transport and
time profile as physical inputs. No full correlated distribution is
fabricated. The integrated budget must not be multiplied by the old
kg-day exposure again.

The decisive result bypasses unknown beam details without assuming them:
for every momentum, the retained-window scalar cross section obeys
sigma_acc<=lambda_X² fN² mN²/(4pi mh⁴)
sum_A xi_A A² mA²/(mX+mA)². This follows by bounding only the admitted
kernel, not by extrapolating coherent physics to all recoil energies.
Give all produced particles the maximum chord of the published xenon
cylinder and readout efficiency one. Use mX>=5Lambda, the existing bare
mass condition and monotonic Lambda⁴/(5Lambda+mA)². Across1–30 GeV,
**N_acc<3.96e-6**; an explicit generous input stress gives5.35e-6.
This is an upper envelope, not a simultaneously attainable threshold
point, an experimental confidence limit or a measured event rate.
Above12.508 GeV the on-shell source is closed; the abs(j)=1 source
already closes at12.2651 GeV.

New information gate: for direct scalar Higgs production and the
diagonal X-only probe, opposite signs have equal priors and identical
full four-dimensional kinematic distributions. Timing, direction and
independent momentum cannot distinguish them: I(j;z)=0 and optimal
error50%. This is source/channel-specific, not an all-orders symmetry
claim for the Phi-holonomy action. Do not transfer the old asymmetric
source's14.057% S2 prior error or2.20% ideal momentum error.

Implementation: new source/bound/endpoint/response/Lorentz checks389/389;
all prior3629 checks remain passing, total4018/4018. Actual NR response
folds are tested for slow particles whose endpoints cross bin edges;
the principal arbitrary-position bound correctly uses efficiency<=1.
No Gaussian optimization, beam Monte Carlo or added coupling is used.
The four-page PDF was compiled twice and every page visually reviewed;
the final build has no layout or unresolved-reference warnings.

Stop: do not refine this rejected chain's beam geometry or sensor noise.
Next options are (1) production/missing-momentum observables that avoid
a second rare interaction; (2) a physical sign-sensitive preparation or
coupling if signed conversion is required. On-shell invisible-Higgs
production measures only an inclusive width, not the tower or signs.
Off-shell sources, astrophysical populations, other targets, electronic/
inelastic and outside-window recoils, and recirculation are untested,
not quietly ruled out or appended to the same model.

Artifacts: [source card](data/HIGGS_SOURCE_CARD.md),
[report](output/g2_higgs_source.md),
[TeX](tex/route_g_higgs_source.tex),
[PDF](output/pdf/route_g_higgs_source.pdf).
Route F and its gates remain independent.

### G2-T — TOF/momentum metrology and signed-mode audit

Status: bounded-done (2026-09-23), **mathematical measurement conditions
only**. Actual source timing, independent X momentum/energy sensing,
directional calibration and signed preparation remain unspecified.
The action, source strength and prior response are unchanged.

The user selected two separate lanes:

1. Different masses: derive m²=p²[(cT/L)²−1], its full logarithmic
   covariance and exact bounded-error intervals. TOF alone fixes m/p;
   two timed hits still do not determine mass without another constraint.
   Unknown t0 and source-time/momentum correlations must enter the joint
   likelihood. No Gaussian or independent-source-time assumption is made.
2. Equal mass, opposite sign: retain the source-specific indistinguishability
   theorem. Direct Higgs production plus the diagonal portal remains
   sign-blind even with perfect timing, direction and momentum.
   The old asymmetric collision is a different preparation.

New practical obstruction: X is neutral; a charged-track momentum
resolution cannot be imported. A missing-momentum pair sum is not an
individual-particle momentum. Existing S2 response is not X TOF, and
electron drift time is not the source-to-interaction flight time.

A same-action alternative is derived, not promoted to hardware.
For known target mass M at rest, q²=ER(ER+2M) and
D=beta q cos(theta)−ER,
E=M ER/D, p=beta E, m=E sqrt(1−beta²).
It requires recoil direction, incident direction, beta, energy and a
target-isotope treatment. At Lambda=5 GeV, E=62.54 GeV, A=132, ER=10 keV,
the angular-only 1% mass-error budget is about0.896 arcsec.
This conditioning result, with other inputs exact, is not a measured
detector resolution or a recommendation to optimize xenon directionality.

The at-rest,10 m illustration gives a139.214 ps j=0/abs(j)=1 time gap.
With bounded dp/p<=.001 and dL<=1 mm, the exact reconstructed-mass boxes
first touch at a timing halfwidth48.672 ps;20 ps boxes are disjoint.
These are illustrative requirements, not sigmas or collider confidence
intervals. The +1/−1 boxes are identical.

Revalidated, not rediscovered: G2-R's undriven source correlation
j=n+l−m, or its full COM kinematic inversion for unobserved l and
m−n nonzero. A real Phi source/tag is still missing; driven events need
the work exchange in the inference. It is not an independent test of
the conservation law used to infer j.

New-interaction guardrail:
|Me+j Mo|²−|Me−j Mo|²=4j Re(Me* Mo).
A sign-odd amplitude by itself, a relative phase pi/2, or orthogonal
final states do not supply the required interference. Any new coupling
must exhibit an observable signed response and its physical reference.
No new interaction or external field is activated in this milestone.

The joint source/response measure includes arbitrary t0–p correlations.
Additional gates lie in[0,1], so the previous arbitrary-momentum,
efficiency-one retained-elastic-event ceiling3.9596863e-6 is unchanged.
No revived LHC-to-xenon experiment, new flux or Gaussian scan is claimed.

Stop/next: mass lane needs a real source tag plus independent momentum/
energy mechanism and a surviving event budget before response forecasts.
Sign lane first needs a physically accessible source correlation, or
an explicitly proposed interaction for separate audit. Do not refine
the rejected Higgs/xenon chain; collider production observables avoiding
a second tiny scattering probability remain a distinct useful lead.
Neither lane is a new Route-F prerequisite.

Artifacts: [input card](data/TOF_READOUT_CARD.md),
[report](output/g2_tof_mode_audit.md),
[TeX](tex/route_g_tof_mode_audit.tex),
[PDF](output/pdf/route_g_tof_mode_audit.pdf).
Verification:97/97 new checks and all4018 prior checks pass, total4115.
The five-page PDF was compiled twice, has no final layout/reference
warnings, and every page was visually reviewed. Algebra tests are not
empirical validation, and none of the synthetic requirements is calibration.

## G3 — Background selection

Status: open, not a prerequisite for the conditional G1 calculation.

Specify what fixes R and alpha, then check a stable stationary background
and whether the field content actually favors that holonomy. A potential
chosen solely to force a desired particle mass is not a derivation.
Do not assume a freely adjustable nonzero alpha is a selected vacuum.

## G4 — Particle-content and Route-F interface

Status: open; no identification presently authorized.

Before comparing with elementary-particle masses, derive spin, charges,
chirality and multiplicities from a single field/boundary content.
The present scalar KK tower cannot be relabeled as the SM fermions.
Any link to Route-F h_D, f_D, Higgs overlaps and Majorana scale must come
from shared normalized modes and one action, not independent fits.
Full string embedding is optional and separate, not a new mandatory gate.

## Explicit rejection conditions

| Test | What fails |
|---|---|
| A passive basis/gauge change alters poles after all sources are transformed | The implementation or proposed interpretation. |
| A free-circle tower requires mode-dependent R, M5 or alpha | This common-action ansatz; do not repair with a separate geometry per particle. |
| Claimed neighboring signed branches violate the common second-difference law | The free-circle mass assignment at its declared accuracy. |
| A pole disappears from both fixed probes despite their nonzero sum rule | The G1 result or the proposed disappearance claim. |
| An interaction needs energy creation or forbidden charge/KK-momentum changes | That transition model. |
| Norm positivity or bounded energy is lost | That proposed completion; negative mode number is not a negative-energy solution. |

G1's failure or success is not a proof or disproof of every possible
energy-spectrum model. A more general internal operator is a new
explicit hypothesis and must earn its extra parameters with predictions.
