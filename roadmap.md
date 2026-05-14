# Roadmap: Another Physics First-Principles Audit

Last updated: 2026-05-07

## Current status

The literal theory in `upload/Another physics.pdf` is not derivable as an unconditional
first-principles theory. The core PDF concepts are valuable as intuition, but several
literal claims conflict with standard first-principles constraints:

- Stable particle mass is the Poincare Casimir `P^mu P_mu = m^2 c^2`; it cannot
  periodically disappear with an internal phase while remaining the same elementary
  particle representation.
- Massive charged particles cannot reach light speed; an electron has `v<c` for every
  finite energy. Magnetic fields are the curvature of a `U(1)` gauge connection, not a
  consequence of electrons moving exactly at `c`.
- Event horizons forbid causal escape from inside a black hole. Jet-like ejection must be
  modeled outside the horizon, for example via disk/magnetosphere dynamics.
- A time-projection explanation of flat rotation curves with finite baryonic mass requires
  `f(r)^2 ~ 1/r`; if `f=dt_obs/dtau=1/sqrt(g_tt)`, this implies `g_tt ~ r`, violating
  asymptotic flatness for an isolated finite-mass source.
- Stable negative-energy or negative-mass elementary particles violate the positive-energy
  condition of relativistic quantum theory unless reinterpreted as holes, quasiparticles,
  tachyonic instabilities, or effective media excitations.

Therefore the unconditional result is a no-go theorem for the literal PDF theory.

## Salvageable concept

The viable direction is not to prove the original document literally, but to rebuild its
usable ideas as a conditional effective theory:

1. Replace `Ei` by a covariant information-current or order-parameter sector.
2. Replace `Ed` by a scalar density functional derived from stress-energy and/or an
   information metric.
3. Keep mass as a Poincare invariant; allow only medium-dependent effective masses or
   visibility weights to vary.
4. Treat time as metric/proper-time structure, not as an arbitrary scalar rescaling.
5. Treat electromagnetism as `U(1)` curvature and gravity as metric dynamics.
6. Treat entanglement as Hilbert-space correlation with no superluminal signalling.

## Candidate first-principles framework

Build an "Energy-Information Effective Unification" (EIEU) action:

```text
S = integral d^4x sqrt(-g) [
    M_Pl^2 R/2
  - 1/4 F^A_{mu nu} F_A^{mu nu}
  + L_matter(psi, Phi, A, g)
  - 1/2 (nabla phi)^2 - V(phi)
  + lambda_1 phi T^mu_mu / M
  + lambda_2 (nabla phi)^2 T^mu_mu / M^4
  + lambda_3 phi F^A_{mu nu} F_A^{mu nu} / M
]
```

Here `phi` is the repaired `Ed` carrier. `Ei` is not a new form of mass; it is encoded in
state counting, entropy, correlations, or effective couplings derived from this action.

## Required theorem targets

- No-go theorem: prove the literal PDF claims fail under Poincare invariance, positive
  energy, gauge invariance, and asymptotically flat GR.
- Conditional theorem: prove that the EIEU action is generally covariant, gauge invariant,
  locally causal, and unitary below its EFT cutoff when the kinetic terms have healthy signs.
- Observable map: derive corrections to
  - gravitational redshift and fifth-force tests,
  - galaxy rotation/lensing consistency,
  - particle mass/width visibility shifts,
  - atomic spectroscopy and equivalence-principle constraints.

## Verification plan

1. Algebraic checks:
   - massive-particle speed bound,
   - flat-curve time-projection asymptotic conflict,
   - positive-energy stability bound.
2. Numerical smoke tests:
   - compute `beta=v/c` for finite-energy electrons;
   - compute required `g_tt` for `f^2~1/r`;
   - show that flat rotation via pure time projection violates asymptotic flatness.
3. Model-building tests:
   - implement the EIEU weak-field limit;
   - compare rotation curve and lensing predictions together;
   - reject any parameter point that fits velocities but fails lensing/redshift bounds.

## Definition of done

The original PDF is closed as a literal theory: refuted.

The research program can continue only as a new conditional EFT/GUT-inspired framework
with explicit action, symmetries, conservation laws, cutoff, and falsifiable observables.

## GUT-from-PSLT audit update

Last updated: 2026-05-05

User request: derive a grand-unification theory from first principles, starting from PSLT
and Another Physics concepts, with unconditional proof/refutation and numerical checks.

### Current theorem status

The unconditional result is again a no-go for a literal derivation:

- PSLT, as written in `upload/main 3.pdf`, explicitly claims an EFT-level conditional
  closed chain, not a full action-level theorem proving the Standard Model family count.
- Another Physics, as written in `upload/Another physics.pdf`, contains useful intuitions
  (`Ei`, `Ed`, spectral visibility, time-density language) but its literal claims conflict
  with Poincare invariance, positive energy, gauge invariance, event-horizon causality,
  and asymptotically flat gravity.
- A true GUT must at minimum specify a compact gauge group, matter representations,
  anomaly cancellation, symmetry breaking, fermion masses, neutrino sector, proton-decay
  operators, gravity interface, and low-energy matching. Neither PSLT nor Another Physics
  alone supplies these as unconditional consequences.

Therefore the strict theorem is:

```text
No local, unitary, Lorentz-invariant, positive-energy quantum field theory can be
unconditionally derived from the literal PSLT+Another-Physics manuscripts as a unique
GUT unless extra action-level axioms are added.
```

### Salvage route: Spectral-Layer GUT EFT

The viable creative direction is a conditional framework:

1. Use PSLT only for a spectral-layer selection functional
   `P_N proportional to B_N g_N (1 - exp(-Gamma_N t_coh))`.
2. Use Another Physics only as a covariant information-density sector:
   `Ed -> phi` and `Ei -> conserved/state-counting information current`.
3. Choose a real GUT gauge group, preferably `SO(10)` for one-generation anomaly-free
   `16` representations and a natural right-handed neutrino.
4. Let PSLT explain why exactly three low-lying chiral layer copies are occupied, while
   the GUT explains gauge charges inside each copy.
5. Couple the information field only through gauge-invariant operators, never by making
   particle rest mass periodically disappear.

Candidate action:

```text
S = int d^4x sqrt(-g) [
    M_Pl^2 R/2
  - 1/4 Tr F_{mu nu} F^{mu nu}
  + sum_N psi_Nbar i gamma^mu D_mu psi_N
  + |D H|^2 - V_GUT(H, Sigma, phi)
  - 1/2 (partial phi)^2 - V_phi(phi)
  + sum_N y_N(phi) psi_N psi_N H
  + L_PSLT[g_N(phi), Gamma_N(phi,D,eta), B_N(phi)]
  + higher-dimensional operators / Lambda
]
```

The central conditional claim to test:

```text
SO(10) fixes anomaly-free gauge unification per occupied layer.
PSLT fixes the number of occupied layers.
Another-Physics/Ed becomes the scalar/information control variable for the occupancy
and visibility functions.
```

### Next verification steps

1. Prove anomaly cancellation for one `SO(10)` spinor `16`.
2. Prove that three copies remain anomaly free.
3. Prove PSLT occupancy is normalizable and gives a finite `R3`.
4. Numerically scan a simple PSLT weight model and show when `R3 > 0.9`.
5. State what would falsify the model:
   - no stable first-three-layer region,
   - proton decay above/below experimental consistency once operators are specified,
   - gauge coupling unification failure after threshold assumptions are fixed,
   - fifth-force or equivalence-principle violation from `phi`.

## Novel element set: topological spectral-layer GUT

Last updated: 2026-05-05

To make the PSLT-to-GUT path less phenomenological, add a topological protection layer:

1. **Topological generation index.** Put the internal spectral layer problem on a compact
   two-center deformation of `CP^1` and twist the internal Dirac/Dolbeault operator by
   the line bundle `O(2)`. Riemann-Roch gives
   `dim H^0(CP^1,O(2)) - dim H^1(CP^1,O(2)) = 3`, and positivity gives `H^1=0`.
   This protects three chiral zero modes under smooth two-center deformations.
2. **PSLT as visibility/gap dynamics, not the source of chirality.** PSLT's
   `W_N = B_N g_N (1-exp(-Gamma_N t))` suppresses higher visible layers and sets
   the occupancy of excited modes. The generation count is topological; the finite
   fourth-layer visibility is dynamical and falsifiable.
3. **Another-Physics `Ed` as a dilaton/information modulus.** Replace time-density
   literalism by a scalar `phi` that controls gauge kinetic functions and spectral geometry:
   `f_G(phi) Tr F^2`, `Omega(x;phi,D,eta)`, and PSLT parameters.
4. **Gauge unification from a critical information point.** At the critical phase-locked
   point `phi=phi_*`, impose `f_1(phi_*)=f_2(phi_*)=f_3(phi_*)`; away from it, threshold
   corrections are calculable as derivatives of `f_a(phi)`.
5. **Mass hierarchy as overlap integrals.** Yukawa matrices become
   `Y_ij = y0 int_X s_i(x) s_j(x) h(x;phi)`, where `s_i` are the three protected
   holomorphic sections of `O(2)`. Moving the two centers changes these overlaps
   exponentially, giving a geometric route to flavor hierarchy.

Verification added in this pass:

- Riemann-Roch/monomial count: `O(2)` on `CP^1` has exactly 3 holomorphic sections.
- SM anomaly sums for one `SO(10)`-compatible generation cancel exactly.
- PSLT smoke scan has broad parameter regions with balanced `P1,P2,P3` and small `P4`.
- One-loop gauge running shows non-SUSY SM does not exactly unify, while supersymmetric
  beta functions are much closer; threshold/dilaton corrections are therefore mandatory
  rather than optional.

## Paper-grade GUT development plan

Last updated: 2026-05-06 20:11 CST

Goal: turn the conditional `Spin(10) + topological spectral-layer + PSLT visibility`
framework into a paper-level theory. The next stage must stop being only a conceptual
proposal and become a chain of derivations with falsifiable numerical checks.

Recommended execution order:

1. **Breaking chain and hypercharge theorem.**
   Use the Pati-Salam route as the primary chain:
   `Spin(10) -> SU(4)_C x SU(2)_L x SU(2)_R -> SU(3)_C x SU(2)_L x U(1)_Y`,
   with
   `Y = T_{3R} + (B-L)/2`.
   This is preferred over the direct `SU(5) x U(1)_chi` chain because it keeps
   `B-L`, the right-handed neutrino, seesaw physics, and dangerous baryon-number
   operators visible at every intermediate step.

   Required derivation:
   `16 -> (4,2,1) + (bar 4,1,2)`, then decompose `4 -> (3,1/3)+(1,-1)`
   and `bar 4 -> (bar 3,-1/3)+(1,+1)` to recover exactly
   `Q, u^c, d^c, L, e^c, nu^c` with Standard Model hypercharges.

   Numerical/algebraic verification:
   write a table-driven checker for all states and verify charge assignments,
   dimensions, and anomaly sums.

2. **Yukawa texture from protected sections.**
   Use the three protected `O(2)` sections on `CP^1` as family wavefunctions.
   In a normalized basis, define
   `Y_ij^(a) = lambda_a int_X s_i(x) s_j(x) h_a(x;phi) dmu_phi(x)`,
   where `a = u,d,e,nu`, and `h_a` are Higgs/localization profiles controlled by the
   information modulus `phi`.

   First tractable toy model:
   approximate the normalized section basis by Bernstein-type functions on `[0,1]`,
   for example `s_0=(1-x)^2`, `s_1=sqrt(2)x(1-x)`, `s_2=x^2`, and use one or two
   Gaussian Higgs centers. The singular values of the overlap matrix become the
   mass hierarchy.

   Required verification:
   scan center position, width, and profile mixing; record whether the model can
   naturally produce hierarchical singular values without fine tuning. A paper-grade
   criterion is stability under small deformations of the centers and widths.

3. **Neutrino seesaw sector.**
   Because `Spin(10)` contains `nu^c` inside each `16`, add a `126_H` or equivalent
   `B-L` breaking operator to generate
   `M_R = y_R v_R`.
   The low-energy light-neutrino matrix is
   `m_nu = -m_D^T M_R^{-1} m_D`.

   Required derivation:
   show which representation or effective operator gives the Majorana mass while
   preserving the breaking chain. Then prove that the same topological family basis
   used for charged Yukawas also gives the neutrino Dirac texture.

   Numerical verification:
   infer the required right-handed neutrino scale from
   `M_R ~ m_D^2 / m_nu`. For example, if `m_D ~ 100 GeV` and `m_nu ~ 0.05 eV`,
   then `M_R ~ 2 x 10^14 GeV`. Scan whether this scale is compatible with the
   intermediate `B-L` breaking scale from the chosen chain.

4. **Proton-decay operators.**
   Separate unavoidable dimension-6 gauge-boson operators from model-dependent
   dimension-5 colored-Higgs operators.

   Gauge-boson operator template:
   `O_6 ~ (g_G^2/M_X^2) (q q q l)`,
   with rough lifetime scaling
   `tau_p ~ M_X^4 / (alpha_G^2 m_p^5)`.

   Model-dependent dangerous operators:
   if a supersymmetric or light-colored-triplet completion is used, include
   dimension-5 structures such as `QQQL/M_T` and `u^c u^c d^c e^c/M_T`.
   The theory must either suppress these by symmetry/localization or explicitly
   show that the triplet threshold is high enough.

   Required verification:
   from the same running and threshold fit that determines `M_G`, compute the
   dimension-6 lifetime estimate and flag whether the model is excluded, marginal,
   or viable after experimental comparison. Then scan the colored-triplet scale
   needed to suppress dimension-5 decay if that sector is present.

5. **Threshold corrections and matching equations.**
   Use one-loop matching as the first paper-grade threshold layer:
   `alpha_i^{-1}(M_Z) = alpha_G^{-1} + (b_i/2pi) log(M_G/M_Z) + Delta_i`,
   with
   `Delta_i = (1/2pi) sum_r b_i^r log(M_G/M_r) + delta_i^phi`.
   Here `delta_i^phi` is the PSLT/Another-Physics contribution through gauge kinetic
   functions `f_i(phi)`.

   Required derivation:
   solve only two independent differences,
   `alpha_1^{-1}-alpha_2^{-1}` and `alpha_2^{-1}-alpha_3^{-1}`, because the common
   shift is absorbed into `alpha_G`. Then express the required threshold vector
   as constraints on heavy multiplet mass splittings and `phi` derivatives.

   Numerical verification:
   compute the required `Delta_i-Delta_j` for the chosen chain; compare them to
   natural threshold sizes `O((1/2pi) log(M_G/M_r))`. If the required correction
   needs extreme mass ratios, the chain is disfavored.

6. **Paper falsifiability table.**
   Each section must end with a failure criterion:
   - breaking chain fails if hypercharges or anomalies do not match exactly;
   - Yukawa texture fails if hierarchy needs unstable fine tuning;
   - seesaw fails if `M_R` is incompatible with the `B-L` threshold;
   - proton decay fails if unavoidable operators are too large;
   - threshold corrections fail if unification requires unnatural splittings;
   - PSLT layer dynamics fails if three-family dominance is not robust.

Immediate next task:

Start with item 1. Build the complete Pati-Salam decomposition table for one `16`,
prove the hypercharge formula state by state, and implement a small algebraic checker.
Only after this succeeds should we fit Yukawa textures, because the texture labels
must correspond to the correctly embedded Standard Model fields.

## Item 1 result: Pati-Salam breaking chain verified

Last updated: 2026-05-05 21:33 CST

Script:
`code/verify_pati_salam_spin10.py`

Outputs:

- `output/pati_salam/pati_salam_report.md`
- `output/pati_salam/pati_salam_decomposition.csv`
- `output/pati_salam/pati_salam_components.csv`
- `output/pati_salam/pati_salam_verification.json`

Verified branching:

```text
Spin(10) 16 -> (4,2,1) + (bar4,1,2)
SU(4)_C -> SU(3)_C x U(1)_{B-L}
SU(2)_R x U(1)_{B-L} -> U(1)_Y
Y = T3R + (B-L)/2
```

Recovered Standard Model left-handed Weyl fields:

```text
Q    : (3,2)_{+1/6}
L    : (1,2)_{-1/2}
u^c  : (bar3,1)_{-2/3}
d^c  : (bar3,1)_{+1/3}
nu^c : (1,1)_0
e^c  : (1,1)_{+1}
```

Exact verification results:

- total left-handed Weyl dimension is `16`;
- all state-by-state hypercharges satisfy `Y = T3R + (B-L)/2`;
- `SU(4)^3` Pati-Salam anomaly cancels between `(4,2,1)` and `(bar4,1,2)`;
- both `SU(2)_L` and `SU(2)_R` have four doublets, so the Witten global anomaly is absent;
- Standard Model anomaly sums are exactly zero:
  `SU(3)^3`, `SU(3)^2-U(1)_Y`, `SU(2)^2-U(1)_Y`, `U(1)_Y^3`, and gravitational-`U(1)_Y`.

Nontrivial elements extracted for the paper:

1. **Hypercharge uniqueness lemma.**
   If `Y = a T3R + b(B-L)/2`, then the observed charges of `Q` and `u^c` force
   `a=b=1`. The Pati-Salam embedding therefore fixes the hypercharge direction
   rather than leaving it as a tunable phenomenological choice.
2. **Threshold-ready normalization.**
   On the full `16`,
   `Tr(T3L^2)=2`, `Tr(Y^2)=10/3`, hence
   `k_Y = Tr(Y^2)/Tr(T3L^2)=5/3`.
   This gives the GUT-normalized coupling `alpha_1=(5/3)alpha_Y`, which must be used
   in the threshold-correction section.
3. **Orthogonal broken Abelian direction.**
   The combination
   `X_PS = 2T3R - 3(B-L)/2`
   satisfies `Tr(Y X_PS)=0`. This gives a clean basis for the heavy neutral gauge
   boson and for later threshold corrections controlled by the `Ed -> phi` modulus.
4. **Layer-factorization check.**
   Since one `16` is anomaly free, any PSLT/topological layer multiplicity that gives
   complete `16` copies preserves anomaly cancellation automatically. Therefore PSLT
   may control the number and visibility of families, but it must not split incomplete
   chiral fragments below the GUT threshold unless a new anomaly-canceling mechanism
   is supplied.

Conclusion for item 1:

The Pati-Salam chain is internally consistent and is now the recommended breaking
chain for the paper. The next item should be the Yukawa texture module, using the
verified field labels above and the `O(2)` family basis.

## Item 2 result: CP1 O(2) Yukawa texture verified as a hierarchy mechanism

Last updated: 2026-05-05 21:40 CST

Script:
`code/scan_cp1_o2_yukawa.py`

Outputs:

- `output/yukawa_o2/cp1_o2_yukawa_report.md`
- `output/yukawa_o2/cp1_o2_yukawa_candidates.csv`
- `output/yukawa_o2/cp1_o2_yukawa_scan.json`

### Mathematical construction

Use the topological family space

```text
H^0(CP1,O(2)) = span{s_0,s_1,s_2}
```

with normalized wavefunctions under the Fubini-Study measure
`dmu_FS = sin(theta)dtheta dphi/(4pi)`:

```text
psi_0 = sqrt(3) cos^2(theta/2)
psi_1 = sqrt(6) exp(i phi) sin(theta/2) cos(theta/2)
psi_2 = sqrt(3) exp(2 i phi) sin^2(theta/2)
```

The exact normalization follows from

```text
int_CP1 z^m bar{z}^n (1+|z|^2)^(-2) dmu_FS = delta_mn m!(2-m)!/3!
```

for `m,n=0,1,2`, hence the coefficients `sqrt(3)`, `sqrt(6)`, `sqrt(3)`.

For each sector
`a in {u,d,e,nu_D}`, define a symmetric Spin(10)-type Yukawa matrix

```text
Y_ij^(a) = lambda_a int_CP1 psi_i(theta,phi) psi_j(theta,phi)
           h_a(theta,phi;phi_Ed) dmu_FS.
```

The first toy Higgs kernel is a PSLT/Ed two-center profile:

```text
h_a = K_1 + A exp(i chi) K_2
      + eps exp(i chi_eps) + q exp(i chi_q) P_2(cos theta),
K_A = exp[kappa_A (n dot n_A - 1)].
```

Interpretation:

- `K_1,K_2` are two phase-locked localization centers on the internal `CP1`;
- `A exp(i chi)` is a relative Berry/Ed phase and amplitude;
- `eps` is an isotropic leakage channel;
- `q P_2(cos theta)` is the leading curvature/quadrupole deformation;
- different `Spin(10)` Higgs sectors correspond to different kernels `h_u,h_d,h_e,h_nu`.

### Nontrivial elements for the paper

1. **Holomorphic moment-matrix principle.**
   Yukawa matrices are not inserted as arbitrary `3x3` matrices; they are degree-two
   holomorphic moment matrices of a scalar kernel on `CP1`.
2. **Rank staircase mechanism.**
   In the point-localized limit, one center gives
   `Y_ij proportional psi_i(x_0) psi_j(x_0)`, hence rank one. Two centers give rank
   at most two. The smallest singular value is generated by finite width plus
   curvature/leakage deformation. This naturally mirrors third, second, first
   generation ordering.
3. **Shared family topology with sector misalignment.**
   All sectors use the same protected `O(2)` family basis, while different Higgs
   kernels produce different symmetric textures. This gives a controlled origin for
   flavor hierarchy and later CKM/PMNS mixing.
4. **CP source from internal two-center phase.**
   The relative phase between two kernels makes `Y` complex symmetric. A later
   `120_H` antisymmetric overlap can be added if CKM/PMNS fitting demands it, but
   it is not needed for the first hierarchy existence test.

### Numerical verification

The scan first imports the item-1 field table and confirms that
`Q,u^c,d^c,L,e^c,nu^c` exist. It then performs exact quadrature on a
`54 x 108 = 5832` point grid. The orthonormality error of the `O(2)` basis is

```text
max |int psi_i^* psi_j dmu - delta_ij| = 1.998e-15.
```

Using `9000` random kernels per sector and stability tests under small geometric
perturbations, representative normalized singular values `[small, mid, large]` are:

```text
up              target [1.0e-5, 7.0e-3, 1]  found [1.233e-5, 6.888e-3, 1]
down            target [1.0e-3, 2.0e-2, 1]  found [9.307e-4, 1.892e-2, 1]
charged lepton  target [3.0e-4, 6.0e-2, 1]  found [2.897e-4, 6.385e-2, 1]
neutrino Dirac  target [1.0e-2, 2.0e-1, 1]  found [9.782e-3, 1.832e-1, 1]
```

Perturbation stability, measured as log10 interquartile width of the small and
middle singular ratios:

```text
up              (0.071, 0.075)
down            (0.058, 0.039)
charged lepton  (0.045, 0.075)
neutrino Dirac  (0.053, 0.023)
```

These widths are well below one order of magnitude, so the examples are not
single-point needle tunings.

Rank-staircase checks:

```text
up:
  single center [1.018e-8, 1.726e-4, 1]
  two centers   [1.205e-5, 6.859e-3, 1]
  full kernel   [1.233e-5, 6.888e-3, 1]

down:
  single center [1.396e-5, 5.022e-3, 1]
  two centers   [1.024e-4, 1.376e-2, 1]
  full kernel   [9.307e-4, 1.892e-2, 1]

charged lepton:
  single center [2.085e-8, 2.578e-4, 1]
  two centers   [2.278e-4, 6.387e-2, 1]
  full kernel   [2.897e-4, 6.385e-2, 1]
```

Conclusion for item 2:

The `CP1/O(2)` spectral-layer family basis can generate realistic-order hierarchical
Yukawa singular values through geometric two-center overlap kernels. This is not yet
a full flavor fit: the next flavor task is to fit CKM/PMNS mixing and decide whether
an antisymmetric `120_H` contribution is needed. For the main roadmap, item 2 is
sufficiently verified to proceed to item 3, the neutrino seesaw sector.

## Item 3 result: type-I seesaw reconstructed and geometrized

Last updated: 2026-05-05 21:46 CST

Script:
`code/verify_seesaw_item3.py`

Outputs:

- `output/seesaw/seesaw_item3_report.md`
- `output/seesaw/seesaw_reconstruction.json`
- `output/seesaw/majorana_veronese_kernel.csv`
- `output/seesaw/majorana_trace_lift_kernel.csv`

### Representation-theoretic derivation

From item 1, the field `nu^c` lies in the `Spin(10)` spinor `16`, specifically in
`(bar4,1,2)` after Pati-Salam branching. Its `B-L` charge is `+1`, so the bilinear
`nu^c_i nu^c_j` carries `B-L=+2`. Therefore a Majorana mass requires a `B-L=-2`
order parameter. In ordinary `Spin(10)` language this is supplied by a
`bar{126}_H` channel, or by the effective operator

```text
(16_i 16_j bar{16}_H bar{16}_H)/Lambda.
```

The neutral-fermion mass matrix is

```text
M_N = [ 0      m_D ]
      [ m_D^T  M_R ].
```

For `||m_D M_R^{-1}|| << 1`, block diagonalization by Schur complement gives

```text
m_light = -m_D M_R^{-1} m_D^T + O(m_D^4/M_R^3).
```

Conversely, for nonsingular `m_D` and `m_light`, the required Majorana matrix is
uniquely

```text
M_R = -m_D^T m_light^{-1} m_D.
```

This inverse-seesaw reconstruction was used with the item-2 `nu_D` texture and the
item-2 charged-lepton left rotation.

### Nontrivial elements for the paper

1. **Schur-complement theorem for the `CP1/O(2)` family basis.**
   The same protected three-family basis used for Dirac textures also supports the
   Majorana sector. The light-neutrino matrix is the Schur complement of the heavy
   `B-L=2` block.
2. **Inverse Majorana reconstruction.**
   Once a benchmark `m_light` is chosen, `M_R` is not arbitrary; it is fixed by
   `M_R = -m_D^T m_light^{-1} m_D`. This gives an immediate falsifiability test:
   the reconstructed `M_R` must be realizable by the allowed internal geometry and
   must have a plausible intermediate/GUT scale.
3. **Scalar Veronese no-go.**
   A pure scalar `CP1/O(2)` Majorana kernel cannot represent a generic complex
   symmetric `M_R`, because
   `psi_1^2 = 2 psi_0 psi_2`. Thus scalar moment matrices obey
   `(M_scalar)_11 = 2(M_scalar)_02`. The reconstructed `M_R` violates this relation.
4. **Majorana trace-lift channel.**
   The missing component is the spin-0 part in
   `Sym^2 H^0(O(2)) = spin-2 plus spin-0`. Add the Majorana-only contact/curvature
   matrix

   ```text
   C_0 = [[0,0,1],[0,-1,0],[1,0,0]]
   ```

   and represent

   ```text
   M_R/M_scale = sum_{a=1}^5 w_a v(x_a)v(x_a)^T + zeta C_0,
   v(x) = (psi_0(x),psi_1(x),psi_2(x)).
   ```

   This is a concrete new theoretical element: Dirac Yukawas can remain scalar
   moment matrices, while the `B-L=2` Majorana sector requires a curvature/contact
   trace-lift.

### Numerical verification

No web lookup was used. The benchmark was an approximate normal-ordering target:

```text
m1 = 1.0e-3 eV
Delta m21^2 = 7.42e-5 eV^2
Delta m31^2 = 2.517e-3 eV^2
sin^2(theta12,theta13,theta23) = (0.3040,0.0222,0.5730)
m_D largest singular value = 100 GeV
```

The item-2 `nu_D` singular ratios were

```text
[1, 0.1832365, 0.0097818].
```

Reconstructed light spectrum and PMNS angles after seesaw:

```text
m1,m2,m3 = (1.000000e-03, 8.671793e-03, 5.017968e-02) eV
sin^2(theta12,theta13,theta23) = (0.304000,0.022200,0.573000)
```

Heavy Majorana Takagi masses:

```text
M1,M2,M3 = (2.599776e10, 9.144789e13, 3.105386e15) GeV
```

Seesaw expansion parameter and residuals:

```text
||Theta||_2 = ||m_D M_R^-1||_2 = 4.393070e-11
seesaw matrix residual = 6.162603e-13
light mass residual = 5.818406e-13
```

Geometric kernel tests:

```text
pure scalar Veronese residual = 4.869104e-02
trace-lift residual = 2.627736e-16
trace-lift basis condition = 8.655981
```

Conclusion for item 3:

The item-2 `nu_D` texture can be connected to a mathematically consistent type-I
seesaw. The reconstructed heavy scale is plausible for a GUT/intermediate
`B-L` sector, and the trace-lift construction gives an exact `CP1/O(2)` geometric
realization of the required `M_R`. The result is still a reconstruction rather than
a prediction; the next stricter neutrino test is to restrict `M_R` to a small
number of PSLT/Ed centers and scan whether the same mass/mixing target remains
reachable without inverse-fitting. For the main roadmap, item 3 is verified enough
to proceed to item 4, proton-decay operators.

## Item 4 result: proton-decay operators separated and bounded

Last updated: 2026-05-05 21:53 CST

Scripts and reports:

- `code/verify_proton_decay_item4.py`
- `output/proton_decay/proton_decay_item4_report.md`
- `output/proton_decay/proton_decay_verification.json`
- `output/proton_decay/dimension6_gauge_lifetimes.csv`
- `output/proton_decay/dimension5_triplet_lifetimes.csv`

TeX draft:

- `paper/gut_framework.tex`
- `paper/gut_framework.pdf`

The TeX draft now records items 1 through 4 in one coherent manuscript skeleton,
so later threshold-correction work can extend the same file without losing the
previous derivations.

### Operator ledger

The two core baryon-violating four-fermion structures are

```text
O_QQQL = epsilon_{abc} (Q^a Q^b)(Q^c L),
O_UUDE = epsilon_{abc} (u^{c,a} u^{c,b})(d^{c,c} e^c),
```

with appropriate Lorentz and `SU(2)_L` contractions. The exact bookkeeping check gives:

```text
O_QQQL: Y=0, B=+1, L=+1, B-L=0, dim=6
O_UUDE: Y=0, B=-1, L=-1, B-L=0, dim=6
```

Thus the dangerous proton-decay operators preserve `B-L`. They are independent
from the item-3 `B-L=2` Majorana trace-lift, which generates neutrino mass but
does not itself force proton decay.

### Dimension-6 gauge exchange

For a heavy gauge boson `X`,

```text
L_eff = - g_G^2 J_X^dagger J_X / M_X^2 + O(M_X^-4),
C6 ~ g_G^2 / M_X^2.
```

The estimate used

```text
Gamma(p -> e+ pi0) ~= K |C6|^2,
K = m_p beta_H^2 (1+D+F)^2 A_R^2 / (64 pi f_pi^2)
  = 1.180182e-03 GeV^5.
```

Representative lifetimes:

```text
alpha_G^-1=24, M_X=1.0e15 GeV -> tau=6.446e31 yr
alpha_G^-1=24, M_X=3.0e15 GeV -> tau=5.222e33 yr
alpha_G^-1=24, M_X=1.0e16 GeV -> tau=6.446e35 yr
alpha_G^-1=24, M_X=2.0e16 GeV -> tau=1.031e37 yr
alpha_G^-1=40, M_X=3.0e15 GeV -> tau=1.450e34 yr
```

For the replaceable benchmark `tau > 1e34 yr`, the required masses are:

```text
alpha_G^-1=24: M_X > 3.529e15 GeV
alpha_G^-1=30: M_X > 3.157e15 GeV
alpha_G^-1=40: M_X > 2.734e15 GeV
```

This is the hard constraint for item 5: threshold corrections cannot lower the
broken gauge-boson mass below this window without reviving dimension-6 decay.

### Dimension-5 colored-triplet operators

In supersymmetric or holomorphic completions with color triplets,

```text
W5 = (C_L/M_T) QQQL + (C_R/M_T) u^c e^c u^c d^c,
C6_dressed ~= (alpha_2/4pi)(m_wino/m_sfermion^2)(S_T y_a y_b/M_T).
```

Here `S_T` is the new proposed `CP1/O(2)` triplet geometric filter:

```text
S_T^{ijkl} =
int psi_i psi_j psi_k psi_l h_T dmu
/ sqrt(int |psi_i psi_j psi_k psi_l|^2 dmu int |h_T|^2 dmu).
```

This filter is the item-4 nontrivial element: doublet kernels may generate
observed Yukawas while the colored-triplet kernel is displaced in the internal
geometry, suppressing first/third generation feed-through.

Using `M_T=1e16 GeV`, `m_wino=1e3 GeV`, `m_sfermion=1e5 GeV`, and the item-2
singular values:

```text
first-generation staircase: yprod=1.652e-10, S_T=1      -> tau=6.388e41 yr
second-generation leakage:  yprod=1.876e-06, S_T=1      -> tau=4.954e33 yr
third-gen unsuppressed:     yprod=1.440e-02, S_T=1      -> tau=8.412e25 yr
third-gen filtered:         yprod=1.440e-02, S_T=1e-5   -> tau=8.412e35 yr
```

For the same benchmark, unsuppressed third-generation feed-through requires

```text
S_T <= 9.172e-5      for 100 TeV sfermions,
S_T <= 9.172e-7      for 10 TeV sfermions.
```

Conclusion for item 4:

Dimension-6 proton decay is controlled by the heavy gauge threshold `M_X` and
cannot be removed by Yukawa geometry. Dimension-5 proton decay is model-dependent:
it is absent if there are no dangerous light triplets, and otherwise requires
either heavy superpartners or a strong `CP1/O(2)` triplet geometric filter. Item 4
is verified enough to proceed to item 5, threshold corrections, where unification
must be solved under the proton lower bound on `M_X`.

## Item 5 result: threshold corrections with proton-safe matching

Last updated: 2026-05-05 21:59 CST

Script and outputs:

- `code/scan_thresholds_item5.py`
- `output/thresholds/threshold_item5_report.md`
- `output/thresholds/threshold_summary.json`
- `output/thresholds/threshold_scan.csv`
- updated `paper/gut_framework.tex`
- updated `paper/gut_framework.pdf`

### Matching theorem

Use the one-loop convention

```text
alpha_i^-1(MZ) = alpha_G^-1 + b_i/(2pi) log(M_G/MZ) + Delta_i.
```

Only threshold differences are observable at this order. Fix the common
unobservable shift by

```text
sum_i Delta_i = 0.
```

For a fixed `M_G`, define

```text
u_i(M_G) = alpha_i^-1(MZ) - b_i/(2pi) log(M_G/MZ).
```

Then the required threshold vector is uniquely

```text
alpha_G^-1 = (u_1+u_2+u_3)/3,
Delta_i = u_i - alpha_G^-1.
```

This turns threshold matching into a two-dimensional vector problem in the
`Delta` difference plane.

### Proton-safe threshold cone

From item 4,

```text
M_X >= M_X^min(alpha_G)
```

where `M_X^min` is fixed by the benchmark `tau(p -> e+ pi0) > 1e34 yr`.
For every threshold point define

```text
rho_X^min = max(1, M_X^min(alpha_G)/M_G).
```

If `rho_X^min >> 1`, then the gauge boson responsible for dimension-6 proton
decay must be far above the nominal matching scale. Such a point may still be
mathematically matchable, but it is not a natural single-scale GUT threshold.
This is the item-5 nontrivial element: unification points must live inside the
proton-safe cone, not merely solve the gauge-coupling equations.

### Numerical inputs

No web lookup was used. Local benchmark values:

```text
alpha_em^-1(MZ) = 127.955
sin^2 theta_W(MZ) = 0.23122
alpha_s(MZ) = 0.1184
alpha_1 = (5/3) alpha_Y
alpha_i^-1(MZ) = (59.021547, 29.585755, 8.445946)
```

### SM one-loop baseline

Beta coefficients:

```text
b_SM = (41/10, -19/6, -7).
```

Pairwise crossings without thresholds:

```text
M12 = 1.031714e13 GeV
M23 = 1.019243e17 GeV
M13 = 2.472388e14 GeV
```

Best point with `M_X=M_G` and proton safety:

```text
M_G = 2.660725e15 GeV
alpha_G^-1 = 42.329763
Delta = (-3.539718, +2.881948, +0.657770)
||Delta||_2 = 4.611712
M_X^min = 2.657381e15 GeV
rho_X^min = 1
```

An illustrative threshold fit using only scalar fragments
`H_C=(3,1,-1/3)`, `Sigma_3=(1,3,0)`, `Sigma_8=(8,1,0)` exactly solves the
linear equations but requires

```text
max |log(M_G/M_r)| = 333.6105.
```

This is physically absurd. Therefore the purely SM one-loop branch is disfavored
unless a concrete intermediate spectrum or a controlled `E_d` gauge-kinetic
threshold supplies the large missing difference vector.

### MSSM-like one-loop baseline

Beta coefficients:

```text
b_MSSM = (33/5, 1, -3).
```

Pairwise crossings:

```text
M12 = 2.010611e16 GeV
M23 = 2.405823e16 GeV
M13 = 2.166712e16 GeV
```

Best proton-safe point:

```text
M_G = 2.150305e16 GeV
alpha_G^-1 = 24.274890
Delta = (-0.016084, +0.043783, -0.027698)
||Delta||_2 = 0.054248
M_X^min = 3.509119e15 GeV
rho_X^min = 1
```

The same illustrative scalar-fragment fit needs only

```text
max |log(M_G/M_r)| = 1.5159,
```

corresponding to masses within a factor `exp(1.5159) ~= 4.55` of `M_G`.
This is a natural threshold scale.

### Novel elements and constraints

1. **Complete-multiplet null direction.**
   Heavy complete GUT multiplets move only the common `alpha_G^-1`; they cannot
   repair relative unification. Only split multiplets or gauge-kinetic moduli
   affect the `Delta` plane.
2. **Proton-safe cone.**
   The scan rejects fake solutions that improve unification only by implicitly
   lowering `M_X` into a proton-decay-dangerous region.
3. **`E_d` gauge-kinetic threshold vector.**
   A controlled modulus contribution can be written
   `f_i(E_d)=f_G+eta_i E_d+...`, giving
   `delta_i^Ed ~= -4pi eta_i E_d`, with `sum_i delta_i^Ed=0`.
   This is allowed only after the heavy spectrum is fixed; it must not be used
   as an unconstrained arbitrary fit vector.

Conclusion for item 5:

The MSSM-like/supersymmetric-threshold branch is viable at one loop and proton
safe. The purely SM one-loop branch requires unnaturally large threshold
differences after imposing the proton bound. The paper-level framework should
therefore proceed with a supersymmetric or SUSY-threshold-inspired unification
branch, while treating the non-SUSY branch as disfavored unless a concrete
intermediate spectrum is supplied.

## Paper Skeleton and Two-Loop Heavy-Spectrum Fit - 2026-05-05

### Paper skeleton status

The current TeX draft has been promoted into a formal paper skeleton:

```text
paper/gut_framework.tex
paper/gut_framework.pdf
paper/refs.bib
```

Structural additions:

1. theorem/proposition/lemma/corollary/definition/assumption/remark environments;
2. a logical-flow TikZ figure for the Spin(10)-to-threshold chain;
3. booktabs tables for Pati-Salam fields, Yukawa singular-value ratios, and the
   two-loop heavy spectrum;
4. BibTeX references for Pati-Salam, early Spin(10), seesaw, and two-loop RGE
   foundations;
5. hidden PDF links and cleaned cross-references.

Verification:

```text
bibtex gut_framework
pdflatex -interaction=nonstopmode -halt-on-error gut_framework.tex
pdflatex -interaction=nonstopmode -halt-on-error gut_framework.tex
```

The latest compile gives a 10-page A4 PDF with no unresolved citations,
undefined references, overfull boxes, or underfull boxes.  Visual pages were
rendered with `pdftoppm` under `tmp/pdfs/gut_framework_visual/`.

### Two-loop RGE fit

The new script

```text
code/two_loop_spectrum_fit.py
```

implements the current gauge-only two-loop refinement:

```text
MZ -> MSUSY: SM two-loop gauge beta functions
MSUSY -> MG: MSSM two-loop gauge beta functions
MSUSY = 1.000e3 GeV
```

The heavy threshold basis used in the complete non-singlet spectrum fit is

```text
H_C + Hbar_C : (2/5, 0, 1)
Sigma_3     : (0, 2, 0)
Sigma_8     : (0, 0, 3)
Delta_i = sum_r b_i^(r) log(MG/Mr)/(2 pi).
```

Best point:

```text
MG = 1.386756e16 GeV
alpha_G^-1 = 25.091863
alpha_i^-1(MG) = (25.096168, 25.036442, 25.180897)
log(MG/M_HC) = +0.067609
log(MG/M_Sigma3) = -0.174110
log(MG/M_Sigma8) = +0.163935
residual_l2 = 6.404746e-15
max_abs_log = 0.174110
```

Fitted spectrum:

```text
M_HC      = 1.296098e16 GeV
M_Sigma3  = 1.650498e16 GeV
M_Sigma8  = 1.177074e16 GeV
M_X       = 1.386756e16 GeV
M_nuR1    = 2.599776e10 GeV
M_nuR2    = 9.144789e13 GeV
M_nuR3    = 3.105386e15 GeV
```

Proton constraints:

```text
M_X^min = 3.451519e15 GeV
rho_X = 1.000000
triplet_filter_required = 1.188726e-04
S_T = 1e-5 is safe
```

Interpretation:

The two-loop gauge-only fit preserves the item-5 MSSM-like branch.  All
non-singlet heavy masses lie within a factor `exp(0.174110) ~= 1.19` of `MG`,
so the threshold spectrum is natural.  The fitted gauge boson mass is well
above the dimension-6 proton lower bound, and the colored-Higgs dimension-5
constraint only requires a geometric triplet filter weaker than the item-4
target `S_T = 1e-5`.

### Next verification targets

1. Add top/bottom/tau and neutrino Yukawa terms to the two-loop gauge running,
   with a scan over `tan beta` and `MSUSY`.
2. Replace the linear heavy-threshold basis by a breaking-chain spectrum model:
   solve for masses from a minimal superpotential rather than treating
   `M_HC`, `M_Sigma3`, and `M_Sigma8` as free logarithms.
3. Propagate the fitted heavy spectrum into the dimension-5 proton operator
   calculation and rescan `S_T`, sfermion mass, and flavor leakage together.
4. Add uncertainty bands for `alpha_s(MZ)`, `sin^2 theta_W`, and the proton
   lifetime bound to turn the present benchmark into a paper-level allowed
   region.

## Yukawa Two-Loop and Breaking-Chain Spectrum Fit - 2026-05-05

### Implemented files

The next refinement is now implemented in

```text
code/scan_yukawa_superpotential_rge.py
output/yukawa_superpotential_rge/yukawa_superpotential_rge_report.md
output/yukawa_superpotential_rge/yukawa_superpotential_rge_summary.json
output/yukawa_superpotential_rge/yukawa_superpotential_rge_scan.csv
```

The paper draft was updated:

```text
paper/gut_framework.tex
paper/gut_framework.pdf
```

### Mathematical refinement

The MSSM two-loop gauge equation now includes the top, bottom, tau, and Dirac
neutrino Yukawa terms:

```text
d alpha_i/dt = b_i alpha_i^2/(2 pi)
             + alpha_i^2/(8 pi^2) [sum_j B_ij alpha_j
               - c_i^t Tr Yu^dag Yu
               - c_i^b Tr Yd^dag Yd
               - c_i^tau Tr Ye^dag Ye
               - c_i^nu Tr Ynu^dag Ynu].
```

The MSSM coefficient vectors are

```text
c^t   = (26/5, 6, 4)
c^b   = (14/5, 6, 4)
c^tau = (18/5, 2, 0)
c^nu  = (6/5, 2, 0)
```

The right-handed-neutrino Yukawa contribution is thresholded at the seesaw
masses already reconstructed in item 3.

### Breaking-chain superpotential spectrum

The heavy threshold logs are no longer treated only as independent free
parameters.  They are mapped into a constrained effective breaking-chain
superpotential:

```text
W_mass = MG lambda_T T Tbar
       + MG lambda_S [(1+2 chi) Tr Sigma_3^2
                      +(1-3 chi) Tr Sigma_8^2]/2 + ...
```

so

```text
M_HC     = lambda_T MG
M_Sigma3 = lambda_S (1+2 chi) MG
M_Sigma8 = lambda_S (1-3 chi) MG
```

Given exact threshold logs `ell_r = log(MG/M_r)`, the analytic map is

```text
kappa_3 = exp(-ell_3)
kappa_8 = exp(-ell_8)
chi = (kappa_3/kappa_8 - 1)/(2 + 3 kappa_3/kappa_8)
lambda_S = kappa_3/(1 + 2 chi)
lambda_T = exp(-ell_T)
```

This makes the threshold fit falsifiable as a superpotential spectrum: a point
can be exact but still rejected or marked tuned if the resulting couplings,
splittings, or cancellations are unacceptable.

### Numerical scan

The scan used

```text
tan beta in {2, 3, 5, 10, 20, 35, 50}
MSUSY in {1e3, 3e3, 1e4, 3e4, 1e5} GeV
MG grid from 10^15.85 to 10^16.75
```

Best proton-safe and dimension-5-safe benchmark:

```text
tan beta = 10
MSUSY = 1.000000e5 GeV
MG = 7.079458e15 GeV
alphaG^-1 = 28.967531
alpha_i^-1(MG) = (28.936255, 29.784595, 28.942088)
yt,yb,ytau,ynu at MG = (0.825810, 0.091706, 0.078759, 0.583294)
lambda_T = 1.634409
lambda_S = 0.404229
chi = -0.405036
log(MG/M_HC) = -0.491282
log(MG/M_Sigma3) = +2.566882
log(MG/M_Sigma8) = +0.110474
residual_l2 = 8.140290e-15
```

Fitted spectrum:

```text
M_HC     = 1.157073e16 GeV
M_Sigma3 = 5.435223e14 GeV
M_Sigma8 = 6.339015e15 GeV
M_X      = 7.079458e15 GeV
M_X^min  = 3.212336e15 GeV
```

Proton checks:

```text
rho_X = 1
S_T_required = 2.017847e-5
tau_dim6 = 2.358936e35 yr
tau_dim5(S_T=1e-5) = 4.071708e34 yr
```

### Interpretation

This is not a fully single-scale natural spectrum.  Across `4445` scanned RGE
points:

```text
proton-safe + dimension-5-safe + perturbative points = 89
safe points with all non-singlet masses within factor 2 of MG = 0
```

The viable point requires an intermediate adjoint-triplet threshold

```text
M_Sigma3/MG = exp(-2.566882) ~= 0.0768,
```

generated by

```text
1 + 2 chi ~= 0.1899.
```

Therefore the updated conditional claim is:

The Yukawa-refined two-loop framework remains viable under the existing proton
constraints only if the breaking-chain superpotential produces a controlled
intermediate `Sigma_3` threshold.  If the theory insists on all non-singlet
heavy states being within a factor `2` of `MG`, this scan gives a no-go under
the current assumptions.

### Next verification targets

1. Derive the `Sigma_3` intermediate threshold from an explicit Spin(10) or
   Pati-Salam breaking superpotential, rather than the effective
   `lambda_S, chi` parameterization.
2. Run an uncertainty scan over `alpha_s(MZ)`, `sin^2 theta_W`, and the proton
   lifetime benchmark to test whether the factor-2 no-go is robust.
3. Replace fixed one-loop Yukawa running by matrix-valued Yukawa RGEs using the
   actual CP1 O(2) textures, including flavor rotations in the dimension-5
   proton operator.
4. Compute whether the intermediate `Sigma_3` threshold has a direct effect on
   neutrino-sector matching or leptogenesis-scale consistency.

## Sigma3 Intermediate Threshold from an Explicit PS Filter - 2026-05-05

### Implemented files

The intermediate adjoint-triplet threshold has now been derived from an explicit
Pati-Salam-stage quadratic superpotential and verified numerically:

```text
code/derive_sigma3_threshold_superpotential.py
output/sigma3_superpotential/sigma3_superpotential_report.md
output/sigma3_superpotential/sigma3_superpotential_summary.json
```

The paper draft was updated with a new derivation:

```text
paper/gut_framework.tex
paper/gut_framework.pdf
```

### Minimal no-go

A pure Spin(10)-singlet quadratic mass

```text
W_univ = (m/2) Tr 45_Sigma^2
```

cannot split the Pati-Salam fragments `(1,3,1)` and `(15,1,1)`.  Its Hessian is
proportional to the identity on the full `45`, so by Schur's lemma it acts as
the same scalar on all Pati-Salam fragments.  Therefore the required

```text
M_Sigma3/M_Sigma8 ~= 0.0857424
```

cannot come from a universal GUT mass.  A non-singlet breaking spurion is
required.

### Explicit Pati-Salam filter superpotential

Use the Spin(10) decompositions

```text
45_Sigma -> (15,1,1) + (1,3,1) + (1,1,3) + (6,2,2)
10_H     -> (6,1,1) + (1,2,2)
```

with

```text
Sigma_L in (1,3,1)  -> Sigma3 = (1,3,0)
Sigma_C in (15,1,1) -> Sigma8 = (8,1,0) + ...
T + Tbar in 10_H    -> H_C + Hbar_C
```

The explicit quadratic superpotential is

```text
W_PS = MG lambda_T T Tbar
     + MG lambda_S/2 [
         Tr Sigma_L^2 + Tr Sigma_C^2
       + chi (2 Tr Sigma_L^2 - 3 Tr Sigma_C^2)
       ].
```

Equivalently, in a Spin(10) completion, this can descend from a schematic
operator

```text
45_Sigma 45_Sigma Theta / M_*
```

where the SM-singlet vev of `Theta` lies in a 54/210-like filter direction and
has Clebsch eigenvalues `+2` on `(1,3,1)` and `-3` on `(15,1,1)`.

The holomorphic mass Hessian in the basis `(H_C, Sigma3, Sigma8)` is

```text
M_Hessian/MG = diag(
  lambda_T,
  lambda_S (1 + 2 chi),
  lambda_S (1 - 3 chi)
).
```

### Solved parameters and numerical verification

Solving against the Yukawa-refined two-loop spectrum gives

```text
lambda_T = 1.634409492
lambda_S = 0.404228616
chi = -0.405035703
1 + 2 chi = 0.189928593
1 - 3 chi = 2.215107110
```

Then

```text
MG = 7.079458e15 GeV
M_HC     = 1.157073e16 GeV
M_Sigma3 = 5.435223e14 GeV
M_Sigma8 = 6.339015e15 GeV
max relative mass error = 4.733e-16
```

The same logs reproduce the two-loop threshold matching:

```text
alphaG^-1 = 28.967530959
threshold residual_l2 = 8.140e-15
tau_dim6 = 2.358936e35 yr
tau_dim5(S_T=1e-5) = 4.071708e34 yr
```

### Nontrivial elements

1. **Clebsch filter.**  The split is controlled by fixed group-theoretic
   eigenvalues `(+2,-3)`, not by an arbitrary independent threshold mass.
2. **Critical adjoint.**  The intermediate `Sigma3` is the near-critical mode
   of the factor `1 + 2 chi`; the exact critical limit gives a massless chiral
   adjoint and an enhanced massless-sector symmetry.
3. **Proton-friendly split.**  The filter lowers only `Sigma3`; the colored
   Higgs remains above `MG`, so the dimension-5 proton constraint remains safe
   for `S_T = 1e-5`.
4. **Falsifiable Spin(10) requirement.**  A future full 54/210 Clebsch
   calculation must reproduce the `(+2,-3)` effective filter or an equivalent
   projection.  Otherwise the current single-spurion explanation fails, and the
   factor-2 single-scale no-go from the previous scan becomes the operative
   conclusion.

### Next verification targets

1. Write the full Spin(10) invariant tensor contraction for
   `45_Sigma 45_Sigma Theta`, choosing explicitly between `54_H` and `210_H`,
   and compute the Clebsch eigenvalues rather than parameterizing them.
2. Check whether the same filter spurion introduces unwanted masses or light
   states in `(1,1,3)` or `(6,2,2)`.
3. Include the intermediate `Sigma3` explicitly in piecewise RGE running rather
   than treating it only as a threshold log at `MG`.
4. Test whether the near-critical `1 + 2 chi` condition is stable under
   plausible Planck-suppressed operators and `E_d` modulus corrections.

## Direct 45-45-54 Clebsch Tensor Test - 2026-05-05

Status: completed as a direct Spin(10)/SO(10)-index calculation and promoted
to a no-go constraint for a lone `54_H`.

### Invariant and first-principles Clebsch calculation

Take the adjoint `45` as an antisymmetric tensor `A_ij=-A_ji` and the `54` as
a symmetric traceless tensor `S_ij=S_ji`, `Tr S=0`.  The unique quadratic
invariant used for the direct test is

```text
W_54 = (eta/2) S_ij A_ik A_jk.
```

For the Pati-Salam preserving direction

```text
S = diag(-2,-2,-2,-2,-2,-2, 3,3,3,3),
```

the Hessian eigenvalue on a plane generator `A_ab` is exactly `s_a+s_b`.
Therefore the decomposition under `SO(6)xSO(4) ~= SU(4)_C x SU(2)_L x SU(2)_R`
gives

```text
(15,1,1):             multiplicity 15, raw Clebsch -4, normalized -4/3
(6,2,2):              multiplicity 24, raw Clebsch +1, normalized +1/3
(1,3,1)+(1,1,3):      multiplicity  6, raw Clebsch +6, normalized +2
```

This proves that a single `54_H` cannot reproduce the previous effective
`(+2,-3)` filter and cannot split `Sigma_L=(1,3,1)` from
`Sigma_R=(1,1,3)`.

### Numerical threshold check

Fitting only the target `Sigma3/Sigma8` masses with the direct `54_H`
Clebsches gives

```text
lambda_54 = 0.567955638
chi_54    = -0.432411471
M_SigmaR  = 5.435223e14 GeV
M_(6,2,2) = 3.441269e15 GeV
```

The target `Sigma3` and `Sigma8` masses are matched, but the literal single-54
Hessian also leaves `Sigma_R` degenerate with `Sigma_L` and leaves a moderately
light `(6,2,2)` fragment.  Keeping these fragments in the one-loop threshold
matching changes the residual from

```text
base residual_l2             = 8.140e-15
literal 54 extra residual_l2 = 4.493821e-01
```

so the lone `54_H` completion is rejected.

### Artifacts

```text
code/compute_54_clebsch.py
output/clebsch_54/clebsch_54_report.md
output/clebsch_54/clebsch_54_summary.json
paper/gut_framework.tex
paper/gut_framework.pdf
```

### Next verification target

Move to a direct `210_H` four-form contraction, or to a controlled `54_H+210_H`
linear combination/lifting sector.  The next test must compute the exact
Clebsch eigenvalues on `(1,3,1)`, `(1,1,3)`, `(15,1,1)`, and `(6,2,2)` and then
rerun the threshold/proton-decay scan with all non-decoupled fragments included.

## Direct 45-45-210 Clebsch and Projector-Lift Test - 2026-05-06

Status: completed as a direct SO(10) four-form contraction.  The result fixes
the `Sigma_L/Sigma_R` chirality problem, but exact matching still requires a
controlled lifting/projector sector for `(6,2,2)`.

### Direct 210_H invariant

Represent `210_H` by an antisymmetric four-form `Phi_ijkl`.  The invariant is

```text
W_210 = (zeta/8) Phi_ijkl A_ij A_kl.
```

For the Pati-Salam preserving D-parity direction

```text
Phi_{7,8,9,10} = +1
```

on the weak SO(4) block, using one-based SO(10) labels, the induced operator
on adjoint two-forms is the SO(4) Hodge star.  The full 45-by-45 matrix has
eigenvalue multiplicities

```text
-1: 3
 0: 39
+1: 3
```

Thus, after orientation choice,

```text
D_210(Sigma_L=(1,3,1)) = +1
D_210(Sigma_R=(1,1,3)) = -1
D_210(Sigma8)          =  0
D_210(X=(6,2,2))       =  0
```

This is the missing nontrivial element: unlike a lone `54_H`, the `210_H`
four-form distinguishes the two weak adjoints while preserving
Pati-Salam gauge symmetry.

### 54_H + 210_H spectrum

Use the normalized linear mass operator

```text
kappa_r = mu + a54 F_54(r) + b210 D_210(r)
F_54 = (2, 2, -4/3, 1/3)
D_210 = (+1, -1, 0, 0)
```

on `(Sigma_L, Sigma_R, Sigma8, X_622)`.  Fixing the verified
`Sigma_L/Sigma8` masses and placing `Sigma_R` at `MG` gives

```text
mu   = 0.752600723
a54  = -0.107106719
b210 = -0.461612714

kappa_Sigma3 = 0.076774572
kappa_Sigma8 = 0.895409681
kappa_SigmaR = 1.000000000
kappa_X_622  = 0.716898484
```

The near-critical triplet is preserved by cancellation
`mu + 2 a54 + b210 = 0.076774572`.

### Spectral-projector lift

The remaining light `(6,2,2)` fragment is selected by the polynomial projector

```text
P_X = -(9/25)(F_54 - 2)(F_54 + 4/3)
```

which equals one on `F_54=1/3` and zero on `F_54=2,-4/3`.  Adding

```text
Delta W_X = MG delta_X A P_X A
delta_X = 0.283101516
```

raises `(6,2,2)` to `MG` without changing `Sigma_L`, `Sigma_R`, or `Sigma8`.
The companion projector

```text
P_R = (D_210^2 - D_210)/2
```

is available if a nearby branch needs a direct `Sigma_R` lift.

### Numerical verification

After the projector lift:

```text
MG        = 7.079458e15 GeV
M_Sigma3  = 5.435223e14 GeV
M_Sigma8  = 6.339015e15 GeV
M_SigmaR  = 7.079458e15 GeV
M_X_622   = 7.079458e15 GeV

extra residual after projector lift = 0.000e+00
base threshold residual             = 8.140e-15
tau_dim6                            = 2.358936e35 yr
tau_dim5(S_T=1e-5)                  = 4.071708e34 yr
```

For comparison, the best linear `54_H+210_H` point without the projector lift
still leaves extra-threshold residual `6.281028e-02`.  Therefore the
projector/lifting sector is required for exact matching.

### Artifacts

```text
code/compute_210_lifting.py
output/clebsch_210/clebsch_210_lifting_report.md
output/clebsch_210/clebsch_210_lifting_summary.json
paper/gut_framework.tex
paper/gut_framework.pdf
```

### Next verification target

Realize the spectral projectors through a renormalizable Spin(10) Higgs
superpotential, for example by introducing explicit heavy mediator multiplets
whose tree-level elimination generates the `P_X` and optional `P_R` operators.
Then rerun the two-loop RGE/proton scan with the resulting physical spectrum
instead of the effective projector masses.

## Renormalizable 45-Mediator Projector Realization - 2026-05-06

Status: completed as an explicit tree-level UV realization of the `P_X`
projector using only mass terms and cubic Spin(10)-invariant couplings.  The
remaining issue is not the existence of the projector, but the small heavy
mediator threshold it induces when `M_med` is finite.

### Superpotential

Introduce two heavy adjoint chiral multiplets `B_45` and `C_45` and write

```text
W/MG = 1/2 A K0 A
     + R B C
     + eta A (F54 - 2) B
     + eta A (F54 + 4/3) C

K0 = mu + a54 F54 + b210 D210
```

where `A=45_Sigma` and `R=M_med/MG`.  The shifted factors are ordinary
bilinear masses plus renormalizable `45-45-54` cubic terms after the `54_H`
vev is inserted.

Solving the mediator F-terms gives the Schur complement

```text
K_eff = K0 - (2 eta^2/R)(F54 - 2)(F54 + 4/3)
      = K0 + [50 eta^2/(9R)] P_X.
```

Therefore the effective `(6,2,2)` projector is generated by a renormalizable
Spin(10) mediator sector rather than assumed by hand.

### Finite mediator benchmark

The full physical mass matrix in a Clebsch sector `(f,d)` is

```text
M/MG = [[mu + a54 f + b210 d, eta(f-2), eta(f+4/3)],
        [eta(f-2),             0,        R],
        [eta(f+4/3),           R,        0]]
```

Solving the exact 3-by-3 eigenvalues at `R=50` gives

```text
mu   = 0.761049382
a54  = -0.108308983
b210 = -0.466795262
eta  = 1.589022391
Schur delta_X = 0.280554684
light fit residual_l2 = 3.980e-14
```

Light eigenvalues:

```text
Sigma_L      kappa = 0.076774572  M = 5.435223e14 GeV
Sigma_R      kappa = 1.000000000  M = 7.079458e15 GeV
Sigma8_block kappa = 0.895409681  M = 6.339015e15 GeV
X_622        kappa = 1.000000000  M = 7.079458e15 GeV
```

The heavy mediator eigenvalues are clustered near `50 MG ~= 3.54e17 GeV`.
Using complete `45` block beta vectors for these two heavy eigenstates gives

```text
projected heavy-threshold residual_l2 = 1.287043e-03
projected vector = [1.1998788e-04, 8.4413082e-04, -9.6411870e-04]
```

Decoupling scan:

```text
R=5    eta=0.474805  residual=1.367372e-02
R=10   eta=0.695423  residual=6.568905e-03
R=20   eta=0.997528  residual=3.238114e-03
R=50   eta=1.589022  residual=1.287043e-03
R=100  eta=2.252401  residual=6.424120e-04
```

This is a controlled UV completion: increasing `R` suppresses the
non-universal heavy threshold, while the couplings remain perturbative through
the displayed range.

### Artifacts

```text
code/construct_uv_projector_mediator.py
output/uv_projector_mediator/uv_projector_mediator_report.md
output/uv_projector_mediator/uv_projector_mediator_summary.json
paper/gut_framework.tex
paper/gut_framework.pdf
```

### Next verification target

Rerun the two-loop gauge/Yukawa/proton fit with the finite mediator threshold
vector included.  In parallel, replace the current simplified `Sigma8_block`
treatment by a complete Goldstone/lifting assignment for the full `(15,1,1)`
block so that the light octet threshold is no longer standing in for a whole
Pati-Salam multiplet.

## Finite Mediator Threshold RGE and `(15,1,1)` Assignment - 2026-05-06

Status: completed for the fixed `R=50` mediator benchmark.  The finite
45-mediator threshold vector has been inserted into the cached two-loop
gauge/Yukawa/proton scan, and the simplified `Sigma8_block` has been replaced
by an explicit `(15,1,1)` Goldstone/lifting assignment.

### Matching theorem now used

The heavy matching equation is

```text
alpha_i^-1(MG) = alpha_G^-1 + Delta_med_i
               + sum_r b_i^r log(MG/M_r)/(2 pi),
r = H_C, Sigma3, Sigma8_octet.
```

Only the traceless projection changes the fitted mass logs:

```text
P (alpha^-1(MG) - Delta_med) = P B ell / (2 pi).
```

For the `R=50` UV mediator completion,

```text
Delta_med = (-9.971399911, -9.970675768, -9.972484018)
P Delta_med = (+1.199878759e-04, +8.441308213e-04, -9.641186972e-04)
||P Delta_med||_2 = 1.287043e-03
```

The large finite threshold is mostly a universal shift of `alphaG_inv`; the
non-universal part is small but has now been included explicitly.

### Historical best post-mediator scan point before the corrected component card

The numbers in this subsection are retained only as provenance for the older
finite-mediator scan.  The current corrected component-Hessian card is the
R=200 result recorded below under
`verify_spin10_component_hessian.py`.

```text
tan beta = 10.00
MSUSY = 1.000000e+05 GeV
MG = 7.079458e+15 GeV
alphaG^-1 = 38.938877
lambda_T = 1.633034289
lambda_S = 0.403573619
chi = -0.404680988
logs = (-0.490439811, +2.564775172, +0.112576116)
residual_l2_after_mediator = 6.153e-15
```

Physical comparison:

```text
M_HC:       1.157073e16 -> 1.156100e16 GeV
M_Sigma3:   5.435223e14 -> 5.446685e14 GeV
M_Sigma8:   6.339015e15 -> 6.325703e15 GeV
tau_dim6:   2.358936e35 -> 4.262457e35 yr
tau_dim5:   4.071708e34 -> 4.064859e34 yr  (S_T=1e-5)
```

Scan count:

```text
Total cached two-loop/Yukawa points = 4445
Proton-safe + dimension-5-safe + perturbative points = 89
Safe single-scale factor-2 points = 0
```

Conclusion: the viable branch survives the finite mediator threshold, but it is
still an intermediate-`Sigma3` branch rather than a single-scale branch.

### Complete `(15,1,1)` assignment

```text
(15,1,1) -> (8,1,0) + (3,1,2/3) + (bar3,1,-2/3) + (1,1,0)
```

Assignments:

```text
(8,1,0):
  physical light octet Sigma8
  b = (0,0,3)
  log(MG/M) = 0.599340563 in the current R=200 component card

(3,1,2/3) + (bar3,1,-2/3):
  SU(4)_C-breaking Goldstone pair, eaten or lifted at MG
  b = (8/5,0,1)
  log(MG/M) = 0

(1,1,0):
  radial singlet at M_Sigma8/2 in the current component card
  b = (0,0,0)
  log(MG/M) = 1.292488, but one-loop SM gauge threshold is zero
```

If the colored triplet pair is incorrectly left at `M_Sigma8`, it contributes
a projected threshold norm `1.090376e-01` in the current R=200 component card, far larger than the finite mediator
non-universal residual.  Therefore the Goldstone/lifting assignment is a
consistency condition, not just notation.

### Artifacts

```text
code/scan_mediator_threshold_rge.py
output/mediator_threshold_rge/mediator_threshold_rge_report.md
output/mediator_threshold_rge/mediator_threshold_rge_summary.json
output/mediator_threshold_rge/mediator_threshold_rge_scan.csv
paper/gut_framework.tex
paper/gut_framework.pdf
```

### Next verification target

Promote the fixed `R=50` benchmark into a full superpotential-spectrum fit:
scan `R`, the Goldstone-sector masses, and the dimension-five dressing/filter
parameters simultaneously.  The next numerical target is to preserve the
post-mediator two-loop fit while minimizing the required `Sigma3` criticality
and keeping the colored Goldstone pair at `MG`.

### Heartbeat continuation - 2026-05-06 01:44 CST

Task status: incomplete.  The fixed `R=50` mediator benchmark is verified, but
the paper-level theory still needs a genuine spectrum fit in which `R`,
Goldstone/lifting masses, and dimension-five dressing parameters are varied
together rather than held fixed.

Current obstacle:

```text
The UV projector sector has been shown to exist at one benchmark point.
It has not yet been converted into a continuous, constrained superpotential
fit with proton-decay and threshold constraints imposed simultaneously.
```

Next nontrivial idea:

```text
Schur-fixed mediator decoupling theorem.

Hold delta_X = 50 eta^2/(9 R) fixed so that integrating out the two 45
mediators preserves the light projector eigenvalue.  Then eta scales as
sqrt(R), the light spectrum is stationary to leading Schur order, and the
non-universal heavy threshold should decay approximately as 1/R because the
leading heavy states form almost-complete 45 blocks whose universal threshold
is absorbed into alpha_G^-1.
```

Numerical check from the existing UV mediator scan:

```text
R       eta          eta/sqrt(R)   ||P Delta_med||   R ||P Delta_med||
5       0.474805474  0.212339463   1.367372e-02     6.836860e-02
10      0.695422832  0.219912009   6.568905e-03     6.568905e-02
20      0.997527783  0.223053993   3.238114e-03     6.476228e-02
50      1.589022391  0.224721702   1.287043e-03     6.435213e-02
100     2.252401149  0.225240115   6.424120e-04     6.424120e-02
```

Interpretation: `eta/sqrt(R)` approaches a constant and
`R ||P Delta_med||` approaches `6.4e-2`, giving a clean numerical signature of
the expected `1/R` decoupling.  This makes the finite mediator threshold a
controlled deformation rather than a new arbitrary threshold vector.

Verification plan for the next active step:

```text
1. Generalize code/scan_mediator_threshold_rge.py so Delta_med is generated
   for each R from the full 3-by-3 mediator mass matrices rather than read
   only from the R=50 JSON benchmark.
2. For each R, rerun the cached two-loop RGE/proton scan and record the best
   viable point, alpha_G shift, M_Sigma3 criticality, tau_dim6, and tau_dim5.
3. Add Goldstone-sector mass variables and require the colored
   (3,1,2/3)+(bar3,1,-2/3) pair to remain at MG within a chosen tolerance.
4. Only after the R-window survives these checks, promote the result into
   paper/gut_framework.tex as a proposition/table.
```

### Heartbeat continuation - 2026-05-06 02:17 CST

Status: partially completed.  The `R`-window part of the previous verification
plan has now been implemented as a dedicated scan, and the result has been
promoted into the TeX draft.  The broader theory is still incomplete because
Goldstone-sector masses and dimension-five dressing/filter parameters are not
yet scanned simultaneously.

New artifact:

```text
code/scan_mediator_r_window.py
output/mediator_r_window/mediator_r_window_report.md
output/mediator_r_window/mediator_r_window_summary.json
output/mediator_r_window/mediator_r_window_scan.csv
paper/gut_framework.tex
```

Mathematical statement tested:

```text
Schur-fixed mediator decoupling:
delta_X = 50 eta^2/(9R) fixed  =>  eta ~ sqrt(R),
and the non-universal heavy threshold obeys ||P Delta_med|| = O(1/R).
```

Numerical result:

```text
fit over R >= 20:
log ||P Delta_med|| = -2.727977 - 1.003607 log R
expected slope = -1

R      eta/sqrt(R)   ||P Delta_med||   R ||P Delta_med||   safe points
20     0.223053993   3.238114e-03     6.476228e-02        90
50     0.224721702   1.287043e-03     6.435213e-02        89
100    0.225240115   6.424120e-04     6.424120e-02        89
200    0.225492161   3.209536e-04     6.419072e-02        89
```

Best decoupled benchmark in the displayed window:

```text
R = 200
eta = 3.188940720
MG = 7.079458e15 GeV
alphaG^-1 = 42.461991
M_Sigma3 = 5.438076e14 GeV
M_Sigma8 = 6.335690e15 GeV
tau_dim6 = 5.068669e35 yr
tau_dim5(S_T=1e-5) = 4.069999e34 yr
safe cached two-loop points = 89
safe single-scale factor-2 points = 0
```

Interpretation: the finite mediator threshold is now shown to be a controlled
decoupling deformation, not an arbitrary threshold dial.  Increasing `R`
suppresses the non-universal threshold while leaving the intermediate-`Sigma3`
branch intact.

Remaining obstacle:

```text
The colored Goldstone/lifting sector is assigned consistently at MG, but its
renormalizable superpotential masses are not yet included as variables in the
same scan as R and the dimension-five proton-decay dressing/filter parameters.
```

Next attempted nontrivial idea:

```text
Goldstone-locking constrained scan.

Introduce a tolerance epsilon_G for the colored
(3,1,2/3)+(bar3,1,-2/3) pair:
|log(MG/M_Goldstone)| <= epsilon_G.
Then scan (R, epsilon_G, S_T) and require simultaneously:
  1. exact projected gauge matching,
  2. tau_dim6 and tau_dim5 above targets,
  3. ||P Delta_med|| below a chosen tolerance,
  4. no accidental colored-pair threshold larger than the mediator residual.
```

Next verification plan:

```text
1. Extend the R-window script with epsilon_G and S_T grids.
2. Add the colored-pair threshold vector (8/5,0,1) log(MG/M_Goldstone)/(2pi)
   before solving the H_C/Sigma3/Sigma8 logs.
3. Determine the largest allowed epsilon_G as a function of R and S_T.
4. Promote the result to TeX only if a nonempty allowed region remains.
```

### Heartbeat continuation - 2026-05-06 04:50 CST

Status: partially completed.  The Goldstone-locking and dimension-five filter
scan has now been implemented and promoted into the TeX draft.  The theory is
still incomplete because the locking relation is currently imposed as a
constraint; it has not yet been derived from an explicit renormalizable
Pati-Salam breaking superpotential.

New artifact:

```text
code/scan_goldstone_locking.py
output/goldstone_locking/goldstone_locking_report.md
output/goldstone_locking/goldstone_locking_summary.json
output/goldstone_locking/goldstone_locking_scan.csv
paper/gut_framework.tex
paper/gut_framework.pdf
```

Mathematical statement tested:

```text
Delta_G = (8/5,0,1) epsilon_G/(2 pi)
epsilon_G = log(MG/M_Goldstone)
locking condition: ||P Delta_G||_2 <= ||P Delta_med(R)||_2
```

The unit conversion is

```text
||P Delta_G||_2 = 1.819292536e-01 |epsilon_G|.
```

Numerical result at `S_T=1e-5`:

```text
R      max |epsilon_G|   safe points   tau_dim5 [yr]
20     1.779875e-02     91            4.038471e34
50     7.074413e-03     90            4.058474e34
100    3.531109e-03     90            4.065098e34
200    1.764167e-03     89            4.068404e34
```

Additional scan result:

```text
S_T = 1e-4 has no safe locked point on the tested grid.
S_T = 3e-5, 1e-5, 3e-6, and 1e-6 all admit locked points.
No safe single-scale factor-2 point appears.
```

Interpretation: stronger mediator decoupling forces stronger Goldstone
locking.  Since `||P Delta_med|| ~ 1/R`, the allowed Goldstone displacement
also scales approximately as `|epsilon_G|max ~ 1/R`.  At `R=200`, the colored
pair must remain within about `0.18%` of `MG`.

Remaining obstacle:

```text
The scan proves a quantitative Goldstone-locking requirement, but the
renormalizable superpotential that enforces
|log(MG/M_Goldstone)| <= epsilon_G,max has not been constructed.
```

Next attempted nontrivial idea:

```text
Goldstone-locking mass theorem.

Construct a Pati-Salam breaking sector with adjoint/singlet mediator fields
whose F-flatness equations make the colored (3,1,2/3)+(bar3,1,-2/3)
directions exact Goldstone modes at the gauge-breaking scale, while assigning
the octet a controlled non-Goldstone mass.  The target is to derive
M_Goldstone/MG = 1 + O(Delta_med) rather than impose it by hand.
```

Next verification plan:

```text
1. Write the minimal PS-breaking superpotential for SU(4)_C -> SU(3)_C x U(1).
2. Compute the Hessian in the (15,1,1) sector and identify exact eaten modes.
3. Add lifting fields for the radial singlet without shifting the colored
   Goldstone eigenvalue beyond epsilon_G,max.
4. Feed the derived masses back into code/scan_goldstone_locking.py and verify
   the locked region survives without imposing epsilon_G by hand.
```

## External Review Triage - 2026-05-06

Status: accepted as a major-revision roadmap update.  The external review is
mostly correct.  The manuscript should be repositioned as a conditional
Spin(10) EFT model-building note, not as an unconditional derivation of a
unique GUT from PSLT/Another Physics.

### P0: two-loop Yukawa normalization audit

The review's most serious technical objection is valid.  The paper currently
prints

```text
d alpha_i/dt = b_i alpha_i^2/(2 pi)
             + alpha_i^2/(8 pi^2)
               [sum_j B_ij alpha_j - c_i^a Tr(Y_a^dag Y_a)].
```

The code currently does the same in
`code/scan_yukawa_superpotential_rge.py`: `yukawa_drag` is built from ordinary
Yukawa squares such as `yt*yt`, not from `alpha_y = y^2/(4 pi)`.

If `Y` denotes ordinary Yukawa matrices, the alpha-form equation should be

```text
d alpha_i/dt = b_i alpha_i^2/(2 pi)
             + alpha_i^2/(8 pi^2)
               [sum_j B_ij alpha_j - c_i^a Tr(Y_a^dag Y_a)/(4 pi)].
```

Action items:

```text
1. Decide and document the convention: ordinary Yukawa y or alpha_y=y^2/(4pi).
2. Patch code/scan_yukawa_superpotential_rge.py and dependent mediator scans.
3. Rerun all two-loop/Yukawa/proton/mediator/Goldstone scans.
4. Update paper/gut_framework.tex tables and conclusions.
5. Treat current post-two-loop benchmark numbers as provisional until this is done.
```

Heartbeat audit result:

```text
artifact:
  code/audit_yukawa_4pi_normalization.py
  output/yukawa_4pi_audit/yukawa_4pi_audit_report.md
  output/yukawa_4pi_audit/yukawa_4pi_audit_summary.json
  output/yukawa_4pi_audit/yukawa_4pi_corrected_scan.csv

corrected convention:
  two-loop alpha_i equation uses Tr(Y^dag Y)/(4 pi)
```

The corrected scan does not kill the branch, but it moves the best spectrum:

```text
quantity                      original              4pi-corrected
MG                            7.079458e15 GeV       7.079458e15 GeV
alphaG^-1                     28.967531             27.867260
M_HC                          1.157073e16 GeV       9.445264e15 GeV
M_Sigma3                      5.435223e14 GeV       9.204404e14 GeV
M_Sigma8                      6.339015e15 GeV       3.891045e15 GeV
tau_dim6                      2.358936e35 yr        2.183141e35 yr
tau_dim5(S_T=1e-5)            4.071708e34 yr        2.913416e34 yr
safe points                   89                    146
safe factor-2 single-scale    0                     0
```

Conclusion: the intermediate-`Sigma3` qualitative branch survives this audit,
but every downstream table based on the old cached Yukawa scan is provisional.
The next concrete step is to regenerate finite mediator, R-window, and
Goldstone-locking scans from the corrected cache rather than from
`output/yukawa_superpotential_rge/yukawa_superpotential_rge_scan.csv`.

Heartbeat downstream replay result:

```text
artifact:
  code/scan_corrected_downstream.py
  output/corrected_downstream/corrected_downstream_report.md
  output/corrected_downstream/corrected_downstream_summary.json
  output/corrected_downstream/corrected_R_window_scan.csv
  output/corrected_downstream/corrected_goldstone_locking_scan.csv
```

Corrected fixed `R=50` mediator result:

```text
safe points = 146
safe factor-2 single-scale points = 0
alphaG^-1 = 37.838606
M_Sigma3 = 9.223815e14 GeV
M_Sigma8 = 3.882874e15 GeV
tau_dim6 = 4.024977e35 yr
tau_dim5(S_T=1e-5) = 2.908515e34 yr
```

Corrected `R`-window result:

```text
fit over R >= 20:
log ||P Delta_med|| = -2.727977 - 1.003607 log R

R      ||P Delta_med||   safe points   alphaG^-1   M_Sigma3 [GeV]   tau_dim6 [yr]
20     3.238114e-03     148           35.518694   9.253640e14      3.546559e35
50     1.287043e-03     146           37.838606   9.223815e14      4.024977e35
100    6.424120e-04     146           39.599014   9.214078e14      4.408206e35
200    3.209536e-04     146           41.361720   9.209235e14      4.809394e35
```

Corrected Goldstone-locking result at `S_T=1e-5`:

```text
R      max |epsilon_G|   safe points   tau_dim5 [yr]
20     1.779875e-02     151           2.889634e34
50     7.074413e-03     148           2.903946e34
100    3.531109e-03     148           2.908686e34
200    1.764167e-03     146           2.911052e34
```

Conclusion: the corrected-cache downstream replay keeps the intermediate
`Sigma3` branch alive and preserves the `1/R` mediator decoupling and
Goldstone-locking structure.  The old provisional downstream numbers should be
replaced by the corrected table above in the paper.  The next P0 item is now
the stronger proton-decay target/channel update and the explicit
PS-breaking superpotential derivation of Goldstone locking.

### P0: update proton-decay targets and channels

The local code and paper still use the benchmark `tau > 1e34 yr`.  The review
claims the modern `p -> e+ pi0` limit is `2.4e34 yr`, implying

```text
M_X_min(new) = (2.4)^(1/4) M_X_min(old) = 1.2457 M_X_min(old).
```

This external number must be verified with current experimental references
before finalizing, but the roadmap should adopt the stronger-bound workflow.

Heartbeat progress, 2026-05-06:

Added a no-web, channel-bound replay script:

```text
code/scan_proton_channel_bounds.py
output/proton_channel_bounds/proton_channel_bounds_report.md
output/proton_channel_bounds/proton_channel_bounds_summary.json
output/proton_channel_bounds/proton_channel_bound_replay.csv
output/proton_channel_bounds/proton_mx_bound_scaling.csv
```

Mathematical scaling used:

```text
Gamma_6 = K r_ch (g_G^2/M_X^2)^2,
M_X^min = [K r_ch g_G^4 tau_bound sec/yr / hbar]^(1/4).

tau_5(S_T) = tau_5(S0) (S0/S_T)^2 / r_ch,
S_T^max = S0 sqrt[tau_5(S0)/(r_ch tau_bound)].
```

For the review-supplied local input `tau(p -> e+ pi0) > 2.4e34 yr`,
the dimension-6 mass bound scales by

```text
(2.4)^(1/4) = 1.2446659546.
```

Numerically:

```text
alphaG^-1=24: M_X^min 3.529158e15 -> 4.392623e15 GeV
alphaG^-1=30: M_X^min 3.156575e15 -> 3.928881e15 GeV
alphaG^-1=40: M_X^min 2.733674e15 -> 3.402511e15 GeV
```

Corrected downstream benchmark checks:

```text
fixed R=50:
  tau_d6 = 4.024977e35 yr, margin over 2.4e34 = 16.771
  tau_d5(S_T=1e-5) = 2.908515e34 yr
  S_T^max for tau_5 > 1e34      = 1.705437e-5
  S_T^max for tau_5 > 2.4e34    = 1.100855e-5
  S_T^max for future tau_5>1e35 = 5.393065e-6

R=200:
  tau_d6 = 4.809394e35 yr, margin over 2.4e34 = 20.039
  tau_d5(S_T=1e-5) = 2.912193e34 yr
  S_T^max for tau_5 > 2.4e34 = 1.101551e-5
```

Interpretation: the corrected branch is comfortably dimension-6 safe under
the stronger `e+ pi0` local input.  The fragile constraint is dimension five:
`S_T=1e-5` barely passes a `2.4e34 yr` stress target, while a future
`1e35 yr` target would force a stronger triplet filter near `5.4e-6`.
The numerical workflow is now in place, but the actual channel limits
(`p -> mu+ pi0`, `p -> K+ nu_bar`, `p -> K0 mu+`) still need verified
experimental references before final-paper use.

Action items:

```text
1. Replace hard-coded TAU_TARGET_YEARS=1e34 in the older scans by the new
   channel-bound input table.
2. Verify and cite p -> e+ pi0, p -> mu+ pi0, p -> K+ nu_bar, p -> K0 mu+
   limits, then update the local channel-bound table.
3. Redo proton-safe cone, M_X lower bounds, and all threshold scans using the
   verified table rather than placeholder/stress targets.
4. Add hadronic matrix element and dressing uncertainties.
```

### P0: derive Goldstone locking from a superpotential

The review is correct that the current Goldstone-locking condition is imposed
as a scan constraint:

```text
||P Delta_G||_2 <= ||P Delta_med(R)||_2.
```

The scan is useful, but it is not yet a full Spin(10)/Pati-Salam model.

Heartbeat progress, 2026-05-06:

Constructed the SU(4)_C Pati-Salam-stage adjoint superpotential

```text
W_C = (m_C/2) Tr Sigma_C^2 + (lambda_C/3) Tr Sigma_C^3,
Sigma_0 = v diag(1,1,1,-3),  v = m_C/(2 lambda_C).
```

New artifacts:

```text
code/derive_ps_goldstone_locking.py
output/ps_goldstone_locking/ps_goldstone_locking_report.md
output/ps_goldstone_locking/ps_goldstone_locking_summary.json
```

The projected F-term on the traceless adjoint is

```text
F_C = m_C Sigma_C
    + lambda_C [Sigma_C^2 - (Tr Sigma_C^2)/4 I_4].
```

At the above vacuum, the F-flat equations reduce to
`m_C - 2 lambda_C v = 0`, so

```text
SU(4)_C -> SU(3)_C x U(1)_{B-L}.
```

Linearizing,

```text
delta F_C(X) = m_C X + lambda_C(
    Sigma_0 X + X Sigma_0 - 1/2 Tr(Sigma_0 X) I_4
).
```

Explicit matrix-basis Hessian check:

```text
sector             multiplicity    Hessian eigenvalue/m_C    residual
octet              8               2                         0
goldstone_3        3               0                         0
goldstone_3bar     3               0                         0
singlet            1              -1                         0
```

Thus the colored `(3,1,2/3)+(bar3,1,-2/3)` pair is an exact chiral Goldstone
multiplet for the broken SU(4)_C generators and is absorbed into the massive
vector multiplets.  It has no independent chiral threshold log:

```text
epsilon_G = log(MG/M_Goldstone) = 0,
Delta_G = 0.
```

Corrected-cache numerical lock:

```text
fixed R=50:
  M_octet = 3.882874e15 GeV
  M_singlet = 1.941437e15 GeV
  epsilon_G derived/max = 0 / 7.074413e-3

R=200:
  M_octet = 3.889004e15 GeV
  M_singlet = 1.944502e15 GeV
  epsilon_G derived/max = 0 / 1.764167e-3
```

Interpretation: the PS-stage SU(4)_C Goldstone lock is now derived rather than
imposed.  The remaining caveat is embedding the cubic PS effective sector into
the final full Spin(10) `54_H/210_H/16_H` breaking superpotential, or stating
explicitly that this is the renormalizable Pati-Salam EFT below the first
Spin(10)-breaking step.

Heartbeat progress, 2026-05-06, second pass:

Added a Spin(10) cubic embedding audit:

```text
code/audit_spin10_cubic_embedding.py
output/spin10_cubic_embedding/spin10_cubic_embedding_report.md
output/spin10_cubic_embedding/spin10_cubic_embedding_summary.json
```

Mathematical result:

```text
For a Spin(10) adjoint represented by a real antisymmetric matrix A,
A^T=-A implies

Tr A^3 = Tr[(A^3)^T] = Tr[(A^T)^3] = Tr[-A^3] = -Tr A^3,

so Tr A^3 = 0.
```

Equivalently, the primitive invariant degrees of `D5=so(10)` are

```text
2, 4, 5, 6, 8,
```

so there is no one-field symmetric cubic invariant for `45^3`.  The
antisymmetric `f_abc Phi^a Phi^b Phi^c` also vanishes for one commuting chiral
superfield.

Numerical audit:

```text
SO(10): 200 random antisymmetric 10x10 matrices
  max |Tr A^3| = 1.954e-14
  mean |Tr A^3| = 4.614e-15

SU(4)_C breaking direction H=diag(1,1,1,-3):
  Tr H = 0
  Tr H^2 = 12
  Tr H^3 = -24
  normalized Tr H^3 = -0.577350
```

Interpretation: the PS-stage Goldstone lock is valid, but the PS cubic is not
the restriction of a primitive one-field Spin(10) `45^3` invariant.  The next
full-Spin(10) target is now precise: compute an explicit `210_H^3` or
`210`-mediated contraction that matches onto the `(15,1,1)` SU(4)_C cubic, or
state the term as a renormalizable Pati-Salam EFT interaction after the first
Spin(10)-breaking step.

Follow-up result:

Implemented the explicit `210_H^3` matching calculation:

```text
code/compute_210_cubic_matching.py
output/210_cubic_matching/210_cubic_matching_report.md
output/210_cubic_matching/210_cubic_matching_summary.json
```

Represent the `210_H` as a four-form `Phi in Lambda^4 R^10` and let it act on
adjoint two-forms by

```text
(D_Phi)_[ij],[kl] = Phi_ijkl.
```

Then

```text
I3(Phi) = Tr_{Lambda^2 R^10}(D_Phi^3)
```

is a Spin(10)-invariant cubic contraction, i.e. a renormalizable `210_H^3`
superpotential term.

Restrict to the color-adjoint branch

```text
Lambda^4 R^6_C ~= Lambda^2 R^6_C ~= (15,1,1),
Phi = *_6 A.
```

For canonical

```text
A = a e12 + b e34 + c e56,
Pf(A) = a b c,
```

the direct contraction gives

```text
I3(*_6 A) = 6 Pf(A).
```

Checks:

```text
(a,b,c)=(1,1,1):       Pf=1,  I3=6,   I3/Pf=6
(a,b,c)=(2,-0.5,3):    Pf=-3, I3=-18, I3/Pf=6
(a,b,c)=(4,1,-2):      Pf=-8, I3=-48, I3/Pf=6

200 random color two-forms:
  mean I3/Pf = 6
  std  I3/Pf = 8.298e-15
```

Normalized local superpotential:

```text
W_loc = 1/2 ||A||^2 - I3(*_6 A)/6.
```

At `A0=e12+e34+e56`, numerical Hessian:

```text
gradient norm = 9.615e-13
eigenvalue clusters = (-1)^1, 0^6, 2^8
max deviation = 1.215e-8
```

Interpretation: the full Spin(10) tensor route exists at the cubic-invariant
level.  The previous PS EFT cubic can be obtained from `210_H^3` restricted to
the `(15,1,1)` branch.  The remaining hard problem is now full multi-field
vacuum alignment: combine this `210_H^3` branch with the existing
`54_H/210_H` projector, finite `45` mediators, Yukawa/flavor sector, and
threshold spectrum without reintroducing unwanted light fragments.

Action items:

```text
1. DONE at PS stage: build explicit SU(4)_C adjoint superpotential.
2. DONE at PS stage: compute the full 15 Hessian and identify 8+3+3bar+1.
3. DONE at PS stage: prove the colored pair is exactly eaten, so Delta_G=0.
4. DONE/no-go: a primitive one-field Spin(10) 45^3 source does not exist.
5. DONE: explicit 210_H^3 invariant matches to the PS cubic on (15,1,1).
6. PARTIAL/DONE at block-operator level: solve the self-consistent
   F54/D210/P_X/I3(210) vacuum-alignment mass matrix, including finite
   mediator thresholds and corrected two-loop/proton replay.  Remaining:
   promote this block-operator alignment to the full component-level
   54_H+210_H+mediator scalar/superpotential vacuum problem.
7. NEXT: decide whether the radial singlet mass M_1=M_Sigma8/2 needs a
   threshold uncertainty entry, even though it is SM neutral at one loop.
```

Follow-up result, 2026-05-06:

Implemented the requested block-operator vacuum-alignment solver:

```text
code/solve_spin10_vacuum_alignment.py
output/spin10_vacuum_alignment/spin10_vacuum_alignment_report.md
output/spin10_vacuum_alignment/spin10_vacuum_alignment_summary.json
output/spin10_vacuum_alignment/spin10_vacuum_alignment_replay.csv
output/spin10_vacuum_alignment/benchmark_R50.json
output/spin10_vacuum_alignment/benchmark_R200.json
```

Mathematical setup:

```text
P_X     = -(9/25)(F54-2)(F54+4/3),
P_color =  (9/50)(F54-2)(F54-1/3),
P_L     =  (D210^2+D210)/2,
P_R     =  (D210^2-D210)/2,
I3(210) = Tr_{Lambda^2}(D_Phi^3), with I3(*_6 A)=6 Pf(A).
```

For a sector with Clebsches `(f,d)=(F54,D210)`, the exact finite-mediator block
mass matrix is

```text
M_fd/MG = [[mu+a54 f+b210 d, eta(f-2), eta(f+4/3)],
          [eta(f-2), 0, R],
          [eta(f+4/3), R, 0]].
```

New element: the script solves a fixed point, not a one-way replay.  For each
`R`, it iterates

```text
(kappa_3,kappa_8) -> mediator mass matrix -> Delta_med
                   -> corrected two-loop/proton replay -> (kappa_3,kappa_8)
```

until convergence.  This removes the mismatch between the old mediator target
spectrum and the 4pi-corrected RGE/proton branch.

Self-consistent corrected-cache results:

```text
R=20:  ||P Delta_med||=5.082702e-03, safe=148,
       M_Sigma3=9.281797e14 GeV, M_Sigma8=3.859002e15 GeV,
       tau_d6=3.549234e35 yr.

R=50:  ||P Delta_med||=2.017822e-03, safe=148,
       M_Sigma3=9.234854e14 GeV, M_Sigma8=3.878242e15 GeV,
       tau_d6=4.026130e35 yr.

R=100: ||P Delta_med||=1.006044e-03, safe=146,
       M_Sigma3=9.219559e14 GeV, M_Sigma8=3.884653e15 GeV,
       tau_d6=4.408812e35 yr.

R=200: ||P Delta_med||=5.022739e-04, safe=146,
       M_Sigma3=9.211965e14 GeV, M_Sigma8=3.887852e15 GeV,
       tau_d6=4.809711e35 yr,
       tau_d5(S_T=1e-5)=2.911503e34 yr.
```

Recommended benchmark card:

```text
R=200,
mu=0.557987598,
a54=0.004787379,
b210=-0.436863941,
eta=3.991775659,
schur_delta_X=0.442618692.
```

The branch survives the self-consistent finite-mediator threshold.  It remains
an intermediate-Sigma3 branch, not a factor-two single-scale branch.

Goldstone/lifting verification:

```text
(15,1,1) -> (8,1,0) + (3,1,2/3) + (bar3,1,-2/3) + (1,1,0).
```

The benchmark keeps only `H_C`, `Sigma_3`, and the octet component of
`Sigma_8` active in the one-loop chiral threshold.  `Sigma_R` and `X_622` are at
`M_G`; the colored pair is eaten by the broken SU(4)_C vector multiplets; the
radial singlet has `M_1=M_Sigma8/2` but zero one-loop SM gauge beta vector.

At the R=200 fixed point, leaving the colored pair at `M_Sigma8` would add a
projected threshold norm

```text
||P Delta_colored_pair||_2 = 1.090376e-01,
```

far above the finite mediator non-universal threshold
`5.022739e-04`.  Thus the largest remaining risk is exactly the full
component-level vacuum alignment proving that the colored pair is eaten/lifted
and that `(6,2,2)` is not re-released.

Next verification task:

```text
1. Promote the block-operator card to an explicit component superpotential:
   W_54 + W_210 + W_med + W_PS-breaking.
2. Compute F-flatness/D-flatness for the combined vev ansatz.
3. Diagonalize the component Hessian at least on the coupled subspace
   {Sigma_L,Sigma_R,Sigma8,X_622,color Goldstones, radial singlet}.
4. Compare the component eigenvalues against benchmark_R200.json.
5. Only after this succeeds, decide whether the radial singlet uncertainty
   deserves a two-loop/gravity-sector note; at one-loop SM gauge matching it is
   currently harmless.
```

Follow-up result, 2026-05-06:

Implemented the first component-Hessian upgrade of `benchmark_R200.json`:

```text
code/verify_spin10_component_hessian.py
output/spin10_component_hessian/spin10_component_hessian_report.md
output/spin10_component_hessian/spin10_component_hessian_summary.json
output/spin10_component_hessian/spin10_component_hessian_blocks.csv
```

This is not yet a full 210-component vacuum solver.  It is the next explicit
component check on the relevant subspace.  The non-color components use the
finite mediator Hessian

```text
M_fd/MG = [[mu+a54 f+b210 d, eta(f-2), eta(f+4/3)],
          [eta(f-2), 0, R],
          [eta(f+4/3), R, 0]].
```

The color `(15,1,1)` component uses the verified `210_H^3` Hessian pattern

```text
radial singlet: -1,
colored Goldstone pair: 0^6,
color octet: +2^8.
```

Benchmark_R200 component eigenvalue verification:

```text
Sigma_L:              mult=3,  light kappa=0.130122468746, rel_err=1.09e-13
Sigma_R:              mult=3,  light kappa=1.000000000000, rel_err=4.13e-14
X_622:                mult=24, light kappa=1.000000000000, rel_err=0
Sigma_8 octet:        mult=8,  light kappa=0.549173662010, rel_err=2.65e-14
colored Goldstone:    mult=6,  light kappa=8.98e-17,       numerical zero/eaten
radial singlet:       mult=1,  light kappa=-0.274586837847,
                      |M|-M_Sigma8/2 relative error=2.49e-08
matching residual:    6.153e-15
```

Threshold audit:

```text
colored pair if unlocked at M_Sigma8: ||P Delta|| = 1.090376e-01
X_622 if unlifted at M_Sigma8:       ||P Delta|| = 1.357954e-01
radial singlet one-loop threshold:   0
```

Conclusion: the R=200 benchmark is reproduced by an explicit component Hessian
on the relevant subspace.  The two dangerous non-SM blocks are genuinely sharp:
`X_622` must remain at `M_G`, and the colored pair must remain an eaten
Goldstone.  The radial singlet is harmless for one-loop gauge matching.

Updated next task:

```text
1. Write the single combined superpotential
   W = W_54 + W_210 + W_med + W_PS + W_10/126
   whose second derivatives reproduce this component Hessian.
2. Compute F-flatness and D-flatness for the combined vev ansatz.
3. Include possible off-block mixings between the 210 color branch and the
   54/mediator sector; prove they vanish by projectors or quantify them.
4. If off-block mixing is nonzero, rerun the component Hessian and RGE/proton
   replay with the mixed eigenvalues.
```

Follow-up result, 2026-05-06:

Implemented the local combined-superpotential flatness/off-block check:

```text
code/verify_combined_superpotential_flatness.py
output/combined_superpotential_flatness/combined_superpotential_flatness_report.md
output/combined_superpotential_flatness/combined_superpotential_flatness_summary.json
```

The restricted local superpotential on the relevant Spin(10)/Pati-Salam
component subspace is

```text
W/MG = W54 + W210 + Wmed + WPS

W54+W210 = 1/2 <A, (1-P_C)(mu+a54 F54+b210 D210) A>
WPS      = m_C [1/2 ||A_C||^2 - I3(*_6 A_C)/6],
           expanded around A0=e12+e34+e56
Wmed     = R <B,C> + eta <A,(F54-2)B> + eta <A,(F54+4/3)C>

P_C = (9/50)(F54-2)(F54-1/3)
P_X = -(9/25)(F54-2)(F54+4/3)
```

Nontrivial new element:

```text
Use P_C to route the full (15,1,1) color branch to WPS, so that the
210_H^3/Pati-Salam cubic controls octet/Goldstone/radial splitting and the
linear 54/210 Clebsch mass controls only the P_C-orthogonal sectors.
```

This prevents double-counting the color mass and makes the off-block Hessian
block-diagonal by projector algebra.

Second nontrivial element:

```text
Schur-flat mediator counter-vev:
  B0/A0 = -eta(F54+4/3)/R,
  C0/A0 = -eta(F54-2)/R.
```

On the color branch `F54=-4/3`, so

```text
B0/A0 = 0,
C0/A0 = 0.066529594317.
```

This cancels the `Wmed` tadpole without changing the Hessian.  Because the
counter-vev is parallel to the color adjoint vacuum, it preserves D-flatness.

Numerical checks:

```text
PS 210H3 gradient norm = 9.615e-13
F_B residual / A0 = 0
F_C residual / A0 = 0
F_A residual / A0 = 0
combined dimensionless F norm = 2.652e-13

adjoint color D proxy = 0
54H diagonal D proxy = 0
210H weak-volume D proxy = 0
210H color-fourform D proxy = 0
max D proxy = 0

restricted Hessian dimension = 135
max abs off-block entry = 0
off-block Frobenius norm = 0
all pass flags = true
```

Current conclusion:

The combined local `W54+W210+Wmed+WPS` card is F-flat, D-flat, and exactly
off-block diagonal on the relevant component subspace, while reproducing the
R=200 component Hessian.  The remaining caveat is now sharply formulated:
this is a local/projector-aligned superpotential card.  A final full model
must show that the same projector orthogonality follows from an untruncated
renormalizable Spin(10) Higgs sector rather than being imposed as the local
alignment basis.

Updated next task:

```text
1. Build the untruncated invariant list for 54_H, 210_H, 45 mediator fields,
   and the PS-breaking color branch.
2. Identify which invariant contractions could mix P_C and 1-P_C sectors.
3. Prove such contractions vanish by representation selection rules, or add
   symmetry/driving fields that forbid them.
4. If any allowed off-block invariant remains, include its coefficient in the
   Hessian and scan the maximum allowed size before X_622 or the colored
   pair becomes dangerous.
```

Follow-up result, 2026-05-06:

Implemented the untruncated renormalizable invariant audit:

```text
code/audit_untruncated_spin10_invariants.py
output/untruncated_spin10_invariants/untruncated_spin10_invariant_audit_report.md
output/untruncated_spin10_invariants/untruncated_spin10_invariant_audit_summary.json
output/untruncated_spin10_invariants/untruncated_spin10_invariant_table.csv
```

Main representation-theory classification:

```text
Forbidden/zero for the present one-copy 54_H and 210_H sector:
  54_H^2 210_H,
  45_i 54_H^2,
  45_i 210_H^2,
  45_i 54_H 210_H,
  45_i^3 and 45_i^2 45_j.

Allowed Higgs vacuum terms:
  54_H^3,
  54_H 210_H^2,
  210_H^3.

Allowed mediator-sector terms that must be scanned unless a symmetry forbids
them:
  m_ij 45_i 45_j,
  sigma_ij 54_H 45_i 45_j,
  tau_ij 210_H 45_i 45_j,
  rho_ABC Tr A[B,C].
```

Numerical sanity checks with random explicit tensors:

```text
Allowed examples:
  Tr S^3              = -6.749921e+01
  S Phi Phi          =  1.409044e+02
  Phi^3 metric       =  4.052987e+02
  S X Y              =  8.792827e+00
  Phi X Y            =  2.094391e+01
  f X Y Z            = -1.486964e+00

Zero examples:
  Tr A^3             = -3.774758e-15
  f A A B            = -7.105427e-15
  A S S moment       = -1.332268e-14
  A Phi Phi moment   =  4.973799e-14
  S S Phi metric     =  9.880985e-15
  A S Phi attempt    =  7.105427e-15
```

Aligned broad-PS off-block audit:

```text
operator                broad off max    spectral norm
F54                     0                 2
D210 weak volume        0                 1
D210 color fourform     0                 2
ad_A0 from fABC         0                 2
```

Conclusion:

The aligned 54/210 Clebsch operators do not mix
`(15,1,1)`, `(6,2,2)`, `(1,3,1)`, `(1,1,3)` broad PS blocks.  The genuinely
dangerous term is not a forbidden off-block contraction, but the allowed
distinct-adjoint invariant

```text
rho_ABC Tr A[B,C].
```

For the aligned color vev it preserves broad PS blocks, but it is an O(1)
component-level deformation with `||ad_A0||=2`.  Using the finite mediator
threshold residual only as a scale marker,

```text
rho * ||ad_A0|| < ||P Delta_med|| = 5.022738709841473e-4
```

would require roughly

```text
rho < 2.511369e-4.
```

This is not a physical exclusion bound, but it shows that `rho_ABC` cannot be
ignored.  The next implementation step should be one of:

```text
Route A: full scan
  Extend the component Hessian by m_ij, sigma_ij, tau_ij, rho_ABC and rerun
  the RGE/proton scan.

Route B: symmetry route
  Add a mediator grading/R-symmetry or driving-field sector that forbids
  rho_ABC and the unwanted entries of m_ij, sigma_ij, tau_ij while preserving
  the projector-generating AA, BC, AB, AC terms.
```

Follow-up result, 2026-05-06:

Implemented Route A:

```text
code/scan_untruncated_invariant_deformations.py
output/untruncated_invariant_deformation_scan/untruncated_invariant_deformation_scan_report.md
output/untruncated_invariant_deformation_scan/untruncated_invariant_deformation_scan_summary.json
output/untruncated_invariant_deformation_scan/untruncated_invariant_deformation_scan.csv
```

The deformation model is

```text
delta W = 1/2 X_i (dm_ij + ds_ij F54 + dt_ij H210) X_j
        + rho Tr A[B,C],
        X_i=(A,B,C).
```

Here `dm,ds,dt` are real symmetric `3x3` matrices with Frobenius-normalized
random directions, and `H210` includes both the weak-volume 210 operator and
the color-fourform PS Hessian with eigenvalues `-1,0,2` on
radial/Goldstone/octet.

Safety definition:

```text
max |kappa_Goldstone| < 1e-6,
||P Delta_X622|| < 5.0227387e-4,
fine-representation off-block max < 1e-10,
and at least one cached corrected two-loop point has matching residual < 1e-3
with d=6 and d=5 proton bounds satisfied.
```

Baseline passes:

```text
kappa3_eff = 0.13012246874638864
kappa8_eff = 0.5491736600061554
Goldstone max |kappa| = 3.076071e-14
X622 projected l2 = 2.806861e-15
fine off-block max = 4.574119e-14
heavy threshold projected l2 = 5.022738709841473e-4
RGE best residual = 5.652822e-10
safe corrected points = 4
tau_d6 = 4.809711e35 yr
tau_d5(S_T=1e-5) = 2.911503e34 yr
```

Random ensemble result, largest amplitude with safe fraction at least `1/2`:

```text
all coefficients together: 1e-6
m_ij only:                 3e-6
sigma_ij only:             1e-6
tau_ij only:               1e-3
rho_ABC only:              >= 1
```

Interpretation:

```text
1. m_ij and sigma_ij are the most dangerous because they directly give the
   colored Goldstone pair a small mass.  The Goldstone criterion fails around
   1e-6--1e-5 unless these entries are correlated or forbidden.
2. tau_ij is much safer, with typical safety through 1e-3 and first median
   RGE failure near 3e-3.  This is because the 210 perturbation preserves the
   protected X622 and Goldstone structure more strongly in the aligned basis.
3. rho_ABC is surprisingly harmless in this aligned Hessian scan through
   O(1): it is an antisymmetric adjoint action inside nearly complete mediator
   subspaces, so it barely changes the protected light eigenvalues or the
   non-universal heavy threshold.  The earlier norm estimate rho<2.5e-4 was
   only a conservative scale marker, not a real bound.
```

Updated next task:

```text
Build Route B using the numerical lesson from Route A:
  impose a minimal mediator grading/R-symmetry that forbids generic m_ij and
  sigma_ij entries at the 1e-6 level, keeps the projector-generating AA, BC,
  AB, AC structure, allows or harmlessly tolerates rho_ABC, and optionally
  allows tau_ij up to about 1e-3.

Then rerun the component Hessian and RGE/proton scan with this symmetry-reduced
coefficient set.
```

Follow-up result, 2026-05-06:

Implemented Route B:

```text
code/verify_routeB_mediator_grading.py
output/routeB_mediator_grading/routeB_mediator_grading_report.md
output/routeB_mediator_grading/routeB_mediator_grading_summary.json
output/routeB_mediator_grading/routeB_reduced_scan.csv
```

First result: there is a field-charge no-go.  A brute-force search over
`Z_N`, `N <= 16`, found no abelian charge assignment that allows

```text
AA, AB, AC, BC, SAA, SAB, SAC
```

while forbidding

```text
BB, CC, SBB, SCC, SBC.
```

Analytic reason:

```text
AA, AB, AC, BC imply q_A=q_B=q_C in any additive abelian grading.
Then SAA, SAB, SAC imply that SBB, SBC, SCC are allowed too.
```

Therefore a single unlabelled `54_H` portal plus abelian charges cannot
protect the projector card.  This is a useful no-go: the model needs
endpoint-labeled mediator links or an equivalent driving-field sector.

Route-B selection rule:

```text
U(1)_M mediator-node grading:
  q(A_45)=0,
  q(B_45)=+1,
  q(C_45)=-1.

Endpoint-labeled portals:
  L_AA^(1,54,210): 0,
  L_AB^(1,54):    -1,
  L_AC^(1,54):    +1,
  L_BC^(1):        0,
  epsilon_tau:     small 210 breaking, |epsilon_tau| <= 1e-3.
```

Allowed projector terms:

```text
L_AA^(1) A A,
L_AA^(54) A A,
L_AA^(210) A A,
L_AB^(1) A B,
L_AB^(54) A B,
L_AC^(1) A C,
L_AC^(54) A C,
L_BC^(1) B C,
rho A[B,C].
```

Forbidden/spurion-suppressed:

```text
BB, CC,
generic BC 54 portal,
generic SBB, SBC, SCC,
generic m_ij and sigma_ij entries not in the link list.
```

Reduced Route-A replay:

```text
scenario              eps_m   eps_sigma  eps_tau  |rho|  safe_fraction
exact_link_grading    0       0          0        1      1.000
routeB_target         1e-6    1e-6       1e-3    1      0.677
routeB_conservative   3e-7    3e-7       1e-3    1      0.802
tau_conservative      1e-6    1e-6       3e-4    1      0.802
tau_stress            0       0          3e-3    1      0.365
```

Median diagnostics:

```text
exact link grading:
  Goldstone kappa = 3.36e-14,
  X622 projected l2 = 2.49e-6,
  RGE residual = 1.55e-5,
  tau_d6 = 4.81e35 yr,
  tau_d5 = 2.91e34 yr.

routeB target:
  Goldstone kappa = 4.93e-7,
  X622 projected l2 = 2.49e-6,
  RGE residual = 3.16e-4,
  tau_d6 = 4.81e35 yr,
  tau_d5 = 2.91e34 yr.
```

Conclusion:

The minimal viable protection is not a pure charge assignment on `A,B,C`, but
a link-graded mediator sector with aligned endpoint-labeled `54` portals.
The exact link-graded card is safe.  The practical benchmark should use either

```text
eps_m=eps_sigma=3e-7, eps_tau=1e-3, |rho|=1
```

or

```text
eps_m=eps_sigma=1e-6, eps_tau=3e-4, |rho|=1
```

depending on whether the manuscript wants stronger Goldstone safety or stronger
RGE margin.  The stress point `eps_tau=3e-3` fails often, confirming the Route-A
tolerance estimate.

Completed link-driving upgrade:

```text
script:
  code/verify_link_driving_superpotential.py
outputs:
  output/routeB_link_driving/routeB_link_driving_summary.json
  output/routeB_link_driving/routeB_link_driving_report.md

W_link:
  1/2 <A,(1-P_C)(l_AA^1 + l_AA^54 F54 + l_AA^210 D210)A>
  + <A,(l_AB^1 + l_AB^54 F54)B>
  + <A,(l_AC^1 + l_AC^54 F54)C>
  + l_BC^1 <B,C> + W_PS,color[A_C] + rho Tr A[B,C].

W_drive:
  X_AA^1(l_AA^1-mu) + X_AA^54(l_AA^54-a54)
  + X_AA^210(l_AA^210-b210)
  + X_AB^n(l_AB^54-eta) + X_AB^r(l_AB^1+2 l_AB^54)
  + X_AC^n(l_AC^54-eta) + X_AC^r(l_AC^1-(4/3)l_AC^54)
  + X_BC^1(l_BC^1-R)
  + perpendicular alignment drivers for the 54 and 210 directions.

fixed ratios:
  l_AB^1/l_AB^54 = -2,
  l_AC^1/l_AC^54 = 4/3.

numerical verification:
  F norm = 0.000e+00,
  D proxy max = 0.000e+00,
  max |H_link-H_routeB| = 0.000e+00,
  Goldstone kappa = 5.65e-14,
  X622 projected l2 = 2.49e-6,
  RGE residual = 1.55e-5,
  safe corrected points = 4,
  tau_d6 = 4.81e35 yr,
  tau_d5(S_T=1e-5) = 2.91e34 yr.
```

Updated next task:

```text
Promote the driving fields themselves into explicit Spin(10) representation
content and check whether their orthogonal alignment modes introduce additional
non-singlet thresholds.  The minimal next card should list the drive-sector
mass eigenvalues, show which drive fields are pure Lagrange multipliers, and
repeat the threshold/proton scan if any drive multiplet remains below M_G.
```

Completed drive-sector representation audit:

```text
script:
  code/audit_drive_sector_spectrum.py
outputs:
  output/drive_sector_spectrum/drive_sector_spectrum_summary.json
  output/drive_sector_spectrum/drive_sector_spectrum_report.md
  output/drive_sector_spectrum/drive_sector_threshold_replay.csv

representation lift:
  AA, AB, AC 54 links:
    L_e^54 + D_e^54 in full 54 + 54_D.
    54 -> (1,1,1) + (20',1,1) + (6,2,2) + (1,3,3).
    Orthogonal components per field: 53.
    Pair beta vector: (24,24,24).

  AA 210 link:
    L_AA^210 + D_AA^210 in full 210 + 210_D.
    210 -> (1,1,1) + (15,1,1) + (15,3,1) + (15,1,3)
           + (6,2,2) + (10,2,2) + (10bar,2,2).
    Orthogonal components per field: 209.
    Pair beta vector: (112,112,112).

  Ratio/norm drivers:
    Spin(10) singlets, no gauge threshold.

orthogonal Hessian:
  W_perp = kappa M_G <D_perp, delta L_perp>,
  M_perp/M_G = [[0,kappa],[kappa,0]],
  masses = |kappa| M_G.

threshold theorem:
  Delta_drive =
    [24 sum_e log(1/kappa_54,e) + 112 log(1/kappa_210)]/(2 pi) * (1,1,1).
  Therefore P Delta_drive = 0 exactly: no non-universal drive threshold enters
  the gauge-matching residual.

RGE/proton replay:
  canonical kappa54=kappa210=1:
    projected threshold = 0,
    safe points = 4,
    alphaG_inv = 41.3631,
    residual = 5.653e-10,
    tau_d6 = 4.810e35 yr,
    tau_d5 = 2.912e34 yr.

  common kappa=0.5:
    universal Delta = 20.2985,
    projected l2 = 6.15e-15,
    safe points = 4,
    alphaG_inv = 21.0646,
    tau_d6 = 1.247e35 yr.

  common kappa=0.3:
    universal Delta = 35.2578,
    projected l2 = 6.15e-15,
    alphaG_inv = 6.105,
    tau_d6 = 1.048e34 yr.

lower-bound estimate for common drive-pair mass:
  kappa > 0.289 for alphaG_inv > 5,
  kappa > 0.299 for tau_d6 > 1e34 yr,
  kappa > 0.334 for tau_d6 > 2.4e34 yr.
```

Updated next task:

```text
The drive-sector orthogonal modes are now safe if they remain complete
Spin(10) pairs near M_G.  The next hard step is to make the alignment source
fields S_54 and Omega_210 dynamical inside the same W_54+W_210+W_med+W_PS
superpotential, then check whether mixing between the original breaking fields
and the link/driver pairs keeps the complete-pair degeneracy or induces
non-complete splittings.
```

Completed dynamic source-mixing audit:

```text
script:
  code/audit_dynamic_source_mixing.py
outputs:
  output/dynamic_source_mixing/dynamic_source_mixing_summary.json
  output/dynamic_source_mixing/dynamic_source_mixing_report.md
  output/dynamic_source_mixing/dynamic_source_mixing_replay.csv

local source-link-driver block:
  W_R = 1/2 m_R S_R^2 + sum_e kappa_e D_e (L_e - ell_e S_R).

benchmark link coefficients:
  54:  ell = [a54, eta, eta]
       = [0.0047873793, 3.9917756590, 3.9917756590],
       fields = [S54, L_AA^54, L_AB^54, L_AC^54,
                 D_AA^54, D_AB^54, D_AC^54].

  210: ell = [b210] = [-0.4368639411],
       fields = [Omega210, L_AA^210, D_AA^210].

fragment beta checks:
  54 fragments:
    (1,1,1), (20',1,1), (6,2,2), (1,3,3)
    sum beta = (12,12,12), P sum = 0.

  210 fragments:
    (1,1,1), (15,1,1), (15,3,1), (15,1,3),
    (6,2,2), (10,2,2)+(10bar,2,2)
    sum beta = (56,56,56), P sum = 1.63e-14.

degenerate-source theorem:
  If m_54 and m_210 are Spin(10)-scalar masses, every PS fragment sees the
  same finite source-link-driver mass matrix.  Therefore every eigenvalue is
  a complete Spin(10) multiplet and the non-universal threshold vanishes.

numerical degenerate-source check:
  54 fragment eigenvalue spread  = 0.000e+00,
  210 fragment eigenvalue spread = 0.000e+00,
  projected threshold l2         = 1.366e-30,
  universal Delta                = 7.987e-15,
  safe points                    = 4,
  alphaG_inv                     = 41.3631,
  tau_d6                         = 4.810e35 yr.

fragment-split stress:
  eps54=eps210=1e-4:
    projected l2 = 4.288e-4,
    safe points = 4,
    residual = 2.789e-4.

  eps54=eps210=1e-3:
    projected l2 = 4.288e-3,
    safe points = 0,
    residual = 1.719e-3,
    d5 lifetime at picked point = 6.52e33 yr.

  eps54=1e-3, eps210=0:
    projected l2 = 7.928e-4,
    safe points = 4.

  eps54=0, eps210=1e-3:
    projected l2 = 3.721e-3,
    safe points = 0.
```

Interpretation:

```text
Making S_54 and Omega_210 dynamical does not by itself break complete-multiplet
degeneracy.  The dangerous object is a non-scalar PS-fragment source Hessian,
especially in the 210 sector.  A conservative model-building requirement is
that any 210 source-fragment splitting induced by W_54+W_210 be absent by
symmetry or suppressed to roughly <= 1e-4 before threshold/proton replay.
```

Updated next task:

```text
Build the source-stabilization sector explicitly.  The target is a
Spin(10)-scalar source Hessian for the non-Goldstone S_54 and Omega_210
orthogonal directions, or a symmetry/projector argument proving that the
dangerous 210 PS-fragment splittings vanish to <= 1e-4.  If this cannot be
done with renormalizable invariants alone, record the obstruction and test the
minimal extra singlet/R-symmetry completion.
```

Completed source-stabilization sector:

```text
script:
  code/construct_source_stabilization_sector.py
outputs:
  output/source_stabilization/source_stabilization_summary.json
  output/source_stabilization/source_stabilization_report.md
  output/source_stabilization/source_split_threshold_scan.csv

nontrivial construction:
  1. Source-field doubling:
     Omega_W is the weak-volume 210 source used by D210.
     Phi_C is a distinct color 210_C field carrying the allowed Phi_C^3
     interaction that produced the PS color cubic.

  2. R-sequestering plus source parity:
     W has R=2.
     S_54, Omega_W, L_54, L_210 have R=0 and P_src=-.
     D_54, D_210 have R=2 and P_src=-.
     X_a have R=2 and P_src=+.
     Phi_C has R=2/3 and P_src=+.
     Optional parity spurion zeta has P_src=-.

allowed:
  D_54,e (L_54,e - ell_e S_54),
  D_210 (L_210 - ell_210 Omega_W),
  X_a f_a(S_54, Omega_W),
  Phi_C^3.

forbidden at tree level:
  Omega_W^3,
  S_54 Omega_W^2,
  Omega_W Phi_C^2.

theorem:
  On the F-flat branch with all R=2 drivers zero, the only terms that would
  give Clebsch-dependent source Hessians are forbidden.  Therefore the
  non-Goldstone source-link-driver Hessian is repeated identically on every
  PS fragment of a full 54 or 210.  Its eigenvalues are complete Spin(10)
  multiplets, and P Delta_source = 0.

residual spurion scan:
  m_54(f)  = m_54  (1 + epsilon54 c_f),
  m_210(f) = m_210 (1 + epsilon210 c_f).

  grid points = 441,
  safe fraction = 0.429,
  max |epsilon210| safe on epsilon54=0 axis = 3.0e-4,
  max |epsilon54| safe on epsilon210=0 axis = 3.0e-3,
  all points safe inside |eps54|,|eps210| <= 1e-4: true,
  all points safe inside |eps54|,|eps210| <= 3e-4: true,
  first unsafe by radius:
    eps54=0, eps210=-5e-4,
    residual=1.154e-3,
    tau_d5=1.069e35 yr.

spurion order:
  lambda_src=0.1 requires dangerous 210 source operators to first appear at
  order >= 5 to ensure epsilon210 <= 1e-4.
  lambda_src=0.2 requires order >= 6.
```

Interpretation:

```text
The sequestered source sector solves the 210 source-splitting problem at tree
level.  If UV source parity leaks, the split is no longer a hidden assumption:
it is now an explicit threshold parameter epsilon210.  For the paper benchmark,
quote the conservative condition |epsilon210| <= 1e-4; the scan suggests some
sign-dependent margin up to about 3e-4 on the epsilon54=0 axis, but 1e-4 is
cleaner and robust.
```

Completed combined Route-B/source charge table:

```text
script:
  code/verify_combined_charge_table.py
outputs:
  output/combined_charge_table/combined_charge_table_summary.json
  output/combined_charge_table/combined_charge_table_report.md

final selection rule:
  G_sel = U(1)_R x U(1)_M x Z3_E,
  R(W)=2,
  e(A)=0, e(B)=e(C)=1 mod 3.

key compatibility repairs:
  1. The old source parity is only a source-sector bookkeeping device.
     Taken literally it would forbid W_link terms L_ij X_i X_j, so it is not
     kept as the combined symmetry.

  2. U(1)_M alone cannot distinguish the neutral AA and BC channels.
     Z3_E separates them:
       e(AA)=0, e(BC)=2,
     allowing L_AA AA and L_BC BC while forbidding L_AA BC and L_BC AA.

  3. A single common S_54 cannot couple to AA, AB and AC links while preserving
     U(1)_M.  The consistent completion uses endpoint source copies
     S_AA^54, S_AB^54, S_AC^54 with identical Spin(10) orbit constraints and
     charges q(S_e^54)=q(L_e^54).

  4. Endpoint source copies leave relative 54-orientation moduli unless one adds
     the allowed drivers Y_AB^54 and Y_AC^54:
       W_rel = <Y_AB, S_AB - v_AB S_AA> + <Y_AC, S_AC - v_AC S_AA>.

  5. The Yukawa/Majorana sector is neutral under U(1)_M x Z3_E, with
     R(16)=1 and R(10_H,120_H,overline126_H,overline16_H)=0.
     Therefore 16 16 10_H, 16 16 120_H, 16 16 overline126_H and
     16 16 overline16_H overline16_H / Lambda are all allowed.

verification:
  terms checked = 58,
  failures = 0,
  W_link = 8/8,
  rho_ABC A[B,C] = 1/1,
  W_drive = 14/14,
  W_src = 10/10,
  W_rel = 4/4,
  Yukawa/Majorana = 4/4,
  forbidden generic/wrong/source terms = 13/13,
  forbidden relative cross terms = 4/4.

forbidden examples:
  BB, CC, BC,
  L_AA BC, L_BC AA,
  S_AA BC, S_AB BB, S_AC CC,
  Omega_AA^3, S_AA Omega_AA^2, Omega_AA Phi_C^2.
```

Completed endpoint source-copy Hessian:

```text
script:
  code/check_endpoint_source_copy_hessian.py
outputs:
  output/endpoint_source_copies/endpoint_source_copy_hessian_summary.json
  output/endpoint_source_copies/endpoint_source_copy_hessian_report.md
  output/endpoint_source_copies/endpoint_source_misalignment_scan.csv

source-copy theorem:
  For each endpoint e=AA,AB,AC,
    W_e = 1/2 m_e ||S_e||^2 + kappa_e <D_e, L_e - ell_e S_e>.
  If m_e and kappa_e are Spin(10)-scalar, the 3x3 source-link-driver matrix is
  repeated on every component of a full 54.  Therefore its eigenvalues are
  complete 54 multiplets and the projected one-loop threshold vanishes.

numerical check:
  aligned F54 operator error = 0.000e+00,
  P beta_54 l2 = 0.000e+00,
  P Delta_source l2 = 2.733e-30,
  aligned endpoint offblock max = 4.574e-14,
  aligned safe points = 4,
  aligned RGE residual = 5.653e-10,
  tau_d6 = 4.810e35 yr,
  tau_d5(S_T=1e-5) = 2.912e34 yr.

relative-orientation result:
  The charge table aligns L_e to its own S_e but does not by itself align
  S_AB and S_AC to S_AA.  Writing
    S_AB = S0 + epsilon_AB T_AB,
    S_AC = S0 + epsilon_AC T_AC,
  generic off-block T changes AB/AC from F54 to F54+epsilon F_T.

  first tested structurally unsafe point:
    epsilon_AB = epsilon_AC = 1e-10,
    offblock max = 2.823e-10,
    RGE residual still = 5.653e-10,
    safe_without_relative_driver = false.

interpretation:
  The problem is structural rather than matching-size: nonzero off-block
  entries invalidate the PS fragment threshold basis even if the cached RGE
  residual stays small.  Therefore keep W_rel as a required part of the
  source-copy completion rather than scanning misalignment as a tolerated
  threshold dial.
```

Completed W_rel drive/source-sector spectrum audit:

```text
updated script:
  code/audit_drive_sector_spectrum.py
updated outputs:
  output/drive_sector_spectrum/drive_sector_spectrum_summary.json
  output/drive_sector_spectrum/drive_sector_spectrum_report.md
  output/drive_sector_spectrum/drive_source_relative_eigenvalues.csv
  output/drive_sector_spectrum/drive_source_relative_replay.csv

new representation content:
  Y_AB^54 and Y_AC^54 are full Spin(10) 54 drivers.  Together with
  S_AA,S_AB,S_AC,L_AA,L_AB,L_AC,D_AA,D_AB,D_AC they form an 11-field
  complete-54 source/link/driver block for each 54 component.

local basis:
  (S_AA, S_AB, S_AC, L_AA, L_AB, L_AC, D_AA, D_AB, D_AC, Y_AB, Y_AC).

superpotential:
  W_rel = <Y_AB, S_AB - v_AB S_AA> + <Y_AC, S_AC - v_AC S_AA>.

54 eigenvalues over MG:
  -3.8003, -3.7870, -1.0019, -0.9981, -0.2112,
   0.08450, 0.2641, 0.99999, 1.9673, 4.7342, 4.7486.

checks:
  min |lambda_54| = 8.450e-2,
  log sum 54 = -1.098612,
  log sum 210 source block = -5.55e-17,
  projected Delta_source+rel l2 = 3.846e-16,
  universal Delta_source+rel = -2.098195.

RGE/proton replay with this universal source block:
  safe points = 4,
  alphaG^-1 = 43.4613,
  residual = 5.653e-10,
  tau_d6 = 5.310e35 yr,
  tau_d5(S_T=1e-5) = 2.912e34 yr.

conclusion:
  The Y_AB,Y_AC non-singlet modes are lifted and repeat over complete 54
  multiplets.  They add no non-universal threshold.  The only effect in the
  canonical normalization is a harmless universal alpha_G shift.
```

Completed source self-coupling / kinetic-leakage audit after W_rel:

```text
new script:
  code/audit_source_leakage_after_wrel.py
new outputs:
  output/source_leakage_after_wrel/source_leakage_after_wrel_summary.json
  output/source_leakage_after_wrel/source_leakage_after_wrel_report.md
  output/source_leakage_after_wrel/source_self_leakage_scan.csv
  output/source_leakage_after_wrel/source_kinetic_leakage_scan.csv

selection-rule result:
  Scalar source masses and scalar source/link Kähler metrics only shift
  complete 54/210 multiplets.  Neutral weak-source self-couplings such as
  Omega_AA^3 and S_AA Omega_AA^2 are forbidden in the strict source sector;
  if an R=2 neutral leakage spurion is admitted, they must be counted as
  epsilon210 threshold deformations.  Phi_C^3 remains safe because Phi_C is
  sequestered from the weak-volume source.

mathematical check:
  For a fragment f,
    M_can(f) = Z_f^{-T/2} M_hol(f) Z_f^{-1/2},
  hence
    sum log |lambda_can(f)| = log |det M_hol(f)| - log det Z_f.
  Therefore a non-scalar Kähler metric regenerates PS-fragment thresholds even
  when W_rel has removed the relative-orientation modulus.

scalar W_rel-completed baseline:
  P Delta_source+rel l2 = 4.710e-16,
  universal Delta = -2.098195,
  log-sum spread 54 = 0,
  log-sum spread 210 = 0.

holomorphic source-Hessian leakage scan:
  points = 441,
  safe fraction = 0.429,
  max |epsilon54| safe at epsilon210=0 = 3.0e-3,
  max |epsilon210| safe at epsilon54=0 = 3.0e-4,
  all points safe inside |epsilon54|,|epsilon210| <= 3e-4,
  first unsafe by radius: epsilon54=0, epsilon210=-5e-4,
  residual = 1.154e-3.

kinetic leakage scan:
  points = 441,
  safe fraction = 0.571,
  max |zeta54| safe at zeta210=0 = 3.0e-3,
  max |zeta210| safe at zeta54=0 = 7.5e-4,
  all points safe inside |zeta54|,|zeta210| <= 3e-4,
  first unsafe by radius: zeta54=0, zeta210=-7.5e-4,
  tau_d5(S_T=1e-5) = 6.80e33 yr.

conclusion:
  W_rel itself does not regenerate non-scalar PS-fragment splitting.  The
  remaining source/Higgs-sector assumption is that source self-couplings and
  source/link kinetic terms are Spin(10)-scalar, or else their 210 leakage is
  suppressed to roughly 1e-4--1e-3.  If not, epsilon210/zeta210 must become
  explicit nuisance threshold parameters in the final fit.
```

Updated next task:

```text
Promote the source-sequestering assumption into a theorem-level symmetry or
spurion statement inside the full Higgs/Kähler sector.  Decide whether the
final paper treats epsilon210 and zeta210 as:
  (A) forbidden/suppressed by a UV selection rule, with the first dangerous
      operator appearing at order >=5 for lambda_src=0.1 or >=6 for
      lambda_src=0.2; or
  (B) explicit nuisance thresholds to be scanned together with the mediator,
      RGE, and proton-decay fit.
```

Heartbeat status 2026-05-06 22:34 Asia/Taipei:

```text
overall status:
  The first-principles GUT program is not complete.  The current manuscript is
  still a conditional Spin(10) EFT branch rather than an unconditional
  derivation from PSLT/Another Physics.  The latest source/Higgs-sector audit
  did, however, close one concrete loophole: W_rel itself is safe, while
  non-scalar source self-couplings or Kähler leakage must be suppressed or
  scanned.

current obstacle:
  The full Higgs/Kähler sector still lacks a theorem-level source-sequestering
  statement.  In present notation the dangerous directions are the non-scalar
  210 leakages epsilon210 and zeta210 generated by Omega_AA^3,
  S_AA Omega_AA^2, or fragment-dependent source/link kinetic metrics.

next nontrivial idea:
  Treat source sequestering as a spurion filtration theorem.  Assign every
  source-only or Kähler leakage operator a filtration degree n under the
  combined Route-B/source grading, and prove that the first non-scalar 210
  invariant occurs only at n >= 5 for lambda_src=0.1, or n >= 6 for
  lambda_src=0.2.  This turns the numerical condition
  |epsilon210|, |zeta210| roughly <= 1e-4 into an operator-order statement.

verification plan:
  1. Extend verify_combined_charge_table.py into a leakage-filtration enumerator
     for holomorphic W terms and Kähler terms up to degree 6.
  2. For every allowed leakage monomial, classify its Spin(10) tensor content:
     scalar complete-multiplet, non-scalar 54, non-scalar 210, or mixed.
  3. Feed the induced epsilon54/epsilon210/zeta54/zeta210 estimates into
     audit_source_leakage_after_wrel.py and confirm the RGE/proton-safe branch.
  4. If any degree <= 4 non-scalar 210 operator survives, abandon pure
     sequestering and promote epsilon210/zeta210 to explicit nuisance
     threshold parameters in the final global scan.
```

Heartbeat status 2026-05-06 23:10 Asia/Taipei:

```text
overall status:
  Still incomplete.  No new theorem-level derivation has closed the remaining
  source-sequestering loophole since the previous heartbeat, and the paper
  remains a conditional Spin(10) EFT rather than an unconditional
  PSLT/Another-Physics derivation.

current obstacle:
  The current bottleneck is unchanged: prove or scan away non-scalar 210
  leakage in the full Higgs/Kähler sector after W_rel.  The concrete danger
  variables remain epsilon210 from holomorphic source self-couplings and
  zeta210 from fragment-dependent source/link kinetic metrics.

next attempted nontrivial idea:
  Build a finite leakage-filtration enumerator.  Instead of merely assigning
  charges by hand, enumerate every holomorphic W monomial and Kähler bilinear
  up to filtration degree 6, then classify its induced Spin(10) tensor as
  scalar, 54-like, 210-like, or mixed.  The desired theorem is that all
  degree <=4 non-scalar 210 leakages are absent in Route-B/source grading.

verification plan:
  Implement code/enumerate_source_leakage_filtration.py using the verified
  charge table; compare the first allowed non-scalar 210 order against the
  numerical requirement from audit_source_leakage_after_wrel.py:
    lambda_src=0.1 requires order >=5,
    lambda_src=0.2 requires order >=6.
  If the enumerator finds a lower-order 210 leakage, update TeX to present the
  result as a nuisance-threshold scan rather than a sequestering theorem.
```

Heartbeat status 2026-05-06 23:41 Asia/Taipei:

```text
overall status:
  Still incomplete.  No new mathematical result has superseded the 23:10
  diagnosis, so no TeX theorem should be added yet.

current obstacle:
  The remaining gap is constructive: the source-sequestering claim must be
  checked by enumeration, not by prose.  The current charge table is verified
  term-by-term, but it has not yet been promoted into an exhaustive
  degree-bounded leakage theorem for W and Kähler terms.

next attempted nontrivial idea:
  Turn the charge table into a graded invariant search:
    monomial charge = sum field charges + spurion charges,
    filtration degree = number of source/leakage spurion insertions,
    dangerous class = monomials whose Hessian or kinetic metric carries a
      non-scalar 210 projector on PS fragments.

verification plan:
  Run a degree <=6 enumeration and require:
    no degree <=4 holomorphic non-scalar 210 source self-coupling,
    no degree <=4 Kähler non-scalar 210 source/link metric,
    all allowed low-degree terms either scalar complete multiplets or already
      included Route-B/projector terms.
  If the check passes, add a source-sequestering proposition to TeX.  If it
  fails, use audit_source_leakage_after_wrel.py to scan the surviving leakage
  coefficient as a nuisance threshold.
```

Referee-style review response recorded 2026-05-07:

```text
overall assessment:
  The review is largely correct.  The new draft is much stronger than the
  earlier one, especially on the corrected 4pi convention, 54/210 projector
  algebra, PS Goldstone locking, and source-stabilization audits.  It is still
  a major-revision working draft rather than a submission-ready paper.

confirmed issues:
  1. Numerical cache mixing:
     paper/gut_framework.tex still contained older Goldstone-locking tolerance
     values |epsilon_G| = 7.074413e-3 at R=50 and 1.764167e-3 at R=200,
     while the self-consistent block-operator cache gives 1.109124e-2 and
     2.760820e-3.  The TeX was patched to the self-consistent values.

  2. UV perturbativity / Landau pole:
     Added code/audit_uv_perturbativity_landau.py.
     For N=1 SUSY Spin(10),
       d alpha_G^{-1}/d log mu = -b10/(2 pi),
       b10 = sum_R T(R) - 24.
     The documented complete drive pairs have
       sum T = 6 T(54) + 2 T(210) = 184,
       b10 = 160.
     This gives
       R=50:  Lambda_LP/MG = 4.42, alpha^{-1}(R MG) = -61.775,
       R=200: Lambda_LP/MG = 5.07, alpha^{-1}(R MG) = -93.557.
     Therefore the high-R mediator/source completion cannot be interpreted as
     a fully propagating perturbative Spin(10) UV completion if those large
     multiplets are active already at MG.

  3. Proton/flavor/reproducibility:
     Existing channel-specific proton replay is useful but remains a toy
     Wilson-coefficient estimate; full p -> K+ nu needs flavor rotations,
     dressing functions, short/long-distance factors, and lattice matrix
     elements.  The seesaw and Yukawa sectors still require full matrices and
     a reproducible CKM/PMNS fit.

new outputs:
  output/uv_perturbativity_landau/uv_landau_audit_summary.json
  output/uv_perturbativity_landau/uv_landau_audit_report.md
  output/uv_perturbativity_landau/uv_landau_audit.csv

paper update:
  Added a UV Perturbativity Above MG subsection with the b10 formula and
  scenario table.  The conclusion is now explicit: complete multiplets are
  projected-threshold silent but not beta-function silent.
```

Updated next task after UV audit:

```text
Resolve the interpretation of the large source/link/driver sector:
  Route U1: declare source/link/driver fields constrained spurions/auxiliary
            fields and state the cutoff EFT range explicitly;
  Route U2: lower R to be below the few-MG Landau scale and redo mediator
            threshold/proton scan;
  Route U3: search for a smaller-representation projector/source realization
            with sum T small enough to reach R=50--200 perturbatively.

The most honest next verification is U3-vs-U1:
  enumerate possible lower-index mediator/projector completions using 45, 16,
  10, and singlet spurions, and compare whether they can reproduce the
  P_X lift without introducing (6,2,2) or colored Goldstone leakage.  If not,
  write the final paper as a cutoff EFT with non-propagating source spurions.
```

Heartbeat status 2026-05-07 00:48 Asia/Taipei:

```text
overall status:
  Still incomplete.  The latest UV Landau audit converted one referee concern
  into a hard theorem-level obstruction: the high-R source/driver completion is
  not a perturbative propagating Spin(10) sector if the large 54/210 multiplets
  are active already at MG.

current obstacle:
  Decide whether the model is a cutoff EFT with constrained/non-propagating
  source spurions, or whether a lower-index projector/source realization can
  replace the 54/210 drive/source completion while preserving the verified
  P_X lift, Sigma_3 threshold, X_622 lifting, and Goldstone lock.

next attempted nontrivial idea:
  Search for a minimal-index projector realization.  The key algebraic target
  is still a polynomial filter equivalent to
    P_X = -(9/25)(F54-2)(F54+4/3),
  but the implementation should try to generate the two roots using smaller
  propagating representations plus spurionic vev insertions, so that sum T is
  small enough to keep Lambda_LP > R MG for R=50--200.

verification plan:
  1. Enumerate candidate mediator/filter blocks built from 45, 16+16bar, 10,
     and singlet spurions, allowing non-propagating source vevs but counting
     only propagating chiral multiplets in b10.
  2. For each candidate, compute effective Clebsch eigenvalues on
     Sigma_L, Sigma_R, Sigma8, and X_622.
  3. Reject candidates that leave X_622 light, move Sigma_R intermediate,
     generate colored-Goldstone mass, or require sum T large enough to put
     Lambda_LP below the mediator scale.
  4. If no candidate passes, update the TeX framing to explicitly state that
     the present construction is a cutoff EFT/spurion completion rather than a
     perturbative Spin(10) UV completion.
```

Heartbeat status 2026-05-07 01:20 Asia/Taipei:

```text
overall status:
  Still incomplete.  The UV perturbativity obstruction remains the active hard
  blocker; no later result has restored a perturbative high-R Spin(10) UV
  interpretation.

current obstacle:
  The verified threshold construction uses large complete 54/210
  source/driver multiplets.  These are harmless for projected one-loop
  threshold differences, but if propagating above MG they drive b10 positive
  enough to Landau-pole below R=50--200.

next attempted nontrivial idea:
  Make the next calculation a falsifiable small-index search rather than a
  prose choice.  Define a candidate score
    score = mismatch(P_X target) + Goldstone leakage penalty
            + X_622 leakage penalty + Landau penalty,
  where the Landau penalty is infinite unless Lambda_LP > R MG for the
  displayed R target.

verification plan:
  Implement a small enumerator for 45-only, 45+singlet, 45+16+16bar, and
  45+10 spurion-assisted filters.  For each, compute the effective diagonal
  filter on {Sigma_L,Sigma_R,Sigma8,X_622}, count propagating Dynkin index,
  and reject any candidate with b10 causing Lambda_LP/MG < R.  If all
  candidates fail, the TeX should be tightened to the cutoff-EFT/spurion
  interpretation rather than leaving this as a maybe.
```

Heartbeat status 2026-05-07 01:55 Asia/Taipei:

```text
overall status:
  Still incomplete.  The current state is a conditional Spin(10) threshold EFT
  with strong local algebraic checks, not a completed first-principles GUT.

current obstacle:
  The UV Landau-pole audit remains decisive for the propagating high-R
  interpretation.  The verified large 54/210 source/driver completion can be
  kept only as a cutoff-EFT/spurion mechanism unless a lower-index realization
  is found.

next attempted nontrivial idea:
  Treat the small-index search as a no-go/proof attempt: test whether the
  required filter can be generated by low-index propagating 45 mediators plus
  non-propagating scalar spurion insertions.  The goal is not merely to fit the
  four target eigenvalues, but to preserve the already verified constraints:
  Sigma_3 intermediate, Sigma_R at MG, X_622 lifted, and colored Goldstone
  exactly locked.

verification plan:
  Build a candidate enumerator with basis functions on the four fragments
  {Sigma_L,Sigma_R,Sigma8,X_622}.  For each candidate:
    - solve least-squares/exact interpolation to the P_X target;
    - compute required propagating Dynkin index and Lambda_LP/MG;
    - reject if the candidate needs non-scalar 54/210 propagating fields;
    - replay the surviving threshold vector through the corrected RGE/proton
      scan.
  If no candidate survives, write the no-go as a proposition and explicitly
  choose the spurion/cutoff-EFT interpretation in the paper.
```

Heartbeat status 2026-05-07 02:26 Asia/Taipei:

```text
overall status:
  Still incomplete.  The active blocker is unchanged: the threshold/model
  algebra is locally strong, but the high-R propagating UV completion is
  obstructed by the one-loop Spin(10) Landau-pole audit.

current obstacle:
  We need a mathematically checked alternative to the large 54/210
  source/driver sector, or else the paper must explicitly state that the
  source/projector sector is a non-propagating spurion/cutoff-EFT completion.

next attempted nontrivial idea:
  Compress the projector into a low-index filter basis.  Instead of starting
  from full 54/210 dynamical tensors, use the already verified four-fragment
  eigenvalue data as the target vector and search over low-index mediator
  blocks whose propagating content is 45-only or 45 plus small reps, with
  source vevs treated as fixed spurions.

verification plan:
  Candidate filters must pass four tests simultaneously:
    1. reproduce P_X on Sigma_L, Sigma_R, Sigma8, X_622 to numerical tolerance;
    2. keep Sigma_R and X_622 lifted near MG while retaining intermediate
       Sigma_3;
    3. leave the colored Goldstone eigenvalue exactly zero/eaten;
    4. satisfy Lambda_LP/MG > R for the intended mediator scale.
  Failure of this enumerator should be recorded as a no-go for perturbative
  high-R propagation and trigger a TeX rewrite toward cutoff-EFT framing.
```

Heartbeat status 2026-05-07 03:00 Asia/Taipei:

```text
overall status:
  Still incomplete, but the UV blocker is now sharper.  A low-index algebraic
  escape exists for the projector filter, yet it requires treating the F54
  direction as a spurion/constrained background rather than using the full
  propagating 54/210 alignment sector.

completed in this heartbeat:
  Added code/enumerate_low_index_projector_filters.py.
  Outputs:
    output/low_index_projector_filters/low_index_projector_filters_summary.json
    output/low_index_projector_filters/low_index_projector_filters_report.md
    output/low_index_projector_filters/low_index_projector_filters.csv

mathematical result:
  On {Sigma_L,Sigma_R,Sigma8,X_622},
    F54 = (2,2,-4/3,1/3),
    P_X = (0,0,0,1).
  The exact polynomial is
    P_X = -(9/25)(F54-2)(F54+4/3)
        = 0.96 + 0.24 F54 - 0.36 F54^2.
  A quadratic F54 spurion filter reproduces the target with max error
  2.22e-16 and propagating sumT=24 (three 45 fields), so b10=0 and R=200
  passes the one-loop Landau audit.

negative result:
  Linear F54 cannot reproduce the target (max error 7.27e-1).
  Including a propagating 210 in the low-index FD basis gives sumT=92,
  b10=68 and fails R=50 and R=200.  The documented large alignment sector
  has sumT=208, b10=184 and fails high-R propagation.

paper update:
  Added a Low-Index Projector Escape subsection after the UV perturbativity
  audit.  The paper now states explicitly that the high-R branch is viable only
  in a constrained-spurion/source-background interpretation unless a genuine
  small-representation alignment dynamics is constructed.

next obstacle:
  The remaining first-principles gap is now the origin of the F54 spurion:
  derive it from a low-index dynamical alignment sector, or declare it a
  background modulus in a cutoff EFT.

verification plan:
  Search for low-index dynamics that fixes the F54 direction without adding
  large propagating 54/210 pairs.  If unavailable, revise the theorem language
  so the final model is explicitly a conditional cutoff EFT with spurionic
  source/projector backgrounds.
```

Heartbeat status 2026-05-07 06:43 Asia/Taipei:

```text
overall status:
  Still incomplete.  The low-index projector audit closed the algebraic
  question for P_X, but it did not derive the F54 background from PSLT/Another
  Physics or from a perturbative low-index Spin(10) alignment sector.

current obstacle:
  The remaining hard gap is the status of the F54 direction.  The model is
  UV-safe at high R if F54 is a fixed spurion/background and only the 45
  mediator triplet propagates, but the current manuscript still lacks a
  first-principles or low-index dynamical origin for that background.

next attempted nontrivial idea:
  Recast the F54 direction as a constrained order parameter rather than a
  propagating multiplet.  Mathematically, treat F54 as a point on the
  traceless-symmetric orbit
    S = diag(-2,-2,-2,-2,-2,-2,3,3,3,3),
  modulo Spin(10), and ask whether a low-index Lagrange-multiplier sector can
  fix only the orbit invariants Tr S^2 and Tr S^3 without introducing a full
  propagating 54 alignment tower.

verification plan:
  1. Write a finite-dimensional orbit-stabilizer check for the F54 background:
     compute stabilizer SO(6)xSO(4), number of Goldstone/orbit directions, and
     invariant constraints needed to fix eigenvalue ratio -2:3.
  2. Count the propagating Dynkin index under two interpretations:
       (A) F54 as background/constrained spurion: count only 45 mediators;
       (B) F54 as propagating 54: count one 54 plus 45 mediators.
  3. Verify both interpretations against the Landau audit and the existing
     low-index projector filter.
  4. If no low-index dynamical stabilizer exists, update TeX theorem language
     to say explicitly: conditional cutoff EFT with a fixed F54 order
     parameter, not a fully perturbative UV-complete Spin(10) GUT.
```

Heartbeat status 2026-05-07 07:15 Asia/Taipei:

```text
overall status:
  Still incomplete.  The current construction has a verified low-index
  projector filter and a verified Landau-pole diagnosis, but it has not yet
  derived the F54 order parameter from first principles.

current obstacle:
  The distinction between a fixed F54 background and a propagating 54_H field
  is now the central theoretical issue.  The former keeps the high-R branch
  UV-safe in the one-loop Spin(10) audit; the latter is only safe if the
  alignment sector is drastically smaller than the previous 54/210 tower.

next attempted nontrivial idea:
  Prove an orbit-stabilizer proposition for the F54 background.  The candidate
  order parameter
    S = diag(-2^6, 3^4)
  has stabilizer SO(6)xSO(4), and the Clebsch map on two-forms gives
    (15,1,1): -4,
    (1,3,1),(1,1,3): 6,
    (6,2,2): 1.
  Normalizing the weak-adjoint value to 2 gives exactly
    F54 = (-4/3, 2, 1/3),
  hence the low-index polynomial P_X follows without needing 210 dynamics.

verification plan:
  Implement an orbit-stabilizer audit script that:
    - computes the stabilizer dimension 15+6=21 and orbit dimension 45-21=24;
    - checks the two-form Clebsch eigenvalues and multiplicities;
    - confirms P_X=(0,0,0,1) on {Sigma_L,Sigma_R,Sigma8,X_622};
    - compares Landau counts for background F54 versus propagating F54.
  If it passes, add a TeX proposition framing F54 as a fixed constrained
  order parameter.  If it fails, collapse the claim to cutoff-EFT spurion
  input with no dynamical alignment claim.
```

Orbit-stabilizer audit completed 2026-05-07 Asia/Taipei:

```text
status:
  Completed the finite-dimensional F54 orbit-stabilizer audit and synchronized
  the TeX draft.

files:
  code/audit_f54_orbit_stabilizer.py
  output/f54_orbit_stabilizer/f54_orbit_stabilizer_summary.json
  output/f54_orbit_stabilizer/f54_orbit_stabilizer_report.md
  output/f54_orbit_stabilizer/f54_two_form_spectrum.csv
  output/f54_orbit_stabilizer/f54_order_parameter_landau.csv

mathematical result:
  For
    S0 = diag(-2,-2,-2,-2,-2,-2,3,3,3,3),
  the checks give
    Tr S0 = 0,
    Tr S0^2 = 60,
    Tr S0^3 = 60,
    S0^2 - S0 - 6 I = 0.
  The stabilizer in so(10) is so(6)+so(4), with dimension 15+6=21.
  Therefore the orbit dimension is 45-21=24.

Clebsch result:
  On adjoint two-forms, F54 is proportional to s_i+s_j.  The script finds
    color-color: multiplicity 15, raw -4, normalized -4/3;
    weak-weak:   multiplicity 6,  raw  6, normalized 2;
    mixed:       multiplicity 24, raw  1, normalized 1/3.
  Thus
    P_X = -(9/25)(F54-2)(F54+4/3)
  gives P_X=(0,0,0,1) on {Sigma_L,Sigma_R,Sigma8,X_622}, with max error
  1.11e-16.

new nontrivial element:
  The same orbit can be generated by a single renormalizable 54_H cubic
  superpotential:
    W54 = m/2 Tr S^2 + lambda/3 Tr S^3.
  The projected F-term is
    F_S = m S + lambda(S^2 - Tr S^2 I/10).
  Since S0^2 - 6 I = S0, S=vS0 is F-flat for v=-m/lambda.
  With m=lambda=1, the numerical F-term residual is exactly 0 in the script.
  D-flatness follows from [S0^dagger,S0]=0.

UV/Landau result:
  fixed F54 background + 3x45:
    sumT=24, b10=0, R=50/200 pass.
  one propagating 54 + 3x45:
    sumT=36, b10=12, R=50/200 pass; at R=200,
    alpha_G^{-1}(R M_G)=31.244.
  previous large 54/210 tower + 3x45:
    sumT=208, b10=184, R=50/200 fail; at R=200,
    alpha_G^{-1}(R M_G)=-113.796.

paper update:
  Added a "F54 Orbit-Stabilizer Audit" subsection after the low-index
  projector escape.  The paper now states that F54 is valid as a fixed
  constrained order parameter, or as one low-index 54_H order parameter, but
  not as the previously enlarged propagating 54/210 source/driver tower.

remaining caveat:
  This closes the internal representation-theory consistency check for F54,
  but it does not derive the S0 orbit uniquely from PSLT/Another Physics.
  The honest claim remains conditional EFT input unless a deeper selection
  principle for the 6+4 traceless-symmetric orbit is found.

next suggested task:
  Do a one-54_H mass-spectrum audit around W54 to see whether the physical
  non-Goldstone 54_H fragments are all at MG or whether their splittings must
  be added to the threshold/proton scan.
```

Single-54_H Hessian audit completed 2026-05-07 Asia/Taipei:

```text
status:
  Completed the full traceless-symmetric 54-dimensional Hessian audit around
  W54 and replayed the resulting threshold vector through the corrected
  two-loop/proton cache.

files:
  code/audit_single_54_hessian_spectrum.py
  output/single_54_hessian/single_54_hessian_summary.json
  output/single_54_hessian/single_54_hessian_report.md
  output/single_54_hessian/single_54_threshold_replay.csv

mathematical derivation:
  With
    W54 = m/2 Tr S^2 + lambda/3 Tr S^3,
    S = v S0 + phi,
    v = -m/lambda,
  and phi traceless symmetric, the Hessian is
    delta F = m phi
      + lambda [v(S0 phi + phi S0) - 2v Tr(S0 phi) I/10].
  The complete numerical 54x54 Hessian gives:
    (20',1,1): multiplicity 20, eigenvalue +5m;
    (1,3,3):  multiplicity 9,  eigenvalue -5m;
    (6,2,2):  multiplicity 24, eigenvalue 0;
    (1,1,1):  multiplicity 1,  eigenvalue -m.
  Matrix rank is 30, zero modes are 24, matching the Spin(10)/PS Goldstone
  orbit tangent.

nontrivial new element:
  The 24 zero modes are exactly eaten Goldstones, so the one-loop chiral
  threshold must be Goldstone-subtracted.  The active physical vector is
    b54_phys = b_(20',1,1) + b_(1,3,3)
             = (16/5,0,8) + (18/5,6,0)
             = (34/5,6,8).
  The eaten (6,2,2) would have b=(26/5,6,4), and
    (34/5,6,8)+(26/5,6,4)=(12,12,12),
  confirming that the full 54 is universal before Goldstone subtraction.

threshold convention:
  Define
    kappa54 = M[(20',1,1),(1,3,3)]/MG = 5|m|/MG.
  Then
    Delta54_i = b54_phys_i log(1/kappa54)/(2 pi),
  and
    ||P Delta54||_2 = 2.265746375e-1 |log kappa54|.
  The radial singlet has M_rad/MG = kappa54/5 and no SM one-loop threshold.

numerical verification:
  The clean benchmark kappa54=1 gives Delta54=0 and preserves the R=200 point:
    alphaG^{-1}=41.363082,
    M_Sigma3=9.211965e14 GeV,
    M_Sigma8=3.887852e15 GeV,
    tau_d6=4.809711e35 yr,
    safe points=146.
  Nonzero kappa54 replay remains viable in the cached scan but shifts the
  fitted thresholds:
    kappa54=0.5: ||P Delta54||=1.570496e-1, safe=129,
                 M_Sigma3=7.126228e14 GeV.
    kappa54=1.1: ||P Delta54||=2.159487e-2, safe=151,
                 M_Sigma3=9.542956e14 GeV.
    kappa54=2.0: ||P Delta54||=1.570496e-1, safe=169,
                 M_Sigma3=1.190817e15 GeV.
    kappa54=5.0: ||P Delta54||=3.646578e-1, safe=196,
                 M_Sigma3=1.671976e15 GeV.

paper update:
  Added a Single-54_H Hessian and Threshold Audit subsection after the
  F54 orbit-stabilizer audit.  The paper now states that the clean benchmark
  should set 5|m|=MG so all non-Goldstone non-singlet 54_H fragments sit at
  MG, while the radial singlet is gauge-neutral.

remaining caveat:
  If kappa54 is not fixed to 1 by the breaking sector, it must be included as
  an explicit heavy-threshold parameter in the global scan.  The next hard
  step is to connect the condition 5|m|=MG to the vector mass normalization
  and check whether it is natural under the PS-breaking/gauge sector.

next suggested task:
  Derive the vector multiplet mass M_V for the same S0 normalization and
  express kappa54=5|m|/M_V in terms of lambda/g.  Then decide whether
  kappa54=1 is a mild coupling relation or a new tuning.
```

Single-54_H vector-mass normalization audit completed 2026-05-07 Asia/Taipei:

```text
status:
  Completed the vector multiplet mass normalization audit for the same S0
  convention as the 54_H Hessian scan and synchronized the TeX draft.

files:
  code/audit_54_vector_mass_normalization.py
  output/single_54_vector_mass/single_54_vector_mass_summary.json
  output/single_54_vector_mass/single_54_vector_mass_report.md
  output/single_54_vector_mass/vector_generator_norms.csv
  output/single_54_vector_mass/vector_coupling_relation.csv
  output/single_54_vector_mass/kappa54_lambda_over_g_scan.csv

normalization:
  K = Tr(S^dagger S)
  Tr_10(T_a T_b)=delta_ab
  g^2 A^2 ||T<S>||^2 = (1/2) M_V^2 A^2
  This is the SUSY chiral normalization used by the superpotential/Hessian
  scripts.  In a purely real (1/2)Tr(DS DS) convention the sqrt(2) below would
  be absorbed.

mathematical result:
  For a mixed broken generator
    T_{i alpha}=(E_{i alpha}-E_{alpha i})/sqrt(2),
  with i in the 6-dimensional -2 eigenspace and alpha in the 4-dimensional
  +3 eigenspace,
    ||[T_{i alpha},S0]|| = 5.
  There are 6*4=24 such broken generators; the remaining 21 generators are the
  so(6)+so(4) stabilizer and have zero tangent norm.

vector mass:
  M_V^2 = 2 g^2 |v|^2 ||T S0||^2 = 50 g^2 |v|^2,
  hence
    M_V = 5 sqrt(2) g |v|.
  Since F-flatness gives |v|=|m|/|lambda| and the non-Goldstone 54_H mass is
    M_54_phys = 5 |m|,
  the clean-threshold ratio is
    kappa54 = M_54_phys/M_V = |lambda|/(sqrt(2) g).

coupling relation:
  kappa54=1 is equivalent to
    |lambda|/g = sqrt(2) = 1.414213562.
  At the R=200 benchmark:
    alphaG^{-1}=41.363082,
    g_G=0.551186391,
    lambda_required=0.779495270.
  This is O(1) and perturbative, not a small-number tuning.

threshold tolerance:
  From the single-54 threshold audit,
    ||P Delta54||_2 = 0.2265746375 |log kappa54|.
  Requiring ||P Delta54||<1e-2 allows
    0.956824 < kappa54 < 1.045124,
    1.353154 < lambda/g < 1.478029.
  Requiring ||P Delta54||<1e-3 tightens this to
    0.995596 < kappa54 < 1.004423,
    1.407986 < lambda/g < 1.420469.

interpretation:
  The kappa54=1 point is a mild gauge-Yukawa locking condition, not a severe
  tuning.  However, ordinary N=1 supersymmetry does not force lambda=sqrt(2)g.
  Therefore the paper has three honest options:
    1. Minimal EFT: keep kappa54 as an explicit heavy-threshold scan parameter.
    2. Boundary locking: impose lambda=sqrt(2)g at M_G as a clean benchmark.
    3. UV target: seek a fixed-line/extended-symmetry reason for lambda/g=sqrt(2).

paper update:
  Added a "Vector Multiplet Mass and lambda/g Locking" subsection.  The paper
  now explicitly states the convention, derives M_V, and frames kappa54=1 as
  a mild but not yet symmetry-protected boundary condition.

next suggested task:
  Choose one of the three interpretations above.  The most publication-honest
  route is to keep kappa54 in the benchmark card and scan, while presenting
  lambda=sqrt(2)g as the clean locked sub-branch.  If we want a stronger claim,
  the next hard problem is a gauge-Yukawa fixed-line or extended-symmetry
  derivation of lambda/g=sqrt(2).
```

Kappa54 global scan and ordinary fixed-line diagnostic completed 2026-05-07 Asia/Taipei:

```text
status:
  Completed the requested upgrade: kappa54 is now a formal benchmark-card and
  global-scan coordinate.  The locked sub-branch lambda=sqrt(2)g is retained
  explicitly, and the ordinary N=1 gauge-Yukawa fixed-line option was tested.

files:
  code/scan_kappa54_global_fixedline.py
  output/kappa54_global_scan/kappa54_global_scan.csv
  output/kappa54_global_scan/gauge_yukawa_fixedline.csv
  output/kappa54_global_scan/kappa54_benchmark_cards.json
  output/kappa54_global_scan/kappa54_global_scan_summary.json
  output/kappa54_global_scan/kappa54_global_fixedline_report.md

explicit scan axis:
  Delta54_i(kappa54)
    = (34/5, 6, 8)_i log(1/kappa54)/(2 pi),
  with
    kappa54 = |lambda|/(sqrt(2) g).
  The locked card is
    kappa54=1, lambda=sqrt(2)g, Delta54=0.

locked benchmark cards:
  R=50:
    lambda=0.814931814, safe points=148,
    M_Sigma3=9.234854e14 GeV,
    tau_d6=4.026130e35 yr.
  R=200:
    lambda=0.779495270, safe points=146,
    M_Sigma3=9.211965e14 GeV,
    tau_d6=4.809711e35 yr.

non-locked diagnostics:
  Score-only best safe grid point:
    R=50, kappa54=5.00, lambda/g=7.071068,
    ||P Delta54||=3.646578e-1, safe points=196.
    This is not a clean/natural branch; lambda is large.
  Perturbative-lambda grid point with lambda<sqrt(4pi):
    R=50, kappa54=3.00, lambda/g=4.242641,
    ||P Delta54||=2.489177e-1, safe points=181.
    This is viable as a broad EFT scan point, but not a clean threshold.
  Clean-threshold grid point with ||P Delta54||<1e-2:
    R=50, kappa54=1.01, lambda/g=1.428356,
    lambda=0.822965, ||P Delta54||=2.254493e-3,
    safe points=148.
    This confirms that clean scan points stay close to the locked branch.

one-loop cubic invariant:
  In an orthonormal 54 basis,
    d_apq d_bpq = (14/5) delta_ab.
  Therefore
    gamma_54 = [(28/5) lambda^2 - 20 g^2]/(16 pi^2),
    beta_lambda/lambda = [(84/5) lambda^2 - 60 g^2]/(16 pi^2).
  With beta_g/g = b10 g^2/(16 pi^2), a one-loop fixed ratio satisfies
    (lambda/g)_*^2 = (60+b10)/(84/5).

ordinary N=1 no-go:
  The target lambda/g=sqrt(2) requires
    b10 = -26.4.
  But ordinary N=1 Spin(10) has
    b10 = sum_R T(R) - 24 >= -24
  for nonnegative-index chiral multiplets.  Hence ordinary one-loop N=1
  gauge-Yukawa flow cannot explain the locked value.

fixed-line numerical audit:
  single 54 only:
    b10=-12, fixed lambda/g=1.690309, fixed kappa54=1.195229.
  54 plus three 45 mediators:
    b10=12, fixed lambda/g=2.070197, fixed kappa54=1.463850.
  54 plus minimal flavor Higgs:
    b10=58, fixed lambda/g=2.650247, fixed kappa54=1.874008.
  54 plus projector plus flavor Higgs:
    b10=82, fixed lambda/g=2.907298, fixed kappa54=2.055770.
  large 54/210 tower:
    b10=184, fixed lambda/g=3.811012, fixed kappa54=2.694792.

paper update:
  Added a "kappa54 Benchmark Cards and Fixed-Line Diagnostic" subsection.
  The paper now separates:
    1. minimal EFT route: scan kappa54 explicitly,
    2. locked-card route: impose lambda=sqrt(2)g as a clean boundary branch,
    3. UV route: derive sqrt(2) using extended symmetry or a nonordinary
       fixed-line mechanism.

next suggested task:
  Do not try to claim ordinary N=1 fixed-line protection.  The next hard
  theoretical task is to build and test one concrete UV locking mechanism that
  evades the b10>=-24 obstruction.  Three candidates to test:
    A. quasi-extended-supersymmetric Higgsing relation,
    B. constrained composite 54 order parameter,
    C. conformal threshold sector with extra anomalous dimensions.
  Each candidate must output a calculable lambda/g prediction and a threshold
  spectrum that can be fed back into the RGE/proton scan.
```

UV-locking mechanism audit completed 2026-05-07 Asia/Taipei:

```text
status:
  Tested three routes for evading the ordinary N=1 b10>=-24 obstruction:
    A. quasi-extended-supersymmetric/link Higgsing,
    B. constrained/composite 54 spectral orbit,
    C. conformal threshold fixed-line with extra anomalous dimensions.
  The recommended route is now B: constrained/composite 54 orbit.

files:
  code/audit_uv_locking_mechanisms.py
  output/uv_locking_mechanisms/uv_locking_mechanism_summary.json
  output/uv_locking_mechanisms/uv_locking_mechanism_report.md
  output/uv_locking_mechanisms/constrained_54_jacobian_spectrum.csv
  output/uv_locking_mechanisms/conformal_required_gamma.csv
  output/uv_locking_mechanisms/mechanism_scorecard.csv

route A: quasi-extended/link Higgsing:
  Idea:
    Put the 54 fragment inside a protected link multiplet so that an
    N=2-like gauge-Higgs relation gives M_link=M_V and hence kappa54=1.
  Group-theory check:
    10 x 10 = 1 + 45 + 54 under diagonal Spin(10),
    T(10 x 10)=20=8+12.
  Threshold check:
    If the whole 1+45+54 link is degenerate at MG, its nonuniversal one-loop
    threshold is zero.
  Status:
    Conditional.  It still needs an explicit Spin(10) link construction that
    proves the 54 fragment, not only adjoint fragments, inherits the extended
    mass locking.

route B: constrained/composite 54 spectral orbit:
  Define the order-parameter orbit by
    S=S^T, Tr S=0, S^2-S-6 I=0.
  The polynomial roots are -2 and 3.  Trace zero forces multiplicities
    (n_-,n_+)=(6,4),
  so every point is Spin(10)-conjugate to
    S0=diag(-2^6,3^4).
  The stabilizer is SO(6)xSO(4), so the orbit dimension is
    45-15-6=24.

  Linearized constraint:
    delta(S^2-S-6I)=S0 deltaS + deltaS S0 - deltaS.
  Block coefficients:
    color-color: -5,
    weak-weak: +5,
    color-weak: 0.
  Therefore all normal color-color and weak-weak symmetric fluctuations are
  removed, modulo the traceless relation:
    [6*7/2 + 4*5/2] - 1 = 30 removed directions.
  The remaining 6*4=24 directions are exactly the Spin(10)/PS orbit tangents
  and are eaten after gauging.

  Numerical verification:
    polynomial residual norm = 0,
    constraint Jacobian rank/nullity = 30/24,
    orbit tangent rank = 24,
    nonuniversal single-54 threshold = 0.
  Interpretation:
    This evades the lambda/g fixed-line obstruction because there is no
    elementary propagating non-Goldstone 54_H threshold to lock.  F54 becomes
    a constrained order parameter.  If the Lagrange multiplier/source sector
    is dynamical rather than auxiliary/composite, its complete multiplet
    spectrum must still be audited.

route C: conformal threshold fixed-line:
  Add an anomalous contribution
    beta_lambda/lambda =
      [(84/5)lambda^2 - 60g^2]/(16pi^2) + 3 gamma_CFT.
  At lambda^2=2g^2, the ordinary coefficient is -26.4, so
    3 gamma_CFT = (b10+26.4) g^2/(16pi^2).
  At R=200, g=0.551186391:
    b10=-12: gamma_CFT=9.234608e-3,
    b10=12:  gamma_CFT=2.462562e-2,
    b10=58:  gamma_CFT=5.412507e-2,
    b10=82:  gamma_CFT=6.951608e-2,
    b10=184: gamma_CFT=1.349279e-1.
  Status:
    Algebraically possible, but unbuilt.  The anomalous dimension must track
    g^2, and the charged conformal-sector thresholds must be fed back into the
    RGE/proton scan.

paper update:
  Added a "UV Locking Mechanism Audit" subsection.  It records the three
  routes and recommends the constrained/composite 54 orbit as the strongest
  current option.

next suggested task:
  Implement route B at action level:
    1. write W_orbit with explicit auxiliary/composite multipliers for
       S^2-S-6I=0 and Tr S=0,
    2. prove F-flatness and D-flatness on the constrained orbit,
    3. audit whether multiplier/source fields are auxiliary, complete
       degenerate multiplets, or new thresholds,
    4. feed any non-auxiliary source spectrum back into the RGE/proton scan.
```

Heartbeat progress check 2026-05-07 08:14 Asia/Taipei:

```text
completion status:
  The global task is not complete.  The PSLT/Another-Physics-to-GUT claim has
  been reduced to a conditional Spin(10) EFT branch, and the current hardest
  open point is now the action-level realization of the constrained/composite
  54 spectral orbit.

current obstacle:
  Route B has a strong algebraic orbit result, but it is not yet an explicit
  action-level sector.  We still must decide whether the multipliers enforcing
    S^2-S-6I=0, Tr S=0
  are auxiliary constraints, complete degenerate propagating multiplets, or
  a composite/source sector with its own threshold vector.

next attempted nontrivial idea:
  Build W_orbit as a Spin(10)-covariant constrained-order-parameter sector:
    W_orbit = <Xi, S^2-S-6I> + xi Tr S
  or its traceless-projected equivalent.  Then compute all F-terms, D-terms,
  and the local Hessian around S0=diag(-2^6,3^4).  The creative point is to
  treat F54 as a spectral projector/orbit variable rather than an elementary
  radial 54_H field, which removes the kappa54 threshold instead of tuning it.

verification plan:
  1. Symbolically derive F_Xi=0, F_xi=0, F_S=0 conditions and identify which
     fields must be auxiliary or source-like.
  2. Numerically construct the Hessian in the 54+multiplier component basis.
  3. Verify that the physical non-Goldstone 54 threshold remains absent, or
     else record the exact source-spectrum beta vector.
  4. Feed any nonzero threshold vector back into the corrected RGE/proton scan.
  5. Synchronize the TeX section and benchmark cards after the Hessian audit.
```

Heartbeat progress check 2026-05-07 08:45 Asia/Taipei:

```text
completion status:
  Still incomplete.  No new proof closes the full conditional Spin(10) model
  yet; the next required deliverable remains the action-level constrained
  orbit sector.

current obstacle:
  The algebraic constrained orbit proves that the normal 54 modes can be
  removed, but the paper has not yet shown an explicit W_orbit Hessian with
  multiplier/source fields classified as auxiliary, complete degenerate
  multiplets, or threshold-generating states.

next attempted nontrivial idea:
  Treat the orbit constraint as a holomorphic matrix equation but project the
  multiplier Xi onto the traceless-symmetric 54 dual, so F_S does not force an
  unwanted scalar identity component.  Compare this projected W_orbit with the
  unprojected <Xi,S^2-S-6I> form to see which version has the cleaner
  auxiliary spectrum.

verification plan:
  Construct the component Hessian for both projected and unprojected
  multiplier choices, count zero modes and massive pairs, then accept only a
  version whose non-Goldstone threshold vector is either zero or a complete
  Spin(10) multiplet contribution that can be fed into the existing
  RGE/proton scripts.
```

Orbit-superpotential Hessian audit completed 2026-05-07 09:16 Asia/Taipei:

```text
status:
  Made concrete progress on route B.  A component Hessian audit now compares
  unprojected, traceless-projected, and normal-bundle multiplier choices for
    W^(2)=<delta Xi, J(delta S)>,
    J(delta S)=S0 deltaS + deltaS S0 - deltaS.

files:
  code/audit_orbit_superpotential_hessian.py
  output/orbit_superpotential_hessian/orbit_superpotential_hessian_summary.json
  output/orbit_superpotential_hessian/orbit_superpotential_hessian_report.md
  output/orbit_superpotential_hessian/orbit_superpotential_hessian_cases.csv

results:
  unprojected symmetric Xi55:
    rank(J)=30, S-nullity=24, Xi-nullity=25,
    Hessian zero modes=49.
    Reject as a dynamical multiplier because it leaves 25 extra massless
    multiplier modes.  It is acceptable only as a purely auxiliary constraint.

  traceless projected Xi54:
    rank(J)=30, S-nullity=24, Xi-nullity=24,
    Hessian zero modes=48.
    Reject as a dynamical multiplier because it leaves 24 extra massless
    multiplier modes.  It is acceptable only if auxiliary or further
    gauge-fixed/quotiented.

  normal-bundle Xi30:
    rank(J)=30, S-nullity=24, Xi-nullity=0,
    Hessian dimension=84, Hessian rank=60, zero modes=24.
    All nonzero singular values are 5 to numerical precision.
    This is the best action-level candidate: it pairs all 30 normal S modes
    and leaves only the 24 Goldstone orbit tangents.

threshold consequence:
  If Xi30 is auxiliary/composite:
    Delta b_orbit = (0,0,0).
  If Xi30 is promoted to a propagating chiral sector:
    Delta b >= 2 b54_phys = (68/5, 12, 16),
  unless the sector is exactly at the matching scale or completed into a
  universal multiplet.  Therefore the clean route requires an auxiliary or
  composite normal-bundle constraint, not a generic dynamical Xi54/Xi55 field.

paper update:
  Added an "Orbit-Superpotential Hessian Audit" subsection to the TeX draft.

next obstacle:
  Formalize the normal-bundle multiplier globally enough to be Spin(10)
  covariant.  Locally it is the image of J at S0; globally it should be the
  conormal bundle to the spectral orbit, or an equivalent auxiliary/composite
  source.  Need prove that using this conormal-bundle constraint is compatible
  with supersymmetric F-flatness and gauge D-flatness.

next verification plan:
  1. Write the conormal-bundle version of W_orbit in invariant language.
  2. Show F_Xi imposes C(S)=0 only in normal directions while F_S is solved
     at Xi=0.
  3. Verify D^a=S^\dagger T^a S vanishes on the real symmetric orbit.
  4. If any conormal source is made propagating, feed its threshold vector
     into the existing RGE/proton scan.
```

Conormal F/D-flatness audit completed 2026-05-07 09:48 Asia/Taipei:

```text
status:
  Advanced route B from local Hessian counting to an explicit conormal-bundle
  F-flatness and D-flatness check over sampled Spin(10) orbit points.

files:
  code/verify_conormal_orbit_flatness.py
  output/conormal_orbit_flatness/conormal_orbit_flatness_summary.json
  output/conormal_orbit_flatness/conormal_orbit_flatness_report.md
  output/conormal_orbit_flatness/conormal_orbit_flatness_samples.csv

mathematical setup:
  C(S)=S^2-S-6I,
  J_S(deltaS)=S deltaS + deltaS S - deltaS.
  The conormal superpotential is locally
    W_conormal=<Xi_N,C(S)>,
    Xi_N in N^*_S O_54 ~= im J_S.
  F-terms:
    F_Xi=P_N C(S),
    F_S=J_S^dagger Xi_N.
  Hence on C(S)=0 with Xi_N=0:
    F_Xi=0, F_S=0.

D-flatness proof:
  The SO(10) action on S is S -> O S O^T, so T.S = TS-ST.
  For real symmetric S and antisymmetric T,
    D_T=<S,T.S>=Tr[S(TS-ST)]
       =Tr(TS^2)-Tr(S^2T)=0.
  Thus every real orbit point S=O S0 O^T is D-flat.

numerical verification:
  Samples: S0 plus five random SO(10)-conjugate points.
  All samples satisfy:
    ||C(S)|| <= 1.672e-14,
    ||D|| <= 1.682e-15,
    rank J_S = 30,
    dim ker J_S = 24,
    tangent rank = 24.

paper update:
  Added a "Conormal F- and D-Flatness Check" subsection to the TeX draft.

current obstacle:
  Flatness is now verified.  The remaining hard issue is UV/model-building:
  Xi_N must be auxiliary/composite, or if propagating it must form a complete
  degenerate threshold sector and be fed into the RGE/proton scan.

next attempted nontrivial idea:
  Construct a minimal auxiliary-source implementation of Xi_N as a conormal
  Lagrange multiplier rather than a chiral propagating 54/55.  If this cannot
  be made precise, define a fallback propagating normal-pair threshold vector
    Delta b >= (68/5,12,16)
  and scan whether the corrected branch survives.

verification plan:
  1. Write the auxiliary-source charge/representation assignment.
  2. Check whether Xi_N has a kinetic term in the effective action.
  3. If no kinetic term, record Delta b=0 and close the kappa54-locking issue.
  4. If kinetic, feed (68/5,12,16) with a mass scale into the threshold/proton
     scripts and determine whether the branch remains viable.
```

Propagating Xi_N threshold fallback scan completed 2026-05-07 10:21 Asia/Taipei:

```text
status:
  Tested the fallback case where Xi_N is not auxiliary/composite but propagates
  as the normal-pair chiral sector found by the Hessian audit.

files:
  code/scan_xin_propagating_threshold.py
  output/xin_propagating_threshold/xin_propagating_threshold_summary.json
  output/xin_propagating_threshold/xin_propagating_threshold_report.md
  output/xin_propagating_threshold/xin_propagating_threshold_scan.csv

threshold:
  b_XiN = 2 b54_phys = (68/5, 12, 16).
  For M_XiN=kappa_XiN MG:
    Delta_i^XiN = b_i^XiN log(1/kappa_XiN)/(2 pi).
  Projected size:
    ||P Delta_XiN||_2 = 0.453149275 |log kappa_XiN|.

clean threshold windows:
  ||P Delta_XiN||<1e-2:
    0.978174 < kappa_XiN < 1.022313.
  ||P Delta_XiN||<1e-3:
    0.997796 < kappa_XiN < 1.002209.

selected scan results:
  R=50:
    kappa=0.50: ||P Delta||=3.140991e-1, safe=107,
                M_Sigma3=5.526432e14 GeV, tau_d6=3.715899e35 yr.
    kappa=1.00: ||P Delta||=0, safe=148,
                M_Sigma3=9.234854e14 GeV, tau_d6=4.026130e35 yr.
    kappa=1.01: ||P Delta||=4.508985e-3, safe=149,
                M_Sigma3=9.303172e14 GeV, tau_d6=4.030674e35 yr.
    kappa=2.00: ||P Delta||=3.140991e-1, safe=192,
                M_Sigma3=1.543175e15 GeV, tau_d6=4.348797e35 yr.
  R=200:
    kappa=0.50: ||P Delta||=3.140991e-1, safe=107,
                M_Sigma3=5.512735e14 GeV, tau_d6=4.470054e35 yr.
    kappa=1.00: ||P Delta||=0, safe=146,
                M_Sigma3=9.211965e14 GeV, tau_d6=4.809711e35 yr.
    kappa=1.01: ||P Delta||=4.508985e-3, safe=148,
                M_Sigma3=9.280114e14 GeV, tau_d6=4.814677e35 yr.
    kappa=2.00: ||P Delta||=3.140991e-1, safe=190,
                M_Sigma3=1.539350e15 GeV, tau_d6=5.161804e35 yr.

interpretation:
  A propagating Xi_N normal pair does not immediately kill the branch; the
  corrected scan can absorb it.  But it becomes a new nonuniversal threshold
  parameter.  For a clean paper-level branch, Xi_N must be auxiliary/composite,
  locked within about 2.2% of MG, or completed into a universal degenerate
  threshold sector.

paper update:
  Added a "Propagating Xi_N Fallback Threshold" subsection.

current obstacle:
  Decide whether to commit to auxiliary/composite Xi_N as an assumption, or
  enlarge the benchmark card with kappa_XiN as a new threshold coordinate.

next verification plan:
  If auxiliary/composite is adopted, write the precise EFT assumption and close
  the kappa54-locking sector.  If not, merge kappa_XiN into the global
  benchmark cards and proton scan.
```

Heartbeat closure update 2026-05-07 10:53 Asia/Taipei:

```text
completion status:
  Full GUT task remains incomplete, but the kappa54-locking subproblem is now
  conditionally closed.

conditional closure statement:
  Under assumptions
    A1. S is a constrained spectral-orbit order parameter:
        S^2-S-6I=0, Tr S=0;
    A2. Xi_N is auxiliary/composite, not a propagating chiral multiplet;
    A3. the 24 orbit tangents are eaten Goldstone multiplets,
  the kappa54 sector contributes no nonuniversal threshold:
    Delta b_54_orbit = 0,
    Delta b_XiN_aux = 0.

locked cards retained:
  R=50:
    safe=148, M_Sigma3=9.234854e14 GeV, tau_d6=4.026130e35 yr.
  R=200:
    safe=146, M_Sigma3=9.211965e14 GeV, tau_d6=4.809711e35 yr.

fallback if A2 is rejected:
  Add kappa_XiN=M_XiN/MG as a new benchmark coordinate with
    b_XiN=(68/5,12,16).
  Clean threshold requires
    0.978174 < kappa_XiN < 1.022313
  for ||P Delta_XiN||<1e-2.

paper update:
  Added "Conditional Closure of the kappa54-Locking Sector" to the TeX draft.

current obstacle:
  The remaining hard blockers have shifted away from kappa54 locking.  The
  most important outstanding paper-level issues are now:
    1. make the CP1/O(2) Yukawa geometry action-level or explicitly downgrade
       it to a texture ansatz,
    2. provide a full reproducible flavor/seesaw benchmark,
    3. complete channel-specific dimension-five proton decay,
    4. decide whether kappa_XiN is an assumption (auxiliary/composite) or a
       benchmark-card threshold axis.

next attempted nontrivial idea:
  Move to the flavor/Yukawa geometry blocker: replace the non-covariant
  int psi_i psi_j h dmu expression by a bundle-covariant O(-4)-valued kernel
  or a Hermitian Toeplitz operator, then rerun singular-value hierarchy scans.

verification plan:
  1. Define the coordinate-invariant Yukawa functional.
  2. Regenerate Yukawa matrices and singular-value ratios.
  3. Check whether the previous seesaw/proton conclusions survive.
  4. Synchronize TeX, roadmap, and benchmark cards.
```

Heartbeat update 2026-05-07 11:25 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the CP1/O(2) Yukawa covariance
  blocker has a concrete action-level resolution.

result:
  The legacy holomorphic texture
    Y_ij = int psi_i psi_j h_a dmu_FS
  is not coordinate-invariant if h_a is an ordinary scalar, because
    psi_i psi_j in H^0(CP1,O(4)).
  It is promoted to a covariant superpotential functional by declaring the
  Higgs functional to be an O(-4)-valued dual density:
    Y_ij = <H_a, psi_i psi_j>,
    H_a in Gamma_sm(O(-4)) tensor Omega^{1,1}.

numerical audit:
  O(2) Gram max error = 1.998e-15.
  O(4) Gram max error = 3.775e-15.
  Product-map reconstruction max error = 8.043e-15.
  Veronese identity psi_1^2=2 psi_0 psi_2 error = 9.930e-16.
  Legacy direct-dual reconstruction errors:
    up=1.319e-16, down=1.667e-16,
    charged_lepton=7.633e-17, neutrino_dirac=5.222e-18.

Toeplitz alternative:
  The Hermitian Toeplitz functional
    T_ij = int psi_i^* h_a psi_j dmu_FS
  is coordinate-invariant for scalar h_a and can generate hierarchies, but it
  naturally belongs to the Kahler/wavefunction sector unless an additional
  holomorphic projection is specified.  Minimal two-center Toeplitz scan:
    up: found (4.636e-4, 1.360e-2), target (1e-5, 7e-3);
    down: found (1.109e-3, 1.995e-2), target (1e-3, 2e-2);
    charged_lepton: found (4.828e-4, 2.956e-2), target (3e-4, 6e-2);
    neutrino_dirac: found (1.015e-2, 1.974e-1), target (1e-2, 2e-1).

paper update:
  Added "Covariant O(2) Yukawa functional audit" to paper/gut_framework.tex.

current obstacle:
  The action-level Yukawa map is now defined, but the benchmark cards still
  need to expose the O(-4) dual coefficients, full Yukawa matrices, basis
  convention, and the downstream seesaw reconstruction.

next attempted nontrivial idea:
  Convert the O(-4) dual-density coefficients into machine-readable benchmark
  cards and rerun the seesaw/flavor reproducibility audit using the exact
  matrices rather than singular-value ratios alone.

verification plan:
  1. Export exact O(-4) dual coefficients and Y_u,Y_d,Y_e,Y_nu matrices.
  2. Regenerate m_D, m_nu, M_R with explicit PMNS/Takagi conventions.
  3. Check whether the existing proton and RGE benchmarks are unchanged.
  4. Synchronize TeX appendix tables with the benchmark JSON/CSV files.
```

Heartbeat update 2026-05-07 12:03 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the Yukawa/seesaw reproducibility
  blocker is substantially reduced.

new script:
  code/export_flavor_benchmark_card.py

new outputs:
  output/flavor_benchmark/flavor_benchmark_card.json
  output/flavor_benchmark/dual_density_coefficients.csv
  output/flavor_benchmark/flavor_benchmark_report.md

result:
  The O(-4) dual-density branch now exports:
    product Clebsch tensor C_ij^m,
    dual coefficients a_m^(a),
    normalized Y_u,Y_d,Y_e,Y_nu matrices,
    U_e, U_PMNS convention,
    m_D, m_nu, M_R,
    trace-lift Majorana fit.

numerical audit:
  O(2) Gram error = 1.998e-15.
  O(4) Gram error = 3.775e-15.
  product-map reconstruction = 8.043e-15.

  Direct-dual normalized reconstruction errors:
    up             = 1.888e-15
    down           = 1.111e-15
    charged_lepton = 1.554e-15
    neutrino_dirac = 1.225e-15

  Exact-card seesaw:
    M1,M2,M3 = 2.599776e10, 9.144789e13, 3.105386e15 GeV.
    theta_norm = 4.393070e-11.
    seesaw residual = 2.778764e-12.
    light-mass residual = 2.012888e-12.
    PMNS sin^2 = (0.304000, 0.022200, 0.573000).
    trace-lift residual = 2.907738e-16.

paper update:
  Added "Flavor benchmark-card reproducibility" to paper/gut_framework.tex.

current obstacle:
  The exact flavor card is now reproducible, but it is still not a full
  phenomenological flavor fit.  It does not yet fit quark masses, charged
  lepton masses, CKM, PMNS, and proton-decay flavor rotations simultaneously.

next attempted nontrivial idea:
  Promote the exact O(-4) Yukawa card from hierarchy benchmark to a small
  complex-coefficient flavor fit: allow controlled sector-dependent
  dual-density perturbations and fit CKM/PMNS while preserving the Veronese
  structure and seesaw trace-lift.

verification plan:
  1. Build a flavor objective using Y_u,Y_d,Y_e,Y_nu from the exported card.
  2. Add unitary rotations and extract CKM/PMNS observables.
  3. Fit or bound quark/lepton mass ratios and mixing angles.
  4. Feed the resulting flavor rotations into channel-specific d=5 proton
     decay, especially p -> K+ nu_bar.
```

Heartbeat update 2026-05-07 12:34 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The exact O(-4) flavor card now has
  a first full-observable audit, and it fails CKM as expected.

new script:
  code/audit_flavor_fit_observables.py

new outputs:
  output/flavor_fit/flavor_observable_audit.json
  output/flavor_fit/flavor_observable_audit_report.md
  output/flavor_fit/ckm_current_vs_target.csv

derivation:
  From exact matrices,
    H_u = Y_u Y_u^dagger = U_u D_u^2 U_u^dagger,
    H_d = Y_d Y_d^dagger = U_d D_d^2 U_d^dagger,
    V_CKM = U_u^dagger U_d.
  Keeping D_d fixed, a CKM-only down-sector target is
    H_d^target = V_CKM^target D_d^2 (V_CKM^target)^dagger
  in the up-mass basis.

numerical audit:
  Mass hierarchy log scores remain small:
    up=8.325e-3,
    down=1.558e-3,
    charged_lepton=9.588e-4,
    neutrino_dirac=1.537e-3.

  Current exact-card CKM:
    |V_us| = 6.702511e-1,
    |V_cb| = 5.175675e-1,
    |V_ub| = 3.726974e-1,
    J      = 3.140935e-2.

  Internal target used only for local audit:
    |V_us| = 2.249985e-1,
    |V_cb| = 4.099972e-2,
    |V_ub| = 3.700000e-3,
    J      = 3.097062e-5.

  CKM magnitude log-score = 5.449970.

  Minimal down-Hermitian deformation with down singular values fixed:
    ||Delta H_d||_F / ||H_d||_F = 8.964034e-1.
    ||Delta H_d||_2 / ||H_d||_2 = 6.339201e-1.
    ||offdiag Delta H_d||_F / y_b^2 = 7.412834e-1.
    max |offdiag Delta H_d| / y_b^2 = 3.950184e-1.

paper update:
  Added "Flavor observable audit" to paper/gut_framework.tex.

current obstacle:
  The O(-4) dual-density card is now a reproducible hierarchy/seesaw
  benchmark, but not a full flavor fit.  CKM is far too large, especially
  |V_ub|.  Proton d=5 cannot be made channel-specific until these flavor
  rotations are fitted.

next attempted nontrivial idea:
  Add a controlled left-misalignment sector.  The cleanest three candidates
  are:
    1. constrained 120_H antisymmetric contribution,
    2. Toeplitz/Kahler canonical-normalization metric,
    3. small O(-4) dual-density perturbations around the exported card.
  Fit CKM first while preserving the successful mass hierarchy and seesaw
  trace-lift, then feed the rotations into p -> K+ nu_bar.

verification plan:
  1. Build a small objective for mass ratios + CKM magnitudes + Jarlskog.
  2. Scan the three misalignment mechanisms separately.
  3. Reject mechanisms that require O(1) destruction of the hierarchy card.
  4. Export the fitted rotations and replay d=5 proton decay channel by channel.
```

Heartbeat update 2026-05-07 13:12 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The first candidate CKM repair
  mechanism, a down-sector effective 120_H antisymmetric perturbation, has
  been tested and rejected in its one-sector form.

new script:
  code/scan_120h_ckm_misalignment.py

new outputs:
  output/flavor_120h_scan/scan_120h_ckm_misalignment.json
  output/flavor_120h_scan/scan_120h_ckm_misalignment_report.md
  output/flavor_120h_scan/scan_120h_ckm_bounds.csv

derivation:
  Add a complex antisymmetric family matrix A to the exact down Yukawa card:
    Y_d -> Y_d + A,  A^T = -A.
  Keep Y_u fixed and compute
    V_CKM = U_u^dagger U_d
  from the left eigensystems of Y_u Y_u^dagger and
    (Y_d+A)(Y_d+A)^dagger.
  Optimize the six real parameters of A against CKM magnitudes, Jarlskog, and
  down mass ratios.

numerical audit:
  Component bounds scanned:
    0.02, 0.05, 0.10, 0.20, 0.50, 1.00.

  Best point:
    |V_us| = 2.936512e-1,
    |V_cb| = 4.390998e-1,
    |V_ub| = 2.746715e-1,
    |J|    = 3.097151e-5.

  Down hierarchy is preserved:
    y_d/y_b = 1.009625e-3,
    y_s/y_b = 2.295305e-2.

  Shift size:
    ||A||_F / ||Y_d||_F = 2.563590e-1,
    ||A||_2 / ||Y_d||_2 = 1.813057e-1.

  CKM magnitude log-score only improves from 5.449970 to 4.573016, so the
  one-sector 120_H correction fails mainly because |V_cb| and |V_ub| remain
  much too large.

paper update:
  Added "Effective 120_H down-sector audit" to paper/gut_framework.tex.

current obstacle:
  A one-sector antisymmetric perturbation can preserve hierarchy and tune J,
  but it cannot make CKM small.  The obstruction appears structural rather
  than a lack of scan range.

next attempted nontrivial idea:
  Do not enlarge the same one-sector scan.  Move to either:
    1. Spin(10) Clebsch-correlated 120_H contributions across u,d,e,nu, or
    2. Toeplitz/Kahler canonical-normalization metric acting on left-handed
       families.
  The second route may be more efficient because it directly rotates the
  shared left-handed Q and L multiplets without destroying holomorphic
  singular hierarchies.

verification plan:
  1. Build a Toeplitz/Kahler metric K_Q on the CP1 O(2) family space.
  2. Canonically normalize Y_u,Y_d via K_Q^{-1/2}.
  3. Scan whether CKM can be made small while preserving mass ratios and the
     exact seesaw card.
  4. If successful, export the rotations for d=5 proton decay; if not, return
     to a full Clebsch-correlated 120_H sector.
```

Heartbeat update 2026-05-07 14:04 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The second candidate CKM repair
  mechanism, an effective Toeplitz/Kahler left-handed Q metric, has been
  tested.  It improves the obstruction but is not sufficient as a Q-only
  solution.

new script:
  code/scan_toeplitz_kahler_ckm.py

new outputs:
  output/flavor_kahler_scan/scan_toeplitz_kahler_ckm.json
  output/flavor_kahler_scan/scan_toeplitz_kahler_ckm_report.md
  output/flavor_kahler_scan/scan_toeplitz_kahler_ckm.csv

derivation:
  In finite H^0(CP1,O(2)), a Hermitian endomorphism can be represented as a
  Berezin-Toeplitz operator with sufficiently general real symbol.  The scan
  uses
    K_Q = exp(G),  G=G^dagger, Tr G=0,
    Y_u,d -> K_Q^{-1/2} Y_u,d.
  CKM is recomputed from the left eigensystems of the canonical matrices.

numerical audit:
  Raw exact-card CKM score was 5.449970.
  Down-only 120_H best CKM score was 4.573016.
  Q-only Toeplitz/Kahler best mass-preserving CKM score is 3.370088.

  Best point:
    |V_us| = 5.778595e-1,
    |V_cb| = 3.501359e-2,
    |V_ub| = 2.271854e-1,
    |J|    = 3.054788e-5.

  Mass ratios:
    y_u/y_t = 7.185015e-6,
    y_c/y_t = 6.244757e-3,
    y_d/y_b = 6.850816e-4,
    y_s/y_b = 1.555073e-2.

  Metric size:
    cond(K_Q) = 3.145789e1,
    span(G)  = 3.448650.

verdict:
  The K_Q route is more promising than down-only 120_H because it suppresses
  |V_cb| to near the target scale, but it leaves |V_us| and |V_ub| much too
  large and requires a sizable metric condition number.  Therefore Q-only
  Kahler normalization is not a complete solution.

paper update:
  Added "Toeplitz/Kahler metric audit" to paper/gut_framework.tex.

current obstacle:
  CKM remains the main blocker before channel-specific d=5 proton decay can be
  meaningfully replayed.  Two simple one-sector mechanisms have now been
  rejected or found insufficient.

next attempted nontrivial idea:
  Move to a correlated flavor model:
    A. K_16 or Pati-Salam correlated Kahler metrics acting simultaneously on
       Q,L,u^c,d^c,e^c,nu^c with seesaw replay; or
    B. full Clebsch-correlated 10_H + overline{126}_H + 120_H Yukawa fit.
  Route A is numerically cheaper and directly tests whether wavefunction
  geometry can handle CKM without destroying PMNS/seesaw.

verification plan:
  1. Implement correlated K_16/K_PS canonical normalization.
  2. Recompute Y_u,Y_d,Y_e,Y_nu, CKM, PMNS, and M_R.
  3. Accept only if mass ratios, CKM, PMNS, and seesaw residual survive.
  4. Export fitted rotations for p -> K+ nu_bar d=5 proton decay.
```

Heartbeat update 2026-05-07 18:10 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The correlated K_16 / Pati-Salam
  Kahler flavor branch has now been tested.  It preserves the inverse-seesaw
  replay but still fails as a phenomenological CKM fit.

new script:
  code/scan_correlated_kahler_flavor.py

new outputs:
  output/flavor_correlated_kahler_scan/scan_correlated_kahler_flavor.json
  output/flavor_correlated_kahler_scan/scan_correlated_kahler_flavor_report.md
  output/flavor_correlated_kahler_scan/scan_correlated_kahler_flavor.csv

derivation:
  K_16 branch:
    K_16 = exp(G),  G=G^dagger, Tr G=0,
    A = exp(-G/2),
    Y_a -> A^T Y_a A for a=u,d,e,nu.

  Pati-Salam branch:
    K_L = exp(G_L), K_R = exp(G_R),
    A_L = exp(-G_L/2), A_R = exp(-G_R/2),
    Y_a -> A_L^T Y_a A_R.

  CKM is recomputed from
    Y_u Y_u^dagger = U_u D_u^2 U_u^dagger,
    Y_d Y_d^dagger = U_d D_d^2 U_d^dagger,
    V_CKM = U_u^dagger U_d.

  The seesaw replay uses the transformed Y_nu and Y_e with
    M_R = - m_D^T m_nu^{-1} m_D.

numerical audit:
  Best K_16 point:
    bound = 0.7,
    |V_us| = 7.8616e-1,
    |V_cb| = 4.6062e-1,
    |V_ub| = 3.2047e-1,
    CKM score = 5.153e0.

  Best Pati-Salam point:
    bound = 1.2,
    |V_us| = 5.600883e-1,
    |V_cb| = 2.336227e-1,
    |V_ub| = 9.250403e-2,
    |J|    = 3.166432e-5,
    CKM score = 2.682300e0,
    mass score = 8.131820e-2,
    max cond(K) = 4.3837e1.

  Seesaw replay at the best Pati-Salam point:
    theta_norm = 5.220052e-11,
    seesaw residual = 7.595655e-12,
    (M1,M2,M3) =
      (1.884656e10, 4.602145e13, 8.512023e15) GeV.

  Internal reconstruction checks pass, but phenomenology viability checks fail.
  In particular |V_ub| is still about twenty-five times too large.

paper update:
  Added "Correlated Kahler flavor audit" to paper/gut_framework.tex.

current obstacle:
  CKM/flavor is now the leading blocker.  One-sector 120_H, Q-only Kahler, and
  correlated K_16 / Pati-Salam Kahler repairs have all been tested and found
  insufficient as complete flavor fits, even though the exact seesaw card is
  stable under the correlated Kahler scan.

next attempted nontrivial idea:
  Stop treating CKM as a wavefunction-only repair.  Build the full
  Clebsch-correlated 10_H + overline{126}_H + 120_H Yukawa fit:
    Y_u = r_u^10 H + r_u^126 F + r_u^120 G,
    Y_d = r_d^10 H + r_d^126 F + r_d^120 G,
    Y_e = r_e^10 H + r_e^126 F + r_e^120 G,
    Y_nu = r_nu^10 H + r_nu^126 F + r_nu^120 G,
  with H,F symmetric and G antisymmetric, constrained by the CP1/O(-4)
  geometric basis and by the existing seesaw benchmark.

verification plan:
  1. Implement a reproducible Clebsch-correlated flavor objective.
  2. Fit quark ratios, charged-lepton ratios, CKM magnitudes, Jarlskog,
     neutrino angles, and inverse-seesaw stability simultaneously.
  3. Export U_u,U_d,U_e,U_nu and all Wilson flavor rotations.
  4. Only then replay channel-specific p -> K+ nu_bar dimension-five proton
     decay.
```

Heartbeat update 2026-05-07 19:10 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The first CP1/O(-4)-restricted
  Clebsch-correlated 10_H + overline{126}_H + 120_H flavor baseline has been
  implemented and scanned.  It is the first flavor branch that can make CKM
  small, but the best CKM point still fails the mass-hierarchy objective.

new script:
  code/scan_clebsch_flavor_fit.py

new outputs:
  output/flavor_clebsch_scan/scan_clebsch_flavor_fit.json
  output/flavor_clebsch_scan/scan_clebsch_flavor_fit_report.md
  output/flavor_clebsch_scan/scan_clebsch_flavor_fit.csv

derivation:
  The renormalizable Spin(10) Yukawa superpotential is modeled as
    W_Y = 16_i (H_ij 10_H + F_ij overline{126}_H + G_ij 120_H) 16_j.

  The CP1/O(-4) geometric restriction is:
    H,F in the five-complex-dimensional symmetric Veronese image,
    G in Lambda^2 C^3.

  The baseline doublet-mixing ansatz is
    Y_u  = H + r_u F + g_u G,
    Y_nu = H - 3 r_u F + g_nu G,
    Y_d  = r_d (H + F + g_d G),
    Y_e  = r_d (H - 3 F + g_e G),
  with complex coefficients.

numerical audit:
  Balanced bound=1.2:
    |V_us| = 2.2147e-1,
    |V_cb| = 1.1981e-1,
    |V_ub| = 3.0609e-2,
    CKM score = 1.059e0,
    mass score = 6.171e-1.

  Best CKM-heavy point, bound=2.4:
    |V_us| = 2.319941e-1,
    |V_cb| = 4.983639e-2,
    |V_ub| = 3.824997e-3,
    |J|    = 4.213567e-5,
    CKM score = 7.570255e-3,
    mass score = 1.306048e0.

  Mass ratios at the best CKM point:
    up:    (1.020173e-5, 4.710678e-2),
    down:  (8.490660e-4, 2.610979e-2),
    lepton:(1.507355e-3, 4.249504e-2),
    nu_D:  (1.961822e-2, 1.789473e-1).

  Seesaw replay at best CKM point:
    theta_norm = 1.113119e-11,
    seesaw residual = 6.071118e-13,
    (M1,M2,M3) =
      (2.624896e11, 1.594168e13, 4.742979e15) GeV.

  Internal reconstruction checks pass, but phenomenology viability checks fail
  because the mass score remains too large.  This is a Pareto tension rather
  than a no-go: the Clebsch branch has enough angular freedom for CKM, but the
  current finite baseline does not yet keep all four Dirac hierarchies.

paper update:
  Added "Clebsch-correlated flavor baseline" to paper/gut_framework.tex.

current obstacle:
  Full flavor fit remains the leading blocker.  The model cannot yet export
  trustworthy Wilson flavor rotations for p -> K+ nu_bar because the CKM-fit
  point does not simultaneously fit the quark/lepton hierarchy targets.

next attempted nontrivial idea:
  Run a Pareto/homotopy Clebsch scan:
    1. start from the CKM-heavy point,
    2. gradually raise the mass-hierarchy weights,
    3. keep CKM score below 0.05,
    4. monitor seesaw rank, theta_norm, heavy M_R, and coefficient naturalness.

verification plan:
  1. Export the best CKM-heavy x-vector as a seed card.
  2. Implement homotopy stages in mass-weight lambda.
  3. Produce a Pareto table CKM-score vs mass-score.
  4. If a viable point is found, export U_u,U_d,U_e,U_nu and begin
     channel-specific d=5 proton decay; if not, enlarge the ansatz only in a
     representation-controlled way.
```

Heartbeat update 2026-05-07 19:37 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The Pareto/homotopy Clebsch scan has
  been completed.  It confirms that the single-120_H Clebsch branch can keep
  CKM small while seesaw remains stable, but the mass score plateaus around
  0.8 and does not reach a joint full-flavor fit.

new script:
  code/scan_clebsch_flavor_homotopy.py

new outputs:
  output/flavor_clebsch_homotopy/scan_clebsch_flavor_homotopy.json
  output/flavor_clebsch_homotopy/scan_clebsch_flavor_homotopy_report.md
  output/flavor_clebsch_homotopy/scan_clebsch_flavor_homotopy.csv

homotopy derivation:
  Start from the CKM-heavy CP1/O(-4)-restricted Clebsch seed and optimize the
  same ansatz
    Y_u  = H + r_u F + g_u G,
    Y_nu = H - 3 r_u F + g_nu G,
    Y_d  = r_d (H + F + g_d G),
    Y_e  = r_d (H - 3 F + g_e G),
  while increasing the mass-hierarchy weight lambda_m.

numerical Pareto table:
  ckm_lock:
    |V_us|=2.1980e-1, |V_cb|=4.8160e-2, |V_ub|=3.7215e-3,
    CKM score=4.996e-3, mass score=9.520e-1.
  mass_1:
    |V_us|=2.1870e-1, |V_cb|=4.7361e-2, |V_ub|=3.7198e-3,
    CKM score=4.081e-3, mass score=8.961e-1.
  mass_2:
    |V_us|=2.1394e-1, |V_cb|=4.7188e-2, |V_ub|=3.7628e-3,
    CKM score=4.259e-3, mass score=8.747e-1.
  mass_4:
    |V_us|=2.0372e-1, |V_cb|=4.7147e-2, |V_ub|=3.8397e-3,
    CKM score=5.801e-3, mass score=8.360e-1.
  mass_8:
    |V_us|=1.982582e-1, |V_cb|=4.731688e-2, |V_ub|=3.848995e-3,
    |J|=3.297653e-5,
    CKM score=7.186513e-3, mass score=7.952350e-1.

best scalarized Pareto point:
  Stage mass_8.
  Mass ratios:
    up:     (1.011333e-5, 4.288794e-2),
    down:   (9.736393e-4, 2.499123e-2),
    lepton: (3.121904e-4, 3.824463e-2),
    nu_D:   (2.238098e-2, 1.699339e-1).

  Seesaw replay:
    theta_norm = 1.108766e-11,
    seesaw residual = 6.388765e-13,
    (M1,M2,M3) =
      (2.907569e11, 1.682372e13, 4.710668e15) GeV.

checks:
  any CKM score < 0.05: yes.
  any mass score < 0.20: no.
  any joint viable point: no.
  best seesaw residual < 1e-10: yes.
  best theta_norm < 1e-8: yes.

paper update:
  Added "Clebsch flavor homotopy obstruction" to paper/gut_framework.tex.

current obstacle:
  The flavor blocker is now localized.  The obstruction is not CKM angular
  freedom and not seesaw rank; it is the restricted single-120_H Clebsch /
  CP1-O(-4) symmetric tensor space.  Proton d=5 Wilson rotations remain
  premature until this joint flavor fit exists.

next attempted nontrivial idea:
  Representation-controlled enlargement:
    A. Add the second independent 120_H doublet Clebsch direction, so the 120_H
       contribution to u,d,e,nu is not forced into one antisymmetric direction.
    B. Alternatively, relax exactly one of H or F by adding a controlled
       orthogonal symmetric tensor outside the Veronese subspace and measure
       the minimal norm required.
  Route A is preferred because it is a genuine Spin(10) Clebsch/doublet-mixing
  extension rather than an ad hoc geometric leak.

verification plan:
  1. Implement a two-120-direction Clebsch scan.
  2. Keep the same targets and viability thresholds.
  3. Compare mass/CKM Pareto front against the single-120 homotopy.
  4. If still not viable, run the orthogonal-symmetric relaxation as a no-go
     diagnostic and quantify the required non-Veronese norm.
```

Heartbeat update 2026-05-07 20:05 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The two-120-direction Clebsch scan is
  complete.  It improves the single-120 Pareto front but still does not produce
  a joint viable full-flavor point, so d=5 proton decay remains premature.

new script:
  code/scan_clebsch_flavor_two120.py

new outputs:
  output/flavor_clebsch_two120/scan_clebsch_flavor_two120.json
  output/flavor_clebsch_two120/scan_clebsch_flavor_two120_report.md
  output/flavor_clebsch_two120/scan_clebsch_flavor_two120.csv

derivation:
  Keep H,F in the CP1/O(-4) Veronese symmetric subspace, but add two
  independent antisymmetric 120-like family directions:
    G_A,G_B in Lambda^2 C^3.

  The tested ansatz is
    Y_u  = H + r_u F + a_u G_A + b_u G_B,
    Y_nu = H - 3 r_u F + a_nu G_A + b_nu G_B,
    Y_d  = r_d (H + F + a_d G_A + b_d G_B),
    Y_e  = r_d (H - 3 F + a_e G_A + b_e G_B).

numerical audit:
  ckm_lock:
    |V_us|=2.2349e-1, |V_cb|=4.5450e-2, |V_ub|=3.7057e-3,
    CKM score=2.012e-3, mass score=7.891e-1.
  mass_2:
    |V_us|=2.2166e-1, |V_cb|=4.4979e-2, |V_ub|=3.7136e-3,
    CKM score=1.663e-3, mass score=7.686e-1.
  mass_4:
    |V_us|=2.2023e-1, |V_cb|=4.5057e-2, |V_ub|=3.7384e-3,
    CKM score=1.786e-3, mass score=7.314e-1.
  mass_8:
    |V_us|=2.1859e-1, |V_cb|=4.5552e-2, |V_ub|=3.8571e-3,
    CKM score=2.574e-3, mass score=6.879e-1.
  mass_12:
    |V_us|=2.161031e-1, |V_cb|=4.616207e-2, |V_ub|=4.059424e-3,
    |J|=3.101923e-5,
    CKM score=4.580687e-3, mass score=6.440132e-1.

best point:
  Stage mass_12.
  Mass ratios:
    up:     (1.011343e-5, 3.458864e-2),
    down:   (9.584436e-4, 2.247952e-2),
    lepton: (3.138696e-4, 3.132297e-2),
    nu_D:   (1.816302e-2, 1.547414e-1).

  Seesaw replay:
    theta_norm = 1.612216e-11,
    seesaw residual = 9.186935e-13,
    (M1,M2,M3) =
      (1.502918e11, 1.758506e13, 4.734118e15) GeV.

checks:
  any CKM score < 0.05: yes.
  any mass score < 0.20: no.
  any joint viable point: no.
  improves mass vs single-120 homotopy: yes, 0.795235 -> 0.644013.
  best seesaw residual < 1e-10: yes.
  best theta_norm < 1e-8: yes.

paper update:
  Added "Two-120_H direction audit" to paper/gut_framework.tex.

current obstacle:
  The flavor blocker survived the minimal antisymmetric enlargement.  The
  second 120-like direction helps but does not remove the hierarchy tension:
  y_c/y_t and y_mu/y_tau remain especially displaced.  The obstruction is now
  likely in the symmetric Veronese restriction or in the too-simple doublet
  mixing ansatz.

next attempted nontrivial idea:
  Orthogonal-symmetric relaxation diagnostic:
    add one symmetric tensor S_perp orthogonal to the five-dimensional Veronese
    image, first to F and then to H, and measure the minimal
    ||S_perp||/||H,F|| required to reach mass score < 0.20 while keeping CKM
    score < 0.05 and the seesaw stable.

verification plan:
  1. Build the 6-dimensional Sym^2(C^3) basis and its 5-dimensional Veronese
     image from product Clebsch C_ij^m.
  2. Construct the normalized orthogonal symmetric direction S_perp.
  3. Run F+epsilon S_perp and H+epsilon S_perp scans seeded by the two-120
     best point.
  4. If |epsilon| is small, interpret it as a controlled geometric correction;
     if large, record a sharper no-go for strict CP1/O(-4) family geometry.
```

### P1: expose Yukawa/flavor benchmark data

The review was correct that the previous holomorphic overlap

```text
Y_ij = int psi_i psi_j h_a dmu_FS
```

was not manifestly bundle-covariant because `psi_i psi_j` is a section of
`O(4)`.  The adopted branch now treats the Higgs functional as an
`O(-4)`-valued dual density.  The remaining issue is reproducibility: the
paper must expose the exact dual coefficients and matrices, not only the
singular-value ratios.

Action items:

```text
1. Export the O(-4) dual-density coefficients a_m^(a).
2. Add Y_u,Y_d,Y_e,Y_nu matrices and basis conventions to the paper or a
   machine-readable benchmark card.
3. Keep the Hermitian Toeplitz construction as an optional Kahler/wavefunction
   branch, not as the default holomorphic superpotential.
4. Rerun seesaw, flavor, and proton checks from the exact matrices.
```

### P1: make seesaw/flavor numerically reproducible

The review is correct that the TeX gives heavy Majorana masses but not enough
data to reproduce them from the paper alone.  The local JSON contains matrices
and phases, but the manuscript does not expose them.

Action items:

```text
1. Add m_D, m_nu, M_R, PMNS phases, basis convention, and Takagi convention.
2. Provide exact matrices in an appendix or machine-readable benchmark card.
3. Expand from singular-value ratios to a full flavor fit:
   quark masses, charged-lepton masses, CKM, neutrino angles, phases.
```

Note: the current TeX already uses `overline{126}_H` and
`overline{16}_H overline{16}_H`; the review's conjugation warning appears to
refer to an older/looser notation, but representation charges should still be
spelled out more explicitly.

### P1: action-level model definition and reproducibility package

The review is correct that the manuscript lacks a complete action/superpotential
and a reproducibility appendix.

Action items:

```text
1. Add W = W_Yukawa + W_Majorana + W_Higgs + W_mediator + W_PS-breaking.
2. State Kähler/gauge-kinetic assumptions or explicitly keep the EFT scope.
3. Add a regeneration appendix: scripts, input JSON/CSV, random seeds, grids,
   thresholds, beta coefficients, tolerances, and benchmark cards.
4. Add modern references for PDG inputs, NuFIT/neutrino data, proton decay,
   lattice matrix elements, SUSY d=5 dressing, and SO(10) model-building.
```

### P2: title and theorem-scope cleanup

The review is correct at the rhetoric level.  The title and theorem language
should avoid claiming a first-principles derivation of a unique microscopic
GUT.

Action items:

```text
1. Rename the paper around "Conditional Spin(10) Effective Framework".
2. Convert the current scope theorem into an assumptions/scope map.
3. Clearly separate PSLT/Another-Physics motivation from Spin(10) EFT inputs.
```

Heartbeat update 2026-05-07 20:19 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the immediate flavor obstruction
  has advanced from a strict no-go to a controlled deformation target.

new script:
  code/scan_clebsch_flavor_symmetric_relaxation.py

new outputs:
  output/flavor_symmetric_relaxation/scan_clebsch_flavor_symmetric_relaxation.json
  output/flavor_symmetric_relaxation/scan_clebsch_flavor_symmetric_relaxation_report.md
  output/flavor_symmetric_relaxation/scan_clebsch_flavor_symmetric_relaxation.csv

derivation:
  Sym^2(C^3) has complex dimension 6, while the CP1/O(-4) Veronese product
  image has rank 5.  With Frobenius pairing <A,B>=Tr(A^dagger B), there is a
  unique unit symmetric direction S_perp orthogonal to the five Clebsch tensors
  C_m:

    <C_m,S_perp> = 0, ||S_perp||_F = 1, S_perp^T = S_perp.

  The scan adds only this one new direction:

    H -> H + epsilon_H S_perp, or F -> F + epsilon_F S_perp,

  while retaining the two independent antisymmetric 120-like directions.

numerical audit:
  Geometric checks:
    rank(Veronese) = 5,
    max_m |<C_m,S_perp>| = 3.403e-16,
    ||S_perp - S_perp^T||_F = 0.

  Representative viable points:
    F branch, eps budget 0.30, mass_40:
      |V_us|=0.164128, |V_cb|=0.045869, |V_ub|=0.004593,
      CKM score=2.995e-2, mass score=1.813e-1,
      |epsilon|=3.740e-4, relative leakage=1.237e-2.

    F branch, eps budget 1.00, mass_12:
      |V_us|=0.221820, |V_cb|=0.041507, |V_ub|=0.003721,
      CKM score=7.249e-5, mass score=1.614e-1,
      |epsilon|=4.288e-4, relative leakage=1.303e-2.

    H branch, eps budget 0.30, mass_40:
      |V_us|=0.178719, |V_cb|=0.044810, |V_ub|=0.004195,
      CKM score=1.447e-2, mass score=1.748e-1,
      |epsilon|=7.553e-4, relative leakage=7.443e-4.

  Best stored point:
    H branch, mass_40, eps budget 0.30.
    theta_norm = 4.340850e-11,
    seesaw residual = 9.305537e-12,
    (M1,M2,M3) =
      (2.446598e10, 2.840222e13, 4.599321e15) GeV.

checks:
  any CKM score < 0.05: yes.
  any mass score < 0.20: yes.
  any joint viable point: yes.
  best seesaw residual < 1e-10: yes.
  best theta_norm < 1e-8: yes.

paper update:
  Added "Orthogonal symmetric relaxation diagnostic" to paper/gut_framework.tex.

current obstacle:
  Strict CP1/O(-4) Veronese flavor geometry is too rigid, but only one
  orthogonal symmetric direction is enough to pass the finite flavor scan.
  This is now a conditional deformation, not yet a first-principles derivation.

next attempted nontrivial idea:
  Derive S_perp covariantly instead of inserting it by hand.  Candidate origins:
    1. curvature/contact correction in the O(-4)-valued dual density;
    2. a controlled higher-dimensional family operator projected onto
       Sym^2(C^3)/Veronese_5;
    3. a second CP1 patch/defect source whose residue is forced into the
       orthogonal complement by a selection rule.

verification plan:
  1. Build an explicit covariant O(-4) contact functional whose projection is
     S_perp and prove coordinate independence.
  2. Rerun the flavor scan with epsilon fixed by that microscopic parameter,
     not free.
  3. Export the viable Yukawa matrices and left/right rotations.
  4. Feed those rotations into p -> K+ nu_bar dimension-5 Wilson coefficients.
```

Heartbeat update 2026-05-07 22:01 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The immediate S_perp-origin blocker is
  partially resolved: the fitted orthogonal symmetric direction is exactly the
  CP1 O(2) second transvectant/contact invariant, up to phase.

new script:
  code/derive_sperp_transvectant.py

new outputs:
  output/flavor_sperp_transvectant/derive_sperp_transvectant.json
  output/flavor_sperp_transvectant/derive_sperp_transvectant_report.md

derivation:
  Use H^0(CP1,O(2)) ~= Sym^2(C^2) with normalized spin-one basis

    u0 = x^2, u1 = sqrt(2) x y, u2 = y^2.

  Then

    Sym^2(Sym^2 C^2) = Sym^4 C^2 + Sym^0 C^2.

  The ordinary O(-4) Veronese product probes the Sym^4 C^2 piece.  The missing
  one-dimensional orthogonal direction is the spin-zero second transvectant:

    K_tr = 1/sqrt(3) [[0,0,-1],[0,1,0],[-1,0,0]],
    T_2(s,s) = a1^2 - 2 a0 a2.

  Pointwise, this is exactly the binary quadratic identity

    u1^2 - 2 u0 u2 = 0,

  so it is invisible to the multiplication/O(-4)-density image but available
  as a contact/trace functional.

numerical audit:
  phase-aligned residual between fitted S_perp and analytic K_tr:
    3.699296e-16.

  max Veronese overlap:
    3.915193e-17.

  pointwise product kernel norm:
    0.

  sampled SU(2) covariance check:
    max_g ||D(g)^T K_tr D(g) - K_tr||_F = 1.097056e-15.

paper update:
  Added "Transvectant origin of S_perp" to paper/gut_framework.tex.

current obstacle:
  S_perp is now covariantly identified, but its small coefficient is still a
  free deformation parameter.  The theory still needs an action-level
  selection rule/source mechanism that generates epsilon_H ~ 7.6e-4 in the
  viable branch without spoiling threshold/proton safety.

next attempted nontrivial idea:
  Promote the transvectant into an explicit family-contact operator:

    W_Yukawa contains
      16_i 16_j H [lambda_4 M_ij[h_{-4}] + lambda_tr K_tr,ij].

  Then either:
    A. derive lambda_tr/lambda_4 from a localized contact/source sector on CP1;
    B. constrain it by a grading/R-symmetry as the leading allowed irrelevant
       family operator;
    C. fit only the scalar coefficient epsilon and export full Yukawa
       rotations for d=5 proton decay.

verification plan:
  1. Add a benchmark card for the viable H-branch transvectant point with
     exact Y_u,Y_d,Y_e,Y_nu and left/right rotations.
  2. Check whether the needed epsilon_H is stable under small changes in the
     CKM/mass weights.
  3. Feed the resulting rotations into the p -> K+ nu_bar dimension-5 Wilson
     coefficient scan.
```

Heartbeat update 2026-05-07 22:39 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The transvectant flavor branch now has
  exact mass-basis rotations and a first dimension-five proton-decay flavor
  leakage proxy, but not yet a full dressed p -> K+ nu_bar Wilson calculation.

new script:
  code/export_transvectant_flavor_rotations.py

new outputs:
  output/flavor_transvectant_rotations/transvectant_flavor_rotations.json
  output/flavor_transvectant_rotations/transvectant_flavor_rotations_report.md
  output/flavor_transvectant_rotations/dimension5_flavor_proxy.csv

derivation:
  For each fitted Yukawa matrix Y_a, compute biunitary rotations

    U_{a,L}^dagger Y_a U_{a,R} = diag(y_{a1},y_{a2},y_{a3})

  with ascending singular values.  Then replay

    V_CKM = U_{u,L}^dagger U_{d,L}

  and export the PMNS convention used by the local seesaw replay.

  For the first d=5 proxy, keep the local item-4 normalization:

    y_t = 0.60, y_b = 0.024,
    M_T = 1e16 GeV, m_wino = 1e3 GeV, m_sfermion = 1e5 GeV,
    alpha_2^{-1} = 25.

  Estimate the dressed coefficient only as

    C6^(5) ~= (alpha_2/4pi)(m_wino/m_sfermion^2)(S_T y_eff/M_T),

  where y_eff includes the fitted CKM leakage.  This is a proxy, not the final
  channel calculation.

numerical audit:
  rotation checks:
    max biunitary off-diagonal residual = 1.145062e-15,
    CKM replay residual = 3.102074e-13.

  d=5 proxy:
    first-generation raw:
      y_eff = 1.445e-10, S_T max for tau>1e34 yr = 9.14e3.
    second-generation raw:
      y_eff = 4.419e-6, S_T max = 2.99e-1.
    second-generation CKM-to-K proxy:
      leakage = 1.756e-1, y_eff = 7.762e-7, S_T max = 1.70.
    third-generation raw:
      y_eff = 1.440e-2, S_T max = 9.17e-5.
    third-generation CKM-to-K proxy:
      |V_td V_ts| = 2.048622e-4,
      y_eff = 2.950015e-6,
      S_T max for tau>1e34 yr = 4.476948e-1,
      S_T max for tau>2.4e34 yr = 2.889858e-1.
    third-generation single-Cabibbo stress:
      leakage = 1.787e-1,
      y_eff = 2.574e-3,
      S_T max for tau>1e34 yr = 5.132e-4.

paper update:
  Added "Transvectant flavor rotations and d=5 proxy" to paper/gut_framework.tex.

current obstacle:
  The old scalar S_T-only estimate is now superseded for the viable flavor
  point, but the real p -> K+ nu_bar Wilson tensor still needs triplet mixing,
  wino/higgsino dressing, and flavor rotations inserted consistently.

next attempted nontrivial idea:
  Build the full LLLL and RRRR dimension-five Wilson tensor in the exported
  mass basis:

    C_L^{abcd} ~ (Y_T^{QQ})_{ij}(Y_T^{QL})_{kl}
      U_Q^{ia} U_Q^{jb} U_Q^{kc} U_L^{ld} / M_T,

  with several triplet-Clebsch hypotheses:
    A. minimal SU(5)-like symmetric triplet,
    B. Spin(10) 10+126bar+120 Clebsch triplet aligned with the fitted doublets,
    C. transvectant-filtered triplet coupling.

verification plan:
  1. Generate C_L and C_R tensors in the mass basis for the three triplet
     hypotheses.
  2. Extract p -> K+ nu_bar and K0 mu+ proxy amplitudes before dressing.
  3. Add wino/higgsino dressing loop functions and compare S_T requirements.
  4. Feed the resulting S_T bound back to the threshold/proton-safe branch.
```

Heartbeat update 2026-05-08 00:18 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The transvectant flavor point now has
  a reproducible four-index mass-basis dimension-five Wilson-tensor ledger, but
  not yet the final channel-by-channel dressed p -> K+ nu_bar calculation with
  triplet mixing and chiral reduction.

new script:
  code/scan_dimension5_wilson_tensors.py

new outputs:
  output/dimension5_wilson_tensors/dimension5_wilson_tensors.json
  output/dimension5_wilson_tensors/dimension5_wilson_tensors_report.md
  output/dimension5_wilson_tensors/dimension5_wilson_proxy.csv

derivation:
  For each triplet-Clebsch hypothesis define

    C_L^{abcd}
      = (Y_QQ)_{ij}(Y_QL)_{kl}
        (U_1)_{ia}(U_2)_{jb}(U_3)_{kc}(U_L)_{ld},

  and audit both doublet assignments

    (U_1,U_2,U_3,U_L) = (U_u,U_u,U_d,U_nu),
    (U_d,U_d,U_u,U_nu).

  The right-handed proxy is

    C_R^{abcd}
      = (Y_UE)_{ij}(Y_UD)_{kl}
        (U_{u,R})_{ia}(U_{e,R})_{jd}
        (U_{u,R})_{kb}(U_{d,R})_{lc},

  stored as (u_a,u_b,d_c,e_d).  The audited triplet hypotheses are:
    minimal down-aligned, lepton-transposed, geometric average,
    antisymmetric up piece, and transvectant-contact QQ.

numerical audit:
  All rows use the local item-4 dressing proxy

    M_T = 1e16 GeV, m_wino = 1e3 GeV, m_sfermion = 1e5 GeV,
    alpha_2^{-1}=25,

  and physical normalization

    y_t=0.60, y_b=0.024, y_tau(proxy)=0.010.

  Most dangerous rows:
    transvectant-contact LLLL up-up-down Knu:
      |C_L| = 2.263424e-3,
      S_T max for tau>2.4e34 yr = 3.766472e-4.
    minimal down-aligned LLLL up-up-down Knu:
      |C_L| = 2.238681e-3,
      S_T max = 3.808101e-4.
    lepton-transposed LLLL Knu:
      |C_L| = 9.041641e-4,
      S_T max = 9.428736e-4.
    antisymmetric up-piece LLLL:
      |C_L| = 9.827458e-6,
      S_T max = 8.674800e-2.
    minimal down-aligned RRRR:
      |C_R| = 1.206154e-6,
      S_T max = 7.068025e-1.

paper update:
  Added "Mass-basis d=5 Wilson-tensor audit" to paper/gut_framework.tex.

current obstacle:
  The scalar S_T proxy is superseded.  Explicit tensor orientation can change
  the leading d=5 flavor entry by more than two orders of magnitude, so a
  publishable proton-safe claim must specify triplet Clebsch orientation,
  triplet mixing, wino/higgsino dressing functions, and chiral/hadronic
  reduction.

next attempted nontrivial idea:
  Replace the common scalar S_T by a small triplet-mixing matrix in the
  10_H + 126bar_H + 120_H + transvectant-contact basis, then scan whether a
  symmetry-compatible null vector can suppress only the LLLL Knu tensor while
  preserving the successful doublet Yukawa and seesaw fit.

verification plan:
  1. Build a triplet-mixing matrix T_AB and compute
       C_L = sum_AB Y_QQ^A (T^{-1})_{AB} Y_QL^B
     and similarly for C_R.
  2. Add wino and higgsino dressing loop functions separately.
  3. Extract p -> K+ nu_bar, K0 mu+, e+ pi0-compatible d=5 amplitudes.
  4. Replay the resulting S_T or T_AB bound through the corrected
     threshold/proton branch.
```

Heartbeat update 2026-05-08 01:14 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The scalar triplet filter has been
  upgraded to a finite inverse triplet-mixing matrix W_AB=(T^{-1})_AB, and the
  first nullspace audit is complete.  It gives a constructive near-null
  direction for LLLL Knu, but not a natural action-level derivation.

new script:
  code/scan_triplet_mixing_nullspace.py

new outputs:
  output/triplet_mixing_nullspace/triplet_mixing_nullspace_summary.json
  output/triplet_mixing_nullspace/triplet_mixing_nullspace_report.md
  output/triplet_mixing_nullspace/triplet_mixing_nullspace_scan.csv

derivation:
  Use the four local triplet-source shapes

    H10 ~ sym(Y_u),
    F126bar ~ sym(Y_d - Y_e^T) / (Y_d - Y_e^T),
    G120 ~ anti(Y_u) / anti(Y_d),
    Ktr ~ CP1 transvectant contact.

  For W_AB=(T^{-1})_AB,

    C_L^{abcd}(W)
      = sum_AB W_AB (Y_QQ^A)_{ij}(Y_QL^B)_{kl}
        (U_1)_{ia}(U_2)_{jb}(U_3)_{kc}(U_L)_{ld}.

  The Knu cancellation problem is a finite SVD problem:

    A_{r,(A,B)} = C_{L,r}^{AB},

  where r runs over all uus nu_i entries for both doublet assignments
  (U_u,U_u,U_d,U_nu) and (U_d,U_d,U_u,U_nu).

numerical audit:
  Diagonal representation mixing:
    rank(A_diag)=4,
    sigma_min(A_diag)=5.359441e-4,
    no exact Knu null.
    Identity max Knu amplitude = 3.296375e-3.
    Least-singular max Knu amplitude = 2.833496e-4.
    Identity overlap of least-singular direction = 4.864907e-2.

  Full 4x4 bipartite mixing:
    rank(A_full)=16,
    sigma_min(A_full)=1.537907e-7,
    no exact all-flavor Knu null, but a sharp near-null exists.
    Identity max Knu amplitude = 3.296375e-3.
    Near-null max Knu amplitude = 9.906694e-8.
    Suppression factor = 3.327422e4.
    Identity overlap = 4.794355e-1,
    angle from identity = 1.070785 rad.
    W rank = 3, so the near-null is a codimension-one triplet inverse-mass
      condition, not a generic perturbation.
    Leading post-null channel = RRRR_uusd_anycharged,
      amplitude = 2.831289e-4,
      S_T max for tau>2.4e34 yr = 3.011040e-3.

paper update:
  Added "Triplet-mixing near-null audit" to paper/gut_framework.tex.

current obstacle:
  The near-null line is mathematically real but not yet natural: it requires a
  rank-3 inverse triplet-mass matrix and shifts the leading obstruction to the
  RRRR proxy.  A real model must derive this rank condition from the
  triplet-sector superpotential or add an independent RRRR filter.

next attempted nontrivial idea:
  Construct an explicit missing-partner/rank-one-lift triplet superpotential

    W_T = T_A M_AB \bar T_B + X p_A T_A + \bar X q_B \bar T_B

  whose Schur complement enforces the near-null W_AB direction while keeping
  doublet Yukawa textures and the successful threshold branch unchanged.

verification plan:
  1. Fit a minimal low-rank mass matrix T_AB whose inverse projects onto the
     near-null W_AB up to controlled epsilon.
  2. Check F-flatness/D-flatness of the triplet-rank sector and whether the
     needed rank deficiency is symmetry protected.
  3. Recompute C_L and C_R with this explicit W_AB, then add separate wino and
     higgsino dressing loop functions.
  4. Feed the resulting d=5 bound back to the corrected RGE/proton scan.
```

Heartbeat update 2026-05-08 02:53 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The triplet near-null now has an
  explicit rank-one-lift superpotential interpretation, but the realization
  reveals a heavy-spectrum hierarchy obstruction and an unresolved RRRR
  proton-decay obstruction.

new script:
  code/construct_triplet_rank_lift.py

new outputs:
  output/triplet_rank_lift/triplet_rank_lift_summary.json
  output/triplet_rank_lift/triplet_rank_lift_report.md
  output/triplet_rank_lift/triplet_rank_lift_scan.csv

derivation:
  Let

    W_near = U diag(s1,s2,s3,0) V^dagger

  be the rank-three inverse triplet propagator from the near-null audit.  A
  finite rank-one lift is

    W_eps = U diag(s1,s2,s3,eps*s1) V^dagger,
    M_eps = V diag(1/s1,1/s2,1/s3,1/(eps*s1)) U^dagger,

  so M_eps^{-1}=W_eps.  This can be written at the action level as

    W_T = T_A (M_r3)_{AB} barT_B
          + M_* (v_4^dagger T)(u_4^dagger barT),

  or with a mediator pair

    W_T = T M_r3 barT + X p_A T_A + q_B barT_B barX + M_X X barX,

  whose Schur complement produces the rank-one update.  For eps>0, M_eps is
  full rank, so the zero-triplet vacuum satisfies F_T=M_eps barT=0 and
  F_barT=M_eps^T T=0 only at T=barT=0; D-terms vanish at the origin.

numerical audit:
  The near-null singular values are

    s = (9.999965e-1, 2.640176e-3, 2.738209e-4, 0),

  treating the fourth singular value as numerical zero.  The finite rank-three
  mass ratios are

    1, 3.787613e2, 3.652009e3,

  plus one infinite/lifted direction.  If M_light=1e16 GeV, the largest finite
  triplet is 3.652009e19 GeV, above the reduced-Planck benchmark used in the
  local audit.

  Representative channel/mass checks:
    rank-three limit:
      max Knu amplitude = 9.906694e-8,
      RRRR amplitude = 2.831289e-4,
      max mass ratio = 3.652009e3.
    eps=1e-4 lift:
      max Knu amplitude = 7.493234e-7,
      RRRR amplitude = 2.831736e-4,
      max mass ratio = 1e4.
    condition cap kappa=100:
      max Knu amplitude = 8.838074e-5,
      RRRR amplitude = 2.884973e-4,
      max mass = 1e18 GeV for M_light=1e16 GeV,
      leading S_T max for tau>2.4e34 yr = 2.955010e-3.

paper update:
  Added "Rank-one lift of the triplet near-null" to paper/gut_framework.tex.

current obstacle:
  The rank-one lift is a valid EFT implementation of the near-null, but not a
  natural UV solution by itself.  The exact near-null wants a triplet mass
  hierarchy above the Planck benchmark if M_light is 1e16 GeV, while the
  Planck-safe condition-capped branch still leaves RRRR as the leading d=5
  obstruction and needs a filter at the 3e-3 level unless the existing
  S_T=1e-5 triplet filter is retained.

next attempted nontrivial idea:
  Build a two-sided missing-partner/selection-rule sector that filters RRRR
  independently of LLLL:

    W_T = T_A M_AB barT_B + X_L p_A T_A + X_R r_A T_A
          + q_B barT_B barX_L + s_B barT_B barX_R

  with one rank condition aligned to the LLLL Knu near-null and a second
  grading that suppresses the u^c e^c u^c d^c source direction.

verification plan:
  1. Extract the dominant RRRR source vector from the rank-lift output.
  2. Solve a joint constrained SVD for W_AB minimizing both Knu and RRRR under
     condition-number caps kappa <= 30, 100, 300.
  3. Translate the resulting two-sided projection into a minimal charge table
     or missing-partner superpotential.
  4. Replay the corrected d=5 amplitudes through the threshold/proton scan.
```

Heartbeat update 2026-05-08 04:09 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  A two-sided triplet filter scan now
  suppresses LLLL Knu and RRRR simultaneously under finite condition-number
  caps, but the preferred direction is still found by SVD rather than derived
  from a symmetry or explicit superpotential grading.

new script:
  code/scan_two_sided_triplet_filter.py

new outputs:
  output/two_sided_triplet_filter/two_sided_triplet_filter_summary.json
  output/two_sided_triplet_filter/two_sided_triplet_filter_report.md
  output/two_sided_triplet_filter/two_sided_triplet_filter_scan.csv

derivation:
  Let A_K be the row matrix for all monitored uus nu_i entries in both doublet
  assignments, and A_R the row matrix for u^c u^c s^c e^c entries.  The scan
  solves

    min_{||W||_F=1} || ( Ahat_K ; omega_R Ahat_R ) vec(W) ||_2,

  where the hats denote Frobenius-normalized blocks.  Each candidate W is then
  regularized by singular-value flooring

    s_i(W) -> max(s_i(W), s_1(W)/kappa_max),

  so that the corresponding triplet hierarchy satisfies

    M_max/M_min <= kappa_max.

numerical audit:
  Block norms before normalization:
    ||A_Knu||_F = 3.725380e-2,
    ||A_RRRR||_F = 2.976942e-2,
    ||A_K0mu||_F = 1.357473e-2.

  Best Knu/RRRR weighted direction for kappa caps 30,100,300 is always
  omega_R = 0.1.

  kappa=30:
    max Knu = 2.481441e-4,
    K0mu = 1.598773e-4,
    RRRR = 2.940061e-4,
    worst(Knu,RRR) = 2.940061e-4,
    M_max = 3.0e17 GeV for M_min=1e16 GeV,
    leading S_T max for tau>2.4e34 yr = 2.899642e-3.

  kappa=100:
    max Knu = 7.455680e-5,
    K0mu = 4.815366e-5,
    RRRR = 7.709749e-5,
    worst(Knu,RRR) = 7.709749e-5,
    M_max = 1.0e18 GeV for M_min=1e16 GeV,
    leading S_T max for tau>2.4e34 yr = 1.105759e-2.

  kappa=300:
    max Knu = 2.410854e-5,
    K0mu = 1.530185e-5,
    RRRR = 2.347488e-5,
    worst(Knu,RRR) = 2.410854e-5,
    M_max = 3.0e18 GeV for M_min=1e16 GeV,
    leading S_T max for tau>2.4e34 yr = 3.536143e-2.

  All-block guard including K0mu:
    kappa=100 gives
      max Knu = 8.224554e-5,
      K0mu = 4.200167e-5,
      RRRR = 9.612328e-5,
      worst all = 9.612328e-5.
    Thus K0mu is not the new leading obstruction.

paper update:
  Added "Two-sided triplet filter" to paper/gut_framework.tex.

current obstacle:
  The two-sided filter works numerically as an EFT direction, especially at
  kappa=100, but it is still selected by an SVD weight omega_R=0.1.  The model
  needs an action-level reason for that weight/direction and a true
  wino/higgsino dressed calculation before claiming proton safety.

next attempted nontrivial idea:
  Translate the omega_R=0.1, kappa=100 filter into a minimal charge-table or
  missing-partner superpotential.  The target is to allow the H10/F126/Ktr
  pairs that dominate the good W_AB direction while forbidding or lifting the
  source combinations that regenerate the RRRR row.

verification plan:
  1. Extract the dominant W_AB entries for the kappa=100 optimum.
  2. Solve a small integer charge-assignment problem for T_A, barT_B, and
     mediator pairs so allowed bilinears reproduce that support pattern.
  3. Recompute W_AB from the allowed/lifted mass matrix and replay the
     two-sided d=5 amplitudes.
  4. Add explicit wino/higgsino dressing functions and feed the resulting
     proton bound back to the corrected RGE threshold scan.
```

Heartbeat update 2026-05-08 04:52 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The d=5 triplet filter has now been
  pushed one step beyond a continuous SVD direction: a minimal cyclic grading
  reproduces the preferred kappa=100, omega_R=0.1 support exactly.  However,
  the exact support is rank deficient, so the remaining task is to realize the
  needed kappa=100 floor from an explicit regulator superpotential rather than
  inserting it phenomenologically.

new script:
  code/derive_triplet_filter_charge_table.py

new outputs:
  output/triplet_filter_charge_table/triplet_filter_charge_table_summary.json
  output/triplet_filter_charge_table/triplet_filter_charge_table_report.md
  output/triplet_filter_charge_table/triplet_filter_charge_scan.csv
  output/triplet_filter_charge_table/triplet_filter_charge_table.csv
  output/triplet_filter_charge_table/triplet_filter_allowed_support.csv

derivation:
  Starting from the omega_R=0.1, kappa=100 two-sided filter, define dominant
  and suppressed entries by

    |W_AB| >= 4.0e-2,     |W_AB| <= 2.0e-2.

  The charge scan searches Z_N gradings, N=2,...,8, with q(H10)=0 gauge-fixed
  and a neutral mass/link spurion.  The best exact support match is a single
  Z_3 grading:

                H10   F126   G120   Ktr
    q(T_A)       0      0      1      0
    q(barT_B)    0      0      1      1

  with allowed entries satisfying

    q(T_A) + q(barT_B) = 0 mod 3.

  Therefore the allowed support is exactly

    (H10, F126, Ktr)_left x (H10, F126)_right,

  while G120-related entries and the right Ktr column are forbidden at the
  neutral-spurion level.

numerical audit:
  False negatives among dominant entries: 0.
  False positives among suppressed entries: 0.

  The support-only matrix is rank deficient.  Replaying the same missing-partner
  support with a kappa=100 universal singular floor gives

    max Knu amplitude = 8.098721968e-5,
    K0mu amplitude    = 6.626369280e-5,
    RRRR amplitude    = 1.066375671e-4,
    leading channel   = RRRR_uusd_anycharged,
    S_T max for tau>2.4e34 yr = 7.994484708e-3,
    final condition number = 100.0,
    M_max = 1.0e18 GeV for M_min=1.0e16 GeV.

paper update:
  Added "Triplet filter charge table" to paper/gut_framework.tex.

current obstacle:
  The Z_3 selection rule explains the sparse support pattern, but the exact
  support is singular.  The one-percent kappa=100 floor must be generated by
  a controlled heavy regulator or small symmetry-breaking spurion without
  regenerating the forbidden G120/right-Ktr entries at order larger than the
  proton-safe leakage budget.  A true wino/higgsino dressed proton-decay
  calculation is still pending.

next attempted nontrivial idea:
  Construct a Z_3-compatible missing-partner regulator superpotential using a
  neutral spurion S_0 for the allowed block plus minimal charged spurions or
  mediator pairs that lift only the null directions.  The desired action-level
  target is a finite mass matrix whose inverse reproduces the support above
  and whose smallest singular value is s_1/100, while all forbidden entries
  remain <= O(1e-2) of the allowed block.

verification plan:
  1. Build the explicit Z_3 regulator mass matrix with spurion charges.
  2. Check rank, singular spectrum, condition number, and leakage into G120 and
     right-Ktr columns.
  3. Replay the mass-basis d=5 Wilson tensors for Knu, K0mu, and RRRR.
  4. If the proxy survives, add wino/higgsino dressing functions and feed the
     resulting proton bound back into the corrected RGE/threshold scan.
```

Heartbeat update 2026-05-08 07:06 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The previous obstacle--the singular
  Z_3 triplet-filter support--has been sharpened into an explicit finite
  regulator candidate.  The kappa=100 floor can be interpreted as a compact
  two-spurion Z_3 missing-partner lift rather than a generic forbidden-entry
  leakage, but the spurion stabilization and radiative filling of the remaining
  charged entries are still open.

new script:
  code/construct_z3_triplet_regulator.py

new outputs:
  output/z3_triplet_regulator/z3_triplet_regulator_summary.json
  output/z3_triplet_regulator/z3_triplet_regulator_report.md
  output/z3_triplet_regulator/z3_triplet_regulator_scan.csv
  output/z3_triplet_regulator/z3_triplet_regulator_matrix_kappa100.csv

derivation:
  Use the charge table

                H10   F126   G120   Ktr
    q(T_A)       0      0      1      0
    q(barT_B)    0      0      1      1

  and decompose every effective triplet inverse-propagator entry by

    r_AB = q(T_A) + q(barT_B) mod 3.

  Neutral entries have r=0.  Entries with r=1 and r=2 require spurions of
  charge 2 and 1, respectively.  The regulator superpotential is

    W_T/M_T = T_A [ S0 N_AB + S2 R1_AB + S1 R2_AB ] barT_B,

  where N has only r=0 entries, R1 only r=1 entries, and R2 only r=2 entries.
  Algebraically, this is the condition-capped singular-value lift

    W_kappa = U diag(max(s_i, s_1/kappa)) V^dagger

  applied to the neutral-support matrix.

numerical audit:
  Neutral-support singular values:

    9.999988e-1, 1.552131e-3, 0, 0.

  For kappa=100 the finite regulator has singular values

    9.998500e-1, 9.998500e-3, 9.998500e-3, 9.998500e-3,

  so det(W_100) = 9.994003e-7 and condition number = 100.0.

  Nonzero charged entries above 1e-10:

    H10-barKtr   residue 1, spurion charge 2, |W|=5.910e-3,
    F126-barKtr  residue 1, spurion charge 2, |W|=5.128e-3,
    Ktr-barKtr   residue 1, spurion charge 2, |W|=6.224e-3,
    G120-barG120 residue 2, spurion charge 1, |W|=9.999e-3.

  Residue Frobenius norms at kappa=100:

    neutral r=0: 9.999000e-1,
    r=1:         9.998500e-3,
    r=2:         9.998500e-3.

  Charged/neutral max-entry ratio = 1.988252e-2.
  Charged/neutral Frobenius ratio = 1.414143e-2.

  d=5 replay at kappa=100:

    max Knu amplitude = 8.098721968e-5,
    K0mu amplitude    = 6.626369280e-5,
    RRRR amplitude    = 1.066375671e-4,
    leading channel   = RRRR_uusd_anycharged,
    S_T max for tau>2.4e34 yr = 7.994484708e-3.

  F/D-flatness at the triplet origin:
    Since det(W_100) != 0, F_T = W_100 barT = 0 and
    F_barT = W_100^T T = 0 force T=barT=0.  Triplet D-terms vanish at the
    origin.  Singlet-spurion stabilization remains separate.

paper update:
  Added "Z3-compatible triplet regulator" to paper/gut_framework.tex.

current obstacle:
  The regulator is now action-level at the effective triplet superpotential
  support level, but S1 and S2 have not yet been stabilized dynamically.
  Higher-dimensional operators or loops could fill the other residue-1 and
  residue-2 entries above the percent budget unless a driver/selection sector
  suppresses them.

next attempted nontrivial idea:
  Build an explicit spurion-stabilization sector for S1 and S2 with F-flat
  driving fields, then audit all Z_3-allowed higher operators up to the first
  dangerous order.  The target is to prove that the four required charged
  entries are generated at O(1e-2), while all other charged entries remain
  <= O(1e-4) in W_AB or are aligned away by a rank-one driver.

verification plan:
  1. Enumerate Z_3-allowed spurion monomials S1^a S2^b up to low order for all
     triplet entries.
  2. Construct a minimal driver superpotential fixing <S1>/Lambda and
     <S2>/Lambda near 1e-2 without creating triplet vevs.
  3. Add the induced higher-order leakage matrix to W_AB and replay Knu, K0mu,
     and RRRR amplitudes.
  4. If stable, propagate the finite triplet spectrum into the full
     threshold/proton scan.
```

Heartbeat update 2026-05-08 07:23 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The charged-spurion magnitudes can be
  stabilized by a standard F-flat singlet driver, but the Z_3 symmetry alone
  does not align the triplet regulator matrix.  A generic degree-one spurion
  sector refills unused charged entries, so a matrix-alignment driver or extra
  shaping symmetry is still needed for a paper-level derivation of the sparse
  support.

new script:
  code/audit_z3_spurion_leakage.py

new outputs:
  output/z3_spurion_leakage/z3_spurion_leakage_summary.json
  output/z3_spurion_leakage/z3_spurion_leakage_report.md
  output/z3_spurion_leakage/z3_spurion_leakage_bounds.csv
  output/z3_spurion_leakage/z3_spurion_monomials.csv
  output/z3_spurion_leakage/z3_triplet_entry_classification.csv

derivation:
  Use q(S1)=1 and q(S2)=2.  A term T_A barT_B S1^a S2^b is invariant when

    q(T_A)+q(barT_B)+a+2b = 0 mod 3.

  Thus degree-one monomials already generate every charged residue:

    S2 generates all r_AB=1 entries,
    S1 generates all r_AB=2 entries.

  A minimal singlet driving sector for the spurion magnitudes is

    W_drive =
      X1(S1^3-v1^3) + X2(S2^3-v2^3) + X12(S1 S2-v12^2),

  with R(Xi)=2 and <Xi>=0.  Then F_X=0 fixes

    S1^3=v1^3,  S2^3=v2^3,  S1 S2=v12^2,

  while F_S1=F_S2=0 at Xi=0.  Consistency requires

    v12^6 = v1^3 v2^3

  up to the chosen Z_3 branch.

dangerous entries:
  The regulator wants four charged entries.  Six additional entries can be
  generated by generic Z_3 spurions:

    H10-barG120, F126-barG120, G120-barH10, G120-barF126,
    Ktr-barG120, G120-barKtr.

numerical audit:
  If every dangerous charged entry obeys |delta W_AB| <= eta, the script uses
  the deterministic row-wise bound

    |C| <= |C0| + eta sum_{dangerous j} |a_j|

  for each monitored proton-decay row C=a.W.

  Bounds:

    eta=0:
      max Knu = 8.098721968e-5,
      K0mu    = 6.626369280e-5,
      RRRR    = 1.066375671e-4.

    eta=1e-4:
      max Knu = 8.286107e-5,
      K0mu    = 6.633639e-5,
      RRRR    = 1.074655e-4.

    eta=1e-3:
      max Knu = 9.972571e-5,
      K0mu    = 7.173965e-5,
      RRRR    = 1.149165e-4.

    eta=1e-2:
      max Knu = 2.710820e-4,
      K0mu    = 1.970695e-4,
      RRRR    = 2.162678e-4.

  Keeping the worst monitored amplitude below 2.0e-4 requires

    eta <= 6.270851e-3.

paper update:
  Added "Z3 spurion leakage audit" to paper/gut_framework.tex.

current obstacle:
  Z_3 plus singlet F-flatness stabilizes the spurion sizes but not the matrix
  direction.  Proton-amplitude safety only requires mild suppression of unused
  charged entries, but a clean model-building derivation of the sparse triplet
  regulator needs a stronger matrix-alignment mechanism.

next attempted nontrivial idea:
  Construct a rank-one/rank-three matrix-alignment driver for the charged
  residue sectors.  The target is:

    R1 support -> right-Ktr column projected onto (H10,F126,Ktr)_left,
    R2 support -> G120-barG120 only,

  while H/F/G/K leakage into the six unused charged entries is either forbidden
  by a second grading or suppressed below the desired structural purity target.

verification plan:
  1. Introduce alignment fields A1,A2 and driving fields Y1,Y2 with charges
     chosen to produce the desired R1 and R2 projectors.
  2. Prove F-flatness of the aligned matrices and count residual moduli.
  3. Add the aligned leakage matrix to W_AB and replay the deterministic
     d=5 bounds for eta=1e-4,1e-3,1e-2.
  4. If stable, feed the finite spectrum and regulator support into the full
     threshold/proton scan.
```

Heartbeat update 2026-05-08 07:36 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The triplet-regulator matrix-support
  alignment obstruction now has a minimal discrete solution: add a Z_2^G parity
  under which only the G120-type triplet source is odd.  This removes the six
  dangerous charged entries to all orders in S1,S2 spurion monomials, but the
  global compatibility of this parity with the full 120_H doublet/Yukawa sector
  still has to be audited.

new script:
  code/find_z3_z2_alignment_symmetry.py

new outputs:
  output/z3_z2_alignment/z3_z2_alignment_summary.json
  output/z3_z2_alignment/z3_z2_alignment_report.md
  output/z3_z2_alignment/z3_z2_alignment_scan.csv
  output/z3_z2_alignment/z3_z2_allowed_monomials.csv

derivation:
  Keep the Z_3 charges from the triplet regulator and add an extra Z_2^G:

                H10   F126   G120   Ktr
    p(T_A)       0      0      1      0
    p(barT_B)    0      0      1      0

  with

    p(S0)=p(S1)=p(S2)=0.

  Since every S1^a S2^b monomial is even under Z_2^G, any entry with exactly
  one G120 leg is forbidden to all orders in these spurions.  Combining this
  with the original Z_3 residue rule gives exactly:

    r=0: (H10,F126,Ktr)_left x (H10,F126)_right through S0,
    r=1: (H10,F126,Ktr)_left x Ktr_right through S2,
    r=2: G120_left x G120_right through S1.

  This is precisely the support of the finite kappa=100 regulator.

numerical/search audit:
  The extra-grading scan over small cyclic groups finds the minimal exact
  solution at Z_2:

    missed desired entries = 0,
    allowed dangerous entries at degree one = 0,
    dangerous entries generated by S1^a S2^b through degree 8 = 0.

  Proton replay is unchanged:

    max Knu amplitude = 8.098721968e-5,
    K0mu amplitude    = 6.626369280e-5,
    RRRR amplitude    = 1.066375671e-4,
    S_T max for tau>2.4e34 yr = 7.994484708e-3.

paper update:
  Added "Z3 x Z2^G triplet alignment" to paper/gut_framework.tex.

current obstacle:
  The triplet support is now symmetry-protected, but Z_2^G must be reconciled
  with the full Spin(10)/Pati-Salam Higgs sector.  If G120 is needed in the
  doublet Yukawa texture, a naive parity on the entire 120_H multiplet may
  remove desired flavor operators.  The parity may need to act only on the
  triplet mediator copy, or the model must split the 120-like triplet filter
  from the physical 120_H Yukawa channel.

next attempted nontrivial idea:
  Perform a compatibility audit of Z_2^G with the Yukawa/Majorana sector:
  list which operators 16_i16_j10_H, 16_i16_j\bar126_H, 16_i16_j120_H,
  triplet-filter links, and Higgs doublet mixings remain allowed under possible
  charge assignments.  Search for the minimal assignment that preserves the
  previously fitted O(2) Yukawa texture while keeping the triplet regulator
  protected.

verification plan:
  1. Build an operator ledger with Z_3 x Z_2^G charges for Yukawa, Majorana,
     triplet-filter, and doublet-mixing terms.
  2. Scan whether a separate mediator-only G120_tr copy is required.
  3. If a consistent charge table exists, replay the Yukawa singular-value
     texture and d=5 proton amplitudes with the protected regulator.
  4. Feed the resulting finite triplet spectrum into the threshold/proton scan.
```

Heartbeat update 2026-05-08 07:46 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The new G120-parity is compatible
  with flavor only if it acts on a mediator-only G_tr copy or is interpreted as
  a post-Spin(10) triplet-sector parity.  It cannot act on the full physical
  120_H multiplet while retaining the Clebsch flavor branch.

new script:
  code/audit_gparity_yukawa_compatibility.py

new outputs:
  output/gparity_yukawa_compatibility/gparity_yukawa_compatibility_summary.json
  output/gparity_yukawa_compatibility/gparity_yukawa_compatibility_report.md
  output/gparity_yukawa_compatibility/gparity_family_parity_scan.csv
  output/gparity_yukawa_compatibility/gparity_operator_ledger.csv
  output/gparity_yukawa_compatibility/gparity_flavor_branch_metrics.csv

mathematical no-go:
  Suppose q(120_H)=1 under the new Z_2^G and all three antisymmetric family
  entries of 16_i 16_j 120_H are required.  Then

    q_i + q_j = 1 mod 2

  for pairs 12,13,23.  Adding these equations gives

    2(q_1+q_2+q_3) = 3 = 1 mod 2,

  which is impossible.  Therefore no Z_2 family-parity assignment can make a
  full physical 120_H odd while retaining the full antisymmetric 120 Yukawa
  tensor.

scan audit:
  The script scans all 128 assignments of

    q(16_i), q(10_H), q(126bar_H), q(120_H), q(16bar_H)

  and checks whether all symmetric 10/126 entries, all three antisymmetric
  120 entries, and all Majorana entries are allowed.

  Results:

    full-flavor assignments with q(120_H)=1: 0,
    full-flavor assignments with q(120_H)=0: 4.

flavor numerical context:
  Existing local flavor audits give:

    no-120 exact O(-4) card:
      |Vus|=0.670251, |Vcb|=0.517568, |Vub|=0.372697,
      CKM score=5.449970, mass score=1.237957e-2.

    single-120 Clebsch:
      |Vus|=0.224273, |Vcb|=0.040579, |Vub|=0.003728,
      CKM score=3.267497e-5, mass score=1.312006.

    two-120 Clebsch:
      |Vus|=0.216103, |Vcb|=0.046162, |Vub|=0.004059,
      CKM score=4.580687e-3, mass score=0.644013.

    two-120 plus S_perp symmetric relaxation:
      |Vus|=0.178719, |Vcb|=0.044810, |Vub|=0.004195,
      CKM score=1.446573e-2, mass score=0.174824,
      phenomenology_viable=True in the stored local criterion.

operator conclusion:
  The viable branch keeps physical 120_H even and introduces a distinct
  mediator-only G_tr copy that is odd under Z_2^G.  Then 16_i16_j120_H remains
  allowed, the triplet regulator is protected, and physical 120_H--G_tr mixing
  is forbidden unless an explicit odd bridge spurion is introduced.

paper update:
  Added "Z2^G compatibility with 120_H Yukawa" to paper/gut_framework.tex.

current obstacle:
  The mediator-only G_tr branch is now the clean symmetry option, but it adds
  extra triplet-sector degrees of freedom whose complete multiplet status and
  threshold contribution must be audited.  If G_tr is not a complete heavy
  Spin(10) multiplet, it may introduce new non-universal thresholds; if it is
  a complete 120-like multiplet, UV perturbativity and doublet partners must be
  revisited.

next attempted nontrivial idea:
  Build a mediator-only G_tr spectrum card: list whether G_tr is a full
  120'_H-like multiplet, a Pati-Salam triplet-only EFT field, or a paired
  vectorlike remnant.  For each option compute threshold vectors, allowed
  mixing terms, and whether the Z_2^G parity remains exact after Spin(10)
  breaking.

verification plan:
  1. Enumerate the SU(5)/Pati-Salam fragments of the mediator-only G_tr option.
  2. Compute its one-loop beta-vector contribution and decide whether it is
     universal or a new threshold parameter.
  3. Replay the corrected RGE/proton scan with this finite regulator spectrum.
  4. If the branch remains safe, update the benchmark card and d=5 proton
     scan to use the mediator-only protected triplet matrix.
```

Heartbeat update 2026-05-08 07:56 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The G-parity branch has now been
  refined: physical 120_H must remain even, while the protected G_tr source is
  mediator-only.  The remaining issue is what representation content this
  mediator-only G_tr actually has.

new script:
  code/audit_gtr_mediator_spectrum.py

new outputs:
  output/gtr_mediator_spectrum/gtr_mediator_spectrum_summary.json
  output/gtr_mediator_spectrum/gtr_mediator_spectrum_report.md
  output/gtr_mediator_spectrum/gtr_threshold_replay.csv
  output/gtr_mediator_spectrum/gtr_lock_windows.csv

mathematical audit:
  For a chiral multiplet (R3,R2)_Y,

    b1 = (3/5) Y^2 dim(R3) dim(R2),
    b2 = T2(R2) dim(R3),
    b3 = T3(R3) dim(R2).

  Therefore a vectorlike triplet pair

    D + Dbar = (3,1,-1/3) + (bar3,1,+1/3)

  has

    b_D = (2/5, 0, 1),
    ||P b_D||_2 = 0.7118052168.

  The inert doublet pair has b_L=(3/5,1,0), so a common 5+5bar pair gives

    b_5 = b_D + b_L = (1,1,1),

  and is exactly projected-threshold silent.  A full 120'_H gives

    b_120 = (28,28,28),

  also projected-threshold silent if degenerate, but it adds T(120)=28 to the
  Spin(10) UV beta function.

numerical threshold windows:
  A triplet-only remnant with mass kappa M_G contributes

    Delta_i = b_i log(1/kappa)/(2 pi),
    ||P Delta_D||_2 = 0.7118052168 |log kappa|/(2 pi).

  Thus the mass-lock windows are:

    ||P Delta|| < 1e-3:
      0.991212 < M_D/M_G < 1.008866.

    ||P Delta|| < 1e-2:
      0.915513 < M_D/M_G < 1.092284.

    ||P Delta|| < 5.4e-2:
      0.620851 < M_D/M_G < 1.610692.

RGE/proton replay:
  The corrected-cache replay does not immediately kill a triplet-only remnant,
  but it makes G_tr a new non-universal threshold parameter.  For the
  triplet-only branch:

    M_D/M_G = 0.1:
      ||P Delta||=0.260854, safe points=114,
      best M_Sigma3=6.009e14 GeV.

    M_D/M_G = 1:
      ||P Delta||=0, safe points=146,
      best M_Sigma3=9.204e14 GeV.

    M_D/M_G = 10:
      ||P Delta||=0.260854, safe points=184,
      best M_Sigma3=1.410e15 GeV.

new nontrivial idea:
  Replace the default "full 120'_H mediator copy" by a sterile Spin(10)
  10'_G completion.  Its 5+5bar contains the needed color triplet pair plus an
  inert doublet partner, so P Delta_Gtr=0 when degenerate.  It raises the UV
  Dynkin index by only T(10)=1 rather than T(120)=28.  The label G_tr then
  denotes the protected triplet source in the regulator matrix, not the
  physical antisymmetric 120_H Yukawa multiplet.

paper update:
  Added "Mediator-only G_tr spectrum audit" to paper/gut_framework.tex.

current obstacle:
  The sterile 10'_G option is now the preferred action-level completion, but
  it requires an explicit projector/selection rule showing that only its
  triplet component enters the d=5 regulator while the doublet partner remains
  inert and degenerate.  Without that rule, the 10'_G doublet may mix with MSSM
  Higgs doublets or create a new mu/doublet-threshold problem.

next attempted nontrivial idea:
  Build a 10'_G missing-partner/projector superpotential.  Use the existing
  F_54 and/or PS adjoint direction to couple the triplet component into the
  regulator matrix while forbidding or lifting doublet mixing.  The target is
  a degenerate 5+5bar threshold for gauge matching plus a triplet-only
  effective d=5 coupling.

verification plan:
  1. Write the 10'_G component mass matrix for D_G, Dbar_G, L_G, Lbar_G.
  2. Check which F_54/PS projectors split D and L and quantify the induced
     non-universal threshold if the split is nonzero.
  3. Add Higgs-doublet mixing operators to the charge ledger and prove they
     are forbidden, or include them in a mu-sector scan.
  4. Replay RGE/proton thresholds with the final 10'_G card.
```

Heartbeat update 2026-05-08 08:04 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The sterile 10'_G branch now has an
  explicit triplet/doublet projector audit and a candidate threshold-safe
  bridge mechanism.  The new blocker is to promote the determinant-locking
  relation to a genuine F-flat superpotential rather than an imposed algebraic
  condition.

new script:
  code/audit_10g_projector_superpotential.py

new outputs:
  output/ten_g_projector/ten_g_projector_summary.json
  output/ten_g_projector/ten_g_projector_report.md
  output/ten_g_projector/ten_g_mass_split_scan.csv
  output/ten_g_projector/ten_g_bridge_scan.csv
  output/ten_g_projector/ten_g_threshold_tolerances.csv

mathematical projector:
  In the Spin(10) vector 10,

    S_54 = diag(-2,-2,-2,-2,-2,-2,3,3,3,3).

  Therefore

    P_D = (3I - S_54)/5,
    P_L = (S_54 + 2I)/5

  are exact projectors:

    ||P_D^2-P_D|| = 0,
    ||P_L^2-P_L|| = 0,
    ||P_D P_L|| = 0,
    rank(P_D)=6,
    rank(P_L)=4.

  Thus singlet+54 insertions can couple only the (6,1,1) color-triplet
  component of 10'_G while leaving the (1,2,2) doublet inert.

mass-splitting warning:
  A generic 54 mass term

    W contains M X Xbar + xi M X S_54 Xbar

  gives

    M_D/M = 1 - 2 xi,
    M_L/M = 1 + 3 xi.

  The threshold scan gives

    |xi| < 1.763854e-3 for ||P Delta|| < 1e-3,
    |xi| < 1.748850e-2 for ||P Delta|| < 1e-2,
    |xi| < 8.938605e-2 for ||P Delta|| < 5.4e-2.

  Therefore the 10'_G mass must be protected as a universal mass; the 54_H
  insertion should be used as a projector in bridge couplings, not as a free
  mass-splitting coefficient.

bridge audit:
  For the physical triplet source A and sterile source G, write

    M_D = [[1, epsilon_R],
           [epsilon_L, m_G]].

  Three cases were checked:

    1. triangular_one_sided:
       epsilon_L=0, m_G=1.
       det M_D=1 and P Delta=0, but (M_D^{-1})_AA=1, so it does not repair
       the physical 120 inverse propagator by itself.

    2. two_sided_equal:
       epsilon_L=epsilon_R=epsilon, m_G=1.
       At epsilon=0.3:
         det M_D=0.91,
         (M_D^{-1})_AA=1.098901,
         ||P Delta||=1.068420e-2,
         safe points=147.

    3. determinant_locked_equal:
       epsilon_L=epsilon_R=epsilon, m_G=1+epsilon^2.
       At epsilon=0.3:
         det M_D=1,
         (M_D^{-1})_AA=1.09,
         ||P Delta||=0,
         safe points=146.

new nontrivial idea:
  Use determinant locking as the 10'_G bridge principle:

    m_G = 1 + epsilon_L epsilon_R.

  Then the bridge modifies the physical inverse propagator,

    (M_D^{-1})_AA = 1 + epsilon_L epsilon_R,

  and generates off-diagonal inverse entries, but keeps the holomorphic
  triplet determinant equal to the complete-pair value.  This protects the
  one-loop gauge threshold while allowing the d=5 Wilson tensor to move.

paper update:
  Added "10'_G triplet projector and determinant locking" to
  paper/gut_framework.tex.

current obstacle:
  Determinant locking is currently an algebraic relation.  To make the model
  action-level, introduce a driving field X_lock or equivalent constrained
  sector enforcing

    m_G - 1 - epsilon_L epsilon_R = 0

  through F-flatness.  The same charge table must forbid 10'_G doublet mixing
  with MSSM Higgs doublets.

next attempted nontrivial idea:
  Construct W_lock + W_bridge explicitly with singlet bridge spurions E_L,E_R
  and mass spurion M_G^D, assign Z3 x Z2^G x R charges, and check:

    F_Xlock = 0 gives determinant locking,
    projected triplet bridge terms are allowed,
    doublet bridge and MSSM mu-mixing terms are forbidden,
    no new non-universal threshold appears.

verification plan:
  1. Enumerate charge assignments for E_L,E_R,M_G^D,X_lock.
  2. Prove the allowed terms generate the determinant-locked 2x2 triplet
     matrix and leave the doublet mass universal.
  3. Scan leakage terms through low operator degree.
  4. Replay d=5 Wilson tensors using the determinant-locked inverse block.
```

Heartbeat update 2026-05-08 08:16 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The image hint has been converted
  into a precise Schur-complement/path-sum determinant-locking theorem.  This
  gives a nontrivial route to move dimension-five Wilson tensors without
  moving one-loop gauge thresholds.

new script:
  code/audit_pathsum_determinant_locking.py

new outputs:
  output/pathsum_determinant_locking/pathsum_determinant_locking_summary.json
  output/pathsum_determinant_locking/pathsum_determinant_locking_report.md
  output/pathsum_determinant_locking/pathsum_determinant_locking_scan.csv

mathematical theorem:
  For a block triplet mass matrix

    M = [[A, U],
         [V, G0 + V A^{-1} U]],

  the Schur complement of A is exactly G0, hence

    det M = det A det G0.

  But the open physical propagator block is shifted:

    (M^{-1})_AA = A^{-1} + A^{-1} U G0^{-1} V A^{-1}.

  Therefore closed-path gauge-threshold data can remain fixed while open-path
  proton-decay Wilson data move.

numerical audit:
  The script tests 1+1, 2+1, and 3+2 block systems for epsilon in
  {0.05,0.10,0.20,0.30,0.50}.

  Across the grid:

    max ||P Delta_locked||_2 = 2.515e-17,
    max inverse-identity error = 2.483e-16,
    max absolute Wilson-proxy relative shift = 0.25.

  In the 3+2 benchmark at epsilon=0.30:

    locked determinant ratio = 1.0000000000000002,
    naive determinant ratio = 0.9450877330,
    locked ||P Delta|| = 2.515e-17,
    naive ||P Delta|| = 6.398e-3,
    locked Wilson-proxy shift = 7.252e-2.

new nontrivial idea:
  Reinterpret the 10'_G determinant lock as a path-sum condition rather than
  an isolated 2x2 trick.  The diagrammatic content is:

    closed loops determine log det M and gauge thresholds,
    open paths determine M^{-1} and proton decay.

  The goal is a determinant-silent open-path deformation of the triplet
  inverse propagator.

paper update:
  Added "Path-sum determinant locking" to paper/gut_framework.tex.

current obstacle:
  The path-sum condition is still an algebraic construction.  A real
  action-level model must promote

    G - G0 - V A^{-1} U = 0

  or its polynomial equivalent to F-flatness.  Since A^{-1} is non-polynomial,
  the actual superpotential must introduce auxiliary mediator fields whose
  F-equations generate the same Schur complement without writing inverse
  matrices by hand.

next attempted nontrivial idea:
  Build a polynomial Schur-locking superpotential with auxiliary fields:

    W_aux = X_A(A Abar - I) + X_G(G - G0 - V Abar U)

  or a finite mediator chain whose integrated-out equations implement
  Abar=A^{-1}.  Then audit whether all fields come in complete multiplets or
  determinant-locked pairs.

verification plan:
  1. Construct the smallest polynomial auxiliary realization of the Schur
     complement identity.
  2. Verify F-flatness gives det M = det A det G0 without non-polynomial
     assumptions.
  3. Feed the resulting shifted inverse block into the existing d=5 Wilson
     tensor scan.
  4. Check whether auxiliary fields add complete multiplet thresholds or
     require their own determinant locks.
```

Heartbeat update 2026-05-08 08:24 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete.  The Schur-complement/path-sum
  determinant lock is no longer just a non-polynomial algebraic condition:
  it now has a polynomial F-flat realization with auxiliary matrices.

new script:
  code/audit_polynomial_schur_locking.py

new outputs:
  output/polynomial_schur_locking/polynomial_schur_locking_summary.json
  output/polynomial_schur_locking/polynomial_schur_locking_report.md
  output/polynomial_schur_locking/polynomial_schur_locking_scan.csv

polynomial superpotential:
  Introduce auxiliary singlet matrices B,C and driver matrices X,Z,Y:

    W_lock = Tr X(A B - I)
           + Tr Z(C - B U)
           + Tr Y(G - G0 - V C).

  The driver F-terms are:

    F_X = A B - I = 0,
    F_Z = C - B U = 0,
    F_Y = G - G0 - V C = 0.

  Therefore, on the invertible-A branch:

    B = A^{-1},
    C = A^{-1} U,
    G = G0 + V A^{-1} U.

  Setting X=Y=Z=0 makes all non-driver F-terms vanish.  The construction uses
  only degree <= 3 interaction terms:

    X A B,  Z B U,  Y V C.

numerical audit:
  The script checks 1+1, 2+1, and 3+2 block examples at epsilon=0.10,0.30,0.50.

  Across the scan:

    max driver F residual = 9.440e-17,
    max non-driver F residual = 0,
    max inverse identity error = 2.483e-16,
    max |log(det M/(det A det G0))| = 2.220e-16.

  For the 3+2, epsilon=0.30 benchmark:

    |log(det M/(det A det G0))| = 2.220e-16,
    Wilson proxy shift = 7.252e-2,
    min singular value of full matrix = 0.8605.

paper update:
  Added "Polynomial Schur-locking superpotential" to paper/gut_framework.tex.

current obstacle:
  The algebraic/F-flat mechanism works.  The remaining blocker is
  representation and grading compatibility: assign Spin(10),
  Z3 x Z2^G, and R charges to A,U,V,G,B,C,X,Y,Z so that:

    1. W_lock is allowed,
    2. doublet bridge and MSSM-Higgs mixing terms are forbidden,
    3. the auxiliary fields are gauge singlets or complete multiplets,
    4. no incomplete threshold is introduced.

next attempted nontrivial idea:
  Build a combined charge ledger for the Schur-locking auxiliaries.  Treat
  A,G as triplet-source mass blocks, U,V as projected triplet bridges, B,C as
  singlet inverse/path auxiliaries, and X,Y,Z as R-charge-2 drivers.  Scan
  finite Z_N x Z2^G x U(1)_R assignments for allowed W_lock terms and forbidden
  dangerous doublet/mu terms.

verification plan:
  1. Enumerate the operator list for W_lock plus dangerous doublet/Higgs
     mixings.
  2. Search small Z_N charge assignments that allow the lock and forbid the
     dangerous terms.
  3. Confirm the allowed bridge terms reproduce the determinant-locked inverse
     block used in the path-sum audit.
  4. If no assignment exists, demote the path-sum lock to a post-Spin(10) EFT
     condition and mark the UV completion as open.
```

Heartbeat update 2026-05-08 08:34 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the polynomial Schur-locking
  auxiliary sector now has a verified finite charge ledger.  The path-sum
  determinant lock is therefore no longer merely an allowed algebraic/F-flat
  pattern; its representative dangerous lower-path and doublet-bridge
  operators can be forbidden by a small grading.

new script:
  code/audit_schur_locking_charge_ledger.py

new outputs:
  output/schur_locking_charge_ledger/charge_table.csv
  output/schur_locking_charge_ledger/term_ledger.csv
  output/schur_locking_charge_ledger/summary.json
  output/schur_locking_charge_ledger/report.md

finite charge ledger:
  The scan finds a minimal solution:

    Z_2^inv x Z_2^path x U(1)_R,    R(W)=2.

  Nonzero charges:

    Z_2^inv:
      A,B,G,G0,U_D,V_D,Y are odd.

    Z_2^path:
      C,U_D,V_D,U_L,V_L,Z are odd.

    U(1)_R:
      X,Y,Z have R=2; the test mass spurion M_mix has R=2;
      lock variables have R=0.

allowed terms:
    XAB, XI_A, ZC, ZBU_D, YG, YG0, YV_D C.

forbidden representative leakage terms:
    XA, XB, ZU_D, ZBU_L, YV_L C, YV_D, YC,
    M_mix H120 U_D, M_mix H120 V_D, M_mix H120 G,
    M_mix H_u U_L, M_mix H_d V_L,
    M_mix U_D U_L, M_mix V_D V_L.

numerical/logic audit:
  required_allowed = True
  forbidden_rejected = True

  The key nontrivial point is the inverse-block grading of A,B.  If B is
  neutral under all finite charges, the desired ZBU_D term is
  charge-indistinguishable from the lower-path term ZU_D.  The minimal
  Z_2^inv grading solves exactly this problem.

paper update:
  Added "Schur-locking charge ledger" to paper/gut_framework.tex.

current obstacle:
  The charge sector is now consistent at the abstract block level.  The next
  hard problem is representation completeness: the 54/210 source, driver, and
  Schur-locking fields must be embedded into explicit Spin(10) multiplets and
  audited component by component.  If non-singlet auxiliary modes are not
  complete multiplets, determinant-locked pairs, or sufficiently heavy, they
  generate new non-universal thresholds and must be included in the
  RGE/proton scan.

next attempted nontrivial idea:
  Build a component spectrum audit for the Schur-locking source/driver sector.
  Treat the charge ledger as fixed, then enumerate the Spin(10) representation
  content of A,B,C,G,U_D,V_D,U_L,V_L,X,Y,Z.  For each PS/SM fragment, classify
  it as:

    1. gauge singlet auxiliary,
    2. complete lifted Spin(10) or SU(5) multiplet,
    3. determinant-locked pair with zero projected threshold,
    4. incomplete threshold requiring RGE/proton replay.

verification plan:
  1. Extend the component Hessian card with the Schur-locking auxiliary fields.
  2. Compute fragment multiplicities and beta-vector contributions.
  3. Verify whether the charge-protected lock adds zero projected threshold.
  4. If not, feed the finite auxiliary threshold vector into the existing
     two-loop RGE and dimension-five proton-decay scans.
```

Heartbeat update 2026-05-08 08:49 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the Schur-locking auxiliary
  component-spectrum audit is now done.  The lock itself no longer appears to
  be the threshold problem.

new script:
  code/audit_schur_auxiliary_component_spectrum.py

new outputs:
  output/schur_auxiliary_component_spectrum/component_classification.csv
  output/schur_auxiliary_component_spectrum/branch_threshold_scan.csv
  output/schur_auxiliary_component_spectrum/clean_windows.csv
  output/schur_auxiliary_component_spectrum/uv_perturbativity_rows.csv
  output/schur_auxiliary_component_spectrum/summary.json
  output/schur_auxiliary_component_spectrum/report.md

minimal Schur completion:
  Treat B,C,X,Y,Z,I_A as SM/Spin(10) singlet auxiliary matrices.
  Embed G,G0,U_D,V_D together with inert U_L,V_L in a degenerate sterile
  10'_G pair.

  Then:

    b_singlet = (0,0,0),
    b_D+Dbar = (2/5,0,1),
    b_L+Lbar = (3/5,1,0),
    b_5+5bar = (1,1,1),
    P b_5+5bar = 0.

  Therefore the minimal Schur-locking auxiliary sector adds zero projected
  one-loop threshold.

numerical audit:
  minimal_projected_threshold_l2 = 0.000e+00
  new SO(10) Dynkin index increment = T(10) = 1
  baseline_plus_10G Landau ratio Lambda_LP/MG = 2.52e+02

fallback triplet-only branch:
  If only the triplet bridge survives as a post-Spin(10) remnant,

    ||P b_D+Dbar||_2 = 0.7118052168.

  Clean-threshold windows:

    ||P Delta|| < 1e-2:
      0.915513 < M_D/M_G < 1.092284

    ||P Delta|| < 1e-3:
      0.991212 < M_D/M_G < 1.008866

UV obstruction sharpened:
  The Schur lock itself can be threshold-safe and UV-mild.  The large
  54/210 source/driver alignment sector is still the hard problem:

    documented complete drive pairs:
      sum T = 184, b10 = 160, Lambda_LP/MG = 5.07

    W_rel source blocks:
      sum T = 300, b10 = 276, Lambda_LP/MG = 2.56

  These do not reach R=50 or R=200 if treated as propagating chiral fields
  at M_G.

paper update:
  Added "Schur auxiliary component spectrum" to paper/gut_framework.tex.

current obstacle:
  The determinant/path Schur lock can now be made charge-consistent,
  F-flat, projected-threshold silent, and UV-mild on the minimal singlet plus
  sterile-10'_G branch.  The remaining obstruction is the large 54/210
  alignment/source sector: it must be interpreted as constrained/composite, or
  replaced by a smaller UV locking mechanism, otherwise the SO(10) coupling
  hits a Landau pole only a few times above M_G.

next attempted nontrivial idea:
  Construct a reduced UV locking mechanism that preserves the verified
  projector data but avoids propagating full 54/210 driver towers.  The most
  concrete route is a constrained/composite F54/210 source sector with only
  orbit-coordinate singlets propagating, while the non-singlet directions are
  nondynamical constraints or confined composites above M_G.

verification plan:
  1. Write a reduced constrained-source effective action with F54 and 210 as
     fixed orbit coordinates plus singlet radial modes.
  2. Compute its local Hessian: only singlet modes may remain propagating.
  3. Recompute SO(10) b10 and require Lambda_LP/MG > 200.
  4. Confirm that the low-energy Clebsch/projector spectrum and the
     Schur-locking charge ledger are unchanged.
```

Heartbeat update 2026-05-08 09:02 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the reduced constrained/composite
  54/210 source-sector branch is now explicitly audited.  This gives a
  concrete way to evade the large 54/210 Landau-pole obstruction, at the price
  of interpreting the projector sources as constrained/composite order
  parameters rather than elementary propagating Higgs towers.

new script:
  code/audit_constrained_54_210_source_sector.py

new outputs:
  output/constrained_54_210_source_sector/summary.json
  output/constrained_54_210_source_sector/report.md
  output/constrained_54_210_source_sector/projector_values.csv
  output/constrained_54_210_source_sector/uv_rows.csv

mathematical audit:
  54 source:
    S0 = diag(-2^6,3^4), S0^2 - S0 - 6I = 0.
    Linearized spectral-constraint rank = 30 on the 54.
    Fixed orbit nullity = 24 = dim SO(10)/(SO(6)xSO(4)).
    If one radial singlet is retained, removed normal modes = 54 - 25 = 29.

  210 source:
    Omega0 = e7 wedge e8 wedge e9 wedge e10.
    SO(10)-orbit rank = 24, same stabilizer SO(6)xSO(4).
    If one radial singlet is retained, removed normal modes = 210 - 25 = 185.
    D210 eigenvalue counts on two-forms = (-1)^3, 0^39, (+1)^3.

projector preservation:
  max projector error = 1.110e-16.

UV audit:
  constrained_54_210_sources_plus_10G:
    sum T = 71, b10 = 47, Lambda_LP/MG = 2.52e2,
    reaches R=50 and R=200.

  elementary_54_210_plus_10G:
    sum T = 139, b10 = 115, Lambda_LP/MG = 9.58,
    fails R=50 and R=200.

paper update:
  Added "Constrained 54/210 source sector" to paper/gut_framework.tex.

current obstacle:
  The model now has a mathematically consistent conditional UV-mild branch:
  constrained/composite F54 and Omega210 sources plus the minimal Schur
  sterile-10G completion.  What remains incomplete is deriving such a
  constrained/composite source sector from a microscopic first-principles
  dynamics, rather than imposing orbit constraints as EFT data.

next attempted nontrivial idea:
  Build a microscopic origin test for the constrained-source branch.  The
  cleanest local target is a Lagrange-multiplier/conormal action whose
  F-terms impose:

    S^2 - S - 6I = 0,
    Omega decomposable and normalized,
    shared SO(10)/(SO(6)xSO(4)) orientation between S and Omega.

  Then audit whether the Lagrange multipliers themselves can be kept
  auxiliary/constrained without reintroducing propagating 54/210 normal modes.

verification plan:
  1. Write explicit conormal F-term constraints for S and Omega.
  2. Compute the combined Jacobian rank and relative-orientation zero modes.
  3. Confirm only one shared orbit plus radial singlets survive.
  4. If multiplier normal modes propagate, reuse the Xi_N threshold fallback
     scan; otherwise keep the constrained branch as the preferred UV-mild EFT.
```

Heartbeat update 2026-05-08 09:12 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the aligned conormal
  Lagrange-multiplier action for the constrained 54/210 source sector is now
  explicitly tested.  The relative orientation problem is solved at the
  linearized F-term level.

new script:
  code/audit_combined_conormal_54_210.py

new outputs:
  output/combined_conormal_54_210/summary.json
  output/combined_conormal_54_210/report.md
  output/combined_conormal_54_210/rank_summary.csv

conormal F-term system:
  Variables:
    S in 54,
    Omega in Lambda^4 R^10.

  Constraints:
    C54(S) = S^2 - S - 6 I = 0,
    A(S,Omega) = (K_S - 12 I) Omega = 0,
    N(Omega) = <Omega,Omega> - 1 = 0.

  Here K_S is the induced action of S on four-forms.  At
  S0=diag(-2^6,3^4), the eigenvalue-12 subspace of K_S is the single
  weak-volume line e7^e8^e9^e10, so A(S,Omega)=0 aligns Omega to the same
  weak four-plane selected by S.

rank audit:
  variables dimension = 54 + 210 = 264
  unprojected constraint rows = 55 + 210 + 1 = 266
  Jacobian rank/nullity = 240 / 24
  shared orbit tangent rank = 24
  relative orientation modes removed = 24

multiplier Hessian audit:
  Full unprojected multipliers:
    zero modes = 264 + 266 - 2*240 = 50.
    This includes 26 extra multiplier/radial modes and is rejected as a
    propagating tower.

  Normal-bundle multipliers:
    multiplier dimension = 240,
    zero modes = 264 + 240 - 2*240 = 24.
    These are exactly the shared SO(10)/(SO(6)xSO(4)) orbit modes.

  If two radial source singlets are deliberately retained:
    physical zero modes = 26.

210 Clebsch check:
  D210 eigenvalue counts on two-forms = (-1)^3, 0^39, (+1)^3.

paper update:
  Added "Aligned conormal action for 54/210 sources" to paper/gut_framework.tex.

current obstacle:
  The constrained 54/210 source sector is now mathematically coherent as a
  conormal EFT: it preserves the projector data, removes relative orientation
  moduli, and avoids the large elementary tower.  The remaining gap is
  microscopic: explain why the conormal multipliers are auxiliary/composite
  normal-bundle objects rather than propagating unprojected chiral fields.

next attempted nontrivial idea:
  Audit the fallback in which the normal-bundle multipliers propagate with a
  finite mass.  Use the existing Xi_N threshold logic, but now with the
  combined 54/210 normal-bundle rank.  If even a modest propagating multiplier
  mass creates unacceptable non-universal threshold or UV cost, then the paper
  must state the constrained branch as an EFT/composite assumption rather than
  a perturbative elementary completion.

verification plan:
  1. Construct the normal-bundle multiplier threshold card for the combined
     rank-240 conormal sector.
  2. Bound its allowed mass window from ||P Delta|| < 1e-2 and 1e-3.
  3. Reuse the corrected RGE/proton scan for representative kappa values.
  4. Decide whether "auxiliary/composite only" should be promoted to an
     explicit assumption in the theorem map.
```

Heartbeat update 2026-05-08 09:22 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the propagating conormal
  multiplier fallback has now been quantified.  The result is strong enough
  to promote "auxiliary/composite normal-bundle multipliers" to an explicit
  preferred-branch assumption.

new script:
  code/scan_combined_conormal_multiplier_threshold.py

new outputs:
  output/combined_conormal_multiplier_threshold/summary.json
  output/combined_conormal_multiplier_threshold/report.md
  output/combined_conormal_multiplier_threshold/combined_conormal_multiplier_threshold_scan.csv

beta-vector derivation:
  The shared orbit removed from both 54 and 210 is (6,2,2), with

    b_622 = (26/5, 6, 4).

  Therefore:

    b_54_normal  = (12,12,12) - b_622
                 = (34/5, 6, 8),

    b_210_normal = (56,56,56) - b_622
                 = (254/5, 50, 52).

  A propagating normal field plus conormal multiplier gives:

    b_Xi,combined = 2 (b_54_normal + b_210_normal)
                  = (576/5, 112, 120).

threshold window:
  coefficient ||P b_Xi/(2 pi)|| = 0.906298550.

  For ||P Delta_Xi|| < 1e-2:
    0.989027 < kappa_Xi < 1.011095.

  For ||P Delta_Xi|| < 1e-3:
    0.998897 < kappa_Xi < 1.001104.

RGE/proton replay:
  At R=200, kappa_Xi=1.01:

    ||P Delta_Xi|| = 9.017970e-03,
    safe points = 150,
    alphaG^-1 = 41.545142,
    M_Sigma3 = 9.348767e14 GeV,
    tau_d6 = 4.852144e35 yr.

paper update:
  Added "Propagating conormal multiplier fallback" to paper/gut_framework.tex.

current obstacle:
  The conormal multiplier fallback is possible only as a tuned threshold
  parameter near M_G.  The clean branch must assume the conormal multipliers
  are auxiliary/composite normal-bundle objects rather than propagating
  incomplete chiral fields.

next attempted nontrivial idea:
  Update the theorem/assumption map of the paper so the conditional branch is
  honest and non-circular:

    1. PSLT/Another Physics does not unconditionally derive a GUT.
    2. The viable branch is a conditional Spin(10) EFT.
    3. Its UV-mild projector sector assumes constrained/composite F54 and
       Omega210 sources with auxiliary/composite conormal multipliers.
    4. If conormal multipliers propagate, the Xi threshold must be scanned and
       mass-locked at the percent/per-mille level.

verification plan:
  1. Patch the introductory theorem/assumption map in TeX.
  2. Add a compact table of preferred branch vs fallback branch.
  3. Ensure no theorem claims elementary renormalizable 54_H+210_H completion.
  4. Recompile and check references/warnings.
```

Heartbeat update 2026-05-08 09:30 Asia/Taipei:

```text
completion status:
  The full GUT task remains incomplete, but the paper's front theorem map is
  now synchronized with the latest technical audits.  The preferred branch is
  explicitly identified as a constrained-source Spin(10) EFT, not an
  unconditional elementary renormalizable 54_H+210_H completion.

paper update:
  Added "Verified conditional branch map" near the start of
  paper/gut_framework.tex.

preferred branch assumptions now explicit:
  A1 Spin(10) EFT with complete 16_i matter multiplets.
  A2 H^0(CP^1,O(2)) protected three-family space.
  A3 Type-I seesaw with trace-lifted Majorana channel.
  A4 Proton bounds and two-loop threshold matching.
  A5 F54 and Omega210 are constrained/composite source order parameters,
     not elementary propagating 54_H+210_H towers.
  A6 Conormal multipliers are auxiliary/composite normal-bundle fields.
  A7 Schur lock uses singlet auxiliaries plus a degenerate sterile 10'_G.

fallback windows listed in the main assumption map:
  propagating conormal multiplier:
    0.989027 < kappa_Xi < 1.011095 for ||P Delta|| < 1e-2.

  triplet-only Schur bridge:
    0.915513 < kappa_D < 1.092284 for ||P Delta|| < 1e-2.

current obstacle:
  The mathematical EFT branch is now internally coherent, but not a completed
  first-principles UV derivation.  The remaining scientific gap is a true
  microscopic origin for the constrained/composite 54/210 source sector and
  auxiliary conormal multipliers.

next attempted nontrivial idea:
  Shift from adding more threshold machinery to writing the final "conditional
  theorem vs no-go theorem" structure:

    1. No-go: PSLT/Another Physics alone do not determine a unique GUT.
    2. Conditional theorem: given A1-A7, the branch is anomaly-safe,
       hypercharge-correct, threshold/proton-safe in the verified scans, and
       UV-mild up to R=200.
    3. Open UV problem: derive A5-A6 from a microscopic composite dynamics.

verification plan:
  1. Add a concise final theorem/summary section.
  2. Ensure every strong claim points to a script/output table.
  3. Mark A5-A6 as assumptions, not derived facts.
  4. Recompile and run warning scan.
```

Heartbeat update 2026-05-08 09:45 Asia/Taipei:

```text
completion status:
  The full first-principles GUT task remains incomplete.  However, the
  conditional theorem/no-go closure step is now complete and synchronized
  between code outputs and the TeX draft.

new verification ledger:
  Added code/audit_conditional_theorem_ledger.py.

  Generated:
    output/conditional_theorem_ledger/claim_ledger.csv
    output/conditional_theorem_ledger/summary.json
    output/conditional_theorem_ledger/report.md

ledger counts:
  FAIL: 1
  NO_GO: 1
  OPEN: 3
  PASS: 3
  PASS_AS_KINEMATICS: 1
  PASS_CONDITIONAL: 4
  TUNED_FALLBACK: 2

main numerical closure benchmark:
  R = 200, kappa_Xi = 1:
    safe points = 146
    M_Sigma3 = 9.211965e14 GeV
    M_Sigma8 = 3.887852e15 GeV
    alphaG^-1 = 41.363082
    tau_d6 = 4.809711e35 yr
    ||P Delta||_2 = 5.022739e-4

paper update:
  Replaced the old "Current Verdict" section in paper/gut_framework.tex with:
    1. No-go theorem for unconditional PSLT/Another Physics -> unique GUT.
    2. Verified conditional Spin(10) branch theorem under A1-A7.
    3. Corollary listing open blockers and tuned fallback windows.

verification:
  python script ran successfully.
  py_compile passed using a sandbox-local PYTHONPYCACHEPREFIX.
  pdflatex compiled twice.
  warning scan found no LaTeX Warning, Error, Undefined, or Overfull entries.

current obstacle:
  The branch is now paper-honest as a conditional EFT, but not a completed
  first-principles UV derivation.  The only remaining blockers are:
    O1 full joint CKM/PMNS/flavor fit,
    O2 channel-specific dressed d=5 proton decay,
    O3 microscopic origin for A5-A6 constrained/composite source and conormal
       multiplier assumptions.

next attempted nontrivial idea:
  Start O3 by testing a microscopic composite constrained-source model:
    - construct a confining or nonlinear-sigma source sector whose moduli are
      exactly the SO(10)/(SO(6)xSO(4)) orbit plus radial singlets;
    - prove that normal-bundle multipliers are auxiliary/composite rather than
      propagating incomplete fields;
    - check anomaly matching, holomorphic constraints, and beta-function cost.

verification plan:
  1. Build a small symbolic model for candidate composite-source constraints.
  2. Compute moduli dimension and stabilizer/orbit tangent rank.
  3. Audit whether any incomplete charged normal modes propagate.
  4. Feed any unavoidable modes into the conditional theorem ledger and
     RGE/proton threshold scan.
```

Heartbeat update 2026-05-08 00:49 UTC / 08:49 Asia/Taipei:

```text
completion status:
  The full first-principles GUT task remains incomplete.  This heartbeat
  converted the user's layer/spiral/dense-volume sketches into a grounded
  A5-A6 research-route audit, not a claimed proof.

new audit:
  Added code/audit_boundary_time_source_hint.py.

  Generated:
    output/boundary_time_source_hint/summary.json
    output/boundary_time_source_hint/report.md
    output/boundary_time_source_hint/rindler_lapse_scan.csv
    output/boundary_time_source_hint/rotating_boundary_stability.csv

mathematical checks:
  Rindler layer:
    N(n)=tanh(n/ell)=n/ell-(n/ell)^3/3+O(n^5).
    Max relative error for n/ell <= 1e-2 is 3.333200e-5.

  constrained source orbits:
    54 source ambient dimension = 54.
    orbit rank = 24.
    radial singlets retained = 1.
    normal modes removed = 29.

    210 source ambient dimension = 210.
    orbit rank = 24.
    radial singlets retained = 1.
    normal modes removed = 185.

  rotating boundary:
    E_co = omega - m Omega.
    For Sigma3, omega/MG = 0.1300156678 at the corrected benchmark.
    The scan brackets the m=1 critical value between stable Omega/MG=0.13
    and unstable Omega/MG=0.15.
    At Omega/MG=0.12, X622 and conormal m=1 margins are both 0.88, so the
    toy boundary rotation can stress Sigma3 without releasing the heavy modes.

paper update:
  Added "Boundary-Time Hint for the A5-A6 Source Sector" to
  paper/gut_framework.tex.  The proposition explicitly states this is a
  candidate mechanism, not a microscopic proof.

verification:
  py_compile passed with sandbox-local PYTHONPYCACHEPREFIX.
  The audit script ran successfully and all six pass flags are true.
  pdflatex compiled twice.
  warning scan found no LaTeX Warning, Error, Undefined, or Overfull entries.

current obstacle:
  A5-A6 are now better motivated by a boundary-time constrained-source ansatz,
  but still not derived from a microscopic action.  The remaining gap is to
  produce a real source-sector dynamics whose low-energy moduli are exactly:
    SO(10)/(SO(6)xSO(4)) orbit tangents + radial singlets,
  with no propagating incomplete normal-bundle modes.

next attempted nontrivial idea:
  Build a nonlinear-sigma/composite boundary source action:
    W_src = Lambda_a C_a(S,Omega) + W_radial + W_boundary,
  where the boundary lapse N(n) explains why normal deformations are
  nondynamical and the rotating interface provides a controlled near-critical
  Sigma3 scale.

verification plan:
  1. Write an explicit constrained-source action with kinetic metric and
     Lagrange/composite fields.
  2. Compute Hessian rank and check that only orbit tangents plus radial
     singlets remain.
  3. Check whether boundary rotation shifts Sigma3 without making X622 or
     conormal modes tachyonic/superradiant.
  4. If new incomplete fields propagate, add them to the threshold/proton scan;
     if not, update the conditional theorem ledger with this as a stronger
     A5-A6 candidate.
```

Heartbeat update 2026-05-08 00:53 UTC / 08:53 Asia/Taipei:

```text
completion status:
  The first-principles GUT task remains incomplete, but A5-A6 have been
  upgraded from a boundary-time hint to a local action-level constrained-source
  ansatz with a verified Hessian.  This is still not a UV microscopic proof.

new audit:
  Added code/audit_boundary_source_action.py.

  Generated:
    output/boundary_source_action/summary.json
    output/boundary_source_action/report.md
    output/boundary_source_action/local_hessian_spectrum.csv

local action:
  W_loc =
    Lambda54.N54^T s
    + Lambda210.N210^T omega
    + Pi.(T54^T s - T210^T omega)
    + mu54/2 (r54.s)^2
    + mu210/2 (r210.omega)^2.

Hessian verification:
  variables = 54 + 210 + 29 + 185 + 24 = 502.
  rank = 478.
  nullity = 24.
  max expected diagonal-orbit zero residual = 0.
  max null normal component = 1.034222e-15.
  max null radial component = 1.199024e-16.
  positive singular values are in [1.0, 1.5] in local units.

boundary/rotation checks:
  At n/ell=1e-2, N^2=tanh^2(1e-2)=9.999333e-5.
  At Omega/MG=0.12:
    Sigma3 margin = 1.001567e-2.
    X622 margin = 0.88.
    conormal margin = 0.88.

paper update:
  Added "Local Boundary Source Action" to paper/gut_framework.tex, including
  the explicit W_src^loc and the Hessian rank/nullity proof.

conditional theorem ledger update:
  Updated code/audit_conditional_theorem_ledger.py and regenerated
  output/conditional_theorem_ledger/*.

  New counts:
    FAIL: 1
    NO_GO: 1
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 5
    TUNED_FALLBACK: 2

  Added ledger row:
    "Local A5-A6 boundary source action leaves only shared orbit modes"
    with rank=478, nullity=24, Nnormal=1.03e-15.

verification:
  py_compile passed for both new/updated scripts with sandbox-local cache.
  audit_boundary_source_action.py ran successfully; all seven pass flags true.
  audit_conditional_theorem_ledger.py regenerated the ledger successfully.
  pdflatex compiled twice.
  warning scan found no LaTeX Warning, Error, Undefined, or Overfull entries.

current obstacle:
  The local A5-A6 action now removes normal modes and relative orientation modes
  correctly.  The remaining hard issue is UV origin: explain why
  Lambda54, Lambda210, and Pi are auxiliary/composite boundary fields rather
  than propagating incomplete charged multiplets.

next attempted nontrivial idea:
  Promote the local boundary action to a microscopic composite/NLSM model:
    - identify preon or sigma-model variables whose gauge-invariant composites
      are S and Omega;
    - derive the constraints C_a(S,Omega)=0 as F-terms or quantum-modified
      moduli constraints;
    - prove anomaly matching and absence of propagating incomplete normal
      multiplets.

verification plan:
  1. Write candidate composite field content and global/gauge charges.
  2. Compute moduli dimension and constraint rank.
  3. Check 't Hooft anomaly matching for any proposed global symmetry.
  4. Compute the induced SO(10) beta-function cost below and above the
     compositeness scale.
  5. If a propagating incomplete sector remains, feed its beta vector into the
     RGE/proton scan and mark the branch as tuned fallback.
```

Heartbeat update 2026-05-08 01:02 UTC / 09:02 Asia/Taipei:

```text
completion status:
  The full first-principles GUT task remains incomplete.  A5-A6 have now been
  strengthened from a local Hessian ansatz to a cutoff NLSM/composite source
  candidate.  This still does not prove the UV microscopic dynamics above the
  compositeness scale.

new audit:
  Added code/audit_nlsm_composite_source_origin.py.

  Generated:
    output/nlsm_composite_source_origin/summary.json
    output/nlsm_composite_source_origin/report.md
    output/nlsm_composite_source_origin/jacobian_singular_values.csv
    output/nlsm_composite_source_origin/anomaly_proxy.csv

candidate parameterization:
  U in SO(10)/(SO(6)xSO(4)).

  S(U,rho54) = rho54 U S0 U^T.
  Omega(U,rho210) = rho210 (U e7) wedge (U e8) wedge (U e9) wedge (U e10).

  The important selection is shared U for both S and Omega.

constraint checks:
  orthogonality error = 6.251720e-16.
  Tr S = 4.884981e-15.
  ||S^2 - rho54 S - 6 rho54^2 I|| = 1.800431e-14.
  ||Omega|| = 0.700000 for rho210=0.7.
  max Plucker residual = 2.602085e-18.

Jacobian/moduli result:
  ambient source dimension = 54 + 210 = 264.
  shared-U source rank = 26.
  independent-U source rank = 50.
  relative-orientation modes removed by shared U = 24.
  shared-U constraint rank = 238.

anomaly/beta proxy:
  coset tangent representation is (6,2,2), real under Pati-Salam.
  radial modes are singlets.
  new chiral anomaly proxy = 0.
  below the compositeness scale there is no elementary 54+210 tower, so the
  preferred EFT Dynkin-index increment from this source sector is 0.
  above the compositeness scale remains unknown and is the next UV problem.

paper update:
  Added "NLSM/Composite Source Candidate" to paper/gut_framework.tex with a
  proposition proving the shared-U rank and constraint count.

conditional theorem ledger update:
  Updated code/audit_conditional_theorem_ledger.py and regenerated
  output/conditional_theorem_ledger/*.

  New counts:
    FAIL: 1
    NO_GO: 1
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 6
    TUNED_FALLBACK: 2

  Added ledger row:
    "Shared-U NLSM source removes relative 54/210 orientation modes"
    with rank=26, constraints=238, Plucker=2.60e-18.

verification:
  py_compile passed for new/updated scripts with sandbox-local cache.
  audit_nlsm_composite_source_origin.py ran successfully; all eight pass flags
  true.
  audit_conditional_theorem_ledger.py regenerated successfully.
  pdflatex compiled twice.
  warning scan found no LaTeX Warning, Error, Undefined, or Overfull entries.

current obstacle:
  The source manifold is now locally and NLSM-consistently constructed, but
  the model still needs a UV microscopic origin for the shared-U sigma model.
  In particular, one must explain what dynamics produces the composite order
  parameters S and Omega, why only a shared U appears, and why no charged
  preon sector creates new anomalies or a Landau pole below the intended
  cutoff.

next attempted nontrivial idea:
  Test a minimal UV completion/no-go ledger for the NLSM source:
    A. purely elementary 54+210: already fails Landau reach;
    B. confining preon composites: check anomaly matching and Dynkin cost;
    C. deconstructed/link-field NLSM: check whether link fields form complete
       multiplets and whether shared-U locking is natural;
    D. accept cutoff NLSM as final conditional EFT and move to O1/O2 flavor
       and d=5 proton decay.

verification plan:
  1. Build a scorecard for A-D with explicit field content, moduli dimension,
     anomaly proxy, and beta-function cost.
  2. Reject candidates that reintroduce incomplete propagating thresholds.
  3. If no UV candidate passes, record a no-go: A5-A6 are conditional NLSM EFT
     assumptions, not first-principles consequences.
  4. If one passes, feed its spectrum into the RGE/proton ledger.
```

Heartbeat update 2026-05-08 01:07 UTC / 09:07 Asia/Taipei:

status:
  The full first-principles GUT task remains incomplete.  The A5-A6 source
  sector is now sharply classified: a shared-U NLSM/constrained-source branch
  passes as a conditional cutoff EFT, but the minimal microscopic routes do
  not close the UV origin problem at the displayed high mediator scales.

new audit:
  Added code/audit_nlsm_uv_completion_scorecard.py.

  Outputs:
    output/nlsm_uv_completion_scorecard/summary.json
    output/nlsm_uv_completion_scorecard/scorecard.csv
    output/nlsm_uv_completion_scorecard/landau_rows.csv
    output/nlsm_uv_completion_scorecard/report.md

  One-loop formula:
    b10 = sum_R T(R) - 24,
    alpha_G^{-1}(R M_G) = alpha_G^{-1}(M_G) - b10 log R/(2 pi).

  R=200 card results:
    A_elementary_54_210:
      sum_T = 139, b10 = 115, Lambda_LP/MG = 9.58,
      alpha_G^{-1}(200 MG) = -55.611, status = FAIL_HIGH_R.
    B_minimal_vector_preon_composite:
      sum_T = 75, b10 = 51, Lambda_LP/MG = 163,
      alpha_G^{-1}(200 MG) = -1.643,
      status = OPEN_FAILS_R200_AS_PERTURBATIVE_PREON.
    C_deconstructed_10x10_link:
      sum_T = 91, b10 = 67, Lambda_LP/MG = 48.4,
      alpha_G^{-1}(200 MG) = -15.135, status = FAIL_HIGH_R.
    D_cutoff_shared_U_NLSM:
      sum_T = 71, b10 = 47, Lambda_LP/MG = 252,
      alpha_G^{-1}(200 MG) = 1.730, status = PASS_CONDITIONAL.

paper update:
  Added "NLSM UV Completion/No-Go Scorecard" to paper/gut_framework.tex.
  Updated the No-Go and Conditional Closure ledger text from 17 to 18 claims.

conditional theorem ledger update:
  Updated code/audit_conditional_theorem_ledger.py to include the new row:
    "Minimal A-D routes give a microscopic high-R origin for A5-A6"
  with status NO_GO.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 6
    TUNED_FALLBACK: 2

current obstacle:
  No tested minimal microscopic route is both first-principles enough and
  perturbative through R=50,200.  The only route that preserves the branch is
  the cutoff shared-U NLSM, so A5-A6 remain conditional assumptions rather
  than derived facts.

next attempted nontrivial idea:
  Stop trying to force A5-A6 into a weakly coupled high-R source sector unless
  a genuinely new mechanism appears.  The next high-impact tasks should return
  to the remaining phenomenological blockers:
    O1. build the full joint CKM/PMNS/flavor fit using the covariant CP1 O(2)
        Yukawa functional and stored benchmark cards;
    O2. compute channel-specific dressed d=5 proton decay, especially
        p -> K+ anti-nu, with the current triplet-filter Wilson tensors;
    O3. keep A5-A6 explicitly labelled as cutoff NLSM/constrained-source
        assumptions in the theorem statement.

verification plan:
  1. Compile the TeX after the new scorecard insertion.
  2. Regenerate the conditional ledger and confirm the 18-claim count.
  3. If the next heartbeat continues, start O2 first because proton decay is
     the most likely branch-killer; then run O1 flavor fit if O2 survives.

verification completed:
  py_compile passed for:
    code/audit_nlsm_uv_completion_scorecard.py
    code/audit_conditional_theorem_ledger.py
  audit_nlsm_uv_completion_scorecard.py regenerated output/nlsm_uv_completion_scorecard/*.
  audit_conditional_theorem_ledger.py regenerated the 18-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice after replacing the one
  remaining \over instance with \frac.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries.

Heartbeat update 2026-05-08 01:48 UTC / 09:48 Asia/Taipei:

status:
  The full first-principles GUT task remains incomplete.  O2 is now partially
  advanced: the dimension-five proton calculation has been upgraded from a
  common scalar lifetime estimate to a channel-specific dressed proxy for
  the monitored K nu, K0 mu, and RRRR rows.  It is still not the final
  paper-level calculation because the true SUSY spectrum, wino/higgsino
  dressing matrices, and channel-specific hadronic reduction are not yet
  derived.

new audit:
  Added code/audit_dressed_dimension5_channels.py.

  Outputs:
    output/dressed_dimension5_channels/summary.json
    output/dressed_dimension5_channels/channel_dressing_rows.csv
    output/dressed_dimension5_channels/report.md

  Formula:
    C5_ch = S_T A_ch / M_T,
    C6_ch = C5_ch (alpha2/4pi) m_wino/m_sfermion^2 d_ch,
    Gamma_ch = K_ch |C6_ch|^2.

  Inputs:
    M_T = M_HC(R=200) = 9.442163e15 GeV,
    S_T(display) = 1e-5,
    baseline dressing = m_sfermion 100 TeV, m_wino 1 TeV, alpha2^{-1}=25.

baseline kappa=100 result:
  For omega_R=0.1, kappa_max=100, the monitored baseline dressed channels all
  pass at S_T=1e-5.

  Worst baseline row:
    channel = RRRR_uusd_anycharged proxy,
    amplitude = 7.709749e-5,
    S_T_max = 1.0440755e-2,
    tau(S_T=1e-5) = 2.616225e40 yr,
    margin over 2.4e34 yr = 1.090094e6.

  All-block kappa=100 worst row:
    channel = RRRR_uusd_anycharged proxy,
    amplitude = 9.612328e-5,
    S_T_max = 8.374205e-3,
    tau(S_T=1e-5) = 1.683055e40 yr,
    margin over 2.4e34 yr = 7.012731e5.

stress result:
  The same kappa=100 branch fails an extreme simultaneous stress with
  m_sfermion = 10 TeV and relative dressing enhancement = 100:
    S_T_max = 1.0440755e-6,
    tau(S_T=1e-5) = 2.616225e32 yr,
    margin = 1.09e-2.

paper update:
  Added "Channel-specific dressed d=5 proxy" to paper/gut_framework.tex.
  Updated the final ledger text from 18 to 19 claims.

conditional theorem ledger update:
  Updated code/audit_conditional_theorem_ledger.py with a new row:
    "Monitored dressed d=5 channels pass at S_T=1e-5 in the baseline spectrum"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 7
    TUNED_FALLBACK: 2

current obstacle:
  O2 is not fully closed.  The dressed channel proxy supports the heavy
  100 TeV baseline branch, but the branch is sensitive to light-sfermion and
  enhanced-dressing stresses.  A true publication-level proton-decay result
  still needs explicit wino/higgsino dressing functions from the SUSY spectrum
  and channel-specific hadronic/chiral factors.

next attempted nontrivial idea:
  Build a spectrum-aware dressing audit:
    1. parameterize wino and higgsino dressing matrices by m_sfermion,
       m_wino, mu_H, tan beta, and flavor-rotation insertions;
    2. use the existing two-sided W_AB triplet matrix to compute C_L and C_R
       channel coefficients separately;
    3. scan the dressing-enhancement plane to draw the true proton-safe
       region in (S_T, m_sfermion, tan beta, mu_H);
    4. if the allowed region requires S_T << 1e-5, feed that stronger bound
       back into the RGE/threshold branch.

verification plan:
  1. Compile the TeX after the new d=5 proxy insertion.
  2. Regenerate the conditional theorem ledger and verify 19 claims.
  3. In the next heartbeat, either start the spectrum-aware dressing scan or,
     if O2 is accepted as a proxy, move to O1 full CKM/PMNS/flavor fit.

verification completed:
  py_compile passed for:
    code/audit_dressed_dimension5_channels.py
    code/audit_conditional_theorem_ledger.py
  audit_dressed_dimension5_channels.py regenerated output/dressed_dimension5_channels/*.
  audit_conditional_theorem_ledger.py regenerated the 19-claim ledger.
  Initial TeX compile exposed two notation issues introduced in the new
  section: one \over use and two undefined \TeV macros.  These were replaced
  by \frac and {\rm TeV}.
  pdflatex then compiled paper/gut_framework.tex successfully twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 01:56 UTC / 09:56 Asia/Taipei:

status:
  The full first-principles GUT task remains incomplete.  O2 is now further
  advanced but not fully closed: the dressed dimension-five proton proxy has
  been upgraded from a fixed relative-dressing stress to a spectrum-aware
  wino/higgsino scan with an explicit coherent RRRR stress parameter xi_H.
  This sharpens the proton-decay risk map, but still does not replace the
  final mass-insertion Wilson-tensor calculation.

new audit:
  Added code/audit_spectrum_aware_d5_dressing.py.

  Outputs:
    output/spectrum_aware_d5_dressing/summary.json
    output/spectrum_aware_d5_dressing/spectrum_scan.csv
    output/spectrum_aware_d5_dressing/channel_scan.csv
    output/spectrum_aware_d5_dressing/report.md

  Formula:
    F(x) = (1 - x + x log x)/(1 - x)^2 with F(1)=1/2,
    D_W = (alpha2/4pi) m_wino/m_sfermion^2 F(m_wino^2/m_sfermion^2),
    D_H = (1/16pi^2) mu_H/m_sfermion^2 y_t y_b (tan beta/10)
          F(mu_H^2/m_sfermion^2),
    D_L = D_W,
    D_R = D_W + xi_H D_H.

  Grid:
    m_sfermion = 10,20,50,100,200,500 TeV,
    m_wino = 0.5,1,3 TeV,
    mu_H = 0.5,1,3,10 TeV,
    tan beta = 5,10,30,50,
    xi_H = 1,3,10,30,100.

key result:
  Preferred filter omega_R=0.1, kappa_max=100:
    total stressed safe points = 1313/1440,
    unstressed xi_H=1 safe points = 288/288,
    baseline point = m_sfermion 100 TeV, m_wino 1 TeV, mu_H 1 TeV,
                     tan beta 10, xi_H 1,
    baseline worst channel = RRRR proxy,
    baseline S_T_max = 9.122911e-3,
    baseline margin tau(S_T=1e-5)/tau_bound = 8.322751e5.

  Coherent-higgsino stress table for the preferred filter:
    xi_H=1:   288/288 safe, worst safe margin 2.885
    xi_H=3:   285/288 safe, worst safe margin 1.330, nearest unsafe 0.839
    xi_H=10:  273/288 safe, worst safe margin 1.139, nearest unsafe 0.832
    xi_H=30:  249/288 safe, worst safe margin 1.120, nearest unsafe 0.860
    xi_H=100: 218/288 safe, worst safe margin 1.002, nearest unsafe 0.921

  Most marginal preferred-filter safe point:
    (m_sfermion, m_wino, mu_H, tan beta, xi_H) = (20 TeV, 1 TeV, 3 TeV, 10, 100),
    S_T_max = 1.001203e-5,
    margin = 1.002408.

  Nearest unsafe preferred-filter point:
    (20 TeV, 3 TeV, 3 TeV, 10, 100),
    S_T_max = 9.594875e-6,
    margin = 0.920616.

paper update:
  Added "Spectrum-aware dressed d=5 proxy" to paper/gut_framework.tex.
  Updated final ledger text from 19 to 20 claims.

conditional theorem ledger update:
  Updated code/audit_conditional_theorem_ledger.py with a new row:
    "Spectrum-aware d=5 dressing scan preserves the unstressed branch"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 8
    TUNED_FALLBACK: 2

current obstacle:
  Exact O2 remains open.  The new scan shows that xi_H=1 is safe across the
  scanned spectrum and that failures require coherent RRRR higgsino stress,
  but the physical xi_H must be derived from actual SUSY flavor rotations,
  squark/slepton mass matrices, wino/higgsino dressing loop functions, and
  channel-specific chiral/lattice hadronic factors.

next attempted nontrivial idea:
  Replace xi_H by explicit mass-insertion tensors:
    1. introduce minimal soft-spectrum cards for squark/slepton mass matrices,
       A-terms, mu_H, and gaugino masses;
    2. construct the wino and higgsino dressing kernels for LLLL and RRRR
       operators as matrices rather than scalar stress factors;
    3. contract those kernels with the two-sided W_AB triplet-filter tensor
       and CKM/PMNS-like rotations;
    4. regenerate channel Wilson coefficients for K+ nu, K0 mu, and e pi;
    5. if the derived effective xi_H exceeds the safe envelope, feed the
       stronger S_T bound back into the threshold/global scan.

verification completed:
  py_compile passed for:
    code/audit_spectrum_aware_d5_dressing.py
    code/audit_conditional_theorem_ledger.py
  audit_spectrum_aware_d5_dressing.py regenerated output/spectrum_aware_d5_dressing/*.
  audit_conditional_theorem_ledger.py regenerated the 20-claim ledger.
  pdflatex compiled paper/gut_framework.tex successfully twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:05 UTC / 10:05 Asia/Taipei:

status:
  The full first-principles GUT task remains incomplete.  O2 advanced again:
  the scalar xi_H stress has now been replaced by an explicit mass-insertion
  Wilson-tensor proxy.  This is a real improvement, but it is still not a
  final d=5 proton-decay calculation because the actual sfermion
  diagonalization, chargino/neutralino dressing, and chiral/lattice reduction
  remain to be derived.

new audit:
  Added code/audit_mass_insertion_d5_dressing.py.

  Outputs:
    output/mass_insertion_d5_dressing/summary.json
    output/mass_insertion_d5_dressing/mass_insertion_scan.csv
    output/mass_insertion_d5_dressing/channel_scan.csv
    output/mass_insertion_d5_dressing/report.md

  Formula:
    C'_{abcd}=R_1{}_{aa'}R_2{}_{bb'}R_3{}_{cc'}R_4{}_{dd'}C_{a'b'c'd'},
    with R_i=1+epsilon Delta_i.

  Insertions:
    zero, up_aligned_LL_MFV, down_aligned_LL_MFV, commutator_LL,
    right_third_split, combined_LL_RR, democratic_RR_stress,
    democratic_all_stress.

  Grid:
    epsilon = 0, 0.01, 0.03, 0.10, 0.30, 1.00;
    spectra = baseline_100TeV, marginal_safe_20TeV, near_unsafe_20TeV.

key result:
  Preferred filter omega_R=0.1, kappa_max=100:
    safe points = 144/144,
    baseline margin = 8.322751e5,
    max xi_eff = 46.316078.

  Worst preferred-filter point:
    scenario = democratic_RR_stress,
    epsilon = 1,
    spectrum = near_unsafe_20TeV,
    S_T_max = 1.927814e-5,
    margin = 3.716466,
    RRRR amplitude ratio = 6.752365.

  All-block filter:
    safe points = 144/144,
    max xi_eff = 43.396436,
    worst margin = 2.676615.

paper update:
  Added "Mass-insertion dressed d=5 proxy" to paper/gut_framework.tex.
  Updated final ledger text from 20 to 21 claims.

conditional theorem ledger update:
  Added row:
    "Explicit mass-insertion d=5 proxy stays within the proton-safe envelope"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 9
    TUNED_FALLBACK: 2

current obstacle:
  Exact O2 remains open.  The mass-insertion proxy is no longer a scalar
  xi_H assumption, but it is still not an actual SUSY spectrum calculation.
  The next missing ingredient is the eigenstate-level sfermion and
  chargino/neutralino dressing map.

next attempted nontrivial idea:
  Build a minimal soft-spectrum card and diagonalize the sfermion/chargino
  blocks:
    1. define soft masses from MFV plus a possible transvectant spurion;
    2. diagonalize squark and slepton mass matrices;
    3. compute loop functions with eigenstate sums;
    4. contract the result with W_AB triplet-filter tensors;
    5. feed the resulting Wilson coefficients into the p -> K nu, K0 mu,
       and e pi lifetime table.

verification completed:
  py_compile passed for:
    code/audit_mass_insertion_d5_dressing.py
    code/audit_conditional_theorem_ledger.py
  audit_mass_insertion_d5_dressing.py regenerated output/mass_insertion_d5_dressing/*.
  audit_conditional_theorem_ledger.py regenerated the 21-claim ledger.
  pdflatex compiled paper/gut_framework.tex successfully twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:13 UTC / 10:13 Asia/Taipei:

status:
  The full first-principles GUT task remains incomplete.  O2 advanced one
  layer beyond mass insertions: the d=5 dressing audit now diagonalizes
  positive soft mass matrices and applies pairwise eigenstate loop kernels to
  the four-index Wilson tensors.  This is still a controlled proxy rather
  than a final chargino/neutralino and chiral/lattice proton-decay calculation.

new audit:
  Added code/audit_eigenstate_d5_dressing.py.

  Outputs:
    output/eigenstate_d5_dressing/summary.json
    output/eigenstate_d5_dressing/eigenstate_scan.csv
    output/eigenstate_d5_dressing/channel_scan.csv
    output/eigenstate_d5_dressing/report.md

  Formula:
    M_f^2 = m0^2 (1 + epsilon Delta_f),
    U_f^\dagger M_f^2 U_f = m0^2 diag(lambda_f).
    For each dressed slot pair p,q:
      K_ab,a'b' = sum_rs U^p_ar U^{p*}_{a'r} U^q_bs U^{q*}_{b's}
                  D_chi(lambda_pr, lambda_qs).
    All six pairings are audited and the most dangerous pairing is retained.

  Grid:
    epsilon = 0, 0.01, 0.03, 0.10, 0.30, 0.60, 0.90;
    spectra = baseline_100TeV, marginal_safe_20TeV, near_unsafe_20TeV;
    soft ansaetze = zero, up_aligned_LL_MFV, down_aligned_LL_MFV,
      commutator_LL, right_third_split, combined_LL_RR,
      democratic_RR_stress, democratic_all_stress.

key result:
  Preferred filter omega_R=0.1, kappa_max=100:
    aligned safe = 126/126,
    positive-grid safe = 168/168,
    invalid positive-definiteness points = 0,
    max amplification versus degenerate soft spectrum = 10.274057.

  Most marginal preferred-filter row:
    scenario = democratic_all_stress,
    epsilon = 0.9,
    spectrum = near_unsafe_20TeV,
    min soft eigenvalue = 0.1,
    worst channel = LLLL_upupdown_Knu,
    worst pair = 23,
    S_T_max = 2.186675e-5,
    margin = 4.781546.

  All-block filter:
    positive-grid safe = 168/168,
    max amplification = 19.735894,
    worst margin = 1.760275.

paper update:
  Added "Soft-eigenstate dressed d=5 proxy" to paper/gut_framework.tex.
  Updated final ledger text from 21 to 22 claims.

conditional theorem ledger update:
  Added row:
    "Positive soft-eigenstate d=5 proxy stays within the proton-safe envelope"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 10
    TUNED_FALLBACK: 2

current obstacle:
  Exact O2 remains open.  The eigenstate proxy respects positive soft spectra
  and pairwise loop enhancement, but it still omits the real chargino and
  neutralino mixing matrices, interference between diagrams, chiral
  reduction, and lattice hadronic matrix elements.

next attempted nontrivial idea:
  Replace pairwise scalar loop kernels by actual MSSM chargino/neutralino
  dressing blocks:
    1. build chargino and neutralino mass matrices from M2, M1, mu_H, tan beta;
    2. diagonalize them and compute eigenstate couplings to sfermions;
    3. sum the dressed Wilson coefficients with diagram signs for LLLL/RRRR;
    4. introduce a minimal chiral reduction matrix for p -> K nu, K0 mu, e pi;
    5. compare the exact-dressing result against the eigenstate proxy envelope.

verification completed:
  py_compile passed for:
    code/audit_eigenstate_d5_dressing.py
    code/audit_conditional_theorem_ledger.py
  audit_eigenstate_d5_dressing.py regenerated output/eigenstate_d5_dressing/*.
  audit_conditional_theorem_ledger.py regenerated the 22-claim ledger.
  pdflatex compiled paper/gut_framework.tex successfully twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:19 UTC / 10:19 Asia/Taipei:

status:
  The full first-principles GUT task remains incomplete.  O2 advanced from
  positive soft-eigenstate dressing to an explicit chargino/neutralino mixing
  proxy.  The preferred omega_R=0.1 triplet filter still survives, while the
  all-block fallback is now visibly less robust.

new audit:
  Added code/audit_mssm_mixing_d5_dressing.py.

  Outputs:
    output/mssm_mixing_d5_dressing/summary.json
    output/mssm_mixing_d5_dressing/mssm_mixing_scan.csv
    output/mssm_mixing_d5_dressing/channel_scan.csv
    output/mssm_mixing_d5_dressing/report.md

  Formula:
    Chargino matrix:
      X_C = [[M2, sqrt(2) mW sin beta],
             [sqrt(2) mW cos beta, mu_H]].
    Neutralino matrix in (B,W0,Hd0,Hu0), with
      M1 = (5/3) tan^2(theta_W) M2.
    The pairwise soft eigenstate kernel coherently sums absolute chargino and
    neutralino contributions.  This is deliberately conservative and does not
    use phase cancellations.

key result:
  Preferred filter omega_R=0.1, kappa_max=100:
    aligned safe = 126/126,
    positive-grid safe = 168/168,
    max amplification versus soft-eigenstate proxy = 2.015148.

  Most marginal preferred-filter row:
    scenario = democratic_all_stress,
    epsilon = 0.9,
    spectrum = near_unsafe_20TeV,
    M1 = 1.503811 TeV,
    m_chargino1 = 2.938222 TeV,
    m_neutralino1 = 1.503212 TeV,
    worst channel = LLLL_upupdown_Knu,
    worst pair = 23,
    S_T_max = 1.085193e-5,
    margin = 1.177643.

  All-block fallback:
    positive-grid safe = 167/168,
    failing row = democratic_all_stress, epsilon=0.9, near_unsafe_20TeV,
    S_T_max = 6.617190e-6,
    margin = 0.437872.
    Interpretation: all-block is no longer equally trustworthy under
    coherent electroweakino mixing; use the omega_R=0.1 preferred branch.

paper update:
  Added "Chargino/neutralino-mixing dressed d=5 proxy" to
  paper/gut_framework.tex.
  Updated final ledger text from 22 to 23 claims.

conditional theorem ledger update:
  Added row:
    "MSSM chargino/neutralino mixing d=5 proxy preserves the preferred filter"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 11
    TUNED_FALLBACK: 2

current obstacle:
  Exact O2 remains open.  The chargino/neutralino eigencontent is now included
  in a coherent-absolute proxy, but the final result still needs signed
  interference-aware Wilson coefficients, channel-specific chiral reduction,
  and lattice/hadronic normalization.

next attempted nontrivial idea:
  Build an interference ledger for the d=5 amplitudes:
    1. keep the complex Wilson tensor entries instead of absolute coherent sums;
    2. assign diagram signs for chargino, neutralino, and higgsino dressing;
    3. define a minimal chiral channel map for K+ nu, K0 mu, and e pi;
    4. scan CP phases in the triplet filter and electroweakino sectors;
    5. classify whether the preferred branch is safe without relying on
       accidental destructive interference.

verification completed:
  py_compile passed for:
    code/audit_mssm_mixing_d5_dressing.py
    code/audit_conditional_theorem_ledger.py
  audit_mssm_mixing_d5_dressing.py regenerated output/mssm_mixing_d5_dressing/*.
  audit_conditional_theorem_ledger.py regenerated the 23-claim ledger.
  pdflatex compiled paper/gut_framework.tex successfully twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:24 UTC / 10:24 Asia/Taipei:

status:
  The full first-principles GUT task remains incomplete.  O2 advanced from
  coherent-absolute MSSM mixing to a signed complex phase-grid interference
  proxy.  The preferred omega_R=0.1 branch remains safe, but now only
  narrowly at the worst constructive phase.

new audit:
  Added code/audit_signed_interference_d5.py.

  Outputs:
    output/signed_interference_d5/summary.json
    output/signed_interference_d5/signed_interference_scan.csv
    output/signed_interference_d5/channel_scan.csv
    output/signed_interference_d5/report.md

  Formula:
    Complex component tensors are retained:
      C = C_chargino + exp(i phi_N) C_neutralino + exp(i phi_H) C_higgsino.
    Proxy relative signs:
      LLLL: C_chargino - exp(i phi_N) C_neutralino.
      RRRR: exp(i phi_N) C_neutralino - exp(i phi_H) C_higgsino.
    Triplet-filter CP deformation:
      W(phi_T)=P(phi_T) W P(phi_T)^dagger,
      P=diag(1,1,exp(i phi_T),exp(2 i phi_T)).
    The scan uses phi_T=0,2pi/3,4pi/3 and
      phi_N, phi_H = 0, pi/2, pi, 3pi/2.
    Minimal e-pi proxies LLLL_uude and RRRR_uude were added.

key result:
  Preferred filter omega_R=0.1, kappa_max=100:
    aligned safe = 126/126,
    positive-grid safe = 168/168,
    max phase spread = 4.794590e2.

  Most marginal preferred-filter row:
    scenario = democratic_all_stress,
    epsilon = 0.9,
    spectrum = near_unsafe_20TeV,
    worst triplet phase = 4pi/3,
    worst channel = LLLL_upupdown_Knu,
    worst pair = 13,
    S_T_max = 1.007201e-5,
    margin = 1.014453.

  All-block fallback:
    positive-grid safe = 167/168,
    failing row remains democratic_all_stress, epsilon=0.9, near_unsafe_20TeV,
    S_T_max = 6.617190e-6,
    margin = 0.437872.

paper update:
  Added "Signed-interference d=5 phase-grid proxy" to
  paper/gut_framework.tex.
  Updated final ledger text from 23 to 24 claims.

conditional theorem ledger update:
  Added row:
    "Signed CP-phase d=5 interference proxy preserves the preferred filter"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 12
    TUNED_FALLBACK: 2

current obstacle:
  Exact O2 remains open.  The preferred branch no longer appears to rely on
  accidental destructive interference in the audited phase grid, but the
  margin is thin.  The remaining missing physics is channel-specific chiral
  reduction and lattice/hadronic normalization for K+ nu, K0 mu, and e pi.

next attempted nontrivial idea:
  Build a chiral/lattice normalization ledger:
    1. separate short-distance Wilson coefficients from hadronic matrix
       elements;
    2. introduce symbolic beta_H, D, F, and channel kinematic factors;
    3. propagate a conservative uncertainty band rather than one scalar
       prefactor;
    4. replay signed-interference worst rows against K+ nu, K0 mu, e pi;
    5. decide whether S_T=1e-5 remains viable or must be tightened.

verification completed:
  py_compile passed for:
    code/audit_signed_interference_d5.py
    code/audit_conditional_theorem_ledger.py
  audit_signed_interference_d5.py regenerated output/signed_interference_d5/*.
  audit_conditional_theorem_ledger.py regenerated the 24-claim ledger.
  pdflatex compiled paper/gut_framework.tex successfully twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:30 UTC / 10:30 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The O2 d=5 proton
  decay audit advanced from a single toy hadronic prefactor to a
  channel-specific chiral/lattice normalization replay.  The central replay
  passes, but the conservative maximum-width envelope fails the previous
  S_T=1e-5 target.  Therefore the d=5 branch is now a tuned fallback unless
  the triplet filter is tightened.

new audit:
  Added code/audit_chiral_lattice_d5_normalization.py.

  Outputs:
    output/chiral_lattice_d5/summary.json
    output/chiral_lattice_d5/chiral_replay.csv
    output/chiral_lattice_d5/report.md

  Formula:
    K_ch = m_p/(64 pi f_pi^2) beta_H^2 A_R^2 C_ch^2 Phi_2.

    C(e pi) = |1 + D + F|.
    C(K nu) = sqrt((1 + (D + 3F)/3)^2 + (2D/3)^2).
    C(K0 mu) = sqrt((1 - D + F)^2 + (2D/3)^2).

  Uncertainty box:
    beta_H = 0.008, 0.012, 0.018 GeV^3.
    A_R = 1.8, 2.4, 3.2.
    D = 0.60, 0.80, 1.00.
    F = 0.35, 0.47, 0.55.

key result:
  Preferred omega_R=0.1, kappa=100 filter:
    central global S_T_max = 1.728394e-5,
    central worst margin = 2.987346,
    central unsafe rows = 0.

    min-width global S_T_max = 3.950075e-5,
    min-width worst margin = 15.603089,
    min-width unsafe rows = 0.

    max-width global S_T_max = 7.693048e-6,
    max-width worst margin = 0.591830,
    max-width unsafe rows = 8.

  Thus S_T=1e-5 is not robust under the conservative hadronic envelope.
  A robust target is:
    S_T <= 7.693048e-6.

  Nearest max-width unsafe preferred-filter row:
    channel = LLLL_downdownup_Knu,
    scenario = democratic_all_stress,
    epsilon = 0.9,
    spectrum = near_unsafe_20TeV,
    triplet_phase = 2pi/3,
    pair = 12,
    S_T_max = 9.752546e-6,
    margin = 0.951122.

  All-block fallback:
    central S_T_max = 1.143347e-5,
    max-width S_T_max = 5.198437e-6,
    max-width worst margin = 0.270238,
    max-width unsafe rows = 8.

paper update:
  Added the "Chiral/lattice normalization replay" proposition to
  paper/gut_framework.tex.
  Updated final ledger text from 24 to 25 claims and tuned fallbacks from
  2 to 3.

conditional theorem ledger update:
  Added row:
    "The d=5 filter target is robust under the conservative chiral/lattice
    envelope"
  with status TUNED_FALLBACK.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 12
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open.  The d=5 channel proof now requires propagating the
  tighter S_T target through the triplet filter / missing-partner regulator
  and checking whether the threshold/RGE branch survives.

next attempted nontrivial idea:
  Tightened-filter replay:
    1. set S_T target to 7.5e-6, or use the exact 7.693048e-6 bound;
    2. rerun the signed-interference and chiral worst rows;
    3. replay the triplet-filter charge/regulator floor and the
       threshold/proton scan;
    4. check whether the R=200 preferred branch and missing-partner
       kappa=100 can naturally supply the tighter filter;
    5. if not, derive an extra suppression source or declare d=5 proton
       safety tuned.

verification completed:
  py_compile passed for:
    code/audit_chiral_lattice_d5_normalization.py
    code/audit_conditional_theorem_ledger.py
  audit_chiral_lattice_d5_normalization.py regenerated
  output/chiral_lattice_d5/*.
  audit_conditional_theorem_ledger.py regenerated the 25-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:35 UTC / 10:35 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The immediate O2
  tightened-filter replay has been completed: the conservative chiral/lattice
  envelope can be closed for the preferred omega_R=0.1,kappa=100 triplet
  filter by replacing the display target S_T=1e-5 with S_T=7.5e-6.  This is a
  filter-specific conditional closure, not a generic d=5 proton-decay proof.

new audit:
  Added code/audit_tightened_triplet_filter_replay.py.

  Outputs:
    output/tightened_triplet_filter/summary.json
    output/tightened_triplet_filter/tightened_replay.csv
    output/tightened_triplet_filter/report.md

  Mathematical replay:
    tau(S_T) / tau_bound = margin(S0) (S0/S_T)^2, with S0=1e-5.
    Equivalently, a row passes iff S_T^target <= S_T,max(row).

  Targets:
    display_1e-5 = 1.0e-5.
    exact_max_width_preferred = 7.693047923553811e-6.
    rounded_7p5e-6 = 7.5e-6.

key result:
  Conservative max-width hadronic envelope:
    preferred omega_R=0.1,kappa=100 at S_T=1e-5:
      worst margin = 0.591830, unsafe rows = 8.
    preferred omega_R=0.1,kappa=100 at exact 7.693047923553811e-6:
      worst margin = 1.000000, unsafe rows = 0.
    preferred omega_R=0.1,kappa=100 at rounded 7.5e-6:
      worst margin = 1.052142, unsafe rows = 0.
    all-block kappa=100 fallback at rounded 7.5e-6:
      worst margin = 0.480422, unsafe rows = 5.

  Interpretation:
    The preferred filter closes the conservative d=5 proxy after tightening,
    but the all-block fallback still fails.  Therefore downstream proton tables
    must use S_T=7.5e-6 for the preferred branch, and cannot claim a generic
    all-block d=5 proof.

paper update:
  Added the "Tightened triplet-filter replay" proposition to
  paper/gut_framework.tex.
  Corrected the C(K nu) chiral-factor formula to include the missing plus sign
  between the direct and indirect squared terms.
  Updated the final ledger text from 25 to 26 claims and passing/conditional
  checks from 16 to 17.
  Updated O2 wording to "full field-basis d=5 proton decay with physical flavor
  and lattice inputs".

conditional theorem ledger update:
  Added row:
    "Tightened triplet filter closes the conservative d=5 chiral/lattice
    envelope"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 13
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open at the publication level.  The proxy now has a
  tightened preferred-filter closure, but it still lacks a full field-basis
  Wilson coefficient calculation with physical flavor rotations, lattice
  matrix elements, and dressed p -> K+ nu / K0 mu / e pi amplitudes.

next attempted nontrivial idea:
  Physical channel Wilson replay:
    1. propagate S_T=7.5e-6 into the RGE/proton benchmark cards;
    2. build a field-basis d=5 Wilson module with explicit CKM/PMNS rotations
       and LLLL/RRRR flavor tensors;
    3. separate short-distance running, dressing matrices, and hadronic matrix
       elements as independent uncertainty knobs;
    4. compare the tightened proxy against physical p -> K+ nu, K0 mu, and
       e+ pi0 amplitudes;
    5. if the physical replay needs S_T below 7.5e-6, either derive a stronger
       triplet-filter symmetry or downgrade d=5 safety to tuned.

verification completed:
  py_compile passed for:
    code/audit_tightened_triplet_filter_replay.py
    code/audit_conditional_theorem_ledger.py
  audit_tightened_triplet_filter_replay.py generated
  output/tightened_triplet_filter/*.
  audit_conditional_theorem_ledger.py regenerated the 26-claim ledger.
  pdflatex compiled paper/gut_framework.tex successfully.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:46 UTC / 10:46 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The tightened
  S_T=7.5e-6 filter has now been propagated into the corrected downstream
  RGE/proton benchmark cards.  The preferred R=200 branch remains alive for
  the current 2.4e34 yr d=5 stress target, but not for a future 1e35 yr
  stress target.

new audit:
  Added code/audit_tightened_downstream_proton.py.

  Outputs:
    output/tightened_downstream_proton/summary.json
    output/tightened_downstream_proton/benchmark_replay.csv
    output/tightened_downstream_proton/report.md

  Mathematical replay:
    tau_d5(S_T) = tau_d5(1e-5) (1e-5/S_T)^2.
    For S_T=7.5e-6, the lifetime scaling factor is 1.777778.

key result:
  Corrected downstream cards at S_T=7.5e-6:
    corrected_base:
      tau_d5 = 5.179406e34 yr,
      margin over 2.4e34 = 2.158086,
      margin over 1e35 = 0.517941,
      triplet_filter_required/S_T = 2.275831.
    fixed_R50:
      tau_d5 = 5.170694e34 yr,
      margin over 2.4e34 = 2.154456,
      margin over 1e35 = 0.517069,
      triplet_filter_required/S_T = 2.273916.
    R_window_R200:
      tau_d5 = 5.177232e34 yr,
      margin over 2.4e34 = 2.157180,
      margin over 1e35 = 0.517723,
      triplet_filter_required/S_T = 2.275353.
    goldstone_locked_R200:
      tau_d5 = 5.175203e34 yr,
      margin over 2.4e34 = 2.156335,
      margin over 1e35 = 0.517520,
      triplet_filter_required/S_T = 2.274907.

  Interpretation:
    The corrected RGE/proton cards remain alive after replacing the display
    S_T=1e-5 by S_T=7.5e-6.  The R=200 branch exceeds the current 2.4e34 yr
    stress target by a factor about 2.16, but would need another suppression
    mechanism to pass a 1e35 yr future stress target.

paper update:
  Added the "Tightened downstream proton replay" proposition to
  paper/gut_framework.tex.
  Updated final ledger text from 26 to 27 claims and passing/conditional
  checks from 17 to 18.

conditional theorem ledger update:
  Added row:
    "Tightened triplet filter remains compatible with downstream RGE/proton
    benchmark cards"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 14
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open.  The tightened proxy is now consistent with the
  threshold/RGE card, but the physical p -> K+ nu / K0 mu / e pi amplitudes
  still need a field-basis Wilson calculation with CKM/PMNS rotations,
  dressing matrices, short-distance running, and lattice/chiral inputs.

next attempted nontrivial idea:
  Physical d=5 Wilson module:
    1. build explicit LLLL and RRRR Wilson tensors in the mass basis;
    2. separate triplet filter, CKM/PMNS rotations, and sfermion/chargino
       dressing into independent matrices;
    3. attach channel-specific hadronic matrix elements as symbolic or scanned
       constants;
    4. reproduce the tightened proxy in the diagonal limit;
    5. then test whether realistic physical mixing demands S_T below 7.5e-6
       or remains within the present R=200 margin.

verification completed:
  py_compile passed for:
    code/audit_tightened_downstream_proton.py
    code/audit_conditional_theorem_ledger.py
  audit_tightened_downstream_proton.py generated
  output/tightened_downstream_proton/*.
  audit_conditional_theorem_ledger.py regenerated the 27-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 02:56 UTC / 10:56 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The next d=5 proton
  step has now moved from scalar proxy margins to a first physical-field-basis
  Wilson replay.  At common 100 TeV dressing, CKM/PMNS flavor rotations do
  not destroy the tightened S_T=7.5e-6 proton-safe margin.

new audit:
  Added code/audit_physical_d5_wilson_replay.py.

  Outputs:
    output/physical_d5_wilson_replay/summary.json
    output/physical_d5_wilson_replay/physical_d5_wilson_replay.csv
    output/physical_d5_wilson_replay/report.md

  Mathematical replay:
    S_T,max^ch = S_T,max^old sqrt(K_old/K_ch).
    tau/tau_bound = (S_T,max^ch / 7.5e-6)^2.
    The replay uses the triplet-mixing Wilson-tensor rows, channel-dependent
    chiral/lattice width prefactors, and the tightened filter target.

key result:
  Conservative max-width chiral/lattice envelope:
    diagonal_rep_mixing_identity:
      worst margin = 2.471619e2,
      extra amplitude headroom = 1.572138e1,
      worst channel = RRRR_uusd_anycharged.
    diagonal_rep_mixing_knu_null:
      worst margin = 9.931894e4,
      extra amplitude headroom = 3.151491e2,
      worst channel = LLLL_upupdown_Knu.
    full_bipartite_mixing_identity:
      worst margin = 2.471619e2,
      extra amplitude headroom = 1.572138e1,
      worst channel = RRRR_uusd_anycharged.
    full_bipartite_mixing_knu_null:
      worst margin = 9.947384e4,
      extra amplitude headroom = 3.153947e2,
      worst channel = RRRR_uusd_anycharged.

  All replay rows have zero unsafe channels at S_T=7.5e-6.

interpretation:
  The bare field-basis CKM/PMNS Wilson rotation is not the limiting d=5
  obstruction.  The identity triplet-mixing branch already has O(10)
  amplitude headroom under the conservative width envelope.  The near-null
  Knu branches have much larger headroom but are rank deficient, so they are
  diagnostics rather than a standalone UV mechanism.

paper update:
  Added the "Field-basis Wilson replay" proposition to paper/gut_framework.tex.
  Updated final ledger text from 27 to 28 claims and passing/conditional
  checks from 18 to 19.

conditional theorem ledger update:
  Added row:
    "Field-basis CKM/PMNS d=5 Wilson replay passes the tightened filter at
    common dressing"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 15
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open.  The physical-field Wilson replay is still tied to a
  common 100 TeV dressing normalization.  A publication-level proof still
  needs full soft-spectrum dressing: chargino/neutralino eigenstates,
  sfermion flavor/eigenstate rotations, short-distance running, lattice/chiral
  hadronic inputs, and an action-level UV origin for the S_T=7.5e-6 filter.

next attempted nontrivial idea:
  Build a full soft-spectrum d=5 Wilson audit:
    1. take the mass-basis LLLL/RRRR Wilson tensors from the physical replay;
    2. replace the common dressing factor by explicit chargino/neutralino and
       higgsino dressing matrices;
    3. diagonalize positive sfermion mass matrices and attach flavor rotations;
    4. keep the chiral/lattice max-width envelope as the conservative
       hadronic stress;
    5. check whether the corrected R=200 branch still passes at
       S_T=7.5e-6, and record the minimum additional suppression needed for a
       future 1e35 yr stress target.

verification completed:
  py_compile passed for:
    code/audit_physical_d5_wilson_replay.py
    code/audit_conditional_theorem_ledger.py
  audit_physical_d5_wilson_replay.py generated
  output/physical_d5_wilson_replay/*.
  audit_conditional_theorem_ledger.py regenerated the 28-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 03:25 UTC / 11:25 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The d=5 proton sector
  has moved one step deeper: the tightened S_T=7.5e-6 filter has now been
  replayed through the existing chargino/neutralino and positive-sfermion
  eigenstate dressing proxy, with the conservative chiral/lattice max-width
  envelope included.

new audit:
  Added code/audit_soft_spectrum_d5_replay.py.

  Outputs:
    output/soft_spectrum_d5_replay/summary.json
    output/soft_spectrum_d5_replay/channel_replay.csv
    output/soft_spectrum_d5_replay/point_replay.csv
    output/soft_spectrum_d5_replay/report.md

  Mathematical replay:
    margin_new = margin_old * (S_T_old/S_T_new)^2 * (K_old/K_channel).
    Here S_T_old=1e-5, S_T_new=7.5e-6, and K_channel is the
    channel-dependent chiral/lattice central, max-width, or min-width
    prefactor.  Equivalently,
    S_T,max^new = S_T,max^old sqrt(K_old/K_channel).

key result:
  Conservative max-width envelope:
    preferred omegaR=0.1,kappa=100 filter:
      current stress safe points = 168/168 positive soft points,
      worst current margin = 1.292080,
      worst current channel = LLLL_upupdown_Knu,
      worst soft point = democratic_all_stress, epsilon=0.9,
      spectrum = near_unsafe_20TeV,
      replayed S_T,max = 8.525229e-6.
    all-block fallback:
      current stress safe points = 167/168 positive soft points,
      worst current margin = 0.480422,
      so the all-block fallback remains unsafe.

  Uniform future 1e35 yr stress:
    preferred filter:
      safe points = 167/168,
      worst future margin = 0.310099,
      extra amplitude suppression needed = 1.795765.
    all-block fallback:
      safe points = 166/168,
      worst future margin = 0.115301,
      extra amplitude suppression needed = 2.944983.

interpretation:
  The preferred branch survives explicit electroweakino and positive-sfermion
  eigenstate dressing at the current stress target, but the surviving margin
  is thin: only a factor 1.292 in lifetime at the worst audited point.  This
  is enough for the current conditional EFT claim, but not enough for a robust
  future-proof proton-decay claim.

paper update:
  Added the "Soft-spectrum d=5 replay" proposition to paper/gut_framework.tex.
  Updated final ledger text from 28 to 29 claims and passing/conditional
  checks from 19 to 20.

conditional theorem ledger update:
  Added row:
    "Soft-spectrum d=5 Wilson replay passes the tightened filter on the
    audited positive grid"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 16
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open.  The most urgent proton-side gap is no longer a
  scalar filter replay but a publication-level d=5 calculation: short-distance
  operator running from M_G to SUSY thresholds, full sfermion flavor matrices,
  chargino/neutralino/higgsino signed amplitudes, and lattice/chiral hadronic
  matrix elements must be combined in one Wilson-coefficient pipeline.
  Separately, the UV/action-level origin of the required S_T=7.5e-6 filter is
  still not first-principles-derived.

next attempted nontrivial idea:
  Build an interference-aware physical d=5 pipeline:
    1. use the soft-spectrum replay rows as the grid;
    2. replace coherent absolute sums by complex signed chargino, neutralino,
       and higgsino amplitudes with CKM/PMNS field-basis tensors;
    3. propagate short-distance running factors per operator class;
    4. attach channel-specific chiral/lattice matrix elements;
    5. scan phases and soft flavor textures to find whether the worst
       margin 1.292 can be made structurally stable, or whether the theory
       requires a stronger symmetry deriving S_T <= 8.5e-6.

verification completed:
  py_compile passed for:
    code/audit_soft_spectrum_d5_replay.py
    code/audit_conditional_theorem_ledger.py
  audit_soft_spectrum_d5_replay.py generated
  output/soft_spectrum_d5_replay/*.
  audit_conditional_theorem_ledger.py regenerated the 29-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 03:37 UTC / 11:37 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The d=5 proton sector
  now has a tightened worst-phase signed-interference replay.  The preferred
  triplet filter still passes the current stress target without relying on
  destructive CP-phase interference, but the margin is now extremely thin.

new audit:
  Added code/audit_interference_tightened_d5_replay.py.

  Outputs:
    output/interference_tightened_d5_replay/summary.json
    output/interference_tightened_d5_replay/channel_replay.csv
    output/interference_tightened_d5_replay/point_replay.csv
    output/interference_tightened_d5_replay/report.md

  Mathematical replay:
    margin_new = margin_old * (S_T_old/S_T_new)^2 * (K_old/K_channel).
    The input rows are already worst over the scanned neutralino and higgsino
    phases, so the replay is a worst-phase pass/fail test, not an optimistic
    destructive-interference test.

key result:
  Conservative max-width envelope:
    preferred omegaR=0.1,kappa=100 filter:
      current stress safe points = 168/168 positive soft points,
      worst current margin = 1.052142,
      worst channel = LLLL_uude_epi,
      worst soft point = democratic_all_stress, epsilon=0.9,
      spectrum = near_unsafe_20TeV,
      worst triplet phase = 2*pi/3,
      replayed S_T,max = 7.693047923553811e-6,
      max phase spread in audited rows = 4.794590e2.
    all-block fallback:
      current stress safe points = 167/168 positive soft points,
      worst current margin = 0.480422,
      so the all-block fallback remains unsafe.

  Uniform future 1e35 yr stress:
    preferred filter:
      safe points = 167/168,
      worst future margin = 0.252514,
      extra amplitude suppression needed = 1.990019.
    all-block fallback:
      safe points = 165/168,
      worst future margin = 0.115301,
      extra amplitude suppression needed = 2.944983.

interpretation:
  This is the sharpest proton-side audit so far.  The branch is still
  conditionally alive for the current stress target because 1.052142 > 1, and
  the pass uses worst scanned phases rather than phase cancellation.  But the
  adopted S_T=7.5e-6 is only about 2.57 percent below the replayed worst-case
  allowance 7.6930479e-6.  A robust publication-level d=5 claim needs either
  a stronger dynamically derived filter, a better channel treatment that
  raises this margin, or a principled reason why the worst epi monitor is too
  conservative.

paper update:
  Added the "Tightened signed-interference replay" proposition to
  paper/gut_framework.tex.
  Updated final ledger text from 29 to 30 claims and passing/conditional
  checks from 20 to 21.

conditional theorem ledger update:
  Added row:
    "Worst-phase signed d=5 replay passes the tightened filter on the audited
    positive grid"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 17
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open.  On the proton side, the bottleneck is now the
  LLLL_uude_epi monitor under the conservative hadronic max-width envelope.
  The surviving margin is too thin to be called robust.  The UV/action-level
  origin of S_T <= 7.5e-6 also remains unproven.

next attempted nontrivial idea:
  Resolve the epi monitor and filter origin together:
    1. replace the minimal epi monitor by a channel-specific p -> e+ pi0
       Wilson contraction with SU(2) and color Clebsch factors;
    2. include short-distance operator running and lattice/chiral matrix
       elements separately for epi, Knu, and K0mu;
    3. test whether the worst epi margin 1.052 is an artifact of the minimal
       monitor or a real proton-decay bound;
    4. in parallel, derive an action-level or charge-symmetry origin for
       S_T <= 7.5e-6, preferably with enough margin for a 1e35 yr stress.

verification completed:
  py_compile passed for:
    code/audit_interference_tightened_d5_replay.py
    code/audit_conditional_theorem_ledger.py
  audit_interference_tightened_d5_replay.py generated
  output/interference_tightened_d5_replay/*.
  audit_conditional_theorem_ledger.py regenerated the 30-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 04:01 UTC / 12:01 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The proton-side
  bottleneck has been sharpened again: the minimal e-pi monitor from the
  worst-phase signed replay is likely an over-conservative proxy once the
  neutral-pion and identical-up projection factors are included.  The current
  preferred branch survives a little more comfortably, but the future
  1e35 yr d=5 stress still fails.

new audit:
  Added code/audit_epi_clebsch_replay.py.

  Outputs:
    output/epi_clebsch_replay/summary.json
    output/epi_clebsch_replay/channel_replay.csv
    output/epi_clebsch_replay/point_replay.csv
    output/epi_clebsch_replay/report.md

  Mathematical replay:
    Only e_pi rows are rescaled.  Knu and K0mu rows are left unchanged.
    For an e_pi amplitude projection a, widths scale as a^2, lifetime margins
    scale as 1/a^2, and S_T,max scales as 1/a.

    The tested representation-level profiles are:
      raw_monitor: a = 1;
      pi0_isospin: a = 1/sqrt(2);
      identical_up_norm: a = 1/sqrt(2);
      pi0_times_identical_up: a = 1/2.

key result:
  Preferred omegaR=0.1,kappa=100 filter under the conservative max-width
  envelope:
    raw monitor:
      current margin = 1.052142,
      bottleneck = LLLL_uude_epi,
      future 1e35 margin = 0.252514,
      future amplitude suppression needed = 1.990019.
    pi0, identical-up, or combined projection:
      current margin = 1.113033,
      bottleneck = LLLL_upupdown_Knu,
      projected S_T,max = 7.912527e-6,
      future 1e35 margin = 0.267128,
      future amplitude suppression needed = 1.934819.

interpretation:
  The raw e-pi monitor is not a stable physical bottleneck unless a full
  channel-specific contraction restores the missing Clebsch strength.  With
  the natural projection sensitivity included, the limiting current channel
  returns to Knu.  This improves the current conditional claim but does not
  close the d=5 sector: the surviving margin remains thin, and the future
  1e35 yr stress still needs about a factor 1.93 in additional amplitude
  suppression.

paper update:
  Added the "e pi Clebsch sensitivity" proposition to paper/gut_framework.tex.
  Replaced raw TeX \over fractions in the new proposition by \frac to keep the
  LaTeX warning scan clean.
  Updated final ledger text from 30 to 31 claims and passing/conditional
  checks from 21 to 22.

conditional theorem ledger update:
  Added row:
    "Neutral-pion and identical-up projection removes the minimal e-pi
    monitor as the current d=5 bottleneck"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 18
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open.  On the proton side, Knu is now the more robust
  bottleneck after the e-pi projection sensitivity replay.  Two things are
  still missing before a publication-level d=5 claim is honest:
    1. a full channel-specific p -> e+ pi0 contraction, to decide whether the
       e-pi projection factor is truly physical in the present operator basis;
    2. a Knu-focused dressed Wilson calculation and/or a stronger dynamical
       origin for S_T <= 7.5e-6 with enough headroom for a 1e35 yr stress.

next attempted nontrivial idea:
  Build a Knu-centered channel-specific audit:
    1. separate the LLLL uud Knu and RRRR Knu Wilson structures by flavor
       rotation, short-distance running, and chiral/lattice matrix elements;
    2. keep the e-pi Clebsch replay as a side sensitivity, not a final
       theorem;
    3. scan whether the existing Z3/missing-partner triplet filter can be
       strengthened from S_T=7.5e-6 to roughly S_T <= 3.9e-6, which would
       meet the 1e35 yr stress if no further channel suppression appears;
    4. synchronize the result with the action-level triplet-filter origin,
       because a purely numerical filter without a symmetry/source mechanism
       remains a conditional EFT assumption.

verification completed:
  py_compile passed for:
    code/audit_epi_clebsch_replay.py
    code/audit_conditional_theorem_ledger.py
  audit_epi_clebsch_replay.py generated output/epi_clebsch_replay/*.
  audit_conditional_theorem_ledger.py regenerated the 31-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice after the final TeX edit.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Heartbeat update 2026-05-08 04:27 UTC / 12:27 Asia/Taipei
==========================================================

status:
  Full first-principles GUT task remains incomplete.  The current turn
  converted the newly identified Knu bottleneck into a numerical target map:
  the preferred branch survives the current proton stress, but future-proof
  1e35 yr d=5 safety needs a quantitatively stronger suppression.

new audit:
  Added code/audit_knu_target_map.py.

  Outputs:
    output/knu_target_map/summary.json
    output/knu_target_map/channel_targets.csv
    output/knu_target_map/report.md

  Mathematical replay:
    The input is output/epi_clebsch_replay/channel_replay.csv with the
    pi0_times_identical_up profile, omegaR=0.1,kappa=100 preferred filter, and
    conservative max-width hadronic envelope.

    For each channel A:
      current amplitude headroom = sqrt(tau_A/tau_now);
      future amplitude suppression needed =
        max((tau_A/1e35 yr)^(-1/2), 1);
      S_T,max for future stress =
        S_T,target * sqrt(tau_A/1e35 yr).

key result:
  The robust preferred-filter bottleneck is:
    channel = LLLL_upupdown_Knu,
    class = Knu,
    current margin = 1.113033,
    future 1e35 margin = 0.267128,
    current amplitude headroom = 1.055004,
    future amplitude suppression needed = 1.934819,
    future-safe S_T target = 3.876331e-6.

  Nearby channels:
    LLLL_downdownup_Knu:
      current margin = 1.339617,
      future margin = 0.321508,
      future amplitude suppression needed = 1.763616.
    LLLL_upupdown_K0mu:
      future margin = 2.075664, already safe under 1e35 yr.
    LLLL_uude_epi after e-pi projection:
      future margin = 1.010056, just future-safe.
    RRRR_uusd_anycharged:
      future margin = 2.948778, safe.

interpretation:
  This audit makes the proton-side design target precise.  The next d=5 proof
  does not need a vague "more suppression" statement: it must either derive
  S_T <= 3.876331e-6 from the triplet-filter/missing-partner sector, or derive
  an equivalent common Knu amplitude scale <= 0.516844.  The latter could come
  from a channel-specific flavor alignment, a dressing zero, or a principled
  chiral/lattice contraction that reduces LLLL_uud,Knu.

paper update:
  Added the "Knu target map" proposition to paper/gut_framework.tex.
  Updated final ledger text from 31 to 32 claims and passing/conditional
  checks from 22 to 23.

conditional theorem ledger update:
  Added row:
    "Knu target map localizes the remaining d=5 future-stress requirement"
  with status PASS_CONDITIONAL.

  New counts:
    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 19
    TUNED_FALLBACK: 3

current obstacle:
  Exact O2 remains open.  For a publication-level conditional EFT paper, the
  active technical bottlenecks remain:
    1. full CKM/PMNS/flavor fit;
    2. full channel-specific d=5 proton decay, now sharply led by
       LLLL_upupdown_Knu;
    3. supplemental reproducibility.
  For a first-principles claim, the additional hard blocker remains the
  microscopic origin of A5-A6 constrained/composite source dynamics.

next attempted nontrivial idea:
  Start the true Knu-centered d=5 calculation:
    1. build explicit LLLL uud Knu Wilson tensors in the mass basis using the
       same flavor matrices that will enter the full flavor fit;
    2. separate short-distance RGE, chargino/neutralino dressing, and
       chiral/lattice matrix elements instead of collapsing them into one
       conservative width factor;
    3. test whether physical neutrino-flavor summation or CKM/PMNS alignment
       gives the required amplitude scale <= 0.516844;
    4. if not, push the triplet-filter sector toward S_T <= 3.876331e-6 and
       check whether the RGE/proton benchmark still survives.

verification completed:
  py_compile passed for:
    code/audit_knu_target_map.py
    code/audit_conditional_theorem_ledger.py
  audit_knu_target_map.py generated output/knu_target_map/*.
  audit_conditional_theorem_ledger.py regenerated the 32-claim ledger.
  pdflatex compiled paper/gut_framework.tex twice after the final TeX edit.
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries after the final compile.

Roadmap note 2026-05-08: mock modular / shadow route
====================================================

decision:
  Worth adding to the roadmap, but only as a secondary exploratory mechanism,
  not as a replacement for the current publication blockers.  It should be
  framed as a possible mathematical origin for PSLT layer counting and for the
  recurring contact/shadow completions, not as a direct proof of Spin(10),
  hypercharge, seesaw, threshold matching, or proton safety.

reason:
  The current paper already contains a CP1 O(2) transvectant/contact sector,
  constrained/composite 54/210 source sectors, and several "visible piece plus
  completion" patterns.  A mock modular or mock Jacobi structure could unify
  these at the PSLT layer-counting level:

    visible holomorphic layer  -> finite-N layer coefficients;
    shadow/completion sector   -> hidden/contact/source compensation;
    wall-crossing sensitivity  -> layer leakage or finite-N correction.

  This is mathematically more disciplined than an arbitrary Cardy or Gaussian
  degeneracy ansatz, provided the polar data and shadow are fixed rather than
  fitted freely.

not allowed to do:
  Do not use mock modularity to modify:
    G_GUT = Spin(10),
    Y = T3R + (B-L)/2,
    anomaly cancellation,
    type-I seesaw algebra,
    the explicit d=5 Wilson/dressing calculation.

  It may motivate a hidden/completion sector, but the action-level GUT checks
  still have to pass independently.

candidate mathematical object:
  Define a vector-valued mock Jacobi layer function

    Zhat_lay(tau,z;E_d) = Z_lay^+(tau,z;E_d) + Z_lay^-(tau,z;E_d),

    Z_lay^+ = sum_{N,l} c^+_{N,l}(E_d) q^{N-Delta} y^l,
    c^+_{N,l}(E_d) in C^3,

  with modular completion

    Zhat_lay |_{k,m} gamma = rho(gamma) Zhat_lay,
    xi_{k,m}(Zhat_lay) = Theta_shadow.

  Then define a PSLT-visible finite-N degeneracy estimator

    g_N^mock(E_d) = sum_l || c^+_{N,l}(E_d) ||_W

  and insert it only into the already existing occupancy functional

    P_N =
      B_N g_N^mock (1-exp(-Gamma_N t))
      / sum_K B_K g_K^mock (1-exp(-Gamma_K t)).

acceptance tests:
  1. finite-N layer stability:
       Compare Cardy, phase-space, and mock g_N choices.  The first three
       visible layers should remain concentrated without inverse-Yukawa or
       manually boosted weights.

  2. parameter discipline:
       The mock object must be fixed by a small amount of polar data plus a
       specified shadow.  If every coefficient is fitted, the route fails.

  3. index versus degeneracy:
       If c_N is a signed index, keep it distinct from a positive degeneracy.
       Do not silently set g_N = |c_N| without a physical estimator argument.

  4. contact-sector check:
       Test whether the shadow/contact interpretation naturally reproduces
       the already derived CP1 transvectant S_perp rather than inventing a new
       unrelated flavor deformation.

  5. proton neutrality:
       The mock route must not be used to claim proton safety until the full
       Knu-centered d=5 calculation is done.  At most it can suggest a
       symmetry/completion origin for a stronger triplet filter.

priority:
  Keep this behind the active hard blockers:
    A. full CKM/PMNS/flavor fit;
    B. full channel-specific d=5 proton decay, led by LLLL_upupdown_Knu;
    C. supplemental reproducibility;
    D. microscopic origin of A5-A6 if claiming first-principles completion.

next attempted use:
  After the Knu-centered d=5 audit has an explicit flavor basis, run a small
  PSLT-only mock-layer experiment:
    code/audit_mock_layer_degeneracy.py
  comparing g_N^Cardy, g_N^phase, and a fixed toy vector-valued mock/Jacobi
  coefficient model.  The output should report R3 = P1+P2+P3 stability and
  whether the mock/shadow route reduces finite-N arbitrariness without adding
  new tunable freedom.

Heartbeat update 2026-05-08 04:44 UTC / 12:44 Asia/Taipei
==========================================================

status:
  The full first-principles GUT task is still incomplete.  The new local result
  is narrower and useful: the Knu target-map bound and the older mass-basis
  d=5 Wilson tensor table have now been compared on a common normalization
  axis.  This shows that the bare flavor/Wilson tensor rotation is not the
  present limiting obstruction by itself.  The remaining d=5 obstruction is
  the bridge from the same Wilson tensor into a full channel-specific
  chargino/neutralino dressing, soft-spectrum, short-distance running, and
  chiral/lattice width calculation.

new audit:
  Added:
    code/audit_knu_wilson_normalization_gap.py

  Generated:
    output/knu_wilson_normalization_gap/summary.json
    output/knu_wilson_normalization_gap/knu_wilson_normalization_gap.csv
    output/knu_wilson_normalization_gap/report.md

mathematical replay:
  Let S_T,A^max be the maximum triplet filter allowed by the older
  mass-basis Wilson tensor audit under its common 100 TeV dressing
  normalization.  Compare it with:

    S_T^current = 7.5e-6,
    S_T^target(10^35 yr) = 3.876330873e-6.

  Define

    H_A^common = S_T,A^max / S_T^current,
    G_A = S_T,A^max / S_T^target(10^35 yr).

  The worst common-dressing row is:

    hypothesis = transvectant_contact_QQ,
    channel = LLLL_upupdown_Knu,
    S_T,A^max = 1.845187e-4,
    H_A^common = 24.60249,
    lifetime headroom = (H_A^common)^2 = 605.2826,
    G_A = 47.60138.

interpretation:
  The old mass-basis Wilson tensor table is optimistic relative to the
  conservative Knu target-map replay by an amplitude factor about 47.6 in its
  worst row.  Therefore the next calculation must not merely rescan bare
  flavor rotations.  It must use one Wilson tensor and carry it through the
  complete d=5 pipeline before deciding whether S_T = 7.5e-6, or the stronger
  S_T <= 3.876331e-6 future target, is really required.

paper and ledger sync:
  Added a "Knu Wilson normalization gap" proposition to paper/gut_framework.tex.
  Updated the conditional theorem ledger from 32 to 33 claims.  The new ledger
  counts are:

    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 20
    TUNED_FALLBACK: 3

current obstacle:
  The active hard blockers remain:
    A. full CKM/PMNS/flavor fit;
    B. full channel-specific d=5 proton decay, especially LLLL_upupdown_Knu;
    C. supplemental reproducibility package;
    D. microscopic first-principles origin of A5-A6 if the paper claims more
       than a conditional EFT branch.

next attempted nontrivial idea:
  Build a unified Knu-centered d=5 pipeline:

    code/audit_full_knu_channel_pipeline.py

  Required inputs:
    1. the same transvectant/contact mass-basis Wilson tensors used in
       output/dimension5_wilson_tensors;
    2. soft-spectrum eigenstates and chargino/neutralino dressing factors;
    3. short-distance and long-distance running factors;
    4. chiral/lattice hadronic matrix element envelope;
    5. neutrino-flavor summation and CKM/PMNS rotations in the same basis as
       the full flavor fit.

  Acceptance tests:
    1. reproduce the older common-dressing Wilson-table lifetime in the
       degenerate-soft, single-K-factor limit;
    2. reproduce the conservative Knu target-map amplitude requirement when
       collapsed to its envelope approximation;
    3. report the physical required suppression as either
         S_T <= 3.876331e-6
       or
         additional Knu amplitude scale <= 0.516844,
       with the bottleneck channel explicitly identified.

verification completed:
  py_compile passed for:
    code/audit_knu_wilson_normalization_gap.py
    code/audit_conditional_theorem_ledger.py

  audit_knu_wilson_normalization_gap.py generated its JSON/CSV/Markdown
  outputs.  audit_conditional_theorem_ledger.py regenerated the 33-claim
  ledger.  paper/gut_framework.tex compiled twice after the final TeX edit,
  and the final warning scan found no LaTeX Warning, Package Warning, Error,
  Undefined, or Overfull entries.

Heartbeat update 2026-05-08 04:50 UTC / 12:50 Asia/Taipei
==========================================================

status:
  The task is still incomplete, but the d=5 bottleneck is now more sharply
  localized.  I added a one-axis Knu pipeline bridge audit that places the
  existing mass-basis Wilson, field-basis Wilson, soft-spectrum dressing, and
  final Knu target-map replays on the same uniform 1e35 yr stress axis.

new audit:
  Added:
    code/audit_full_knu_channel_pipeline.py

  Generated:
    output/full_knu_channel_pipeline/summary.json
    output/full_knu_channel_pipeline/pipeline_bridge.csv
    output/full_knu_channel_pipeline/report.md

mathematical replay:
  Define the future-stress margin

    M_1e35 = M_current * (2.4e34 / 1.0e35) = 0.24 M_current

  for Knu-type rows whose local current bound is normalized at 2.4e34 yr, and
  define amplitude headroom

    h = sqrt(M_1e35).

  The audited headroom ladder is:

    mass_basis_common_wilson:
      channel = LLLL_upupdown_Knu
      M_1e35 = 6.052826e2
      h = 24.60249

    field_basis_common_dressing_max_width:
      channel = RRRR_uusd_anycharged
      M_1e35 = 5.931887e1
      h = 7.701874
      loss from previous = 3.194351

    soft_spectrum_preferred_max_width:
      channel = LLLL_upupdown_Knu
      M_1e35 = 3.100993e-1
      h = 0.5568656
      loss from previous = 13.83076

    epi_projected_knu_target:
      channel = LLLL_upupdown_Knu
      M_1e35 = 2.671278e-1
      h = 0.5168441
      loss from previous = 1.077434

key conclusion:
  The total future-headroom loss from common Wilson normalization to the final
  projected Knu target map is 47.60138, matching the previous normalization-gap
  factor.  The largest single contraction is not the final projection step; it
  is the move from common field-basis dressing to the explicit preferred soft
  spectrum:

    7.701874 / 0.5568656 = 13.83076.

  Therefore the next proton-side work should focus on the soft-spectrum
  dressing kernel and its normalization to short-distance/chiral/lattice width,
  not on another bare flavor tensor rescan.

paper and ledger sync:
  Added a "Knu pipeline headroom ladder" proposition to paper/gut_framework.tex.
  Added a corresponding conditional ledger row:

    Full Knu pipeline bridge localizes the d=5 headroom loss.

  The conditional theorem ledger now has 34 claims:

    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 21
    TUNED_FALLBACK: 3

current obstacle:
  Full channel-specific d=5 proton decay is still open.  The next missing
  object is a single calculation that starts from the same Wilson tensor and
  computes the physical LLLL_upupdown_Knu width with chargino/neutralino
  dressing matrices, soft eigenstates, short-distance running, neutrino-flavor
  summation, and chiral/lattice matrix elements without changing
  normalization conventions midstream.

next attempted nontrivial idea:
  Build the first true physical Knu-width module, not just a bridge:

    code/audit_full_knu_width.py

  It should satisfy these replay checks before being trusted:
    1. in the degenerate 100 TeV dressing limit, reproduce the common Wilson
       table's Knu row;
    2. with the audited soft-spectrum point, reproduce the soft-spectrum
       preferred max-width margin 0.3100993 at 1e35 yr;
    3. with the e-pi/identical-up projection profile, reproduce the final Knu
       target margin 0.2671278 and amplitude headroom 0.5168441;
    4. report whether a structural suppression can replace the fallback
       requirement S_T <= 3.876331e-6.

verification completed:
  py_compile passed for:
    code/audit_full_knu_channel_pipeline.py
    code/audit_conditional_theorem_ledger.py

  audit_full_knu_channel_pipeline.py generated its JSON/CSV/Markdown outputs.
  audit_conditional_theorem_ledger.py regenerated the 34-claim ledger.  The
  TeX draft compiled twice after adding the new proposition and ledger count,
  and the final warning scan found no LaTeX Warning, Package Warning, Error,
  Undefined, or Overfull entries.

Heartbeat update 2026-05-08 04:56 UTC / 12:56 Asia/Taipei
==========================================================

status:
  The task is still incomplete, but one important d=5 ambiguity has been
  removed.  I built the first calibrated Knu-width audit and verified that the
  preferred soft-spectrum bottleneck and the final projected Knu target row
  use the same effective width normalization.

new audit:
  Added:
    code/audit_full_knu_width.py

  Generated:
    output/full_knu_width/summary.json
    output/full_knu_width/knu_width_calibration.csv
    output/full_knu_width/report.md

mathematical replay:
  The audited width formula is

    Gamma(Knu) = K_dyn K_had (S_T |A_dress|)^2.

  For the preferred max-width Knu bottleneck:

    soft_spectrum_width:
      pair = 23
      |A_dress| = 2.361106e-11
      M_1e35 = 3.100993e-1
      K_dyn = 5.377683e-1 yr^-1
      S_T^max(1e35 yr) = 4.176492e-6

    projected_worst_phase_width:
      pair = 13
      triplet_phase = 4.188790
      |A_dress| = 2.543937e-11
      M_1e35 = 2.671278e-1
      K_dyn = 5.377683e-1 yr^-1
      S_T^max(1e35 yr) = 3.876331e-6

  Cross-check:

    M_projected / M_soft = 0.8614268,
    (A_soft / A_projected)^2 = 0.8614268,
    mismatch = 1.110223e-16,
    max relative K_dyn spread = 0.

key conclusion:
  The remaining d=5 Knu future-stress gap is not a hidden normalization
  mismatch between local scripts.  It is a physical requirement inside a
  single calibrated width formula: either an extra Knu amplitude suppression
  of 1.934819 or a stronger triplet filter

    S_T <= 3.876331e-6.

  The hard next task is therefore no longer "align the replay conventions";
  it is "derive a structural amplitude suppression in the same calibrated
  width formula" or accept the stronger triplet filter as a tuned fallback.

paper and ledger sync:
  Added a "Calibrated Knu width formula" proposition to
  paper/gut_framework.tex.

  Added a corresponding conditional ledger row:

    Calibrated Knu width formula aligns soft and projected target rows.

  The conditional theorem ledger now has 35 claims:

    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 22
    TUNED_FALLBACK: 3

current obstacle:
  The d=5 obstruction is now a true model-building obstruction: find a
  symmetry, triplet-mixing near-null, phase structure, or soft-spectrum
  alignment that gives a real Knu amplitude factor <= 0.516844 inside the
  calibrated width formula.  If no such structural suppression exists, the
  branch must use S_T <= 3.876331e-6 for a 1e35 yr future stress and then
  rerun threshold/RGE/proton consistency.

next attempted nontrivial idea:
  Convert the triplet-mixing near-null direction into the calibrated width
  formula rather than treating it as a separate Wilson proxy.  Proposed script:

    code/audit_knu_nearnull_width.py

  It should:
    1. read the near-null W_AB direction from the existing triplet-mixing
       nullspace audit;
    2. compute A_dress for LLLL_upupdown_Knu after the same soft-spectrum
       dressing used above;
    3. insert it into Gamma = K_dyn K_had (S_T |A_dress|)^2;
    4. test whether the required amplitude scale 0.516844 is achieved without
       rank-deficient or proton-dangerous side effects in K0mu and RRRR.

verification completed:
  py_compile passed for:
    code/audit_full_knu_width.py
    code/audit_conditional_theorem_ledger.py

  audit_full_knu_width.py generated JSON/CSV/Markdown outputs.  The conditional
  theorem ledger regenerated with 35 claims.  The TeX draft compiled twice
  after the new proposition and ledger update, and the final warning scan found
  no LaTeX Warning, Package Warning, Error, Undefined, or Overfull entries.

Heartbeat update 2026-05-08 05:04 UTC / 13:04 Asia/Taipei
==========================================================

status:
  The task is still incomplete, but the calibrated Knu bottleneck now has a
  concrete candidate solution mechanism.  I inserted the triplet-mixing
  near-null direction into the calibrated Knu width formula by a channel-wise
  Wilson-ratio transfer.

new audit:
  Added:
    code/audit_knu_nearnull_width.py

  Generated:
    output/knu_nearnull_width/summary.json
    output/knu_nearnull_width/knu_nearnull_width.csv
    output/knu_nearnull_width/report.md

mathematical replay:
  Use the transfer rule

    A_near_dress(c) =
      A_identity_dress(c) |C_near(c)| / |C_identity(c)|,

  then insert into the already calibrated formula

    Gamma(c) = K_dyn K_had(c) (S_T |A_near_dress(c)|)^2.

  For S_T = 7.5e-6 and the conservative max-width envelope:

    LLLL_upupdown_Knu:
      |C_near|/|C_id| = 3.182422e-5
      |A_near_dress| = 8.095881e-16
      M_1e35 = 2.637567e8
      S_T^max(1e35 yr) = 1.218044e-1

    LLLL_downdownup_Knu:
      |C_near|/|C_id| = 2.786993e-5
      |A_near_dress| = 6.462579e-16
      M_1e35 = 4.139234e8
      S_T^max(1e35 yr) = 1.525883e-1

    LLLL_upupdown_K0mu:
      |C_near|/|C_id| = 1.391926e-5
      |A_near_dress| = 2.764838e-16
      M_1e35 = 1.071334e10
      S_T^max(1e35 yr) = 7.762894e-1

    RRRR_uusd_anycharged:
      |C_near|/|C_id| = 4.984669e-2
      |A_near_dress| = 3.349118e-13
      M_1e35 = 1.541243e3
      S_T^max(1e35 yr) = 2.944400e-4

key conclusion:
  The near-null direction easily beats the required calibrated target

    h_required = 0.5168441.

  The worst future channel becomes RRRR_uusd_anycharged with margin
  1.541243e3, so the previous Knu amplitude suppression requirement 1.934819
  is solved at the proxy level.  However, the same audit shows:

    rank(W) = 3,
    condition(W) = infinity,
    |<W,W_id>| = 0.479435.

  Therefore this is not a final proton-decay proof.  It is a strong conditional
  branch: proton safety follows if an action-level triplet sector can generate
  this codimension-one inverse propagator, or a finite-rank lift that preserves
  the displayed suppression.

paper and ledger sync:
  Added a "Near-null inserted into the calibrated width" proposition to
  paper/gut_framework.tex.

  Added a corresponding conditional ledger row:

    Triplet near-null solves the calibrated Knu width target as a conditional branch.

  The conditional theorem ledger now has 36 claims:

    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 23
    TUNED_FALLBACK: 3

current obstacle:
  The d=5 obstruction has moved from "find an O(1) Knu suppression" to
  "derive the near-null triplet inverse propagator from a consistent
  action-level mechanism."  The most direct obstruction is that the audited
  near-null has rank 3 and infinite condition number.

next attempted nontrivial idea:
  Return to the rank-one lift / triplet-sector superpotential and test whether
  a finite epsilon lift preserves the calibrated width margins without
  generating super-Planckian triplet eigenvalues or dangerous threshold
  splittings.  Proposed script:

    code/audit_nearnull_finite_lift_width.py

  It should scan epsilon in W_epsilon, compute:
    1. min singular value and triplet mass eigenvalue range;
    2. calibrated width margins for LLLL_upupdown_Knu, LLLL_downdownup_Knu,
       K0mu, and RRRR;
    3. threshold beta-vector leakage from the lifted triplet sector;
    4. whether there is an epsilon window with
         M_triplet < M_Pl,
         M_1e35 > 1 for all channels,
         and projected threshold residual below the existing R=200 tolerance.

verification completed:
  py_compile passed for:
    code/audit_knu_nearnull_width.py
    code/audit_conditional_theorem_ledger.py

  audit_knu_nearnull_width.py generated JSON/CSV/Markdown outputs.  The
  conditional theorem ledger regenerated with 36 claims.  The TeX draft
  compiled twice after the new proposition and ledger update, and the final
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries.

Heartbeat update 2026-05-08 05:41 UTC / 13:41 Asia/Taipei
==========================================================

status:
  The task is still incomplete.  The finite-lift audit now shows that the
  rank-deficient near-null can be regularized into Planck-safe, proton-safe
  finite cards, but only if the triplet-lift sector is also threshold-locked.
  Without such locking, the mass spread gives an unacceptably large
  conservative threshold proxy.

new audit:
  Added:
    code/audit_nearnull_finite_lift_width.py

  Generated:
    output/nearnull_finite_lift_width/summary.json
    output/nearnull_finite_lift_width/lift_summary.csv
    output/nearnull_finite_lift_width/channel_widths.csv
    output/nearnull_finite_lift_width/report.md

mathematical replay:
  For every finite-lift card, use the same calibrated transfer rule as the
  near-null audit:

    A_lift_dress(c) =
      A_identity_dress(c) |C_lift(c)| / |C_identity(c)|,

    Gamma(c) = K_dyn K_had(c) (S_T |A_lift_dress(c)|)^2.

  The finite-lift card is then checked against:

    M_1e35(c) = tau(c) / 1e35 yr > 1,
    M_max < Mbar_Pl = 2.435e18 GeV,
    Delta_split_proxy = log(M_max/M_min)/(2*pi)
      < ||P Delta||_R200 = 5.022739e-4.

  The threshold proxy is deliberately conservative: it assumes a single
  uncancelled unit beta-vector.  A complete-multiplet or locked-triplet lift
  can cancel it, but an arbitrary split triplet sector cannot.

numerical result:
  For S_T = 7.5e-6 and M_min = 1e16 GeV:

    condition_cap_10:
      Mmax/Mmin = 10
      Mmax = 1.0e17 GeV
      worst future margin = 3.222499
      Planck safe = yes
      Delta_split_proxy = 3.664678e-1

    condition_cap_30:
      Mmax/Mmin = 30
      Mmax = 3.0e17 GeV
      worst future margin = 2.918314e1
      Planck safe = yes
      Delta_split_proxy = 5.413174e-1

    condition_cap_100:
      Mmax/Mmin = 100
      Mmax = 1.0e18 GeV
      worst future margin = 3.313948e2
      Planck safe = yes
      Delta_split_proxy = 7.329356e-1

    condition_cap_300:
      Mmax/Mmin = 300
      Mmax = 3.0e18 GeV
      worst future margin = 1.608078e3
      Planck safe = no
      Delta_split_proxy = 9.077852e-1

key conclusion:
  Three audited finite-lift cards are both proton-safe at the 1e35 yr stress
  target and below the reduced-Planck benchmark:

    condition_cap_10, condition_cap_30, condition_cap_100.

  The best Planck/proton row is:

    condition_cap_100:
      worst future margin = 3.313948e2,
      Mmax = 1.0e18 GeV.

  The lowest-hierarchy Planck/proton row is:

    condition_cap_10:
      worst future margin = 3.222499,
      Mmax/Mmin = 10.

  However, no Planck/proton-safe row passes the naive uncancelled threshold
  proxy.  Therefore the next real obstruction is not the width formula; it is
  the action-level construction of a threshold-degenerate complete-multiplet or
  Goldstone/projector-locked triplet lift.

paper and ledger sync:
  Added "Finite near-null lift in the calibrated width formula" to
  paper/gut_framework.tex.

  Added a corresponding conditional ledger row:

    Finite triplet near-null lifts can be proton-safe and Planck-safe only with
    threshold locking.

  The conditional theorem ledger now has 37 claims:

    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 24
    TUNED_FALLBACK: 3

current obstacle:
  The finite near-null lift removes the rank-deficiency objection at the proxy
  level, but a generic split triplet mass spectrum produces a projected
  threshold proxy much larger than the R=200 tolerance.  The lift must be
  generated without leaving an uncancelled non-universal heavy threshold.

next attempted nontrivial idea:
  Build a threshold-locked Schur/clockwork triplet-lift sector in which the
  effective inverse propagator is near-null, but the propagating eigenstates
  remain in complete or pairwise-degenerate multiplets.  The creative route is:

    1. introduce degenerate mediator pairs whose Schur complement generates the
       condition-capped inverse propagator;
    2. keep the physical mediator spectrum complete-multiplet degenerate so the
       non-universal beta-vector cancels;
    3. allow the near-null to live in the effective inverse propagator rather
       than in large physical mass splittings.

  Proposed next script:

    code/construct_threshold_locked_triplet_lift.py

verification plan:
  The next audit should compute:

    1. full component Hessian eigenvalues of the locked triplet/mediator sector;
    2. effective inverse propagator W_eff and its distance to condition_cap_10,
       condition_cap_30, and condition_cap_100;
    3. exact projected threshold vector using the real fragment spectrum;
    4. calibrated d=5 width margins for Knu, K0mu, and RRRR;
    5. whether at least one row satisfies
         Mmax < 2.435e18 GeV,
         min M_1e35 > 1,
         ||P Delta||_2 < 5.022739e-4.

secondary roadmap note:
  The mock modular / mock Jacobi idea is worth keeping, but as a secondary
  PSLT-layer item rather than the current hard blocker.  It may organize the
  layer-counting function g_N and explain shadow/contact sectors, but it does
  not replace the action-level flavor fit or the locked d=5 proton-decay
  calculation.  Add it later as a controlled "mock-completed layer degeneracy"
  audit only after the triplet-lift threshold issue is resolved.

verification completed:
  py_compile passed for:
    code/audit_nearnull_finite_lift_width.py
    code/audit_conditional_theorem_ledger.py

  audit_nearnull_finite_lift_width.py generated JSON/CSV/Markdown outputs.  The
  conditional theorem ledger regenerated with 37 claims.  The TeX draft
  compiled twice after the new proposition and ledger update, and the final
  warning scan found no LaTeX Warning, Package Warning, Error, Undefined, or
  Overfull entries.

Heartbeat update 2026-05-08 11:39 UTC / 19:39 Asia/Taipei
==========================================================

status:
  The task is still incomplete, but the previous threshold-locking bottleneck
  has a concrete algebraic solution.  A Julia/Sz.-Nagy unitary dilation can
  realize the near-null inverse triplet propagator as the visible subblock of
  the inverse of a larger exactly degenerate heavy mass matrix.

new audit:
  Added:
    code/construct_threshold_locked_triplet_lift.py

  Generated:
    output/threshold_locked_triplet_lift/summary.json
    output/threshold_locked_triplet_lift/locked_lift_summary.csv
    output/threshold_locked_triplet_lift/report.md

mathematical construction:
  For any target inverse propagator W with ||W||_2 <= 1, define

    D_L = (I - W W^\dagger)^{1/2},
    D_R = (I - W^\dagger W)^{1/2},

  and the Julia/Sz.-Nagy dilation

    U(W) = [[W, D_L],
            [D_R, -W^\dagger]].

  Then U(W) is unitary.  If

    M_full = M_lock U(W)^\dagger,

  then all physical singular values of M_full are exactly M_lock, while

    P M_full^{-1} P = W / M_lock.

  Thus the near-null can be an interference subblock of a larger inverse mass
  matrix, rather than a literal split physical spectrum.

numerical result:
  The audit used

    M_lock = M_Sigma8 = 3.887852e15 GeV

  from the current R=200 benchmark.  It checked the exact rank-three near-null
  and the condition-capped rows:

    rank3_limit:
      ||W||_2 = 9.999965e-1
      ||U^\dagger U - I||_2 = 3.346e-14
      full singular spread / M_lock = 3.364e-14
      worst future margin = 1.541243e3
      projected threshold = 0

    condition_cap_10:
      ||U^\dagger U - I||_2 = 5.303e-15
      full singular spread / M_lock = 3.331e-16
      worst future margin = 3.222499
      projected threshold = 0

    condition_cap_30:
      ||U^\dagger U - I||_2 = 3.710e-14
      full singular spread / M_lock = 3.697e-14
      worst future margin = 2.918314e1
      projected threshold = 0

    condition_cap_100:
      ||U^\dagger U - I||_2 = 5.762e-14
      full singular spread / M_lock = 5.729e-14
      worst future margin = 3.313948e2
      projected threshold = 0

key conclusion:
  All four rows pass the conditional algebraic checks:

    1. calibrated d=5 future margin > 1,
    2. M_lock below the reduced Planck benchmark,
    3. exact degenerate physical singular spectrum,
    4. projected threshold zero under complete-degenerate multiplet embedding.

  This removes the rank-deficiency and split-threshold obstruction at the
  algebraic EFT level.  The exact rank-three near-null is now viable as a
  visible inverse-propagator subblock of an 8x8 degenerate triplet/mediator mass
  matrix.

paper and ledger sync:
  Added "Threshold-locked unitary dilation of the triplet lift" to
  paper/gut_framework.tex.

  Added a corresponding conditional ledger row:

    Threshold-locked unitary dilation removes the triplet rank and
    split-threshold obstruction algebraically.

  The conditional theorem ledger now has 38 claims:

    FAIL: 1
    NO_GO: 2
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 25
    TUNED_FALLBACK: 3

current obstacle:
  The near-null triplet sector now has an algebraic threshold-locked EFT
  realization, but it is not yet derived from a Spin(10) representation table,
  charge grading, or driving-sector superpotential.  The construction requires
  four mediator pairs and a unitary 8x8 block; the next hard problem is to show
  that this block can be enforced by allowed Spin(10)-invariant renormalizable
  couplings or by a controlled constrained/composite sector.

next attempted nontrivial idea:
  Build a "unitary-link" mediator sector that enforces the Julia dilation by
  F-term orthonormality constraints.  The proposed creative route is:

    1. introduce four visible triplet sources T_A and four mediator copies X_A;
    2. promote the 8x8 unitary matrix U to a constrained link field L;
    3. use driving fields to impose L^\dagger L = I and fix the visible
       subblock P L P = W_target;
    4. embed the T/X fields into complete degenerate Spin(10) multiplets so
       the threshold vector remains zero;
    5. check whether the required constraints can be written as renormalizable
       singlet/adjoint/54/210 contractions or must remain a constrained
       sigma-model ansatz.

  Proposed next script:

    code/audit_unitary_link_spin10_embedding.py

verification plan:
  The next audit should compute:

    1. representation content and multiplicity of the visible and mediator
       triplet copies;
    2. all allowed quadratic and cubic Spin(10)-invariant couplings involving
       the proposed unitary link/drivers;
    3. F-term residuals for L^\dagger L = I and P L P = W_target;
    4. non-singlet driver/source spectrum and whether it forms complete
       degenerate multiplets;
    5. exact projected threshold vector and calibrated d=5 width margins after
       including the link/driver spectrum.

remaining global blockers:
  The first-principles GUT is still not complete.  The open blockers remain:

    O1. full joint CKM/PMNS/flavor fit;
    O2. complete field-basis d=5 proton decay calculation with physical soft
        spectrum, flavor rotations, and lattice/chiral inputs;
    O3. microscopic first-principles origin of the constrained source/link
        sectors.

verification completed:
  py_compile passed for:
    code/construct_threshold_locked_triplet_lift.py
    code/audit_conditional_theorem_ledger.py

  construct_threshold_locked_triplet_lift.py generated JSON/CSV/Markdown
  outputs.  The conditional theorem ledger regenerated with 38 claims.  The
  TeX draft compiled twice after the new proposition and ledger update, and
  the final warning scan found no LaTeX Warning, Package Warning, Error,
  Undefined, or Overfull entries.

Heartbeat update 2026-05-08 11:45 UTC / 19:45 Asia/Taipei
==========================================================

status:
  The task is still incomplete.  The unitary dilation is a valid
  threshold-locked EFT block, but the newest audit shows that it is not yet an
  ordinary perturbative elementary Spin(10) completion.

new audit:
  Added:
    code/audit_unitary_link_spin10_embedding.py

  Generated:
    output/unitary_link_spin10_embedding/summary.json
    output/unitary_link_spin10_embedding/matrix_obstructions.csv
    output/unitary_link_spin10_embedding/scenario_scorecard.csv
    output/unitary_link_spin10_embedding/report.md

mathematical checks:
  A single species of real Spin(10) 10-plets has a superpotential mass

    W = 10_a M_ab 10_b,

  so the copy-space mass matrix must be symmetric.  The audited unitary
  dilation blocks fail this:

    rank3_limit:
      ||U - U^T||_F / ||U||_F = 3.007790e-1
      ||U^T U - I||_2 = 8.149241e-1

    condition_cap_10:
      ||U - U^T||_F / ||U||_F = 3.298594e-1
      ||U^T U - I||_2 = 8.204971e-1

    condition_cap_30:
      ||U - U^T||_F / ||U||_F = 3.043147e-1
      ||U^T U - I||_2 = 8.146548e-1

    condition_cap_100:
      ||U - U^T||_F / ||U||_F = 3.011472e-1
      ||U^T U - I||_2 = 8.146781e-1

  Therefore:

    1. single-species complete 10 mass matrices cannot realize the audited U;
    2. holomorphic F-term orthogonality L^T L = I is also not the same as
       unitarity L^\dagger L = I;
    3. exact unitarity must be a Kahler/D-term/constrained-link datum, or a
       composite/NLSM condition.

representation scorecard:
  Using the current R=200 benchmark and baseline constrained-source branch:

    single_species_8_complete_10s:
      extra Dynkin T = 8
      alpha_G^-1(200 M_G) = -5.015855
      projected threshold = 0
      status = FAIL_MATRIX_AND_R200

    two_species_16_complete_10s:
      extra Dynkin T = 16
      alpha_G^-1(200 M_G) = -11.761882
      projected threshold = 0
      status = FAIL_R200

    post_spin10_8_su5_complete_5pairs:
      extra Spin(10) T = 0 in the post-breaking EFT bookkeeping
      alpha_G^-1(200 M_G) = 1.730172
      projected threshold = 1.922963e-16
      status = PASS_CONDITIONAL_POST_GUT_EFT

    constrained_composite_unitary_link:
      extra Spin(10) T = 0 below the compositeness scale
      alpha_G^-1(200 M_G) = 1.730172
      projected threshold = 0
      status = PASS_CONDITIONAL_CONSTRAINED

    triplet_only_8_pairs:
      alpha_G^-1(200 M_G) = 1.730172
      projected threshold = 5.431815e-1
      status = FAIL_THRESHOLD

key conclusion:
  The naive elementary Spin(10) 10-copy completion is ruled out for the current
  R=200 branch:

    single complete-10 species fails the symmetric mass requirement;
    two complete-10 species allow arbitrary U but destroy perturbative reach;
    triplet-only remnants fail the threshold bound.

  The surviving interpretations are explicitly conditional:

    A. post-Spin(10) SU(5)-complete EFT link;
    B. constrained/composite unitary-link sector.

paper and ledger sync:
  Added "Unitary-link Spin(10) embedding scorecard" to paper/gut_framework.tex.

  Added a corresponding no-go ledger row:

    Perturbative elementary Spin(10) 10-copy unitary link closes the triplet
    sector.

  The conditional theorem ledger now has 39 claims:

    FAIL: 1
    NO_GO: 3
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 25
    TUNED_FALLBACK: 3

current obstacle:
  The triplet near-null can be made threshold-locked algebraically, but not by
  the simplest elementary perturbative Spin(10) 10-copy sector.  The next hard
  problem is to choose between a post-Spin(10) SU(5)-complete EFT link and a
  constrained/composite unitary-link origin, then propagate that choice into
  the paper's assumptions A5-A7.

next attempted nontrivial idea:
  Build a constrained/composite unitary-link action analogous to the existing
  shared-U source NLSM.  The idea is to treat the 8x8 unitary link as a
  Stiefel/unitary sigma-model coordinate with:

    L^\dagger L = I,
    P L P = W_target,

  and to couple it to complete post-GUT SU(5) or constrained Spin(10)
  multiplet sources.  This would make the unitary condition a geometric
  constraint, not a holomorphic F-term.

  Proposed next script:

    code/audit_unitary_link_nlsm_action.py

verification plan:
  The next audit should compute:

    1. dimension of U(8), the fixed-subblock constraint rank, and residual
       moduli;
    2. whether the constrained manifold is nonempty for rank3_limit and
       condition_cap_10/30/100;
    3. tangent/normal decomposition and whether normal multipliers are gauge
       singlets or complete multiplets;
    4. threshold vector and UV perturbativity in the constrained bookkeeping;
    5. calibrated d=5 width margins after enforcing the exact constrained
       unitary link.

verification completed:
  py_compile passed for:
    code/audit_unitary_link_spin10_embedding.py
    code/audit_conditional_theorem_ledger.py

  audit_unitary_link_spin10_embedding.py generated JSON/CSV/Markdown outputs.
  The conditional theorem ledger regenerated with 39 claims.  The TeX draft
  compiled twice after the new no-go proposition and ledger update, and the
  final warning scan found no LaTeX Warning, Package Warning, Error, Undefined,
  or Overfull entries.

## 2026-05-08 21:01 Asia/Taipei heartbeat -- constrained unitary-link NLSM audit

status:
  The first-principles GUT derivation remains incomplete.  The present update
  closes the immediate "can the threshold-locked unitary link be a controlled
  constrained EFT target?" check, but it does not yet derive that target from
  microscopic Spin(10) dynamics.

new audit:
  Added:

    code/audit_unitary_link_nlsm_action.py

  Outputs:

    output/unitary_link_nlsm_action/summary.json
    output/unitary_link_nlsm_action/target_geometry.csv
    output/unitary_link_nlsm_action/scenario_scorecard.csv
    output/unitary_link_nlsm_action/report.md

mathematical construction:
  Treat the locked triplet link as a constrained unitary target rather than an
  elementary holomorphic mass matrix:

    L in U(8),        P L P = W_target.

  The tangent space is:

    T_L U(8) = { L K : K^\dagger = -K }.

  The fixed 4x4 block condition gives a real linear map:

    K -> (L K)_{4x4},

  from 64 real anti-Hermitian directions to 32 real block constraints.

numerical results:
  For all audited near-null targets:

    all_targets_nonempty = true
    all_fixed_block_constraints_full_rank = true
    fixed-block rank = 32
    residual real moduli = 32
    min contraction defect = 7.045509e-06

  The calibrated future d=5 proton-safety margins are:

    rank3_limit:        1.541243e+03
    condition_cap_10:   3.222499e+00
    condition_cap_30:   2.918314e+01
    condition_cap_100:  3.313948e+02

  Scenario scorecard:

    constrained_unitary_link:
      projected threshold = 0
      alpha_G^-1(200 M_G) = 1.730172
      min future margin = 3.222499
      status = PASS_CONDITIONAL

    post_spin10_su5_complete_8_fivepairs:
      projected threshold = 1.922963e-16
      alpha_G^-1(200 M_G) = 1.730172
      min future margin = 3.222499
      status = PASS_CONDITIONAL

interpretation:
  The constrained NLSM route is nonempty, full-rank, and threshold-silent in
  the preferred bookkeeping.  It provides a controlled conditional EFT origin
  for the unitary link, but it is not an unconditional microscopic Spin(10)
  derivation.  In particular, unitarity is being imposed as a Kahler/D-term or
  composite sigma-model datum, not as a holomorphic F-term consequence.

paper and ledger sync:
  Added "Constrained unitary-link NLSM action" to paper/gut_framework.tex.

  Added one corresponding ledger row:

    Constrained unitary-link NLSM gives a nonempty threshold-silent triplet-link
    origin.

  The conditional theorem ledger now has 40 claims:

    FAIL: 1
    NO_GO: 3
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 26
    TUNED_FALLBACK: 3

current obstacle:
  The A5-A6 unitary-link sector now has a mathematically controlled conditional
  EFT realization, but three publication-level blockers remain:

    O1. full joint CKM/PMNS/flavor fit;
    O2. full field-basis d=5 proton decay with physical flavor rotations and
        channel-specific hadronic inputs;
    O3. a microscopic first-principles Spin(10) origin for the constrained
        unitary-link/source sector.

next attempted nontrivial idea:
  Freeze the constrained unitary-link branch as the current conditional A5-A7
  completion and move the next numerical audit to the coupled flavor/proton
  bottleneck.  The useful next script should build a reproducible locked-link
  benchmark card that contains explicit Yukawa matrices, PMNS/CKM conventions,
  triplet couplings, soft dressing choices, and the channel-by-channel d=5
  Wilson tensors:

    code/audit_locked_link_full_flavor_d5_card.py

verification plan:
  The next audit should check:

    1. exact input matrices Y_u, Y_d, Y_e, Y_nu and M_R in one convention;
    2. quark/lepton masses, CKM, PMNS, and light-neutrino observables from the
       same matrices;
    3. C_{5L}^{ijkl}, C_{5R}^{ijkl} before and after flavor rotations;
    4. wino/higgsino dressing for p -> K^+ nu_bar, e^+ pi^0, mu^+ pi^0, and
       K^0 mu^+ channels;
    5. whether the constrained unitary link keeps the d=5 safety margin after
       the full flavor basis is used.

verification completed:
  py_compile passed for:

    code/audit_unitary_link_nlsm_action.py
    code/audit_conditional_theorem_ledger.py

  using a local PYTHONPYCACHEPREFIX to avoid writing outside the workspace.
  The TeX draft compiled twice after the new proposition and closure update.
  The final warning scan found no LaTeX Warning, Package Warning, Error,
  Undefined, or Overfull entries.

## 2026-05-08 21:17 Asia/Taipei heartbeat -- locked-link full flavor/d5 card

status:
  The first-principles GUT derivation remains incomplete.  This update turns
  the current full flavor + d=5 proton bottleneck into one reproducibility
  card, but it does not close the CKM fit or final Knu proton-stability proof.

new audit:
  Added:

    code/audit_locked_link_full_flavor_d5_card.py

  Outputs:

    output/locked_link_full_flavor_d5_card/locked_link_full_flavor_d5_card.json
    output/locked_link_full_flavor_d5_card/summary.json
    output/locked_link_full_flavor_d5_card/component_status.csv
    output/locked_link_full_flavor_d5_card/input_manifest.csv
    output/locked_link_full_flavor_d5_card/report.md

mathematical content:
  The card records the exact local matrices and the d=5 conventions on one
  axis.  The operator definitions are:

    C_5L^{abcd}
      = sum_AB (Y_QQ^A)_{ij} W_AB (Y_QL^B)_{kl}
        U_Q^{ia} U_Q^{jb} U_Q^{kc} U_L^{ld},

    C_5R^{abcd}
      = sum_AB (Y_UE^A)_{ij} W_AB (Y_UD^B)_{kl}
        U_u^{ia} U_e^{jd} U_u^{kb} U_d^{lc}.

  The Knu width normalization is:

    Gamma = K_dyn K_had (S_T |A_dress|)^2.

numerical results:
  Component status:

    exact_CP1_O2_yukawa_and_seesaw_card: PASS
    joint_CKM_PMNS_flavor_fit: OPEN
    locked_unitary_link_threshold: PASS_CONDITIONAL
    mass_and_field_basis_d5_replay: PASS_CONDITIONAL
    calibrated_Knu_width_formula: PASS_CONDITIONAL
    future_stress_d5_proton_safety: OPEN

  Key numbers:

    exact inputs available = true
    PMNS reconstructed = true
    CKM fit completed = false
    locked link threshold silent = true
    d5 current bound passes = true
    d5 future 1e35 yr stress passes = false
    future Knu amplitude suppression needed = 1.934819

  The card also writes SHA256 hashes and file sizes for all source JSON/CSV
  inputs, so later runs can tell whether a benchmark changed.

paper and ledger sync:
  Added "Locked-link flavor and d=5 reproducibility card" to
  paper/gut_framework.tex.

  Added one corresponding ledger row:

    Locked-link flavor and d=5 reproducibility card localizes the remaining
    bottleneck.

  The conditional theorem ledger now has 41 claims:

    FAIL: 1
    NO_GO: 3
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 27
    TUNED_FALLBACK: 3

current obstacle:
  The reproducibility layer is now much cleaner, but the two physics bottlenecks
  are unchanged:

    O1. CKM/full flavor is not fit by the current exact CP1 O(2) matrices.
    O2. The conservative Knu future-stress proton target needs an extra
        amplitude suppression factor 1.934819, or equivalently a stronger
        triplet filter S_T <= 3.876330873e-6.

  O3 also remains: the constrained unitary-link/source sector is conditional,
  not yet derived from microscopic first-principles Spin(10) dynamics.

next attempted nontrivial idea:
  Start a joint correction scan that attacks O1 and O2 together instead of
  separately.  The most promising next move is a constrained shadow/contact
  deformation:

    Y_a -> Y_a + epsilon_a S_perp + delta_a Delta_mock

  where S_perp is the already derived transvectant singlet/contact direction
  and Delta_mock is a small mock/shadow-compatible finite-N correction to the
  visible holomorphic Yukawa sector.  The scan should penalize CKM error and
  the Knu dressed amplitude simultaneously, while keeping mass ratios, PMNS,
  seesaw perturbativity, and locked-link threshold silence.

proposed next script:

    code/scan_shadow_contact_flavor_d5_joint.py

verification plan:
  The joint scan should output:

    1. deformed Y_u, Y_d, Y_e, Y_nu matrices and exact CKM/PMNS observables;
    2. seesaw replay and heavy M_R eigenvalues;
    3. updated C_5L, C_5R tensors in the same field basis;
    4. Knu future margin and required S_T after the deformation;
    5. a Pareto table showing whether any point improves both CKM score and
       Knu suppression without spoiling masses or threshold locking.

## 2026-05-08 21:23 Asia/Taipei heartbeat -- shadow/contact joint flavor-Knu scan

status:
  The first-principles GUT derivation remains incomplete.  This update tests
  whether the same transvectant/contact direction can jointly improve the CKM
  flavor obstruction and the Knu d=5 proton bottleneck.

new audit:
  Added:

    code/scan_shadow_contact_flavor_d5_joint.py

  Outputs:

    output/shadow_contact_flavor_d5_joint/summary.json
    output/shadow_contact_flavor_d5_joint/report.md
    output/shadow_contact_flavor_d5_joint/shadow_contact_joint_summary_rows.csv
    output/shadow_contact_flavor_d5_joint/shadow_contact_joint_scan_top_ckm.csv

mathematical construction:
  Starting from the current transvectant flavor point, scan:

    Y_a -> Y_a + epsilon_a K,

  with the CP1 O(2) second transvectant/contact matrix:

    K = (1/sqrt(3)) [[0,0,-1],[0,1,0],[-1,0,0]].

  For each sample the audit recomputes:

    1. CKM magnitudes Vus, Vcb, Vub and J proxy;
    2. mass-ratio scores for up, down, charged-lepton, and neutrino-Dirac
       sectors;
    3. a seesaw replay through the existing convention;
    4. mass-basis Knu Wilson proxy amplitude;
    5. inferred future 1e35 yr margin by rescaling the calibrated Knu target
       margin from output/knu_target_map.

numerical results:
  Samples:

    total = 50005
    modes = down_only, up_down, down_lepton_tied, all_independent
    budgets = 0, 0.015, 0.03, 0.06, 0.10, 0.16

  Representative rows:

    baseline:
      CKM score = 1.446573e-2
      mass score = 1.748238e-1
      Knu amplitude ratio = 1
      inferred future margin = 2.671278e-1

    best_ckm:
      mode = up_down
      CKM score = 4.362117e-4
      mass score = 5.754102
      Knu amplitude ratio = 6.220927e-1
      inferred future margin = 6.902541e-1

    best_knu:
      mode = all_independent
      CKM score = 1.192497
      mass score = 6.746977
      Knu amplitude ratio = 1.747133e-2
      inferred future margin = 8.751189e2

    best_balanced:
      mode = down_only
      CKM score = 5.436502e-3
      mass score = 1.858465e-1
      Knu amplitude ratio = 9.780730e-1
      inferred future margin = 2.792393e-1

  Counts:

    joint_improver_count = 205
    future_safe_count = 1217
    future_safe_with_CKM_gain_count = 0

interpretation:
  The contact/shadow direction is a real Pareto direction: there are 205 points
  that improve both CKM magnitude and the Knu proxy while keeping the mass
  score within 25 percent of the transvectant baseline.  However, no sampled
  point is simultaneously future-safe and CKM-improved.  Therefore a single
  K-line shadow/contact deformation is useful but insufficient.

paper and ledger sync:
  Added "Shadow/contact flavor-Knu Pareto direction" to paper/gut_framework.tex.

  Added one corresponding ledger row:

    Shadow/contact deformation has a joint CKM-Knu Pareto direction.

  The conditional theorem ledger now has 42 claims:

    FAIL: 1
    NO_GO: 3
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 28
    TUNED_FALLBACK: 3

current obstacle:
  A single spin-zero contact/transvectant direction is too small a deformation
  space.  It shows a joint gradient but does not hit the simultaneous target:

    CKM gain + mass tolerance + Knu future 1e35 yr safety.

next attempted nontrivial idea:
  Enlarge the shadow sector from one K-line to the smallest mathematically
  controlled two-direction completion:

    span{K, D_K},

  where D_K should be a covariant finite-N/mock-shadow partner constrained to
  preserve the O(2) mass hierarchy and seesaw replay.  Operationally, test
  two candidates:

    A. K plus a Veronese-image finite-N correction with fixed modular/shadow
       phase relation;
    B. K plus a triplet-only alignment deformation W_AB in the locked-link
       kernel, leaving doublet Yukawas almost unchanged.

proposed next script:

    code/scan_two_direction_shadow_triplet_joint.py

verification plan:
  The next scan should require:

    1. CKM magnitude score below 1e-3 and J score monitored;
    2. mass score no worse than 25 percent above baseline, or a separate
       Pareto front if this is impossible;
    3. inferred Knu future margin above 1;
    4. explicit updated C_5L/C_5R tensor row for the bottleneck channel;
    5. no change to the locked-link threshold vector.

## 2026-05-08 21:31 Asia/Taipei heartbeat -- two-direction shadow/triplet audit

status:
  The first-principles GUT derivation remains incomplete.  This heartbeat
  tests the roadmap's proposed two-direction completion after the one-line
  shadow/contact scan.

new audit:
  Added:

    code/scan_two_direction_shadow_triplet_joint.py

  Outputs:

    output/two_direction_shadow_triplet_joint/summary.json
    output/two_direction_shadow_triplet_joint/report.md
    output/two_direction_shadow_triplet_joint/two_direction_summary_rows.csv
    output/two_direction_shadow_triplet_joint/two_direction_scan_top.csv

mathematical construction:
  The scan separates the doublet flavor shadow deformation from a
  triplet-only colored-Higgs deformation:

    Y_a -> Y_a + epsilon_a K,
    Y_QL^T -> Y_QL^T + delta_T.

  Here K is again the CP1 O(2) spin-zero transvectant:

    K = (1/sqrt(3)) [[0,0,-1],[0,1,0],[-1,0,0]].

  Two triplet-only directions were tested:

    1. contact_K:
       delta_T proportional to K, scaled to the down-type triplet size.

    2. linear_cancel:
       for a selected dangerous Knu entry,

         A_abcd = sum_kl (Y_QL)_kl L_kl^abcd,

       use the Frobenius-minimal diagnostic direction

         D_abcd = conjugate(L^abcd)/||L^abcd||_F,
         zeta_abcd = -A_abcd/||L^abcd||_F.

       This proves whether a triplet-only functional direction exists, but it
       is not a natural mechanism unless later derived from a symmetry or
       mediator Hessian.

numerical results:
  Scan size:

    flavor_samples = 36005
    kept_flavor_cards = 483
    triplet_rows = 38640

  Closure target:

    CKM score < 1e-3
    mass score <= 1.25 * baseline mass score
    inferred Knu future 1e35 yr margin > 1

  Counts:

    local_closure_count = 0
    natural_local_closure_count = 0
    contact_closure_count = 0
    linear_cancel_closure_count = 0
    ckm_future_rows_count = 0

  Representative rows:

    baseline:
      CKM score = 1.446573e-2
      mass score = 1.748238e-1
      Knu amplitude ratio = 1
      inferred future margin = 2.671278e-1

    best_ckm:
      flavor mode = up_down
      triplet mode = contact_K
      CKM score = 8.344852e-4
      mass score = 4.422271
      Knu amplitude ratio = 7.854299e-1
      inferred future margin = 4.330163e-1

    best_knu:
      flavor mode = all_independent
      triplet mode = linear_cancel
      CKM score = 5.137168e-3
      mass score = 6.896499
      Knu amplitude ratio = 4.170101e-1
      inferred future margin = 1.536123
      triplet_delta_over_down_norm = 4.150564e-1

interpretation:
  The minimal two-direction completion is a useful negative result.  It shows
  that triplet-only linear cancellation can make Knu future-safe in isolation,
  and that K contact can improve CKM, but the two goals still do not coincide
  with acceptable mass ratios.  This route should now be recorded as a local
  no-go, not repeated as the next main path.

paper and ledger sync:
  Added the two-direction negative result to paper/gut_framework.tex.

  Added one corresponding no-go ledger row:

    Minimal two-direction K plus triplet-only shadow scan closes flavor and
    d=5 together.

  The conditional theorem ledger now has 43 claims:

    FAIL: 1
    NO_GO: 4
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 28
    TUNED_FALLBACK: 3

current obstacle:
  The bottleneck is no longer just "find a proton filter."  The obstruction is
  joint alignment: directions that improve CKM/masses are not the same
  directions that suppress the physical Knu tensor.

next attempted nontrivial idea:
  Stop scanning ad hoc matrix directions and derive the next deformation from
  an operator-level flavor ansatz.  The most promising route is a full
  CP1-covariant two-kernel Yukawa sector:

    Y_a = lambda_a V[h_a] + rho_a T[h_a],

  where V is the O(-4) Veronese functional and T is the transvectant/contact
  shadow functional.  The triplet couplings must be derived from the same
  operator with explicit doublet-triplet Clebsch phases rather than assigned
  by a posteriori delta_T.

proposed next script:

    code/fit_two_kernel_flavor_then_d5.py

verification plan:
  1. Fit Y_u, Y_d, Y_e, Y_nu from a shared two-kernel operator basis, not
     independent matrix perturbations.
  2. Require CKM, J, charged-lepton/down/up mass ratios, PMNS replay, and
     seesaw perturbativity simultaneously.
  3. Derive the colored-triplet Wilson tensors from the same fitted operator
     coefficients and Clebsch phases.
  4. Only after the flavor fit passes, replay Knu future margin and threshold
     silence.

## 2026-05-08 21:39 Asia/Taipei heartbeat -- operator-level two-kernel flavor and d5 audit

status:
  The first-principles GUT derivation remains incomplete.  This heartbeat
  implements the previous roadmap instruction: stop scanning ad hoc matrix
  perturbations and derive both doublet flavor and colored-triplet tensors
  from one CP1 two-kernel operator ansatz.

new audit:
  Added:

    code/fit_two_kernel_flavor_then_d5.py

  Outputs:

    output/two_kernel_flavor_then_d5/summary.json
    output/two_kernel_flavor_then_d5/report.md
    output/two_kernel_flavor_then_d5/two_kernel_summary_rows.csv

mathematical construction:
  The symmetric flavor tensors are:

    H = V[h_H] + epsilon_H K,
    F = V[h_F] + epsilon_F K,

  where V is the CP1 O(-4) Veronese dual-density product map and

    K = (1/sqrt(3)) [[0,0,-1],[0,1,0],[-1,0,0]]

  is the spin-zero second transvectant/contact functional.  The doublet
  Yukawas are generated by the same Spin(10)-Clebsch ansatz as the existing
  transvectant flavor point, with two antisymmetric 120-like tensors G_A,G_B.

  For d=5 proton decay, no arbitrary triplet matrix delta_T is allowed.  The
  QQ and QL triplet tensors are finite Clebsch-phase reuses of the same
  H,F,G_A,G_B and doublet-mixing coefficients.  The scanned finite profiles
  include:

    doublet_like
    ql_F_minus
    both_F_minus
    anti_flip_ql
    F_quarter_phase

  A convention bug was caught and fixed during this audit: the old numerical
  S_perp differs from the analytic K by a global phase.  The seed epsilon must
  be multiplied by <K,S_perp> to reproduce the previous transvectant point.

numerical results:
  Rows:

    total = 7
    strict_closure_count = 0
    loose_closure_count = 0

  Closure definitions:

    strict:
      CKM score < 1e-3,
      mass score <= 1.25 * seed mass score,
      Knu future 1e35 yr margin > 1.

    loose:
      CKM score < 5e-2,
      mass score <= 0.2,
      Knu future margin > 1.

  Representative rows:

    seed:
      CKM score = 1.446573e-2
      mass score = 1.748238e-1
      Knu future margin = 2.856121e-1
      best finite triplet profile = F_quarter_phase

    best_ckm:
      profile_hint = anti_flip_ql
      best finite triplet profile = F_quarter_phase
      CKM score = 1.041311e-3
      mass score = 1.614449e-1
      Knu future margin = 2.986321e-1

    best_knu:
      profile_hint = both_F_minus
      best finite triplet profile = F_quarter_phase
      CKM score = 1.453203e-3
      mass score = 1.602662e-1
      Knu future margin = 3.200720e-1

interpretation:
  This is the sharpest negative result so far.  The CP1 two-kernel ansatz
  almost reaches the strict CKM target and keeps mass ratios acceptable, but
  the derived finite-Clebsch triplet tensors still fall short of the required
  Knu future margin by roughly a factor sqrt(1/0.320072) in amplitude.  The
  problem is no longer arbitrary matrix fitting; it is a missing operator-level
  triplet suppression mechanism.

paper and ledger sync:
  Added the two-kernel operator-level negative result to paper/gut_framework.tex.

  Added one corresponding no-go ledger row:

    Operator-level two-kernel finite-Clebsch ansatz closes flavor and d=5 together.

  The conditional theorem ledger now has 44 claims:

    FAIL: 1
    NO_GO: 5
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 28
    TUNED_FALLBACK: 3

current obstacle:
  A CP1-covariant two-kernel flavor geometry can nearly repair CKM and masses,
  but the same operator data does not naturally suppress the dangerous
  Knu triplet tensor enough.  Therefore the remaining bottleneck is a genuine
  doublet-triplet operator splitting, not just a family-geometry fit.

next attempted nontrivial idea:
  Derive a doublet-triplet selective zero from representation theory instead
  of finite phase scans.  Concretely, build the SU(5)/Pati-Salam component
  contractions for the 10_H, overline{126}_H, and 120_H triplet couplings and
  search for a Clebsch null in the Knu channel that leaves the doublet Yukawa
  fit intact.

proposed next script:

    code/audit_triplet_clebsch_nullspace_from_components.py

verification plan:
  1. Enumerate QQT and QLT component Clebsch coefficients for 10, 126bar, and
     120 channels in the current field basis.
  2. Build the linear map from operator coefficients to the Knu Wilson tensor.
  3. Compute its nullspace under the constraint that the doublet Yukawa
     matrices remain fixed to the two-kernel fit.
  4. If a nontrivial null exists, replay CKM/mass/seesaw/Knu margin; if not,
     record a stronger no-go for this minimal Higgs flavor sector.

## 2026-05-08 21:54 Asia/Taipei heartbeat -- triplet Clebsch component nullspace

status:
  The first-principles GUT derivation remains incomplete.  This heartbeat
  executes the previous roadmap item: keep the successful CP1 two-kernel
  doublet data fixed and ask whether the triplet Clebsch/component space has
  a genuine Knu null rather than another finite phase scan.

new audit:
  Added:

    code/audit_triplet_clebsch_nullspace_from_components.py

  Outputs:

    output/triplet_clebsch_nullspace/summary.json
    output/triplet_clebsch_nullspace/report.md
    output/triplet_clebsch_nullspace/component_nullspace_rows.csv

mathematical construction:
  With the two-kernel flavor fit fixed, expand the triplet couplings as

    Y_QQ^T = sum_I p_I B_Q^I,    B_Q = (H, r_u F, a_u G_A, b_u G_B),
    Y_QL^T = sum_J q_J B_L^J,    B_L = r_d (H, F, a_d G_A, b_d G_B).

  For the nine monitored Knu entries, the bilinear Wilson map is

    C_a = sum_{I,J} p_I M_a^{IJ} q_J,

  where a runs over the three (u,u,s) permutations times three neutrino
  flavors.  The audit computes fixed-side nullspaces and an explicit
  rank-one finite-norm null in p_I q_J.

numerical results:
  Linear algebra:

    fixed QQ nullity in QL coefficients = 0
    fixed QL nullity in QQ coefficients = 2
    unrestricted p tensor q nullity = 10
    rank-one null residual = 0
    rank-one null is antisymmetric QQ only = True

  Representative rows:

    natural_doublet_like:
      future Knu margin = 2.671278e-1
      amplitude ratio = 1.0

    fixed_QQ_smallest_QL_singular:
      future Knu margin = 3.979722e-1
      amplitude ratio = 8.192816e-1
      QL deformation ratio = 1.611938

    fixed_QL_smallest_QQ_singular:
      future Knu margin = 2.671278e29
      amplitude ratio = 0
      QQ deformation ratio = 1.0

    rank_one_bilinear_null_search:
      future Knu margin = 2.671278e29
      residual = 0

interpretation:
  The new nontrivial element is not another tuned scalar filter.  The null is

    p = (0, 0, 2, 0),    q = (1, 1, 1, 1).

  Thus QQT is pure antisymmetric 120-like, while QLT remains natural.  The
  monitored LLLL Knu operator depends on sym(Y_QQ^T), so this component-level
  choice kills the dangerous tensor exactly in the audited map.  This gives a
  concrete action-level target: derive a triplet mass-matrix / Clebsch sector
  that sends QQT into the antisymmetric 120-like direction while preserving
  the doublet flavor fit and a natural QLT side.

important caveat:
  This is not yet a full d=5 proton decay proof.  The audit closes the
  monitored LLLL Knu tensor at component level.  It still has to be checked
  against RRRR operators, complete 120_H component matching, color-triplet
  mass mixing, and the full dressed proton width pipeline.

paper and ledger sync:
  Added the component nullspace result to paper/gut_framework.tex.

  Added one corresponding conditional ledger row:

    Component triplet Clebsch nullspace contains an antisymmetric QQ Knu null.

  The conditional theorem ledger now has 45 claims:

    FAIL: 1
    NO_GO: 5
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 29
    TUNED_FALLBACK: 3

  paper/gut_framework.tex compiled twice after the new paragraph/operator
  definitions, and the final log scan found no LaTeX warnings, package
  warnings, overfull boxes, undefined references, or errors.

current obstacle:
  The remaining bottleneck is no longer whether a component-level Knu null
  exists; it does.  The bottleneck is whether a realistic Spin(10)/PS Higgs
  triplet sector can realize that null without breaking the successful
  doublet flavor fit, without regenerating RRRR proton decay, and without
  introducing nonuniversal threshold damage.

next attempted nontrivial idea:
  Build the explicit 120_H triplet consistency audit.  The key question is
  whether the SO(10)/Pati-Salam component contractions and triplet mass matrix
  can make QQT purely antisymmetric while keeping QLT natural and keeping
  doublets fixed.

proposed next script:

    code/audit_120_triplet_rrrr_and_mass_matrix.py

verification plan:
  1. Derive the 120_H triplet component couplings to QQ, QL, u^c e^c, and
     u^c d^c in the current field basis.
  2. Check whether the pure antisymmetric QQT LLLL null is compatible with
     the same SO(10) Clebsch assignments that generate the doublet Yukawas.
  3. Add the RRRR Wilson tensors and replay all monitored channels.
  4. Insert the resulting color-triplet mass eigenvalues into the threshold
     and proton scans; if the null requires split incomplete multiplets,
     record it as a tuned fallback rather than a natural closure.

## 2026-05-08 22:36 Asia/Taipei heartbeat -- 120-like RRRR and finite mass-matrix audit

status:
  The first-principles GUT derivation remains incomplete, but the d=5
  proton-side obstruction moved again.  The previous LLLL Knu null was tested
  against RRRR operators and a finite triplet mass-matrix lift.

new audit:
  Added:

    code/audit_120_triplet_rrrr_and_mass_matrix.py

  Outputs:

    output/triplet_120_rrrr_mass_matrix/summary.json
    output/triplet_120_rrrr_mass_matrix/report.md
    output/triplet_120_rrrr_mass_matrix/triplet_120_rows.csv

mathematical construction:
  Use the locked component bookkeeping:

    Y_QQ^T and Y_UE^T use the same p_I B_10^I,
    Y_QL^T and Y_UD^T use the same q_J B_5^J.

  The LLLL and RRRR maps are

    C_L^a = sum_{I,J} p_I q_J M_{L,a}^{IJ},
    C_R^b = sum_{I,J} p_I q_J M_{R,b}^{IJ}.

  LLLL uses sym(B_10^I) on QQ.  RRRR uses the full B_10^I on the UE side in
  this component proxy.  The finite mass-matrix lift is

    W_kappa = U diag(1, 1/kappa, 1/kappa, 1/kappa) V^\dagger,

  where U,V complete the exact rank-one p,q directions to unitary bases.

numerical results:
  Linear algebra:

    fixed antisymmetric p RRRR q-nullity = 2
    fixed antisymmetric p RRRR residual = 1.019781e-25
    joint W rank = 10
    joint W nullity = 6
    joint rank-one LLLL+RRRR residual = 1.019781e-25
    exact joint rank-one null = True

  Physical p,q rows under the common ST=7.5e-6 item-4 width proxy:

    natural_locked:
      LLLL amplitude = 2.173374e-3
      RRRR amplitude = 1.934137e-6
      worst 1e35 margin = 6.564794e2

    antisymmetric_120_QQT_locked_natural_QLT:
      LLLL amplitude = 0
      RRRR amplitude = 4.122672e-4
      worst 1e35 margin = 1.824452e4

    antisymmetric_120_QQT_with_best_RRRR_q:
      LLLL amplitude = 0
      RRRR amplitude = 6.564897e-19
      worst 1e35 margin = 7.195058e33
      q deformation ratio = 1.190939

  The exact joint null is

    p = (0,0,2,0),
    q ~= (0,0,0,1.62686 + 1.16332 i),

  with only ~1e-5 leakage in the first three q entries.  In words: use the
  G_A antisymmetric source for the 10-10 triplet coupling and the orthogonal
  G_B-like source for the 10-5bar triplet coupling.

finite mass-matrix leakage:
  Equal orthogonal singular-value lifts reintroduce leakage at O(1/kappa).
  For kappa=100:

    LLLL leakage ratio = 5.951391e-5
    RRRR leakage ratio = 1.522444e-2

  For kappa=1000:

    LLLL leakage ratio = 5.951391e-6
    RRRR leakage ratio = 1.522502e-3

interpretation:
  This is the strongest proton-side structural clue so far.  The component
  algebra admits a simultaneous LLLL Knu and RRRR uusd null, but only as a
  crossed two-120 triplet projector: 10-10 selects G_A and 10-5bar selects a
  G_B-like orthogonal antisymmetric source.  The finite mass-matrix lift is
  not free: unless the orthogonal triplet directions are symmetry/projector
  suppressed, they leak back into RRRR at the percent level for kappa=100.

paper and ledger sync:
  Added the crossed 120A/120B result to paper/gut_framework.tex.

  Added one corresponding conditional ledger row:

    Crossed 120A/120B triplet projector has a component-level LLLL+RRRR null.

  The conditional theorem ledger now has 46 claims:

    FAIL: 1
    NO_GO: 5
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 30
    TUNED_FALLBACK: 3

  paper/gut_framework.tex compiled twice after the crossed-120 update, and
  the final log scan found no LaTeX warnings, package warnings, overfull
  boxes, undefined references, or errors.

current obstacle:
  The bottleneck is now action-level origin, not component-level existence.
  We need a renormalizable or constrained Spin(10)/Pati-Salam triplet
  superpotential that realizes the crossed G_A/G_B projector, keeps doublet
  flavor unchanged, and makes the finite-lift leakage threshold-safe.

next attempted nontrivial idea:
  Build a crossed-two-120 triplet projector sector.  The natural candidate is
  a 2x2 120_A/120_B color-triplet mass matrix plus a grading that allows
  10-10 to couple only to 120_A and 10-5bar only to 120_B in the dangerous
  triplet channel, while the light doublet Yukawa fit still sees the original
  two-120 combination.

proposed next script:

    code/construct_crossed_120_triplet_projector.py

verification plan:
  1. Write the minimal 120_A,120_B triplet source/mass matrix W_AB that
     realizes p=(0,0,2,0), q~G_B as an inverse-propagator subblock.
  2. Enumerate the required charge/grading selection rules for triplet and
     doublet components separately.
  3. Feed the finite W_kappa leakage ratios into the threshold and dressed
     proton scans for kappa=30,100,300,1000.
  4. If the selection rules force incomplete light multiplets, record the
     mechanism as a tuned EFT fallback; if they are complete/constrained,
     promote it to a conditional closure candidate.

## 2026-05-09 00:12 Asia/Taipei heartbeat -- crossed 120 projector construction audit

status:
  The first-principles GUT derivation remains incomplete.  The component
  crossed-120 proton null is now promoted into an explicit projector audit.
  The result is mixed but useful: it is a strong post-breaking
  PS/constrained EFT mechanism, and simultaneously a no-go for a field-only
  unbroken Spin(10) grading origin.

new audit:
  Added:

    code/construct_crossed_120_triplet_projector.py

  Outputs:

    output/crossed_120_triplet_projector/summary.json
    output/crossed_120_triplet_projector/report.md
    output/crossed_120_triplet_projector/finite_leakage_replay.csv
    output/crossed_120_triplet_projector/dilation_rows.csv

mathematical construction:
  The target visible inverse triplet block is

    W_vis = |120_A><120_B|.

  This realizes the previous component null:

    10-10 triplet source = G_A,
    10-5bar triplet source = G_B-like.

  The algebraic threshold-locked realization uses the Julia-Sz.-Nagy
  dilation:

    U(W) = [[W, sqrt(I-W W^\dagger)],
            [sqrt(I-W^\dagger W), -W^\dagger]],

    P (M_lock U^\dagger)^(-1) P = W / M_lock.

field-only Spin(10) no-go:
  If both Yukawa terms

    16 16 120_A,    16 16 120_B

  are invariant under a field-only abelian grading, then

    q_A = q_B = -2 q_16.

  Therefore

    q_A+q_A = q_A+q_B = q_B+q_A = q_B+q_B.

  Any field-only symmetry that allows one 120_i 120_j mass entry allows all
  four.  The brute-force search over Z_N, N<=24, found zero crossed-only
  witnesses.  This means the crossed projector cannot be honestly advertised
  as an unbroken Spin(10) field-charge consequence.

conditional PS/constrained result:
  As a post-Spin(10)-breaking Pati-Salam/constrained triplet projector, the
  construction is viable:

    field_only_unbroken_spin10_projector_possible = False
    ps_eft_or_constrained_projector_possible = True
    all_finite_leakage_rows_future_safe = True
    all_unitary_dilations_clean = True

  The common locked mass is

    M_lock = 3.887852e15 GeV.

  The dilated complete-multiplet construction has projected threshold l2 = 0
  for all audited kappas.

finite leakage replay:
  Replaying the finite W_kappa leakage against the 1e35 yr stress target:

    kappa=30:
      LLLL margin = 1.668117e10
      RRRR margin = 3.218681e11

    kappa=100:
      LLLL margin = 1.853464e11
      RRRR margin = 3.576291e12

    kappa=300:
      LLLL margin = 1.668117e12
      RRRR margin = 3.218608e13

  Direct split realizations are below the reduced Planck benchmark for
  kappa=30,100,300, but not for 1000 or 10000.  The unitary-dilated complete
  realization remains degenerate and threshold silent for all audited kappas.

paper and ledger sync:
  Added the crossed-projector construction and field-only Spin(10) no-go to
  paper/gut_framework.tex.

  Recompiled paper/gut_framework.tex twice with pdflatex after the insertion.
  The refreshed paper/gut_framework.pdf is 100 pages, and the log scan for
  LaTeX warnings, package warnings, overfull boxes, undefined references, and
  errors is clean.

  Added two ledger rows:

    Field-only unbroken Spin(10) grading enforces the crossed 120 triplet projector.
    Post-breaking constrained crossed-120 projector is threshold-silent and d=5 safe in the replay.

  The conditional theorem ledger now has 48 claims:

    FAIL: 1
    NO_GO: 6
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 31
    TUNED_FALLBACK: 3

current obstacle:
  The crossed 120 mechanism is not blocked by proton algebra anymore.  The
  bottleneck is now deriving a post-breaking/constrained projector from an
  explicit Higgs/source action without damaging doublet flavor or introducing
  incomplete nonuniversal thresholds.

next attempted nontrivial idea:
  Build the PS-level source action that realizes the crossed projector after
  Spin(10) breaking.  It should use PS projectors to distinguish triplets from
  doublets, enforce W_vis=|120_A><120_B| in the triplet block, leave the
  doublet Yukawa combination unchanged, and place all mediator/source
  remnants into complete or constrained threshold-silent multiplets.

proposed next script:

    code/construct_ps_crossed_120_source_action.py

verification plan:
  1. Define PS fragment labels of 120_A and 120_B relevant to color triplets
     and electroweak doublets.
  2. Write a constrained source/action that projects only the triplet inverse
     block to |A><B| while leaving the doublet Yukawa fit untouched.
  3. Compute F-flatness at zero triplet vev and the mass eigenvalues of
     triplet, doublet, and mediator fragments.
  4. Feed any non-complete finite fragments into the threshold and proton
     scans; if the projected threshold is nonzero, classify the mechanism as a
     tuned PS EFT fallback rather than a conditional closure candidate.

## 2026-05-09 00:29 Asia/Taipei heartbeat -- PS crossed 120 source action

status:
  The first-principles GUT derivation remains incomplete, but the crossed
  120A/120B branch now has a concrete post-breaking Pati-Salam component
  source action.  The new result sharply separates the viable constrained
  source interpretation from the nonviable literal triplet-only mediator
  interpretation.

new audit:
  Added:

    code/construct_ps_crossed_120_source_action.py

  Outputs:

    output/ps_crossed_120_source_action/summary.json
    output/ps_crossed_120_source_action/report.md
    output/ps_crossed_120_source_action/offblock_mixing_scan.csv
    output/ps_crossed_120_source_action/threshold_interpretation_scan.csv

mathematical construction:
  Use the Pati-Salam decomposition

    120 -> (1,2,2) + (15,2,2) + (6,1,1) + (10,1,1) + (10bar,1,1).

  The doublet fragments are in (1,2,2)+(15,2,2).  The color-triplet
  fragments sit inside the remaining PS pieces.  The component source
  ansatz is

    W_src =
      M_lock bar(Phi_T)^T U_T(W_kappa)^dagger Phi_T
      + M_lock bar(D)^T D
      + rho M_lock bar(Phi_T)^T C D,

  where the rho term is only a deterministic off-block leakage probe.  For
  kappa=100,

    W_kappa = [[0,1],[1e-2,0]].

numerical validation:
  At rho=0:

    F_norm_at_zero = 0
    triplet_inverse_subblock_residual_fro = 0
    triplet_unitary_residual_2norm = 6.245109e-19
    triplet_mass_singular_spread_over_Mlock = 3.330669e-16

  The doublet block is exactly diagonal.  In the off-block leakage scan:

    rho=1e-3:
      triplet inverse residual = 7.071068e-9
      doublet inverse residual = 5.000250e-7
      proxy worst 1e35 yr margin = 1.853023e11

    rho=1e-2:
      proton proxy remains safe, but doublet residual = 5.000250e-5,
      so it fails the 1e-6 doublet-preservation criterion.

threshold interpretation:
  With M_lock=3.887852e15 GeV and M_G=7.079458e15 GeV:

    constrained_nonpropagating_source:
      projected_l2 = 0
      within R=200 budget = True

    two_complete_5_plus_5bar_mediator_pairs:
      projected_l2 = 4.807407e-17
      within R=200 budget = True

    two_literal_triplet_mediator_pairs:
      projected_l2 = 1.357954e-1
      budget ratio = 270.36
      within R=200 budget = False

  Therefore the PS component source action is conditionally viable only as a
  constrained/auxiliary source or as a completed inert 5+5bar mediator
  package.  A literal triplet-only propagating PS mediator completion is a
  threshold no-go.

paper and ledger sync:
  Added the PS source-action audit to paper/gut_framework.tex.

  Recompiled paper/gut_framework.tex twice after the insertion.  The refreshed
  paper/gut_framework.pdf is 101 pages, and the log scan for LaTeX warnings,
  package warnings, overfull boxes, undefined references, and errors is clean.

  Added two ledger rows:

    Post-breaking PS crossed-120 source action realizes the projector without doublet mixing.
    Literal triplet-only crossed-120 mediator completion is threshold silent.

  The conditional theorem ledger now has 50 claims:

    FAIL: 1
    NO_GO: 7
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 32
    TUNED_FALLBACK: 3

current obstacle:
  The crossed 120 source action is now a valid PS-level constrained closure
  candidate, but it is still not a first-principles microscopic origin.  The
  remaining hard issue is deriving the constrained/auxiliary source, or the
  inert completed 5+5bar partner package, from a holomorphic UV action without
  reintroducing Landau-pole or doublet-flavor problems.

next attempted nontrivial idea:
  Build the completed inert-doublet partner action and compare it to a
  genuinely constrained auxiliary-source action.  The useful distinction is:
  a completed 5+5bar mediator is threshold silent but may affect UV
  perturbativity and doublet flavor; a constrained source is threshold silent
  but still needs a microscopic origin.

proposed next script:

    code/construct_completed_120_partner_action.py

verification plan:
  1. Add inert doublet partners for the two Julia defect triplet pairs and
     enforce a common M_lock.
  2. Check F-flatness, doublet off-block mixing, and whether the inert
     doublets can be symmetry-decoupled from physical Yukawa operators.
  3. Compute UV beta contribution of the completed package and Landau-pole
     ratio above M_G.
  4. Re-run the threshold and proton replay with any residual noncomplete
     splitting; if UV perturbativity fails, keep only the constrained-source
     branch as the conditional theorem path.

## 2026-05-09 00:39 Asia/Taipei heartbeat -- completed inert 5+5bar partner audit

status:
  The first-principles GUT derivation remains incomplete.  The literal
  triplet-only crossed-120 mediator no-go now has a minimal propagating repair:
  add inert doublet partners to make two mass-locked 5+5bar packages.  This
  repairs the one-loop nonuniversal threshold, but introduces a per-mille
  mass-locking condition and still does not derive the crossed projector.

new audit:
  Added:

    code/construct_completed_120_partner_action.py

  Outputs:

    output/completed_120_partner_action/summary.json
    output/completed_120_partner_action/report.md
    output/completed_120_partner_action/split_threshold_scan.csv
    output/completed_120_partner_action/uv_beta_scan.csv

mathematical construction:
  The completion is

    W_comp = W_T[crossed source] + M_L sum_a L_a Lbar_a,

  with two inert doublet pairs completing the two triplet mediator pairs.
  Exact threshold cancellation requires

    M_L = M_T = M_lock.

  At exact locking the two 5+5bar pairs give the universal threshold

    2 (b_T+b_D) log(M_G/M_lock)/(2 pi)
      = (0.190776, 0.190776, 0.190776),

  hence projected_l2 = 4.807407e-17.

mass-locking validation:
  Let xi = M_D/M_T.  The nonuniversal residual is

    ||P Delta(xi)||_2 = |log xi| ||P(2 b_D)||_2/(2 pi).

  For the R=200 budget 5.022739e-4:

    max |log xi| = 2.216814e-3
    xi window = [0.997785641, 1.002219273]
    percent window = 0.221927%

  Scan result:

    xi=0.998,0.999,1.000,1.001,1.002 pass.
    xi=0.995 and xi=1.005 fail.

inert grading:
  A post-GUT Z2_inert can decouple the added doublet partners:

    L_inert Lbar_inert: allowed
    16_i 16_j L_inert: forbidden
    H_phys_doublet Lbar_inert: forbidden
    Phi_T crossed-source mass block: allowed

  This protects the physical doublet Yukawa fit from inert partner pollution.

UV beta audit:
  If embedded as two 10_G copies:

    two_10G package alone: b10 = -22, no one-loop Landau pole.
    minimal Yukawa sector + two 10_G: b10 = 48,
      Lambda_LP/M_G = 224.621 > 200.

  The completed package is therefore UV-benign for the R=200 target in this
  restricted audit, unlike the previous large 54/210 drive-sector problem.

paper and ledger sync:
  Added the completed partner audit to paper/gut_framework.tex.

  Recompiled paper/gut_framework.tex twice after the insertion.  The refreshed
  paper/gut_framework.pdf is 101 pages, and the log scan for LaTeX warnings,
  package warnings, overfull boxes, undefined references, and errors is clean.

  Added two ledger rows:

    Inert completed 5+5bar partner package repairs the crossed-120 threshold.
    Inert completed 5+5bar partners can be freely split from the triplet mass.

  The conditional theorem ledger now has 52 claims:

    FAIL: 1
    NO_GO: 8
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 33
    TUNED_FALLBACK: 3

current obstacle:
  The threshold issue for a propagating crossed-120 completion can be repaired
  by mass-locked inert 5+5bar partners.  The remaining bottleneck is now
  microscopic: derive the crossed triplet source block and the inert partner
  mass-locking from a holomorphic post-GUT or Spin(10)-descended action,
  rather than imposing them as component-level constraints.

next attempted nontrivial idea:
  Search for a single holomorphic post-GUT link action that simultaneously
  generates the crossed triplet Julia block and locks the inert doublet partner
  masses.  The key creative element should be a common spurion or Lagrange
  multiplier enforcing both W_T and M_L=M_T, so the per-mille mass locking is
  no longer an independent tuning.

proposed next script:

    code/derive_crossed_120_link_locking_action.py

verification plan:
  1. Write a holomorphic link/driver superpotential with a common singlet
     spurion S_lock coupling to triplet Julia fields and inert doublet
     partners.
  2. Solve F-flatness to check whether W_T=|A><B|+epsilon|B><A| and M_L=M_T
     follow from the same equations.
  3. Compute Hessian eigenvalues and off-block doublet mixing.
  4. Recompute split_threshold_scan under small allowed deformations of the
     lock equations; if the F-term lock only fixes masses to O(1%), classify
     the completed partner branch as tuned fallback.

## 2026-05-09 00:48 Asia/Taipei heartbeat -- crossed 120 link-locking action

status:
  The first-principles GUT derivation remains incomplete.  The common-spurion
  link action can derive the crossed inverse block and a shared mass scale, but
  pure holomorphic F-term constraints are not enough to guarantee the unitary
  mass lock needed for threshold safety.  A D-term/Kahler/NLSM or composite
  unitary constraint is still required.

new audit:
  Added:

    code/derive_crossed_120_link_locking_action.py

  Outputs:

    output/crossed_120_link_locking_action/summary.json
    output/crossed_120_link_locking_action/report.md
    output/crossed_120_link_locking_action/holomorphic_moduli_deformation_scan.csv

link action:
  The tested holomorphic action is

    W =
      Tr X(AB-I) + Tr Y(PBP-W_kappa) + Z(S_lock-M_lock)
      + S_lock [bar(Phi_T)^T A Phi_T + L_inert Lbar_inert].

  The constraints AB=I and PBP=W_kappa derive the crossed inverse triplet
  block.  The common S_lock ties the triplet-link and inert-doublet mass
  scales.

reference-point validation:
  At A0=U(W_kappa)^dagger and B0=U(W_kappa):

    constraint_residual = 8.831917e-19
    unitarity_residual_2norm = 6.245109e-19
    zero_field_matter_F_norm = 0

  With driver fields X=Y=Z=0, all nonconstraint F-terms vanish at the
  zero-field point.

constraint geometry:
  The linearized holomorphic constraint system has

    real variables = 64
    real equations = 40
    real rank = 40
    real nullity = 24

  This proves the holomorphic constraints fix the inverse block but leave
  residual moduli.

moduli deformation scan:
  Exact deformations with AB=I and PBP=W_kappa preserved at ~1e-16 give:

    amplitude=1e-3:
      max |log singular| = 4.764397e-4
      inside mass-lock window = True

    amplitude=3e-3:
      max |log singular| = 1.429323e-3
      inside mass-lock window = True

    amplitude=1e-2:
      max |log singular| = 4.764439e-3
      inside mass-lock window = False

  Since the completed partner branch allows only

    max |log xi| = 2.216814e-3,

  generic holomorphic moduli can break threshold-safe mass locking.

paper and ledger sync:
  Added the link-locking audit to paper/gut_framework.tex.

  Recompiled paper/gut_framework.tex twice after the insertion.  The refreshed
  paper/gut_framework.pdf is 101 pages, and the log scan for LaTeX warnings,
  package warnings, overfull boxes, undefined references, and errors is clean.

  Added two ledger rows:

    Common-spurion crossed-120 link action derives the inverse block and shared mass scale.
    Pure holomorphic crossed-120 link F-terms are sufficient for threshold-safe mass locking.

  The conditional theorem ledger now has 54 claims:

    FAIL: 1
    NO_GO: 9
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 34
    TUNED_FALLBACK: 3

current obstacle:
  The current crossed-120 completion has a precise conditional form:
  holomorphic link constraints plus a common S_lock plus inert completed
  5+5bar partners.  The remaining obstacle is deriving the unitary constraint
  on the link from a microscopic D-term/Kahler/composite mechanism rather than
  imposing it as an NLSM constraint.

next attempted nontrivial idea:
  Build a D-term/Kahler quotient or composite meson model whose moduli space
  contains a unitary link variable with PBP=W_kappa.  The goal is to replace
  the imposed NLSM unitarity condition by an orbit-stabilizer or strong-sector
  quotient that leaves only unitary link moduli and keeps the inert 5+5bar
  mass lock.

proposed next script:

    code/construct_unitary_link_dterm_quotient.py

verification plan:
  1. Define a toy U(4)_L x U(4)_R link field with D-flat constraint
     L^\dagger L = I up to gauge quotient.
  2. Add holomorphic drivers for P L P = W_kappa and common S_lock.
  3. Compute tangent rank after quotienting gauge/orbit directions and verify
     that nonunitary holomorphic deformations are removed.
  4. Re-run mass-lock and threshold scans under the remaining quotient moduli.
     If any residual physical modulus violates max |log xi|<2.216814e-3,
     classify the branch as an NLSM assumption rather than a derived action.

## 2026-05-09 00:52 Asia/Taipei heartbeat -- unitary link D-term quotient audit

status:
  The first-principles GUT derivation remains incomplete, but the previous
  "pure F-term is insufficient" obstruction has a local resolution.  A
  D-term/Kahler quotient imposing B^\dagger B=I removes the dangerous
  nonunitary holomorphic moduli and leaves only unitary completion moduli,
  which do not split the singular values or violate the completed-partner
  mass-locking window.

new audit:
  Added:

    code/construct_unitary_link_dterm_quotient.py

  Outputs:

    output/unitary_link_dterm_quotient/summary.json
    output/unitary_link_dterm_quotient/report.md
    output/unitary_link_dterm_quotient/finite_unitary_completion_scan.csv

D-term/Kahler quotient:
  Impose the moment-map constraint

    mu = B^\dagger B - I = 0

  and then the holomorphic fixed-block driver

    P B P = W_kappa.

  On the 32 real components of B:

    D-term real equations = 16
    D-term linear rank = 16
    D-flat tangent real dimension = 16

  On this D-flat tangent space:

    fixed-block real equations = 8
    fixed-block rank = 7
    residual unitary moduli = 9

  The rank is 7 rather than 8 because W_kappa lies on the contraction
  boundary ||W_kappa||_2=1.

finite unitary completion validation:
  Random finite defect-space unitary completions with angles from 1e-6 to 1
  obey:

    fixed_block_residual_fro = 0
    max unitarity_residual_2norm = 1.784572e-15
    max |log singular| <= 6.661338e-16

  Therefore all sampled completions remain inside the completed-partner
  mass-lock window

    max |log xi| = 2.216814e-3.

paper and ledger sync:
  Added the D-term/Kahler quotient audit to paper/gut_framework.tex.

  Recompiled paper/gut_framework.tex twice after the insertion.  The refreshed
  paper/gut_framework.pdf is 102 pages, and the log scan for LaTeX warnings,
  package warnings, overfull boxes, undefined references, and errors is clean.

  Added one ledger row:

    D-term/Kahler quotient removes dangerous crossed-120 link moduli.

  The conditional theorem ledger now has 55 claims:

    FAIL: 1
    NO_GO: 9
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 35
    TUNED_FALLBACK: 3

current obstacle:
  The unitary link is now locally derived as a D-term/Kahler quotient ansatz,
  but not yet as a microscopic first-principles sector.  The remaining
  obstacle is to realize this quotient from a concrete strong/composite sector
  or an explicit gauged linear sigma model with controlled UV beta function.

next attempted nontrivial idea:
  Promote the quotient to a microscopic GLSM/composite model: introduce
  bifundamental fields Q,\tilde Q under a hidden U(4) gauge group whose
  meson B=Q\tilde Q/f^2 becomes the unitary link after D-flatness and
  quotienting.  Then check whether the hidden gauge coupling and extra matter
  can remain perturbative or safely confined below the GUT cutoff.

proposed next script:

    code/construct_composite_unitary_link_glsm.py

verification plan:
  1. Define hidden-gauge charges for Q,\tilde Q and the visible source fields.
  2. Derive the D-flat quotient and identify the meson coordinate B.
  3. Check that PBP=W_kappa can be imposed by holomorphic drivers without
     reintroducing nonunitary moduli.
  4. Compute hidden-sector beta functions and any SM/GUT threshold
     contribution.  If the hidden sector must be strongly coupled below M_G,
     classify it as a composite EFT assumption rather than an elementary UV
     completion.

## 2026-05-09 01:03 Asia/Taipei heartbeat -- composite unitary-link hidden GLSM audit

status:
  The first-principles GUT derivation remains incomplete, but the unitary-link
  origin is now sharper.  The D-term/Kahler quotient has a concrete hidden
  U(4)_H composite/GLSM caricature whose meson reproduces the locked unitary
  link with zero visible threshold.  This removes the immediate visible
  threshold and UV-Landau obstruction for the link variable, but it does not
  yet prove the full nonperturbative hidden dynamics.

new audit:
  Added:

    code/construct_composite_unitary_link_glsm.py

  Outputs:

    output/composite_unitary_link_glsm/summary.json
    output/composite_unitary_link_glsm/report.md
    output/composite_unitary_link_glsm/finite_meson_samples.csv
    output/composite_unitary_link_glsm/hidden_beta_scan.csv

mathematical construction:
  Introduce a visible-singlet hidden U(4)_H sector with Q and Qtilde and define
  the link meson

    B = Qtilde Q / f^2.

  On the locked-radial branch

    Q = f I_4,
    Qtilde = f B,
    B in U(4),

  the hidden moment maps vanish:

    Q Q^dagger - Qtilde^dagger Qtilde = 0,
    Q^dagger Q - Qtilde Qtilde^dagger = 0.

  The same holomorphic driver then imposes

    P B P = W_kappa.

  Dimension bookkeeping:

    raw Q,Qtilde real components = 64
    hidden D-term rank = 16
    hidden gauge orbit rank = 16
    radial-lock rank = 16
    unitary-link real dimension = 16
    fixed-block rank on unitary tangent = 7
    residual unitary moduli = 9

numerical verification:
  Finite meson samples using the audited unitary completions give:

    all_finite_meson_samples_dflat_and_locked = True
    max D residual = 1.778812e-15
    max unitarity residual = 1.778803e-15
    max |log singular| = 8.881784e-16
    completed-partner lock window = 2.216814e-3
    visible threshold vector = (0,0,0)

hidden beta audit:
  With convention

    d alpha_H^{-1}/d log mu = - b_H/(2 pi),
    b_H = N_f - 3 N_c,  N_c=4,

  the minimal N_f=4 route has

    b_H = -8,
    phase = Nf=Nc quantum-deformed/constrained meson.

  It is UV-Landau safe to R=200, but should be interpreted as a constrained
  composite meson branch rather than a weakly coupled UV proof.

  The N_f=8 route has

    b_H = -4,
    phase = inside conformal window.

  This is the cleaner microscopic direction if hidden spectators are allowed.

paper and ledger sync:
  Added the hidden U(4)_H meson audit to paper/gut_framework.tex in the
  crossed-120/unitary-link proof chain.

  Regenerated the conditional theorem ledger.  It now has 56 claims:

    FAIL: 1
    NO_GO: 9
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 36
    TUNED_FALLBACK: 3

  New ledger row:

    Hidden U(4) composite GLSM can realize the unitary crossed-120 link locally.

  The TeX draft was compiled twice from the paper directory.  The refreshed
  paper/gut_framework.pdf has 103 pages.  A log scan for Warning, undefined,
  Overfull, Underfull, Error, and Fatal returned no matches.

current obstacle:
  O3 is reduced but not closed.  The unitary link has a threshold-silent
  composite meson realization, yet the actual nonperturbative hidden-sector
  dynamics, radial stabilization, and compatibility with the rest of A5--A7
  still need to be derived rather than assumed.

next attempted nontrivial idea:
  Decide whether to make the hidden link branch an N_f=N_c quantum-deformed
  meson model or an N_f=8 conformal-window model.  The sharper route is:

    1. write a Seiberg-style effective superpotential with meson B and
       baryon/determinant constraint for N_f=N_c;
    2. add radial/source drivers and check F-flatness plus D-flatness;
    3. verify that the determinant/quantum constraint does not force
       det B away from the unitary completions;
    4. if N_f=4 fails, switch to N_f=8 and audit spectator decoupling,
       conformal anomalous dimensions, and induced visible thresholds.

proposed next script:

    code/construct_quantum_deformed_link_moduli.py

verification plan:
  1. Use B=Qtilde Q/f^2 with det B and baryon variables for the N_f=N_c hidden
     branch.
  2. Check whether the unitary completions satisfy the quantum-deformed
     constraint after allowed baryon vevs.
  3. Compute the Hessian of the meson/baryon/radial driver sector and verify
     that no nonunitary light moduli reappear.
  4. Record whether this closes O3 conditionally, or whether the model must
     move to the N_f=8 conformal-window branch.

## 2026-05-09 01:15 Asia/Taipei heartbeat -- quantum-deformed hidden-link moduli audit

status:
  The task is still incomplete.  The N_f=N_c hidden-SQCD quantum deformation is
  compatible with the required unitary crossed-120 link, but it does not derive
  the unitary lock by itself.  This is a useful split result: it validates the
  hidden composite branch as a consistent moduli-space embedding, while adding
  a no-go against using the quantum deformation alone to close O3.

new audit:
  Added:

    code/construct_quantum_deformed_link_moduli.py

  Outputs:

    output/quantum_deformed_link_moduli/summary.json
    output/quantum_deformed_link_moduli/report.md
    output/quantum_deformed_link_moduli/finite_quantum_samples.csv

mathematical construction:
  For hidden SU(4) SQCD with N_f=N_c=4, use the gauge-invariant moduli

    M = Qtilde Q,
    Baryon = calB,
    anti-Baryon = calBtilde,

  with quantum-deformed constraint

    det M - calB calBtilde = Lambda_H^8.

  In dimensionless variables

    m = M/f^2,
    beta = calB/f^4,
    betatilde = calBtilde/f^4,

  the tested equations are

    det m - beta betatilde = (Lambda_H/f)^8,
    P m P = W_kappa.

  For each sampled unitary link m, choose

    beta = betatilde = sqrt(det m - (Lambda_H/f)^8).

numerical verification:
  The samples use Lambda_H/f in {0.3, 0.7, 1.0} and the same unitary
  completion angle grid as the D-term quotient audit.  Results:

    all_unitary_completions_satisfy_quantum_constraint_with_baryons = True
    max quantum constraint residual = 2.220446e-16
    max baryon magnitude = 1.402261
    visible threshold vector = (0,0,0)

linearized rank audit:
  Holomorphic fixed-block plus quantum-deformed constraint:

    complex variables = 18
    complex constraints = 5
    complex rank = 5
    complex nullity = 13
    real nullity = 26

  With the D-term/radial unitary lock retained and baryon variations included:

    real variables = 20
    real rank = 9
    real nullity = 11

  Therefore the quantum deformation removes no dangerous nonunitary directions
  by itself; the D-term/radial unitary lock removes 15 extra real moduli in
  this linearized comparison.

paper and ledger sync:
  Added the quantum-deformed hidden-link moduli audit to
  paper/gut_framework.tex.

  Added two ledger rows:

    Nf=Nc quantum-deformed hidden moduli are compatible with the unitary link.
    Nf=Nc quantum deformation alone enforces the unitary crossed-120 link.

  The first is PASS_CONDITIONAL; the second is NO_GO.  The conditional theorem
  ledger now has 58 claims:

    FAIL: 1
    NO_GO: 10
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 37
    TUNED_FALLBACK: 3

  The TeX draft was compiled twice from the paper directory.  The refreshed
  paper/gut_framework.pdf has 103 pages.  A log scan for Warning, undefined,
  Overfull, Underfull, Error, and Fatal returned no matches.

current obstacle:
  O3 remains open in its strong form.  The quantum-deformed hidden moduli space
  allows the unitary link, but unitarity still comes from the D-term/Kahler or
  radial-lock sector rather than from the holomorphic quantum constraint.

next attempted nontrivial idea:
  Make the radial/unitarity lock itself dynamical inside the hidden model.
  The least artificial route is a hidden vectorlike pair plus singlet/radial
  driver sector that fixes Q^\dagger Q and Qtilde Qtilde^\dagger without
  adding visible thresholds.  If that cannot be made holomorphic enough, move
  to the N_f=8 conformal-window branch and use anomalous dimensions to make
  the unitary-link constraint an IR-attractive fixed manifold.

proposed next script:

    code/construct_hidden_radial_lock_sector.py

verification plan:
  1. Write the minimal hidden-source/radial driver variables and their charge
     assignments.
  2. Linearize F- and D-flatness around Q=fI, Qtilde=fB and compute whether
     radial/nonunitary modes get masses.
  3. Check that the remaining light modes are only unitary link moduli plus
     hidden singlets.
  4. Confirm again that the induced visible threshold vector is zero and that
     the hidden beta function remains safe to R=200.

## 2026-05-09 01:22 Asia/Taipei heartbeat -- hidden radial-lock D-term sector

status:
  The task is still incomplete as an unconditional first-principles GUT, but
  O3 has been narrowed again.  The missing unitary lock can be realized by a
  hidden radial/D-term moose at the finite-rank component level.  This is the
  cleanest current A5-A6 closure: the remaining caveat is the fully propagating
  endpoint gauge/vectorlike completion, not the link geometry itself.

new audit:
  Added:

    code/construct_hidden_radial_lock_sector.py

  Outputs:

    output/hidden_radial_lock_sector/summary.json
    output/hidden_radial_lock_sector/report.md
    output/hidden_radial_lock_sector/finite_radial_lock_samples.csv
    output/hidden_radial_lock_sector/rank_audit.csv
    output/hidden_radial_lock_sector/charge_table.csv

mathematical construction:
  Use hidden bifundamental matrices Q and Qtilde with link meson

    B = Qtilde Q / f^2.

  The radial/hidden moment maps are

    mu_L = Q Q^dagger - f^2 I_4,
    mu_R = Qtilde^dagger Qtilde - f^2 I_4,
    mu_H = Q^dagger Q - Qtilde Qtilde^dagger.

  On the branch

    Q = f I_4,
    Qtilde = f B,
    B in U(4),

  all moment maps vanish.  The hidden U(4)_H quotient removes the internal
  frame, leaving B as the physical constrained source/link coordinate.

rank and Hessian audit:
  Linearizing on the 64 real components of (Q,Qtilde):

    radial D-term rank = 32
    hidden moment rank = 16
    combined D rank = 32
    combined D nullity before quotient = 32
    hidden gauge orbit rank = 16
    residual after radial D and hidden quotient = 16

  The residual 16 real dimensions are exactly dim_R U(4).  Imposing the
  visible fixed block gives:

    fixed-block rank on residual link = 7
    residual unitary completion moduli = 9

finite numerical verification:
  For the same unitary-completion angle grid:

    all_finite_samples_Dflat_unitary_and_locked = True
    max D residual = 2.220797e-15
    max |log singular| = 1.221245e-15
    completed-partner lock window = 2.216814e-3
    visible threshold vector = (0,0,0)
    hidden SU(4)_H beta b = -8

paper and ledger sync:
  Added the radial D-term moose audit to paper/gut_framework.tex.

  Added one ledger row:

    Hidden radial D-term moose dynamically realizes the crossed-120 unitary lock.

  The conditional theorem ledger now has 59 claims:

    FAIL: 1
    NO_GO: 10
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 38
    TUNED_FALLBACK: 3

  The TeX draft was compiled twice from the paper directory.  The refreshed
  paper/gut_framework.pdf has 104 pages.  A log scan for Warning, undefined,
  Overfull, Underfull, Error, and Fatal returned no matches.

current obstacle:
  O3 is now reduced to a UV-completion caveat: if endpoint copy groups are made
  fully propagating gauge groups, one must add vectorlike/anomaly-completing
  hidden singlet matter and verify that it remains threshold-silent and safe
  to R=200.  The local D-term geometry itself passes.

next attempted nontrivial idea:
  Build the vectorlike endpoint-gauge completion card.  The goal is not to add
  new visible matter, but to show that any endpoint-gauge anomaly completion is
  entirely visible-singlet, vectorlike, and either asymptotically free or
  safely below the R=200 Landau target.

proposed next script:

    code/audit_endpoint_vectorlike_completion.py

verification plan:
  1. List the endpoint U(4)_L x U(4)_R x U(4)_H charges for Q,Qtilde and the
     required conjugate/vectorlike partners.
  2. Compute cubic anomaly proxies for each hidden copy gauge factor.
  3. Compute one-loop hidden beta coefficients for the minimal vectorlike
     completion.
  4. Confirm visible Spin(10) threshold vector remains zero, or classify any
     nonzero visible charge as a no-go.

## 2026-05-09 01:31 Asia/Taipei heartbeat -- endpoint vectorlike completion audit

status:
  The task remains incomplete as a full first-principles GUT, but the A5-A6
  hidden unitary-link sector is now conditionally closed at the finite-rank
  action-bookkeeping level.  Raw endpoint gauging is anomalous; the minimal
  vectorlike endpoint completion cancels the anomalies, is one-loop beta-safe
  to R=200, and remains visible-threshold silent.

new audit:
  Added:

    code/audit_endpoint_vectorlike_completion.py

  Outputs:

    output/endpoint_vectorlike_completion/summary.json
    output/endpoint_vectorlike_completion/report.md
    output/endpoint_vectorlike_completion/field_content.csv
    output/endpoint_vectorlike_completion/anomaly_audit.csv
    output/endpoint_vectorlike_completion/beta_audit.csv
    output/endpoint_vectorlike_completion/visible_threshold_audit.csv

field content:
  Raw hidden radial-link fields:

    Q       ~ (4_L, bar4_H, 1)
    Qtilde  ~ (1, 4_H, bar4_R)

  Vectorlike completion:

    Qc        ~ (bar4_L, 4_H, 1)
    Qtilde_c  ~ (1, bar4_H, 4_R)

  All fields are visible Spin(10) singlets.

anomaly audit:
  Raw endpoint gauging gives cubic anomaly proxies:

    A_L = +4
    A_H = 0
    A_R = -4

  Therefore raw endpoint gauging is a no-go.  With Qc and Qtilde_c included:

    A_L = A_H = A_R = 0.

beta and threshold audit:
  Using

    b = sum_R T(R) - 3 C_2(SU(4)),
    T(4)=1/2,
    C_2(adj)=4,

  the vectorlike-completed hidden beta coefficients are:

    b_L = -8
    b_H = -4
    b_R = -8

  With convention d alpha^{-1}/d log mu = -b/(2pi), there is no one-loop UV
  Landau obstruction before R=200.  Since all endpoint-completion fields are
  visible singlets:

    visible threshold vector = (0,0,0).

paper and ledger sync:
  Added the endpoint vectorlike completion audit to paper/gut_framework.tex.

  Added two ledger rows:

    Raw endpoint gauging of the hidden radial link is anomaly-free.
    Vectorlike endpoint completion is anomaly-free, beta-safe, and threshold-silent.

  The first is NO_GO; the second is PASS_CONDITIONAL.  The conditional theorem
  ledger now has 61 claims:

    FAIL: 1
    NO_GO: 11
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 39
    TUNED_FALLBACK: 3

  The TeX draft was compiled twice from the paper directory.  The refreshed
  PDF has 104 pages, and the log scan found no Warning/undefined/Overfull/
  Underfull/Error/Fatal lines.

current obstacle:
  The A5-A6 hidden unitary-link sector is now conditionally consistent at the
  local action/rank/anomaly/beta/threshold level.  The main remaining blockers
  are therefore again phenomenological and reproducibility-side:

    O1. full joint CKM/PMNS/flavor fit;
    O2. channel-specific dressed d=5 proton decay with physical flavor,
        soft spectrum, and lattice/chiral inputs.

  A residual theoretical caveat remains: the global nonperturbative hidden
  dynamics beyond the finite-rank D-term/vectorlike completion card is not a
  complete UV theorem.

next attempted nontrivial idea:
  Return to the phenomenological bottleneck and build a reproducible full
  flavor + d=5 pipeline around the now-stabilized hidden-link sector.  The
  most efficient route is not another toy Wilson proxy, but a single benchmark
  card carrying Y_u,Y_d,Y_e,Y_nu, CKM/PMNS observables, triplet Wilson tensors,
  soft dressing inputs, and Knu/e pi widths.

proposed next script:

    code/build_full_flavor_d5_pipeline.py

verification plan:
  1. Load the current best two-kernel/shadow/contact flavor cards and exact
     unitary-link triplet source.
  2. Fit or score all SM charged masses, CKM, PMNS, and neutrino splittings in
     one table.
  3. Construct C_5L and C_5R in the same field basis, rotate to mass basis,
     and apply a documented soft-dressing card.
  4. Export a supplemental reproducibility package: matrices, Wilson tensors,
     channel widths, thresholds, and all input constants.

## 2026-05-09 02:02 Asia/Taipei heartbeat -- full flavor plus dressed d=5 closure pipeline

status:
  The task remains incomplete as an unconditional GUT.  The hidden-link
  endpoint/vectorlike sector is locally closed, so the bottleneck is now
  decisively phenomenological: full flavor plus dressed dimension-five proton
  decay.

new closure pipeline:
  Added:

    code/build_full_flavor_d5_pipeline.py

  Outputs:

    output/full_flavor_d5_pipeline/summary.json
    output/full_flavor_d5_pipeline/report.md
    output/full_flavor_d5_pipeline/candidate_scoreboard.csv
    output/full_flavor_d5_pipeline/component_status.csv
    output/full_flavor_d5_pipeline/input_manifest.csv

mathematical bridge:
  The pipeline uses the same mass-basis Wilson definitions as the paper:

    C5L^{abcd} = sum_AB (Y_QQ^A)_{ij} W_AB (Y_QL^B)_{kl}
                 U_Q^{ia} U_Q^{jb} U_Q^{kc} U_L^{ld}

    C6_ch = D_ch C5_ch,
    Gamma_ch = K_ch |C6_ch|^2.

  Hence a lifetime margin mu=tau/tau_target below one requires a common
  amplitude suppression mu^{-1/2}.  This gives a direct numerical bridge
  between the flavor/Wilson score and the proton lifetime target.

numerical result:
  The exact reproducibility card passes exact-input and PMNS checks, and the
  current Knu bound passes with margin 1.113033.  However:

    baseline CKM score = 5.449970e+00
    operator strict closures = 0
    operator loose closures = 0
    best operator CKM score = 1.041311e-03
    best operator future Knu margin = 3.200720e-01
    best operator suppression still needed = 1.767568
    publication-level closure = False

  Therefore the current branch has a reproducible no-closure certificate for
  simultaneous full flavor and future-stress dressed d=5 proton safety.  This
  is not a refutation of the conditional EFT branch; it is a clean diagnosis
  that the current CP1 two-kernel/finite-Clebsch profiles are not enough.

paper and ledger sync:
  Added a new TeX proposition, "Full flavor plus dressed d=5 closure
  pipeline", recording the suppression rule and candidate scoreboard.

  Added a theorem-ledger row:

    Full flavor plus dressed d=5 pipeline is publication-complete.

  Its status is NO_GO.  The conditional theorem ledger now has 62 claims:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 39
    TUNED_FALLBACK: 3

  The TeX draft compiled twice from the paper directory after this update.  The
  refreshed PDF has 105 pages, and the log scan found no Warning/undefined/
  Overfull/Underfull/Error/Fatal lines.

current obstacle:
  The bottleneck is now sharper than "do flavor and proton decay": one needs a
  triplet-sector mechanism that lowers the Knu Wilson amplitude by at least a
  factor of roughly 1.77 at the best operator-level flavor point, while also
  nudging the CKM score below 1e-3 and preserving mass/seesaw quality.

next attempted nontrivial idea:
  Use the already discovered antisymmetric 120-like QQ Knu null and crossed
  120A/120B LLLL+RRRR null as the next constructive direction.  The next useful
  script should not add arbitrary triplet matrices; it should build a
  flavor-compatible 10+126bar+120+120' Clebsch card whose triplet subblock lies
  near the component nullspace while its doublet subblock keeps the best CKM
  point.

proposed next script:

    code/fit_clebsch_null_flavor_d5.py

verification plan:
  1. Start from the best operator CKM row in
     output/full_flavor_d5_pipeline/candidate_scoreboard.csv.
  2. Parameterize triplet Clebsch angles around the antisymmetric QQ and crossed
     120A/120B null directions, with the doublet Yukawa sector held to the
     existing flavor score.
  3. Recompute C5L and C5R in the same mass basis; require both LLLL Knu and
     RRRR uusd channels to improve.
  4. Accept only rows with CKM score < 1e-3, mass score <= 0.2, seesaw residual
     < 1e-10, and future Knu margin > 1.

## 2026-05-09 08:06 Asia/Taipei heartbeat -- Clebsch-null flavor/d5 bridge

status:
  The task remains incomplete as a first-principles GUT.  The full flavor plus
  dressed d=5 closure pipeline gave a no-closure result for finite
  CP1/two-kernel Clebsch profiles, so this heartbeat tested the more
  structured crossed-120 nullspace target.

new bridge script:
  Added:

    code/fit_clebsch_null_flavor_d5.py

  Outputs:

    output/clebsch_null_flavor_d5/summary.json
    output/clebsch_null_flavor_d5/report.md
    output/clebsch_null_flavor_d5/closure_rows.csv

mathematical idea:
  Use the component-level identity

    Y_QQ^T pure antisymmetric 120-like => sym(Y_QQ^T)=0

  so the LLLL Knu tensor vanishes.  Then choose the 10-5bar triplet source in
  the crossed 120_B-like direction to null the monitored RRRR uusd tensor.  A
  finite mass lift

    W_kappa = U diag(1,1/kappa,1/kappa,1/kappa) V^dagger

  leaks only at O(1/kappa).

numerical result:
  Ordinary finite doublet-reused Clebsch profiles still fail future d=5:

    ordinary best CKM future Knu margin = 2.986321e-01
    ordinary source-consistent future Knu margin = 3.200720e-01

  Accepting the crossed 120A/120B post-Spin(10)/constrained projector gives:

    kappa=30 future Knu margin  = 1.668117e+10
    kappa=100 future Knu margin = 1.853464e+11
    loose conditional closures = 2
    strict conditional closures = 0

  Therefore the crossed projector conditionally closes the dressed d=5 gap,
  but does not close the paper.  The source-consistent CKM score remains

    chi_CKM^2 = 1.453203e-03,

  so it misses the strict 1e-3 target by a factor 1.453203.

selection-rule caveat:
  The crossed projector is not enforceable by a field-only unbroken Spin(10)
  grading once both 16 16 120_A and 16 16 120_B Yukawa terms are allowed.  It
  must be treated as a post-Spin(10)-breaking Pati-Salam/constrained
  triplet-sector projector unless an explicit stronger Higgs/source symmetry is
  constructed.

paper and ledger sync:
  Added a new TeX proposition, "Clebsch-null flavor/d=5 bridge".

  Added a theorem-ledger row:

    Crossed-120 Clebsch-null bridge closes the d=5 gap without closing strict flavor.

  Its status is PASS_CONDITIONAL.  The conditional theorem ledger now has 63
  claims:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 40
    TUNED_FALLBACK: 3

  TeX compiled cleanly after the sync pass:

    pdflatex x2 from paper/
    PDF pages: 106
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean

current obstacle:
  The proton-side gap is now conditionally closed if the crossed 120 projector
  is accepted, so the sharp remaining task is to lower the source-consistent
  CKM score by about 31.2 percent, from 1.453203e-3 to below 1e-3, without
  leaving the crossed projector source label or spoiling the mass score.

next attempted nontrivial idea:
  Fit a small doublet-only deformation around the source-consistent
  d5_both_F_minus_0 row while freezing the crossed triplet projector.  This
  should be framed as a 10+126bar+120+120' doublet-Higgs mixing correction,
  not as a free triplet Wilson retuning.

proposed next script:

    code/fit_source_consistent_ckm_with_crossed120.py

verification plan:
  1. Start from the d5_both_F_minus_0 matrix card.
  2. Allow only doublet-sector Clebsch/mixing perturbations that preserve the
     crossed 120A/120B triplet source assignment.
  3. Require CKM score < 1e-3, mass score <= 0.2, seesaw residual < 1e-10.
  4. Reuse the crossed-projector kappa=30 leakage row for d=5 and verify that
     the projected threshold remains zero in the complete-degenerate branch.

## 2026-05-09 09:11 Asia/Taipei heartbeat -- Source-consistent CKM repair with crossed-120 frozen

status:
  The full first-principles GUT task is still incomplete, but the previous
  strict flavor-plus-d5 numerical bottleneck has been locally closed under the
  already stated crossed-120 conditional assumption.

new script:

    code/fit_source_consistent_ckm_with_crossed120.py

outputs:

    output/source_consistent_ckm_crossed120/summary.json
    output/source_consistent_ckm_crossed120/report.md
    output/source_consistent_ckm_crossed120/scan_rows.csv

mathematical idea:
  Freeze the source-consistent CP1/operator tensors H,F,G_A,G_B from the
  d5_both_F_minus_0 card and freeze the crossed 120_A/120_B triplet projector
  at kappa=30.  Vary only the doublet-Higgs mixing coefficients

    r_u,r_d,a_u,b_u,a_d,b_d,a_e,b_e,a_nu,b_nu.

  This tests whether the CKM miss is a doublet-mixing artifact rather than a
  proton-side Wilson obstruction.  The d5 margin is not retuned.

numerical result:
  The source seed was

    CKM score  = 1.453203e-03
    mass score = 1.602662e-01
    crossed d5 future margin = 1.668117e10
    seesaw residual = 5.859e-12.

  The best mass-valid local doublet refit gives

    CKM score  = 7.634719e-04
    mass score = 1.490597e-01
    crossed d5 future margin = 1.668117e10
    seesaw residual = 5.174e-12
    ||delta x||/||x0|| = 4.995e-02.

  There are 3 sampled strict conditional closures satisfying

    CKM < 1e-3,
    mass score < 0.2,
    seesaw residual < 1e-10,
    d5 future margin > 1.

interpretation:
  This is a local conditional closure, not a publication-level or
  first-principles completion.  It still relies on the crossed 120_A/120_B
  projector, which is forbidden by a field-only unbroken Spin(10) grading and
  must be derived from a post-Spin(10)-breaking constrained/source mechanism.

paper and ledger sync:
  Added a TeX proposition,

    Source-consistent CKM repair with crossed 120 frozen.

  Added a theorem-ledger row:

    Source-consistent doublet-only CKM repair closes strict flavor with frozen
    crossed-120 d=5 projector.

  The theorem ledger now has 64 claims:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 41
    TUNED_FALLBACK: 3

  TeX compiled cleanly after the sync pass:

    pdflatex x2 from paper/
    PDF pages: 106
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean

current obstacle:
  The numerical flavor+d5 bottleneck is now conditionally closed.  The sharp
  remaining physics bottleneck is action-level derivation of the crossed
  120_A/120_B triplet projector, compatible with the doublet-only mixing refit,
  threshold silence, and the existing source/constrained-sector UV logic.

next attempted nontrivial idea:
  Construct an explicit post-Spin(10) Pati-Salam constrained-source
  superpotential that enforces crossed triplet pairing while leaving the
  doublet sector free to realize the new local CKM refit.  The candidate
  should separate triplet source labels from doublet Higgs mixing by a
  source/link grading rather than by a forbidden field-only Spin(10) charge.

proposed next script:

    code/audit_crossed120_action_level_source_symmetry.py

verification plan:
  1. Define PS-fragment source charges for 120_A and 120_B triplet fragments,
     doublet fragments, and link/driver fields.
  2. Require the crossed triplet mass block allowed, same-source triplet blocks
     forbidden or lifted, and all doublet-mixing operators needed by the
     CKM-refit card allowed.
  3. Build the component Hessian and verify the triplet projector rank,
     threshold vector, and finite-leakage d5 margin.
  4. Replay the local CKM-refit Yukawa matrices with the derived source
     selection rules and check that the theorem-ledger PASS_CONDITIONAL row is
     not merely imposed by hand.

## 2026-05-09 09:41 Asia/Taipei heartbeat -- crossed-120 action-level source symmetry

status:
  The task is still incomplete as a first-principles GUT derivation.  This
  heartbeat upgrades the crossed 120_A/120_B flavor+d5 repair from a naked
  imposed projector to a conditional post-Spin(10)-breaking Pati-Salam
  source/spurion selection rule.

new artifact:
  Added

    code/audit_crossed120_action_level_source_symmetry.py

  with outputs

    output/crossed120_action_level_source_symmetry/summary.json
    output/crossed120_action_level_source_symmetry/report.md
    output/crossed120_action_level_source_symmetry/operator_ledger.csv
    output/crossed120_action_level_source_symmetry/charge_table.csv.

mathematical idea:
  Use a post-breaking source grading

    U(1)_R x Z4_src x Z2_frag x Z2_inert

  acting on PS fragments/sources rather than on a single unbroken Spin(10)
  multiplet.  The allowed operator classes are

    Y_QQ_A B_QQ T_A,
    Y_QL_B B_QL T_B,
    Y_D_A B_16_16 D_A,
    Y_D_B B_16_16 D_B,
    S_AB Tbar_A T_B,
    S_BA Tbar_B T_A,
    S_lock L Lbar.

  The forbidden classes include same-source triplet bare masses, wrong
  crossed Yukawas, inert Yukawa leakage, and triplet-doublet off-block mixing:

    Tbar_A T_A,
    Tbar_B T_B,
    Y_QQ_B B_QQ T_B,
    Y_QL_A B_QL T_A,
    B_16_16 L,
    S_lock D_A Lbar,
    S_lock T_A D_A.

numerical verification:
  The operator ledger passes all fourteen audited classes:

    operator_ledger_passes = true.

  The local source-consistent CKM+d5 closure remains intact:

    chi2_CKM  = 7.634719e-04
    chi2_mass = 1.490597e-01
    M_Knu^{10^35 yr} margin = 1.668117e10
    seesaw residual = 5.174e-12.

  The same audit also keeps the previous hard warnings visible:

    literal triplet-only mediator threshold overshoot at R=200 = 270.361...
    inert 5+5bar completion lock window |log(M_D/M_T)| <= 2.216814e-3
    equivalent mass-lock percent window = 0.221927 percent
    holomorphic crossed-link nullity = 24 real moduli.

interpretation:
  This is a real advance, but still conditional.  It shows that a PS
  fragment/source EFT can make the crossed triplet projector, physical
  doublet Yukawa/mixing repair, and inert partner lift mutually compatible.
  It does not yet prove a microscopic unbroken Spin(10) selection rule.
  The remaining constraint is the NLSM/D-term/composite unitarity mechanism
  needed to remove the 24 real holomorphic link moduli without releasing new
  non-universal thresholds.

paper and ledger sync:
  Added a TeX proposition,

    Post-Spin(10) source symmetry for the crossed 120 branch.

  Added the theorem-ledger row,

    Post-Spin(10) source symmetry makes the crossed-120 flavor+d=5 closure
    action-level compatible.

  The theorem ledger now has 65 claims:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 42
    TUNED_FALLBACK: 3.

  TeX compiled cleanly after the sync pass:

    pdflatex x2 from paper/
    PDF pages: 107
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean.

current obstacle:
  The immediate bottleneck is no longer merely "invent a crossed 120
  projector"; it is to derive the unitary/constrained source link dynamics
  that fixes AB=1 and PBP=W_kappa without leaving the 24 real relative
  orientation moduli as physical flat directions or threshold-carrying modes.
  The broader publication bottleneck remains the same: full flavor/PMNS fit,
  channel-specific d=5 proton decay, and supplemental reproducibility.

next attempted nontrivial idea:
  Build the crossed 120 link as a constrained or composite unitary sigma-model
  variable rather than as a generic holomorphic matrix.  The intended
  mechanism is that the holomorphic equations AB=1 and PBP=W_kappa set the
  algebraic projector, while D-terms/Kahler curvature/composite unitarity
  remove the 24 real moduli and keep the inert 5+5bar completion in complete,
  threshold-silent multiplets.

proposed next script:

    code/audit_crossed120_unitary_link_stabilization.py

verification plan:
  1. Start from the crossed-link solution of AB=1 and PBP=W_kappa.
  2. Add explicit D-term or NLSM constraints such as A^\dagger A=1, or a
     composite meson constraint with the same tangent-space projection.
  3. Compute the tangent Hessian and verify that all 24 real relative
     orientation modes are lifted or gauge/equivalence directions.
  4. Re-run the crossed source-symmetry ledger, the CKM+d5 local closure, the
     inert-partner mass-lock audit, and the conditional theorem ledger.
  5. If any lifted mode carries a non-universal SM beta vector, feed it back
     into the two-loop RGE/proton scan instead of declaring it harmless.

## 2026-05-09 10:18 Asia/Taipei heartbeat -- reconcile crossed-link unitarity bottleneck

status:
  The first-principles GUT task remains incomplete.  However, this heartbeat
  corrects the immediate bottleneck listed in the previous 09:41 note: the
  dangerous nonunitary crossed-link moduli are already addressed by the stored
  D-term/Kahler, composite GLSM, hidden radial-lock, and endpoint-vectorlike
  audits.  The active bottleneck should therefore shift to publication-level
  closure of full flavor/PMNS, channel-specific d=5 proton decay, and a single
  supplemental reproducibility pipeline.

checks rerun:
  Re-ran

    code/construct_unitary_link_dterm_quotient.py
    code/construct_composite_unitary_link_glsm.py
    code/construct_hidden_radial_lock_sector.py
    code/audit_conditional_theorem_ledger.py

  and inspected

    output/unitary_link_dterm_quotient/summary.json
    output/hidden_radial_lock_sector/summary.json
    output/endpoint_vectorlike_completion/summary.json
    output/full_flavor_d5_pipeline/summary.json.

mathematical reconciliation:
  The holomorphic equations AB=1 and PBP=W_kappa leave 24 real nonunitary
  moduli.  Imposing the D-term/Kahler constraint

    B^\dagger B = I_4

  linearizes to

    B^\dagger delta B + delta B^\dagger B = 0,

  so the tangent space is anti-Hermitian and has real dimension 16.  The
  fixed visible block condition P delta B P=0 has rank 7 on this D-flat
  tangent space, leaving 9 real unitary completion moduli.  These residual
  moduli are harmless for the completed-partner threshold lock because B is
  unitary: all singular values of the mass link remain equal to one.

numerical verification:
  D-term/Kahler quotient:

    real variables in B = 32
    D-term linear rank = 16
    D-flat tangent dimension = 16
    fixed-block rank on D-flat tangent = 7
    residual unitary completion moduli = 9
    max finite-completion unitarity residual = 1.785e-15
    max |log singular value| in finite samples = 6.661e-16
    completed-partner lock window = 2.216814e-3

  Hidden radial-lock realization:

    real variables in Q,Qtilde = 64
    radial D-term rank = 32
    hidden gauge orbit rank = 16
    residual unitary link dimension = 16
    residual after fixed block = 9
    max D residual = 2.221e-15
    max |log singular value| = 1.221e-15
    visible threshold vector = (0,0,0).

  Endpoint/vectorlike completion:

    raw endpoint gauging is anomalous in SU(4)_L and SU(4)_R;
    adding Qc and Qtilde_c cancels all hidden cubic anomaly proxies;
    b_L=-8, b_H=-4, b_R=-8 for the vectorlike completion;
    vectorlike completion is one-loop safe to R=200;
    visible threshold vector remains (0,0,0).

  Conditional theorem ledger remains unchanged:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 42
    TUNED_FALLBACK: 3.

TeX sync:
  Updated the post-Spin(10) crossed-120 source-symmetry proposition in
  paper/gut_framework.tex so it no longer leaves the 24-moduli statement as
  the active bottleneck.  It now points to the later unitary-link audits:
  the D-term quotient removes the nonunitary directions, the fixed-block
  driver has rank 7 on the 16-real-dimensional unitary link, and finite
  completions keep the completed-partner mass splitting at machine precision.

  Recompiled

    pdflatex -interaction=nonstopmode -halt-on-error gut_framework.tex

  successfully.  The PDF has 107 pages.  A log scan for
  Warning/undefined/Overfull/Underfull/Error/Fatal is clean.

current obstacle:
  The crossed-link unitarity issue is now locally conditionally controlled.
  What is not finished is the global phenomenology/reproducibility closure.
  The older full_flavor_d5_pipeline still reports

    baseline CKM completed = false
    d5 current bound passes = true
    d5 future 1e35 passes = false
    publication_level_complete = false.

  Later source-consistent crossed-120 refits give a strong local closure, but
  the paper still needs a single unified pipeline that merges:

    full mass/CKM/PMNS observables,
    the source-consistent crossed-120 flavor card,
    channel-specific d=5 Wilson tensors and dressing,
    d=6 proton bounds,
    threshold/RGE cards,
    exact input manifests and hashes.

next attempted nontrivial idea:
  Promote the local source-consistent crossed-120 refit into the full
  publication pipeline rather than treating it as a separate local audit.
  The creative step is to make one benchmark card the single source of truth:
  it should contain the Yukawa matrices, PMNS/CKM rotations, crossed triplet
  inverse block, unitary-link/radial-lock card, inert partner lock, and all
  channel-specific proton-decay constants.

proposed next scripts:

    code/build_publication_closure_card.py
    code/audit_publication_flavor_d5_reproducibility.py

verification plan:
  1. Merge output/source_consistent_ckm_crossed120 with the exact matrix card
     in output/locked_link_full_flavor_d5_card.
  2. Recompute masses, CKM, PMNS, seesaw residual, and Majorana spectrum from
     only that merged card.
  3. Recompute C5L and C5R channel-by-channel for K+ nu, e+ pi0, mu+ pi0,
     and K0 mu+ using the same rotations and dressing inputs.
  4. Require present proton limits and the 1e35 yr stress target to be clearly
     labelled as separate gates.
  5. Emit a manifest with hashes for every input and output used by the paper.
  6. Update the theorem ledger so the remaining OPEN items are either closed
     by this card or explicitly downgraded to future-work assumptions.

## 2026-05-09 10:48 Asia/Taipei heartbeat -- publication closure candidate card

status:
  The first-principles GUT task is still incomplete.  This heartbeat closes a
  narrower bookkeeping bottleneck: the source-consistent crossed-120 local
  flavor+d=5 branch is now represented by a single reproducible publication
  candidate card.  It passes the local gates, but it is not yet a
  publication-final package.

new artifacts:
  Added

    code/build_publication_closure_card.py
    code/audit_publication_flavor_d5_reproducibility.py

  with outputs

    output/publication_closure_card/summary.json
    output/publication_closure_card/publication_closure_card.json
    output/publication_closure_card/report.md
    output/publication_closure_card/input_manifest.csv
    output/publication_flavor_d5_reproducibility/summary.json
    output/publication_flavor_d5_reproducibility/report.md
    output/publication_flavor_d5_reproducibility/gate_status.csv
    output/publication_flavor_d5_reproducibility/recompute_differences.csv.

mathematical setup:
  The card selects the source-consistent crossed-120 row

    all10_radius_0p30_ckmheavy_1.

  It recomputes

    V_CKM = U_u^\dagger U_d

  from the left rotations of Y_u Y_u^\dagger and Y_d Y_d^\dagger, recomputes
  sector singular-value ratios, and replays the inverse type-I seesaw

    M_R = -m_D^T m_nu^{-1} m_D

  against the fixed PMNS benchmark.  The d=5 proton gate is the finite-lift
  crossed-120 kappa=30 worst Knu margin against the 1e35 yr stress target.

numerical verification:
  From the card's own Yukawa matrices:

    CKM score       = 7.634719e-04
    mass score      = 1.490597e-01
    seesaw residual = 5.173996e-12
    future d5 margin at 1e35 yr = 1.668117e10.

  CKM observables:

    |V_us| = 0.2349862183
    |V_cb| = 0.0429226120
    |V_ub| = 0.0036712892
    J      = 3.579814639e-05.

  Heavy Majorana masses from the replay:

    M1 = 2.380447592e10 GeV
    M2 = 3.220379683e13 GeV
    M3 = 3.926531008e15 GeV.

  Recompute differences from the stored card are exactly zero at printed
  precision for

    CKM score, mass score, seesaw residual, |V_us|, |V_cb|, |V_ub|.

local gates:
  The reproducibility audit returns true for

    card recomputes within tolerance,
    strict CKM,
    loose mass,
    seesaw residual,
    future 1e35 d5 margin,
    post-Spin(10) source symmetry,
    D-term unitary lock,
    hidden radial lock,
    endpoint vectorlike safety.

paper and theorem-ledger sync:
  Added a TeX proposition:

    Single-card source-consistent flavor and d=5 candidate.

  Updated the theorem ledger with a PASS_CONDITIONAL row:

    Single-card source-consistent crossed-120 flavor+d=5 candidate is reproducible.

  The theorem ledger now has 66 claims:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 43
    TUNED_FALLBACK: 3.

  TeX compiled successfully:

    gut_framework.pdf has 107 pages
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean.

interpretation:
  The earlier full_flavor_d5_pipeline remains useful as a legacy obstruction
  record, but it no longer describes the strongest local branch.  The current
  strongest branch passes local source-consistent CKM, mass, seesaw, source
  symmetry, unitary-link, endpoint, and d=5 stress gates from one card.

  This is not an unconditional or publication-final theorem.  The crossed
  projector is still a post-Spin(10) source-sector mechanism, not a field-only
  unbroken Spin(10) derivation.

current obstacle:
  The remaining active bottleneck is final publication closure:

    1. replace internal flavor targets by a final cited input table;
    2. regenerate all channel-specific d=5 widths from this exact card;
    3. emit final tables for K+nu, e+pi0, mu+pi0, and K0mu+;
    4. lock TeX tables to generated JSON/CSV artifact hashes.

next attempted nontrivial idea:
  Build the channel-specific proton-decay replay directly from

    output/publication_closure_card/publication_closure_card.json

  rather than from older Wilson proxy outputs.  The new replay should treat
  the card as the single source of truth for Yukawas, rotations, crossed
  inverse block, and seesaw data.

proposed next script:

    code/audit_publication_channel_d5_tables.py

verification plan:
  1. Read the publication closure card and recompute all left/right flavor
     rotations.
  2. Build C5L and C5R tensors for K+nu, e+pi0, mu+pi0, and K0mu+ using the
     crossed-120 inverse block in the same card.
  3. Apply the existing dressing constants and width normalization, but record
     them explicitly in the manifest.
  4. Emit channel-by-channel lifetime margins for current bounds and the 1e35
     yr stress target.
  5. If any channel fails, identify whether the failure is a flavor rotation,
     dressing, or triplet-filter issue rather than hiding it in a scalar Knu
     proxy.

## 2026-05-09 11:24 Asia/Taipei heartbeat -- publication-card d=5 channel scaffold

status:
  The full first-principles GUT task remains incomplete.  This heartbeat makes
  the proton-side publication bottleneck sharper: the single publication card
  can now generate channel-level d=5 rows, but the generated table is still a
  proxy/class-level scaffold, not an exact dressed proton-decay proof.

new artifacts:
  Updated

    code/build_publication_closure_card.py

  so the card also embeds the source-consistent triplet tensor data
  \(H,F,G_A,G_B\) from the two-kernel source model.  Added

    code/audit_publication_channel_d5_tables.py

  with outputs

    output/publication_channel_d5_tables/summary.json
    output/publication_channel_d5_tables/channel_table.csv
    output/publication_channel_d5_tables/report.md
    output/publication_channel_d5_tables/input_manifest.csv.

mathematical setup:
  The channel scaffold rebuilds the flavor rotations from the publication
  card and evaluates

    C_L^{abcd}
      = Y_QQ^{ij} Y_QL^{kl}
        U_{q1}^{ia} U_{q2}^{jb} U_{q3}^{kc} U_l^{ld},

    C_R^{abcd}
      = Y_UE^{ij} Y_UD^{kl}
        U_u^{ia} U_e^{jd} U_u^{kb} U_d^{lc}.

  It compares two hypotheses:

    finite_Clebsch_profile_F_quarter_phase,
    crossed_pure_120_source_proxy.

  The crossed pure 120 source has exact antisymmetric residuals

    ||G_A + G_A^T||_F = ||G_B + G_B^T||_F = 0,

  with normalized largest singular values

    sigma_max(Y_QQ) = 0.6,
    sigma_max(Y_QL) = 0.024.

numerical verification:
  The scaffold emits 14 rows covering

    K+ nu_bar, K0 mu+, e+ pi0, mu+ pi0

  across LLLL and RRRR monitors.  The crossed kappa=30 class-level finite-lift
  margins are

    LLLL margin at 1e35 yr = 1.668117e10,
    RRRR margin at 1e35 yr = 3.218681e11.

  Representative raw central proxy amplitudes are

    crossed pure LLLL K+nu      |C| = 0, tau(ST=1) = infinity,
    finite Clebsch LLLL K+nu    |C| = 1.942324e-3, tau(ST=1) = 1.380314e28 yr,
    finite Clebsch RRRR e+pi0   |C| = 5.398852e-8, tau(ST=1) = 6.239844e36 yr,
    finite Clebsch RRRR mu+pi0  |C| = 2.394907e-7, tau(ST=1) = 3.171023e35 yr.

paper and theorem-ledger sync:
  Added a TeX proposition:

    Publication-card d=5 channel scaffold.

  Regenerated the theorem ledger with a new PASS_CONDITIONAL row:

    Publication-card d=5 channel table scaffold is generated.

  The theorem ledger now has 67 claims:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 44
    TUNED_FALLBACK: 3.

  TeX compiled successfully with pdflatex x2:

    paper/gut_framework.pdf has 108 pages,
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean.

current obstacle:
  The d=5 publication item is still open.  The scaffold explicitly reports

    exact_channel_wilson_complete = false.

  Missing before a final theorem:

    1. triplet mass-eigenstate mixing matrices,
    2. short-distance and long-distance SUSY dressing matrices by channel,
    3. final C5L and C5R tensor arrays after the crossed inverse block,
    4. channel-specific chiral/lattice reduction operators with cited
       uncertainty propagation.

next attempted nontrivial idea:
  Promote the crossed inverse block from a class-level leakage ratio to an
  explicit triplet mass-eigenstate map.  In practice, extend the publication
  card with a tensor-valued `triplet_eigenstate_card`:

    M_T, U_T, V_T, Lambda_T^{-1},
    C5L^{abcd}, C5R^{abcd}

  so every channel table row is computed from the same diagonalized triplet
  sector rather than inheriting a class-level margin.

verification plan:
  1. Build the triplet mass matrix in the source basis
     (120_A,120_B,completed partners).
  2. Diagonalize it and export the inverse block
     \(Y_L^T M_T^{-1} Y_R\) as exact four-index tensors.
  3. Replay K+nu, K0mu, e+pi0, and mu+pi0 with explicit wino/higgsino
     dressing factors.
  4. Replace the `PROXY_ONLY_CLASS_LEVEL_PASS` labels by exact channel gates,
     or keep the OPEN item if any required tensor remains implicit.

## 2026-05-09 13:38 Asia/Taipei heartbeat -- triplet eigenstate C5 tensor card

status:
  The full first-principles GUT task remains incomplete, but the d=5 proton
  bottleneck advanced from a class-level leakage table to an explicit
  source-basis tensor card.  The crossed triplet inverse block is now
  diagonalized and contracted into exported C5L/C5R proxy arrays.

important convention correction:
  The publication CKM repair is doublet-only, while the crossed triplet
  projector is frozen from the source-consistent triplet sector.  Therefore
  `audit_publication_channel_d5_tables.py` was corrected to rebuild triplet
  component bases from the triplet-source `model_data.mixing_coefficients`,
  not from the repaired doublet `selected_row.mixing_coefficients`.

new artifact:
  Added

    code/build_publication_triplet_eigenstate_card.py

  with outputs

    output/publication_triplet_eigenstate_card/summary.json
    output/publication_triplet_eigenstate_card/triplet_eigenstate_card.json
    output/publication_triplet_eigenstate_card/channel_table.csv
    output/publication_triplet_eigenstate_card/report.md
    output/publication_triplet_eigenstate_card/input_manifest.csv.

mathematical setup:
  In the source-consistent component basis

    B_10 = (H, r_u F, a_u G_A, b_u G_B),
    B_5bar = r_d (H, F, a_d G_A, b_d G_B),

  the crossed component-null direction is

    p = (0,0,2,0),
    q approximately (-3.09e-7,
                     4.49e-7 + 6.41e-7 i,
                     8.86e-6 + 5.41e-6 i,
                     1.163 + 1.627 i),

  with residual 1.019781e-25.  Completing p and q to unitary frames U,V gives

    W_kappa = U diag(1, 1/kappa, 1/kappa, 1/kappa) V^\dagger,
    kappa = 30.

  The mass singular values are proportional to

    M_lock (1, kappa, kappa, kappa),
    M_lock = 3.887851789e15 GeV,

  so the finite block condition number is exactly 30.

  The exported tensors are

    C5L^{abcd}
      = sum_{I,J} W_{IJ} (B_10^I)^{ij} (B_5bar^J)^{kl}
        U_u^{ia} U_u^{jb} U_d^{kc} U_l^{ld},

    C5R^{abcd}
      = sum_{I,J} W_{IJ} (B_10^I)^{ij} (B_5bar^J)^{kl}
        U_u^{ia} U_e^{jd} U_u^{kb} U_d^{lc}.

numerical verification:
  The rank-one block keeps monitored LLLL channels exactly null and gives
  minimum unfiltered RRRR lifetime

    2.139208e42 yr.

  The finite kappa=30 block has worst unfiltered row

    LLLL mu+ pi0:
      |C5| = 3.949138e-5,
      tau(ST=1) = 1.762753e30 yr,
      S_T^max(1e35 yr) = 4.198515e-3.

  Other limiting finite rows:

    LLLL K+ nu_bar:
      |C5| = 4.156717e-5,
      tau(ST=1) = 4.555538e30 yr,
      S_T^max(1e35 yr) = 6.749472e-3.

    RRRR mu+ pi0:
      |C5| = 6.793427e-7,
      tau(ST=1) = 5.956870e33 yr,
      S_T^max(1e35 yr) = 2.440670e-1.

    RRRR e+ pi0:
      |C5| = 1.314489e-7,
      tau(ST=1) = 1.591044e35 yr,
      S_T^max(1e35 yr) = 1.261366.

  With the already-used reference triplet filter

    S_T = 1e-5,

  the worst finite-block margin becomes

    min tau(ST=1e-5) / 1e35 yr = 1.762753e5.

paper and theorem-ledger sync:
  Added a TeX proposition:

    Publication triplet eigenstate C5 card.

  Updated the conditional theorem ledger with a new PASS_CONDITIONAL row:

    Publication-card triplet eigenstate card exports source-basis C5 tensors.

  The theorem ledger now has 68 claims:

    FAIL: 1
    NO_GO: 12
    OPEN: 3
    PASS: 3
    PASS_AS_KINEMATICS: 1
    PASS_CONDITIONAL: 45
    TUNED_FALLBACK: 3.

  TeX compiled successfully with pdflatex x2:

    paper/gut_framework.pdf has 109 pages,
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean.

current obstacle:
  The d=5 publication item is still open.  We now have explicit source-basis
  C5L/C5R tensors, but not a final physical decay calculation.  Missing:

    1. chargino/neutralino and sfermion dressing matrices acting on these exact
       exported tensors,
    2. short- and long-distance running factors by operator,
    3. channel-specific chiral reductions with uncertainty propagation,
    4. a final present-bound and 1e35 yr stress table sourced only from the
       new eigenstate card.

next attempted nontrivial idea:
  Build `audit_publication_dressed_c5_from_eigenstate_card.py`: read the
  exported C5L/C5R arrays, apply the existing MSSM dressing proxy as an
  explicit matrix operation rather than a scalar channel prefactor, and emit
  final-like channel lifetimes with a clear `DRESSED_PROXY` status.

verification plan:
  1. Parse `triplet_eigenstate_card.json` and reconstruct C5L/C5R exactly.
  2. Apply wino/higgsino dressing matrices to C5L/C5R in the same flavor
     basis used by the publication card.
  3. Replay K+nu, K0mu, e+pi0, and mu+pi0 with central and envelope chiral
     constants.
  4. Compare to the present channel bounds and to the uniform 1e35 yr stress.
  5. Keep the OPEN d=5 item unless every row is sourced from the eigenstate
     card and no scalar leakage shortcut remains.

## 2026-05-09 14:05 Asia/Taipei note -- PSLT/Another Physics index-theorem hint

status:
  A useful external suggestion was reviewed: do not present PSLT/Another
  Physics as an unconditional derivation of a unique GUT.  Instead, isolate the
  strongest rigorous theorem that the documents can support:

    two-center divisor R = p_+ + p_- on CP1
    + protected Spin^c Dolbeault carrier
    + no extra unpaired SM-charged kernels
    => exactly three chiral SM-family copies.

mathematical core:
  The proposed clean theorem is

    L_R = O(R) = O(2),
    D_R = sqrt(2)(dbar_L + dbar_L^\dagger),
    ker D_R^+ = H^0(CP1,O(2)),
    ker D_R^- = H^1(CP1,O(2)).

  Since

    h^0(CP1,O(2)) = 3,
    h^1(CP1,O(2)) = h^0(CP1,O(-4)) = 0,

  the protected internal index is exactly

    ind D_R = 3.

interpretation:
  This is a strong first-principles-style conditional theorem for the family
  count.  It does not prove Spin(10), threshold viability, Yukawa texture,
  seesaw, or d=5 proton safety.  It should be used as a front-end family-index
  theorem, not as a replacement for the existing conditional Spin(10) EFT
  model-building chain.

caveats to preserve:
  1. The theorem assumes CP1, the two-center divisor, and the Spin^c carrier.
     These are axioms/inputs, not derived from empty first principles.
  2. The "ordinary spin gives two modes" contrast is useful but should be
     stated as a convention check: the three-family result depends on the
     canonical Spin^c/Dolbeault carrier, not on an arbitrary spin twist.
  3. The visibility formula, g_N, B_N, and rank-2 kinetics remain EFT/spectral
     layer data; they affect formation/visibility, not the protected family
     count.
  4. The anomaly calculation fixes one SM family representation up to
     normalization and assumptions about Yukawa invariance; it does not by
     itself determine family number.

roadmap action:
  Add a future theorem/proposition section before the Spin(10) branch:

    "Conditional Spin^c family-index theorem from the two-center divisor."

  It should explicitly separate:

    PASS_CONDITIONAL: three-family protected kernel from O(2) on CP1;
    OPEN/ASSUMPTION: derivation of CP1, two-center divisor, and exclusivity
    from deeper PSLT dynamics;
    separate EFT layer: Omega, V_eff, WKB/rank-2 kinetics, B_N visibility.

## 2026-05-09 15:07 Asia/Taipei heartbeat -- index theorem integrated, publication bottleneck unchanged

status:
  The overall first-principles GUT task is still incomplete.  The new
  PSLT/Another Physics insight is useful and should be promoted to a formal
  front-end theorem, but it does not close the remaining publication-level
  Spin(10) EFT checks.

current obstacle:
  The family-count subproblem is now conditionally clean:

    two-center CP1 divisor + protected Spin^c Dolbeault carrier
    + no extra SM-charged unpaired kernels
    => h^0(CP1,O(2)) = 3 and h^1(CP1,O(2)) = 0.

  The remaining blockers are downstream:

    1. insert the conditional Spin^c family-index theorem into TeX with the
       assumptions visibly labeled;
    2. keep CP1/two-center/exclusivity as OPEN inputs rather than hidden
       conclusions;
    3. continue the exact channel-specific d=5 proton-decay calculation from
       the triplet eigenstate C5 tensor card;
    4. keep full flavor and supplemental reproducibility tied to the same
       publication card.

next attempted nontrivial idea:
  Treat the Spin^c theorem as a "protected kernel layer" in the paper: the
  PSLT geometry supplies a conditional index theorem for three copies, while
  the crossed-120/source-consistent Spin(10) EFT supplies the conditional
  phenomenology.  This prevents the paper from overclaiming and gives a sharp
  interface between geometry and model-building.

verification plan:
  1. Add a TeX proposition before the Spin(10) branch:
       "Conditional Spin^c family-index theorem from the two-center divisor."
  2. Include the exact Riemann-Roch/Serre-duality check:
       ker D_R^+ = H^0(CP1,O(2)) = 3,
       ker D_R^- = H^1(CP1,O(2)) = H^0(CP1,O(-4)) = 0.
  3. Add an ordinary-spin contrast:
       K_CP1^{1/2} tensor O(2) = O(1), so the ordinary spin twist would give
       two sections, not three.
  4. Mark the theorem PASS_CONDITIONAL and mark CP1/two-center/Spin^c/no-exotic
     exclusivity as explicit assumptions.
  5. Resume the d=5 publication path with the exact eigenstate C5 tensors and
     no scalar leakage shortcut.

sync result:
  Added the TeX proposition "Conditional Spin^c family-index theorem" to
  paper/gut_framework.tex in the O(2) family-geometry section.  Recompiled
  paper/gut_framework.tex with pdflatex successfully; paper/gut_framework.pdf
  remains 109 pages and the log scan reports no Warning/undefined/Overfull/
  Underfull/Error/Fatal lines.

## 2026-05-09 15:39 Asia/Taipei heartbeat -- dressed C5 replay from eigenstate card

status:
  The task remains incomplete.  The recommended next step was executed: a new
  dressed d=5 replay now starts directly from the publication triplet
  eigenstate C5 tensor card, so the scalar leakage shortcut is no longer the
  active bottleneck.

new artifact:
  Added

    code/audit_publication_dressed_c5_from_eigenstate_card.py

  and generated

    output/publication_dressed_c5_from_eigenstate_card/summary.json
    output/publication_dressed_c5_from_eigenstate_card/report.md
    output/publication_dressed_c5_from_eigenstate_card/dressed_channel_rows.csv
    output/publication_dressed_c5_from_eigenstate_card/dressed_point_rows.csv

mathematical/numerical check:
  The replay reads the exported C5L/C5R tensors from
  output/publication_triplet_eigenstate_card/triplet_eigenstate_card.json and
  applies the existing chargino/neutralino plus positive sfermion-eigenstate
  dressing kernel,

    K^{pq}_{ab,a'b'} =
      sum_{r,s} U^p_{ar} U^{p*}_{a'r} U^q_{bs} U^{q*}_{b's}
      D_{pq}(m_chi,m_sfermion_pr,m_sfermion_qs),

  then evaluates

    C6_proxy = S_T * max |K.C5| / M_T

  for K+ nubar, K0 mu+, e+ pi0, and mu+ pi0 under central/max/min chiral
  width normalizations.

result:
  The rank-one inverse block is extremely safe on the audited stress grid:

    rank_one central unsafe 1e35 rows = 0,
    rank_one max-width unsafe 1e35 rows = 0.

  The finite kappa=30 block is not safe at the reference filter S_T=1e-5 on
  the full soft/chiral stress grid:

    finite central unsafe 1e35 rows = 25 / 7056,
    finite max-width unsafe 1e35 rows = 74 / 7056,
    finite min-width unsafe 1e35 rows = 8 / 7056.

  Worst finite central row:

    channel = p_to_muplus_pi0,
    operator = LLLL,
    scenario = democratic_all_stress,
    epsilon = 0.9,
    spectrum = near_unsafe_20TeV,
    pair = 01,
    selected index = 0001,
    tau(S_T=1e-5) = 5.314364e33 yr,
    margin versus 1e35 yr = 5.314364e-2,
    S_T^max(1e35 yr) = 2.305290e-6.

  Worst finite max-width row has the same channel/operator/soft point and
  requires

    S_T^max(1e35 yr) = 1.026080e-6.

interpretation:
  The publication d=5 obstruction has changed character.  It is no longer
  "missing exact C5 tensors"; those tensors now exist and have a direct dressed
  replay.  The obstruction is a real model-building one: the finite kappa=30
  crossed triplet block plus S_T=1e-5 fails under the democratic-all-stress
  near-unsafe 20 TeV soft spectrum.

TeX/theorem-ledger sync:
  Added a TeX proposition "Dressed replay from the publication C5 eigenstate
  card" after the triplet eigenstate C5 proposition.

  Updated code/audit_conditional_theorem_ledger.py with a new OPEN row:

    Publication-card eigenstate C5 tensors have a dressed channel replay.

  Regenerated output/conditional_theorem_ledger.  New status counts:

    FAIL: 1,
    NO_GO: 12,
    OPEN: 4,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 45,
    TUNED_FALLBACK: 3.

  Recompiled paper/gut_framework.tex successfully after adding the proposition:

    paper/gut_framework.pdf has 110 pages,
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean.

current obstacle:
  Decide how to close the finite-block d=5 failure without hiding it:

    A. tighten the triplet filter to S_T <= 1.0e-6 and rerun threshold/RGE;
    B. impose a physically motivated soft-alignment condition excluding the
       democratic_all_stress near_unsafe_20TeV corner;
    C. promote the rank-one/near-rank-one crossed triplet inverse block to an
       action-level completion, since rank_one is safe but finite kappa=30
       leaks through LLLL mu+ pi0.

next attempted nontrivial idea:
  Prefer Route C first: derive a constrained rank-one or near-rank-one
  crossed-triplet inverse block from the source/link action, then quantify how
  small the three orthogonal singular values must be.  Numerically this is a
  scan in kappa with the new dressed C5 replay, not a scalar leakage scan.

verification plan:
  1. Add a kappa scan to audit_publication_dressed_c5_from_eigenstate_card.py
     or a companion script, replacing finite kappa=30 by kappa in
     {30,50,100,300,1000,3000}.
  2. For each kappa, rerun the dressed eigenstate-card replay and record the
     global S_T^max under central and max-width normalizations.
  3. Identify the smallest kappa for which S_T=1e-5 passes the full stress
     grid, and the smallest kappa for which max-width 1e35 passes.
  4. If kappa must be too large to realize naturally, switch to Route B and
     formalize a soft-alignment assumption with an explicit excluded stress
     norm.
  5. Recompile TeX and regenerate the conditional theorem ledger after the
     kappa/soft-alignment decision.

## 2026-05-09 15:47 Asia/Taipei heartbeat -- dressed C5 kappa scan

status:
  The task remains incomplete, but the d=5 obstruction is now quantified as a
  rank-one-locking requirement.

new artifact:
  Added

    code/scan_publication_dressed_c5_kappa.py

  and generated

    output/publication_dressed_c5_kappa_scan/summary.json
    output/publication_dressed_c5_kappa_scan/report.md
    output/publication_dressed_c5_kappa_scan/kappa_scan.csv

derivation:
  Since C5 is linear in the inverse triplet block, the arbitrary-kappa tensor
  is reconstructed exactly from the two exported endpoints:

    C5(kappa) = C5_rank_one + (30/kappa) (C5_finite30 - C5_rank_one).

  This avoids rebuilding the flavor sector and tests only how much the three
  orthogonal triplet inverse-propagator singular values must be suppressed.

numerical result:
  At S_T=1e-5, the dressed eigenstate-card scan gives:

    kappa = 30:
      central unsafe 1e35 rows = 25,
      max-width unsafe 1e35 rows = 74,
      central global S_T^max = 2.305290e-6,
      max-width global S_T^max = 1.026080e-6.

    kappa = 100:
      central unsafe 1e35 rows = 1,
      max-width unsafe 1e35 rows = 18,
      central global S_T^max = 7.684301e-6,
      max-width global S_T^max = 3.420268e-6.

    kappa = 300:
      central unsafe 1e35 rows = 0,
      max-width unsafe 1e35 rows = 0,
      central global S_T^max = 2.305290e-5,
      max-width global S_T^max = 1.026080e-5.

  Therefore the smallest scanned passing value is

    kappa_min = 300.

interpretation:
  The finite-block d=5 problem can be recast as a precise rank-one-locking
  requirement:

    orthogonal inverse-block singular values <= 1/300

  relative to the protected crossed triplet direction.  This is stronger than
  the previous kappa=30 finite-lift card by a factor of 10, but it is no longer
  an unknown proton-decay gap.

TeX sync:
  Added the proposition "Rank-one locking requirement from the dressed C5 scan"
  to paper/gut_framework.tex after the dressed replay proposition.

  Recompiled paper/gut_framework.tex twice successfully:

    paper/gut_framework.pdf has 111 pages,
    log scan for Warning/undefined/Overfull/Underfull/Error/Fatal: clean.

current obstacle:
  Need an action-level reason for kappa >= 300.  Without it, the branch is a
  tuned fallback.  Three possible closures remain:

    A. derive kappa >= 300 from the crossed triplet source/link superpotential;
    B. use a soft-alignment assumption to remove the democratic_all_stress
       near_unsafe_20TeV corner;
    C. tighten the triplet filter to S_T <= 1e-6 and re-run threshold/RGE.

next attempted nontrivial idea:
  Route A is most structural: construct a rank-one-locking driver where three
  orthogonal triplet eigenstates receive masses enhanced by a small hierarchy
  parameter h^{-1}, with h <= 1/300, while the protected crossed direction
  remains at M_lock.  The same driver must be complete/threshold-silent or its
  threshold vector must be scanned.

verification plan:
  1. Build a minimal source-link superpotential for the triplet inverse block
     with one protected singular direction and three orthogonal lifted
     directions.
  2. Derive the Hessian eigenvalues and show kappa = M_orth/M_parallel.
  3. Check whether kappa >= 300 can arise from an order-one coupling ratio,
     a seesaw-like Schur complement, or a symmetry-protected small parameter.
  4. Feed any new non-complete fields into the threshold/RGE/proton scan.
  5. If Route A is too tuned, switch to Route B and quantify the minimal
     soft-alignment cut needed to restore S_T=1e-5.

## 2026-05-09 16:05 Asia/Taipei heartbeat -- finite-link clockwork locking route

status:
  The task remains incomplete.  The publication-card dressed d=5 bottleneck is
  now sharpened one step further: kappa >= 300 can be generated by a finite
  source/link hierarchy, but only as a conditional completion.

new artifact:
  Added

    code/audit_rank_one_clockwork_locking.py

  and generated

    output/rank_one_clockwork_locking/summary.json
    output/rank_one_clockwork_locking/report.md
    output/rank_one_clockwork_locking/clockwork_scan.csv

derivation:
  Replace a bare 300:1 singlet ratio by a nearest-neighbour link chain

    F_a = X_a - q X_{a+1},        a=0,...,n-1.

  The constraint matrix has a null vector

    v_a proportional to q^{-a}.

  If the protected crossed direction couples to the first site while the
  orthogonal leakage directions couple to the endpoint, then

    epsilon = |v_n/v_0| = q^{-n},
    kappa = epsilon^{-1} = q^n.

  This converts the rank-one-locking requirement into a finite graded-chain
  problem rather than a single large coupling.

numerical result:
  The dressed C5 scan requires kappa >= 300.  The smallest passing chains by q
  include:

    q=2.0, n=9:  kappa=512,  max-width 1e35 margin at S_T=1e-5 = 3.066621.
    q=2.5, n=7:  kappa=610.352, margin = 4.357930.
    q=3.0, n=6:  kappa=729,  margin = 6.216919.
    q=4.0, n=5:  kappa=1024, margin = 12.266483.

  The preferred moderate-coupling benchmark is

    q=3, n=6, kappa=729, epsilon=1.371742e-3.

  The nonzero singular values of the link constraint matrix for this card are

    s_min = 2.143405,
    s_max = 3.925024,
    condition = 1.831209,

  so the link constraint itself does not create an exponentially light tower.

UV/threshold warning:
  If the three orthogonal chains are implemented literally as propagating
  complete 10+10bar pairs, the q=3,n=6 card adds roughly 21 vectorlike pairs,
  giving delta b10 = 42 on top of the minimal Yukawa/Higgs sector.  The one-loop
  Landau-pole ratio estimate then falls to

    Lambda_LP/M_G = 12.641502,

  which does not reach R=200.  Therefore the clockwork branch must be realized
  as complete/degenerate threshold-silent links below the cutoff, or preferably
  as constrained/composite source data, not as a large elementary propagating
  tower.

interpretation:
  Route A is not dead.  A finite graded-chain mechanism gives an order-one
  explanation for kappa >= 300:

    kappa = 3^6 = 729.

  However, this is only a conditional source/link completion.  The next hard
  step is to embed the clockwork link into the existing unitary-link D-term or
  composite/constrained source sector so that the added fields do not reopen the
  UV Landau-pole and threshold problems.

next attempted nontrivial idea:
  Build a Spin(10)-compatible constrained clockwork link:

    W_link = sum_a Lambda_a (X_a - q X_{a+1})

  for the three orthogonal triplet directions, with the protected crossed
  direction left at the original M_lock.  The driver/link fields should either
  be SM singlet constrained fields or come in complete degenerate multiplets.

verification plan:
  1. Write the component Hessian for the q=3,n=6 constrained chain and verify
     that only the intended protected endpoint mode remains light.
  2. Check whether the unitary-link D-term quotient can impose the chain
     without 24-real-moduli leakage.
  3. Compute the actual threshold vector for a literal complete-multiplet
     realization; if nonzero, feed it into the two-loop RGE/proton scan.
  4. If constrained/composite realization is chosen, state explicitly that the
     clockwork chain is a boundary/source sector rather than an elementary
     perturbative Spin(10) tower.

## 2026-05-09 16:57 Asia/Taipei heartbeat -- constrained clockwork Hessian audit

status:
  The task remains incomplete, but Route A has gained a component-level
  Hessian check.  The q=3,n=6 clockwork lock can be made source-constrained
  without introducing an exponentially light tower.

new artifact:
  Added

    code/audit_constrained_clockwork_source_hessian.py

  and generated

    output/constrained_clockwork_source_hessian/summary.json
    output/constrained_clockwork_source_hessian/report.md
    output/constrained_clockwork_source_hessian/hessian_spectrum.csv

derivation:
  For one orthogonal crossed-triplet direction, take

    W_link = M sum_{a=0}^{n-1} Lambda_a (X_a - q X_{a+1}).

  The constraint matrix C has rank n and a one-dimensional null vector

    v_a proportional q^{-a}.

  For q=3,n=6,

    ||C v|| = 2.220463e-16,
    |v_6/v_0| = 1.371742e-3,
    kappa = 729.

  The raw chiral Hessian

    H_raw = [[0, C^T], [C, 0]]

  has one intended zero mode per chain.  Its nonzero singular masses are

    2.143405 <= |m|/M_lock <= 3.925024.

  Adding a boundary/source driver

    W_bd = M Xi (v . X)

  stacks the row v^T under C.  The closed Hessian becomes full rank with

    zero modes = 0,
    min |m|/M_lock = 1.000000,
    max |m|/M_lock = 3.925024.

numerical result:
  The induced effective inverse block has singular values

    (1, 0.001371742, 0.001371742, 0.001371742),

  hence condition number 729.  Reusing the dressed C5 scaling gives the
  max-width 1e35 margin at S_T=1e-5:

    margin = 6.216919.

interpretation:
  As a constrained/composite source sector, the nonuniversal threshold vector
  is

    (0,0,0).

  As a complete degenerate link sector, the differential one-loop threshold is
  also zero.  However, as a literal propagating visible tower, the earlier UV
  audit still applies: Lambda_LP/M_G ~= 12.641502, below R=200.  Therefore the
  clockwork lock is now a consistent conditional source mechanism, not a
  perturbative elementary Spin(10) tower.

current obstacle:
  The remaining hard gap is microscopic origin: derive or explicitly model the
  constrained/composite source sector that produces W_link and W_bd while
  staying compatible with the crossed-120 unitary-link D-term quotient and the
  existing Yukawa/Majorana charge table.

next attempted nontrivial idea:
  Merge the constrained clockwork chain with the existing unitary-link D-term
  quotient.  The clockwork should generate the small endpoint block W_kappa,
  while the unitary dilation keeps the physical crossed-triplet mass matrix
  singular values locked at M_lock, so threshold safety and proton suppression
  are both maintained.

verification plan:
  1. Construct an explicit combined constraint system:
       C_clock X = 0,
       P B P = W_kappa[X],
       B^\dagger B = I.
  2. Linearize F-terms plus D-terms around q=3,n=6 and verify that the only
     remaining moduli are the known 9 real unitary-completion moduli.
  3. Confirm that the full triplet mass singular values remain locked to
     M_lock after the clockwork-generated W_kappa is inserted.
  4. If the combined quotient adds propagating charged fields, compute their
     threshold vector and feed it into the two-loop RGE/proton scan.

## 2026-05-09 21:00 Asia/Taipei heartbeat -- clockwork plus unitary-link quotient

status:
  The task remains incomplete, but the clockwork rank-one lock now passes the
  combined unitary-link D-term quotient check.  This is a stronger conditional
  Route A result: proton suppression and physical mass locking can coexist in
  one finite-dimensional constrained-source construction.

new artifact:
  Added

    code/audit_clockwork_unitary_link_quotient.py

  and generated

    output/clockwork_unitary_link_quotient/summary.json
    output/clockwork_unitary_link_quotient/report.md
    output/clockwork_unitary_link_quotient/unitary_completion_scan.csv

derivation:
  Use the constrained clockwork result

    epsilon = 3^-6 = 1.371742e-3,
    kappa = 729,

  to build the visible inverse block

    W_kappa = diag(1, epsilon, epsilon, epsilon).

  Since ||W_kappa||=1, it sits on the contraction boundary but still admits a
  Julia/Halmos unitary dilation:

    B = [[W, sqrt(1-W W^\dagger)],
         [sqrt(1-W^\dagger W), -W^\dagger]].

  This realizes the combined constraints

    P B P = W_kappa,
    B^\dagger B = I,
    M_T = M_lock B^\dagger.

numerical result:
  For the q=3,n=6 card:

    dim B = 8,
    ||B^\dagger B - I||_2 = 2.221326e-16,
    ||PBP - W||_F = 0,
    max |log singular(M_T/M_lock)| = 2.220446e-16.

  The D-term linearization has

    rank = 64,
    D-flat tangent dimension = 64.

  The fixed-block F-term on this D-flat tangent has

    rank = 31,
    residual unitary moduli = 33.

  This is larger than the older 2x2 crossed-link quotient's 9 real moduli,
  because the publication-card block is now 4x4.  The enlarged moduli are
  harmless in the tested finite-dimensional quotient: all sampled completions
  keep

    min singular(M_T/M_lock) = max singular(M_T/M_lock) = 1

  up to 1e-15 numerical precision.

threshold/proton interpretation:
  The clockwork source still gives the dressed max-width 1e35 margin

    6.216919 at S_T=1e-5.

  The physical triplet singular spread is

    2.220446e-16,

  so the complete-degenerate nonuniversal threshold vector is

    (0,0,0).

  As before, a literal propagating visible 10+10bar clockwork tower is
  UV-disfavored; this branch must remain constrained/composite or
  cutoff-local complete-degenerate.

current obstacle:
  The finite-dimensional constrained-source algebra now works.  The remaining
  first-principles gap is microscopic origin: derive the boundary/composite
  sector that enforces W_link, W_boundary, and B^\dagger B=I without adding a
  large elementary Spin(10) tower.

next attempted nontrivial idea:
  Build a boundary/strong-sector origin for the combined quotient, for example
  a hidden vectorlike endpoint completion whose meson matrix is B and whose
  quantum-deformed constraint imposes unitarity or fixed determinant, while
  visible Spin(10) sees only the constrained source block.

verification plan:
  1. Write a microscopic hidden-sector ansatz with Q, Qtilde endpoint fields
     whose meson B=Q Qtilde/Lambda^2 carries the unitary-link data.
  2. Check anomaly cancellation and one-loop beta functions for the hidden
     endpoint gauge factors.
  3. Derive whether the meson constraint can impose B^\dagger B=I or only
     det B = constant; if only determinant is fixed, quantify the residual
     nonunitary leakage.
  4. Feed any non-singlet or non-complete endpoint spectrum back into
     threshold/RGE/proton scans.

## 2026-05-09 21:08 Asia/Taipei heartbeat -- TeX sync check after quotient audit

status:
  The task is still incomplete as a first-principles GUT derivation.  The
  conditional finite-dimensional branch is now locally synchronized: the
  clockwork-plus-unitary-link quotient audit is in the roadmap, the TeX paper
  contains the corresponding proposition, and a final pdflatex pass produced
  paper/gut_framework.pdf with 114 pages.

local verification:
  A log scan after the final compile found no undefined references, no
  Overfull/Underfull boxes, no LaTeX errors, and no remaining rerun warning
  beyond the package-identification line for rerunfilecheck.

current obstacle:
  The remaining hard gap is unchanged and microscopic: derive a boundary,
  composite, or strong-sector origin for W_link, W_boundary, and the unitary
  constraint B^\dagger B=I without introducing a large elementary Spin(10)
  tower or a non-universal threshold spectrum.

next attempted nontrivial idea:
  Start from the hidden endpoint/meson route.  Treat B as a constrained meson
  matrix of vectorlike hidden fields, then test whether the strong constraint
  can enforce the full unitary-link condition or only a determinant/radial
  condition.  If the latter happens, the residual angular leakage becomes the
  next scan parameter.

verification plan:
  1. Build the hidden endpoint ansatz and list all global/gauge charges.
  2. Check hidden gauge anomalies and one-loop beta functions.
  3. Derive the meson constraint algebraically and compare it to
     B^\dagger B=I and PBP=W_kappa.
  4. If nonunitary residual modes remain, quantify their effect on the triplet
     singular spread, threshold vector, and dressed d=5 proton margins.

## 2026-05-09 21:18 Asia/Taipei heartbeat -- 8x8 hidden endpoint meson audit

status:
  The task remains incomplete as an unconditional first-principles GUT
  derivation, but the microscopic-origin question is now sharper.  I upgraded
  the older hidden endpoint/meson tests from the crossed-link 2x2 fixed block
  to the current clockwork publication card with an 8x8 unitary link and a
  4x4 visible block.

files:
  Added code/audit_clockwork_hidden_endpoint_meson.py.

outputs:
  output/clockwork_hidden_endpoint_meson/summary.json
  output/clockwork_hidden_endpoint_meson/report.md
  output/clockwork_hidden_endpoint_meson/finite_endpoint_samples.csv
  output/clockwork_hidden_endpoint_meson/hidden_beta_scan.csv

mathematical result:
  On the hidden endpoint branch

    Q = f I_8,
    Qtilde = f B,
    B = Qtilde Q/f^2,
    P B P = W_kappa,
    B^\dagger B = I,

  the radial/hidden D-term rank is 128.  The hidden U(8)_H gauge orbit has
  rank 64, leaving a 64-real-dimensional unitary meson link.  The visible
  4x4 fixed block has rank 31 on this link, so the residual unitary-completion
  moduli count is

    64 - 31 = 33,

  exactly matching the direct clockwork unitary quotient.

numerical result:
  All finite endpoint samples pass:

    max ||D||_2 = 3.081163e-15,
    max ||B^\dagger B-I||_2 = 3.017862e-15,
    max quantum residual = 1.227317e-16,
    max |log singular(B)| = 9.992007e-16,
    visible threshold vector = (0,0,0).

negative result:
  The N_f=N_c=8 quantum-deformed meson constraint

    det m - beta beta_tilde = (Lambda_H/f)^16

  is compatible with the unitary branch, but it does not imply the branch.
  The holomorphic fixed-block plus determinant system has

    66 complex variables,
    17 complex rank,
    49 complex nullity,

  i.e. 98 real holomorphic flat directions.  With the unitary/radial lock and
  baryons included, the real nullity is 35.  Thus the radial/D-term/Kahler
  lock removes 63 real nonunitary directions and remains a necessary
  conditional ingredient.

current obstacle:
  We have not derived the radial/unitary lock from a microscopic holomorphic
  superpotential.  The hidden meson route says that determinant-type strong
  constraints are not enough.  A boundary, D-term, Kähler quotient, or
  conformal-sector mechanism must be made explicit.

next attempted nontrivial idea:
  Promote the radial lock to a concrete hidden-gauge moose or conformal
  endpoint sector with vectorlike anomaly cancellation.  The target is not to
  change the visible threshold vector, but to prove that the D-term/Kahler
  potential naturally fixes Q Q^\dagger and Qtilde^\dagger Qtilde near f^2 I
  while leaving only the unitary meson B.

verification plan:
  1. Construct a charge table for the hidden U(8)_H endpoint sector and any
     spectator flavors needed for N_f=16 conformal-window control.
  2. Check hidden gauge and mixed global anomalies.
  3. Estimate the hidden beta function and the allowed matching window for
     N_f=8 and N_f=16.
  4. Derive the radial mass matrix and show whether nonunitary meson modes
     are lifted above M_lock without visible Spin(10) thresholds.

## 2026-05-09 21:28 Asia/Taipei heartbeat -- 8x8 endpoint vectorlike completion

status:
  The task remains incomplete as an unconditional first-principles GUT
  derivation, but the latest hidden endpoint route now has an 8x8
  anomaly/beta bookkeeping card rather than relying on the older 4x4 crossed
  link result.

files:
  Added code/audit_clockwork_endpoint_vectorlike_completion.py.

outputs:
  output/clockwork_endpoint_vectorlike_completion/summary.json
  output/clockwork_endpoint_vectorlike_completion/report.md
  output/clockwork_endpoint_vectorlike_completion/field_content.csv
  output/clockwork_endpoint_vectorlike_completion/anomaly_audit.csv
  output/clockwork_endpoint_vectorlike_completion/beta_audit.csv
  output/clockwork_endpoint_vectorlike_completion/visible_threshold_audit.csv

mathematical result:
  For endpoint groups

    SU(8)_L x SU(8)_H x SU(8)_R

  the raw fields

    Q       ~ (8, 8bar, 1),
    Qtilde  ~ (1, 8, 8bar)

  have cubic anomaly proxies

    A_L = 8,
    A_H = 0,
    A_R = -8.

  Hence raw endpoint gauging is chiral and not acceptable.  Adding

    Qc       ~ (8bar, 8, 1),
    Qtilde_c ~ (1, 8bar, 8)

  cancels all hidden cubic anomaly proxies.

numerical/beta result:
  For the vectorlike-completed field content, using

    b = sum_R T(R) - 3N,   N=8,

  the one-loop hidden beta coefficients are

    b_L = -16,
    b_H = -8,
    b_R = -16.

  Therefore there is no one-loop UV Landau pole before R=200 for the tested
  matching couplings.  The bookkeeping identifies SU(8)_H as equivalent
  N_f=16, inside the SQCD conformal-window range, while SU(8)_L and SU(8)_R
  sit at the N_f=N_c edge.

threshold result:
  All endpoint-completion fields are visible Spin(10) singlets, so

    Delta b_visible = (0,0,0).

current obstacle:
  Endpoint anomaly and one-loop beta bookkeeping no longer blocks the 8x8
  hidden meson route.  The remaining hard issue is still dynamical: derive the
  radial/Kahler potential or D-term sector that fixes

    Q Q^\dagger = f^2 I,
    Qtilde^\dagger Qtilde = f^2 I

  rather than imposing it as a finite-rank moment-map constraint.

next attempted nontrivial idea:
  Build a radial-driver/mass Hessian audit for the vectorlike endpoint
  completion: include Qc and Qtilde_c, write the most economical hidden
  singlet/radial driver terms, and verify that nonunitary radial meson modes
  are lifted while the unitary B moduli and visible threshold vector remain
  harmless.

verification plan:
  1. Propose W_radial or a D-term potential for Q,Qtilde,Qc,Qtilde_c.
  2. Linearize around Q=fI, Qtilde=fB and conjugate partners at the vectorlike
     locked point.
  3. Count massive nonunitary modes, hidden gauge orbits, and residual unitary
     B moduli.
  4. Confirm no visible Spin(10)-charged thresholds are introduced.

## 2026-05-09 22:25 Asia/Taipei heartbeat -- radial-driver Hessian for 8x8 endpoint lock

status:
  The task remains incomplete as an unconditional first-principles GUT
  derivation.  However, the local hidden endpoint route is now closed at the
  finite-rank D-term/Kahler Hessian level: the radial/unitary lock is no
  longer just a rank-count slogan.

files:
  Added code/audit_clockwork_radial_driver_hessian.py.

outputs:
  output/clockwork_radial_driver_hessian/summary.json
  output/clockwork_radial_driver_hessian/report.md
  output/clockwork_radial_driver_hessian/finite_radial_driver_samples.csv

tested potential:
  Around

    Q = I_8,
    T = B,
    Qc = Tc = 0,
    B^\dagger B = I,
    PBP = W_kappa,

  the tested D-term/Kahler driver is

    V = m_D^2 ||mu_rad/H(Q,T)||^2
        + m_F^2 ||P(TQ)P - W_kappa||^2
        + m_c^2 (||Qc||^2 + ||Tc||^2),

  with

    mu_rad/H = (Q Q^\dagger-I, T^\dagger T-I, Q^\dagger Q-TT^\dagger).

rank result:
  Active Q,T real variables = 256.

    D-term rank on active fields = 128,
    D-term plus fixed-block rank = 159,
    full rank with vectorlike partners = 415,
    full raw nullity = 97,
    hidden gauge orbit rank = 64,
    physical nullity after quotient = 33.

  The decomposition is therefore:

    97 raw flat directions = 64 hidden gauge directions + 33 physical unitary
    completion moduli.

  The 33 physical moduli exactly match the previous unitary quotient audit.

numerical result:
  finite samples give

    max ||B^\dagger B-I||_2 = 3.403932e-15,
    max |log singular(B)| = 1.554312e-15,
    vectorlike partner residual = 0,
    visible threshold vector = (0,0,0).

current obstacle:
  The hidden endpoint sector now passes as a local D-term/Kahler finite-rank
  completion.  What is still not first-principles is the microscopic origin of
  this D-term/Kahler potential itself: we have not derived why the strong or
  boundary dynamics must generate precisely these radial moment-map terms.

next attempted nontrivial idea:
  Either:
    A. formulate the endpoint lock as a hidden supersymmetric gauge quotient
       with FI/radial data and show the above moment maps are its actual
       D-flat equations; or
    B. demote this to an explicit conditional D-term/Kahler postulate and move
       effort back to the remaining phenomenology blockers: full flavor fit
       and fully dressed d=5 proton decay.

verification plan:
  1. If choosing A, write the hidden gauge/Kahler action and derive the moment
     maps rather than imposing them.
  2. Check whether FI/radial parameters are technically natural under the
     hidden beta functions already audited.
  3. If choosing B, mark the endpoint lock as a conditional EFT assumption in
     the final theorem map and prioritize full CKM/PMNS plus d=5 proton decay.

## 2026-05-09 23:31 Asia/Taipei heartbeat -- hidden gauge quotient origin of radial lock

status:
  The task remains incomplete as an unconditional first-principles GUT
  derivation.  But the endpoint radial/unitary equations have now been derived
  as hidden gauge-quotient moment maps, rather than only written as an
  ad hoc real potential.

files:
  Added code/audit_clockwork_hidden_gauge_quotient_origin.py.

outputs:
  output/clockwork_hidden_gauge_quotient_origin/summary.json
  output/clockwork_hidden_gauge_quotient_origin/report.md
  output/clockwork_hidden_gauge_quotient_origin/u1_charge_table.csv
  output/clockwork_hidden_gauge_quotient_origin/u1_trace_audit.csv

moment-map derivation:
  Gauge

    U(8)_L x U(8)_H x U(8)_R

  with endpoint charges

    Q  : ( 1,-1, 0),
    T  : ( 0, 1,-1),
    Qc : (-1, 1, 0),
    Tc : ( 0,-1, 1)

  for the abelian factors.  The D-term maps are

    mu_L = Q Qdag - Qcdag Qc - I,
    mu_H = Qdag Q - T Tdag - Qc Qcdag + Tcdag Tc,
    mu_R = Tdag T - Tc Tcdag - I.

  On the branch

    Q=I_8, T=B, Qc=Tc=0, Bdag B=I,

  these become exactly the radial lock equations used in the Hessian:

    Q Qdag=I,
    Tdag T=I,
    Qdag Q=T Tdag.

rank result:
  Linearized moment-map rank = 128.
  Moment plus fixed-block rank = 159.
  Fixed-block increment rank = 31.

  These exactly match the previous radial-driver Hessian on the active
  Q,T fields.

FI/anomaly result:
  The vectorlike U(1) charge traces vanish:

    Tr Q_U1L = Tr Q_U1H = Tr Q_U1R = 0,

  and the cubic U(1) proxies also vanish.  Thus the FI/radial levels are
  one-loop trace-safe in this bookkeeping.  The nonabelian vectorlike anomaly
  and beta card from the previous heartbeat remains valid.

threshold result:
  All fields are visible Spin(10) singlets, so the direct visible threshold
  vector remains

    (0,0,0).

current obstacle:
  The local endpoint lock is now a hidden gauge quotient with trace-safe FI
  data.  What remains open is not the local D-term algebra; it is whether this
  hidden Kahler/FI sector has a deeper UV origin or should be declared as a
  conditional EFT assumption.

next attempted nontrivial idea:
  Stop pushing this sector unless the goal is a fully UV-complete hidden
  theory.  The rational next move is to update the theorem map: mark the
  endpoint D-term/Kahler quotient as a conditional assumption, then return to
  the two remaining publication blockers:

    1. full CKM/PMNS/flavor fit,
    2. fully dressed d=5 proton decay tables.

verification plan:
  1. Add the hidden gauge quotient result to the conditional theorem ledger.
  2. Make sure the paper states this is conditional on the hidden Kahler/FI
     sector, not a holomorphic theorem.
  3. Resume either the full flavor fit or the dressed d=5 proton decay scan.

## 2026-05-10 00:29 Asia/Taipei heartbeat -- clockwork-rescued dressed C5 replay

status:
  The task remains incomplete as an unconditional first-principles GUT
  derivation.  But the publication-card dimension-five obstruction has moved:
  the old finite kappa=30 eigenstate-card replay is now a genuine NO-GO row,
  while the q=3,n=6 clockwork-rescued branch passes the exact local dressed
  C5 stress replay.

files:
  Added code/audit_publication_dressed_c5_clockwork_rescue.py.
  Updated code/audit_conditional_theorem_ledger.py.

outputs:
  output/publication_dressed_c5_clockwork_rescue/summary.json
  output/publication_dressed_c5_clockwork_rescue/report.md
  output/publication_dressed_c5_clockwork_rescue/dressed_channel_rows_kappa729.csv

mathematical bridge:
  The previous finite eigenstate-card replay used

    W_30 = diag(1, 1/30, 1/30, 1/30)

  and failed because the dressed Wilson tensors scale linearly in the inverse
  crossed-triplet block.  The finite clockwork branch instead generates

    W_kappa = diag(1, epsilon, epsilon, epsilon),
    epsilon = q^{-n},
    kappa = q^n.

  For the preferred card

    q = 3,
    n = 6,
    kappa = 729,
    epsilon = 1.371742112482853e-3.

  Since C5(kappa) is recomputed directly as

    C5(kappa) = C5_rank_one + (30/kappa)(C5_finite30-C5_rank_one),

  this is not a new scalar leakage shortcut; every channel row still starts
  from the exported publication eigenstate C5L/C5R tensors.

numerical result:
  Old finite kappa=30 replay:

    central unsafe 1e35 rows   = 25,
    central S_T^max            = 2.305290e-6,
    max-width unsafe 1e35 rows = 74,
    max-width S_T^max          = 1.026080e-6.

  Exact kappa=729 replay:

    central rows               = 7056,
    central unsafe 1e35 rows   = 0,
    central S_T^max            = 5.601856e-5,
    central worst margin       = 31.3807879,

    max-width rows             = 7056,
    max-width unsafe 1e35 rows = 0,
    max-width S_T^max          = 2.493375e-5,
    max-width worst margin     = 6.21691896,

    min-width rows             = 7056,
    min-width unsafe 1e35 rows = 0,
    min-width S_T^max          = 1.304227e-4,
    min-width worst margin     = 170.100788.

  The same worst channel remains p -> mu+ pi0 under the
  democratic_all_stress / near_unsafe_20TeV row, but it is now safely above
  the 1e35 yr stress target at S_T=1e-5.

ledger update:
  The theorem ledger now records:

    finite kappa=30 publication-card eigenstate C5 replay is sufficient:
      NO_GO,

    clockwork-rescued publication-card dressed C5 replay passes the full
    stress grid:
      PASS_CONDITIONAL.

  Updated status counts:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 53,
    TUNED_FALLBACK: 3.

current obstacle:
  The local d=5 proton proxy is now conditionally safe, but it is not yet a
  paper-grade proton-decay theorem.  The remaining d=5 obstacle is packaging:
  regenerate the final channel-specific tables from the kappa=729 card with an
  explicit hadronic/chiral input manifest, convention ledger, and exact
  reproduction script.  In parallel, the full CKM/PMNS/flavor fit remains open.

next attempted nontrivial idea:
  Promote the kappa=729 rescue from a local proxy pass to a publication card:
  make a single machine-readable "clockwork-rescued flavor+d5 card" that
  contains Yukawa matrices, seesaw matrices, CKM/PMNS observables, C5L/C5R
  tensors, dressed channel tables, chiral constants, and all hashes.  This
  should prevent the old problem where flavor, triplet tensors, and proton
  tables lived in separate caches.

verification plan:
  1. Generate the final d=5 channel table package from
     output/publication_dressed_c5_clockwork_rescue/summary.json and the exact
     dressed channel CSV.
  2. Record central/min/max hadronic-width cases, present bounds, the 1e35 yr
     stress target, and the worst row for every channel/operator pair.
  3. Merge the card with the existing publication flavor reproducibility
     output and check whether the full CKM/PMNS/flavor observable table can be
     regenerated from the same manifest.
  4. Only then consider changing the remaining d=5 OPEN item to a
     PASS_CONDITIONAL theorem-row.

## 2026-05-10 01:36 Asia/Taipei heartbeat -- single clockwork-rescued publication manifest

status:
  The task remains incomplete as an unconditional first-principles GUT
  derivation.  But the cache-synchronization problem between flavor, seesaw,
  hidden clockwork, and dressed proton decay is now locally closed: all of
  those data have been merged into one hash-locked manifest.

files:
  Added code/build_clockwork_rescued_publication_card.py.
  Updated code/audit_conditional_theorem_ledger.py.

outputs:
  output/clockwork_rescued_publication_card/clockwork_rescued_publication_card.json
  output/clockwork_rescued_publication_card/report.md
  output/clockwork_rescued_publication_card/channel_worst_rows.csv
  output/clockwork_rescued_publication_card/flavor_observables.csv
  output/clockwork_rescued_publication_card/input_manifest.csv

manifest content:
  The new manifest binds together:

    publication_closure_card/publication_closure_card.json,
    publication_flavor_d5_reproducibility/summary.json,
    publication_triplet_eigenstate_card/triplet_eigenstate_card.json,
    publication_dressed_c5_clockwork_rescue/summary.json,
    publication_dressed_c5_clockwork_rescue/dressed_channel_rows_kappa729.csv,
    clockwork_hidden_gauge_quotient_origin/summary.json,
    conditional_theorem_ledger/summary.json.

  Each input is recorded with size and sha256 hash.

flavor result:
  From the single manifest:

    CKM score        = 7.634719e-4,
    mass score       = 1.490597e-1,
    seesaw residual  = 5.173996e-12,

    M_R = (2.380448e10, 3.220380e13, 3.926531e15) GeV.

  The flavor card recomputes exactly from its embedded Yukawa matrices and
  passes the local strict CKM, loose mass, and seesaw gates.

d5 result:
  The exact kappa=729 dressed channel CSV has 21168 rows.  Compressing by
  normalization/operator/channel gives 21 worst-row records.  Every worst row
  passes the 1e35 yr stress target at S_T=1e-5.

  The global worst row is:

    normalization = max_width,
    operator      = LLLL,
    channel       = p -> mu+ pi0,
    scenario      = democratic_all_stress,
    epsilon       = 0.9,
    spectrum      = near_unsafe_20TeV,
    margin        = 6.216918955.

  The weakest K+ nubar row is also safe:

    max_width LLLL p -> K+ nubar margin = 14.8805947.

ledger update:
  Added a PASS_CONDITIONAL row:

    Single clockwork-rescued flavor+d5 manifest is locally complete.

  Updated counts:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 54,
    TUNED_FALLBACK: 3.

current obstacle:
  Locally, flavor + seesaw + dressed d=5 + hidden clockwork are now in one
  reproducible card.  The remaining obstacle is publication-grade external
  closure: replace internal target values and chiral/hadronic constants by a
  final cited input table, then decide whether the d=5 OPEN row can be
  honestly relabelled PASS_CONDITIONAL.  The full CKM/PMNS wording also needs
  tightening: we currently have a local CKM/seesaw fit, not a full modern
  global-data flavor fit.

next attempted nontrivial idea:
  Build a final "input convention ledger" that contains every numerical
  target and every width prefactor used by the publication card, with a local
  provenance hash.  Since this thread forbids web search, this should be
  presented as a no-web internal convention card, not as a final literature
  refresh.

verification plan:
  1. Extract all CKM, mass-ratio, seesaw, d=6, d=5 width, and proton-bound
     constants into a single JSON/CSV input ledger.
  2. Re-run the clockwork publication card using only that input ledger.
  3. Check that the regenerated flavor_observables.csv and
     channel_worst_rows.csv are bitwise or numerically identical.
  4. Then the remaining d=5 OPEN item can be narrowed to "external citation
     refresh only"; if the user later permits web lookup, replace the local
     input ledger with cited current values.

## 2026-05-10 heartbeat update: no-web input convention ledger

status:
  Incomplete as an unconditional GUT derivation.  The branch remains a
  conditional Spin(10)/Pati-Salam EFT.  The local flavor + seesaw +
  clockwork-rescued dressed d=5 card is now internally convention-locked.

completed in this heartbeat:
  Added `code/build_no_web_input_convention_ledger.py` and generated

    output/no_web_input_convention_ledger/input_convention_ledger.json,
    output/no_web_input_convention_ledger/summary.json,
    output/no_web_input_convention_ledger/constant_rows.csv,
    output/no_web_input_convention_ledger/report.md.

  The ledger records 78 local numerical convention rows: flavor targets,
  flavor gates, seesaw benchmark inputs, d=5 width prefactors, channel bounds,
  display filter, clockwork kappa, lifetime constants, soft-spectrum points,
  MSSM mixing inputs, RGE inputs, and legacy local proton constants.  Its
  constant-row digest is

    f3a354855dc9983fe62c3b9b7434b516e0ddd9f88828844bf2af5381c6f4db04.

  Numerical replay:

    dressed d=5 rows checked        = 21168,
    width-prefactor mismatches      = 0,
    present-bound mismatches        = 0,
    display-filter mismatches       = 0,
    kappa mismatches                = 0,
    max C6 relative error           = 0,
    max tau relative error          = 0,
    max present-margin relative err = 0,
    max 1e35-margin relative err    = 0.

ledger update:
  Added a PASS_CONDITIONAL row:

    No-web input convention ledger reproduces the dressed d=5 card arithmetic.

  Updated counts:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 55,
    TUNED_FALLBACK: 3.

current obstacle:
  The publication d=5 bottleneck is no longer an internal arithmetic or cache
  mismatch.  It is now a convention-to-literature replacement problem: the
  no-web local targets and width constants must eventually be replaced by a
  cited final input table, then the same ledger and dressed-row replay must be
  regenerated.  Full CKM/PMNS remains a local CKM/seesaw fit rather than a
  full global-data flavor fit.

next attempted nontrivial idea:
  Split the remaining "full flavor + full d=5 + supplemental reproducibility"
  blocker into two tracks:

    A. a no-web internal closure track, where every table in the TeX is
       sourced from the new ledger and manifest hashes;
    B. a later citation-refresh track, only if web/literature lookup is
       allowed, replacing the no-web convention values while preserving the
       same code-level replay.

verification plan:
  1. Insert the no-web input convention ledger into the TeX as a proposition
     and appendix pointer.
  2. Compile and verify the PDF/log.
  3. Next heartbeat: make the TeX tables explicitly point to the generated
     JSON/CSV hashes, then decide whether the d=5 OPEN item should be renamed
     "external input refresh remains open" instead of "channel-specific d=5
     proton decay is completed".

## 2026-05-10 heartbeat update: d=5 open item narrowed

status:
  Incomplete as an unconditional first-principles GUT.  The no-web local
  clockwork-rescued flavor+d=5 package is internally arithmetic-consistent,
  but the final theory still has three open blockers.

completed in this heartbeat:
  The theorem ledger OPEN item was renamed from

    Channel-specific d=5 proton decay is completed

  to the narrower and more accurate

    External cited d=5 input refresh is completed.

  This records the fact that the exact clockwork-rescued dressed replay already
  passes locally:

    local dressed rows   = 21168,
    worst 1e35 margin    = 6.216918955,
    no-web ledger replay = true.

  The TeX no-web input convention proposition now includes short SHA prefixes
  for the key machine artifacts:

    constant_rows.csv                         f3a354855dc9,
    clockwork_rescued_publication_card.json   92adc5b90ba1,
    dressed_channel_rows_kappa729.csv         e9ab804c9577,
    proton_decay_verification.json            8a52caa56c40.

ledger update:
  Counts are unchanged, but the meaning of the d=5 OPEN row is sharper:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 55,
    TUNED_FALLBACK: 3.

  Current OPEN blockers are:

    1. Full CKM/PMNS/flavor fit is completed.
    2. External cited d=5 input refresh is completed.
    3. A5-A6 have a microscopic first-principles origin.

current obstacle:
  The d=5 proton-decay calculation is locally closed only as a no-web
  convention-locked proxy.  It still needs an externally cited hadronic/chiral
  input table and final flavor target table before it can be called
  publication-final.  The larger first-principles issue remains A5-A6: the
  hidden clockwork/Kahler/FI quotient is conditional, not derived from PSLT
  alone.

next attempted nontrivial idea:
  Build a "no-web table provenance appendix" that makes every numerical table
  in the TeX traceable to one generated JSON/CSV artifact, then audit whether
  the full-flavor blocker can be similarly narrowed into (i) local matrix
  reproducibility and (ii) external target refresh.  This keeps the mathematical
  core honest while reducing publication risk.

verification plan:
  1. Recompile the TeX after the artifact-hash table insertion.
  2. Check the LaTeX log for warnings/errors/overfull boxes.
  3. Next: add a generated table-provenance appendix or CSV-to-TeX manifest,
     then verify that no major numerical claim in the paper is orphaned from a
     machine-readable output.

## 2026-05-10 heartbeat update: table provenance manifest

status:
  Still incomplete as an unconditional first-principles GUT.  The local
  clockwork-rescued flavor+d=5 package is now convention-locked and
  table-provenance locked.  The remaining blockers are unchanged in count:
  full CKM/PMNS flavor fit, external cited d=5 input refresh, and A5-A6
  microscopic origin.

completed in this heartbeat:
  Removed a provenance cycle from `build_clockwork_rescued_publication_card.py`:
  the publication card no longer hashes the conditional theorem ledger as an
  input, because the theorem ledger itself reads the publication card.  The
  publication card now manages only its local flavor+d5+hidden-clockwork
  package, while the theorem ledger owns the OPEN blocker list.

  Regenerated:

    output/clockwork_rescued_publication_card/clockwork_rescued_publication_card.json,
    output/no_web_input_convention_ledger/summary.json,
    output/conditional_theorem_ledger/summary.json.

  The new publication-card hash prefix is

    563d883ce1f6.

  Added `code/build_no_web_table_provenance_manifest.py` and generated

    output/no_web_table_provenance_manifest/provenance_rows.csv,
    output/no_web_table_provenance_manifest/summary.json,
    output/no_web_table_provenance_manifest/report.md.

numerical verification:
  The provenance manifest contains 14 rows.  All artifacts exist.  For the
  11 rows whose numerical key is required to appear directly in the TeX, all
  11 are found:

    required TeX keys found = 11/11,
    all required TeX keys found = true.

  Covered printed keys include:

    convention digest f3a354855dc9...,
    publication card hash 563d883ce1f6,
    d5 row CSV hash e9ab804c9577,
    proton convention hash 8a52caa56c40,
    CKM score 7.634719e-4,
    mass score 1.490597e-1,
    seesaw residual 5.173996e-12,
    d5 rows 21168,
    worst d5 margin 6.216918955,
    weakest K+ nubar margin 14.8805947,
    clockwork kappa 729.

TeX sync:
  Added a "No-web table provenance manifest" proposition to the paper.  The
  TeX no-web input ledger hash table now uses the current publication-card
  prefix 563d883ce1f6.

current obstacle:
  No major local closure-package numerical claim is now orphaned from generated
  artifacts.  The remaining risk is physical rather than bookkeeping:
  replace local no-web inputs by cited external values, then rerun the same
  provenance and d=5 ledger; and separately derive or demote the A5-A6 hidden
  clockwork/Kahler/FI quotient.

next attempted nontrivial idea:
  Start attacking the full-flavor blocker in the same style: build a
  flavor-target provenance ledger that separates local matrix reproducibility
  from external target replacement, then check whether PMNS/CKM target updates
  can be swapped without invalidating the source-consistent crossed-120
  branch.

verification plan:
  1. Compile the updated TeX and check the log.
  2. On the next heartbeat, generate a flavor-target provenance ledger:
     local targets, fitted observables, residuals, matrices, and external
     placeholders separated into distinct machine-readable rows.
  3. If the no-web target ledger is clean, update the OPEN blocker wording
     from "Full CKM/PMNS/flavor fit is completed" to the sharper
     "External cited flavor target refresh and PMNS completion are completed".

## 2026-05-10 12:08 UTC heartbeat: flavor-target provenance split

status:
  Still incomplete as an unconditional first-principles GUT.  The current
  result remains a verified conditional Spin(10)/Pati-Salam EFT branch.  This
  heartbeat narrows the flavor blocker: the local source-consistent flavor card
  is reproducible, while external target refresh and PMNS completion remain
  open.

completed in this heartbeat:
  Added `code/build_no_web_flavor_target_provenance.py` and generated

    output/no_web_flavor_target_provenance/target_rows.csv,
    output/no_web_flavor_target_provenance/seesaw_rows.csv,
    output/no_web_flavor_target_provenance/matrix_manifest.csv,
    output/no_web_flavor_target_provenance/summary.json,
    output/no_web_flavor_target_provenance/report.md.

  Updated `code/audit_conditional_theorem_ledger.py` so the previous OPEN
  blocker

    Full CKM/PMNS/flavor fit is completed

  is now the sharper

    External cited flavor target refresh and PMNS completion are completed.

  Updated `code/build_no_web_table_provenance_manifest.py` to include the new
  flavor-target ledger and its all-Yukawa matrix digest.  Added a new TeX
  proposition, "No-web flavor-target provenance ledger", and updated the
  no-web table provenance manifest proposition.

numerical verification:
  The flavor-target provenance ledger contains

    target rows = 12,
    seesaw rows = 15,
    matrix-manifest rows = 5.

  The replayed local scores exactly reproduce the stored publication flavor
  card:

    CKM score diff = 0.000e+00,
    mass score diff = 0.000e+00.

  The local source-consistent gates remain true:

    local targets match no-web ledger = true,
    card recomputes within tolerance = true,
    local flavor gates pass = true.

  The ledger intentionally keeps the two publication-level items false:

    external_target_refresh_done = false,
    pmns_completion_done = false.

  The combined digest of the embedded Y_u,Y_d,Y_e,Y_nu matrices begins

    dccb2b3d8002.

  Regenerated the conditional theorem ledger:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 55,
    TUNED_FALLBACK: 3.

  The current OPEN blockers are now:

    1. External cited flavor target refresh and PMNS completion are completed.
    2. External cited d=5 input refresh is completed.
    3. A5-A6 have a microscopic first-principles origin.

  Regenerated the no-web table provenance manifest.  It now has 16 rows, with
  all 12 required TeX keys found:

    required TeX keys found = 12/12,
    provenance_manifest_complete = true.

current obstacle:
  There is no longer a local hidden cache/provenance mismatch in the flavor
  card.  The remaining flavor obstacle is physical and publication-facing:
  choose final cited CKM, quark/lepton mass, neutrino oscillation, phase, and
  uncertainty inputs; then rerun the same source-consistent flavor/d5 card and
  complete PMNS observables in the common field basis.

next attempted nontrivial idea:
  Build an external-input placeholder manifest schema with columns for value,
  uncertainty, scale, scheme, running convention, and citation key.  Then run a
  no-web dry replay using the current local values through that schema; when
  references are allowed, only the manifest rows need replacement.

verification plan:
  1. Compile the TeX twice and check the log.
  2. Add a PMNS observable replay table sourced from the same embedded
     Y_e,Y_nu,m_nu matrices, including angles, phase convention, and
     mass-splitting residuals.
  3. Replace the dry local target rows by cited external target rows and rerun
     `build_no_web_flavor_target_provenance.py`, the conditional theorem
     ledger, and the table provenance manifest.

## 2026-05-10 16:22 UTC heartbeat: PMNS benchmark replay ledger

status:
  Still incomplete as an unconditional first-principles GUT.  The branch is
  still a conditional Spin(10)/Pati-Salam EFT construction.  This heartbeat
  narrows the flavor blocker further by separating an exact local PMNS
  benchmark replay from the still-open source-consistent/global PMNS fit.

completed in this heartbeat:
  Added `code/build_no_web_pmns_benchmark_replay.py` and generated

    output/no_web_pmns_benchmark_replay/observable_rows.csv,
    output/no_web_pmns_benchmark_replay/matrix_rows.csv,
    output/no_web_pmns_benchmark_replay/summary.json,
    output/no_web_pmns_benchmark_replay/report.md.

  Updated `code/audit_conditional_theorem_ledger.py` so the flavor OPEN item
  now records both facts:

    pmns_bench = true,
    pmns_pub = false.

  Updated `code/build_no_web_table_provenance_manifest.py` and the TeX paper
  so the PMNS benchmark pair digest is part of the no-web provenance manifest.

numerical verification:
  The PMNS benchmark replay contains

    observable rows = 13,
    matrix rows = 5.

  The PMNS pair digest begins

    f8c599abcf19.

  The local exact CP1/O(2) seesaw benchmark reconstructs

    sin^2 theta12 = 0.30400000000005034,
    sin^2 theta13 = 0.022200000000032458,
    sin^2 theta23 = 0.5729999999999305,

  against the local target

    (0.304, 0.0222, 0.573).

  Residuals:

    max PMNS angle residual = 6.950e-14,
    max light-mass residual = 1.025e-13 eV,
    max mass-splitting residual = 1.029e-14 eV^2,
    seesaw matrix residual = 2.779e-12.

  The local PMNS benchmark replay passes the 1e-10 gate.

  Regenerated the conditional theorem ledger.  Counts remain:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 55,
    TUNED_FALLBACK: 3.

  Regenerated the no-web table provenance manifest.  It now has 18 rows, with
  all 13 required TeX keys found:

    required TeX keys found = 13/13,
    provenance_manifest_complete = true.

TeX sync:
  Added a proposition "No-web PMNS benchmark replay".  Updated the no-web
  table provenance proposition from 16/12 rows to 18/13 rows and included the
  PMNS digest.

current obstacle:
  The local PMNS benchmark itself is reproducible, but it is not yet the final
  publication flavor theorem.  The remaining task is to attach a PMNS replay to
  the source-consistent crossed-120 publication card after final cited
  CKM/mass/neutrino target replacement.

next attempted nontrivial idea:
  Build a source-consistent PMNS replay card that uses the current
  publication-closure Y_e and Y_nu matrices, then explicitly states whether
  the PMNS target is inherited from the old exact CP1/O(2) benchmark or must be
  refit.  This will tell us whether the source-consistent CKM repair preserved
  PMNS or only preserved the quark/mass/d5 side.

verification plan:
  1. Extract U_e from the publication-closure charged-lepton matrix.
  2. Construct the PMNS target/replay convention in the same basis and compare
     angles, Jarlskog, masses, and splittings.
  3. If the source-consistent card fails PMNS, add a constrained lepton-sector
     deformation or demote the flavor closure to CKM/mass/seesaw-only until a
     joint PMNS refit is performed.

## 2026-05-10 17:22 UTC heartbeat: source-consistent PMNS compatibility

status:
  Still incomplete as an unconditional first-principles GUT.  The branch remains
  a conditional Spin(10)/Pati-Salam EFT construction.  This heartbeat checks
  whether the source-consistent crossed-120 publication card actually remains
  PMNS-compatible after the CKM repair.

completed in this heartbeat:
  Added `code/build_source_consistent_pmns_replay.py` and generated

    output/source_consistent_pmns_replay/observable_rows.csv,
    output/source_consistent_pmns_replay/matrix_manifest.csv,
    output/source_consistent_pmns_replay/summary.json,
    output/source_consistent_pmns_replay/report.md.

  Updated `code/audit_conditional_theorem_ledger.py` so the flavor OPEN item
  now records:

    pmns_bench = true,
    source_pmns = true,
    pmns_pred = false.

  Updated `code/build_no_web_table_provenance_manifest.py` and the TeX paper to
  include the source-consistent PMNS compatibility digest.

mathematical derivation:
  For the publication card's embedded Y_e and Y_nu, compute U_e from

    U_e^\dagger Y_e Y_e^\dagger U_e = diagonal.

  Fix the same local PMNS benchmark U_PMNS and define

    U_nu = U_e U_PMNS,
    m_nu = U_nu^* diag(m_1,m_2,m_3) U_nu^\dagger.

  With m_D = 100 GeV * Y_nu, reconstruct

    M_R = -m_D^T m_nu^{-1} m_D.

  The replay then checks

    m_nu^reco = -m_D M_R^{-1} m_D^T

  and compares PMNS angles, light masses, splittings, and the matrix residual.
  This is a compatibility proof, not a predictive Majorana texture derivation.

numerical verification:
  The source-consistent PMNS replay contains

    observable rows = 9,
    matrix rows = 4.

  The source PMNS pair digest begins

    e33901f5bb79.

  Residuals:

    max PMNS angle residual = 1.110e-16,
    max light-mass residual = 2.014e-13 eV,
    max mass-splitting residual = 2.022e-14 eV^2,
    seesaw matrix residual = 5.174e-12,
    theta_norm = 4.359e-11.

  Heavy Majorana singular values:

    M_R = (2.3804475922e10, 3.2203796830e13, 3.9265310084e15) GeV,
    condition number = 1.6494927346e5.

  The source-consistent PMNS compatibility gate passes.  The predictive PMNS
  fit gate remains false because M_R is reconstructed from the chosen PMNS
  target.

  Regenerated the conditional theorem ledger.  Counts remain:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 55,
    TUNED_FALLBACK: 3.

  Regenerated the no-web table provenance manifest.  It now has 20 rows, with
  all 14 required TeX keys found:

    required TeX keys found = 14/14,
    provenance_manifest_complete = true.

current obstacle:
  The CKM-repaired source-consistent card is PMNS-compatible, but only because
  the inverse seesaw reconstructs M_R.  The remaining flavor obstacle is now
  sharply the lack of a predictive or externally refit Majorana/PMNS sector:
  either derive M_R from a constrained Majorana texture, or install cited final
  CKM/mass/neutrino targets and perform a genuine joint source-consistent fit.

next attempted nontrivial idea:
  Build a Majorana-texture rank audit: restrict M_R to the same CP1/O(2)
  Veronese plus trace/contact/shadow basis used elsewhere, then solve for the
  best approximation to the reconstructed source-consistent M_R.  If the
  residual is small and stable under target variation, PMNS can be promoted
  from compatibility to a constrained texture fit.

verification plan:
  1. Decompose the reconstructed M_R into Veronese, trace-lift, and contact
     components.
  2. Report the relative residual, condition number, and sensitivity of PMNS
     angles to coefficient perturbations.
  3. If the residual is too large, mark the Majorana sector as an explicit
     free inverse-seesaw input in the paper rather than a derived PSLT/CP1
     prediction.

## 2026-05-10 18:24 Taipei heartbeat: source Majorana texture rank audit closed

status:
  The task is still incomplete as an unconditional first-principles GUT
  derivation, but the immediate Majorana-texture ambiguity is now sharply
  quantified.  The source-consistent PMNS replay remains compatible only after
  inverse-seesaw reconstruction of M_R; the reconstructed M_R is not a
  Veronese-only CP1/O(2) texture.

files:
  - `code/audit_source_majorana_texture_rank.py`
  - `output/source_majorana_texture_rank/summary.json`
  - `output/source_majorana_texture_rank/report.md`
  - `code/audit_conditional_theorem_ledger.py`
  - `code/build_no_web_table_provenance_manifest.py`
  - `paper/gut_framework.tex`

mathematical derivation:
  Starting from the source-consistent replay, reconstruct

    M_R = -m_D^T m_nu^{-1} m_D.

  Normalize by the largest singular value M_* and decompose the complex
  symmetric matrix into

    M_R / M_* = M_V + zeta K + R_perp,

  where M_V lies in the five-complex-dimensional CP1/O(2) Veronese Majorana
  subspace obeying

    (M_V)_{11} = 2 (M_V)_{02},

  and

    K = (1/sqrt(3)) [[0,0,1],[0,-1,0],[1,0,0]]

  is the orthogonal trace/contact direction.  The projection is just the
  Frobenius-orthogonal split of Sym^2(Sym^2 C^2) into the Veronese component
  plus the spin-zero contact component.

numerical verification:
  - largest singular scale: `3.926531e15 GeV`
  - condition number: `1.649493e5`
  - Veronese-only relative residual: `1.304275e-01`
  - contact fraction: `1.304275e-01`
  - Veronese+contact residual: `2.687332e-17`
  - post-projection Veronese constraint residual: `1.570092e-16`
  - contact-lift essential: `true`
  - predictive Majorana texture closed: `false`

ledger/provenance sync:
  Regenerated the conditional theorem ledger.  Counts are now:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 56,
    TUNED_FALLBACK: 3.

  Regenerated the no-web table provenance manifest after adding the source
  Majorana rank artifact.  It now has:

    rows = 22,
    required TeX keys found = 15/15,
    provenance_manifest_complete = true.

current obstacle:
  The Majorana/PMNS sector is still not predictive.  The contact component is
  O(0.13), so it cannot be dismissed as a small numerical artifact.  It must
  either be derived from the same shadow/contact/source mechanism used in the
  flavor and triplet sectors, or be declared as an explicit conditional
  inverse-seesaw input.

next attempted nontrivial idea:
  Build a constrained Majorana source action whose F-terms force exactly the
  Veronese-plus-contact split above.  The useful target is not Veronese-only;
  it is a one-contact-lift theorem with a small number of coefficients that
  predicts the PMNS target rather than reconstructing it.

verification plan:
  1. Parameterize M_R by the five Veronese coefficients plus one contact
     coefficient and scan whether PMNS angles and splittings can be fitted
     without arbitrary inverse reconstruction.
  2. Add stability checks under small target perturbations and under the
     crossed-120 flavor source deformations.
  3. Feed the resulting PMNS-compatible coefficient card back into the d=5
     proton scan, because the same right-handed rotations enter the triplet
     Wilson tensors.

## 2026-05-10 19:39 Taipei heartbeat: Majorana contact sensitivity scan

status:
  The task remains incomplete as a first-principles GUT theorem.  The current
  Majorana/PMNS obstruction is now sharper: the source-consistent contact lift
  is not merely nonzero; its complex coefficient must be locked extremely
  accurately if the present PMNS benchmark is to be a prediction rather than an
  inverse-seesaw reconstruction.

files:
  - `code/scan_majorana_contact_sensitivity.py`
  - `output/majorana_contact_sensitivity/summary.json`
  - `output/majorana_contact_sensitivity/report.md`
  - `output/majorana_contact_sensitivity/real_scale_scan.csv`
  - `output/majorana_contact_sensitivity/phase_ring_scan.csv`
  - `paper/gut_framework.tex`
  - `code/audit_conditional_theorem_ledger.py`
  - `code/build_no_web_table_provenance_manifest.py`

mathematical derivation:
  With the source-consistent inverse-seesaw matrix decomposed as

    M_R / M_* = M_V + zeta K,

  scan only the contact source coefficient:

    M_R(s) = M_* (M_V + s zeta K),

  and the phase ring:

    M_R(delta) = M_* (M_V + exp(i delta) zeta K).

  For each point compute

    m_nu(s) = -m_D M_R(s)^{-1} m_D^T,

  Takagi-diagonalize `m_nu`, form

    U_PMNS(s) = U_e^\dagger U_nu(s),

  and compare the three PMNS angles and the two mass splittings to the local
  benchmark.  The first derivative follows directly from

    d m_nu / ds =
      m_D M_R^{-1} (dM_R/ds) M_R^{-1} m_D^T.

numerical verification:
  At Veronese-only `s=0`:

    max PMNS-angle residual = 5.869949e-01,
    max relative splitting residual = 9.974501e-01.

  At the reconstructed point `s=1, delta=0`:

    max PMNS-angle residual = 9.007795e-13,
    max relative splitting residual = 2.122693e-11.

  Under the deliberately loose local gate

    max angle residual <= 3.0e-2,
    max splitting relative residual <= 3.0e-1,

  the adaptive contact windows are

    scale half-width Delta s = 3.208798e-05,
    phase half-width Delta delta = 5.424119e-05 rad.

  Tight gate half-widths are even smaller:

    Delta s_tight = 9.937927e-07,
    Delta delta_tight = 1.799257e-06 rad.

  Local derivatives at `s=1`:

    d sin^2(theta12)/ds = -1.053664e2,
    d sin^2(theta13)/ds =  1.371521e2,
    d sin^2(theta23)/ds = -2.742749e2,
    d Delta m21^2/ds = -1.141987e-1 eV^2,
    d Delta m31^2/ds = -5.363911e1 eV^2.

ledger/provenance sync:
  Regenerated the conditional theorem ledger.  Counts are now:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 57,
    TUNED_FALLBACK: 3.

  Regenerated the no-web table provenance manifest after adding contact
  sensitivity keys:

    rows = 24,
    required TeX keys found = 17/17,
    provenance_manifest_complete = true.

current obstacle:
  The contact coefficient is a source-locking datum.  The PMNS result cannot
  be honestly promoted to a predictive theorem until the complex coefficient
  `zeta` is derived from an F-term/D-term/source equation, or until the paper
  explicitly labels it as an inverse-seesaw input.

next attempted nontrivial idea:
  Construct a minimal Majorana source-locking sector with one driving field
  imposing

    M_{11} - 2 M_{02} = -sqrt(3) zeta

  and one clockwork/hidden-radial equation fixing both `|zeta|` and `arg zeta`.
  The goal is not to remove the contact term, but to derive its complex value
  from the same shadow/contact mechanism already used elsewhere.

verification plan:
  1. Write an explicit holomorphic source superpotential for the contact
     coefficient and compute its F-flat Hessian.
  2. Check whether the source-locking modes form complete visible-threshold
     multiplets or hidden singlets.
  3. Re-run the PMNS and d=5 scans using the locked coefficient, including
     small source-sector perturbations to quantify stability.

## 2026-05-10 20:13 Taipei heartbeat: affine Majorana source-locking fallback

status:
  The task is still incomplete.  A minimal source-locking superpotential can
  lock the reconstructed contact coefficient without flat directions, but it
  does not predict the coefficient.  This is therefore a controlled conditional
  fallback, not a first-principles Majorana/PMNS theorem.

files:
  - `code/audit_majorana_source_locking_sector.py`
  - `output/majorana_source_locking_sector/summary.json`
  - `output/majorana_source_locking_sector/report.md`
  - `paper/gut_framework.tex`
  - `code/audit_conditional_theorem_ledger.py`
  - `code/build_no_web_table_provenance_manifest.py`

mathematical derivation:
  Introduce hidden-singlet chiral fields

    Z, P, Q, A, B, C

  and the renormalizable holomorphic source superpotential

    W_zeta = A (Z - P Q) + B (P - p0) + C (Q - q0),

  with

    p0 q0 = zeta.

  The F-terms are

    F_Z = A,
    F_P = -A Q + B,
    F_Q = -A P + C,
    F_A = Z - P Q,
    F_B = P - p0,
    F_C = Q - q0.

  At the vacuum

    Z = zeta, P = p0, Q = q0, A = B = C = 0,

  all F-terms vanish.  The holomorphic Hessian in field order
  `(Z,P,Q,A,B,C)` has nonzero entries

    W_ZA = 1,
    W_PA = -q0,
    W_QA = -p0,
    W_PB = 1,
    W_QC = 1.

  This tests whether the source lock has local moduli or visible threshold
  cost.  It does not derive zeta because p0 and q0 are target inputs.

numerical verification:
  Target contact coefficient:

    zeta = 0.1076472949 + 0.0736514853 i,
    |zeta| = 1.304319e-01,
    arg(zeta) = 6.000380e-01 rad.

  Balanced factorization:

    p0 = sqrt(|zeta|),
    q0 = sqrt(|zeta|) exp(i arg zeta).

  Results:

    |p0 q0 - zeta| = 3.103168e-17,
    F-term norm at vacuum = 3.103168e-17,
    Hessian rank = 6/6,
    Hessian singular values =
      (1.287467, 1.287467, 1.000000, 1.000000, 0.776719, 0.776719),
    visible threshold vector = (0,0,0).

  Inherited precision requirement from the previous PMNS sensitivity scan:

    loose scale half-width = 3.208798e-05,
    loose phase half-width = 5.424119e-05 rad,
    tight scale half-width = 9.937927e-07,
    tight phase half-width = 1.799257e-06 rad.

ledger/provenance sync:
  Regenerated the conditional theorem ledger:

    FAIL: 1,
    NO_GO: 13,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 57,
    TUNED_FALLBACK: 4.

  Regenerated the no-web table provenance manifest:

    rows = 26,
    required TeX keys found = 19/19,
    provenance_manifest_complete = true.

current obstacle:
  The affine source sector locks zeta after p0 and q0 are supplied.  It does
  not explain why the hidden sector should choose those constants to 1e-5
  precision.  The remaining first-principles problem is therefore the origin
  of the hidden source values, not the algebraic existence of a lock.

next attempted nontrivial idea:
  Replace the affine constants p0 and q0 by a hidden quotient/clockwork or
  modular-shadow extremum whose vacuum equations determine |zeta| and arg(zeta).
  The concrete target is a small hidden sector whose F/D equations yield
  p0 q0 = zeta while keeping all added fields hidden singlets or complete
  visible multiplets.

verification plan:
  1. Promote p0 and q0 to mesons of a hidden U(1) or U(N) quotient.
  2. Compute the combined F/D Hessian and check for phase/radial moduli.
  3. Quantify whether Planck-suppressed or soft-breaking perturbations move
     zeta by more than the PMNS tolerance windows.
  4. If no non-tuned source appears, keep the Majorana contact coefficient as
     an explicit conditional EFT datum in the final theorem map.

## 2026-05-10 21:15 Taipei heartbeat: hidden quotient no-go for Majorana zeta

status:
  The task remains incomplete.  The obvious hidden quotient upgrade of the
  affine source constants does not predict the PMNS-sensitive Majorana contact
  coefficient.  This is a useful no-go: the bottleneck has moved from "can zeta
  be locked?" to "what microscopic hidden dynamics chooses zeta?"

files:
  - `code/audit_majorana_hidden_quotient_origin.py`
  - `output/majorana_hidden_quotient_origin/summary.json`
  - `output/majorana_hidden_quotient_origin/report.md`
  - `paper/gut_framework.tex`
  - `code/audit_conditional_theorem_ledger.py`
  - `code/build_no_web_table_provenance_manifest.py`

mathematical derivation:
  Candidate A uses a hidden U(1) quotient with fields P,Q of charges (+1,-1):

    V_D = 1/2 (|P|^2 - |Q|^2)^2.

  At the balanced target |P|=|Q|=sqrt(|zeta|), D-flatness fixes only the radial
  difference.  The radial sum and the gauge-invariant product magnitude/phase
  remain moduli.

  Candidate B adds a holomorphic product constraint:

    W_X = X(PQ - Z),

  with scalar potential

    V = |PQ-Z|^2 + |XQ|^2 + |XP|^2 + |X|^2 + V_D.

  At Z=PQ=zeta and X=0, this enforces Z=PQ but transfers the modulus into Z
  rather than selecting its complex value.

numerical verification:
  Target inherited from the affine source-locking audit:

    |zeta| = 1.304319e-01,
    arg(zeta) = 6.000380e-01 rad.

  Required PMNS tolerance inherited from contact sensitivity:

    loose scale half-width = 3.208798e-05,
    loose phase half-width = 5.424119e-05 rad.

  Hessian audit:

    D-only hidden U(1): rank = 1/4,
    flat real directions before quotient = 3,
    predicts |zeta| = false,
    predicts arg(zeta) = false.

    D plus product constraint: rank = 5/8,
    flat real directions before quotient = 3,
    predicts |zeta| = false,
    predicts arg(zeta) = false.

    affine source fallback: F norm = 3.103168e-17,
    Hessian rank = 6/6,
    singular-value floor = 7.767189e-01,
    visible threshold vector = (0,0,0),
    but predicts |zeta| and arg(zeta) = false.

ledger/provenance sync:
  Regenerated the conditional theorem ledger:

    FAIL: 1,
    NO_GO: 14,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 57,
    TUNED_FALLBACK: 4.

  Regenerated the no-web table provenance manifest:

    rows = 28,
    required TeX keys found = 21/21,
    provenance_manifest_complete = true.

current obstacle:
  Pure D-flat quotient and the simple product constraint cannot determine
  zeta.  The affine source sector only hides zeta in p0,q0.  Therefore the
  remaining first-principles problem is a genuine hidden dynamics that fixes
  the product magnitude and phase to the PMNS tolerance window.

next attempted nontrivial idea:
  Test whether a discrete clockwork/monomial hidden sector can reduce the
  arbitrary complex input to integer charges plus one real scale, e.g.

    W = X(P^n Q^m - Lambda^{n+m}) + phase-locking term,

  and quantify whether the residual phase/scale still requires a spurion.  If
  yes, mark the Majorana contact coefficient as an irreducible conditional EFT
  datum.

verification plan:
  1. Enumerate small integer monomial constraints (n,m <= 6) and count how many
     independent dimensionful/phase inputs remain.
  2. Compute Hessian ranks and leftover moduli for each monomial clockwork.
  3. Compare the natural discretization precision to the PMNS windows
     Delta s ~ 3.2e-5 and Delta phase ~ 5.4e-5 rad.
  4. Promote only mechanisms that fix both modulus and phase without inserting
     a target complex spurion.

## 2026-05-10 22:29 Taipei heartbeat: monomial clockwork no-go for Majorana zeta

status:
  The task remains incomplete.  The discrete monomial/clockwork upgrade of the
  affine Majorana source does not predict the PMNS-sensitive contact
  coefficient in the small-order regime.  This is now a second no-go after the
  minimal hidden quotient: zeta can be locked, but its complex value is still a
  conditional input unless a stronger hidden dynamics is supplied.

files:
  - `code/audit_majorana_monomial_clockwork_origin.py`
  - `output/majorana_monomial_clockwork_origin/summary.json`
  - `output/majorana_monomial_clockwork_origin/report.md`
  - `code/audit_conditional_theorem_ledger.py`
  - `code/build_no_web_table_provenance_manifest.py`
  - `paper/gut_framework.tex`
  - `roadmap.md`

mathematical derivation:
  The test tries to replace the affine source constants p0,q0 by an integer
  monomial or root-of-unity phase.  In the real-coefficient version the phase
  is restricted to

    theta = 2 pi k / N.

  The target extracted from the inverse-seesaw Majorana contact is

    arg(zeta) = 6.000380e-01 rad.

  The loose PMNS-compatible phase window from the contact-sensitivity scan is

    |theta - arg(zeta)| <= 5.424119e-05 rad.

  Therefore small denominators predict zeta only if one of their discrete
  phases falls inside this very narrow interval.  A complex monomial
  coefficient would trivially fit the phase, but it would be exactly the
  spurion we were trying to remove.  A real monomial scale can fit the
  magnitude only by supplying a continuous hidden scale with relative tolerance

    Delta |zeta| / |zeta| = 3.208798e-05

  for the direct zeta-scale case, and tighter if zeta is generated as a power.

numerical verification:
  Target and tolerance:

    |zeta| = 1.304319e-01,
    arg(zeta) = 6.000380e-01 rad,
    loose scale half-width = 3.208798e-05,
    loose phase half-width = 5.424119e-05 rad.

  Small-denominator scan:

    best N <= 6: N = 6, k = 1,
    phase residual = 4.471595e-01 rad,
    residual / tolerance = 8.243910e3,
    passes loose phase = false.

  Extended denominator check:

    first acceptable denominator = 178,
    k = 17,
    phase residual = 4.147512e-05 rad.

  Verdict:

    small_monomial_clockwork_predicts_phase = false,
    small_monomial_clockwork_predicts_scale = false,
    complex_spurion_needed_for_small_orders = true,
    continuous_scale_tuning_required = true.

ledger/provenance sync:
  Regenerated the conditional theorem ledger after adding the monomial
  clockwork no-go row:

    FAIL: 1,
    NO_GO: 15,
    OPEN: 3,
    PASS: 3,
    PASS_AS_KINEMATICS: 1,
    PASS_CONDITIONAL: 57,
    TUNED_FALLBACK: 4.

  Regenerated the no-web table provenance manifest:

    rows = 30,
    required TeX keys found = 23/23,
    provenance_manifest_complete = true.

  Updated `paper/gut_framework.tex` with a new proposition
  "Monomial clockwork no-go for zeta" and updated the no-web table provenance
  proposition from 28/21 rows to 30/23 rows.  Recompiled
  `paper/gut_framework.pdf` twice.  The log check reports only the package
  line for rerunfilecheck and no unresolved warnings/errors.

current obstacle:
  The Majorana contact coefficient zeta is now the dominant first-principles
  obstruction.  Veronese-only Majorana fails; contact works only if zeta is
  fixed at the 1e-5 level.  The affine source locks it but imports p0,q0; the
  D/product quotient leaves moduli; small monomial/clockwork phases miss the
  target by thousands of tolerance widths.  The honest paper-level status is
  therefore a conditional Spin(10)/PSLT-inspired EFT with zeta as an explicit
  source datum, unless a non-minimal hidden sector is introduced and justified.

next attempted nontrivial idea:
  Stop spending cycles on small-order mechanisms.  There are only two credible
  routes left:

    Route 1: demote zeta formally to an irreducible conditional EFT datum and
    close the theorem ledger as a conditional model-building paper.

    Route 2: audit a high-denominator or modular-shadow hidden sector whose
    denominator/phase is not hand-picked, and include its UV/threshold cost in
    the same ledger.  The first real root that matches the loose phase window
    appears at N=178, so this route must justify why such a large denominator
    is natural rather than disguised tuning.

verification plan:
  1. Add a final "conditional theorem boundary" proposition: PSLT + Spin(10)
     EFT + source zeta implies the current branch; PSLT alone does not imply
     zeta or a unique GUT.
  2. If Route 2 is pursued, scan whether a modular-shadow or high-denominator
     clockwork sector predicts N=178-like phase locking without adding a
     complex spurion, and compute its added Dynkin/threshold/hidden-field cost.
  3. Otherwise, freeze zeta as a benchmark-card datum and move the remaining
     work to external cited input refresh and paper polishing.

## 2026-05-10 23:29 Taipei heartbeat: conditional theorem boundary certificate

status:
  The task remains incomplete as an unconditional first-principles GUT
  derivation, but the logical boundary is now explicit and machine-readable.
  The strongest supported result is a locally verified conditional Spin(10)
  EFT branch, not a PSLT-only theorem.  The dominant first-principles
  obstruction remains the microscopic origin of the Majorana contact
  coefficient zeta.

files:
  - `code/build_conditional_theorem_boundary.py`
  - `output/conditional_theorem_boundary/summary.json`
  - `output/conditional_theorem_boundary/report.md`
  - `code/build_no_web_table_provenance_manifest.py`
  - `output/no_web_table_provenance_manifest/summary.json`
  - `paper/gut_framework.tex`
  - `roadmap.md`

mathematical theorem boundary:
  The locally supported implication is:

    Spin(10) EFT
    + CP1/O(2) family kinematics
    + constrained/composite 54/210 threshold sector
    + source-fixed Majorana zeta
    + clockwork-rescued local d=5 proxy
      => current internally checked benchmark branch.

  The stronger implication is not established:

    PSLT/Another Physics alone
      => unique GUT with predicted zeta.

  In theorem language, zeta is not derived by the current first-principles
  inputs.  It is a conditional source datum unless a new hidden dynamics fixes
  both its magnitude and phase without importing a complex spurion.

numerical verification:
  The boundary certificate reads the current claim ledger:

    claim rows = 84,
    PASS_CONDITIONAL = 57,
    NO_GO = 15,
    TUNED_FALLBACK = 4,
    OPEN = 3.

  It returns:

    unconditional_first_principles_gut_derived = false,
    conditional_spin10_eft_branch_verified_locally = true,
    majorana_zeta_is_current_conditional_datum = true.

  The zeta-specific numerical obstruction is:

    contact fraction = 1.304275e-01,
    loose scale half-width = 3.208798e-05,
    loose phase half-width = 5.424119e-05 rad,
    best small monomial phase residual = 4.471595e-01 rad,
    first acceptable root denominator = 178.

provenance sync:
  Added the theorem-boundary artifact to the no-web table provenance manifest.
  To avoid a hash cycle, the theorem-boundary certificate records only the
  provenance manifest's completeness boolean, not its hash.

  Regenerated the no-web table provenance manifest:

    rows = 32,
    required TeX keys found = 25/25,
    provenance_manifest_complete = true.

  Updated `paper/gut_framework.tex` with a "Conditional theorem boundary"
  proposition and updated the no-web provenance proposition from 30/23 rows to
  32/25 rows.

current obstacle:
  The framework is now logically organized enough that further hidden-sector
  invention should be optional, not automatic.  Without a principled origin for
  zeta, the correct paper posture is "conditional EFT model-building branch".

next attempted nontrivial idea:
  If continuing physics development, only pursue a high-denominator or
  modular-shadow hidden sector if it predicts the phase denominator and scale
  without a target spurion and with acceptable threshold/UV cost.  Otherwise,
  freeze zeta as an explicit benchmark-card datum and spend the next work cycle
  on paper hygiene: theorem-scope wording, cited input refresh, and final
  reproducibility appendix.

verification plan:
  1. Recompile the paper and check for TeX warnings/errors after the boundary
     proposition.
  2. If Route 2 is chosen later, build a separate high-denominator/modular
     audit with the same ledger discipline.
  3. If Route 1 is chosen, mark the long-running heartbeat complete after the
     final conditional-EFT summary is written.

## 2026-05-12 Taipei rewrite: face-projection / Rubik covariance paper pass

status:
  Backed up the prior manuscript and rewrote the main theory framing around
  the "element labels as projection faces" idea.  This is a conceptual rewrite
  of the paper's spine, not a new numerical fit.  The strict conclusion remains
  unchanged: the framework is a locally verified conditional Spin(10) EFT
  branch, not an unconditional PSLT-only GUT proof.

backup:
  - `backups/2026-05-12_rubik_rewrite/gut_framework.tex`
  - `backups/2026-05-12_rubik_rewrite/gut_framework.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `roadmap.md`

new framing:
  The core principle is now stated as a face-projection/Rubik covariance
  principle:

    hidden state x in S,
    visible label e_theta(x) = pi_theta(x),
    e_theta(g x) = e_{g^{-1} theta}(x),

  with legal rotations required to preserve anomaly cancellation, index,
  chirality, charge normalization, and representation dimension.

  The physical translation is:

    one family = one Spin(10) half-spinor 16,
    six visible faces = Q, L, u^c, d^c, nu^c, e^c,
    three families = H^0(CP1,O(2)) protected sections.

mathematical additions to TeX:
  - Changed title to:
      A Face-Projection Spin(10) Effective Framework:
      Rubik Covariance, O(2) Family Geometry, Majorana Contact, and
      Proton-Safe Thresholds.
  - Rewrote the abstract to state the projection principle and conditional
    theorem boundary.
  - Rewrote the introduction so low-energy multiplets are projections rather
    than primitive labels.
  - Added Definition: Face-projection or Rubik covariance principle.
  - Added Proposition: Minimal six-face family object.
  - Added Proposition: Index-protected multiplicity.
  - Updated A1-A3 in the branch map:
      A1 = Spin(10) face-projection EFT with complete 16_i,
      A2 = H^0(CP1,O(2)) protects three copies of the 16 cubie,
      A3 = type-I seesaw with source-fixed Majorana contact.
  - Updated the no-go section's ledger counts to the current 84-row ledger.
  - Added final corollary: Face-projection status of the theory.

numerical/theorem status preserved:
  - claim rows = 84,
  - NO_GO = 15,
  - OPEN = 3,
  - PASS_CONDITIONAL = 57,
  - zeta contact fraction = 1.304275e-01,
  - zeta loose phase half-width = 5.424119e-05 rad,
  - small monomial/clockwork best residual = 4.471595e-01 rad,
  - first acceptable root denominator = 178.

verification:
  - Recompiled `paper/gut_framework.pdf` twice.
  - `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun" gut_framework.log`
    reports only the package line for `rerunfilecheck`.
  - Regenerated no-web table provenance manifest:
      rows = 32,
      required TeX keys found = 25/25,
      provenance_manifest_complete = true.

current obstacle:
  The rewrite improves the conceptual route from "periodic-table labels" to
  "high-dimensional projection faces", and it gives a cleaner motivation for
  Spin(10) as the minimal six-face family object.  It still does not derive
  the Majorana contact zeta.  Therefore the honest final status remains:

    face-projection motivation + conditional Spin(10) EFT
      != unconditional PSLT-only GUT proof.

next recommendation:
  Treat this as the new master narrative.  The next useful work should be a
  short "paper polish" pass: make all later section headings and terminology
  consistently use face-projection language, then decide whether to freeze zeta
  as a benchmark-card datum or start a separate high-denominator/modular-shadow
  hidden-sector project.

## 2026-05-13 Taipei trim: lean face-projection main paper

status:
  The main TeX has now been fully rewritten as a short face-projection theory
  note.  The long appendix/audit/calculation material has been removed from the
  main paper because it is not yet clear which pieces genuinely continue the
  face-projection architecture.  The old long version is preserved in backup,
  so this is a main-text trim rather than a loss of project history.

backup:
  - `backups/2026-05-13_trim_appendix/gut_framework.tex`
  - `backups/2026-05-13_trim_appendix/gut_framework.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `roadmap.md`

new main-paper structure:
  - Title:
      A Face-Projection Spin(10) Effective Framework: Rubik Covariance,
      O(2) Family Geometry, and the Conditional GUT Boundary.
  - Abstract now states the core boundary:
      PSLT/Another Physics do not yet prove a unique GUT unconditionally;
      the surviving result is a conditional Spin(10) face-projection EFT.
  - Section 1 defines the face-projection/Rubik covariance principle:
      e_theta(x)=pi_theta(x),
      e_theta(gx)=e_{g^{-1}theta}(x),
      with legal rotations preserving anomaly cancellation, index, chirality,
      charge normalization, and representation dimension.
  - Section 2 proves the six visible faces of one Spin(10) 16:
      Q, L, u^c, d^c, nu^c, e^c,
      and checks k_Y=5/3.
  - Section 3 proves the protected three-family carrier:
      H^0(CP1,O(2)) has dimension 3 and H^1(CP1,O(2))=0.
  - Section 4 keeps only the conditional EFT branch data A1-A7 and the
      current zeta boundary.
  - Section 5 states the theorem boundary:
      face-projection motivation + conditional Spin(10) EFT
        != unconditional PSLT-only GUT proof.
  - Section 6 explicitly defers the removed detailed calculations until their
      relation to the new architecture is decided.

verification:
  - Recompiled `paper/gut_framework.pdf`.
  - Current PDF length is 6 pages.
  - TeX log check reports no warnings/errors except the normal package line for
      `rerunfilecheck`.
  - Important caveat: the older no-web table provenance manifest was generated
      for the long TeX.  It should not be treated as a verification of the lean
      six-page TeX unless the manifest script is updated for lean-paper mode.

current obstacle:
  The conceptual architecture is now cleaner: low-energy particle labels are
  projections of a high-dimensional Spin(10) half-spin object, and three copies
  are protected by the O(2) index on CP1.  The hard physics obstacle remains the
  same: the Majorana contact parameter zeta and any hidden-sector origin for it
  are not derived from first principles.

next recommendation:
  Freeze this lean draft as the master narrative.  Next work should not restore
  the large appendices wholesale.  Instead, add back only one short conceptual
  appendix at a time, starting with a D5 half-spin weight-hypercube appendix
  that makes the Rubik/face-projection picture mathematically explicit.  After
  that, decide whether zeta is a benchmark-card datum or a separate hidden-sector
  derivation project.

## 2026-05-13 Taipei review-tighten pass

status:
  Applied the reviewer-style tightening pass to the lean face-projection draft.
  The result keeps the same conceptual architecture but lowers overstrong
  theorem language and makes the mathematical objects sharper.

backup:
  - `backups/2026-05-13_review_tighten/gut_framework.before_review_tighten.tex`
  - `backups/2026-05-13_review_tighten/gut_framework.before_review_tighten.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `roadmap.md`

changes made:
  - Replaced the loose state-to-label face map with a representation-theoretic
    restriction functor:
      Pi_theta(R)=Res^G_{H_theta} R,
      H_{g theta}=g H_theta g^{-1}.
  - Renamed the Spin(10) claim as a "standard minimal six-face realization"
    rather than a uniqueness theorem.
  - Tightened the assumptions:
      simple compact GUT candidate,
      one irreducible complex chiral representation,
      exactly one SM family plus nu^c,
      no chiral SM exotics,
      minimal representation dimension.
  - Made the hypercharge trace convention explicit:
      trace over one left-handed 16,
      included the nu^c zero hypercharge contribution,
      stated Tr_2 T_3^2=1/2.
  - Promoted the Spin^c/Dolbeault relation into the family theorem:
      D_R=sqrt(2)(dbar_L+dbar_L^\dagger),
      ker D_R^+=H^0(CP1,O(2)),
      ker D_R^-=H^1(CP1,O(2)).
  - Moved the no-extra-unpaired-SM-charged-sector condition into the theorem
    assumptions, not only the remark.
  - Added the missing transvectant lemma:
      Sym^2(Sym^2 C^2)=Sym^4 C^2 + Sym^0 C^2,
      and K_tr is the unique SL(2)-invariant contact direction up to phase.
  - Recast the zeta statement as a benchmark datum, not a proposition:
      M_R/M_*=M_V+zeta K_tr,
      zeta=0.1076472949+0.0736514853 i,
      contact fraction=1.304275e-01.
  - Recast the zeta-origin failure as a limited negative scan for simple
    hidden-quotient/small-monomial-clockwork classes, not a general no-go.
  - Changed the final theorem into a conditional benchmark boundary:
      D_face + D_index + D_zeta + D_companion_audits
        leadsto one conditional benchmark branch,
      while PSLT/Another Physics alone still does not imply a unique GUT with
      predicted zeta.
  - Deferred full flavor, full d=5 proton decay, threshold, and source-sector
    numerical checks to companion artifacts instead of implying the lean note
    contains the full audit package.

verification:
  - Recompiled `paper/gut_framework.pdf`.
  - Current PDF length is 7 pages.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

current obstacle:
  The lean note is now suitable as a conditional face-projection model-note
  skeleton.  It still does not derive zeta from first principles, and it still
  does not contain full flavor/proton/threshold/source-sector reproducibility.

next recommendation:
  Do not restore the old calculation appendices yet.  The next targeted
  addition should be a short D5 half-spin weight-hypercube appendix, because it
  directly supports the new face-projection/Rubik-covariance narrative without
  reopening the unresolved full phenomenology package.

## 2026-05-13 Taipei Route-B hidden zeta addition

status:
  Added a conservative Route-B hidden-sector candidate for the Majorana contact
  coefficient.  This is not presented as a PSLT-only derivation of zeta.
  Instead it is a new conditional hidden-sector theorem: a visible-singlet
  sterile transvectant messenger generates the already geometrically fixed
  contact tensor K_tr, while a hidden phase/radial sector supplies the complex
  coefficient.

backup:
  - `backups/2026-05-13_routeB_hidden_zeta/gut_framework.before_routeB_hidden_zeta.tex`
  - `backups/2026-05-13_routeB_hidden_zeta/gut_framework.before_routeB_hidden_zeta.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `code/verify_hidden_zeta_origin.py`
  - `output/hidden_zeta_origin/routeB_hidden_zeta_verification.json`
  - `roadmap.md`

TeX additions:
  - Added section:
      Route B: Hidden Transvectant Origin for zeta.
  - Added theorem:
      W_R/M_* = 1/2 N^T M_V N + lambda X^T N
                - 1/2 X^T K_tr^{-1} X
      integrates out X to give
      W_R^eff/M_* = 1/2 N^T (M_V + lambda^2 K_tr) N.
      Hence lambda^2=zeta generates the benchmark contact.
  - Added target square-root coupling:
      zeta = 0.1076472949 + 0.0736514853 i,
      |zeta| = 0.13043190325293763,
      arg zeta = 0.600038020318215,
      sqrt(zeta)=0.34502115743308964+0.10673473744038917 i,
      |sqrt(zeta)|=0.3611535729477664.
  - Added Z_178 phase option:
      theta_178 = 2*pi*17/178 = 0.6000794956295111,
      phase error = 4.1475311296e-05 rad,
      direct zeta_178 error = 5.4097037899e-06.
      This is inside the earlier loose phase window but is not an exact
      prediction unless an axion shift or continuous phase is added.
  - Added instanton/spurion option:
      zeta = A exp(-S_H+i theta_H),
      for A=1, S_H = -log|zeta| = 2.0369040025655094.
  - Added theorem:
      if X is visible singlet and Majorana-only, then delta b_visible=(0,0,0),
      Dirac Yukawas are unchanged, and d=5 colored-triplet Wilson tensors are
      unchanged.  If zeta equals the benchmark value, the seesaw replay is
      unchanged.
  - Updated theorem boundary:
      D_face + D_index + (D_zeta or D_hidden_zeta) + companion audits
        leadsto one conditional benchmark branch.
      Route B improves the conditional explanation of zeta but is still not a
      PSLT-only prediction.

verification:
  - Added and ran `code/verify_hidden_zeta_origin.py`.
  - Wrote JSON output to:
      `output/hidden_zeta_origin/routeB_hidden_zeta_verification.json`.
  - Numerical checks:
      sqrt_zeta_squared_error = 1.3877787807814457e-17,
      K_tr symmetric error = 0,
      K_tr inverse check = 3.8459253727671276e-16,
      messenger_contact_error = 1.2018516789897274e-17,
      phase_error = 4.147531129616855e-05,
      zeta_178_error = 5.409703789992232e-06,
      delta_b_visible = [0,0,0].
  - Recompiled `paper/gut_framework.pdf` twice.
  - Current PDF length is 9 pages.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

current obstacle:
  Route B derives the tensor insertion zeta K_tr after adding a hidden sterile
  messenger and phase/radial data.  The remaining open issue is no longer the
  contact tensor direction; it is whether the hidden phase/radial sector can be
  made microscopic, symmetry-protected, and compatible with the full flavor,
  proton, threshold, and UV audits.

next recommendation:
  Audit the Route-B hidden sector as a standalone module before restoring old
  phenomenology appendices.  The minimal next check is to promote the
  Majorana-only selection rule into an explicit charge table/superpotential and
  verify that it forbids X insertions into colored-triplet and Dirac-Yukawa
  operators while allowing the sterile transvectant messenger term.

## 2026-05-14 Taipei review of Route-B hidden zeta gaps

status:
  Reviewed the 2026-05-13 Route-B hidden-zeta addition.  The mechanism is
  algebraically sound as a conditional hidden-sector completion: integrating
  out a visible-singlet sterile transvectant messenger can generate
  zeta K_tr, and the visible-singlet/Majorana-only assumptions protect the
  one-loop visible threshold, Dirac Yukawas, and d=5 triplet Wilson tensors at
  the matching level.  It is still not a PSLT-only derivation.

review verdict:
  Correct and should be kept:
    - Face projection as representation restriction Res^G_{H_theta}R.
    - Spin^c CP1/O(2) index theorem for three protected families.
    - K_tr as the unique SL(2)-invariant second-transvectant contact direction.
    - Schur-complement derivation:
        W_R/M_* = 1/2 N^T M_V N + lambda X^T N
                  - 1/2 X^T K_tr^{-1} X
        => M_R/M_* = M_V + lambda^2 K_tr.
    - Treating zeta as conditional hidden-sector data unless a microscopic
      phase/radial sector is supplied.
    - Z_178 as a candidate/approximation, not a prediction.

must-fix items to add to theory:
  1. Representation precision for X:
     - Define V=H^0(CP1,O(2)).
     - Use N in V after B-L breaking.
     - Use X in V^vee, not ambiguously in V.
     - Treat K_tr in Sym^2 V^vee and K_tr^{-1} in Sym^2 V.
     - Rewrite the messenger theorem coordinate-free:
         X(N) is the natural pairing,
         X K_tr^{-1} X is a scalar.

  2. Post-B-L projector and Majorana-only selection rule:
     - Define N_i=P_{nu^c}(16_i) after B-L breaking.
     - State explicitly that Route B is a post-B-L-breaking EFT matching
       theorem unless a full Spin(10)-covariant source representation is added.
     - Add a minimal hidden charge/projector table allowing the Majorana
       messenger terms but forbidding X insertions into Dirac Yukawa and
       colored-triplet operators.
     - Required allowed terms:
         lambda X(N),
         X K_tr^{-1} X,
         N^T M_V N.
     - Required forbidden structures:
         X times Dirac-Yukawa operators,
         X times colored-triplet QQ/QL/UE/UD operators,
         X leakage into the triplet source sector.

  3. Canonical-normalization/radiative-stability audit:
     - The theorem as written is a tree-level matching theorem after canonical
       normalization, not yet a complete radiative-stability theorem.
     - Estimate:
         delta Z_N ~ |lambda|^2/(16*pi^2) log(M_*/M_X),
         |lambda|^2/(16*pi^2) ~ 8.26e-4.
     - This is small but larger than the loose PMNS phase/splitting windows, so
       the paper must either:
         a. define the benchmark in the already canonical-normalized matching
            basis, or
         b. run an explicit post-matching seesaw replay with delta Z_N included.

  4. Minimal numerical table for Route B:
     - Add a compact table to Section 5 using output from
       `code/verify_hidden_zeta_origin.py`:
         |zeta| = 0.13043190325293763,
         arg zeta = 0.600038020318215,
         sqrt(zeta)=0.34502115743308964+0.10673473744038917 i,
         |lambda|=0.3611535729477664,
         messenger_contact_error=1.2018516789897274e-17,
         K_tr inverse check=3.8459253727671276e-16,
         Z_178 phase error=4.147531129616855e-05,
         zeta_178_error=5.409703789992232e-06,
         instanton action=2.0369040025655094.

priority:
  P0:
    - Fix the X in V^vee representation statement and rewrite Theorem 5.1 in
      coordinate-free form.
    - Add the post-B-L projector P_{nu^c}(16) statement.
  P1:
    - Add explicit Majorana-only charge/projector table.
    - Add canonical-normalization caveat or audit formula.
  P2:
    - Add the compact Route-B numerical verification table.
    - Decide whether Z_178 remains a side observation or is moved to an
      appendix/remark.

next recommendation:
  Implement P0 first in the TeX.  The current Route-B mechanism is defendable,
  but the X in V versus V^vee issue is a genuine mathematical precision gap.
  After that, build the selection-rule table as a standalone hidden-sector
  module before restoring full flavor/proton/threshold appendices.

## 2026-05-14 Taipei P0 Route-B precision fix

status:
  Completed the P0 precision fix for the Route-B hidden transvectant messenger.
  The messenger theorem is now coordinate-free and explicitly post-B-L.  The
  previous ambiguity "X as a family vector" is removed.

backup:
  - `backups/2026-05-14_p0_routeB_precision/gut_framework.before_p0_precision.tex`
  - `backups/2026-05-14_p0_routeB_precision/gut_framework.before_p0_precision.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `roadmap.md`

TeX changes:
  - Added coordinate-free type information:
      V=H^0(CP1,O(2)),
      K_tr in Sym^2 V^vee,
      K_tr^{-1} in Sym^2 V,
      K_tr^2=(1/3)I and K_tr^{-1}=3K_tr in the displayed orthonormal basis.
  - Rewrote the Route-B theorem as:
      N_i=P_{nu^c}(16_i), N in V,
      X in V^vee,
      W_R/M_* = 1/2 M_V(N,N) + lambda X(N)
                - 1/2 K_tr^{-1}(X,X),
      with M_V,K_tr in Sym^2 V^vee and K_tr^{-1} in Sym^2 V.
  - Preserved matrix notation only as a chosen-basis representation:
      1/2 N^T M_V N + lambda X^T N - 1/2 X^T K_tr^{-1}X.
  - Rewrote the proof using maps:
      K_tr^{-1}: V^vee -> V,
      K_tr: V -> V^vee,
      F_X=0 gives X=lambda K_tr N,
      W_eff/M_* = 1/2 (M_V+lambda^2 K_tr)(N,N).
  - Added a "Scope of the projector" remark:
      N_i=P_{nu^c}(16_i) is a post-B-L EFT variable.
      A Spin(10)-singlet messenger cannot select nu^c before breaking without
      an additional Majorana source/projector.  Route B therefore assumes a
      post-breaking projector, e.g. from a (B-L)=-2 source sector.
  - Updated the Majorana-only theorem to explicitly use the same post-B-L
      projector N_i=P_{nu^c}(16_i).

verification:
  - Recompiled `paper/gut_framework.pdf`.
  - Current PDF length is 10 pages.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

remaining items:
  P1:
    - Build an explicit Majorana-only hidden charge/projector table.
    - Add canonical-normalization/radiative-stability caveat or audit formula.
  P2:
    - Add compact Route-B numerical verification table in the TeX.
    - Decide whether the Z_178 phase option remains in the main text or moves
      to a remark/appendix.

next recommendation:
  Start P1 with the selection-rule table.  The clean target is a minimal charge
  assignment that allows X(N) and X K_tr^{-1} X but forbids X insertions into
  Dirac Yukawa and colored-triplet operators.

## 2026-05-14 Taipei P1 Route-B selection and canonical audit

status:
  Completed P1 for the Route-B hidden transvectant messenger.  The paper now
  contains an explicit post-B-L R-selection rule that keeps X Majorana-only,
  plus a canonical-normalization caveat for loop-level stability.

backup:
  - `backups/2026-05-14_p1_routeB_selection/gut_framework.before_p1_selection.tex`
  - `backups/2026-05-14_p1_routeB_selection/gut_framework.before_p1_selection.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `code/verify_hidden_zeta_origin.py`
  - `output/hidden_zeta_origin/routeB_hidden_zeta_verification.json`
  - `roadmap.md`

TeX changes:
  - Added Lemma: Minimal post-B-L R-selection rule.
  - Selection rule:
      R(W)=2,
      R(16_i components, including N_i)=1,
      R(H,T,Tbar)=0,
      R(X in V^vee)=1,
      R(lambda,M_V,K_tr^{-1})=0.
  - Allowed by R-charge:
      X(N): 1+1=2,
      K_tr^{-1}(X,X): 2,
      M_V(N,N): 2,
      baseline Dirac/triplet Yukawa terms 16_i16_jH and 16_i16_jT: 2.
  - Forbidden:
      Any Dirac or colored-triplet operator with at least one extra X has
      R=2+m with m>=1, so it is forbidden while the R-selection rule is exact.
  - Updated the Majorana-only theorem to cite the R-selection lemma.
  - Added Observation: Canonical basis caveat.
      The theorem is tree-level matching in a canonical-normalized EFT basis.
      A microscopic messenger interval can generate
        delta Z_N ~ |lambda|^2/(16 pi^2) log(M_*/M_X).
      For the benchmark,
        |lambda|^2/(16 pi^2)=8.259697e-04.
      If non-negligible, the companion seesaw replay must use
        Y_nuD -> Y_nuD (1 - delta Z_N/2),
        M_R -> M_R - 1/2(delta Z_N^T M_R + M_R delta Z_N).

verification:
  - Updated and reran `code/verify_hidden_zeta_origin.py`.
  - Added JSON field:
      one_loop_lambda_factor = 0.000825969676394408.
  - Recompiled `paper/gut_framework.pdf`.
  - Current PDF length is 11 pages.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

remaining issue:
  P1 now supplies a clean selection rule at the holomorphic matching level.
  It is not yet a full microscopic UV completion.  If R-breaking spurions are
  introduced later, they must be checked not to regenerate X insertions into
  Dirac or colored-triplet operators.

next recommendation:
  Move to P2: add a compact numerical Route-B verification table to the TeX,
  using the JSON output, and decide whether the Z_178 phase option should stay
  in the main text or be demoted to a remark/appendix.

## 2026-05-14 Taipei P2 Route-B numerical table and phase demotion

status:
  Completed P2 for the Route-B hidden transvectant messenger.  The lean paper
  now contains a compact numerical verification table for the hidden-zeta
  mechanism, and the Z_178 phase option has been demoted to a side diagnostic
  rather than a premise of the benchmark branch.

backup:
  - `backups/2026-05-14_p2_routeB_table/gut_framework.before_p2_table.tex`
  - `backups/2026-05-14_p2_routeB_table/gut_framework.before_p2_table.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `roadmap.md`

TeX changes:
  - Renamed the Z_178 item from a phase option to a phase diagnostic.
  - Added a sentence clarifying that the Route-B theorem uses only
      lambda^2 = zeta,
    while the Z_178 value is retained only as a side diagnostic.
  - Added Table `tab:routeb-zeta-checks` with the compact numerical checks:
      |zeta_tar| = 0.13043190325293763,
      arg zeta_tar = 0.600038020318215,
      sqrt(zeta_tar) = 0.345021 + 0.106735 i,
      |lambda| = 0.3611535729477664,
      ||lambda^2 K_tr - zeta K_tr||_F = 1.20185e-17,
      ||K_tr(3K_tr)-I||_F = 3.84593e-16,
      2pi(17/178)-arg zeta = 4.14753e-05 rad,
      |zeta_178-zeta_tar| = 5.40970e-06,
      -log|zeta| = 2.0369040025655094,
      |lambda|^2/(16pi^2) = 8.25970e-04,
      Delta b_vis^X = (0,0,0).
  - The table caption explicitly states that these are the lean algebraic
    Route-B checks, not a substitute for the deferred full flavor/proton/
    threshold companion audits.

verification:
  - Recompiled `paper/gut_framework.pdf` twice after adding the table.
  - Current PDF length is 11 pages.
  - `pdfinfo` reports file size 303283 bytes.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

current theory status after P0-P2:
  Route B is now a clean conditional matching mechanism:
    face projection + Spin^c O(2) index + unique transvectant K_tr
    + post-B-L hidden messenger X in V^vee
    + Majorana-only R-selection
    + canonical-basis caveat
    + compact numerical checks
    => M_R/M_* = M_V + zeta K_tr
  without claiming a PSLT-only first-principles derivation.

remaining hard problems:
  - The radial/phase origin of lambda itself is still conditional.  Z_178 is only
    a diagnostic; an axion, instanton, or explicit hidden potential is needed for
    a dynamical origin.
  - The lean paper still defers full flavor fitting, full d=5 proton-decay
    Wilson tensors, and threshold/proton companion audits.
  - If R-breaking spurions are added later, they must be checked not to
    regenerate X insertions into Dirac or colored-triplet sectors.

next recommendation:
  Choose one of two clean next branches:
    A. Hidden-origin branch: write an explicit phase/radial stabilization
       potential for lambda and verify it preserves the Majorana-only selection
       rule.
    B. Representation branch: add a short D5 half-spin weight-hypercube appendix
       explaining the six SM faces as the Spin(10) 16 branching, while keeping it
       motivational rather than a uniqueness theorem.

## 2026-05-14 Taipei review of post-P2 wording boundaries

status:
  The external 2026-05-13-style assessment is broadly correct.  The current
  draft is defensible as a conditional face-projection Spin(10) EFT note with
  a Route-B hidden-transvectant origin for zeta, but it still must not be
  presented as a PSLT-only first-principles GUT proof, a full flavor fit, or a
  complete proton-decay/threshold closure.

confirmed strengths:
  - Spin(10) is now framed as a standard minimal single-object six-face
    realization, not as an overstrong uniqueness theorem.
  - The CP1 O(2) Spin^c index theorem is the strongest proved part:
      dim H^0(CP1,O(2))=3 and H^1(CP1,O(2))=0
    under the stated branch assumptions.
  - K_tr is now fixed as the unique SL(2)-invariant second-transvectant
    contact direction, not an arbitrary fitted tensor.
  - Route B gives a coordinate-free Schur-complement mechanism with
      N in V, X in V^vee, lambda^2=zeta,
    producing M_R/M_* = M_V + zeta K_tr.
  - The post-B-L R-selection rule and canonical-basis caveat correctly limit
    the claim to tree-level/canonical matching plus visible threshold silence.
  - The compact Route-B numerical table is useful and correctly scoped.

P3 wording fixes to implement in TeX:
  P3.1 Abstract:
    - Update the statement "central unresolved representation-level datum is
      zeta".
    - More precise wording:
        K_tr is fixed by the unique second transvectant and Route B can
        generate zeta K_tr by Schur complement; the remaining open datum is
        the hidden phase/radial sector supplying lambda^2=zeta, plus companion
        flavor/proton/threshold/source audits.
  P3.2 Theorem 5.7 title and first sentence:
    - Rename from "Majorana-only hidden contact does not disturb visible
      audits" to a narrower title such as:
        "Tree-level canonical matching silence of Majorana-only hidden contact".
    - Add explicitly:
        In a canonical-normalized post-B-L EFT at the matching scale...
      so the theorem is not read as an all-orders radiative-stability claim.
  P3.3 R-selection caveat:
    - Add a remark after Lemma 5.6:
        The U(1)_R is used as a holomorphic EFT matching selection rule, not
        yet as an anomaly-free microscopic symmetry.
    - If upgraded to a UV symmetry later, audit continuous/discrete R anomalies
      or introduce spurions whose R-breaking insertions do not regenerate
      X-decorated Dirac or colored-triplet operators.
  P3.4 Projector visibility:
    - Make the statement near Theorem 5.1/Scope remark more explicit:
        P_{nu^c} is an input of the post-B-L EFT matching, not derived by
        face projection alone.
  P3.5 Deferred-audit ledger:
    - Add a small checklist before the conclusion/deferred section:
        proved in lean note:
          Spin(10) face skeleton, O(2) index, K_tr uniqueness,
          Route-B Schur complement;
        conditional inputs:
          C=CP1, R=p_+ + p_-, Spin^c carrier, no unpaired exotics,
          P_{nu^c}, R-selection rule, hidden phase/radial sector;
        deferred audits:
          full flavor fit, full d=5 proton-decay Wilson tensors,
          threshold/proton companion audits, source-sector UV origin;
        not claimed:
          PSLT-only unique GUT proof.

items not worth changing now:
  - Do not promote Z_178 from diagnostic to prediction unless a dynamical
    reason for N=178 and m=17 is derived.
  - Do not restore the long flavor/proton/threshold appendices until their
    assumptions are compatible with the lean Route-B skeleton.

recommended next action:
  Implement P3.1--P3.5 as a wording-boundary pass before choosing between the
  hidden-origin branch and the D5 half-spin representation appendix.  This is
  a low-risk pass and will make the draft harder to attack for overclaiming.

## 2026-05-14 Taipei P3 wording-boundary pass completed

status:
  Completed the P3 wording-boundary pass.  The paper now more accurately
  states what Route B does and does not solve: K_tr is fixed, zeta K_tr can be
  generated by hidden transvectant matching, but the phase/radial origin of
  lambda^2=zeta and the full companion phenomenology remain open.

backup:
  - `backups/2026-05-14_p3_wording_boundary/gut_framework.before_p3_wording.tex`
  - `backups/2026-05-14_p3_wording_boundary/gut_framework.before_p3_wording.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `roadmap.md`

TeX changes:
  - Abstract:
      Replaced the older wording that zeta itself is the central unresolved
      representation-level datum.  The abstract now says K_tr is fixed by the
      unique SL(2)-invariant second transvectant, Route B generates
      zeta K_tr when lambda^2=zeta, and the remaining open datum is the hidden
      phase/radial sector that supplies this complex coefficient.
  - Projector scope:
      Made explicit that P_{nu^c} is an input of the post-B-L EFT matching
      problem, not a result of face projection alone.
  - R-selection rule:
      Added a remark that the U(1)_R assignment is a holomorphic EFT matching
      selection rule, not yet an anomaly-free microscopic symmetry.  Any UV
      implementation must use an anomaly-free continuous/discrete R symmetry
      or controlled spurions that do not regenerate X-decorated Dirac or
      colored-triplet operators.
  - Theorem 5.7:
      Renamed it to a narrower tree-level canonical matching theorem:
        Tree-level canonical matching silence of Majorana-only hidden contact.
      The first sentence now specifies a canonical-normalized post-B-L EFT at
      the matching scale.
  - Status ledger:
      Added a lean-note status ledger separating:
        proved in this lean note:
          Spin(10) face skeleton, Spin^c(CP1,O(2)) index,
          K_tr uniqueness, Route-B Schur complement;
        conditional inputs:
          C=CP1, R=p_+ + p_-, protected Spin^c carrier,
          no unpaired SM-charged exotics, P_{nu^c},
          R-selection rule, hidden phase/radial sector;
        deferred audits:
          full flavor fit, d=5 proton-decay Wilson tensors,
          threshold/proton companion audits, source-sector UV origin;
        not claimed:
          PSLT-only unique GUT proof.

verification:
  - Recompiled `paper/gut_framework.pdf` after the wording pass.
  - Current PDF length is 12 pages.
  - `pdfinfo` reports file size 305960 bytes.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

current theory status:
  The draft is now cleanly positioned as:
    conditional face-projection Spin(10) EFT note
    + Route-B hidden-transvectant origin for zeta K_tr
    + explicit non-claim of PSLT-only first-principles GUT proof.

next recommendation:
  The next nontrivial scientific step should be one of:
    A. Hidden phase/radial origin:
       build an explicit potential or nonperturbative spurion model for
       lambda, then verify it preserves the Majorana-only selection rule.
    B. D5 half-spin representation appendix:
       add the weight-hypercube derivation showing how the Spin(10) 16
       branches into the six Standard-Model faces, while keeping it as
       motivation/minimality rather than a full uniqueness theorem.

## 2026-05-14 Taipei Route-A main-paper boundary strategy

status:
  Adopted the conservative paper strategy suggested by the latest assessment:
  Route A should be the main paper, while Route B remains a conditional
  candidate mechanism or companion-paper seed.  Also fixed the theorem-boundary
  wording so companion audits cannot be read as optional when Route B is used.

backup:
  - `backups/2026-05-14_routeA_boundary_strategy/gut_framework.before_routeA_boundary.tex`
  - `backups/2026-05-14_routeA_boundary_strategy/gut_framework.before_routeA_boundary.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `roadmap.md`

TeX changes:
  - Theorem Boundary:
      Replaced the ambiguous wording
        "A1--A7 and companion audits, or Route-B hidden messenger..."
      with:
        face-projection data + O(2) index data + companion benchmark audits
        + either D_zeta or D_hidden_zeta.
      This makes companion audits common to both the benchmark-zeta path and
      the Route-B hidden-zeta path.
  - Theorem Boundary proof:
      Added that A1--A7 and benchmark-level checks are required whether zeta
      is supplied as a benchmark datum or generated through Route B.
      Route B improves the conditional explanation of zeta but does not
      replace companion phenomenological audits.
  - Lean-note status ledger:
      Rephrased the first row as "proved/constructed" so Route-B Schur
      complement is not conflated with the strictly proved index/representation
      parts.
      Converted the ledger to a compact ragged-right table to avoid overfull
      display-array layout.
  - Deferred calculations:
      Replaced the binary "next paper-level decision" wording with an explicit
      conservative strategy:
        Route A = main paper: conditional face-projection model note.
        Route B = candidate mechanism or companion: hidden origin for zeta.
      Added that promoting Route B to the main line would require restoring only
      flavor/proton/threshold/source calculations that survive the R-selection
      rule and canonical-normalization replay.

verification:
  - Recompiled `paper/gut_framework.pdf`.
  - Current PDF length is 12 pages.
  - `pdfinfo` reports file size 314689 bytes.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

current strategic state:
  Route A is now the recommended main-paper form:
    conditional face-projection Spin(10) EFT note with O(2) index and K_tr
    uniqueness.
  Route B is preserved as a valuable conditional mechanism:
    hidden transvectant origin for zeta K_tr,
  but it should not carry the main paper until hidden phase/radial dynamics,
  R-symmetry/spurion consistency, delta Z_N replay, and companion audits are
  restored.

next recommendation:
  For the Route-A main-paper path, the most useful next addition is the short
  D5 half-spin weight-hypercube appendix.  It strengthens the six-face
  representation story without dragging the draft back into unfinished
  hidden-sector phenomenology.

## 2026-05-15 Taipei D5 half-spin weight-hypercube appendix

status:
  Completed the short Route-A appendix deriving the six visible Standard-Model
  faces directly from the D5 half-spin weight hypercube.  This strengthens the
  representation-restriction skeleton of the main paper without reopening
  hidden-sector phenomenology.

backup:
  - `backups/2026-05-15_d5_half_spin_appendix/gut_framework.before_d5_appendix.tex`
  - `backups/2026-05-15_d5_half_spin_appendix/gut_framework.before_d5_appendix.pdf`

files changed:
  - `paper/gut_framework.tex`
  - `paper/gut_framework.pdf`
  - `code/verify_d5_half_spin_hypercube.py`
  - `output/d5_half_spin/d5_half_spin_weights.csv`
  - `output/d5_half_spin/d5_half_spin_summary.json`
  - `output/d5_half_spin/d5_half_spin_report.md`
  - `roadmap.md`

TeX changes:
  - Added Appendix A, "D5 Half-Spin Weight Hypercube".
  - Wrote the D5 root and half-spin convention:
      roots +/-e_a +/- e_b,
      half-spin weights 1/2 sum_a s_a e_a with an even number of minus signs.
  - Split the signs under the Pati-Salam face:
      D5 -> D3 + D2,
      Spin(6) x Spin(4) ~= SU(4)_C x SU(2)_L x SU(2)_R.
  - Proved the parity split:
      even color parity + even weak parity gives (4,2,1) with 8 weights;
      odd color parity + odd weak parity gives (bar4,1,2) with 8 weights.
      Hence 16 -> (4,2,1) + (bar4,1,2).
  - Added the SU(4)_C convention for B-L and the weak-sign convention for
      T_{3L}, T_{3R}.
  - Applied Y = T_{3R} + (B-L)/2 and listed the six visible faces:
      Q:      (3,2),      Y=1/6,   multiplicity 6
      L:      (1,2),      Y=-1/2,  multiplicity 2
      u^c:    (bar3,1),   Y=-2/3,  multiplicity 3
      d^c:    (bar3,1),   Y=1/3,   multiplicity 3
      nu^c:   (1,1),      Y=0,     multiplicity 1
      e^c:    (1,1),      Y=1,     multiplicity 1
  - Explicitly stated that the appendix is a weight-level explanation of the
      standard Spin(10) half-spinor branch, not a new uniqueness theorem.

verification:
  - Added and ran `code/verify_d5_half_spin_hypercube.py`.
  - Generated reproducibility artifacts in `output/d5_half_spin/`.
  - Script summary:
      number of half-spin weights = 16;
      Pati-Salam counts = (4,2,1):8 and (bar4,1,2):8;
      field multiplicities = Q:6, L:2, u^c:3, d^c:3, nu^c:1, e^c:1;
      all expected multiplicities match;
      Tr Y^2 = 10/3, Tr T3L^2 = 2, k_Y = 5/3.
  - Recompiled `paper/gut_framework.pdf`.
  - Current PDF length is 13 pages.
  - `pdfinfo` reports file size 335560 bytes.
  - Log scan:
      `rg "Warning|undefined|Overfull|Underfull|Error|Fatal|Rerun|does not exist"`
      reports only the normal `rerunfilecheck` package line.

current strategic state:
  Route A is stronger and still clean:
    face projection is supported by an explicit D5 half-spin weight appendix,
    the three-family count remains protected by the O(2) index theorem,
    and Route B remains only a conditional hidden-zeta mechanism or companion
    direction.

next recommendation:
  The next natural Route-A polish pass is not more hidden-sector physics, but
  final-paper hygiene:
    tighten the introduction/conclusion around the theorem boundary,
    make sure citations cover Spin(10), Pati-Salam, Riemann-Roch/Spin^c, and
    representation branching,
    and add a minimal reproducibility manifest pointing to the D5 verifier and
    any surviving companion-audit artifacts.
