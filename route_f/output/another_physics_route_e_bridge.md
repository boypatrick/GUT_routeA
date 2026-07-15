# Another Physics / Route-E Bridge Verification

Status: **mechanical_separation_checks_pass_no_physics_promotion** — `27/27` checks pass.

This is a fail-closed algebraic/numerical card.  A green result means only
that the displayed identities are reproducible.  It does **not** identify the
two theories and `physics_promotion_allowed=false`.

## Separation theorems

1. On the assumed carrier, `H^0(CP1,O(2))` has basis
   `(X^2,sqrt(2)XY,Y^2)`, dimension `3`, and Cartan weights `(+2,0,-2)`.
2. A nonzero quadratic section `aX^2+bXY+cY^2` has discriminant
   `Delta=b^2-4ac`: `Delta != 0` gives two distinct projective zeros, while
   `Delta=0` gives a double zero.  The identically zero section is a separate
   zero orbit and is outside this divisor statement.  This separates the
   regular-semisimple and nilpotent cases; it does not make every degree-two
   pattern the same physical object.
3. The pairing loaded from Route-E obeys
   `rho(A)^T K_tr+K_tr rho(A)=0` for `A=H,E,F` and `K_tr^2=I/3`.
   In Route-E spherical coordinates it also obeys the exact identity
   `B_Kill(x,x)=2 Delta(p)=2 sqrt(3) K_tr(x,x)`, while the standard
   projective matrix satisfies `A^2=Delta I/4`.  Thus the contact-null cone is
   the double-zero locus.  With the displayed generator normalization the CP1
   moment map assigns `+1/-1` to the two fixed points; their opposition, not
   that numerical magnitude, is invariant.  These are charge/orientation
   labels, not energies.
4. For `zeta=0.1076472949+0.0736514853i`,
   `arg(zeta)=0.600038020318215`,
   `cos^2(arg zeta)=0.681143440291929`, and the weight-one lift gives
   `cos^2(arg zeta/2)=0.912657073213319`.  Under complex `y`,
   `zeta'=(y^2/|y|^2)zeta`; this is weight-two covariance, not phase invariance.
5. For the explicitly chosen sextic Q-ball action, the existence window is
   `0.866025403784<omega<1 GeV`.  At `Q=10^6`, the leading thin-wall card gives
   `R=12.842415686 fm` and `E/Q=0.872678753982 GeV`; the independent
   finite-Q step-profile minimization gives `E/Q=0.872667433681 GeV < m`.
   Both are controlled approximations, not an exact radial-ODE solution.
6. The displayed oscillatory equation has continuous solution families
   `x^2-y^2=2 pi n` and `x^2+y^2=pi/2+2 pi n`.  For the circular family,
   `Delta r=2 pi/(r_(n+1)+r_n)~pi/r_n`; this proves a shrinking interference
   pitch, not atomic-orbital quantization.
7. The AP-E3 successor card proves an exact constrained-cell `k=+2` theorem
   and an anomaly-consistent mixed-WZW intermediate candidate (`26/26`), while
   the AP-E4 card solves the declared canonical Spin-c spectrum (`22/22`).
   The new nonlinear proxy (`18/18`) passes its vacuum, `B=1`, radial, and
   charged-scalar necessary gates but explicitly leaves lattice QC2D and the
   full Hessian open.  The global-UV card defines the bosonic non-extendible
   WZW phase and proves a pure-`SU(3)` no-go for the original charge-one
   scalar; its two spin-torsion signs and an equivalent all-scale embedding
   remain open.  The moduli-SQM card (`27/27`) independently derives a
   physical tangent fermion, but neither the mother-model canonical vacuum
   line nor its CPT antiparticle polarization.  It is also not the same
   charged-two-colour soliton, so the spatially transgressed AP-E3 `O(2)`
   line, gauge descent, and Route-E portal remain non-derived.

## Bridge axioms retained as non-derived

- **BA1_phase_identification** (`not_derived`): Identify an Another-Physics visibility phase delta with either arg(zeta) or the weight-one lift arg(zeta)/2.  These choices predict different cos^2 values and are not selected by Route-E.
- **BA2_bubble_identification** (`not_derived`): Identify the proposed energy bubble with the charge-Q soliton of the displayed complex-scalar action.  Q-ball existence does not derive this map.
- **BA3_two_center_identification** (`not_derived`): Identify a quadratic-phase circle/hyperbola morphology with the two-zero divisor of a CP1/O(2) section.  Shared degree-two algebra is insufficient to identify their state spaces or dynamics.
- **BA4_UV_embedding** (`open`): Embed the Q-ball scalar and its conserved U(1) charge in the Route-E Spin(10)/family sector while preserving gauge, anomaly, and threshold gates.

## Mechanical checks

- [PASS] `S0_AP_E3_AP_E4` — AP-E3 exact-cell and intermediate mixed-WZW audit is green but non-promoting: AP-E3 UV=26/26; exactly-two=True; k+2=True; full UV=False; portal=False
- [PASS] `S0_AP_E3_AP_E4` — AP-E4 canonical Spin-c spectrum is mathematically green but its physical origin is open: AP-E4=22/22; math=True; physical tangent fermion=False; physics closure=False
- [PASS] `S0_AP_E3_AP_E4` — charged two-colour nonlinear proxy closes only radial classical-EFT gates: proxy=18/18; radial=True; full3d=False; lattice=False
- [PASS] `S0_AP_E3_AP_E4` — non-extendible bosonic WZW is defined while spin torsion and pure-SU(3) equivalence stay open: global=38/38; char=True; U2 bundle=False; torsion=False; SU3 q=1=False; q=2 equivalent=False
- [PASS] `S0_AP_E3_AP_E4` — intrinsic N=2 SQM derives a tangent fermion but not the charged-two-colour composition: SQM=27/27; tangent=True; charged embedding=False; portal=False
- [PASS] `S1_CP1_O2` — H0(CP1,O(2)) has dimension three: h^0(O(d))=d+1 for d>=0; d=2 gives 3
- [PASS] `S1_CP1_O2` — Cartan weights are (+2,0,-2): weights=[2, 0, -2] in basis ['X^2', 'sqrt(2) X Y', 'Y^2']
- [PASS] `S2_quadratic_divisor` — nonzero discriminant gives two distinct projective zeros: Delta=4.0, roots=[(-1-0j), (1-0j)], residual=0.00e+00
- [PASS] `S2_quadratic_divisor` — a nonzero section with zero discriminant gives one double zero: Delta=0.0, double root=-0.0
- [PASS] `S3_Ktr` — locally reconstructed pairing matches the Route-E convention card: max matrix difference=0.00e+00
- [PASS] `S3_Ktr` — K_tr squared equals I/3: max|K_tr^2-I/3|=1.11e-16
- [PASS] `S3_Ktr` — K_tr is invariant under the full sl2 generator triple: residuals={"E": 0.0, "F": 0.0, "H": 0.0}
- [PASS] `S3_Ktr` — Route-E Killing norm equals twice the binary-quadratic discriminant: regular Delta=-0.120000, B(x,x)=-0.240000, |B-2Delta|=0.00e+00; null Delta=2.22e-16, null B=0.00e+00
- [PASS] `S3_Ktr` — standard projective matrix gives A^2=Delta I/4 and canonical B=2 Delta: max|A^2-Delta I/4|=0.00e+00; |4 Tr(A^2)-2Delta|=0.00e+00
- [PASS] `S3_Ktr` — the normalized CP1 moment map labels the two fixed points by plus/minus one: mu={"equator": 0.0, "p_minus": -1.0, "p_plus": 1.0}
- [PASS] `S4_zeta` — benchmark magnitude and argument are reproduced: canonical-card error=0.00e+00; |zeta|=0.1304319032529376, arg(zeta)=0.600038020318215
- [PASS] `S4_zeta` — both explicit cos-squared phase conventions are reproducible: cos^2(arg zeta)=0.6811434402919294; cos^2(arg zeta/2)=0.9126570732133188
- [PASS] `S4_zeta` — complex rescaling obeys weight-two covariance: |error|=2.80e-17, modulus error=2.78e-17, phase error=0.00e+00
- [PASS] `S4_zeta` — weight-two covariance is not complex-phase invariance: |zeta'-zeta|=0.094332 for arg(y)=0.37
- [PASS] `S5_Qball` — chosen sextic potential has a nonempty Q-ball frequency window: omega0^2=0.750000; window 0.866025403784<omega<1.0; factorization residual=5.03e-17
- [PASS] `S5_Qball` — analytic wall tension agrees with independent quadrature: sigma_analytic=0.125000000000, sigma_numeric=0.125000000000
- [PASS] `S5_Qball` — leading thin-wall benchmark reproduces the quoted radius and E/Q: R=65.081904460 GeV^-1=12.842415686 fm; E/Q=0.872678753982 GeV
- [PASS] `S5_Qball` — finite-Q step-profile variational cross-check stays inside the window and below free-particle threshold: R_var=64.971265436 GeV^-1, omega_var=0.870457184321, E/Q=0.872667433681 GeV
- [PASS] `S6_quadratic_phase` — circle family x^2+y^2=pi/2+2pi n solves f(x)=f(y): maximum sampled residual=1.40e-14
- [PASS] `S6_quadratic_phase` — hyperbola family x^2-y^2=2pi n solves f(x)=f(y): maximum sampled residual=2.61e-14
- [PASS] `S6_quadratic_phase` — radial pitch obeys the exact rationalized identity and pi/r asymptotic: n=1000000, Delta r=1.253313667348e-03, identity error=2.53e-14, asymptotic relative error=2.50e-07
- [PASS] `S7_boundary` — every cross-theory identification remains an explicit bridge axiom: 4 bridge axioms retained as non-derived

## Promotion boundary

- `physics_promotion_allowed=false`.
- A physical bridge requires an action-level phase map, a gauge/anomaly-safe
  embedding compatible with the exactly-two dressing, a full gauge-theory
  soliton/fluctuation analysis, a same-mother-model SQM/WZW composition, a
  degree-one portal, and reruns of Route-E's flavor, threshold, proton, and
  cosmology gates.
