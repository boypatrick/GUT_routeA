# AP-E4 same-soliton Callias/descent audit

- Passed: **27/27**
- Conditional sufficient composition criterion established: **True**
- Same-model physical tangent fermion derived: **False**
- Same-model Callias c1 computed: **False**
- Gauge-basic WZW descent proved: **False**
- Degree-one Route-E portal built: **False**
- Physics promotion: **False**

## Quantitative regression

- Discrete Berry Chern number of the conditional tangent/Callias line: `1.999999999999981`.
- Finite template/control kernel residual: `4.561e-16`; gap: `1.000000000000`.
- Conditional Callias boundary data `(p,q,epsilon)=(1,2,+1)` give rank `1` and determinant c1 `2`.
- B=+1 / B=-1 fixed-polarization kernels: `{'positive': 3, 'negative': 0}` / `{'positive': 0, 'negative': 3}`; naive O(-2) control: `{'positive': 0, 'negative': 1}`.
- Factorized WZW pullback integer `n B d = 2`.

## Theorem-level conclusion

At low energy, a one-multiplet tangent N=2 SQM requires normalizable CP1
modes, a uniform gap, a Fredholm tangent kernel/connection, N=2 Ward
identities, and no extra zero modes.  Exact microscopic supercharges are one
sufficient way to enforce those Ward identities, but are not logically necessary
because accidental/emergent worldline supersymmetry is possible.  The conservative
same-model audit conditions below require either microscopic supercharges or an
independent derivation of all emergent Ward identities.  The current bosonic proxy contains none of the
required supersymmetry transformations or Yukawa mass family, so the result is
underdetermined rather than proved.

The degree-five differential character can always be fiber-integrated along
the spatial three-cycle on the framed configuration space.  Descent to the
gauge quotient additionally requires an equivariant differential refinement,
basic curvature, and trivial vertical/stabilizer holonomy.  Those data are not
present.  The exact pullback criterion is the integer integral over
Sigma3 x CP1; the factorized conditional map gives +2, but the present proxy
has not identified its scalar S2 with the same soliton CP1.

## Checks

- [PASS] `S0_provenance` — all same-soliton audit sources exist and are hashed: hashed=3/3
- [PASS] `S1_tangent_geometry` — Fubini--Study kinetic metric is chart invariant: max residual=4.441e-15
- [PASS] `S1_tangent_geometry` — psi^v=(-w^-2)psi^w has the tangent-bundle norm: max residual=2.665e-15
- [PASS] `S1_tangent_geometry` — Levi-Civita connection obeys the CP1 chart law: max residual=7.569e-14
- [PASS] `S1_tangent_geometry` — discrete Berry regression gives c1(TCP1)=+2: c1=1.999999999999981; triangles=2208
- [PASS] `S1_tangent_geometry` — Berry-Chern sequence returns k for k=0,1,2,4: sequence=[0.0, 0.999999999999988, 1.9999999999999807, 3.999999999999958]
- [PASS] `S2_callias_family` — finite H=1-|u_2><u_2| template/control has an exact one-dimensional kernel: kernel=4.561e-16; projector=4.564e-16
- [PASS] `S2_callias_family` — finite template/control has a uniform unit nonzero gap: gap=0.999999999999999
- [PASS] `S2_callias_family` — finite template/control kernel Berry line is O(+2): c1=1.999999999999998
- [PASS] `S2_callias_family` — conditional positive-mass eigenbundle (p,q)=(1,2) gives rank one and c1=+2: result={'rank': 1, 'determinant_c1': 2}
- [PASS] `S2_callias_family` — q=0 negative control retains pointwise index but has trivial determinant line: result={'rank': 1, 'determinant_c1': 0}
- [PASS] `S2_callias_family` — CPT-doubled q=+2 and q=-2 control cancels determinant c1: c1 sum=0
- [PASS] `S3_sqm_algebra` — finite B=+1 and B=-1 blocks satisfy Q^2=0 and {Gamma,D}=0: max residual=0.000e+00
- [PASS] `S3_sqm_algebra` — O(2) has three positive-chirality Dolbeault ground states: (h0,h1,index)=(3,0,3)
- [PASS] `S3_sqm_algebra` — fixed-polarization CPT/Serre duality requires K tensor E^vee=O(-4): degree=-4; (h0,h1,index)=(0,3,-3)
- [PASS] `S3_sqm_algebra` — Serre-dual massive spectra and gaps agree between B=+1 and B=-1: gap+=2.412090756622; gap-=2.412090756622
- [PASS] `S3_sqm_algebra` — naive dual O(-2) is a one-state negative control, not the CPT triplet: (h0,h1)=(0,1)
- [PASS] `S3_sqm_algebra` — the conditional tangent line and Callias determinant have matching c1=2: c1(TCP1)=c1(det Ind D)=2 in the declared conditional template
- [PASS] `S4_differential_transgression` — fiber integration along Sigma3 lowers differential degree 5 to degree 2: 5-3=2
- [PASS] `S4_differential_transgression` — factorized same-sphere map gives c1=n B d=+2 for (n,B,d)=(2,1,1): integer=2
- [PASS] `S4_differential_transgression` — orientation reversal B=+1 to B=-1 reverses the raw WZW line: integer=-2
- [PASS] `S4_differential_transgression` — degree-zero orientation-map control pulls the WZW line back trivially: integer=0
- [PASS] `S4_differential_transgression` — level-one control produces c1=1 and therefore cannot supply O(2): integer=1
- [PASS] `S4_differential_transgression` — descent audit enumerates local and global necessary conditions: conditions=['equivariant_differential_refinement', 'curvature_horizontal_and_invariant', 'vertical_and_large_gauge_holonomy_trivial', 'stabilizer_character_trivial', 'class_in_image_of_pullback_from_quotient']
- [PASS] `S4_differential_transgression` — gauge descent fails closed when no condition is supplied by the mother model: all same-model descent evidence flags are false
- [PASS] `S5_underdetermination` — same-model audit records all eight conservative closure conditions: conditions=8
- [PASS] `S5_underdetermination` — present charged two-colour proxy supplies none of the CP1/SUSY/Callias closure conditions: bosonic potential and radial Hessian do not determine a Yukawa family

## Fail-closed composition gate

- Specify one charged-two-colour mother action and derive either microscopic supercharges or complete emergent N=2 Ward identities, together with its Yukawa operator.
- Prove a normalizable CP1 family on the existing B=1 soliton and a uniform noncollective gap.
- Compute the actual Callias asymptotic positive-mass eigenbundle; do not substitute the finite matrix template/control.
- Derive the anti-canonical vacuum-line factor K=O(-2) in the B=-1 regulator, rather than using the naive dual O(-2).
- Construct an equivariant degree-five differential character and prove trivial large-gauge/stabilizer holonomy.
- Evaluate the same-soliton integral and obtain +2 before composing AP-E3 with AP-E4.
- Do not build the degree-one Route-E portal until this gate is closed.
