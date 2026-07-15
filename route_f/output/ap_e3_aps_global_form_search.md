# AP-E3 APS generators and semisimple-global-form search

- Status: `pass_with_open_aps_threshold_and_dynamics_gates`
- Mechanical checks: `58/58`
- Best simple candidate in scanned set: `Sp(4)=USp(4)=Spin(5), simply connected`
- Minimal product backup: `SU(2)c x SU(2)H, no diagonal quotient`
- Torsion pair selected by charged-QC2D UV: `false`
- Exactly-two operator realized at group-theory level: `true`
- Full UV closed: `false`
- Physics promotion allowed: `false`
- Degree-one portal allowed: `false`

## APS / Dai-Freed result

The ordinary product Dirac spectra on `G3=S1_R x S3` and
`G2=T2_odd x S2` are paired and give the plain-regulator phase `+1`.
Decorating the transverse inverse image by the mod-two circle index or the
Arf phase gives independent regulator stacks:

- `(epsilon_3,epsilon_2)=(0,0)` -> `(G3,G2)=(+1,+1)`
- `(epsilon_3,epsilon_2)=(0,1)` -> `(G3,G2)=(+1,-1)`
- `(epsilon_3,epsilon_2)=(1,0)` -> `(G3,G2)=(-1,+1)`
- `(epsilon_3,epsilon_2)=(1,1)` -> `(G3,G2)=(-1,-1)`

Thus the two generator phases are computable after a regulator is chosen,
but the present charged-QC2D Lagrangian does not choose that regulator.

## Candidate classification

- `SU(3)` / `simply connected`: `passes=false` - unbroken subgroup is (SU2 x U1)/Z2; 1_(+1) forbidden
- `SU(2)c x SU(2)H` / `direct product`: `passes=true` - (2,2) and (1,2) contain 2_(+1) and 1_(+1)
- `SO(4)` / `(SU2 x SU2)/Z2`: `passes=false` - (1,2) scalar does not descend
- `Sp(4)=Spin(5)` / `simply connected`: `passes=true` - 5 contains 2_(+1); 4 contains 1_(+1)
- `SO(5)` / `Sp4/Z2`: `passes=false` - spinor 4 scalar does not descend
- `SU(4)` / `simply connected`: `passes=true` - 4 and bar4 furnish equal-magnitude doublet/singlet charges
- `SU(4)/Z2` / `central quotient`: `passes=false` - 4 and bar4 do not descend

The selected `Sp(4)` cover has

```text
5  -> 2_(+1) + 2_(-1) + 1_0
4  -> 1_(+1) + 2_0 + 1_(-1)
10 -> 1_(+2) + 2_(+1) + 3_0 + 1_0 + 2_(-1) + 1_(-2)
b0(Sp4)=15/2 > 0
```

Its tree-level adjoint split is explicit, but radiative protection and the
heavy-threshold eta phase are not established.

## Checks

- [PASS] `provenance` - all declared sources exist: TeX, bibliography, and verifier are present and nonempty
- [PASS] `provenance` - source hashes are well formed: hashed=3/3; unique=3
- [PASS] `provenance` - claim-boundary anchors in TeX: anchors=5/5
- [PASS] `bordism` - stable splitting: Omega_4^Spin(S3 x S2)=Z + Z2 + Z2
- [PASS] `dirac_product` - periodic circle KO zero mode: ker D_S1 periodic=1 mod 2
- [PASS] `dirac_product` - antiperiodic circle negative control: ker D_S1 antiperiodic=0
- [PASS] `aps` - periodic circle mod-two spectral flow: crossings=1
- [PASS] `aps` - antiperiodic circle spectral-flow control: crossings=0
- [PASS] `dirac_product` - odd torus Arf index: dim_C ker D_T2^+(RR)=1
- [PASS] `dirac_product` - three even torus spin structures: zero counts={'RR': 1, 'RA': 0, 'AR': 0, 'AA': 0}
- [PASS] `dirac_product` - ambient S1 x S3 product gap: minimum |lambda|=1.500000000000
- [PASS] `dirac_product` - ambient T2 x S2 product gap: minimum |lambda|=1.000000000000
- [PASS] `aps` - plain product eta symmetry: truncated sign-sum residual=0.000e+00
- [PASS] `aps` - four independent torsion regulators: all Hom(Z2+Z2,U(1)) characters realized
- [PASS] `aps` - epsilon3 spectator toggle: stacking the codimension-three defect regulator flips only G3
- [PASS] `aps` - epsilon2 spectator toggle: stacking the Arf defect regulator flips only G2
- [PASS] `aps` - local curvature blind to torsion stack: the four choices share the same degree-five differential curvature
- [PASS] `global_form` - diagonal quotient admits odd charged doublet: (-1)^(2j+q)=+1 for 2_(+1)
- [PASS] `global_form` - diagonal quotient excludes charge-one singlet: (-1)^(2j+q)=-1 for 1_(+1)
- [PASS] `global_form` - unquotiented form admits both required fields: H=SU(2)c x U(1)H has no diagonal-centre constraint
- [PASS] `operator` - exactly-two neutrality: Q[qq(phi^dagger)^2]=2 Xq-2 Xphi=0
- [PASS] `candidate_su3` - SU3 obstruction reproduced: 3 -> 2_(+1)+1_(-2); every singlet has even charge
- [PASS] `candidate_product` - SU2 x SU2 branching: (2,2)->2_(+1)+2_(-1), (1,2)->1_(+1)+1_(-1)
- [PASS] `candidate_product` - direct product centre is faithful: no nontrivial central pair is identified with the identity
- [PASS] `candidate_product` - SO4 quotient negative control: the diagonal quotient removes the required (1,2) scalar
- [PASS] `candidate_sp4` - Sp4 fundamental branching: 4 -> [(1, 1), (2, 0), (1, -1)]
- [PASS] `candidate_sp4` - Sp4 vector branching from wedge square: 5 -> [(2, 1), (1, 0), (2, -1)]
- [PASS] `candidate_sp4` - Sp4 adjoint branching from symmetric square: 10 -> [(1, 2), (2, 1), (3, 0), (1, 0), (2, -1), (1, -2)]
- [PASS] `candidate_sp4` - Spin5-cover embedding has trivial kernel: (-1c,-1H) maps to the nontrivial Spin5 centre, not to identity
- [PASS] `candidate_sp4` - SO5 quotient excludes spinor scalar: 4 is a Spin5 spinor and does not descend to SO5
- [PASS] `candidate_su4` - SU4 fundamental branches contain both charges: 4 -> 2_(+1)+2 copies of 1_(-1); bar4 is conjugate
- [PASS] `candidate_su4` - SU4 block embedding kernel: minimum nontrivial central-image distance=2.828427
- [PASS] `candidate_su4` - fundamental SU(N) equal-charge uniqueness: 2a+(N-2)b=0 and a=-b!=0 imply N=4
- [PASS] `candidate_su4` - SU4 quotient negative control: quotient by {-I4} forbids 4 and bar4, so the candidate needs SU4 itself
- [PASS] `threshold` - product fermion split: M(+1)=0.0, M(-1)=-6.0
- [PASS] `threshold` - product scalar split: m2(+1)=-1.0, m2(-1)=3.0
- [PASS] `threshold` - Sp4 fermion partner split: masses={1: 0.0, -1: -6.0, 0: -3.0}
- [PASS] `threshold` - Sp4 scalar partner split: mass-squared={1: -1.0, -1: 3.0, 0: 1.0}
- [PASS] `threshold` - SU4 adjoint split: the first adjoint separates upper and lower two-dimensional blocks
- [PASS] `anomaly` - product deep Witten parities: left-Weyl doublets: colour=8, H=8
- [PASS] `anomaly` - low colour Witten parity: light left-Weyl colour doublets=4
- [PASS] `anomaly` - Sp4 vectorlike local anomaly: two Dirac 5s have cancelling left/right anomaly phases
- [PASS] `anomaly` - Sp4 pi4 parity for displayed Dirac matter: chiral copies of each real 5=4
- [PASS] `anomaly` - SU4 vectorlike cubic anomaly: A(4)+A(bar4)=0 per Dirac flavour
- [PASS] `running` - product colour asymptotic freedom: b0_c=14/3
- [PASS] `running` - product H asymptotic freedom: b0_H=4
- [PASS] `running` - Sp4 asymptotic freedom: b0_Sp4=15/2
- [PASS] `running` - SU4 asymptotic freedom: b0_SU4=35/3
- [PASS] `orientation` - light-field level magnitude and sign: Nc Xq=2*(+1)=+2 before heavy-threshold counterterms
- [PASS] `claim_boundary` - charged-QC2D epsilon3 remains open: the codimension-three spectator toggles the same local data
- [PASS] `claim_boundary` - charged-QC2D epsilon2 remains open: the Arf spectator toggles the same local data
- [PASS] `claim_boundary` - torsion pair is not inferred from curvature: four regulator choices remain until the microscopic determinant is fixed
- [PASS] `claim_boundary` - Sp4 full gauge-bordism audit remains open: displayed vectorlike/Witten checks do not replace Omega_5^Spin(BSp4) evaluation with all backgrounds
- [PASS] `claim_boundary` - radiative threshold protection remains open: M(+1)=0 uses m=-yV and must be symmetry-protected or retuned
- [PASS] `claim_boundary` - heavy eta threshold remains open: the sign of the light +2 coefficient is not yet an all-scale Dai-Freed match
- [PASS] `claim_boundary` - nonperturbative Sp4 route remains open: no Sp4 strong-phase or B=1 soliton computation is supplied here
- [PASS] `claim_boundary` - unit-cell rule is not inferred from operator neutrality: qq(phi^dagger)^2 is neutral/triplet, but the microscopic exactly-two Gauss rule is not rederived
- [PASS] `claim_boundary` - no premature Route-E promotion: group theory closes first; portal waits for APS, thresholds, and dynamics

## Remaining blockers

- Specify the actual charged-QC2D Euclidean Dirac/Yukawa operator and Pauli-Villars mass signs, then evaluate both generator determinant phases; the product defect models prove possible values, not microscopic selection.
- For Sp(4), compute the complete Dai-Freed gauge-bordism phase with every retained fermion and Higgs background; vectorlike and pi4 parity checks are necessary but deliberately not called a full audit.
- Protect or retune the Sp(4) adjoint split M(+1)=0, and calculate heavy-fermion threshold counterterms fixing whether the infrared +2 orientation survives all-scale matching.
- Demonstrate the desired charged two-colour phase and a stable B=1 soliton in the Sp(4)-broken theory, including monopole-induced violation of the emergent magnetic one-form symmetry.
- Only after APS, threshold, faithful-bundle, and nonperturbative gates close may the degree-one Route-E portal be constructed.

## Source manifest

- `route_f/tex/ap_e3_aps_global_form_search.tex` - `8fa56b8473ee6bdc93ae809da3f33ab0c756acf63f93c37df16f9288d04f954d` (31713 bytes)
- `route_f/tex/ap_e3_aps_global_form_search.bib` - `2f274a4c4fcd6447d90d28a664eae58af317eb20ae976bdf56ec2e0399323cf0` (3444 bytes)
- `route_f/code/verify_ap_e3_aps_global_form_search.py` - `b99502b2b797502604abb4ed22f78dc0d9b16889317ee459209e051e97c36b17` (33702 bytes)
