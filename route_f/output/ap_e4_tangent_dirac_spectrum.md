# AP-E4 tangent-valued Dirac audit

- Status: `ap_e4_canonical_spinc_mathematics_complete_physical_origin_and_anomalies_open`
- Checks: `22/22`
- Tangent bundle derived: `true`
- Canonical Spin-c mathematics complete: `true`
- Three chiral zero modes: `true`
- Full partner spectrum and gap: `true`
- Physical fermionic tangent mode derived: `false`
- AP-E4 physics closed: `false`
- Physics promotion allowed: `false`

## Canonical Spin-c result

```text
E=T^(1,0)CP1=O(2)
D_T^c=sqrt(2)(dbar_T+dbar_T^dagger)
ker D_+=H0(O(2))=C^3, ker D_-=H1(O(2))=0
lambda_(n,+/-)=+/- sqrt(n(n+3))/R
mult per sign=2n+3, n>=1
R=1/2 -> gap=4; first massive multiplicity=5 per sign
```

Ordinary spin Dirac twisted only by `T=O(2)` instead gives two zero modes and
gap `2 sqrt(3)`.  Canonical Spin-c with `O(2)` is equivalent to ordinary spin
twisted by `O(3)`; the half-canonical shift is physical input, not notation.

## Remaining blockers

- The horizontal fluctuation Q_z delta z is a commuting sigma-model scalar; it does not derive a fermionic Clifford module or a Dirac operator on target CP1.
- Choose and derive either N=2 moduli-space/SQM quantization or an explicit higher-dimensional product containing CP1.  If CP1 is only field-value space, the target spectrum is not a Kaluza-Klein particle spectrum.
- Explain microscopically why the canonical Spin-c structure is selected.  Relative to the unique ordinary spin structure it supplies a half-canonical O(1) shift; ordinary spin twisted only by T has two zero modes.
- Prove that H0(TCP1) represents matter families rather than automorphism, gauge, or collective-coordinate redundancies.
- For a six-dimensional Weyl realization, cancel the full I8 gauge/gravitational anomaly (or exhibit Green-Schwarz factorization) and audit the Spin-c determinant flux.
- Show that all four-dimensional representation-weighted anomalies cancel and exclude conjugate bundles or additional sectors with compensating zero modes.
- Construct the Route-E portal below the gap 4 v/R_normalization and show its orientation degree is +1 without closing the gap or mixing in the five-fold first massive level.

## Mechanical checks

- [PASS] `E40_provenance` - all AP-E4 critical sources exist and are hashed: hashed=5/5
- [PASS] `E41_geometry` - the AP-E1 Fubini-Study convention is a round sphere of radius one half: R=0.5; K=4.0; Scal=8.0
- [PASS] `E41_geometry` - the tangent transition has degree two, hence T^(1,0)CP1=O(2): transition winding=2.000000000000
- [PASS] `E41_geometry` - the tangent Chern connection has c1=2: flux sequence=[2.0082484079079745, 2.0020576482854167, 2.000514134394606, 2.000128516254437, 2.0000321279797606]
- [PASS] `E41_geometry` - Q_z delta z is the horizontal tangent-valued sigma-model fluctuation: horizontal=1.097e-17; norm=2.082e-17
- [PASS] `E42_index` - canonical Dolbeault Spin-c twisted by T=O(2) has (h0,h1,index)=(3,0,3): h0=3; h1=0; index=3
- [PASS] `E42_index` - Serre duality independently gives H1(O(2))=H0(O(-4))*=0: h1(O2)=0; h0(O-4)=0
- [PASS] `E42_index` - Riemann-Roch/Todd integration gives chi(O(q))=q+1: integral ch(O(2)) Td(T)=1+2=3
- [PASS] `E42_index` - negative control: ordinary spin Dirac twisted only by T=O(2) has two, not three, zero modes: effective Dolbeault line=O(1); h0=2; index=2
- [PASS] `E42_index` - canonical Spin-c O(2) equals ordinary spin twist O(3), exposing the half-canonical shift: S^+=O(-1); S^+ tensor O(3)=O(2)
- [PASS] `E42_index` - nearby bundle controls q=1 and q=3 give two and four canonical zero modes: q=1 -> 2; q=3 -> 4
- [PASS] `E43_spectrum` - the monopole-Casimir identity yields lambda_n^2=n(n+q+1)/R^2 exactly: levels checked=10
- [PASS] `E43_spectrum` - independent finite SU(2) matrices reproduce every tested D^2 level and degeneracy: Casimir=1.421e-14; shifted spectrum=1.421e-14
- [PASS] `E43_spectrum` - the truncated spectral operator is Hermitian and anticommutes with chirality: Hermitian=0.000e+00; anticommutator=0.000e+00
- [PASS] `E43_spectrum` - the kernel contains exactly three positive-chirality states: zero count=3; kernel chirality trace=3.0
- [PASS] `E43_spectrum` - every nonzero level is paired at plus/minus lambda as implied by chirality anticommutation: paired nonzero states per sign=45
- [PASS] `E43_spectrum` - at R=1/2 the canonical first gap is four with multiplicity five per sign: gap=4.000000000000; multiplicity per sign=5
- [PASS] `E43_spectrum` - ordinary spin twist T=O(2) has the distinct gap 2 sqrt(3) and first multiplicity four: ordinary gap=3.464101615138; multiplicity=4
- [PASS] `E43_spectrum` - radius drift control: changing R from one half to one halves the canonical gap: gap(R=1)=2.000000; gap(R=1/2)=4.000000
- [PASS] `E44_fail_closed` - a chirality-breaking mass lifts the protected zero modes: mass=0.3; remaining zeros=0
- [PASS] `E44_fail_closed` - nonintegral monopole degree is rejected as a global line bundle: q=5/2 rejected; q=2 accepted
- [PASS] `E44_fail_closed` - absence of opposite chirality applies only to the kernel, not the massive tower: kernel=3 positive; massive states=90

## Source hashes

- `route_f/tex/ap_e4_tangent_dirac_spectrum.tex` - `8a212bc5db3fb0c1c176a89ed6c77943f914972c5a6f78d4a2f31db737d63307` (22766 bytes)
- `route_f/tex/ap_e4_tangent_dirac_spectrum.bib` - `63ef8f67f0b30b4465d4e0a4fa911e163b3ee6358a373effa421bb6de97e7774` (1254 bytes)
- `route_f/tex/ap_e1_projective_doublet_action.tex` - `cd4e36772688a48f97c09c9e1251b4b61256d2441ff9b7247e3a7d27d0014ca0` (37998 bytes)
- `route_f/LITERATURE_2023_2026.md` - `71f24c2d5725b32d50b5ff8073919b86b4eb69c2092d9bc52929b97034956a41` (15797 bytes)
- `route_f/code/verify_ap_e4_tangent_dirac_spectrum.py` - `ceaee9ca07bd1eaeff4a10e4092917c17d2b4b24ff15252bb95e05a2d288dc1c` (20996 bytes)
