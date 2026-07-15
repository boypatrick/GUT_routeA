# AP-E3 microscopic origin of projective level magnitude two

- Status: `ap_e3_hund_pair_derives_abs_k2_conditionally_uv_and_sign_open`
- Checks: `27/27`
- Candidate-level magnitude derivation complete: `true`
- Level magnitude: `2`
- Signed level for the declared aligned action: `-2`
- Positive signed Route-E level selected: `false`
- AP-E3 physics closure: `false`
- Physics promotion allowed: `false`

## Conditional first-principles result

For two orbitals with one doublet per orbital and ferromagnetic Hund locking,

```text
H_micro = sum_r [U N_r(N_r-1)/2 - mu N_r]
          - J_H S_1.S_2 - h n.(S_1+S_2),
C = mu + h/2 + J_H/4 < U - J_H/8.
```

the low-energy state is the symmetric pair

```text
|s;2> = |s> tensor |s>,
-i <s;2|d|s;2> = 2 A,
L_ket = nu_2^* O_CP2(-1) = O_CP1(-2),
L_pre = dual(L_ket) = O_CP1(2).
```

Thus the declared dimer derives `|k_eff|=2`, a three-state spin-one carrier,
and a singlet gap `Delta_singlet=J_H`.  Relative to `S_k=hbar*k*integral(A)`,
the displayed aligned microscopic action gives `k=-2`; reversing the
orientation gives `k=+2`.  This is not a derivation of the two-orbital UV field
content or of the signed Route-E chirality.

Numerical Chern sequence: `[2.0020576482854167, 2.000128516254437, 2.000008031927202, 2.000000501994128]`.
Hund spectrum at `J_H=1`: `[-0.25, -0.25, -0.25, 0.75]`.
Interacting-plateau sufficient-condition margin: `2.025`.

## Alternative branches

- **fermion_determinant:** r filled same-orientation eigenlines induce k=sum C_a; r=2 can give |k|=2, but flavor count, filling, and signs must be independently fixed
- **mixed_WZW:** a quantized mixed WZW coefficient can reduce to k=n_c B and is nonrenormalization-protected, but an n_c=2 UV completion with the correct global coset is not present in this repository
- **tangent_Dirac:** T_CP1=O(2) gives a three-section Spin-c/Dolbeault index only after the physical mode is proved tangent-valued; this is deferred to AP-E4

## Remaining blockers

- derive why the UV field content contains exactly two orbitals/constituents
- prove odd or single-constituent sectors are absent rather than merely heavy by assumption
- embed the dimer and its Mott/Hund interactions into an anomaly-consistent four-dimensional field theory
- derive the portal to Route-E families and keep it below charge, singlet, and orientation gaps
- show the spin-one triplet is chiral family data rather than an internal vector multiplet
- derive the physical orientation/chirality that selects the required signed Route-E level rather than only |k|=2

## Mechanical checks

- [PASS] `M0_provenance` - all AP-E3 input sources exist and are hashed: hashed=6/6
- [PASS] `M0_provenance` - AP-E2 exact regression is green but explicitly leaves the Berry level underived: AP-E2=30/30; promotion=False
- [PASS] `M1_Mott` - the bare Bose-Hubbard window 0<mu<U locks exactly one constituent on each orbital: U=4.0; mu=1.5; one-orbital ground=1; pair=(1, 1)
- [PASS] `M1_Mott` - the bare analytic charge gap min(mu,U-mu) matches the enumerated spectrum: analytic=1.500000000000; enumerated=1.500000000000
- [PASS] `M1_Mott` - bare-window limits and the interacting false-inference counterexample are detected: mu<0 ground n=0; mu>U ground n=2; counterexample E11=-2.955000, E22=-5.860000
- [PASS] `M1_Mott` - a proved sufficient inequality selects the unique interacting (1,1) plateau: C=1.850000 < U-J_H/8=3.875000 (margin=2.025000); ground=(1, 1); full gap=1.850000000000; N>=8 lower bound=31.200000; ground energy=-3.450000
- [PASS] `M2_Hund` - the ferromagnetic Hund Hamiltonian is Hermitian and globally SU(2)-invariant: Hermitian residual=0.00e+00; SU2 residual=0.00e+00
- [PASS] `M2_Hund` - exact diagonalization gives a triply degenerate triplet and one singlet: eigenvalues=[-0.25, -0.25, -0.25, 0.75]
- [PASS] `M2_Hund` - the three low states have S^2=2 and the high state has S^2=0: S2 labels=[2.0, 1.9999999999999996, 2.0, 0.0]
- [PASS] `M2_Hund` - the spin-one projector equals the symmetric-exchange projector: max projector residual=0.00e+00; rank trace=3.0
- [PASS] `M2_Hund` - the singlet exclusion gap is exactly Delta_singlet=J_H: gap=1.000000000000; J_H=1.000000000000
- [PASS] `M2_Hund` - reversing the Hund sign selects a singlet and destroys the |k|=2 triplet: antiferromagnetic spectrum=[-0.75, 0.25, 0.25, 0.25]
- [PASS] `M2_Hund` - a weak orientation field gives a unique gapped ground-state line over every n in S2: 64 directions; max spectral residual=8.88e-16; min gap=0.200000000000
- [PASS] `M3_Veronese` - the symmetric tensor product is the normalized quadratic Veronese state: samples=256; max norm/map residual=9.99e-16; max singlet leakage=0.00e+00
- [PASS] `M3_Veronese` - the pair transforms with phase weight two: samples=256; max residual=3.51e-16
- [PASS] `M3_Veronese` - the pair connection doubles while the declared microscopic kinetic action has signed k=-2: samples=256; connection residual=2.24e-16; action-sign residual=2.24e-16
- [PASS] `M3_Veronese` - the Veronese pullback quantum metric is twice the single-doublet metric: samples=256; max metric residual=2.22e-16
- [PASS] `M3_Veronese` - the ket transition has winding -2 and its dual prequantum line has Chern number +2: single winding=-1.000000000000; pair winding=-2.000000000000
- [PASS] `M3_Veronese` - direct curvature quadrature converges to dual-prequantum first Chern number two: cells=[20, 80, 320, 1280]; c1=[2.0020576482854167, 2.000128516254437, 2.000008031927202, 2.000000501994128]; final error=5.02e-07
- [PASS] `M4_representation` - the locked pair is the three-state spin-one representation: dimension=3; max Casimir residual=4.44e-16; commutator residual=2.22e-16
- [PASS] `M4_representation` - large-gauge invariance quantizes k but does not uniquely select signed k=+2: levels 1,2,3 all invariant; max phase residual=2.20e-15
- [PASS] `M5_failure_modes` - opposite Berry orientations cancel to k=0 rather than |k|=2: samples=128; max cancellation residual=1.53e-16
- [PASS] `M5_failure_modes` - a permitted single constituent would retain an unwanted |k|=1 doublet: single Hilbert dimension=2 at |k|=1; pair dimension=3 at |k|=2; singleton_exclusion_derived=False
- [PASS] `M5_failure_modes` - without Hund locking the target is CP1xCP1 rather than its diagonal CP1: unlocked product Hilbert dimension=4; symmetric locked subspace dimension=3
- [PASS] `M6_EFT_gate` - an explicit benchmark satisfies portal, temperature, and drive scales below every retained gap: min gap=0.200000; max perturbation=0.030000; ratio=0.150000
- [PASS] `M6_EFT_gate` - a portal above the smallest gap is rejected by the scale-separation gate: bad portal=0.220000; min gap=0.200000; hierarchy pass=False
- [PASS] `M7_boundary` - the declared dimer derives |k|=2 while the signed Route-E orientation and UV selection remain open: magnitude-two proof complete; aligned action k=-2; signed orientation, UV field content, and Route-E identification remain open

## Source hashes

- `route_f/tex/ap_e1_projective_doublet_action.tex` - `cd4e36772688a48f97c09c9e1251b4b61256d2441ff9b7247e3a7d27d0014ca0` (37998 bytes)
- `route_f/output/ap_e2_discriminant_regression.json` - `c360753719e4fb3d1f75ccb3b04d8e7cfe89a37d173f6d3130989a7fd9bc25a8` (11893 bytes)
- `route_f/tex/ap_e3_level_two_microscopic_origin.tex` - `a649f72721cca74e3370ae8ee71bcd6da5808c97e493aefc01cf36d685d9ea94` (32635 bytes)
- `route_f/tex/ap_e3_level_two_microscopic.bib` - `101a75d92fcd026e8be1d266ac9e527e84768317fb170b9c8f8b76a4ae812a63` (3550 bytes)
- `route_f/LITERATURE_2023_2026.md` - `71f24c2d5725b32d50b5ff8073919b86b4eb69c2092d9bc52929b97034956a41` (15797 bytes)
- `route_f/code/verify_ap_e3_level_two_microscopic.py` - `8263b14e644d372150fd6db5abc1edda636eb391c69fbf4c5ff1b9a57bcb1c7a` (32248 bytes)
