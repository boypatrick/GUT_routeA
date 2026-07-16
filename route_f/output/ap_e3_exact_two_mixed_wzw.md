# AP-E3 exactly-two and mixed-WZW audit

- Status: `ap_e3_exact_cell_and_anomaly_consistent_mixed_wzw_candidate_physics_portal_open`
- Checks: `26/26`
- Exactly-two constrained-cell theorem: `true`
- Physical Pauli orientation gives k=+2: `true`
- Dynamical gauge anomalies cancel: `true`
- Mixed-WZW intermediate UV candidate: `true`
- Full all-scale/nonperturbative UV closure: `false`
- Route-E portal closed: `false`
- Physics promotion allowed: `false`

## Derived chain

```text
G_r=N_r-1=0 (r=1,2) -> H_phys=C2 tensor C2
-J_H S1.S2 + h n.(S1+S2), J_H,h>0 -> unique m=-1 triplet
Pauli electron sign -> anti-aligned line Q tensor Q = O(2)
i hbar <Omega_-|d Omega_-> = +2 hbar A_+ -> k=+2
SU(2)c x U(1)g, nc=2, Xq=1 -> kappa_L=-kappa_R=2
S_mix=2 pi hbar 2 integral(omega3 wedge omega2)
B=+1 worldline -> k=nB=+2
```

Numerical Hund-Zeeman spectrum: `[-0.45, -0.25, -0.04999999999999999, 0.75]`.
Numerical Chern sequence: `[2.0082484079079745, 2.0020576482854167, 2.000514134394606, 2.000128516254437, 2.0000321279797606]`.

## Remaining blockers

- Prove nonperturbatively that the charged two-colour theory condenses in the mesonic SU(2)_L x SU(2)_R -> SU(2)_V channel and gaps every Pauli-Guersey/diquark direction.
- Keep confinement at or above the Higgs-vector scale, or otherwise show that approximate Pauli-Guersey restoration cannot unwind the B=1 soliton.
- Audit compact-U(1) monopoles, the discrete axial anomaly, the spin-bordism sign, and an explicit all-scale embedding beyond the abelian Landau pole.
- Construct the differential-character/Cech-bordism definition on non-extendible four-dimensional sectors; the present test proves only extension independence when an extension exists.
- Compute the soliton Hessian, Finkelstein-Rubinstein statistics, collective-coordinate inertia, and the LLL/adiabatic gap before retaining only the three O(2) states.
- Derive a degree +1 Route-E portal map; neither the constrained cell nor the mixed-WZW sector yet identifies its triplet with three chiral Spin(10) families.
- The complete-cell boundary theorem covers positive parent interactions only; arbitrary intercell couplings may generate emergent projective edge modes even though literal odd-occupancy states are absent.

## Mechanical checks

- [PASS] `U0_provenance` - all exact-two/mixed-WZW critical sources exist and are hashed: hashed=4/4
- [PASS] `U1_exact_cell` - two compact Gauss laws implement the exact Fourier projector P_1 P_1: basis=36; projector residual=3.538e-17
- [PASS] `U1_exact_cell` - the constrained physical cell has rank four and exactly total occupation two: rank=4
- [PASS] `U1_exact_cell` - separate Gauss laws remove literal odd, singleton, and charge-transfer sectors: forbidden=(odd totals),(2,0),(0,2)
- [PASS] `U1_exact_cell` - negative control in the N_r<=2 bosonic audit: one total-number Gauss law is strictly weaker: rank(P_total=2)=10; transfer survives=True
- [PASS] `U1_exact_cell` - negative control in the N_r<=2 bosonic audit: even parity neither selects N=2 nor removes charge transfer: rank(P_even)=20
- [PASS] `U1_exact_cell` - integer background charges pass 0+1d large-gauge invariance: integer phase=(1+2.4492935982947064e-16j); half-integer phase=(-1-1.2246467991473532e-16j)
- [PASS] `U2_orientation` - Hund plus Pauli-Zeeman Hamiltonian has the exact triplet/singlet spectrum: spectrum=[-0.45, -0.25, -0.04999999999999999, 0.75]
- [PASS] `U2_orientation` - the unique local gap is min(h,J_H+h)=h for J_H,h>0: gap=0.200000000000
- [PASS] `U2_orientation` - random-direction diagonalization selects the anti-aligned tensor-square ground line: projector=1.110e-15; eigen=2.720e-16
- [PASS] `U2_orientation` - the microscopic kinetic term gives signed action connection +2 A_+: max residual=1.624e-16
- [PASS] `U2_orientation` - the anti-aligned pair quantum metric is twice the one-doublet metric: max residual=8.327e-17
- [PASS] `U2_orientation` - pair curvature quadrature converges to c1=+2 with second-order refinement: sequence=[2.0082484079079745, 2.0020576482854167, 2.000514134394606, 2.000128516254437, 2.0000321279797606]
- [PASS] `U2_orientation` - the quotient-line pair transition has winding +2: winding=2.000000000000
- [PASS] `U3_anomaly` - all local dynamical gauge-anomaly sums vanish: U1^3=0; grav^2-U1=0; SU2^2-U1=0.0; SU2-U1^2=0.0
- [PASS] `U3_anomaly` - the SU(2)_c Witten anomaly cancels because the Weyl-doublet count is even: left Weyl doublets=4
- [PASS] `U3_anomaly` - the mixed flavour/one-form Postnikov classes are kappa_L=-kappa_R=2: kappa_L=2; kappa_R=-2
- [PASS] `U3_anomaly` - both SU(2)_L and SU(2)_R background Witten anomalies have even colour multiplicity: copies per flavour SU(2)=2
- [PASS] `U3_anomaly` - gauged U(1)_g exactly removes off-diagonal Pauli-Guersey generators: off-diagonal generators change U(1)_g charge by 2 and do not commute with the gauge generator
- [PASS] `U4_WZW` - the five-dimensional mixed WZW exponent is extension-independent for tested extendible integral classes: tested=200; n=2
- [PASS] `U4_WZW` - a positive unit baryon reduces the mixed WZW term to signed worldline k=+2: k(B=+1)=2; k(B=-1)=-2
- [PASS] `U4_WZW` - the minimal local colour singlet requires two quarks and two conjugate Higgs dressings: qq charge=2; dressing charge=-2; total=0
- [PASS] `U4_WZW` - the local Higgs dressing is Sym^2(conjugate doublet), hence a triplet of degree two: dim Sym^2(C^2)=3; associated projective line is O(2)
- [PASS] `U4_WZW` - a coloured singleton is excluded from the local gauge-invariant operator algebra: one SU(2)_c fundamental is not a colour singlet; the first local baryon uses epsilon_ab q^a q^b
- [PASS] `U5_boundaries` - the running diagnostic detects asymptotically free colour but an abelian Landau pole: b_U1=6.000000; b_SU2=6.000000
- [PASS] `U5_boundaries` - pure two-colour QCD is retained as a negative control, not used for the completion: without U(1)_g, SU(4)->Sp(4) gives S^5 and pi_3(S^5)=0; the charged theory needs a separate phase audit

## Source hashes

- `route_f/tex/ap_e3_exact_two_mixed_wzw_uv.tex` - `6f1a2490ddfd909c3f5a5b72c0466cd93b4097b0aa6e852dd715923697faa781` (32225 bytes)
- `route_f/tex/ap_e3_exact_two_mixed_wzw_uv.bib` - `707d2671e0fcee6074c6d96fef5101445b5e702b9c4457ace5161a7ff55b49cf` (2631 bytes)
- `route_f/LITERATURE_2023_2026.md` - `81bbc6917cf9a62cc0944bf999823a9f0dbade0534412b0a63818eea00825d63` (29765 bytes)
- `route_f/code/verify_ap_e3_exact_two_mixed_wzw.py` - `0e9844ab5f00526f8a828d8b9f3bc64330ac1a330b4b5248b6b4b2cc102d6a8c` (22905 bytes)
