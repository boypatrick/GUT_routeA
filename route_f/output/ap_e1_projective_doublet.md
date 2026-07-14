# AP-E1 projective-doublet audit

- Status: `ap_e1_conditionally_proved_with_explicit_o2_and_stability_blockers`
- Checks: `30/30`
- Physics promotion allowed: `false`

## Logical result

- **local_U1:** proved under one nonzero charge-one complex doublet, fixed norm, and an auxiliary local U(1): S3/U(1)=CP1 with the induced FS action
- **global_U1:** refuted as a pointwise CP1 quotient: the local Hopf-fiber phase remains physical
- **fixed_global_charge:** proved for the collective ansatz: Routh reduction gives a monopole-coupled T*CP1 rotor, not by itself H0(CP1,O(k))
- **holomorphic_sector:** conditional on first-order reduction or a controlled LLL/Kahler polarization: for k>=0, H_k=H0(CP1,O(k)), dim=k+1; for k<0 the conjugate anti-holomorphic sector has dim=|k|+1
- **route_E_O2:** not derived: k=2, the LLL gate, and the identification with T_CP1 and three chiral families remain separate assumptions

## Required branch decision

- **A_charge_two:** set the physical fixed charge to k=2; obtains the triplet but abandons the macroscopic Q=10^6 thin-wall benchmark
- **B_separate_level_two_internal_sector:** retain the macroscopic Q-ball charge but introduce an independent quantized level-two Wess--Zumino/Schwinger-boson sector
- **forbidden_shortcut:** identify Q=10^6 with k=2 or infer O(2) merely from CP1

## Non-derived inputs

- why the microscopic order parameter is exactly one complex doublet
- why the radial mode is frozen and what fixes its scale
- why the relevant prequantum level is k=2
- why the second-order rotor may be projected to its LLL
- why the resulting spin-one triplet is a family carrier
- how a stable four-dimensional bubble survives radial, gauge, and portal fluctuations

## Mechanical checks

- `PASS` **G1_Hopf_quotient / normalized doublet defines a rank-one Hermitian projector** — max|P^2-P|=5.56e-17, max|P^dag-P|=5.55e-17, TrP=1.0000000000000000, |detP|=5.80e-18
- `PASS` **G1_Hopf_quotient / projector is invariant under the common U(1) phase** — max|P(e^(i alpha)z)-P(z)|=5.57e-17
- `PASS` **G1_Hopf_quotient / Pauli Hopf vector lies on S2 and reconstructs the projector** — n^2=0.9999999999999996, max|P-(I+n.sigma)/2|=1.11e-16
- `PASS` **G2_induced_geometry / flat C2 metric splits into Hopf-fiber plus Fubini--Study metric** — |ds|^2=0.0183320320132069, A^2=0.0000379346403124, g_FS=0.0182940973728944, residual=3.47e-18
- `PASS` **G2_induced_geometry / Fubini--Study metric agrees on north and south charts** — north=0.0182940973728944, south=0.0182940973728944
- `PASS` **G2_induced_geometry / Hopf connection transforms as A to A+d alpha while horizontal norm is invariant** — A'-A-alpha_dot=0.00e+00, horizontal residual=3.47e-18
- `PASS` **G2_induced_geometry / Hopf curvature has first Chern number one** — integral F=6.283185468670607, integral F/(2pi)=1.000000025702
- `PASS` **G3_local_global_split / polar decomposition of the canonical doublet kinetic term is exact** — direct=1.1595327812435738, polar=1.1595327812435741
- `PASS` **G3_local_global_split / auxiliary local gauge field removes only the Hopf-fiber square** — global kinetic minus fiber=0.1057699414076649, local quotient=0.1057699414076649
- `PASS` **G4_fixed_charge / direct Routh transform equals the monopole-coupled CP1 expression** — direct=-0.2210913783783783, closed=-0.2210913783783784, residual=1.11e-16
- `PASS` **G4_fixed_charge / Legendre transform gives the minimally coupled monopole Hamiltonian** — direct=1.5639753783783781, closed=1.5639753783783783
- `PASS` **G4_fixed_charge / fixed-Q profile energy retains nonnegative orientation kinetic energy** — E-E_min=0.7355969999999999, (I/2)g*qdot^2=0.7355969999999999
- `PASS` **G5_holomorphic_quantization / C2 moment-map quotient is CP1 only at positive level; k=0 is a point** — k>0 regular-level gradient norm=2.828427; mu^-1(0) radius=0.0
- `PASS` **G5_holomorphic_quantization / corrected Branch-B symplectic sign imposes N=k and a positive level-k endpoint phase** — dL/d(hbar*c)=k-N=0.300000; delta L/(hbar*alpha_dot)=2.000000; reduced flux level=+2
- `PASS` **G5_holomorphic_quantization / Sym^0(C2) has k+1 states and the correct Cartan weights** — k=0, dimension=1, weights=[0]
- `PASS` **G5_holomorphic_quantization / Sym^1(C2) has k+1 states and the correct Cartan weights** — k=1, dimension=2, weights=[1, -1]
- `PASS` **G5_holomorphic_quantization / Sym^2(C2) has k+1 states and the correct Cartan weights** — k=2, dimension=3, weights=[2, 0, -2]
- `PASS` **G5_holomorphic_quantization / Sym^3(C2) has k+1 states and the correct Cartan weights** — k=3, dimension=4, weights=[3, 1, -1, -3]
- `PASS` **G5_holomorphic_quantization / Sym^4(C2) has k+1 states and the correct Cartan weights** — k=4, dimension=5, weights=[4, 2, 0, -2, -4]
- `PASS` **G5_holomorphic_quantization / Sym^5(C2) has k+1 states and the correct Cartan weights** — k=5, dimension=6, weights=[5, 3, 1, -1, -3, -5]
- `PASS` **G5_holomorphic_quantization / Sym^6(C2) has k+1 states and the correct Cartan weights** — k=6, dimension=7, weights=[6, 4, 2, 0, -2, -4, -6]
- `PASS` **G5_holomorphic_quantization / k=2 Schwinger sector is the spin-one triplet** — max commutator error=2.22e-16, max|J^2-2I|=4.44e-16
- `PASS` **G5_holomorphic_quantization / second-order k=2 rotor has a three-state LLL but infinitely many higher levels** — first four levels=[{"degeneracy": 3, "energy_times_I_over_hbar2": 2.0, "j": 1.0, "n": 0}, {"degeneracy": 5, "energy_times_I_over_hbar2": 10.0, "j": 2.0, "n": 1}, {"degeneracy": 7, "energy_times_I_over_hbar2": 22.0, "j": 3.0, "n": 2}, {"degeneracy": 9, "energy_times_I_over_hbar2": 38.0, "j": 4.0, "n": 3}], gap*I/hbar^2=8.0
- `PASS` **G5_holomorphic_quantization / monopole-rotor spectrum depends on |k| and is invariant under flux reversal** — k=+2 levels=[(1.0, 3, 2.0), (2.0, 5, 10.0), (3.0, 7, 22.0), (4.0, 9, 38.0)]; k=-2 levels=[(1.0, 3, 2.0), (2.0, 5, 10.0), (3.0, 7, 22.0), (4.0, 9, 38.0)]
- `PASS` **G6_route_e_gate / prior Q-ball angular frequency is Q divided by collective inertia** — I=1148821.582511051325, Q/I=0.870457184321, stored omega=0.870457184321
- `PASS` **G6_route_e_gate / macroscopic Q=10^6 sector is not the k=2 triplet** — if k=Q/hbar=1000000, dim H0(O(k))=1000001, not 3
- `PASS` **G7_Derrick / three-dimensional size dilation has E2~lambda, E4~lambda^-1, E0~lambda^3** — lambda=1.73, numerical ratios={"E0": 5.177716999999997, "E2": 1.729999999999994, "E4": 0.578034682080926}, expected={"E0": 5.1777169999999995, "E2": 1.73, "E4": 0.5780346820809249}, max relative error=3.44e-15
- `PASS` **G7_Derrick / positive E2 plus E0 cannot have a finite-size Derrick stationary point** — dE/dlambda at lambda=1 is E2+3E0=5.400000
- `PASS` **G7_Derrick / an E4 stabilizer can balance E2 at finite size** — lambda_star=sqrt(E4/E2)=2.000000000000, E''=3.000000000000
- `PASS` **G8_source_regression / TeX contains the corrected action signs, charge domain, and signed-spectrum formulas** — missing critical fragments=[] ; carriage_return_present=False
