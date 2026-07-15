# AP-E4 CP1 moduli-space N=2 SQM audit

- Passed: **27/27**
- Fail-closed physics promotion: **False**
- Intrinsic physical tangent fermion: **True**
- Mother-model canonical vacuum line: **False**
- CPT antibaryon polarization closed: **False**
- Charged-two-colour embedding: **False**
- AP-E4 worldline completion: **False**
- Mother-model L2 coefficient surrogate: `1.400009382609`
- Twisted k=2 Dirac gap: `3.380605690754`
- Twisted k=2 SQM gap: `5.714247418179`

## Main conclusion

The BPS vortex supplies CP1 orientational zero modes, a finite L2 metric,
and a tangent fermion.  Quantizing that tangent fermion supplies the
Clifford algebra, but the mother model does **not** select the canonical
Spin-c vacuum line used in the spectral card and does not supply an additional
T=O(2) coefficient bundle.  The declared canonical Spin-c
re-quantization has one untwisted ground state; the source-selected
half-form ordering has none.  Three positive-chirality ground states require a separately
derived E=O(2) flux/Fermi bundle.
The AP-E3 level-two WZW term can be that coefficient line only if
the same B=1 CP1 collective-coordinate embedding is proved and the
degree-five WZW character is first transgressed along the spatial
three-cycle, with the resulting line pulling back as O(+2).  The independent Chern--Simons vortex used here
does not establish that composability statement.

## Checks

- [PASS] `M0_provenance` — all moduli-SQM critical sources exist and are hashed: hashed=3/3
- [PASS] `M1_moduli_metric` — projected orientational variations are horizontal CP1 zero modes: max horizontal residual=3.817e-16
- [PASS] `M1_moduli_metric` — the L2 metric has the Fubini--Study shape |dw|^2/(1+|w|^2)^2: max metric residual=8.882e-16
- [PASS] `M1_moduli_metric` — the Fubini--Study kinetic norm agrees in the w and v=1/w charts: max chart residual=7.105e-15
- [PASS] `M1_moduli_metric` — the eliminated-rho near-BPS profile integral is positive, finite, and converged: I0 sequence=[0.22281841361727756, 0.2228184136172765, 0.22281841361727456, 0.22281841361726343]; C(e=1)=1.400009382609
- [PASS] `M1_moduli_metric` — the profile surrogate is not promoted to a solved microscopic vortex: profile_surrogate_is_microscopic_proof=false
- [PASS] `M2_tangent_fermion` — psi transforms as a holomorphic tangent vector under v=1/w: max D_t covariance residual=9.607e-14
- [PASS] `M2_tangent_fermion` — the tangent-fermion kinetic norm is chart invariant: max norm residual=6.217e-15
- [PASS] `M3_dolbeault` — Riemann--Roch and Serre duality give index(O(k))=k+1: k range=-6..6
- [PASS] `M3_dolbeault` — untwisted canonical Dolbeault SQM has one, not three, ground state: (h0,h1) for O(0)=(1,0)
- [PASS] `M3_dolbeault` — an independent O(2) coefficient line gives three positive-chirality zero modes: (h0,h1,index) for O(2)=(3,0,3)
- [PASS] `M3_dolbeault` — orientation reversal requires O(-4) for three negative-chirality zero modes: (h0,h1,index) for O(-4)=(0,3,-3)
- [PASS] `M3_dolbeault` — the de Rham completion has two CP1 vacua and is also not a three-state mechanism: b0+b2=2
- [PASS] `M3_flux` — the candidate magnetic line has integral flux c1(O(2))=2: curvature sequence=[2.0082484079079745, 2.0020576482854167, 2.000514134394606, 2.000128516254437, 2.0000321279797606]
- [PASS] `M3_flux` — the north--south transition of O(2) has winding two: winding=2.000000000000
- [PASS] `M3_flux` — canonical Spin-c twisting by O(2) has determinant line O(6): c1(det W_E)=2+2k=6
- [PASS] `M4_superalgebra` — the finite Dolbeault block has nilpotent Q and Q-dagger: Q^2=(Qdag)^2=0
- [PASS] `M4_superalgebra` — D is Hermitian, odd, and satisfies D^2={Q,Qdag}=2H: Hermitian/odd/superalgebra residuals below tolerance
- [PASS] `M4_spectrum` — the k=2 kernel contains exactly three positive-chirality zero modes: kernel=3; chiralities=[1.0, 1.0, 1.0]
- [PASS] `M4_spectrum` — the first twisted Dirac and SQM gaps scale as 4/sqrt(C) and 8/C: Delta_D=3.380605690754; Delta_H=5.714247418179; C=1.400009382609
- [PASS] `M4_spectrum` — the first massive level has five states at each Dirac sign: multiplicity(+gap)=multiplicity(-gap)=5
- [PASS] `M5_fail_closed` — a tangent fermion is not an automatic T-valued coefficient twist: ordinary_tangent_fermion_implies_E_equals_O(0), not E_equals_T
- [PASS] `M5_fail_closed` — the unmodified Chern--Simons vortex mother model has no derived O(2) Berry line: mother_model_level_two_flux_derived=false
- [PASS] `M5_composability` — the intrinsic vortex SQM is not silently identified with the charged-two-colour soliton: charged_two_colour_soliton_sqm_embedding_closed=false
- [PASS] `M5_composability` — AP-E3 k=2 can twist AP-E4 only after the same B=1 CP1 moduli map and spatially transgressed WZW line pullback are proved: same_B1_CP1_map=false; signed_WZW_pullback_to_O2=false
- [PASS] `M5_fail_closed` — the product compactification route is not used to obtain the SQM result: product_compactification_used=false
- [PASS] `M5_fail_closed` — the backup product compactification has not passed its six-dimensional anomaly gate: I8_factorization_or_cancellation=false; global_bordism_audit=false

## Fail-closed gates

- The analytic radial profile is a convergence surrogate, not a solved bulk vortex.
- The unmodified mother model has no induced level-two Berry line.
- The source-selected half-form ordering has no normalizable zero mode; the one-state untwisted count belongs to a declared canonical re-quantization.
- Its CP1 vortex is not the charged-two-colour AP-E3 soliton.
- AP-E3/AP-E4 composability needs the same B=1 moduli map and a spatial transgression whose signed line pulls back as O(2).
- The transgressed line must descend through the gauge quotient, and the B=-1 sector needs a derived CPT/anti-canonical polarization map.
- A bulk Callias/determinant-line or supersymmetric WZ derivation of O(2) remains open.
- The product-compactification route remains a backup; its six-dimensional local, global, and Green--Schwarz anomaly gate is open.
- No zero mode is identified with a chiral Route-E family.
