# AP-E7 common Sp(4) APS/PV regulator audit

- Status: mechanical_pass_physics_fail_closed
- Checks: 42/42 passed
- Common continuum quadratic APS/PV regulator: specified
- Physical fermion determinant phase: zero spectral cut with gapped endpoints
- Conditional Higgsed parity prediction: (+1,+1)
- Higgsed finite overlap Dirac-Yukawa/5D kernel: not defined; eta pair not computed
- Original deep-Sp(4) target eta pair: not selected
- Pure gauge/gravity heavy phase: exactly +1
- All-scale mixed-WZW k=+2: not matched
- Portal: not started

## Exact formulas

- \(I_h=\int[-c_2(V_c)+x^2-p_1(TM)/8]\) and \(N_f=2\) gives \(\exp(2\pi iI_h)=1\).
- Physical sign: \(\tau_F=\exp(i\pi N_f\,\mathrm{SF}_0)\); arbitrary \(\alpha\) is used only for the Fredholm index.
- \(J\widehat\Sigma J^{-1}=-\widehat\Sigma\), so fixed-adjoint-mass Kramers reality is not claimed.
- \(P_\pm=(h^2\pm h)/2,\ P_0=1-h^2\), with ranks \((2,2,1)\).
- Exact separated-projector obstruction: \(C(0)=cI_5\) has \(\mathrm{Tr}(HC(0))=0\), but \(\mathrm{Tr}(HP_+)=2\).
- Complete-irrep mixed trace: \(\kappa_{\rm full}=+2-2=0\).

## Gate ledger

| gate | value |
|---|---:|
| euclidean_gamma_reference_reality_data_specified | true |
| fixed_adjoint_mass_kramers_reality_proven | false |
| common_quadratic_aps_pv_regulator_specified | true |
| five_dimensional_aps_domain_specified | true |
| physical_fermion_phase_zero_cut_fixed | true |
| conditional_higgsed_gw_source_specified | true |
| finite_higgsed_overlap_dirac_yukawa_operator_defined | false |
| five_dimensional_higgsed_target_kernel_defined | false |
| minimal_higgsed_overlap_eta_pair_computed | false |
| conditional_higgsed_parity_predicts_plus_plus | true |
| pure_unbroken_gauge_gravity_heavy_phase_computed | true |
| pure_heavy_phase_on_G3_computed | true |
| pure_heavy_phase_on_G2_computed | true |
| exact_separated_rank_two_projector_obstruction_proven | true |
| full_nonperturbative_sp4_measure_constructed | false |
| original_action_target_mass_family_lift_defined | false |
| original_deep_sp4_action_uniquely_selects_target_pair | false |
| target_mapping_torus_eta_pair_selected | false |
| mixed_wzw_heavy_determinant_computed | false |
| k_plus_two_all_scale_matched | false |
| sp4_aps_pv_lane_closed | false |
| degree_one_route_e_portal_started | false |
| physics_promotion_authorized | false |

## Deterministic checks

| group | check | pass | detail |
|---|---|---:|---|
| provenance | all declared sources exist | yes | TeX, bibliography, and verifier are nonempty |
| provenance | source hashes are unique sha256 values | yes | hashed=3/3; unique=3 |
| provenance | all theorem anchors are present | yes | anchors=9/9 |
| provenance | physical fermion phase fixes the zero cut with gapped endpoints | yes | zero is excluded from both endpoint spectra and SF_0 defines the phase |
| provenance | arbitrary APS cut is restricted to Fredholm-index bookkeeping | yes | general alpha remains available for the APS index, not the physical sign |
| provenance | lateral Agmon cut and finite log are explicit | yes | upper lateral cut and renormalized finite part found |
| clifford | Hermitian Euclidean gamma matrices | yes | max residual=0.000e+00 |
| clifford | Euclidean Clifford algebra | yes | max residual=0.000e+00 |
| clifford | gamma5 involution and anticommutation | yes | max residual=0.000e+00 |
| reality | B gamma-star equals gamma B | yes | max residual=0.000e+00 |
| reality | quaternionic reality squares to minus one | yes | norm(BB*+1)=0.000e+00 |
| reality | adjoint Hermitian mass is odd under vector conjugation | yes | Hermitian residual=0.000e+00; K5 Sigmahat K5^-1 + Sigmahat=0.000e+00 |
| reality | anti-linear map sends the plus-Sigma mass to minus-Sigma | yes | family residual=0.000e+00; fixed-background residual=1.188e+00 |
| reality | finite scalar-reference gamma5-hermitian kernel | yes | K residual=0.000e+00; H residual=0.000e+00 |
| reality | scalar-reference two-flavour determinant square is nonnegative | yes | det(K)^2=0.725468539966 |
| quadratic_complex | ghost R-dagger R is symmetric positive | yes | min eigenvalue=2.826726e+00 |
| quadratic_complex | gauge-fixed boson toy is symmetric positive | yes | min eigenvalue=4.288672e-01 |
| pv | first three PV moments cancel | yes | moments={0: 0, 1: 0, 2: 0, 3: -6} |
| pv | large-lambda log coefficient is minus two | yes | lambda^3 sum(c log(1+a/lambda))=-1.999100213 |
| pv | two-flavour fermion amplitude exponents are integral | yes | Nf*c/2=[1, -3, 3, -1] |
| aps | gapped scalar endpoints give the physical zero-cut sign | yes | SF_0=-1; one-flavour sign=-1 |
| aps | auxiliary cut changes raw spectral flow and is not the physical sign | yes | SF_0=-1; SF_2=0 |
| aps | periodic circle path has one crossing | yes | crossings=1 |
| aps | antiperiodic circle path has no crossing | yes | crossings=0 |
| target | four abstract target torsion characters | yes | characters={(0, 0): (1, 1), (0, 1): (1, -1), (1, 0): (-1, 1), (1, 1): (-1, -1)} |
| projectors | P plus P minus P zero are orthogonal projectors | yes | idempotent=0.000e+00; orthogonal=0.000e+00; sum=0.000e+00 |
| projectors | projector ranks are two two one | yes | ranks=(2, 2, 1) |
| extendibility | Schur scalar has zero H trace | yes | Tr(H)=0.0; Tr(H cI)=0.000e+00 |
| extendibility | Higgsed light selector has H trace two | yes | Tr(H P+)=2.000000000000 |
| extendibility | Schur idempotent ranks exclude rank two | yes | ranks={0.0: 0, 1.0: 5} |
| conditional_higgsed_source | conditional Higgsed G3 copy parity is even | yes | rank(P+) mod 2=0 |
| conditional_higgsed_source | conditional Higgsed G2 copy parity is even | yes | light=4; heavy=4; parity=0 |
| conditional_higgsed_source | copy-counting parity prediction is derived as plus plus | yes | derived pair=((-1)^0,(-1)^0)=(1, 1) |
| heavy_index | charged plus neutral gravitational coefficient | yes | -1/12-1/24=-1/8 |
| heavy_index | symbolic degree-four Chern expansion gives minus c2 plus x squared | yes | (2-c2)(1-x+x^2/2) [degree four]=(-1)c2+(1)x^2 |
| heavy_index | finite spin characteristic grid spot-checks heavy-index integrality | yes | grid points=567; all denominators one |
| heavy_index | declared generator data directly evaluate to the pure pair plus plus | yes | indices={'G3': Fraction(0, 1), 'G2': Fraction(0, 1)}; pair=(1, 1) |
| mixed_wzw | complete five charge trace vanishes | yes | charges=[1, 1, -1, -1, 0]; sum=0 |
| mixed_wzw | light heavy split is plus two minus two | yes | light=2; heavy=-2; total=0 |
| declared_scope | derived and conditional-source gates match the declared true set | yes | true gates=['common_quadratic_aps_pv_regulator_specified', 'conditional_higgsed_gw_source_specified', 'conditional_higgsed_parity_predicts_plus_plus', 'euclidean_gamma_reference_reality_data_specified', 'exact_separated_rank_two_projector_obstruction_proven', 'five_dimensional_aps_domain_specified', 'physical_fermion_phase_zero_cut_fixed', 'pure_heavy_phase_on_G2_computed', 'pure_heavy_phase_on_G3_computed', 'pure_unbroken_gauge_gravity_heavy_phase_computed'] |
| declared_scope | unconstructed operators, physics promotion, and portal stay false | yes | false gates=['degree_one_route_e_portal_started', 'finite_higgsed_overlap_dirac_yukawa_operator_defined', 'five_dimensional_higgsed_target_kernel_defined', 'fixed_adjoint_mass_kramers_reality_proven', 'full_nonperturbative_sp4_measure_constructed', 'k_plus_two_all_scale_matched', 'minimal_higgsed_overlap_eta_pair_computed', 'mixed_wzw_heavy_determinant_computed', 'original_action_target_mass_family_lift_defined', 'original_deep_sp4_action_uniquely_selects_target_pair', 'physics_promotion_authorized', 'sp4_aps_pv_lane_closed', 'target_mapping_torus_eta_pair_selected'] |
| declared_scope | Sp4 APS/PV lane remains fail closed | yes | computed lane closure=False |

## Source manifest

| path | sha256 |
|---|---|
| route_f/tex/ap_e7_common_sp4_aps_pv_regulator.tex | 5f4012b72a335895ed946250a693420915ad1cdc959e9e33de426d5154e367a2 |
| route_f/tex/ap_e7_common_sp4_aps_pv_regulator.bib | e1c251769526c57a21d51cb33d788956e3529190e92dc03900c62420503e111d |
| route_f/code/verify_ap_e7_common_sp4_aps_pv_regulator.py | 9834cfdf8a251656a0219ae791f7365ea4dc6303b83f8d7b35da245e8d26c4a3 |
