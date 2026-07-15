# AP-E completion-frontier promotion gate

- Mechanical checks: **20/20**
- 4D lattice/full-Hessian lane closed: **False**
- APS/semisimple all-scale lane closed: **False**
- Same-soliton SQM/Callias/descent lane closed: **False**
- Any pre-portal lane closed: **False**
- Degree-one portal authorized: **False**
- Degree-one portal constructed: **False**
- Physics promotion allowed: **False**

## Decision

No pre-portal lane is closed.  The degree-one Route-E portal is therefore deliberately not constructed; only its theorem-level prerequisites are recorded for the next audit.

## Portal theorem prerequisites

- `at_least_one_preportal_lane_closed`: `false`
- `same_moduli_space_identified`: `false`
- `uniform_gap_below_portal_scale`: `false`
- `gauge_basic_wzw_line_o_plus_2`: `false`
- `physical_tangent_or_callias_line_matched`: `false`
- `physical_cpt_anti_canonical_map`: `false`
- `anomaly_free_portal_operator`: `false`
- `degree_one_map_and_orientation_proven`: `false`

## Checks

- [PASS] `provenance` - all three audit source sets exist and are hashed: hashed=10/10
- [PASS] `source_cards` - all three mechanical cards are green: lattice=24/24; APS=58/58; SQM=27/27
- [PASS] `lattice` - sampled B=1 background has a controlled continuum extrapolation: B0=0.999981468235
- [PASS] `lattice` - gauge and ghost determinant xi dependence cancels up to the predicted constant: residual=1.421e-14
- [PASS] `lattice` - declared meson and diquark blocks are positive with a tachyonic control: meson=0.621297848; diquark=0.694069173; control=-0.349308274
- [PASS] `lattice` - Wilson-Dirac Schur and one-loop curvature controls agree: sigma_min=0.702570766; Schur=6.301e-16; curvature=1.316e-08
- [PASS] `lattice` - Wilson and tree-Symanzik scaling controls have orders two and four: pW=1.99860957; pS=3.99627578
- [PASS] `lattice_boundary` - finite diagnostic is not promoted to a dynamical 4D/full-Hessian result: all six microscopic/full-closure flags remain false
- [PASS] `aps` - both reduced bordism generators and reference product spectra are explicit: plain phases={'G2_T2_odd_times_S2': 1, 'G3_S1_periodic_times_S3': 1, 'reason': 'paired nonzero ambient spectrum gives eta=0'}
- [PASS] `aps` - defect regulators realize all four torsion characters: characters=[(-1, -1), (-1, 1), (1, -1), (1, 1)]
- [PASS] `aps` - charged-QC2D does not yet select either torsion character: reference/defect calculations establish regulator dependence, not UV selection
- [PASS] `semisimple` - the simply connected Sp(4) candidate passes the scanned algebraic screens: candidate=Sp(4)=USp(4)=Spin(5), simply connected; b0=15/2
- [PASS] `aps_boundary` - Sp(4) remains a tree-level candidate rather than an all-scale UV completion: APS selection, gauge bordism, protected thresholds, eta matching, and strong dynamics are open
- [PASS] `same_soliton` - conditional Callias boundary class gives rank one and determinant c1=+2: template={'determinant_c1': 2, 'rank': 1}
- [PASS] `same_soliton` - Serre-dual anti-canonical CPT sector has three states while O(-2) has one: plus={'negative': 0, 'positive': 3}; minus={'negative': 3, 'positive': 0}; naive={'negative': 1, 'positive': 0}
- [PASS] `same_soliton` - spatial differential pushforward has degree two and conditional Chern number +2: degree=5-3=2; c1=2
- [PASS] `same_soliton` - gauge-basic descent conditions are explicit and all remain unsupplied: conditions=['class_in_image_of_pullback_from_quotient', 'curvature_horizontal_and_invariant', 'equivariant_differential_refinement', 'stabilizer_character_trivial', 'vertical_and_large_gauge_holonomy_trivial']
- [PASS] `same_soliton_boundary` - conditional template is not substituted for a same-mother-model derivation: all eleven same-model closure flags remain false
- [PASS] `portal` - no lane closure means no degree-one portal construction: routes={'four_dimensional_lattice_and_full_hessian': False, 'aps_and_semisimple_uv': False, 'same_soliton_sqm_callias_descent': False}
- [PASS] `promotion` - all source cards and the integrated frontier remain non-promoting: physics_promotion_allowed=false
