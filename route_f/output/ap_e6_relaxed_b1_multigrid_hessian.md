# AP-E6 relaxed B=1 multigrid sparse-Hessian audit

- Status: `mechanical_pass_physics_fail_closed`
- Checks: `20/20`
- Same relaxed B=1 solution computed: `true`
- Multigrid/multivolume converged: `false`
- Aggregate same-grid sparse projected Hessian complete: `false`
- Exact n+s sparse second variation complete: `true`
- Projected-Hessian stability gate: `false`
- Discrete stationarity achieved: `false`
- Lattice topology converged to B=1: `false`
- Off-shell projected-curvature diagnostic: `true`
- Full gauge-meson-ghost-fermion Hessian complete: `false`
- 4D importance sampling performed: `false`
- Lane closed: `false`

This card is a deterministic classical calculation.  Exact sparse second
differentiation refers only to the declared `n+s` energy.  The charge-two
diquark and periodic free-gauge/ghost operators are separate algebraic audits,
not one same-grid assembled Hessian.  The continuum-relaxed BVP is
sampled, not re-relaxed, on each cubic grid: the cubic spectra are therefore
off-shell curvature diagnostics, not a physical stability Hessian.  The card
does not include a fermion determinant, dynamical colour links, HMC, or a
quantum continuum extrapolation.

## Coupled BVP volume sequence

| R | nodes | E/(4pi f/e) | s(0) | B | virial/E |
|---:|---:|---:|---:|---:|---:|
| 8 | 521 | 6.1284788541 | 0.917136444 | 0.999999999832 | 3.665e-04 |
| 10 | 575 | 6.1282475647 | 0.917129875 | 1.000000000027 | 5.457e-05 |
| 12 | 565 | 6.1282195127 | 0.917129078 | 0.999999999717 | 8.168e-06 |
| 16 | 580 | 6.1282155251 | 0.917128965 | 0.999999999686 | 1.796e-07 |

## Fixed-volume grid sequence

| N | L | a | B_lattice | E_lattice | residual density | projected lambda0 | Delta lambda0 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 9 | 8.0 | 1.000000 | 0.168017735 | 42.81849799 | 3.477e-01 | -10.79221928 | 1.21908897 |
| 13 | 8.0 | 0.666667 | 0.440012354 | 56.70681532 | 4.302e-01 | -23.42271592 | 1.22357411 |
| 17 | 8.0 | 0.500000 | 0.625354259 | 64.24791193 | 4.217e-01 | -28.44285305 | 1.22488418 |

## Fixed-spacing volume sequence

| N | L | a | B_lattice | projected lambda0 | Delta lambda0 |
|---:|---:|---:|---:|---:|---:|
| 9 | 4.0 | 0.500000 | 0.621121457 | -28.41847154 | 2.45369874 |
| 13 | 6.0 | 0.500000 | 0.625215361 | -28.44004659 | 1.54065100 |
| 17 | 8.0 | 0.500000 | 0.625354259 | -28.44285305 | 1.22488418 |

## Negative controls

- Strong-core diquark eigenvalue: `[-1.6342233282777296, 1.1811017338409013, 1.181101733840908, 1.4300837867734544, 2.0814974451051644]`
- Projector with one omitted collective direction has leakage: `1.000000000000004`
- No-Skyrme Derrick derivative: `3.3082997659535196`

## Mechanical checks

- [PASS] `P0_provenance` - script, independent derivation, and bibliography exist and are hashed: hashed=3/3
- [PASS] `P0_provenance` - the public same-solution API produces canonical r-F-s arrays with a SHA-256 checksum: points=4097; checksum=81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59
- [PASS] `P1_relaxation` - the coupled F-s BVP converges from the regular origin series: max residual=2.956e-07
- [PASS] `P1_relaxation` - the relaxed solution has B=+1 and the breathing amplitude remains positive: B=0.999999999686; min s=0.917128965
- [PASS] `P1_relaxation` - the generalized Derrick virial identity is satisfied: virial/E=1.796e-07; curvature=8.08175706
- [PASS] `P2_volume` - R=12 to R=16 stabilises energy and the coupled core field: delta E=3.988e-06; delta s0=1.133e-07
- [PASS] `P3_multigrid` - the fixed-volume B estimator improves monotonically with resolution: errors=[0.8319822650099871, 0.5599876462859643, 0.37464574146737817]
- [PASS] `P3_multigrid` - the exact off-shell projector detects a negative meson-breathing curvature while the diquark block stays positive: projected negative=[-10.792219279276424, -23.42271591904794, -28.44285305075885]; diquark positive=[1.2190889714700144, 1.2235741122210273, 1.224884179556587]
- [PASS] `P3_multigrid` - the final refinement exposes nonconvergence of the off-shell curvature while the diquark gap remains controlled: delta projected=5.020e+00; delta diquark=1.310e-03
- [PASS] `P3_multigrid` - the negative mode is strongly localised and is not misreported as a continuum collective zero mode: inverse participation ratios=[0.5079310927762034, 0.30844271552165403, 0.1421571026828438]
- [PASS] `P4_multivolume` - fixed-spacing volume scan reproduces the negative projected mode and positive diquark gap: projected negative=[-28.418471544202234, -28.440046592199526, -28.44285305075885]; diquark positive=[2.453698743285602, 1.540651003481715, 1.224884179556587]
- [PASS] `P4_multivolume` - the negative mode stabilises and the unbound diquark level has the expected 1/L^2 threshold extrapolation: delta projected=2.806e-03; Delta extrapolate=0.813709480; fit rms=1.278e-03
- [PASS] `P5_projector` - translation-isorotation projector is idempotent and horizontal: idem=6.791e-17; orth=6.224e-15; eig overlap=2.275e-16
- [PASS] `P5_projector` - the exact n+s sparse Hessian is symmetric: dofs=13500; nnz=395538; residual=6.538e-17
- [PASS] `P5_projector` - an independent exponential-map energy difference reproduces v^T H_tan v: best relative residual=1.023e-09
- [PASS] `P6_gauge_ghost` - toron and constant-ghost projectors are idempotent and orthogonal: max residual=8.687e-16
- [PASS] `P6_gauge_ghost` - free gauge/ghost determinant ratio obeys the xi cancellation law: max xi residual=1.137e-13
- [PASS] `P7_negative_controls` - a strong attractive diquark core creates a tachyon: lambda0=-1.634223328
- [PASS] `P7_negative_controls` - omitting one collective vector is detected by unit leakage: leakage=1.000000000000
- [PASS] `P7_negative_controls` - removing the Skyrme term destroys Derrick stationarity: dE/dscale=3.308299766

## Fail-closed boundary

- `sparse_projected_hessian_complete_in_declared_sector`: `false`
- `all_blocks_assembled_on_one_same_grid_and_boundary_condition`: `false`
- `discrete_stationarity_achieved`: `false`
- `lattice_topology_converged_to_B1`: `false`
- `physical_projected_hessian_at_stationary_discrete_solution`: `false`
- `full_three_dimensional_continuum_stationarity_proven`: `false`
- `multigrid_multivolume_converged`: `false`
- `full_gauge_meson_ghost_fermion_hessian_complete`: `false`
- `fermion_determinant_included`: `false`
- `dynamical_colour_links_included`: `false`
- `dynamical_4d_importance_sampling_performed`: `false`
- `HMC_or_Monte_Carlo_performed`: `false`
- `continuum_quantum_phase_proven`: `false`
- `microscopic_charged_QC2D_phase_proven`: `false`
- `finite_amplitude_global_unwinding_excluded`: `false`
- `degree_one_route_E_portal_allowed`: `false`
- `physics_promotion_allowed`: `false`
- `lane_closed`: `false`
