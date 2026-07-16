# AP-E7 discrete re-relaxation and same-grid operator audit

- Status: `mechanical_pass_physics_fail_closed`
- Checks: `16/16`
- Finite-site configuration space path-connected: `true`
- Exact unconstrained lattice B sector exists: `false`
- Unconstrained minimizers unwind to B=0: `true`
- Admissible interior B=1 stationary point found: `false`
- Artificial 3x3x3-core constrained solver converged: `true`
- Exact n+s sparse Hessian in the artificial constrained sector: `true`
- Same-grid declared proxy blocks assembled: `true`
- Physical SU(2) gauge-ghost blocks assembled: `false`
- Physical aggregate super-Hessian complete: `false`
- 4D dynamics performed: `false`
- Physics promotion allowed: `false`
- Portal start allowed: `false`
- Lane closed: `false`

## Unconstrained same-action minimization

| N | L | a | E_final | free gradient density | B_geom | admissibility margin |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 6.0 | 1.000 | 0.00000000 | 1.990e-09 | 0 | 1.000e+00 |
| 9 | 6.0 | 0.750 | 0.00000000 | 3.913e-09 | 0 | 1.000e+00 |
| 11 | 6.0 | 0.600 | 0.00000000 | 3.170e-09 | 0 | 1.000e+00 |
| 9 | 8.0 | 1.000 | 0.00000000 | 1.844e-08 | 0 | 1.000e+00 |
| 11 | 10.0 | 1.000 | 0.00000000 | 2.135e-09 | 0 | 1.000e+00 |

## Admissible-component guard at N=7, L=6

| floor | accepted steps | E_final | gradient density | final margin | B_geom | termination |
|---:|---:|---:|---:|---:|---:|:---|
| 0.100 | 21 | 28.63959591 | 4.729e-01 | 0.100000002 | 1 | admissibility_boundary_blocked_descent |
| 0.050 | 25 | 26.04401742 | 4.377e-01 | 0.050000002 | 1 | admissibility_boundary_blocked_descent |
| 0.020 | 28 | 24.72976823 | 4.178e-01 | 0.020000002 | 1 | admissibility_boundary_blocked_descent |

## Artificial fixed-core control only

| N | L | a | iters | E | free gradient density | B_geom | B_derivative | restricted n+s lambda0 | Delta lambda0 | sigma_min(D_W) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 6.0 | 1.000 | 19 | 43.53504850 | 1.454e-08 | 1 | 0.166471 | 2.2643404 | 1.5246250 | 1.5052234 |
| 9 | 6.0 | 0.750 | 28 | 52.59825336 | 3.333e-08 | 1 | 0.347786 | 1.8043338 | 1.5386160 | 1.4631199 |
| 11 | 6.0 | 0.600 | 43 | 57.29959805 | 6.107e-08 | 1 | 0.477536 | 1.5002079 | 1.5508433 | 1.4323985 |
| 9 | 8.0 | 1.000 | 27 | 42.92064261 | 9.461e-09 | 1 | 0.166379 | 1.1610706 | 1.2193475 | 1.4466224 |
| 11 | 10.0 | 1.000 | 33 | 42.78844632 | 1.023e-08 | 1 | 0.166440 | 0.7562822 | 1.0759295 | 1.4213190 |

For fixed boundary, the finite-site space is a finite product of connected
spaces.  It therefore has no exact pi_3 sector.  The regular-value degree is
locally constant only after excluding tetrahedra whose convex hull meets the
origin.  Unconstrained minimization crosses precisely that dislocation set.
Every accepted guard endpoint has B=1, but each search reaches its imposed
admissibility floor with a nonzero gradient.  The fixed 3x3x3 core is only an artificial
solver and sparse-operator control, never a genuine unanchored B=1 solution.
The Wilson/Yukawa normal operator is not the bosonic second variation of
`-log det D_W`; the physical aggregate and four-dimensional gates remain false.

## Mechanical checks

- [PASS] `P0_provenance` - AP-E7 source, imported AP-E6 action, TeX, and bibliography are present and hashed: hashed=4/4
- [PASS] `P0_provenance` - the AP-E6 canonical continuum profile checksum is reproduced: checksum=81104f59a5f3f4337739fc2a217cafc8da91581c9f1fb40b4fdba6ce0f948d59
- [PASS] `P1_relaxation` - the independent chart gradient agrees with a centred directional derivative: relative residual=7.313e-11
- [PASS] `P1_geometry` - a regular tetrahedron whose equal-weight convex combination is zero is rejected as non-admissible: margin=2.212e-16; admissible=False
- [PASS] `P2_configuration_space` - the fixed-boundary finite-site configuration space is path-connected: C_N=(S^3)^(N_int) x R^(N_int); finite products of path-connected spaces are path-connected
- [PASS] `P2_unwind` - unconstrained same-action minimization converges but unwinds every tested B=1 sample to B=0: max gradient=1.844e-08; final B rows=[0, 0, 0, 0, 0]
- [PASS] `P2_unwind` - pinning only the centre to n=-1 does not define a B=1 lattice sector: final B=[0, 0, 0]
- [PASS] `P3_admissible_guard` - three regular targets give B=1 at every accepted endpoint-guard iterate: B rows=[[1, 1, 1], [1, 1, 1], [1, 1, 1]]
- [PASS] `P3_admissible_guard` - each guarded descent reaches its admissibility floor before an interior stationary point: (floor,margin,gradient)=[(0.1, 0.10000000187060631, 0.4729201941100107), (0.05, 0.05000000222274066, 0.4376742955004592), (0.02, 0.020000002121115898, 0.4178457722209829)]
- [PASS] `P4_artificial_control` - each fixed-3x3x3-core control converges in its explicitly restricted sector: max free gradient=6.107e-08
- [PASS] `P4_artificial_control` - the artificial controls retain the three-target B=1 preimage count: B rows=[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]
- [PASS] `P5_sparse_hessian` - each artificial-sector exact tangent n+s Hessian is symmetric and positive: lowest restricted values=[2.2643404361070756, 1.8043338467253038, 1.5002079240515442, 1.1610706047982788, 0.7562821876616538]
- [PASS] `P5_sparse_hessian` - a core-compatible geodesic second difference reproduces every analytic n+s Hessian: best residuals=[6.741557463411512e-10, 1.055946123290509e-09, 2.44358253583415e-09, 1.3313563321569516e-09, 9.750525479804493e-10]
- [PASS] `P5_same_grid_blocks` - diquark, gauge, ghost, and declared Wilson/Yukawa normal operators are positive on every identical grid: minima Delta=1.07593, gauge=0.293661, sigma(D)=1.42132
- [PASS] `P5_same_grid_blocks` - every declared proxy has the advertised same-grid multiplicity: multiplicities=(Delta 2M, abelian spatial gauge 3M, abelian ghost M, fermion 8M)
- [PASS] `P5_determinism` - an independent repeat of the artificial N=7,L=6 control is byte-stable after canonical JSON serialization: first=27c5e4f06d5a27dc8866bd8bad32c8205a0d5f9d7ada10a9eb1dfecb6c4ca6c8; repeat=27c5e4f06d5a27dc8866bd8bad32c8205a0d5f9d7ada10a9eb1dfecb6c4ca6c8

## Fail-closed boundary

- `genuine_unanchored_B1_lattice_stationary_point_found`: `false`
- `full_unanchored_discrete_stationarity_achieved`: `false`
- `admissibility_independent_continuum_B1_solution_constructed`: `false`
- `continuous_line_search_segments_certified_admissible`: `false`
- `physical_su2_gauge_ghost_blocks_assembled`: `false`
- `fermion_determinant_bosonic_second_variation_computed`: `false`
- `interacting_gauge_meson_cross_blocks_computed`: `false`
- `brst_superdeterminant_on_the_relaxed_soliton_computed`: `false`
- `physical_aggregate_superhessian_complete`: `false`
- `projected_stability_continuum_extrapolated`: `false`
- `four_dimensional_dynamics_performed`: `false`
- `importance_sampling_or_hmc_performed`: `false`
- `quantum_continuum_limit_proven`: `false`
