# AP-E3 charged two-colour nonlinear-EFT proxy

- Status: `pass`
- Model class: `charged-two-colour nonlinear EFT proxy; not lattice QC2D`
- Checks: `18/18`
- Nonlinear EFT proxy pass: `true`
- Radial Hessian necessary gate: `true`
- Full three-dimensional Hessian closed: `false`
- Lattice QC2D closed: `false`
- Physics promotion allowed: `false`

This is a deterministic classical low-energy proxy, **not** a lattice or
continuum nonperturbative solution of the microscopic charged two-colour gauge
theory.  The word nonlinear refers to solving the full classical hedgehog BVP,
not to quantum nonperturbative closure.

## Homogeneous charged phase scan

| mu_lift^2 | phase | sigma | |d| | lowest vacuum Hessian |
|---:|---|---:|---:|---:|
| -0.40 | diquark_condensed_U1g_broken | 0.625000 | 0.796477 | 0.000000e+00 |
| -0.30 | diquark_condensed_U1g_broken | 0.833333 | 0.560258 | 0.000000e+00 |
| -0.26 | diquark_condensed_U1g_broken | 0.961538 | 0.277688 | 0.000000e+00 |
| -0.25 | critical | 1.000000 | 0.000000 | 0.000000e+00 |
| -0.10 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 1.500000e-01 |
| 0.00 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500000e-01 |
| 0.20 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500000e-01 |
| 0.56 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500000e-01 |
| 1.00 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500000e-01 |

The analytic transition is `mu_lift^2=-m_pi^2=-0.250000`.
At the benchmark, `(m_pi,m_sigma,m_Delta)=[0.5, 3.5, 0.9]`.

## Declared charged-coordinate chart

The chart uses the explicitly non-derived matching ansatz
`mu_lift^2/Lambda_c^2 = mu_PG,0^2 + c_g (e_g v_g/Lambda_c)^2`.
It is an EFT sensitivity scan, not a microscopic prediction.

| e_g | v_g/Lambda_c | matched mu_lift^2 | phase | sigma | |d| | min Hessian |
|---:|---:|---:|---|---:|---:|---:|
| 0.00 | 0.50 | -0.4000 | diquark_condensed_U1g_broken | 0.625000 | 0.796477 | 0.000e+00 |
| 0.00 | 1.00 | -0.4000 | diquark_condensed_U1g_broken | 0.625000 | 0.796477 | 0.000e+00 |
| 0.00 | 2.00 | -0.4000 | diquark_condensed_U1g_broken | 0.625000 | 0.796477 | 0.000e+00 |
| 0.00 | 3.00 | -0.4000 | diquark_condensed_U1g_broken | 0.625000 | 0.796477 | 0.000e+00 |
| 0.25 | 0.50 | -0.3963 | diquark_condensed_U1g_broken | 0.630915 | 0.791405 | 0.000e+00 |
| 0.25 | 1.00 | -0.3850 | diquark_condensed_U1g_broken | 0.649351 | 0.775141 | 0.000e+00 |
| 0.25 | 2.00 | -0.3400 | diquark_condensed_U1g_broken | 0.735294 | 0.688725 | 4.441e-08 |
| 0.25 | 3.00 | -0.2650 | diquark_condensed_U1g_broken | 0.943396 | 0.335416 | 0.000e+00 |
| 0.50 | 0.50 | -0.3850 | diquark_condensed_U1g_broken | 0.649351 | 0.775141 | 0.000e+00 |
| 0.50 | 1.00 | -0.3400 | diquark_condensed_U1g_broken | 0.735294 | 0.688725 | 4.441e-08 |
| 0.50 | 2.00 | -0.1600 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 9.000e-02 |
| 0.50 | 3.00 | 0.1400 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500e-01 |
| 0.75 | 0.50 | -0.3663 | diquark_condensed_U1g_broken | 0.682594 | 0.743936 | 0.000e+00 |
| 0.75 | 1.00 | -0.2650 | diquark_condensed_U1g_broken | 0.943396 | 0.335416 | 0.000e+00 |
| 0.75 | 2.00 | 0.1400 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500e-01 |
| 0.75 | 3.00 | 0.8150 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500e-01 |
| 1.00 | 0.50 | -0.3400 | diquark_condensed_U1g_broken | 0.735294 | 0.688725 | 4.441e-08 |
| 1.00 | 1.00 | -0.1600 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 9.000e-02 |
| 1.00 | 2.00 | 0.5600 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500e-01 |
| 1.00 | 3.00 | 1.7600 | mesonic_U1g_unbroken | 1.000000 | 0.000000 | 2.500e-01 |

## B=1 benchmark

- beta: `0.5`
- alpha: `0.9`
- chi: `-2.0`
- origin slope: `2.244905536362364`
- numerical B: `0.9999999995940666`
- B residual: `4.059333980066526e-10`
- dimensionless energy integral: `6.135914192963568`
- Derrick relative residual: `2.5223696174580667e-08`
- lowest radial squared frequency `omega^2/(ef)^2`: `0.31140056712190023`
- lowest charged-diquark squared frequency `omega^2/(ef)^2`: `0.6965061007409071`
- unstable negative-control squared frequency: `-1.4305158104185836`

The radial eigenvalue is a finite-box necessary stability test.  It is not a
claim that every nonradial meson, vector, gauge, radial-sigma, or fermion-
determinant channel has been diagonalised.  Positivity in the charged scalar
`l=0` channel does imply positivity for `l>0` within this scalar block because
the centrifugal term is nonnegative.

## Volume convergence

| L | intervals | E | B | lambda_F,min | lambda_Delta,min |
|---:|---:|---:|---:|---:|---:|
| 14 | 700 | 6.135914701 | 1.000000000 | 0.349739800 | 0.696547693 |
| 18 | 900 | 6.135914193 | 1.000000000 | 0.311400833 | 0.696489297 |
| 22 | 1100 | 6.135914184 | 1.000000000 | 0.291384645 | 0.696485381 |

## Grid convergence certificate

- Observed radial order: `1.9951686188672073`
- Observed diquark order: `2.000718983825481`
- Radial Richardson extrapolate: `0.3114002243846333`
- Diquark Richardson extrapolate: `0.6965277010072412`
- Radial fine-grid error estimate: `8.568422943389109e-08`
- Diquark fine-grid error estimate: `5.4000665677955695e-06`

## Mechanical checks

- [PASS] `P0_provenance` - proxy source, independent TeX derivation, and bibliography exist and are hashed: hashed=3/3
- [PASS] `P1_vacuum` - the declared mesonic vacuum solves the exact stationarity equation: residual=2.220e-16
- [PASS] `P1_vacuum` - finite-difference and analytic six-field vacuum Hessians agree: max residual=3.267e-09
- [PASS] `P1_vacuum` - benchmark pion, sigma, and doubly charged diquark gaps are positive: gaps=(0.500000,3.500000,0.900000)
- [PASS] `P1_vacuum` - analytic scan finds diquark onset exactly at mu_lift^2=-m_pi^2: critical=-0.250000; critical Hessian min=0.000e+00
- [PASS] `P1_vacuum` - the declared (e_g,v_g/Lambda_c) chart resolves both sides of the charged phase gate: phases=['diquark_condensed_U1g_broken', 'mesonic_U1g_unbroken']; matching=declared_EFT_ansatz_not_UV_derived
- [PASS] `P2_soliton` - nonlinear hedgehog BVP converges with the regular origin series: nodes=579; max rms=2.955e-07
- [PASS] `P2_soliton` - integrated baryon density agrees with the B=1 boundary formula: B_numeric=0.999999999594
- [PASS] `P2_soliton` - Derrick virial identity E4=E2+3E0 is numerically satisfied: relative residual=2.522e-08
- [PASS] `P2_soliton` - scale Hessian is positive at the Skyrme saddle: d2E/dlambda2=8.145879695
- [PASS] `P2_soliton` - negative control: deleting E4 removes the finite-scale Derrick stationary point: d(E2*lambda+E0*lambda^3)/dlambda|1=3.319202687
- [PASS] `P3_hessian` - all computed spherically symmetric meson Hessian eigenvalues are positive: lowest five=[0.31140056712190023, 0.4257651809692827, 0.5821246592778764, 0.7752867401972081, 1.024793300584382]
- [PASS] `P3_hessian` - the doubly charged l=0 soliton mode is positive and bound below its continuum threshold: lambda0=0.696506101; continuum=0.810000000; free-box=0.840461742
- [PASS] `P3_hessian` - centrifugal ordering makes l=0 the most dangerous charged scalar channel: lambda_l0=0.696506101; lambda_l1=0.872004739
- [PASS] `P3_hessian` - negative control: a stronger attractive core coupling produces a charged tachyon: chi=-5; lambda0=-1.430515810
- [PASS] `P4_convergence` - the N=600,1200,2400 sequence shows second-order Hessian convergence and a final change below 2e-5: delta_radial=2.571e-07, p=1.9952; delta_diquark=1.620e-05, p=2.0007
- [PASS] `P4_convergence` - the final box enlargement stabilises soliton energy and the localized diquark eigenvalue: delta_E=8.679e-09; delta_diquark=3.916e-06
- [PASS] `P4_convergence` - the lowest radial finite-box state decreases toward the pion continuum beta^2 without crossing zero: sequence=[0.3497397998444644, 0.31140083333049695, 0.29138464523810853]; beta^2=0.250000000

## Remaining blockers

- Replace the classical EFT proxy by a controlled lattice or other nonperturbative calculation of the charged SU(2)c gauge theory.
- Compute the full three-dimensional coupled meson, radial-sigma, vector, gauge, and ghost Hessian, including collective zero modes.
- Include the fermion determinant and quantum corrections; the present Hessian is tree-level classical.
- Prove that the microscopic charged theory dynamically selects the mesonic vacuum and match mu_lift_sq and chi from UV observables.
- Exclude global escape/unwinding paths in the faithful full target; positivity of local Hessian blocks is only a necessary condition.
- No Route-E family portal follows from this proxy scan.
