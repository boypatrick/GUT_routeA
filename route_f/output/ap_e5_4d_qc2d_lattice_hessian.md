# AP-E5 4D charged-QC2D lattice/Hessian diagnostic

**Status:** PASS (24/24)

This is a deterministic finite-volume background-field diagnostic. It is not a Monte Carlo or dynamical-lattice result.

## Numerical certificate

- B=1 extrapolation: `0.999981468235` (absolute error `1.853e-05`).
- Gauge/ghost xi-cancellation residual: `1.421e-14`.
- Meson / diquark minimum eigenvalues: `0.621297848377` / `0.694069172619`.
- Diquark negative-control minimum: `-0.349308273811`.
- Wilson-Dirac minimum singular value: `0.702570765924`.
- Schur determinant / solve residuals: `4.263e-14` / `6.301e-16`.
- Fermion log-det curvature identity relative error: `1.316e-08`.
- Observed Wilson / tree-Symanzik orders: `1.99860957` / `3.99627578`.

## Fail-closed boundary

- `importance_sampling_or_Monte_Carlo_performed`: `false`
- `dynamical_4d_lattice_QC2D_phase_established`: `false`
- `fermion_determinant_positivity_proven`: `false`
- `complete_nonlinear_gauge_meson_ghost_fermion_Hessian_computed`: `false`
- `Skyrme_Jacobi_translation_isorotation_zero_modes_resolved`: `false`
- `one_loop_renormalised_continuum_limit_established`: `false`
- `Finkelstein_Rubinstein_constraint_derived`: `false`
- `finite_amplitude_global_stability_established`: `false`
- `degree_one_route_E_portal_allowed`: `false`
- `physics_promotion_allowed`: `false`

## Checks

- [x] **provenance / all declared source files exist** — route_f/code/verify_ap_e5_4d_qc2d_lattice_hessian.py, route_f/tex/ap_e5_4d_qc2d_lattice_hessian.tex, route_f/tex/ap_e5_4d_qc2d_lattice_hessian.bib
- [x] **provenance / deterministic seed recorded** — numpy Generator seed=20260716
- [x] **topology / finite-lattice B estimator approaches one monotonically** — [0.7804621958442137, 0.8957980638679076, 0.9400008778824941, 0.972881794380486, 0.9846555997802079]
- [x] **topology / B=1 continuum extrapolation** — B0=0.999981468235
- [x] **topology / B estimator has second-order convergence** — p=1.97945840
- [x] **gauge_ghost / four gauge toron zero modes per generator** — [4, 4, 4, 4]
- [x] **gauge_ghost / one constant ghost zero mode per generator** — [1, 1, 1, 1]
- [x] **gauge_ghost / all nonzero gauge and ghost modes positive** — prime spectra only
- [x] **gauge_ghost / xi dependence is the predicted background-independent constant** — max residual=1.421e-14
- [x] **boson_hessian / declared meson tangent block is positive** — lambda_min=0.621297848377
- [x] **boson_hessian / charge-two diquark block is positive** — lambda_min=0.694069172619
- [x] **boson_hessian / Delta=0 removes classical gauge-diquark mixing** — K_A_Delta=0 at the benchmark background
- [x] **negative_control / over-attractive core coupling creates a negative diquark mode** — lambda_min=-0.349308273811; count=1
- [x] **fermion / Euclidean Clifford algebra** — error=0.000e+00
- [x] **fermion / finite Wilson operator is nonsingular** — sigma_min=0.702570765924
- [x] **fermion / even-odd Schur determinant factorisation** — logabs=4.263e-14; phase=5.194e-17
- [x] **fermion / even-odd Schur solve agrees with direct solve** — relative residual=6.301e-16
- [x] **fermion / fermion log-determinant curvature trace identity** — relative error=1.316e-08
- [x] **scaling / Wilson dispersion has second-order scaling** — p=1.99860957
- [x] **scaling / tree-level Symanzik dispersion has fourth-order scaling** — p=3.99627578
- [x] **fail_closed / no Monte Carlo claim** — deterministic background-field diagonalisation only
- [x] **fail_closed / determinant positivity remains open** — no antiunitary/global-form proof is supplied
- [x] **fail_closed / FR and global stability remain open** — local Gaussian spectra cannot close configuration-space topology
- [x] **fail_closed / portal and physics promotion remain false** — necessary finite-lattice diagnostic only
