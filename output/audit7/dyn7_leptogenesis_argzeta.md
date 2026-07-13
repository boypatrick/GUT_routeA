# DYN-7: Historical Unflavored Leptogenesis Diagnostic

`13/13` checks passed.

## Chain (order of magnitude, SUSY loop function, BDP-type washout)

Takagi -> h = Y_nu V -> epsilon_1 (Davidson-Ibarra respected point-wise)
-> m_tilde, kappa -> eta_B = -0.96e-2 eps1 kappa.
Old benchmark: M_1 = 2.38e+10 GeV,
eps1 = +2.11e-07, m_tilde = 4.52e-02 eV,
kappa = 6.00e-03, eta_B = -1.21e-11.
Refreshed central: M_1 = 2.49e+10 GeV,
eta_B = -8.98e-12 (observed +6.1e-10).

## Regression statistics (DYN-4b priors, 4000 samples; not a posterior)

P(eta_B > 0) = 0.502; P(|eta_B| within half a decade of 6.1e-10) =
0.007; P(success) = 0.004.

## Diagnostic only: no flavored phase constraint

Unconditioned (DYN-4b): 68% [0.4669, 0.7535],
95% [0.1383, 1.0942] rad.
Leptogenesis-conditioned: 68% [0.2231, 0.7767],
95% [0.1241, 1.1358] rad.
The old benchmark phase 0.6000 is
INSIDE the conditioned 95% window.
Success drivers: cos(a21) corr -0.02, log10 m1 +0.08,
log10 M_1 +0.03.

These numbers are retained only to regression-test the historical unflavored
pipeline.  They cannot cut the contact phase in the tau-resolved regime.

## Boundary

Status: `blocked_missing_branch_thermal_inputs`.  The missing inputs are
tan(beta), branch-local charged-lepton/spectator rates, and a flavored kinetic
solution.  See `audit7_dyn7_flavor_regime_gate.py`.
