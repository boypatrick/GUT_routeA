# DYN-4a: Seesaw Replay, Target Refresh, and the zeta Posterior

`19/19` checks passed.

## Provenance closed

The pruned benchmark artifacts were restored by the author under
route_E/output/ (sha256 recorded).  Using ONLY the recipe preserved in HEAD
code, this audit regenerates the paper's reproducibility anchors from the
archival Yukawas: zeta = 0.1076472949 + 0.0736514853 i and contact
fraction 0.1304275 -- both match the printed digits.  The
forward seesaw closes on the old targets to 1e-9, and the heavy Majorana
spectrum (3.93e+15, 3.22e+13, 2.38e+10 GeV) is exported
for DYN-7.

## Target refresh (NuFit-6.0, arXiv:2410.05380, NO)

Old hardcoded benchmark vs IC24-with-SK-atm: the sin^2 theta_23 octant flip
dominates (+6.9 sigma); other shifts are mild.  Both
NuFit variants are recorded in the ledger.

## zeta posterior (conditional on the benchmark Dirac sector)

Benchmark chi^2 vs NuFit-6.0 = 47.5 (pull table in ledger,
theta_23-driven).  Scanning M_R(zeta') = M_*(M_V + zeta' c_hat):
best zeta' = 0.107647 + 0.073651 i with
chi^2_min = 47.46; profile 68% intervals
|zeta| in [0.13043, 0.13043],
arg zeta in [0.60004, 0.60004] rad;
distance from the benchmark 0.00000.
Residual pulls at the best point are in the ledger: with the Dirac sector
frozen, the zeta direction can only partially absorb the octant flip --
the remainder is DYN-4b's covariant refit problem.
delta_CP (sin branch) at the best point: -36.0 deg vs NuFit
212 +/- 33.5 deg (prediction,
not fitted).

## Windows

Historical loose windows REGRESS to the paper anchors:
Delta_s = 3.208798e-05 (anchor 3.208798e-5),
Delta_delta = 5.424119e-05 rad (anchor 5.424119e-5).
Refreshed Delta-chi^2 <= 1 windows at the new best point:
Delta_s = 1.008e-06, Delta_delta = 1.506e-06 rad.

## Boundary

Posterior conditional on benchmark Y_nu, Y_e (full covariant refit =
DYN-4b); NO assumed, m1 = 1e-3 eV fixed; v_u = 100 GeV archival convention;
a fit, not a prediction -- Theorem VII.4 stands.
