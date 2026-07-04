# Route-E Dependency-Closure Audit

`24/24` checks passed.

- CI-1 removed: w0 computed for every simple algebra of rank 4-15 plus
  E6/E7/E8/F4; the self-conjugacy rules are machine-derived.
- Enumeration closed over ALL rank 4-15 simple algebras: the complex
  dim-15/16 list is unchanged; so(9)'s 16 and so(15)'s 15 are real
  near-misses (fail chirality, not dimension).
- CI-4 removed: E6 min = 27, E7 = 56, F4 = 26, E8 = 248 by scan;
  rank-16 fundamentals are 17/33/32/32 (> 16).
- CI-2 removed: SU(15)/SU(16) fundamental cubic traces are nonzero
  (exact integers).
- CI-3 removed: the chiral 16 of Spin(10) built from five fermionic
  modes has vanishing symmetrized cubic trace over all 45^3 triples;
  its Cartan weights reproduce the demicube.

Boundary: pi_4(Spin(10)) = 0 (Witten global anomaly) is cited, not
computed; Frobenius-Schur types are not computed; rank >= 17 minimal
dimensions are textbook-cited; zeta is not derived; no PSLT-only
unconditional GUT proof is claimed.
