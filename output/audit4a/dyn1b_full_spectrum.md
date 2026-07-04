# DYN-1b: Complete CMSGUT Heavy Spectrum -- Derivation Record

`16/16` checks passed.

## How the spectrum was derived (two independent computations)

1. **Decomposition side (derived here).**  A PS -> SM engine from first
   principles: SU(4) -> SU(3) x U(1)_B-L (4 = 3_1/3 + 1_-1; 6 = 3_-2/3 +
   3b_2/3; 10 = 6_2/3 + 3_-2/3 + 1_-2; 15 = 8 + 3_4/3 + 3b_-4/3 + 1) and
   SU(2)_R -> T3R, with Y_SM = T3R + (B-L)/2.  Applied to
   210 = (1,1,1)+(15,1,1)+(15,1,3)+(15,3,1)+(6,2,2)+(10,2,2)+(10b,2,2),
   126 = (6,1,1)+(10,3,1)+(10b,1,3)+(15,2,2), 126bar (conj pattern),
   10 = (1,2,2)+(6,1,1): total 472 chiral states.  The 45 gives 12 unbroken
   + 33 broken gauge states.
2. **Transcription side (AG hep-ph/0405074, sha 4022df72...).**  Table I
   (lines 2439-2531): 21 unmixed masses, 19 letter types; mixed pure chiral
   R/h/t; mixed chiral-gauge G/E/F/J/X with gauge masses m_lambda (items
   i-v, lines 736-830); Table 2 (lines 2540-2567): threshold indices.

**Completeness gate:** the sector-census multiset EQUALS the decomposition
multiset (all 26 types, 472 states, no sector missing or double-counted),
and the gaugino sectors equal 45 minus the SM adjoint (33).

## Key verifications

- Threshold indices {S3,S2,S1} recomputed from Dynkin indices and
  Y_SM = Y_AG/2 match AG Table 2 for all 26 types, including the printed
  S_W = 4S1-9.6S2+5.6S3 and S_X = 5S1+3S2-8S3 columns.
- R[8,1,0]: AG's closed-form eigenvalues match the 2x2 matrix.
- The five m_lambda formulas equal the gaugino-column norms of the
  transcribed blocks to 1e-12 (super-Higgs bookkeeping closed).
- At the benchmark (x = 0.1, M_H det-tuned) exactly one massless level
  exists: the MSSM Higgs pair.  The 33 Goldstones are eaten.
- Vacuum-consistency insight: on the branch M + eta(p+3a-6w) = 0, so
  m_A = 24 eta w exactly.
- Accidental zero-mass loci of unmixed sectors inside the SM window are
  scanned and disclosed in the ledger (DYN-2 must avoid or handle them).

## Output

`output/audit4a/heavy_spectrum.json` now covers ALL 26 sector types
(57 states): 21 unmixed Dirac/Majorana levels, R eigenvalues,
4 doublet + 5 triplet levels, all mixed-block levels, and the five massive
gauge supermultiplets.  Masses are benchmark units of m; GeV normalization,
vector-multiplet threshold bookkeeping, and the AG Table-2 cross-check of
per-threshold b-vectors are DYN-2's job.
