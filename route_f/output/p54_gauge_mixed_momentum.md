# Same-prescription radial gauge / mixed momentum kernel

Strict-Landau one-loop background kernel, not a gauge-independent pole or all-orders stability proof.
Checks: 99/99.

| Euclidean p2 | Bosonic generalized eigenvalues including gauge momentum |
|---:|---|
| 0.001 | -0.758574603, 0.0910321031, 1.13817219 |
| 0.01 | -0.726254818, 0.101339356, 1.73030098 |
| 0.05 | -0.589136529, 0.142533661, 2.18892505 |
| 0.1 | -0.42162596, 0.193296044, 2.42971285 |
| 0.5 | 0.596221906, 0.79553823, 3.3617128 |
| 1.0 | 1.09860159, 2.10101216, 4.19405963 |

The gauge-only local kinetic insertion has generalized eigenvalues [-0.14757559836612907, 0.0, 0.00034481327091172434].
Adding it to the previous hard scalar kinetic insertion gives [0.008995184178964154, 0.2733736138424651, 2.4665928385677534]. The soft scalar nonlocal part and unfitted fermions are not included in this local diagnostic; it is not a physical residue.

Finite-xi checks hold the original finite mass CT fixed and report residual tadpoles; convergence to Landau is not a Nielsen identity certificate.
The MS dimensional rational terms +2a (vector seagull), -2 (vector bubble), and -1/8 (mixed slope) are retained and independently checked.
The exact radial Majorana function can be added with three specified masses; no fitted values are invented.

## Still open

- fixed-bare finite-xi Nielsen/BRST pole test and physical background/quantum identification
- perturbative control, resummation/higher-loop remainder and fitted Majorana inputs
- full Wilson/box, lower gauge/ghost/Yukawa matching and finite CHN seesaw
