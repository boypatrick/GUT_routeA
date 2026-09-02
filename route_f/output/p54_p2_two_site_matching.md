# Hierarchical P54PQ two-site matching

Status: **11/11 checks passed**.

The calculation uses exact PS Casimir projectors and spectral traces, not a nearest-mass assignment of broken-phase states.

- PS-breaking parent: `Sigma126:(10-pair,1,3)`;
- scalar `lambda_I(3,2,1) = [31.213843839581514, 4.930380657631324e-32, 51.86695627747886]`;
- scalar `lambda_U(4,L,R) = [-8.128234982577338, -12.047442539256684, -2.5218265047839274]`;
- total `lambda_I = [31.213843839581497, 4.930380657631324e-32, 51.86695627747884]`;
- total `lambda_U = [-8.461621820713413, -12.547522796460797, -3.021906761988041]`;
- thresholded `MI=1.079859e+14 GeV`, `MU=8.534751e+14 GeV`, `alphaU^-1=39.7014114`;
- resulting `MI/MU=0.126524923` versus the input vacuum ratio `0.1265` (replay factor `1.0002`).
- covariance sigma in `(log10 MI, log10 MU, alphaU^-1)`: experiment `[0.007939101281430294, 0.025592134547277202, 0.06129503867460127]`, thresholds `[0.039516897509813166, 0.028477412750478995, 0.24056951564284654]`, total `[0.04030650714162524, 0.03828733978284576, 0.24825546040873836]`.

## PS parent census

| parent | real dimension | C4 | C2L | C2R |
|---|---:|---:|---:|---:|
| Phi54:(1,1,1) | 1 | 2.46494e-31 | 4.02182e-17 | 2.70399e-18 |
| Phi54:(6,2,2) | 24 | 2.5 | 0.75 | 0.75 |
| Phi54:(20prime,1,1) | 20 | 6 | 7.65791e-32 | 7.65791e-32 |
| Phi54:(1,3,3) | 9 | 9.29763e-31 | 2 | 2 |
| Sigma126:(6,1,1) | 12 | 2.5 | 0 | 0 |
| Sigma126:(10-pair,3,1) | 60 | 4.5 | 2 | 5.05265e-34 |
| Sigma126:(15,2,2) | 120 | 4 | 0.75 | 0.75 |
| Sigma126:(10-pair,1,3) | 60 | 4.5 | -6.9032e-34 | 2 |
| phi10:(6,1,1) | 12 | 2.5 | 0 | 0 |
| phi10:(1,2,2) | 8 | 0 | 0.75 | 0.75 |
| S:(1,1,1) | 2 | 0 | 0 | 0 |

## Checks

| check | result |
|---|---|
| PS parent projectors cover 328 real coordinates | PASS |
| fixed-point scalar ledger has 35 SM sectors | PASS |
| fixed-point Hessian has 38 zeros and no tachyon | PASS |
| four nonsymmetry zeros are one doublet | PASS |
| all PS Casimirs have recognized labels | PASS |
| 126 VEV lies in a unique PS parent | PASS |
| vector site census is 9 plus 24 | PASS |
| two-site unification residual closes | PASS |
| one replay keeps the hierarchy within a factor two | PASS |
| gU threshold iteration converges | PASS |
| two-site threshold covariance is positive semidefinite | PASS |
