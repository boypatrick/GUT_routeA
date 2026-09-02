# P54PQ-v2 P2 running, threshold covariance, and PQ repair

Verifier status: **23/23 checks passed**.

P2 matching gate closed with Yukawa nuisance: **True**.
Gauge-independent pole mass predicted: **False**.

## Two-loop baseline

- no-threshold current-input solution: `MI=4.609432e+13 GeV`, `MU=1.181261e+15 GeV`, `alphaU^-1=37.601980`;
- preferred `MI/MU=3.902129e-02` and `ln(MU/MI)=3.243648`;
- the superseded P1 diagnostic used `sigma/omega=0.350000`, a no-threshold mismatch factor `8.969`;
- the parent-resolved threshold fixed point is `sigma/omega=0.1265000` and returns `MI/MU=0.1265249`;
- fixed-point scales: `MI=1.079859e+14 GeV`, `MU=8.534751e+14 GeV`, `alphaU^-1=39.7014114`.

The former hierarchy blocker is closed by re-deriving the radial quadratic masses at the fixed point; the quartic/cubic P1 couplings are unchanged.

## Actual P1 threshold ledger

- rows: `35`; positive physical rows: `29`; positive real dimension: `290`;
- complete scalar trace indices `(3,2,1)=[84.0, 83.99999999999999, 83.99999999999997]`;
- positive-mode trace indices `(3,2,1)=[79.0, 76.99999999999999, 75.39999999999999]`;
- scalar mass band: `[0.116278,1.682012] omega`;
- collapsed `mu=omega` scalar lambda: `[-56.89359583551077, -56.33639848841842, -53.167052880163546]`.

The exact PS parent projectors assign the breaking `(10-pair,1,3)` sector to MI and the complementary parent space to MU. Scalar and vector mass-block Jacobians are propagated on both sites.

Two-site total sigma in `(log10 MI,log10 MU,alphaU^-1)` is `[0.04030650714162524, 0.03828733978284576, 0.24825546040873836]`; the threshold-only part is `[0.039516897509813166, 0.028477412750478995, 0.24056951564284654]`.

## Loop-corrected doublet retuning

The fixed-point light mode has 10_H weight `0.999997474290` and therefore

`delta xi02 = 1.000002525717 h^T Pi_D h/omega^2`.

The nearest heavy-doublet gap is `0.030557030 omega^2`.

## Bosonic Coleman-Weinberg matching

In `background-field Landau gauge (xi_gf=0)` and `MSbar` at `muU=gU*omega`, the hard-mode result is `kappa_B=-2.023665109501e-02 omega^2` and `delta xi02,B=-2.023670220706e-02`.

The heavy-Yukawa projection is not set to zero: `eta_Y(mu)=h^T Pi_heavy-Y(0;mu)h/omega^2`, and the matching condition is `delta xi02=[kappa_B/omega^2+eta_Y]/w10`. This is a zero-momentum matching curvature, not a pole mass.

## PQ repair

The baseline remains `N_DW=3`.  The minimal candidate extension adds two `10_F` Weyl fields with `qPQ=+2` and the mass operator `S 10_F 10_F`.  It changes `Nhat=-4.0` and gives physical `N_DW=1` while preserving one-loop differential unification.

The mass-dependent two-loop replay (the extension is not added to the baseline) gives:

| y_F | M_F/MU | MI/MU | alphaU^-1 |
|---:|---:|---:|---:|
| 0.5 | 0.156270 | 0.125482 | 39.279799 |
| 1 | 0.313169 | 0.125873 | 39.438279 |
| 2 | 0.627590 | 0.126263 | 39.596062 |

## Blocker and next action

there is no remaining bosonic blocker to opening P3; a fully predictive pole mass remains unavailable until P3 fits the heavy-Yukawa nuisance and supplies momentum-dependent matching

begin P3 with eta_Y(mu) profiled as a matching nuisance and test the corrected heavy-doublet block against the 0.0305570303 omega^2 gap

## Verification

| Group | Check | Result |
|---|---|---|
| input | P1 spectrum contains exactly 35 SM-Casimir rows | PASS |
| input | P1 spectrum covers all 328 real coordinates | PASS |
| threshold | physical positive scalar threshold census has 290 real dimensions | PASS |
| threshold | complete Spin(10) scalar trace indices are universal and equal to 84 | PASS |
| rge | two-loop solver closes both unification differences | PASS |
| rge | historical no-threshold scale reproduction is within five percent | PASS |
| hierarchy | verifier detects whether P1 sigma/omega supports the RGE interval | PASS |
| hierarchy | actual scalar band overlap with both vector bands is explicitly detected | PASS |
| covariance | threshold covariance is positive semidefinite | PASS |
| covariance | all positive P1 rows enter the collapsed regression Jacobian | PASS |
| doublet | Hellmann-Feynman tuning slope is nonzero | PASS |
| doublet | nearest heavy-doublet gap is positive | PASS |
| PQ | baseline physical domain-wall number remains three | PASS |
| PQ | two charge-two 10_F fields reduce physical N_DW to one | PASS |
| PQ | preferred repair preserves one-loop differential running | PASS |
| PQ | F10 benchmark thresholds lie inside the PS interval | PASS |
| PQ | F10 two-loop benchmark replays close unification | PASS |
| two-site | parent-resolved fixed-point verifier passes | PASS |
| two-site | fixed-point ratio closes below two per mille | PASS |
| covariance | two-site threshold covariance is positive semidefinite | PASS |
| doublet | one-loop light-doublet renormalization condition is reimposed | PASS |
| CW | bosonic hard-mode CW verifier passes | PASS |
| CW | heavy-Yukawa projection is an explicit nonzero-default nuisance | PASS |
