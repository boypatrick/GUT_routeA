# P54PQ-v2 bosonic Coleman-Weinberg matching Hessian

Status: **14/14 checks passed**.

Declared scheme: background-field Landau gauge, MSbar, zero-momentum hard-mode matching at `muU=gU*omega`. This is not advertised as a gauge-independent pole mass.

- hard trace: `290` scalars and `33` vectors;
- EFT trace removed: `38` scalar directions and `12` vectors;
- `gU=0.562602889985`, `muU/omega=0.562602889985`;
- scalar kappa: `-1.825257059514e-02 omega^2`;
- vector kappa: `-1.984080499865e-03 omega^2`;
- bosonic kappa: `-2.023665109501e-02 omega^2`;
- bosonic retuning: `delta xi02,B=-2.023670220706e-02`.

The full real light-doublet Hessian is proportional to `I4` by SM invariance. The numerical second direction and mixed-direction checks give

- scalar diagonal mismatch: `2.291e-13`;
- scalar mixed ratio: `7.695e-14`;
- vector isotropy residual: `2.940e-18`.

## Scale variation

| mu/muU | scalar kappa | vector kappa | total kappa | delta xi02,B |
|---:|---:|---:|---:|---:|
| 0.5 | 1.07814599e-02 | 6.03714528e-03 | 1.68186052e-02 | 1.68186477e-02 |
| 1.0 | -1.82525706e-02 | -1.98408050e-03 | -2.02366511e-02 | -2.02367022e-02 |
| 2.0 | -4.72866011e-02 | -1.00053063e-02 | -5.72919074e-02 | -5.72920521e-02 |

## Matching nuisance

The unfitted fermionic term is not set to zero:

`eta_Y(mu)=h^T Pi_heavy-Y(0;mu) h/omega^2`,

`delta xi02(mu)=[kappa_B(mu)+eta_Y(mu)]/w10`.

The nearest heavy-doublet gap remains a separate norm inequality. Thus the bosonic matching coefficient is now numerical, while the heavy-Yukawa projection remains an explicit P3 matching nuisance.

## Checks

- [x] tree background has 38 soft modes and no tachyon
- [x] hard scalar cluster contains 290 positive modes
- [x] hard vector cluster contains 33 positive modes
- [x] light subspace is one real four-plane
- [x] quartic Hessian first derivative is step stable
- [x] quartic Hessian second derivative is step stable
- [x] soft scalar block has no linear mass splitting
- [x] Frechet trace Hessian matches direct finite difference
- [x] scalar CW curvature agrees on two doublet directions
- [x] scalar CW mixed curvature vanishes
- [x] vector CW Hessian is SM-isotropic
- [x] Landau ghost hard contribution vanishes
- [x] bosonic xi02 counterterm retunes all four light modes
- [x] heavy-Yukawa term is retained as a matching nuisance
