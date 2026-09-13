# Bounded P54 scalar parameter reselection

**Engineering probes only; none is automatically a viable model or a replacement matching input.**

Checks: 36/36.

| Scalar multiplier | sigma/omega | Lowest hard radial eigenvalue | Max relative radial insertion | Max hard scalar kinetic correction | Pass 0.3 radial-control margin |
|---:|---:|---:|---:|---:|:---:|
| 0.03 | 0.3 | -0.040837813 | 18.86484 | 0.090097461 | False |
| 0.03 | 0.5 | -0.044064027 | 7.9561084 | 0.092964388 | False |
| 0.1 | 0.3 | -0.19859462 | 27.000032 | 0.30032487 | False |
| 0.1 | 0.5 | 0.009713662 | 1.4322992 | 0.30988129 | False |

The exact quartic-action polynomial jets avoid new Hessian evaluations. All accepted points would still require full nonradial, gauge/momentum, light-doublet and threshold verification. No experimental mass was used. Parameter cards, tree classifications and all finite tadpole shifts are retained in JSON.

## Missing

- full nonradial one-loop stability and complete gauge/fermion momentum consistency
- one-loop light-doublet retuning and actual new upper PS saddle/matching
- two-loop gauge running and full thresholds at any reselected point
