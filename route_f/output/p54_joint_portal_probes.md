# Two coarse joint-portal engineering probes

Explicit second-stage design, not an optimizer or fit. No default benchmark changes.

Checks: 36/36; new Hessians: 0.

| sigma/omega | min hard radial mass2 | rho_H | rho_Z | Next radial audit only |
|---:|---:|---:|---:|:---:|
| 0.3 | -0.0031774696 | 1.17695117 | 1.13159926 | False |
| 0.5 | 0.0229910277 | 0.538986352 | 1.89903168 | False |

The quartic parameter lambda0 is set to 0.2 (not multiplied by 0.2).
lambda2/lambda4/lambda4p/chi4 are multiplied by 0.01; alpha/beta/chi2 by 0.1.
All other independent scalar couplings are retained; three stationary masses and xi02 are solved anew.

No additional full Hessian evaluation is used: cache preflight plus exact invariant continuation only.
Passing a radial gate would still leave full mixed BFB, global vacua, nonradial stability, gauge/ghost/fermion momentum, loop doublet tuning, Wilson/box, lower matching and finite seesaw open.
- joint_portals_sigma_0.3: One or more independent radial positivity/control gates failed
- joint_portals_sigma_0.5: One or more independent radial positivity/control gates failed
