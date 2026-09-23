# Route G2 — candidate circle-uniform probe coupling tests

This is a synthetic EFT normalization test, not a physical-scale calibration, particle fit, or detector qualification. No existing G2 scattering kernel is changed. Every number uses one arbitrary engineering mass unit; none is a measured Higgs, nucleon or target input.

## Candidate interaction and selection rule

$$\mathcal L_{\rm portal}=-H^\dagger H\left[\lambda_\Phi\sum_n|\phi_n|^2+\lambda_X\sum_l|x_l|^2\right],\qquad H^\dagger H=\frac{(v+h)^2}{2}.$$

$$\int_0^{2\pi}R\,d\theta\,u_m^*u_n=\delta_{mn},\quad \Delta M_S^2=\lambda_Sv^2/2,\quad hS^*S:\ -i\lambda_Sv.$$

The normalized uniform-circle overlap is the identity. A different, brane-local probe has matrix u_m*(theta0)u_n(theta0), rank one with off-diagonal entries; this alternative is not adopted. Its matrix is rescaled by 2 pi R only for a dimensionless rank comparison.

A common mass shift preserves the signed-label second difference 2/R^2. A neutral tower remains exactly degenerate under j -> -j, and this uniform scalar probe has identical masses and cross sections for both signs. It cannot serve as a direct sign tag.

If both species' portals are nonzero, an extra tree h-mediated diagonal Phi-X elastic amplitude interferes with the original q=0 elastic amplitude. Old total rates and branch fractions must then be recomputed. The lambdaPhi=0, X-only option avoids this additional tree amplitude, but is not a radiative-stability proof.

## Complex-scalar Higgs widths

$$\Gamma(h\to S_nS_n^*)=\frac{\lambda_S^2v^2}{16\pi m_h}\sqrt{1-\frac{4m_n^2}{m_h^2}},\qquad 2m_n<m_h.$$

Each signed label is one complex field. Its distinguishable particle-antiparticle pair is already included: do not multiply by another two. Different +n and -n fields are separate channels when kinematically open. An independent radial two-body phase-space integration verifies the width normalization, and a neutral light tower checks Gamma_total=Gamma_0+2 Gamma_1. The finite signed-label bound is followed by the exact threshold; no infinite degeneracy factor is inserted.

Synthetic card: R=1, alpha=0.25, bare M^2=0.25, lambda=0.02, v=2, mh=4, mN=0.8, fN=0.3, stated EFT cutoff=8. After the common shift, physical M^2=0.29.

| signed n | physical scalar mass | beta | partial width |
|---:|---:|---:|---:|
| -2 | 1.83098334 | 0.402336923 | 3.20169551e-06 |
| -1 | 0.923309266 | 0.88705975 | 7.0589972e-06 |
| 0 | 0.593717104 | 0.954921463 | 7.59902355e-06 |
| 1 | 1.36106576 | 0.732717544 | 5.83078095e-06 |

Bsum=2.97703568; total width=2.36904972e-05. The synthetic neutral heavy tower has no open modes at mh=4, hence no bound from this partial-width channel. The listed cutoff is an assumed scope condition, not an independently derived UV cutoff.

## Explicit low-energy nucleon matching

$$\mathcal L_{hNN}=-\frac{f_Nm_N}{v}h\bar NN,\qquad\overline{|\mathcal M|^2}=\frac{\lambda_S^2f_N^2m_N^2(4m_N^2-t)}{(m_h^2-t)^2}.$$

$$\frac{d\sigma}{dt}=\frac{\overline{|\mathcal M|^2}}{16\pi\mathcal K(s,m_S^2,m_N^2)},\quad -4p_{\rm COM}^2\leq t\leq0,\quad\sigma_{\rm NR}=\frac{\lambda_S^2f_N^2m_N^4}{4\pi m_h^4(m_S+m_N)^2}.$$

The initial nucleon spin is averaged and the final spin summed. Explicit Dirac-matrix traces verify 4mN^2-t. Direct t integration, a rescaled integration, and an analytic antiderivative agree. This is a fixed synthetic pointlike fN matching ansatz; momentum-dependent form factors and nuclear composition are absent.

| scalar speed in nucleon-rest frame | integrated sigma | relative departure from NR |
|---:|---:|---:|
| 0.3 | 2.31898513e-09 | 0.0172691 |
| 0.1 | 2.355399044e-09 | 0.00183776 |
| 0.03 | 2.359347228e-09 | 0.000164617 |
| 0.01 | 2.359692538e-09 | 1.82832e-05 |
| 0.003 | 2.359731799e-09 | 1.64541e-06 |
| 0.001 | 2.35973525e-09 | 1.82823e-07 |

NR limit=2.359735681e-09 in inverse-square engineering mass units. Tests also cover zero coupling, lambda-squared scaling at fixed physical masses, increasing-mediator suppression, and exact neutral +/-j degeneracy.

## Width versus interaction compatibility

For $B=\sum_{\rm open}\sqrt{1-4m_n^2/m_h^2}>0$ and a hypothetical upper allowance $\Gamma_{\rm allow}$ on this entire tower partial width:

$$\lambda_S^2\leq\frac{16\pi m_h\Gamma_{\rm allow}}{v^2B},\qquad \sigma_{SN}^{\rm NR}\leq\frac{4f_N^2m_N^4\Gamma_{\rm allow}}{v^2m_h^3(m_S+m_N)^2B}.$$

Using only the synthetic allowance Gamma_allow=0.0001 gives lambda^2_max=0.00168844071 and sigma_max=9.96068449e-09. Substitution saturates both algebraic bounds. If B=0, this decay channel supplies no such limit. These are not experimental bounds.

The elimination and lambda-squared tests keep the physical spectrum fixed. Holding bare masses fixed instead changes the common portal mass shift, beta sum and thresholds, making the constraint implicit. Interpreting a partial-width allowance as an invisible-width constraint additionally requires the final states really to be invisible under the relevant selection.

## A necessary X-only fixed-mass matching gate

$$M_{D,\rm bare}^2=25\Lambda^2-\frac{\lambda_Xv^2}{2},\qquad \lambda_X\leq\frac{50\Lambda^2}{v^2}.$$

This condition follows only for the stated action with no X self-quartic or other stabilizing X interaction. At H=Phi=0, the neutral j=0 X mode may be constant, so its gradients and portal terms vanish: V(X)-V(0)=MD_bare^2 |X|^2. A negative coefficient runs to minus infinity as |X| grows. At equality this restricted direction is flat, not strictly isolated.

Synthetic checks use Lambda=0.1,1,10, v=2 and lambdaX at 0,0.5,1,1.1 times the bound. All preserve the frozen physical MD^2=25 Lambda^2; the last factor gives a negative bare coefficient and fails this necessary gate. Large algebraic coupling values are not weak-coupling or EFT recommendations.

This constant-X proof must not be generalized to the charged Phi field: nontrivial circle holonomy can supply an unavoidable covariant-gradient term, so a negative Phi bare mass coefficient alone need not imply the same runaway. The X bound itself does not certify the full vacuum, loop stability, perturbativity, or detector performance.

## A per-pass probability is not a resolution model

$$P(\text{at least one interaction})=1-e^{-n_Td\sigma}.$$

This expression is tested only as an independent dilute-target first-interaction law with synthetic number density and path length. It does not determine momentum resolution, sign-classification error, trigger acceptance, luminosity, or a real detector efficiency. A larger coupling is therefore not itself a readout specification.

**113/113 checks pass.** Existing G2, G2-S and driven verifier hashes are unchanged. No GeV map, physical hardware number, cosmology, dark-matter abundance, experimental exclusion, or UV completion has been added.

Reproduce: python3 route_g/code/verify_g2_probe_coupling.py. JSON retains overlap matrices, all widths, spin traces, velocity checks, hypothetical compatibility bounds and source hashes.
