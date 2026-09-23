# Route G2-D — explicitly powered scattering sidebands

A bounded tree-level engineering calculation. The original G2 masses, recoil field, and exact finite internal packet are retained; no measured masses or optimization enter.

## One additional physical assumption

$$g(t)=g[1+\epsilon\cos(\Omega t)],\qquad g=0.05,\quad\epsilon=0.5,\quad\Omega=0.2.$$

The classical pump is homogeneous in ordinary space and on the circle. Its preferred frame is the incoming center-of-momentum frame. It supplies energy, not ordinary or compact momentum. Both zero-VEV quadratic mass spectra stay unchanged at tree level.

$$g(t)=\sum_qg_qe^{-iq\Omega t},\quad g_0=g,\quad g_{\pm1}=g\epsilon/2;\qquad n+l=m+j,\quad E_{\rm out}=E_n+E_l+q\Omega.$$

$$K_{lmq}=\frac{|g_q|^2p_f}{16\pi E_nE_lE_{\rm out}},\qquad p_f=\frac{\sqrt{[E_{\rm out}^2-(m_m+\mu_j)^2][E_{\rm out}^2-(m_m-\mu_j)^2]}}{2E_{\rm out}}.$$

Set K=0 for a closed threshold. The independent radial delta-function integral gives the two-body phase space p_f/(4 pi E_out), with incoming flux normalization 1/(4 E_n E_l). This yields K directly; the background breaks Lorentz invariance, so no invariant conserved-s cross section is assumed.

## Complete first-order rate budget

The exact G2 packet has l=-4,...,4, w=1.2 and theta0=0.7; incoming n=1 and COM p=0.2. Every open (l,m,q) for q=-1,0,+1 is enumerated using |m+alpha|<R E_out and verified against a wider integer interval. JSON retains all joint outputs and rejected candidates.

| q | open channels | weighted K | conditional share | matter work per event |
|---:|---:|---:|---:|---:|
| -1 | 14 | 7.71283296279e-08 | 0.0366294884092 | -0.2 |
| 0 | 25 | 1.85661444734e-06 | 0.88173616241 | 0 |
| 1 | 27 | 1.71892135709e-07 | 0.081634349181 | 0.2 |

Total common-overlap rate coefficient: **2.10563491267e-06**. Mode-conversion share (m != n): **0.873817014233**. Elastic share (m=n,q=0): **0.104094953228**. Driven mode-preserving share (m=n,q != 0): **0.0220880325391**.

These are fractions conditional on counted scattering, not absolute conversion probabilities. The unscattered identity amplitude is excluded. In particular, q != 0 with m=n is not elastic.

$$\overline W_{\rm event}=\frac{\sum_{lmq}|c_l|^2K_{lmq}q\Omega}{\sum_{lmq}|c_l|^2K_{lmq}},\qquad C_W=\sum_{lmq}|c_l|^2K_{lmq}q\Omega.$$

Mean pump work absorbed per scattering: **0.00900097215435** (mass units). Work-rate coefficient C_W: **1.89527612162e-08** (inverse-mass units). C_W is **not power**: an independently specified density/overlap or luminosity is required.

## Target m=0, |j|=1: six lines, not one recoil line

| incoming l | outgoing j | q | outgoing momentum | weighted K |
|---:|---:|---:|---:|---:|
| 0 | 1 | -1 | 0.828215831746 | 2.03863176326e-08 |
| 0 | 1 | 0 | 1.02175457609 | 3.8975944446e-07 |
| 0 | 1 | 1 | 1.20257368551 | 2.77974898862e-08 |
| -2 | -1 | -1 | 1.18924246709 | 6.37966406817e-09 |
| -2 | -1 | 0 | 1.36191317264 | 1.13431642387e-07 |
| -2 | -1 | 1 | 1.52787996799 | 7.72454517115e-09 |

Closest opposite-sign recoil lines: j=+1,q=1 at p=1.20257368551, and j=-1,q=-1 at p=1.18924246709. Gap: **0.0133312184222**. A pump can create near-overlapping outputs even while increasing the rate.

If an independent ideal sideband label q were supplied, the smallest same-q opposite-sign gap would be 0.325306282482. This is a conditional information advantage, not a detector or quantized-pump measurement provided by this calculation.

## Predeclared frequency sensitivity, not optimization

| Omega | total K | conversion share | mean pump work/event | closest opposite-sign momentum gap |
|---:|---:|---:|---:|---:|
| 0.1 | 2.09094796035e-06 | 0.878721396283 | 0.00298367338856 | 0.163086735506 |
| 0.2 | 2.10563491267e-06 | 0.873817014233 | 0.00900097215435 | 0.0133312184222 |
| 0.4 | 2.12712454507e-06 | 0.866964456376 | 0.0303815562826 | 0.0127614291903 |

## Long-time boundary and checks

For a centered top-hat time envelope only, $A_T(\Delta E)\propto\sum_qg_qF_T(\Delta E-q\Omega)$, $F_T(x)=2\sin(xT/2)/x=T\,\mathrm{sinc}(xT/2)$ with sinc defined as sin(x)/x. The script independently integrates this Fourier transform. It does not turn it into a finite-collision probability.

The rate sum uses the long-time/period-averaged limit with Omega*T much greater than one. At finite T the sideband amplitudes interfere. Omega=0 cannot be substituted into an incoherent sum: coincident amplitudes must be combined. A finite 4D collision envelope, timing, and acceptance are absent.

**565/565 checks pass.** Independent phase-space/root integration, Fourier coefficients, zero-drive G2 regression, energy plus pump work, compact/gauge conservation, full enumeration, packet/share normalization, finite-time Fourier transform, and positive instantaneous quartic coefficients are checked.

The analytic quartic minimum is kappa5_min=0.157079632679>0. This proves nonnegativity of that term at every time, not conservation of the matter Hamiltonian.

q beyond +-1 is absent only at first order for this vertex. A homogeneous pump cannot directly shift KK momentum or convert a lone particle's mode at tree level; the recoil field remains necessary. The unmodeled work reservoir, finite-time collisions, quantized pump, radiative corrections, EFT cutoff, and UV completion remain outside the result.

Reproduce: python3 route_g/code/verify_g2_driven.py. No changes to the undriven G2 verifier are needed. JSON records both source hashes and every joint output.
