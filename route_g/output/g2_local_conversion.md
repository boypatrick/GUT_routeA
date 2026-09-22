# Route G G2 — local dynamical mode conversion

A bounded tree-level engineering calculation, not an empirical particle model. No measured masses or fitted parameters enter.

## Shared action and interaction

Add one positive-kinetic complex bulk detector field X, neutral under Phi's external U(1), and the nonnegative potential +kappa5 |Phi|^2 |X|^2. Both VEVs are zero, so the G1 quadratic masses remain unchanged at tree level. The local 5D interaction is not a pinned defect or an arbitrarily selected projector.

$$m_n^2=M_5^2+(n+\alpha)^2/R^2,\quad \mu_l^2=M_D^2+l^2/R^2,\quad g=\kappa_5/(2\pi R).$$

$$\Phi_n+X_l\longrightarrow\Phi_m+X_j,\qquad \mathcal M=-g\,\delta_{n+l,m+j}.$$

The action conserves both complex fields' net U(1) charges, four-momentum, and total compact momentum. This channel has one particle of each species before and after; total particle-plus-antiparticle number is not an exact symmetry of the full field theory. Changing Phi's internal mode is balanced by detector recoil.

One declared card: R=1, M5=0.5, alpha=0.25, MD=5, g=0.05, initial COM momentum p=0.2; n=1 and target m=0. All quantities use a common arbitrary mass unit. kappa5 has dimension -1, g is dimensionless, sigma and K have dimension -2.

$$p_f=\frac{\sqrt{[s-(m_m+\mu_j)^2][s-(m_m-\mu_j)^2]}}{2\sqrt{s}},\quad \sigma=\frac{g^2}{16\pi s}\frac{p_f}{p_i},\quad K=\sigma\left(\frac{p_i}{E_\Phi}+\frac{p_i}{E_X}\right).$$

Set sigma=K=0 when the final threshold is closed. These expressions are for distinct scalars, with no identical-particle factor.

## Explicit nonzero process

| Quantity | Engineering value |
|---|---:|
| initial Phi mass m1 | 1.34629120178 |
| final Phi mass m0 | 0.559016994375 |
| initial X mass mu0 | 5 |
| final X mass mu1 | 5.09901951359 |
| sqrt(s) | 6.36506416016 |
| final momentum | 1.02175457609 |
| sigma | 6.27164081328e-06 |
| K = sigma v_rel | 1.17224310906e-06 |
| Phi energy change | -0.196384509024 |
| detector energy change | 0.196384509024 |

This is genuine scattering from one Phi tree-level mass eigenstate to another; it is not a passive basis change. The same tree spectrum supplies both states. Interacting loop pole shifts are not computed. Phi's lost energy is detector recoil/energy, not missing energy.

The inverse channel Phi0 + X1 -> Phi1 + X0 at the *same small p=0.2* has sqrt(s)=5.69665743322, below threshold 6.34629120178, so it is closed. At the same s as the forward reaction, inverse detailed balance is satisfied. These are different energy preparations, not a violation of reversibility.

For l=0 the analytic rest-threshold condition is MD > [((n-m)/R)^2-(mn-mm)^2]/[2(mn-mm)] = 0.241465628349; the chosen MD=5 is well above it. Direct threshold checks at 0.9, 1.0 and 1.1 times this critical mass agree. Heavy recoil costs sqrt(MD^2+((n-m)/R)^2)-MD, which approaches (n-m)^2/(2 MD R^2), not the massless emitted-mode cost |n-m|/R.

As p_initial approaches zero, the exothermic cross section scales as 1/p_initial, but sigma*v_rel has the finite limit 1.16914471395e-06. The smallest tested momentum 0.0002 agrees to relative error 2.79225e-09. A cross section is not a probability.

## A finite localized dynamical detector state

$$c_l=Z^{-1/2}\exp[-l^2/(4w^2)-il\theta_0],\quad l=-4,\ldots,4,\quad w=1.2,\quad\theta_0=0.7.$$

This is an exactly normalized finite prepared state; no approximation to an infinite Gaussian is claimed. Each component has the same 3-momentum magnitude but different s_l. The free packet spreads. Its canonical angular one-particle wavefunction is not a covariant relativistic position density.

| time | norm | mean detector energy | magnitude of first circular moment | mean angle |
|---:|---:|---:|---:|---:|
| 0 | 1 | 5.14217672149 | 0.916405360076 | 0.7 |
| 1 | 1 | 5.14217672149 | 0.89416520353 | 0.7 |
| 5 | 1 | 5.14217672149 | 0.480418720087 | 0.7 |

For the fixed target m=0, recoil j=l+1 distinguishes incoming components. After tracing recoil, the inclusive rate kernel is diagonal in l. A pure localized packet and its momentum-dephased mixture have equal inclusive target rates. The illustrative rate matrices retain different off-diagonal entries, but do not determine physical recoil coherence after unobserved momenta or timing are traced. Translating this packet cannot change that rate against a delocalized incident Phi_n mode.

$$K_{\rm packet}=\sum_l|c_l|^2K_l,\qquad A_{jl}=\sqrt{K_l}\,\delta_{j,l+1},\qquad \mathrm{tr}(A\rho A^\dagger)=\sum_lK_l\rho_{ll}.$$

This is a **common-overlap dilute rate coefficient**, not a universal wavepacket cross section, scattering probability, or a G1 source-visibility weight. The rate map A is not a single-energy S-matrix: the components have different s_l. Conditional incoming and outgoing compact momenta must both be weighted by the target event rate.

| initial l | probability weight | sqrt(s_l) | target K_l (m=0) | all open final m | all-channel K_l |
|---:|---:|---:|---:|---|---:|
| -4 | 0.0012853809 | 7.76731271 | 1.190752033e-06 | -2, -1, 0, 1 | 3.880268915e-06 |
| -3 | 0.014608604 | 7.195446617 | 1.30890424e-06 | -2, -1, 0, 1 | 3.782260199e-06 |
| -2 | 0.082907187 | 6.749943193 | 1.368176232e-06 | -1, 0, 1 | 2.876737317e-06 |
| -1 | 0.23495369 | 6.464006088 | 1.327634702e-06 | -1, 0, 1 | 2.561761111e-06 |
| 0 | 0.33249028 | 6.36506416 | 1.172243109e-06 | -1, 0, 1 | 1.807486204e-06 |
| 1 | 0.23495369 | 6.464006088 | 9.358529204e-07 | 0, 1 | 1.157416715e-06 |
| 2 | 0.082907187 | 6.749943193 | 6.822115734e-07 | 0, 1 | 8.831313177e-07 |
| 3 | 0.014608604 | 7.195446617 | 4.614465027e-07 | 0, 1 | 6.355343718e-07 |
| 4 | 0.0012853809 | 7.76731271 | 2.909258195e-07 | 0, 1 | 4.378003781e-07 |

Weighted target coefficient: **1.11933296351e-06**. All-channel coefficient: **1.85661444734e-06**. Target fraction of all scattering rates: **0.602889288681**. These fractions depend on this prepared state and common-overlap prescription; the chosen m=0 channel is not exclusive. The total includes elastic scattering m=n, not the probability of remaining in that KK sector: the latter also includes the unscattered identity amplitude.

Conditional mean incoming detector l=-0.221167768855; conditional mean outgoing j=0.778832231145; shift=1=n-m. Comparing outgoing mean with the *unconditioned* incident mean would incorrectly mix recoil with event-selection bias.

Every open channel is enumerated by the rigorous finite candidate bound |m+alpha| < R sqrt(s_l), followed by the exact two-body threshold. JSON also records rejected candidates and every open channel's energies, cross section, independent radial phase-space check, and inverse rate.

## Verification and scope

**262/262 numerical/code checks passed.** This certifies this implementation, not the truth of the model as a particle theory.

Checks cover normalized circle overlap and forbidden channels; independent radial-root/phase-space versus Kallen kinematics; complete open-channel enumeration; energy and compact-momentum conservation; inverse detailed balance at common s; zero coupling; large-gauge relabelling; exact finite-packet normalization, free spreading and energy; inclusive recoil trace, dephasing, translation and conditional momentum shift.

The positivity statement is analytic: all three potential coefficients are nonnegative. Sampled potential values are regression checks, not a numerical proof for all field amplitudes.

Limits: 5D EFT and tree amplitudes only; the EFT cutoff must exceed the prepared scattering energies, and no loop matching/naturalness assertion is made. Radius/holonomy stabilization, localized-state preparation, spin/chirality, three families, the full CERN-paper claims, and a Route-F or Standard Model bridge remain open. This does not derive the G1 paired probe operator from the new detector field.

Reproduce: python3 route_g/code/verify_g2_local_conversion.py. No plot, lattice scan, per-particle fitting, or external work source is used. JSON retains the source SHA-256 and every check.

Normalization references: [PDG Kinematics, equations 49.27--49.33](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-kinematics.pdf); [Tong QFT, sections 3.4 and 3.6](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S3.html). Full action and packet derivation: ../tex/route_g_local_conversion.tex.
