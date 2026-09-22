# Route G — Scheme A: common spectrum and internal modes

Created: 2026-09-22. Separate research branch, explicitly requested by the user.

Route G 回到「同一內部結構，產生不同能譜與可見性」的構想。它不把 Spin(10) 的表示維度
當成額外時空，也不把轉換座標誤認為轉換粒子質量。
Route F 的四維 P54 模型與先前計算保留不動；Route G 尚未導出或取代它。

## G2 — Local dynamical conversion (2026-09-23)

G2 adds one neutral complex detector field X and the positive local
potential kappa5 |Phi|² |X|². The detector is localized by a quantum
wavepacket, not a pinned external defect. Its recoil makes the process

\[
\phi_n+X_l\longrightarrow\phi_m+X_{l+n-m}
\]

possible, with one tree amplitude M=-g, g=kappa5/(2 pi R).
Total four-dimensional energy/momentum, compact momentum and both
charges are conserved. The outgoing excitation occupies another existing
mass branch; no coordinate transformation changes a mass eigenvalue.
This action is an explicit minimal candidate, not uniquely derived from
the motivating paper.

On the unchanged G1 quarter-holonomy card (R=1, M5=0.5, alpha=0.25),
add MD=5, g=0.05 and incoming COM momentum 0.2. The component
Phi_1 + X_0 -> Phi_0 + X_1 has outgoing momentum 1.02175458 and
tree cross section 6.27164081e-6 in inverse-square common energy units.
The Phi rest-energy decrease 0.78727421 pays for detector internal
recoil 0.09901951; the remaining 0.68825469 becomes ordinary kinetic
energy. No external power source is present.

The finite detector preparation is localized around theta=0.7 with
coefficients c_l proportional to exp[-l²/(4*1.2²)] exp[-i l theta0],
l=-4,...,4. Different components have different collision energies.
Tracing over recoil removes cross terms for a fixed incoming/outgoing
Phi mode. A common-overlap averaged sigma*v is a rate coefficient,
not an absolute encounter probability; the latter still needs specified
four-dimensional envelopes and geometry.

The result is **tree-level EFT only**. The classical vacuum is stable
and the G1 quadratic tower is unchanged, but interacting loop masses,
counterterms and higher operators have not been calculated. G1's exact
free response is not promoted to an all-orders interacting identity.
No measured particle data were fitted.
The G2 verifier passes **262/262 implementation checks**; the G1
regression still passes 331/331. The tree conversion milestone is
complete within this scope, not a closure of the full physical theory.

## First bounded milestone: G1

The initial prototype is a positive-energy complex scalar on fixed
four-dimensional Minkowski spacetime times one circle. One shared
compactification radius R, bulk mass M5 and background holonomy alpha
determine the entire tower:

\[
m_n^2=M_5^2+(n+\alpha)^2/R^2,\qquad n\in\mathbb Z.
\]

Two probes at opposite points of the same circle have fixed interactions,
including the parallel transporter between them. They see the SAME poles
with different residues. Their normalized relative weights are

\[
p_{\pm,n}=\frac{1\pm\cos[\pi(n+\alpha)]}{2},\qquad
p_{+,n}+p_{-,n}=1.
\]

These are probe weights, not decay probabilities. At a degenerate mass,
residues are summed over the complete eigenspace.
In particular, at alpha=0 a degenerate sine standing wave can be dark
to BOTH probes: the nonzero branch sum does not ensure visibility of
every superposition. The mass pole remains visible through the other
combination.

**The useful distinction:** changing a basis leaves masses and the
correctly transformed response invariant; changing a physical holonomy
or radius changes the tower. A mode dark to one probe still exists.
No measured particle masses are used to select the example parameters.

## Files and reproduction

- [Research plan and rejection gates](ROADMAP.md)
- [G2 local action, recoil and packet derivation](tex/route_g_local_conversion.tex)
- [G2 typeset derivation](output/pdf/route_g_local_conversion.pdf)
- [G2 numerical report](output/g2_local_conversion.md)
- [G2 executable](code/verify_g2_local_conversion.py) and [saved tests](output/g2_local_conversion.json)
- [Action and full short derivation](tex/route_g_scheme_a.tex)
- [Typeset derivation](output/pdf/route_g_scheme_a.pdf)
- [Reproducible numerical checks](code/verify_g1_circle_spectrum.py)
- [Numerical report](output/g1_circle_spectrum.md)
- [Saved spectra, residues, tests and source hash](output/g1_circle_spectrum.json)
- [Sources and interpretation boundary](SOURCES.md)

~~~sh
python3 route_g/code/verify_g1_circle_spectrum.py
python3 route_g/code/verify_g2_local_conversion.py
~~~

G1 uses NumPy; G2 uses NumPy and SciPy.
The four engineering cards are R=1, M5=0.5 with alpha=0, 0.25, 0.5,
and R=2, M5=0.5, alpha=0.25, in common arbitrary energy units.
One spectral computation is checked independently against a link-covariant
finite-difference operator and explicit gauge/basis transformations.
This is not a new AP-E soliton lattice scan.
The independently rerun verifier passes **331/331 implementation checks**.
Halving the grid spacing gives low-spectrum error ratios 3.98--3.99,
consistent with the derived second-order convergence; the truncated
responses satisfy their analytic tail bounds. These are code checks,
not experimental evidence for an extra dimension.

## What is and is not being claimed

| Claim | Status |
|---|---|
| One action produces a correlated tower and probe couplings | Analytically derived; numerical verification in G1. |
| Coordinates alone change the particle masses | Rejected by basis covariance. |
| Darkness in one channel means the mode ceases to exist | Rejected for the two-probe prototype. |
| A physical change of alpha or R can change the spectrum | Demonstrated between the prescribed static cards; symmetries can still relate equal spectra. |
| Every state must be visible in at least one of these two probes | False in a degenerate eigenspace: a joint dark standing wave exists at alpha=0. |
| The background dynamically selects alpha and R | Not yet modeled. |
| A collision physically changes the occupied Phi mode | G2 gives an allowed tree-level process with dynamical detector recoil. G1 alone does not. |
| A complete finite-encounter probability has been numerically computed | Not claimed: cross sections and a scoped common-overlap rate are supplied. |
| Electron, muon, quarks or three families are identified with this tower | Not claimed; the present field is scalar, not a chiral fermion. |
| The original paper's Ei/Ep/En have been derived as these fields | Not claimed. |
| Negative energy, emergent time, superluminal travel, black holes or galaxy curves are explained | Not included. |
| This is a complete superstring compactification or a Spin(10) embedding | Not claimed. |

R here is a physical circle radius, not an identification with the symbol
R or the dimensionality discussion in the source paper. Alpha is a
modeling hypothesis motivated by internal phases, not a measured Ei phase.
The ordinary background spacetime is assumed in G1, not derived as emergent.

## Connection to Route F

The previous Route-F low-order flavor/seesaw tension remains unchanged.
Only a later explicit map from shared geometry and normalized modes to
the P54 fields, couplings and thresholds could connect the two programs.
That map would be a new tested hypothesis, not an inference from similar
pictures. No full matching or physical fit is promoted by G1.

G2 implements the dynamical-detector option. Its heavy finite recoil
bypasses the previously closed one-body emission process by supplying
an incoming collision partner, not by violating the old inequality.
The most useful next physical test is recoil-resolved conversion:
check j-l=n-m together with the predicted outgoing ordinary momentum.
Specified finite-encounter packets are an option when an absolute
probability is needed. An energy-accounted driven background stays a
backup; no drive has been added.
