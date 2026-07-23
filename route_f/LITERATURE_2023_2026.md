# 2023--2026 Primary-Literature Matrix

Checked on 2026-07-23.  These papers are inputs for Route-F design, not proof
that the present model inherits their results.

## AP-E15 defect-support boundary

AP-E15 adds no transferred theorem and no new literature-dependent claim.
Its defect-measure decomposition, compact fibre-phase convexity, and
Lebesgue-a.e. weighted-cutoff decay are derived directly.  The AP-E14 sources
remain relevant only as structural context for the Hopf split.

The literature boundary is now sharper: a fixed-map theorem valid
Lebesgue-a.e. neither controls a possibly singular relaxation measure nor
interchanges the shrinking-radius limit with a recovery-sequence limit.
None of the sources below is being cited as proving the required
sequence-uniform `lim_r limsup_j` cutoff decay for `mu`-a.e. points, local
recovery for the present complete-minor relaxation, or isolation of its
degree-one minimizer.

## AP-E14 WZ/Morrey and Hopf-transfer boundary (superseded by AP-E15 for the defect mainline)

AP-E14 uses the Hopf split as an exact in-repository algebraic change of
variables, not as a transfer of a regularity theorem.  Auckly--Kapitanski's
Sobolev Pontryagin--Hopf theory supports the meaning of weak lifts and
invariants.  Guerra--Lamy--Zemas (2026) proves constraint-preserving Hopf
lifting/approximation in `W1,2 x L2`, and their companion work proves a
special global-minimality theorem for the Faddeev--Skyrme Hopf map.  Neither
controls the present `L2` second- and third-complete-minor graph norm or the
relaxed `S3` minimizer's linear Morrey oscillation.

Bousquet--Ponce--Van Schaftingen (2025) gives powerful generic topological
screening and smooth approximation results for ordinary Sobolev maps.  It
does not provide the complete-minor endpoint recovery needed to identify the
present relaxation with the classical integral.  The scoped 2023--2026
search found no primary theorem that removes AP-E14's rank-one vertical--
horizontal defect measure or proves the required `O(r)` target oscillation.

| Source | Relevant result | Transfer boundary |
|---|---|---|
| Auckly and Kapitanski, [The Pontryagin--Hopf invariants for Sobolev maps](https://doi.org/10.1142/S0219199710003752), Commun. Contemp. Math. 12 (2010) 121 | Defines Sobolev Pontryagin--Hopf invariants and develops weak lifting machinery. | Supports the invariant/lift framework, not `L2` complete-minor density or Morrey regularity of this relaxed minimizer. |
| Bousquet, Ponce, and Van Schaftingen, [Generic topological screening and approximation of Sobolev maps](https://arxiv.org/abs/2501.18149) (2025) | Establishes broad topological screening and approximation mechanisms in Sobolev mapping spaces. | Ordinary Sobolev approximation does not control simultaneous strong `L2` convergence of `Du`, `M2`, and `M3`. |
| Guerra, Lamy, and Zemas, [Global Minimality of the Hopf Map in the Faddeev--Skyrme Model with Large Coupling Constant](https://doi.org/10.1007/s00220-026-05702-5), Commun. Math. Phys. 407, 171 (2026) | Proves global minimality/uniqueness in a special Hopf sector. | Different target, functional, homogeneous background, and coupling regime; it is not a regularity theorem for Route-E's relaxed `S3` map. |
| Guerra, Lamy, and Zemas, [Some lifting and approximation properties for maps in `W1,2(B3;S2)`](https://arxiv.org/abs/2605.14507) (2026) | Preserves the exact curvature constraint in smooth Hopf approximation. | Gives first-order base and `L2` gauge control, not the positive shear/sextic endpoint bounds exposed by AP-E14. |

## AP-E13 annular/reverse-Hölder transfer boundary (superseded by AP-E14 for the Hopf test)

Luckhaus's foundational manifold-valued interpolation controls ordinary
first-order Sobolev energy; Mironescu--Van Schaftingen explains why critical
manifold trace extension has separate analytic and topological conditions.
Neither controls the present `L2` third-minor filling cost.  Auckly--Kapitanski
explicitly identified uniformly energy-bounded smooth approximation as an
unresolved point even for the usual Skyrme energy.  Esteban--Müller instead
establish integer degree directly in the finite-energy class.  These sources
support the fail-closed boundary, but the AP-E13 Wess--Zumino cap no-go and
the scale-compatible geodesic-cone estimates are derived in-repository.

The 2026 Guerra--Lamy--Zemas result proves global minimality of the Hopf map
for a special large-coupling Faddeev--Skyrme problem.  Its target, invariant,
and homogeneous setting differ from the present `S3` complete-minor graph
problem, so it cannot be transferred as endpoint density or reverse Hölder.
Their separate 2026 lifting paper is more structurally suggestive: it
characterizes `W1,2` Hopf liftability by exactness of the pulled-back area form
and constructs smooth approximations preserving that constraint.  It still
controls only `W1,2 x L2` in the base-map/gauge variables, not the AP-E13
`L2` complete-minor graph norm.  A useful next experiment is nevertheless to
rewrite the `S3` target field as its `S2` Hopf projection plus the vertical
gauge one-form and test whether the Maurer--Cartan cutoff remainder becomes an
exact/div--curl pairing.

| Source | Relevant result | Transfer boundary |
|---|---|---|
| Luckhaus, [Partial Hölder Continuity for Minima of Certain Energies among Maps into a Riemannian Manifold](https://iumj.org/article/3325/), Indiana Univ. Math. J. 37 (1988) 349 | Supplies the classical annular interpolation architecture for manifold-valued Sobolev minimizers. | Controls first-order Sobolev cost, not Wess--Zumino flux or `L2` second/third complete minors. |
| Mironescu and Van Schaftingen, [Trace theory for Sobolev mappings into a manifold](https://doi.org/10.5802/afst.1675), Ann. Fac. Sci. Toulouse 30 (2021) 281 | Separates analytic and topological trace-extension obstructions. | Does not prove the endpoint complete-minor estimates needed here. |
| Auckly and Kapitanski, [Holonomy and Skyrme's model](https://doi.org/10.1007/s00220-003-0901-x), Commun. Math. Phys. 240 (2003) 97 | Develops the finite-energy Skyrme variational class and flags the missing uniformly energy-bounded smooth approximation. | The standard Skyrme action lacks the present sextic `L2` third-minor term; it is a warning, not a no-go theorem for AP-E13. |
| Esteban and Müller, [Sobolev maps with integer degree and applications to Skyrme's problem](https://doi.org/10.1098/rspa.1992.0014), Proc. R. Soc. A 436 (1992) 197 | Proves integer degree and applies it to existence in Skyrme's problem. | Direct degree theory does not yield strong graph-norm density. |
| Guerra, Lamy, and Zemas, [Global Minimality of the Hopf Map in the Faddeev–Skyrme Model with Large Coupling Constant](https://doi.org/10.1007/s00220-026-05702-5), Commun. Math. Phys. 407, 171 (2026) | Establishes special global minimality and uniqueness modulo symmetries in a large-coupling Hopf sector. | It is an `S3 -> S2` Hopf problem, not the `S3` target endpoint replacement or selected relaxed `B=1` regularity theorem. |
| Guerra, Lamy, and Zemas, [Some lifting and approximation properties for maps in `W1,2(B3;S2)`](https://arxiv.org/abs/2605.14507) (2026) | Characterizes `W1,2` Hopf lifts through exactness of `u*omega_S2` and gives smooth approximation preserving `d eta=u*omega_S2`. | The convergence is first-order for the base map and `L2` for the gauge one-form; it does not control `M2/M3` in the AP-E13 graph norm. |

## AP-E12 transfer boundary

The closest current result is Gmeineder--Kristensen's 2024 smooth partial-
regularity theorem for relaxed strongly quasiconvex `(p,q)` functionals.  Its
range is `q<min(np/(n-1),p+1)`, which becomes `q<3` for `(n,p)=(3,2)`.
AP-E11's continuum density reaches pointwise sixth-order growth, so this
theorem cannot be cited as partial regularity for the present `(2,6)` model.
It also treats Euclidean vector-valued maps, not the `S3` constraint plus a
fixed-degree complete-minor graph relaxation.

Cerf--Mariconda's 2024 review is used only to keep the Lavrentiev boundary
visible: existence in a relaxed class does not by itself identify the smooth
problem.  The endpoint-minor statements instead rely on foundational primary
work by Malý (strong approximation with strict exponent loss), Mucci (`L1`
area/minor approximation), and Giaquinta--Modica--Souček (Cartesian-current
closure).  The scoped 2023--2026 search found no primary theorem proving
fixed-degree `L2` density of all first/second/third minors in dimension three.

| Source | Main result relevant to AP-E12 | Transfer boundary |
|---|---|---|
| Gmeineder and Kristensen, [Quasiconvex Functionals of `(p,q)`-Growth and the Partial Regularity of Relaxed Minimizers](https://doi.org/10.1007/s00205-024-02013-8), Arch. Ration. Mech. Anal. 248, 80 (2024) | Proves smooth partial regularity for relaxed strongly quasiconvex integrals in the essentially maximal range `q<min(np/(n-1),p+1)`. | For `(n,p)=(3,2)` the theorem requires `q<3`, not the present `q=6`; the target manifold and degree graph closure are also absent. |
| Cerf and Mariconda, [The Lavrentiev Phenomenon](https://arxiv.org/abs/2404.02901), Am. Math. Monthly 131 (2024) 881 | Reviews why a relaxed minimizer and smooth approximation problem can have different values and why recovery must be proved. | Conceptual warning only; it supplies no complete-minor endpoint density, sphere-valued recovery, or regularity theorem for AP-E12. |

## AP-E11 transfer boundary

AP-E11 uses the finite-element exterior-calculus principle that cochains,
discrete differentiation, and Hodge products must form one compatible
complex.  The recent shifted-product work of Ptackova and Velho supplies
useful evidence that exact Leibniz identities can be retained in discrete
calculus.  Its constructions are two-dimensional and polygonal, however;
they do not prove the present three-dimensional tetrahedral Skyrme action,
its normalized-affine degree identity, or its Gamma limit.  Those statements
are derived independently in the AP-E11 ledger.

The continuum proof uses classical polyconvex lower semicontinuity and weak
minor compactness.  Briani--Cicalese--Kreutz provides a recent example of a
topology-sensitive discrete/continuum proof architecture, but only for a
two-dimensional `S2` regime.  Bethuel's density theorem likewise warns that
topology constrains Sobolev approximation; it does not supply the still
missing fixed-degree density theorem in the complete `M_2/M_3` graph norm.

| Source | Main result relevant to AP-E11 | Transfer boundary |
|---|---|---|
| Ptackova and Velho, [A simple and complete discrete exterior calculus on general polygonal meshes](https://arxiv.org/abs/2401.15436) (2024) | Builds a DEC framework on general polygonal meshes with discrete differential-form operations. | Structural motivation only: the setting is two-dimensional and does not give the tetrahedral Hodge/cell/degree theorem used here. |
| Ptackova, [Discrete exterior calculus with a discrete Leibniz rule](https://arxiv.org/abs/2504.14275) (2025) | Constructs shifted discrete products satisfying an exact Leibniz identity in its declared complex. | Confirms the design principle, not this three-dimensional Alexander--Whitney/Whitney realization or its continuum limit. |

## AP-E10 transfer boundary

The AP-E10 counterexamples are derived in-repo and do not depend on a
literature no-go theorem.  They identify why standard weak-Jacobian and
div--curl results cannot be invoked: the AP-E7 unshifted one-corner products
are not a compatible discrete differential/cup-product complex.  Briane and
Casado-Diaz give the continuum differential-constraint hypotheses; Arnold,
Falk, and Winther give the discrete subcomplex/cochain-projection repair
principle.  Neither source proves the repaired Skyrme model in advance.

For the finite-`R` barrier, the general discrete-variational literature
motivates a periodic cell formula, but AP-E10 evaluates the model-specific
formula directly: direction-wise cycle balance plus Jensen forces a zero
corrector.  This turns AP-E9's raw nonproportionality warning into an exact
homogenized nonuniversality theorem.  Briani--Cicalese--Kreutz remains the
closest recent topology-sensitive Gamma-convergence architecture, but its
two-dimensional `S2` charge-concentration regime is deliberately not
transferred to the three-dimensional `S3` stencil.

| Source | Main result relevant to Route F | Required response |
|---|---|---|
| Jarkovska, Malinsky, Susic, [The trouble with the minimal renormalizable SO(10) GUT](https://arxiv.org/abs/2304.14227) (2023) | The minimal non-SUSY `45 + 126 + 10_C` model has a low-scale/proton problem and difficulty supporting a perturbative SM-like Higgs doublet. | Do not select this field content by minimality alone; require a full light-doublet and perturbativity gate. |
| Haba, Shimizu, Yamada, [Neutrino Mass in Non-Supersymmetric SO(10) GUT](https://arxiv.org/abs/2304.06263) (2023) | A `54 + 126 + 10_C` two-step branch performs gauge-unification and fermion/neutrino fits, including both mass orderings. | Use as the independent F-54 comparison branch and reproduce its scale/matching assumptions. |
| Djouadi, Fonseca, Ouyang, Raidal, [Non-supersymmetric SO(10) models with Gauge and Yukawa coupling unification](https://arxiv.org/abs/2212.11315), revised/published 2023 | Two-loop gauge/Yukawa running for Pati--Salam and left-right intermediate groups shows strong proton and threshold sensitivity. | Route F must use two-loop running and explicit matching; one-loop central values are triage only. |
| Fu et al., [Testing Realistic SO(10) SUSY GUTs with Proton Decay and Gravitational Waves](https://arxiv.org/abs/2308.05799) (2023) | Demonstrates the appropriate joint scale of analysis: two-loop RGEs, thresholds, flavor, leptogenesis, dark matter, proton decay, and gravitational waves. | Use as a benchmark for scope and uncertainty discipline; do not transfer its SUSY conclusions to the non-SUSY branch. |
| Kaladharan and Saad, [Fermion mass, Axion dark matter, and Leptogenesis in SO(10) GUT](https://arxiv.org/abs/2308.04497) (2023) | Shows how vectorlike fermion/scalar extensions can relieve minimal Yukawa-fit tension and simultaneously affect axion/leptogenesis predictions. | If the Route-F minimal flavor fit fails, pre-register this as an extension rather than adding ad hoc texture parameters. |
| Babu, Di Bari, Fong, Saad, [Leptogenesis in SO(10) with Minimal Yukawa sector](https://arxiv.org/abs/2409.03840) (2024) | A `10 + 120 + overline{126}` Yukawa sector fits masses/mixings, gives `N_2`-dominated leptogenesis, and checks unification/proton limits. | Use as the primary comparison for the F-210 flavor sector; reproduce full RG and density-matrix/flavored assumptions before claiming agreement. |
| Bertucci et al., [Positivity Bounds on Massive Vectors](https://arxiv.org/abs/2402.13327) (2024) | Massive-vector positivity requires nontrivial crossing functions, sum rules, and null constraints. | Replace Route-C constructed PSD residues with amplitude-extracted residues and crossing/dispersion tests. |
| Sato, [Majorana Fermion Zero Modes and Anomalies in String and M Theories](https://arxiv.org/abs/2404.18138) (2024 thesis) | Gives a conceptual nonperturbative zero-mode-parity/anomaly criterion for pointlike solitons and brane configurations. | Use only as an anomaly cross-check; it is not E3/M5 charged-zero-mode cohomology or a Pfaffian calculation. |
| Super-Kamiokande, [Search for proton decay via p to e+ eta and mu+ eta](https://arxiv.org/abs/2409.19633) (2024) | Updated channel-specific null search and lifetime limits illustrate that a theory needs a channel table, not one generic lifetime. | Refresh all experimental inputs and publish branching fractions and correlated theory errors by channel. |
| LEGEND Collaboration, [First Results on the Search for Lepton Number Violating Neutrinoless Double Beta Decay with LEGEND-200](https://arxiv.org/abs/2505.10440) (2025) | No signal; reports a `76Ge` half-life limit and nuclear-matrix-element-dependent effective-mass range. | Treat absence of a signal as a likelihood, not proof of Dirac neutrinos; propagate nuclear matrix element uncertainty. |
| Zheng et al., [A lab-based test of the gravitational redshift with a miniature clock network](https://doi.org/10.1038/s41467-023-40629-8) (2023) | Directly resolves gravitational redshift over a laboratory-scale height difference. | Use clock redshift as an external constraint on any repaired `E_d` carrier; local density alone may not replace metric/proper-time structure. |
| Cafasso et al., [Quantum time and the time-dilation induced interaction transfer mechanism](https://doi.org/10.1103/PhysRevD.110.106014) (2024) | Develops a relational quantum-time mechanism with time-dilation-induced transfer in a finite-dimensional setting. | Keep relational time as a parallel interface; it does not derive Route E or make a local energy density equal to time. |
| Mistele et al., [Indefinitely Flat Circular Velocities and the Baryonic Tully--Fisher Relation from Weak Lensing](https://doi.org/10.3847/2041-8213/ad54b0) (2024) | Weak-lensing data probe the persistence of flat circular velocities beyond resolved rotation curves. | Any gravity repair must fit rotation and lensing jointly; a pure clock-rescaling ansatz is not sufficient. |
| Hamada et al., [Q-Balls in the presence of attractive force](https://doi.org/10.1007/JHEP08(2024)242) (2024) | Extends Q-ball stability analysis when an additional attractive force is present. | Use as a modern stability/gravity interface for AP-E5, not as evidence that the present bubble-to-Q-ball identification is already realized. |
| Bykov and Smilga, [Monopole harmonics on complex projective spaces](https://arxiv.org/abs/2302.11691), SciPost Phys. 15, 195 (2023) | Gives spectra and homogeneous-coordinate eigenfunctions for charged quantum mechanics on projective space; states form `SU(n)` multiplets and connect to twisted Dolbeault/Dirac complexes. | Use to audit the AP-E1 monopole spectrum. Keep all Landau levels unless an LLL/holomorphic projection and its gap are separately derived. |
| Aoki et al., [Eta Invariant of Massive Wilson Dirac Operator and the Index](https://arxiv.org/abs/2501.02873) (2025) | Relates the eta invariant of a massive Wilson--Dirac operator to the continuum index under the stated lattice and continuum-limit hypotheses. | Use as a regulator bridge for AP-E7, not as permission to omit the actual target Yukawa kernel, GW chiral-projector placement, endpoint orientation, or admissibility/locality hypotheses of this model. |
| van den Dungen, [Generalised Dirac--Schrodinger operators and the Callias Theorem](https://arxiv.org/abs/2312.17600), Forum Math. Sigma 13, e11 (2025) | Proves a generalized Callias theorem in which the index of a Dirac--Schrodinger operator is computed on a suitable compact hypersurface. | Use it to justify reducing a uniformly Fredholm same-soliton family to its asymptotic positive-mass eigenbundle. It does not provide the charged-QC2D Yukawa map, its uniform gap, or the family determinant class. |
| Jia and Yi, [Discrete Gauge Anomalies and Instantons](https://arxiv.org/abs/2405.09007), Phys. Rev. D 110, 025011 (2024) | Re-examines the `pi_4(Sp(k))=Z2` anomaly, its instanton zero-mode relation, and the distinction between Euclidean and Lorentzian zero-mode roles. | For a semisimple `Sp(4)` candidate, retain the Witten-parity test as necessary but do not replace a full Dai--Freed/bordism and heavy-threshold calculation by that parity count. |
| Saito and Tachikawa, [Cancelling mod-2 anomalies by Green--Schwarz mechanism with B-mu-nu](https://arxiv.org/abs/2411.09223), SciPost Phys. 19, 017 (2025) | Treats the four-dimensional `Sp(n)` Witten anomaly as the nontrivial character of the relevant spin-bordism class and shows that it cannot be removed by the proposed Green--Schwarz mechanism; its appendix gives a direct low-degree `BSp(n)` bordism analysis. | Use it to audit `Omega_5^Spin(BSp(4))` and the mod-two representation character.  Keep this gauge-bordism calculation separate from the two target-space torsion signs on `S3 x S2`, which still require the specified Dirac--Yukawa mass map and regulator. |
| Davighi and Lohitsiri, [WZW terms without anomalies: generalised symmetries in chiral Lagrangians](https://arxiv.org/abs/2407.20340), SciPost Phys. 17, 168 (2024) | Constructs a weakly coupled UV completion of a specific 4d sigma model on `SU(N) x S2`; a 2-group matches the mixed WZW coefficient `n=n_c X_q`, and the result is tree-level exact.  It also sketches `SU(n_c+1)` semi-simple deep-UV unification, while leaving extra-scalar/lepton decoupling to a future explicit model.  Nonminimal scalar charge leaves a discrete gauge sector and changes normalized flux bookkeeping. | AP-E3 specializes the direct-product intermediate construction to charged two-color matter and obtains `n=2`.  The 2026-07-16 audit now defines the bosonic non-extendible sector and performs the proposed `SU(3)` branching/decoupling exercise.  It finds an exact parity no-go for the original charge-one singlet scalar; the closest charge-two branch decouples but changes the `O(2)` dressing to `O(1)`.  Its covering-`U(1)` normalization is not promoted to faithful `U(2)` until correlated centre-flux bundles and determinant lines are matched. |
| Saito, [Wess-Zumino-Witten Terms of Sp QCD by Bordism Theory](https://arxiv.org/abs/2404.06185), JHEP 10 (2024) 099 | Uses invertible-field-theory/bordism methods to test perturbative and nonperturbative consistency of ungauged and gauged 4d WZW terms on spin manifolds. | Treat ordinary anomaly sums as necessary but not sufficient.  For the special target `S3 x S2`, retain both reduced spin-bordism `Z2` factors: the new audit defines their geometric inverse-image invariants but leaves their APS/Dai--Freed signs to the microscopic regulator. |
| Pai, Akiyama, and Todo, [Two-color lattice QCD in (1+1) dimensions with Grassmann tensor renormalization group](https://arxiv.org/abs/2501.18918) (2025) | Computes number density, chiral condensate, and diquark condensate with a Grassmann tensor-network formulation of two-colour lattice gauge theory. | Use as a recent methodological control for measuring competing condensates without a sign problem.  Its dimension, discretization, and matter content do not validate the present four-dimensional charged model; a dedicated 4d continuum-extrapolated calculation is still required. |
| Suenaga et al., [Mass spectrum of spin-one hadrons in dense two-color QCD](https://arxiv.org/abs/2312.17017), Phys. Rev. D 109, 074031 (2024) | Builds a two-color effective theory respecting Pauli--Guersey `SU(4)` and shows that meson, diquark, vector, and axial-vector directions can all matter dynamically. | Pure two-color QCD is a negative control: `SU(4)/Sp(4)=S5` has no stable `pi_3` Skyrmion.  In the charged model, explicitly measure the Pauli--Guersey/diquark gaps rather than inferring them from gauge-anomaly cancellation. |
| Mazitov and Katanin, [Local magnetic moment formation and Kondo screening in the presence of Hund exchange](https://arxiv.org/abs/2403.20036) (2024) | Finds distinct local-moment, quasiparticle, screening, and Mott regimes in a two-band Hubbard model with Hund exchange. | A finite-`U` plateau is not an exact Hilbert-space projection.  Use it only as the old AP-E3 EFT candidate; the new exactly-two theorem must be stated with exact Gauss constraints. |
| Yao et al., [Lieb-Schultz-Mattis Theorem for 1D Quantum Magnets with Antiunitary Translation and Inversion Symmetries](https://arxiv.org/abs/2307.09843), Phys. Rev. Lett. 133, 136705 (2024) | Sharpens constraints on projective unit-cell representations and symmetric gapped phases in one-dimensional magnets. | Complete-cell Gauss constraints remove literal odd constituents, but arbitrary intercell interactions may still generate projective edge modes.  Restrict the proved boundary theorem to the positive-parent cone until a gapped homotopy is supplied. |
| Zhang et al., [Large-N SU(4) Schwinger boson theory for coupled-dimer antiferromagnets](https://arxiv.org/abs/2409.04627) (2024; revised 2025) | Uses `SU(4)` coherent states of spin-1/2 dimers and a Schwinger-boson large-N expansion to describe coupled-dimer phases and their transition. | Treat this as a boundary check that “paired/dimer” does not automatically mean the diagonal `CP1`, a ferromagnetic triplet, or `k=2`; it is not a literature derivation of the proposed Mott/Hund candidate. |
| Lee et al., [Controllable tunability of a Chern number within the electronic-nuclear spin system in diamond](https://www.nature.com/articles/s41534-023-00732-6), npj Quantum Information 9, 66 (2023) | Demonstrates controllable Chern-number transitions from 0 to 3 in an electronic-nuclear NV-spin system and maps the control problem to an interacting three-qubit system. | Higher or level-two Chern sectors are physically realizable in engineered systems, but existence and tunability do not select `k=2` for Route E.  The microscopic Hamiltonian and exactly-two rule must still be derived. |
| Heeck and Sokhashvili, [Q-balls in polynomial potentials](https://arxiv.org/abs/2211.00021), Phys. Rev. D 107, 016006 (2023) | Studies polynomial Q-balls numerically and analytically and shows how single-field effective potentials can arise from renormalizable multifield models. | Use for AP-E5 radial/profile and scale-separation tests; integrating out heavy fields does not automatically preserve a CP1 orientation sector. |
| Espinosa, Heeck, Sokhashvili, [The Tunneling Potential Approach to Q-Balls](https://arxiv.org/abs/2307.05667), Phys. Rev. D 108, 056019 (2023) | Recasts the Q-ball radial ODE as a field-space variational problem analogous to a three-dimensional bounce and constructs approximate/exact examples. | Use as a profile-solver cross-check after the action is fixed; it does not derive the potential or projective doublet. |
| Cuomo et al., [Numerical tests of the large charge expansion](https://arxiv.org/abs/2305.00499), JHEP 05 (2024) 161 | Lattice measurements in the three-dimensional critical `O(2)` model agree with large-charge superfluid EFT predictions over the tested charges. | Treat fixed-charge EFT as a testable expansion, not a universal derivation of AP-E1 coefficients or a justification for identifying macroscopic Q with Chern level two. |
| Sheu and Shifman, [More on Classical Stability of Hopf-like Solitons of the Toroidal-Twisted type](https://arxiv.org/abs/2605.00757) (2026 preprint) | In two-flavor scalar QED, a low-energy Faddeev--Hopf description and large-size numerical/vorton analysis find locally stable twisted configurations in the studied regime. | This is the closest recent doublet-to-Hopf bridge, but it is a preprint and not a global stability theorem under all perturbations or a derivation of Route-E `O(2)`. |
| Balakrishnan, Dandoloff, Saxena, [Exact Hopfion Vortices in a 3D Heisenberg Ferromagnet](https://arxiv.org/abs/2202.07195), Phys. Lett. A 480, 128975 (2023) | Exact Hopf-labelled static configurations are stabilized because spatial inhomogeneity supplies a length scale that changes Derrick scaling. | Use as evidence that topology alone is insufficient and a stabilizing scale is essential; do not transfer the condensed-matter Hamiltonian to a Lorentz-invariant UV theory. |
| Briani, Cicalese, Kreutz, [Energy Concentration in a Two-Dimensional Magnetic Skyrmion Model](https://arxiv.org/abs/2503.15151), Commun. Math. Phys. 407, 28 (2026) | Proves weak-star compactness of quantized topological-charge measures and a discrete/continuum Gamma limit for a two-dimensional `S2` baby-Skyrme-type scaling. | Use its compactness--liminf--recovery architecture as the closest rigorous modern analogue.  It does not transfer to the three-dimensional `S3` AP-E8 action or its logarithmic edge barrier. |
| Alicandro, Braides, Cicalese, Solci, [Discrete Variational Problems with Interfaces](https://doi.org/10.1017/9781009298766), Cambridge University Press (2023) | Develops discrete variational and discrete-to-continuum methods, including cell-formula and interface architectures for lattice energies. | Use the general homogenization architecture only.  AP-E10 separately enumerates the actual bond graphs and proves the zero corrector; the book does not imply five-/six-tet equivalence. |
| Bartels, Bohnlein, Palus, Sander, [Benchmarking numerical algorithms for harmonic maps into the sphere](https://doi.org/10.1515/jnma-2024-0149), J. Numer. Math. (2025) | Compares nodal and projected sphere-valued finite-element discretizations and their Riemannian solvers; the underlying Dirichlet setting has a controlled discrete-to-continuum theory. | Use its manifold-optimization and interpolation checks, but do not import them across the Skyrme quartic, fixed-degree sector, or topology barrier. |
| Gaster, Loustau, Monsaingeon, [Computing harmonic maps between Riemannian manifolds](https://doi.org/10.4153/S0008414X22000074), Can. J. Math. 75 (2023) 531 | Gives conditions on weighted triangulations, including almost-asymptotically-Laplacian weights, under which discrete harmonic maps approach smooth ones. | Alternative triangulation scans are only a necessary diagnostic.  A theorem also needs shape regularity, weight consistency, and a same-action stability argument. |
| Bartels, Palus, Wang, [Quasi-optimal error estimates for the approximation of stable harmonic maps](https://doi.org/10.1137/22M1524497), SIAM J. Numer. Anal. 61 (2023) 1819 | Gives local finite-element error control for geometrically stable sphere-valued harmonic maps under regularity and stability hypotheses. | Use only after a regulator-stable continuum stationary background and its local stability are independently established; it cannot prove their existence. |
| Cicalese, Reggiani, Solombrino, [From discrete to continuum in the helical XY-model](https://arxiv.org/abs/2412.15994) (2024) | Identifies a critical parameter-growth threshold at which a spin-system Gamma limit changes which chirality transitions survive. | Supports an explicit regime split in `gamma(a),epsilon(a)`; it contains neither the AP-E8 logarithmic barrier nor a three-dimensional degree theorem. |

## Foundational sources used through the 2026-07-20 checkpoint

These sources fall outside the three-year search window but are needed for
the global definitions and the independent moduli-space construction.  They
are not substitutes for recent model-specific evidence.

| Source | Exact role in the checkpoint | Boundary retained |
|---|---|---|
| Malý, [$L^p$-approximation of Jacobians](https://dml.cz/dmlcz/118445), Comment. Math. Univ. Carol. 32 (1991) 659 | Gives strong smooth approximation of Cartesian maps in the all-minor topology after a strict exponent reduction `q<p`. | It does not claim the endpoint `q=p=2`, uniform endpoint bounds, sphere-valued recovery, or fixed trace. |
| Mucci, [A characterization of graphs which can be approximated in area by smooth graphs](https://doi.org/10.1007/s100970000025), JEMS 3 (2001) 1 | Characterizes strong graph/area approximation using `W1,1` and `L1` convergence of every minor, with a slicing criterion. | Its `L1` topology is weaker than the AP-E12 `L2` complete-minor graph norm. |
| Giaquinta, Modica, and Souček, [Cartesian currents and variational problems for mappings into spheres](https://www.numdam.org/item/ASNSP_1989_4_16_3_393_0/), Ann. SNS Pisa 16 (1989) 393 | Supplies the sphere-valued Cartesian-current closure and topology-sensitive variational setting. | Closure language is not an endpoint smooth-density theorem for the present exponents. |
| Davighi, Gripaios, and Randal-Williams, [Differential cohomology and topological actions in physics](https://arxiv.org/abs/2011.05768), Commun. Math. Phys. 396 (2022) 1205 | Supplies the differential-character formulation of Wess--Zumino and related topological actions directly on cycles, without demanding a chosen spacetime extension. | Differential cohomology fixes the bosonic holonomy at fixed integral curvature here; it does not choose the two spin-bordism torsion phases of a fermionic UV regulator. |
| Aldrovandi and Schaposnik, [Non-Abelian vortices in Chern--Simons theories and their induced effective theory](https://arxiv.org/abs/hep-th/0702209), Phys. Rev. D 76 (2007) 045010 | Provides an explicit supersymmetric non-Abelian-vortex mother model with internal projective moduli and a low-energy sigma-model description. | This independently realizes a physical moduli fermion, but it is not the charged two-colour `B=1` soliton and supplies no automatic AP-E3 `O(2)` coefficient line. |
| Witten, [Supersymmetry and Morse theory](https://doi.org/10.4310/jdg/1214437492), J. Differential Geom. 17 (1982) 661 | Identifies supersymmetric quantum mechanics with the differential-form complex and its cohomological ground states. | The ordinary tangent fermion produces the Clifford/Dolbeault module; an extra coefficient line must be independently derived before using `H0(O(2))`. |
| Begun et al., [Study of two color QCD on large lattices](https://arxiv.org/abs/2203.04909), Phys. Rev. D 105, 114505 (2022) | Simulates two-flavour QC2D with staggered fermions on `40^4` and `32^4` lattices at `a=0.048 fm`, focusing on confinement and deconfinement at finite chemical potential. | This is a genuine four-dimensional large-lattice control, but it has neither the extra charged-scalar sector nor a `B=1` gauge--meson--ghost--fermion Hessian. The new `2^4` Route-F calculation must remain labelled a deterministic background-field diagnostic. |
| Witten and Yonekura, [Anomaly Inflow and the eta-Invariant](https://arxiv.org/abs/1909.08775) (2019/2021), and Dai--Freed, [Eta-Invariants and Determinant Lines](https://arxiv.org/abs/hep-th/9405012) (1994) | Express nonperturbative anomaly phases through odd-dimensional eta invariants and determinant lines; after perturbative anomaly cancellation the eta phase is bordism invariant. | A chirality-paired four-dimensional product spectrum is only a reference control. The charged-QC2D torsion signs require a specified microscopic regulator or mapping-torus problem; defect mod-two characters show the available choices but do not select one. |
| Yonekura, [Dai--Freed theorem and topological phases of matter](https://arxiv.org/abs/1607.01873), JHEP 09 (2016) 022 | Gives a practical determinant-line and eta-invariant formulation of fermionic phases with boundaries and families. | AP-E7 fixes the physical sign at zero spectral cut with gapped endpoints.  An arbitrary cut requires its spectral-section/end-point orientation and cannot be inserted into the determinant sign by spectral-flow counting alone. |
| Fukaya, Onogi, and Yamaguchi, [Atiyah--Patodi--Singer index from the domain-wall fermion Dirac operator](https://arxiv.org/abs/1710.03379), Phys. Rev. D 96, 125004 (2017), and Fukaya et al., [Mod-two APS index and domain-wall fermion](https://arxiv.org/abs/2012.03543), Commun. Math. Phys. 393 (2022) 1189 | Realize ordinary and mod-two APS indices through domain-wall fermions and clarify how a finite regulator carries the eta phase. | They motivate the AP-E7 finite-regulator direction but do not define this model's missing target-dependent overlap/Yukawa operator or select its two mapping-torus characters. |
| Berg and Luescher, [Definition and statistical distributions of a topological number in the lattice O(3) sigma-model](https://doi.org/10.1016/0550-3213(81)90568-X), Nucl. Phys. B 190 (1981) 412 | Defines geometric lattice degree using spherical simplices and exposes the role of exceptional/dislocation configurations. | AP-E7 uses the analogous Freudenthal regular-value construction only on the admissible subset.  The unrestricted finite-site configuration space is connected, so a naive optimizer may unwind through the exceptional locus. |
| Ward, [Stable topological Skyrmions on the 2-D lattice](https://arxiv.org/abs/hep-th/9502048), Lett. Math. Phys. 35 (1995) 385 | Demonstrates that a lattice action can be designed so topology-changing configurations carry an energetic obstruction and stable discrete solitons become possible. | AP-E8 adopts the design lesson, not Ward's two-dimensional action.  Its own three-dimensional Freudenthal barrier, finite-grid existence theorem, and continuum caveats are derived independently. |
| Schramm and Svetitsky, [Topology and metastability in the lattice Skyrme model](https://arxiv.org/abs/hep-lat/0008003), Phys. Rev. D 62 (2000) 114020 | Uses spherical tetrahedra and a topology-sensitive lattice construction for Skyrmions, emphasizing exceptional configurations and metastability. | AP-E8 uses normalized-affine degree plus a unique-edge subtracted logarithmic barrier.  It does not claim an exact spherical-volume implementation or inherit continuum/regulator independence from this paper. |
| Briane and Casado-Diaz, [A new div-curl result](https://arxiv.org/abs/1501.02152), J. Math. Pures Appl. 106 (2016) 377 | Proves div--curl and weak-Jacobian results under explicit integrability and differential-compactness hypotheses. | The AP-E7 forward products fail before these hypotheses can be used: an exact three-periodic sequence has a nonzero weak minor defect.  A replacement must supply compatible discrete curl/div identities. |
| Arnold, Falk, Winther, [Finite element exterior calculus](https://arxiv.org/abs/0906.4325), Bull. AMS 47 (2010) 281 | Shows that stable differential-form discretizations are organized by a Hilbert subcomplex and bounded cochain projection. | Use this as the structural blueprint for the next Skyrme action: primal cochains, compatible cup product, and positive Hodge star.  It does not choose the model's weights or prove degree equality. |

## Literature-driven design decisions

1. Keep both F-210 and F-54 until the same two-loop/threshold/proton pipeline
   compares them.
2. Use `10 + 120 + overline{126}` as a motivated flavor candidate, but count
   parameters and demand out-of-fit predictions.
3. Add massive-vector crossing/positivity to Route C only after full amplitudes
   are available.
4. Keep D1/D2 conditional: recent zero-mode consistency results strengthen the
   list of gates but do not construct the missing global geometry.
5. Refresh experimental likelihoods at execution time.  Numerical limits are
   versioned inputs, not permanent theorem constants.
6. For AP-E1, keep four claims separate: projective quotient, monopole/line-
   bundle quantization, fixed-charge Q-ball dynamics, and Hopf stabilization.
   No recent primary source derives the full chain in one theory.
7. Do not equate a second-order fixed-charge rotor with
   `H0(CP1,O(k))`; derive a first-order polarization or show every retained
   coupling is below the computed LLL gap.
8. Do not identify the macroscopic Q-ball charge with Route E's level two.
   Either study the true charge-two sector or introduce an independently
   quantized level-two internal sector.
9. AP-E2 is an exact regression theorem, not literature-dependent physics:
   keep its `30/30` discriminant/Killing/transvectant gate green while retaining
   `physics_promotion_allowed=false`.
10. Retain the finite-`U` AP-E3 Mott/Hund construction as the EFT predecessor,
    not the exact rule.  The new `G_r=N_r-1=0` constrained-cell theorem removes
    literal singleton and charge-transfer sectors, while the electron Pauli
    sign derives the anti-aligned `Q^2=O(2)` branch and signed `k=+2`.  Its
    complete-cell positive-parent boundary theorem does not cover arbitrary
    Haldane-like intercell phases.
11. The charged two-color mixed-WZW branch is now an anomaly-consistent
    intermediate UV candidate with `n=2` and `k=nB=+2` for `B=1`.  Its
    bosonic non-extendible action now has a differential-character definition,
    but its two reduced spin-bordism signs, compact-monopole regime, and
    microscopic strong phase remain open.
12. AP-E4's canonical Spin-c mathematics is complete, including the full
    partner tower and gap.  Keep the ordinary-spin tangent result as a
    mandatory negative control: `T=O(2)` alone gives two zero modes; three
    requires the canonical half-canonical shift or ordinary twist `O(3)`.
    Neither operator count is a family theorem without a fermionic UV origin,
    anomaly cancellation, and a Route-E portal.
13. Treat the pure-`SU(3)` deep-UV proposal as refuted for the original field
    normalization: `[SU(2) x U(1)]/Z2` representation parity forbids a
    charge-one colour singlet.  The explicit charge-two scalar variant is a
    useful decoupling control, not an equivalent completion, because its
    minimal baryon dressing is `O(1)` rather than `O(2)`.
14. Keep “physical tangent fermion” and “three twisted ground states” as two
    separate assertions.  The moduli-space `N=2` construction proves the
    first; the second requires an independently derived `O(2)` coefficient
    line and a proof that it is pulled back from AP-E3 on the same soliton.
15. The nonlinear charged-two-colour program is an auditable EFT proxy.  Its
    positive vacuum/radial blocks motivate a dedicated gauge-theory
    calculation but cannot be cited as lattice-QC2D or full-Hessian closure.
16. For the charge-two `SU(3)` near-miss, keep the nonminimal-charge WZW
    normalization conditional: the faithful group is `U(2)`, not a direct
    product, so the `U(1)` flux is correlated with the `SU(2)` centre class.
    Prove the global bundle/determinant-line match before calling its level an
    all-bundle theorem.

17. Treat the new '2^4' charged-QC2D matrix calculation as a background-field
    regression only. Genuine large four-dimensional QC2D ensembles exist,
    but neither the 2022 large-lattice work nor the recent 1+1-dimensional
    tensor-network calculation contains the added charged-scalar sector or
    the requested coupled 'B=1' Hessian.
18. Prefer the simply connected 'Sp(4)=Spin(5)' cover for the next
    semisimple audit, with direct 'SU(2)c x SU(2)H' as the minimal backup.
    Keep two complex scalar '4' copies explicit. Do not quotient by the
    centre, and do not call the tree-level branch an all-scale completion
    before its heavy eta threshold, gauge bordism, and radiative protection
    are computed.
19. Interpret the four product/defect mod-two character values as a proof of
    regulator non-uniqueness. A paired four-dimensional product spectrum is
    not the odd-dimensional mapping-torus eta invariant of the charged UV
    theory.
20. Keep three same-soliton assertions independent: a tangent SQM,
    'det Ind D=O(+2)', and gauge-basic WZW transgression. Fixed-polarization
    CPT further requires the anti-canonical 'O(-4)' sector. The conditional
    Callias theorem and differential fiber integration do not supply the
    missing physical Yukawa family, equivariant descent, or CPT regulator.
21. Construct no degree-one Route-E portal until at least one complete
    pre-portal lane closes. The current integrated gate evaluates all three
    lane booleans to false.
22. For the simply connected compact group written physically as
    `Sp(4)=USp(4)` (mathematical `Sp(2)`), keep the exact result
    `Omega_5^Spin(BSp(2))=Z2` separate from the two target-space torsion
    characters.  The gauge-bordism generator is detected by a fundamental
    `4`, whereas `5` and `10` have even instanton index; this does not select
    either reduced `Omega_4^Spin(S3 x S2)` sign.
23. Do not rename Pauli--Villars moment cancellation, gamma-five reality, or
    local theta-angle periodicity as a complete Euclidean/APS regulator.
    Closure requires the regulator kernel and mass interpolation, boundary
    eta prescription, local counterterm scheme, regulated bosonic/ghost
    sectors, and the full heavy determinant including gravitational phases.
24. Distinguish a continuum radial stationary point from a stationary point
    of its sampled cubic action.  A sparse analytic second derivative can be
    algebraically exact while still being an off-shell curvature operator;
    topology, projected gradient, spacing, and volume limits must all close
    before its low eigenvalues are called a physical soliton Hessian.
25. For the present hedgehog, use the actual localized orbit
    `SU(2)/Z2=SO(3)`, not an assumed `CP1`.  A spatially uniform scalar-vacuum
    `CP1` has norm proportional to volume and is not a localized modulus.
    Therefore reformulating the target on the `SO(3)`/`SU(2)` collective
    coordinate and computing its torsion/Finkelstein--Rubinstein line is now
    preferable to importing `O(2)` from a spectator family.
26. Retain the two-band obstruction as a UV-design test.  On
    `S2_infinity x CP1`, a nondegenerate two-level mass in a trivial rank-two
    bundle has `c1(E_+)^2=0`, so it cannot realize `c1=x+2y`.  Any future
    construction must use at least three bands, a nontrivial ambient bundle,
    or an infinite-dimensional operator family and then prove its uniform
    Fredholm gap and regulator threshold.
27. In APS determinant-sign formulae, use the physical zero cut with gapped
    endpoints, or include the full spectral-section/end-point orientation.
    `SF_alpha` by itself depends on `alpha` and is not an invariant physical
    sign.
28. Keep the AP-E7 regulator claim at the determinant-prescription level
    unless every half-integral power is realized by a valid local first-order
    or Pfaffian PV system.  Likewise, do not infer fixed-background Kramers
    symmetry for an adjoint Higgs mass from a scalar-mass gamma-matrix test.
29. The finite-site no-sector theorem is now a design constraint: naive
    fixed-boundary `(S3)^M` has one connected component.  A future exact
    lattice `B` sector must declare an admissibility barrier or a
    topology-preserving spherical-volume action and demonstrate regulator-
    independent continuum observables.
30. On the physical `SO(3)` orbit, distinguish classification from selection.
    `H^2(SO(3);Z)=Z2` gives two line candidates and a topological Real lift;
    only the actual Yukawa mapping-torus mod-two index and microscopic CPT
    regulator can select the Pfaffian line.
31. The explicit rank-three `O(1,2)` mass family is the minimal bundle-level
    solution of the `c1^2` obstruction.  Keep it mathematical until a localized
    normalizable modulus, charged-two-colour three-band Yukawa representation,
    and same-operator uniform Fredholm gap are derived.
32. AP-E8 now implements the admissibility alternative required by item 29.
    The finite-grid sector and no-pin stationary-family statements are true,
    but the new barrier occupies roughly one quarter of the displayed coarse-
    grid energy and has strong gamma dependence.  Treat its compact-supported
    `O(a^2)` regression as a consistency check on fixed smooth fields only.
    Require equicoercivity or a Gamma-limit/counterexample, alternative
    triangulations, and joint `(a,L,gamma,epsilon)` extrapolation before any
    continuum or regulator-independent interpretation.
33. AP-E9 closes only the first layer of item 32.  On a fixed box the base
    coordinate Dirichlet term proves strong-`L2` equicoercivity.  With
    `d=1-epsilon`, `w=gamma*a`, and `R=gamma*a^2/d^2`, uniform relative
    edge interiority requires `inf w>0`, whereas smooth barrier disappearance
    requires `R->0`; for fixed epsilon and `gamma=c*a^(-p)` this gives the
    unique window `1<=p<2`.  The barrier alone still cannot close degree when
    `R->0`; for `inf R>0` it gives `W1,4` degree compactness.  At critical
    `p=2`, the raw smooth-sampling quartics are mesh-dependent, but the actual
    checkerboard Gamma density requires an unevaluated periodic homogenized
    cell problem and must not be identified with the raw `Q_T`.
    The 23-case production scan also fails the predeclared volume-radius,
    barrier-slope, cross-triangulation, and `p=0` versus `p=1` regulator gates.
    Keep the same-action Hessian and determinant variation embargoed until a
    discrete Skyrme-minor Gamma theorem, the required barrier cell formula,
    and a regulator-stable background are all available.
34. AP-E10 disproves compensated compactness for the old one-corner forward
    minor stencil and proves finite-`R` five-/six-tetrahedron cell
    inequivalence.  Do not use continuum div--curl literature to repair a
    product which is not a compatible discrete differential algebra.
35. AP-E11 replaces that stencil by a compatible complete-minor cochain
    action and proves its mesh-independent cell formula and Gamma convergence
    to the fixed-degree lower-semicontinuous relaxation.  Keep equality with
    the unrelaxed classical functional as an independent graph-density gate.
36. AP-E12 proves graph-norm continuity of degree, exact strong
    quasiconvexity, and a no-go for the direct radial Malý dipole mechanism.
    Do not infer endpoint density from exponent-losing Cartesian
    approximation, and do not import the 2024 relaxed partial-regularity
    theorem across its `q<3` boundary to this `(2,6)` sphere-constrained
    problem.  The next analytic object is a target-valued annular replacement
    with `L2` minor equiintegrability or a model-specific reverse-Hölder gain.
37. AP-E13 disproves that annular replacement under shell-energy control
    alone: a shrinking geodesic-cap trace has vanishing shell energy but a
    non-vanishing Wess--Zumino-forced cubic-minor filling cost.  Use the proved
    geodesic-cone replacement only under linear-in-radius target oscillation.
    Exact strong quasiconvexity does not bypass the blocker because its
    highest cutoff remainder requires that oscillation or prior higher
    integrability.  Seek WZ/Morrey decay or a Noether/Maurer--Cartan
    compensation before any Hessian construction.
38. AP-E14 closes fixed-map WZ-flux decay by `L2` absolute continuity, but
    disproves both its sufficiency for linear Morrey oscillation and a
    universal exact div--curl collapse of the full Hopf cutoff remainder.
    Seek an inner-variation formula for the lower-semicontinuous relaxation,
    retain its stress defect measure explicitly, and eliminate the rank-one
    vertical--horizontal defect before attempting epsilon regularity.  Keep
    Hessian, determinant, and portal closed.

The scoped recent search did not identify a 2023--2026 primary paper that
already supplies the exact global `Spin(10)` F-theory geometry, quantized flux,
massless hypercharge, exotic lifting, and E3/M5 Pfaffian needed here.  Route F
should therefore retain the relevant foundational literature rather than use a
recent but only adjacent zero-mode paper as a substitute.

Likewise, no 2023--2026 primary paper found in the AP-E1 search derives in one
model the chain `4D doublet -> CP1 -> O(2) -> Q-ball -> stable Hopf soliton ->
three chiral families`.  Each arrow remains an independent theorem or gate in
the roadmap.

The 2026-07-15 follow-up changes the intermediate conclusion but not the final
promotion boundary.  The repository now proves the constrained-cell chain
`G_1=G_2=0 -> exactly two -> anti-aligned Hund triplet -> signed k=+2 -> O(2)`
and constructs an anomaly-consistent mixed-WZW intermediate candidate.  It
also solves the canonical AP-E4 spectrum.  No checked primary source or
in-repository theorem yet supplies the remaining chain `nonperturbative strong
phase + all-scale UV + physical Spin-c fermion + anomaly cancellation +
degree-one Route-E portal -> three chiral families`.

The 2026-07-17 AP-E8 follow-up changes the finite-grid conclusion but still
does not alter that promotion boundary.  Seven unanchored admissible
stationary representatives and an exact fixed-grid minimizer theorem now
replace AP-E7's finite-grid negative result.  No 2023--2026 primary source in
the scoped matrix supplies the missing equicoercive continuum limit, physical
charged-QC2D determinant Hessian, or Route-E portal for this action.
