# Fixed-VEV fermionic tadpole and doublet interface

Status: executable synthetic-matrix regression. No flavor fit or actual Spin(10) Clebsch assignment is claimed.

For $M_R=\sigma f_M$, $X=M_R^\dagger M_R$, at fixed renormalization scale $\mu$:

The local code variable `f` denotes $f_M$, not the unit-bidoublet $f_D$. The independently verified canonical translation is $f_M=2\sqrt6 f_D$. The arbitrary synthetic $Y_a$ are not assigned actual Spin(10) Clebsch tensors.

$$V_F(0)=-\frac{\mathrm{Tr}[X^2(\log(X/\mu^2)-3/2)]}{32\pi^2},\quad \frac{dX}{d\sigma}=\frac{2X}{\sigma}.$$

Differentiating the spectral trace gives $t_{\sigma,F}=-\mathrm{Tr}[X^2(\log(X/\mu^2)-1)]/(8\pi^2\sigma)$.

The invariant finite counterterm has radial potential $-\delta\nu_F^2\sigma_{\rm field}^2/2$, hence $\delta\nu_F^2=t_{\sigma,F}/\sigma$. Its canonical doublet Hessian is $-\delta\nu_F^2P_{126}/2$.

$$\Delta M_{D,F}^{2,\mathrm{fixed\ VEV}}=\Pi_F-\frac{\delta\nu_F^2}{2}P_{126}.$$

The same counterterm value is held fixed during the finite-difference field variation. The scale is not varied with the radial field.

Synthetic example: t_sigma,F=8.85324236242e-05, delta_nu2,F=0.000699861056318 in omega=1 units.

Three nontrivial complex 3x3 families, including a degenerate Takagi pair, test the radial derivative, finite-counterterm cancellation, all eight real doublet directions, family unitaries, and copy unitaries with P126 rotated consistently.

Actual bosonic matrix plus synthetic fermion assembly: delta_xi02=-0.0202318898299; heavy gap=0.0302613550298; Schur residual=7.898e-16. This is an assembly regression, not a physical fit point.

Physical fit gates: historical PS census/matching failed; repaired logarithmic identities pass; absolute spinor normalization is independently checked. Finite PS thresholds, complete P54 PS Yukawa betas, common four-copy phases, and type-I plus type-II matching remain open.

Checks: **11/11**. Full matrices and residuals are in `p54_fermion_tadpole.json`.
