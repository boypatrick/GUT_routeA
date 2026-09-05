# P54 Spin(10) Yukawa intertwiner audit

explicit Clifford and Yukawa intertwiners verified; normalization translation required

Checks: 20/20.

Selected chirality: -1; SM matter orientation: -1.

For the literal P1 U126 and C Gamma[5] contraction, the selected chiral matter has conjugate SM hypercharges. Standard matter charges require conjugating scalar irrep and spinor chirality together (or an equivalent global embedding conjugation). Earlier geometric overlap keys cannot be assigned physical u/d labels without this dictionary.

| Matrix contraction | Absolute coefficient |
|---|---:|
| raw_majorana_per_sigma | 5.65685424949 |
| factor_for_MR_equal_sigma_f | 0.176776695297 |
| raw_vector_down_doublet | 1.41421356237 |
| factor_for_unit_vector_Dirac_h | 0.707106781187 |
| raw_126_down_doublet | 1.15470053838 |
| 126_down_doublet_if_MR_equal_sigma_f | 0.204124145232 |
| 126_lepton_doublet_if_MR_equal_sigma_f | 0.612372435696 |
| canonical_MR_over_Dirac126_ratio | 4.89897948557 |
| raw_LL_triplet_coefficient_magnitude | 5.65685424949 |
| LL_triplet_if_MR_equal_sigma_f | 1 |

Signed lepton/down relative Clebsch: {'re': -3.0000000000000013, 'im': 1.9783251369615097e-16}.

The raw 1/5! action-card contraction on the actual unit P1 omega is not normalized to MR=sigma f.
The output gives the conversion explicitly. It does not insert a phenomenological convention silently.

Remaining gates:

- Translate the action-card h,f and four VEV/overlap symbols to these canonically normalized matrix intertwiners before fitting.
- Build all actual heavy/light fermion matrices on the P1/P2 scalar background with one fixed normalization and re-evaluate the CW projection.
- The minus-three Clebsch is verified; defining both unit doublet coupling f and MR=sigma f requires a nontrivial conversion of the canonical singlet sigma or of the doublet VEV.
- The absolute LL-triplet coefficient is verified; carry its complex phase with the same fermion/scalar dictionary into the actual type-II source response.
- No global flavor/seesaw fit or two-loop PS Yukawa transport is supplied by this algebraic audit.
