# Physics-first two-matrix LO flavor screen

no_LO_shape_witness_in_this_bounded_search

**Conditional common-scale flavor screen, not a matched P54 physical fit.**

Comparison scale: 5.0822675e+13 GeV; sigma=MI/g4=8.9080652e+13 GeV.
Nine charged masses and four CKM parameters are checked after one-loop SM target transport. Three PMNS angles and the solar/atmospheric ratio use frozen NuFIT 6.0 normal-ordering inputs without neutrino RG transport.
The charged 10%/30% widths and CKM 10%/0.12-rad widths are declared LO discrepancy filters, not experimental uncertainties or adjustable threshold shifts.

Bounded search objective: 2.3766912; numerical checks 18/18.

- charged_masses_within_declared_LO_widths: False
- three_CKM_sines_within_10_percent: True
- CKM_delta_within_0p12_rad: True
- three_PMNS_sin2_inside_frozen_NuFIT_3sigma_NO: True
- delta_m21_over_delta_m31_inside_3sigma_rectangle: True

Normalized-overlap sigma interval (central atmospheric scale, ||h_D||,||f_D||<=1): [703992701730.2725, 18038702229413.59] GeV.
Baseline sigma belongs to this witness family: False.
Shape-minimum down Yukawa residual: -31.971%, compared with a predeclared 30% width. The formal gate is not relaxed after seeing this near miss.
Even removing norm caps leaves the normalized-overlap supremum sigma=2.0631235e+13 GeV for THIS shape texture.
A failure excludes only this texture/overlap profile under these approximations, not all two-matrix models.

## One bounded fixed-sigma refinement

Exactly one refinement from the best shape start used at most 1000 evaluations. Flavor objective=112.373; joint LO screen pass=False.
Its sigma interval is [13948252662004.67, 89151884191705.23] GeV; the baseline is inside, but the charged/CKM/neutrino gates must also pass.
Largest charged discrepancy: mu +80.285% versus the fixed 10% width.
At actual sigma, MN/MI=[0.0007219339373983908, 0.13726928115378156, 1.3445219878018533]. Sterile threshold/EFT assignment remains uncomputed, including a state above MI if indicated.
No extra starts, wider tolerances, or scalar retuning are used to turn the result green. Bounded failure is not a global no-go.

## Exact algebraic interface

The matrices are Hprime=a h_D and Fprime=d f_D. The common relations are Yd=Hprime+Fprime, Ye=Hprime-3Fprime, Yu=r(Hprime+s Fprime), and Ynu=r(Hprime-3s Fprime).
Choose a=cos(theta)/sqrt(1+r^2), d=sin(theta)/sqrt(1+|r s|^2), b=r a, e=r s d. This gives unit norm exactly and MR=(i 2 sqrt(6)) sigma Fprime/d in the stored common phase convention.
Writing H=||Hprime||_2, F=||Fprime||_2, the cap-compatible interval is F/fcap <= d <= sqrt([1-(1+r^2)H^2/hcap^2]/[1+|r s|^2]). It exists iff (1+r^2)H^2/hcap^2+(1+|r s|^2)F^2/fcap^2 <= 1.
For Q=Ynu Fprime^(-1) Ynu^T with singular values q_i, matching the atmospheric central value gives sigma/d = 174^2 * 10^9 * sqrt(q_3^2-q_1^2)/(2 sqrt(6) sqrt(Delta m31^2)). This is the only neutrino scale adjustment.
Caps on Dirac-unit matrices are screening conventions, not a perturbativity theorem. Raw h/f and Majorana f_M norms, a stricter Majorana-component cap, and a looser Dirac-unit cap sensitivity are all retained in JSON without refitting.
The four normalized overlaps have NOT been obtained from a scalar stationary point or mass eigenvector.

## Still untested

- MU-to-MI Pati-Salam Yukawa running
- neutrino-Yukawa feedback and sequential sterile thresholds
- finite upper/lower gauge/ghost/Yukawa matching, Wilson/box and C_HN seesaw matching
- actual scalar action realization of overlaps and suppression of type II
- full perturbative control of collective Yukawa tensors
- leptonic Dirac CP phase and Majorana phases
- absolute lightest-neutrino mass bounds and cosmology likelihood
- neutrinoless double beta decay likelihood
- baryon asymmetry/leptogenesis
- proton decay and gauge unification in this subtask
- precision electroweak/Higgs observables

## Primary sources

- [Table II MSbar charged masses and gauge inputs at 173.1 GeV; Appendix B provides numerical starts ONLY, originally a different 2HDM fit](https://arxiv.org/html/2109.04050v2)
- [frozen PDG 2024 standard CKM central inputs already encoded in the legacy transport helper](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-ckm-matrix.pdf)
- [NuFIT 6.0 IC24 with SK, normal ordering: three mixing angles and two mass-squared differences; frozen 2024 input, not claimed to be the latest fit](https://www.nu-fit.org/sites/default/files/v60.tbl-parameters.pdf)
- [primary NuFIT 6.0 methodology and conditional ordering choice](https://arxiv.org/abs/2410.05380)
