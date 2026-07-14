# DYN-8: fail-closed falsifiability collection

30/30 mechanical/disclosure checks pass.  The physics main line remains open; this collector derives nothing new; zeta is NOT derived.

## Branch map

- **kinematic_core**: SUSY-agnostic, machine-verified (44+24), untouched by all dynamics results
- **susy_minimal_slice**: EXCLUDED (slice family, not model class): unification x proton decay + obstruction map
- **nonsusy_variant**: PRELIMINARY: PS 210-compatible marginal (4.3e33, Hyper-K window, TRIPLE-TESTED against the computed 210 and 126bar spectra); G_LR alive (tau ~ 1e36); the 210-only LR vacuum has rare sampled local minima near vev ratios 0.6--0.7, so the 45-adjoint route is an alternative rather than a necessity; among the two unpromoted D3 mechanisms, only the scale-decoupled instanton benchmark survives the price card; this is not an exhaustive source classification; the normalized contact content is branch-independent by the invariance theorem (M_* branch-local)
- **string_conditional**: three pricing cards, unpromoted, must not be cited as evidence

## Kill bank (branch-tagged, all re-verified)

- **K1** [kinematic core] fourth sequential chiral family kills the a-priori bound N_fam <= 3 (self-carried-index principle) (`route_E/output/route_e_first_principles.json`; physics-promotable=`true`)
- **K2** [conditional Majorana-card diagnostic] Dirac neutrinos remove this Majorana contact target; within the stated nuisance-prior card, all sampled contact fractions exceed 0.01 (not a data posterior) (`output/audit1/dyn4b_unconditional_zeta.json`; physics-promotable=`false`)
- **K3** [SUSY completion (FIRED)] unification x proton decay excludes the slice family; obstruction map recorded (`output/audit2 + output/audit3 ledgers`; physics-promotable=`true`)
- **K4** [non-SUSY PS chain] tau(p -> e+ pi0) = 4.3e33 vs 2.4e34: dead or inside the Hyper-K discovery window (`output/audit9/dyn9_nonsusy_intermediate.json`; physics-promotable=`true`)
- **K5** [non-SUSY branch] e+ pi0 dominance; G_LR at 1e36 partially probeable; K+ nubar dominance would disfavor the branch (`output/audit9/dyn9_nonsusy_intermediate.json`; physics-promotable=`true`)
- **K6** [SUSY slice (historical diagnostic)] unflavored regression gives a small success fraction, but DYN-7F places M1 in the tau-resolved regime; probabilities and arg-zeta windows are not publishable (`output/audit7/dyn7_leptogenesis_argzeta.json`; physics-promotable=`false`)
- **K7** [hidden-messenger extension (invalid model diagnostic)] historical odd-R ceiling replays arithmetically, but DYN-5V shows that the displayed bilinear has no asserted one-loop anomalous dimension and allows X L H_u (`output/audit5/dyn5_messenger_one_loop.json`; physics-promotable=`false`)
- **K8** [string-conditional] for the fixed archival tower, N2 and N3 lie ABOVE the intermediate gauge scale (39x, 4.8e3x on PS); the instanton benchmark permits this coexistence but does not require it for a generic refitted tower (`route_d/output/d3_instanton_majorana_pricing.json`; physics-promotable=`false`)
- **K9** [string-conditional (invalid numerical diagnostic)] the historical DYN-5 ceilings replay to a finite zero-mode selectivity gap, but DYN-5V invalidates those ceilings, so no numerical Delta S bound is physics-promotable (`route_d/output/d4_stueckelberg_protection.json`; physics-promotable=`false`)
- **K10** [string-conditional] preliminary SUSY-breaking toy scan: G_LR window log10 M_SS in [10.30, 15.55]; PS window empty under the stated one-loop/ESH/degenerate-spectrum assumptions (`route_d/output/d5_susy_breaking_bridge.json`; physics-promotable=`false`)
- **K11** [benchmark bookkeeping] Z_178 diagnostic retired: misses the refreshed card by 4.1e-3 rad vs the 5.4e-5 window (76x) (`output/audit1/dyn4b_unconditional_zeta.json`; physics-promotable=`true`)
- **K4r** [non-SUSY PS chain (refresh)] K4 TRIPLE-TESTED: the Hyper-K window survives the computed 210 spectrum (tree positivity 39/250, tau = 4.25e33 all percentiles) AND the computed 126bar spectrum (45/150 viable, same tau) (`output/audit9/dyn9b1b + dyn9b1c ledgers`; physics-promotable=`true`)
- **K5r** [non-SUSY LR chain (refresh, CORRECTED by the ratio scan)] the 210-only LR vacuum IS a tree-level local minimum in rare coupling regions at vev ratios near 0.6-0.7 (8/8800 sampled; He-Meljanac 1986 reproduced at claim level); the earlier fixed-ratio negative (0/249 at a/p = 0.8) was an artifact of the fixed ratio; the epsilon invariant neither rescues nor kills any sampled point; the adjoint route is an ALTERNATIVE, not a necessity (`output/audit9/dyn9b1d_lr_ratio_scan.json (corrects dyn9b1b/1c)`; physics-promotable=`true`)
- **K6r** [non-SUSY branch (refresh)] K6 is replaced by a branch-tagged unflavored regression: ~6x harder (suppression 0.165, 0/4000 target-band hits without boost, largest scanned hit fraction near x60); no likelihood or thermal kinetics make these viability probabilities, and absence of a gravitino ceiling is not a reheating prediction (`output/audit9/dyn9b3_nonsusy_leptogenesis.json`; physics-promotable=`false`)
- **K8r** [string-conditional (refresh)] K8 remains an unpromoted conditional diagnostic: the fixed archival-kernel ceiling and optional top-like ansatz tension, together with the type-II estimate, disfavor the fixed-kernel renormalizable benchmark, while the leptogenesis comparison lacks a global flavor fit and flavored thermal kinetics; D3 itself remains unpromoted (`output/audit9/dyn9b2_nonsusy_flavor_refit.json`; physics-promotable=`false`)
- **K2r** [conditional kernel-refit diagnostic] K2 is refined by the kernel refit: the theta_23 tension is a frozen-anchor property (exact Majorana absorption at unchanged kernels; a 5% Y_e perturbation absorbs it alone) and contact essentiality is lifted to the perturbation level (cf never below 0.01 over 1000 refits with eps up to 0.3) (`output/audit1/dyn4c_kernel_dirac_refit.json`; physics-promotable=`false`)
- **K11r** [benchmark bookkeeping (refresh)] positive-real rescaling theorem: within the fixed-kernel y>0 family, normalized contact content is exactly invariant and only M_* moves; for complex y, zeta is phase-covariant rather than invariant (`output/audit9/dyn9b2_nonsusy_flavor_refit.json`; physics-promotable=`false`)

## Blocking physics work

- DYN-5 has no valid interacting messenger action and permits X L H_u
- DYN-7 and DYN-9b-3 require branch-local flavored kinetics
- DYN-9b-2 has not performed a global non-SUSY Spin(10) flavor fit
- RE-SC3/4/5 remain unpromoted pricing cards; RE-SC4 inherits an invalid DYN-5 numerical ceiling

Upstream: 304 checks green across 21 source ledgers.

## Checks

- [PASS] every source ledger loads and reports all_pass
- [PASS] no source ledger claims a derived zeta (boundary theorem respected program-wide)
- [PASS] validity guards fail closed: DYN-5 is invalid pending an interacting messenger action and DYN-7/9b-3 require flavored thermal inputs
- [PASS] BRANCH kinematic core: 44/44 first-principles + 24/24 dependency-closure checks green; DYN-9 records the core as SUSY-agnostic
- [PASS] BRANCH SUSY slice EXCLUDED: d=5 lifetime short by 7.7 orders at the unification-compatible point AND zero living points in the joint rescue scan
- [PASS] BRANCH non-SUSY PRELIMINARY LEDGER: G_LR has M_I 1e9.4, M_X 1e16.3 and the central d=6 estimate near 1e36; PS is a 210-compatible marginal central estimate at 4.3e33
- [PASS] BRANCH non-SUSY: perturbativity gate makes the archival contact card slice-local (refit REQUIRED); worst PS factor 379x over 4 pi
- [PASS] BRANCH string-conditional: all three cards unpromoted with explicit negative flags
- [PASS] K1 [kinematic core]: fourth sequential chiral family kills the a-priori bound N_fam <...
- [PASS] K2 [conditional Majorana-card diagnostic]: Dirac neutrinos remove this Majorana contact target; within the ...
- [PASS] K3 [SUSY completion (FIRED)]: unification x proton decay excludes the slice family; obstructio...
- [PASS] K4 [non-SUSY PS chain]: tau(p -> e+ pi0) = 4.3e33 vs 2.4e34: dead or inside the Hyper-K ...
- [PASS] K5 [non-SUSY branch]: e+ pi0 dominance; G_LR at 1e36 partially probeable; K+ nubar dom...
- [PASS] K6 [SUSY slice (historical diagnostic)]: unflavored regression gives a small success fraction, but DYN-7F...
- [PASS] K7 [hidden-messenger extension (invalid model diagnostic)]: historical odd-R ceiling replays arithmetically, but DYN-5V show...
- [PASS] K8 [string-conditional]: for the fixed archival tower, N2 and N3 lie ABOVE the intermedia...
- [PASS] K9 [string-conditional (invalid numerical diagnostic)]: the historical DYN-5 ceilings replay to a finite zero-mode selec...
- [PASS] K10 [string-conditional]: preliminary SUSY-breaking toy scan: G_LR window log10 M_SS in [1...
- [PASS] K11 [benchmark bookkeeping]: Z_178 diagnostic retired: misses the refreshed card by 4.1e-3 ra...
- [PASS] K4r [non-SUSY PS chain (refresh)]: K4 TRIPLE-TESTED: the Hyper-K window survives the computed 210 s...
- [PASS] K5r [non-SUSY LR chain (refresh, CORRECTED by the ratio scan)]: the 210-only LR vacuum IS a tree-level local minimum in rare cou...
- [PASS] K6r [non-SUSY branch (refresh)]: K6 is replaced by a branch-tagged unflavored regression: ~6x har...
- [PASS] K8r [string-conditional (refresh)]: K8 remains an unpromoted conditional diagnostic: the fixed archi...
- [PASS] K2r [conditional kernel-refit diagnostic]: K2 is refined by the kernel refit: the theta_23 tension is a fro...
- [PASS] K11r [benchmark bookkeeping (refresh)]: positive-real rescaling theorem: within the fixed-kernel y>0 fam...
- [PASS] reconstruction paper contains the Falsifiability section with branch map and all eleven criteria K1-K11
- [PASS] reconstruction paper status ledger updated (mechanical evidence + still-open list) and manifest extended with the collector
- [PASS] main paper status ledger updated (deferred -> executed with findings) and Z_178 retirement note added
- [PASS] self-containment discipline: the reconstruction paper uses working aliases only inside the terminology remark and literal file paths (spot check: no 'DYN-' outside texttt/remark context)
- [PASS] 9b-4 REFRESH carried by both papers: the reconstruction paper's branch map records the executed re-derivation (vacuum / source / contact stages), K4 cites three spectrum treatments, K6 is branch-tagged, K8 is explicitly unpromoted, and the invariance theorem replaces the bookkeeping remark; the main paper's status note records the executed re-derivation
