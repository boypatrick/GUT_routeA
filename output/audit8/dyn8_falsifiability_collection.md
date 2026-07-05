# DYN-8: falsifiability collection

22/22 checks pass.  Closes the dynamics main line; derives nothing new; zeta NOT derived.

## Branch map

- **kinematic_core**: SUSY-agnostic, machine-verified (44+24), untouched by all dynamics results
- **susy_minimal_slice**: EXCLUDED (slice family, not model class): unification x proton decay + obstruction map
- **nonsusy_variant**: ALIVE: G_LR alive (tau ~ 1e36), PS 210-compatible marginal (4.3e33, Hyper-K window); flavor refit REQUIRED (contact card slice-local)
- **string_conditional**: three pricing cards, unpromoted, must not be cited as evidence

## Kill bank (branch-tagged, all re-verified)

- **K1** [kinematic core] fourth sequential chiral family kills the a-priori bound N_fam <= 3 (self-carried-index principle) (`route_E/output/route_e_first_principles.json`)
- **K2** [kinematic core target] Dirac neutrinos remove the Majorana contact target; essentiality posterior-wide (`output/audit1/dyn4b_unconditional_zeta.json`)
- **K3** [SUSY completion (FIRED)] unification x proton decay excludes the slice family; obstruction map recorded (`output/audit2 + output/audit3 ledgers`)
- **K4** [non-SUSY PS chain] tau(p -> e+ pi0) = 4.3e33 vs 2.4e34: dead or inside the Hyper-K discovery window (`output/audit9/dyn9_nonsusy_intermediate.json`)
- **K5** [non-SUSY branch] e+ pi0 dominance; G_LR at 1e36 partially probeable; K+ nubar dominance would disfavor the branch (`output/audit9/dyn9_nonsusy_intermediate.json`)
- **K6** [SUSY slice (soft)] thermal unflavored leptogenesis fails by ~1.5 orders, wrong central sign; x3-10 boost -> 34-45% (`output/audit7/dyn7_leptogenesis_argzeta.json`)
- **K7** [hidden-messenger extension] odd-R spurion ceiling eps < 4.4e-4 (refreshed windows); light-sector silence structural (`output/audit5/dyn5_messenger_one_loop.json`)
- **K8** [string-conditional] instanton source requires N2, N3 ABOVE the intermediate gauge scale (39x, 4.8e3x on PS): coexistence correlation falsifiable (`route_d/output/d3_instanton_majorana_pricing.json`)
- **K9** [string-conditional] zero-mode selectivity gap Delta S >= 6.6 (refreshed) / 3.2 (loose) at M_s = 2e16 (`route_d/output/d4_stueckelberg_protection.json`)
- **K10** [string-conditional] SUSY-breaking bridge: G_LR window log10 M_SS in [10.30, 15.55]; PS bridge EMPTY (structural) (`route_d/output/d5_susy_breaking_bridge.json`)
- **K11** [benchmark bookkeeping] Z_178 diagnostic retired: misses the refreshed card by 4.1e-3 rad vs the 5.4e-5 window (76x) (`output/audit1/dyn4b_unconditional_zeta.json`)

Upstream: 210 checks green across 12 source ledgers.

## Checks

- [PASS] every source ledger loads and reports all_pass
- [PASS] no source ledger claims a derived zeta (boundary theorem respected program-wide)
- [PASS] BRANCH kinematic core: 44/44 first-principles + 24/24 dependency-closure checks green; DYN-9 records the core as SUSY-agnostic
- [PASS] BRANCH SUSY slice EXCLUDED: d=5 lifetime short by 7.7 orders at the unification-compatible point AND zero living points in the joint rescue scan
- [PASS] BRANCH non-SUSY ALIVE: G_LR (M_I 1e9.4, M_X 1e16.3, tau ~ 1e36) alive; PS 210-compatible marginal at 4.3e33
- [PASS] BRANCH non-SUSY: perturbativity gate makes the archival contact card slice-local (refit REQUIRED); worst PS factor 379x over 4 pi
- [PASS] BRANCH string-conditional: all three cards unpromoted with explicit negative flags
- [PASS] K1 [kinematic core]: fourth sequential chiral family kills the a-priori bound N_fam <...
- [PASS] K2 [kinematic core target]: Dirac neutrinos remove the Majorana contact target; essentiality...
- [PASS] K3 [SUSY completion (FIRED)]: unification x proton decay excludes the slice family; obstructio...
- [PASS] K4 [non-SUSY PS chain]: tau(p -> e+ pi0) = 4.3e33 vs 2.4e34: dead or inside the Hyper-K ...
- [PASS] K5 [non-SUSY branch]: e+ pi0 dominance; G_LR at 1e36 partially probeable; K+ nubar dom...
- [PASS] K6 [SUSY slice (soft)]: thermal unflavored leptogenesis fails by ~1.5 orders, wrong cent...
- [PASS] K7 [hidden-messenger extension]: odd-R spurion ceiling eps < 4.4e-4 (refreshed windows); light-se...
- [PASS] K8 [string-conditional]: instanton source requires N2, N3 ABOVE the intermediate gauge sc...
- [PASS] K9 [string-conditional]: zero-mode selectivity gap Delta S >= 6.6 (refreshed) / 3.2 (loos...
- [PASS] K10 [string-conditional]: SUSY-breaking bridge: G_LR window log10 M_SS in [10.30, 15.55]; ...
- [PASS] K11 [benchmark bookkeeping]: Z_178 diagnostic retired: misses the refreshed card by 4.1e-3 ra...
- [PASS] reconstruction paper contains the Falsifiability section with branch map and all eleven criteria K1-K11
- [PASS] reconstruction paper status ledger updated (executed audits + still-open list) and manifest extended with the collector
- [PASS] main paper status ledger updated (deferred -> executed with findings) and Z_178 retirement note added
- [PASS] self-containment discipline: the reconstruction paper uses working aliases only inside the terminology remark and literal file paths (spot check: no 'DYN-' outside texttt/remark context)
