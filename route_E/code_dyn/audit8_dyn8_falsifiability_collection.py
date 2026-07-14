#!/usr/bin/env python3
"""DYN-8: fail-closed falsifiability/status collection.

Collector audit: loads every dynamics-lane and string-card ledger,
re-verifies each number quoted in the reconstruction paper's "Dynamics
Audits and Falsifiability" section against its source ledger
(chain-of-custody), checks that both tex files carry the update, and
writes the collected branch map + kill bank as a ledger of its own.

This audit derives nothing new and does not close the physics program.  It
collects the bookkeeping chain DYN-0 -> (1 || 4a) -> 2 -> 3 -> 5 -> 8 with
the DYN-7/9 lanes and the D3-D5 string-conditional cards, while consuming the
DYN-5V and DYN-7F validity guards.  ``all_pass`` below means that the
chain-of-custody and disclosure checks pass; it is not a promotion verdict.
zeta is NOT derived anywhere.
"""

import cmath
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from route_e_paths import AUDIT_OUTPUT, REPO_ROOT, ROUTE_E_ROOT

ROOT = REPO_ROOT
OUT = AUDIT_OUTPUT / "audit8"
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


def load(rel):
    return json.loads((ROOT / rel).read_text())


L = {
    "core": load("route_E/output/route_e_first_principles.json"),
    "closure": load("route_E/output/route_e_dependency_closure.json"),
    "dyn2": load("output/audit3/dyn2_thresholds_unification.json"),
    "dyn2b": load("output/audit3/dyn2b_rescue_scan.json"),
    "dyn3": load("output/audit2/dyn3_proton_d5_kill_criterion.json"),
    "dyn4b": load("output/audit1/dyn4b_unconditional_zeta.json"),
    "dyn5": load("output/audit5/dyn5_messenger_one_loop.json"),
    "dyn5v": load("output/audit5/dyn5_model_validity.json"),
    "dyn7": load("output/audit7/dyn7_leptogenesis_argzeta.json"),
    "dyn7f": load("output/audit7/dyn7_flavor_regime_gate.json"),
    "dyn9": load("output/audit9/dyn9_nonsusy_intermediate.json"),
    "d3": load("route_d/output/d3_instanton_majorana_pricing.json"),
    "d4": load("route_d/output/d4_stueckelberg_protection.json"),
    "d5": load("route_d/output/d5_susy_breaking_bridge.json"),
    # DYN-9b series (the 9b-4 refresh)
    "9b1": load("output/audit9/dyn9b1_nonsusy_vacuum_thresholds.json"),
    "9b1b": load("output/audit9/dyn9b1b_210_quartic_descent.json"),
    "9b1c": load("output/audit9/dyn9b1c_eps_and_126_quartics.json"),
    "9b2": load("output/audit9/dyn9b2_nonsusy_flavor_refit.json"),
    "9b3": load("output/audit9/dyn9b3_nonsusy_leptogenesis.json"),
    "9b1d": load("output/audit9/dyn9b1d_lr_ratio_scan.json"),
    "4c": load("output/audit1/dyn4c_kernel_dirac_refit.json"),
}

print("== DYN-8 section 1: chain of custody ==")
check("every source ledger loads and reports all_pass",
      all(v.get("all_pass") for v in L.values()),
      f"{len(L)} ledgers, "
      f"{sum(v['checks_passed'] for v in L.values())} upstream checks")
check("no source ledger claims a derived zeta (boundary theorem "
      "respected program-wide)",
      all(v.get("zeta_value_derived") is False for v in L.values()
          if "zeta_value_derived" in v))
check("validity guards fail closed: DYN-5 is invalid pending an interacting "
      "messenger action and DYN-7/9b-3 require flavored thermal inputs",
      L["dyn5v"].get("claim_status") == "invalid_pending_rederivation"
      and L["dyn5v"].get("physics_claim_valid") is False
      and L["dyn7f"].get("physics_claims_closed") is False
      and L["dyn7f"]["statuses"]["DYN-7"]["success_probability_publishable"]
      is False
      and L["dyn7f"]["statuses"]["DYN-9b-3"]
      ["success_probability_publishable"] is False)

print("== DYN-8 section 2: branch map re-verified ==")

check("BRANCH kinematic core: 44/44 first-principles + 24/24 "
      "dependency-closure checks green; DYN-9 records the core as "
      "SUSY-agnostic",
      L["core"]["checks_passed"] == 44
      and L["closure"]["checks_passed"] == 24)

k3a = L["dyn3"]["kill_criterion"]
check("BRANCH SUSY slice EXCLUDED: d=5 lifetime short by 7.7 orders at "
      "the unification-compatible point AND zero living points in the "
      "joint rescue scan",
      abs(k3a["orders_short"] - 7.677891420646367) < 1e-9
      and L["dyn2b"]["living_points"] == [],
      f"orders_short = {k3a['orders_short']:.2f}, living points = 0")

g = L["dyn9"]["chains"]["G_LR"]
p = L["dyn9"]["chains"]["PS"]
check("BRANCH non-SUSY PRELIMINARY LEDGER: G_LR has M_I 1e9.4, "
      "M_X 1e16.3 and the central d=6 estimate near 1e36; PS is a "
      "210-compatible marginal central estimate at 4.3e33",
      abs(g["log10_MI"] - 9.43) < 0.01 and abs(g["log10_MX"] - 16.32) < 0.01
      and g["tau_d6_years"] > 2.4e34
      and abs(p["tau_d6_years"] / 4.3e33 - 1) < 0.05,
      f"tau: G_LR {g['tau_d6_years']:.1e}, PS {p['tau_d6_years']:.1e}")
check("BRANCH non-SUSY: perturbativity gate makes the archival contact "
      "card slice-local (refit REQUIRED); worst PS factor 379x over "
      "4 pi",
      L["dyn9"]["archival_zeta_card_transplantable_to_nonsusy_chains"]
      is False
      and abs(max(p["seesaw_f_needed"]) / (4 * math.pi) / 379 - 1) < 0.02)

check("BRANCH string-conditional: all three cards unpromoted with "
      "explicit negative flags",
      all(L[c]["promoted_to_paper"] is False for c in ("d3", "d4", "d5")))

print("== DYN-8 section 3: kill bank re-verified ==")

KILL_BANK = []


def bank(tag, branch, statement, ledger, ok, detail="", promotable=True):
    KILL_BANK.append({"id": tag, "branch": branch, "statement": statement,
                      "source_ledger": ledger, "verified": bool(ok),
                      "physics_promotable": bool(promotable)})
    check(f"{tag} [{branch}]: {statement[:64]}...", ok, detail)


bank("K1", "kinematic core", "fourth sequential chiral family kills the "
     "a-priori bound N_fam <= 3 (self-carried-index principle)",
     "route_E/output/route_e_first_principles.json",
     L["core"]["all_pass"], "theorem-level; LHC Higgs data exclude SM4")

ce = L["dyn4b"]["contact_essentiality"]
bank("K2", "conditional Majorana-card diagnostic", "Dirac neutrinos remove "
     "this Majorana contact target; within the stated nuisance-prior card, "
     "all sampled contact fractions exceed 0.01 (not a data posterior)",
     "output/audit1/dyn4b_unconditional_zeta.json",
     ce["P(cf > 0.01)"] == 1.0,
     f"conditional sample fraction = {ce['P(cf > 0.01)']}",
     promotable=False)

bank("K3", "SUSY completion (FIRED)", "unification x proton decay "
     "excludes the slice family; obstruction map recorded",
     "output/audit2 + output/audit3 ledgers",
     k3a["orders_short"] > 7 and L["dyn2b"]["living_points"] == [],
     "criterion has fired; branch dead")

bank("K4", "non-SUSY PS chain", "tau(p -> e+ pi0) = 4.3e33 vs 2.4e34: "
     "dead or inside the Hyper-K discovery window",
     "output/audit9/dyn9_nonsusy_intermediate.json",
     p["marginal_within_10x"] and not p["alive_now"],
     f"tau = {p['tau_d6_years']:.2e}")

bank("K5", "non-SUSY branch", "e+ pi0 dominance; G_LR at 1e36 partially "
     "probeable; K+ nubar dominance would disfavor the branch",
     "output/audit9/dyn9_nonsusy_intermediate.json",
     g["tau_d6_years"] > 1e35, f"tau(G_LR) = {g['tau_d6_years']:.1e}")

ps7 = L["dyn7"]["posterior_statistics"]
bs7 = L["dyn7"]["boost_sensitivity"]
bank("K6", "SUSY slice (historical diagnostic)", "unflavored regression "
     "gives a small success fraction, but DYN-7F places M1 in the "
     "tau-resolved regime; probabilities and arg-zeta windows are not "
     "publishable",
     "output/audit7/dyn7_leptogenesis_argzeta.json",
     abs(ps7["P_success"] - 0.0035) < 1e-6
     and abs(bs7["P_success_x3"] - 0.337) < 0.01
     and abs(bs7["P_success_x10"] - 0.4545) < 0.01,
     f"P = {ps7['P_success']*100:.2f}%, x3 -> {bs7['P_success_x3']*100:.0f}%",
     promotable=False)

c5 = L["dyn5"]["derivation_log"]["S4_r_selection"]["ceilings_order_estimates"]
bank("K7", "hidden-messenger extension (invalid model diagnostic)",
     "historical odd-R ceiling replays arithmetically, but DYN-5V shows "
     "that the displayed bilinear has no asserted one-loop anomalous "
     "dimension and allows X L H_u",
     "output/audit5/dyn5_messenger_one_loop.json",
     c5["eps_odd_tight"] > 0
     and L["dyn5"]["one_loop_silence_of_light_sector"]
     and L["dyn5v"].get("claim_status") == "invalid_pending_rederivation"
     and L["dyn5v"].get("physics_claim_valid") is False,
     f"eps_odd_tight = {c5['eps_odd_tight']:.2e}", promotable=False)

cx = L["d3"]["derivation_log"]["S4_price_card"]
bank("K8", "string-conditional", "for the fixed archival tower, N2 and "
     "N3 lie ABOVE the intermediate gauge scale (39x, 4.8e3x on PS); "
     "the instanton benchmark permits this coexistence but does not "
     "require it for a generic refitted tower",
     "route_d/output/d3_instanton_majorana_pricing.json",
     "N_2, N_3 > M_I" in " ".join(cx["buys"]),
     "verified against the unpromoted D3 price card", promotable=False)

g4 = L["d4"]["derivation_log"]["S3_closure"]["selectivity_gap_by_Ms"]["2e16"]
bank("K9", "string-conditional (invalid numerical diagnostic)",
     "the historical DYN-5 ceilings replay to a finite zero-mode "
     "selectivity gap, but DYN-5V invalidates those ceilings, so no "
     "numerical Delta S bound is physics-promotable",
     "route_d/output/d4_stueckelberg_protection.json",
     math.isfinite(g4["dS_refreshed"])
     and math.isfinite(g4["dS_loose"])
     and L["d4"].get("numerical_selectivity_gap_valid") is False
     and L["d4"].get("physics_promotion_allowed") is False,
     f"Delta S = {g4['dS_refreshed']:.2f} / {g4['dS_loose']:.2f}; "
     "blocked because the ceiling is inherited from invalid DYN-5",
     promotable=False)

w5 = L["d5"]["windows"]
bank("K10", "string-conditional", "preliminary SUSY-breaking toy scan: "
     "G_LR window log10 M_SS in [10.30, 15.55]; PS window empty under "
     "the stated one-loop/ESH/degenerate-spectrum assumptions",
     "route_d/output/d5_susy_breaking_bridge.json",
     w5["PS"] is None and abs(w5["G_LR"]["log10_MSS_min"] - 10.30) < 0.01
     and abs(w5["G_LR"]["log10_MSS_max"] - 15.55) < 0.01,
     f"G_LR window = [{w5['G_LR']['log10_MSS_min']:.2f}, "
     f"{w5['G_LR']['log10_MSS_max']:.2f}] (unpromoted toy scan)",
     promotable=False)

zc = L["dyn4b"]["refreshed_central_card"]["zeta"]
miss = abs(2 * math.pi * 17 / 178 - cmath.phase(complex(zc[0], zc[1])))
bank("K11", "benchmark bookkeeping", "Z_178 diagnostic retired: misses "
     "the refreshed card by 4.1e-3 rad vs the 5.4e-5 window (76x)",
     "output/audit1/dyn4b_unconditional_zeta.json",
     abs(miss - 4.098e-3) < 1e-5 and miss / 5.424119e-5 > 70,
     f"miss = {miss:.3e} rad = {miss/5.424119e-5:.0f}x window")

print("== DYN-8 section 4 (9b-4 REFRESH): the non-SUSY series ==")

f1b = L["9b1b"]["derivation_log"]["S4_verdict"]["fractions"]
k4c = L["9b1c"]["derivation_log"]["S4_K4_retest"]
bank("K4r", "non-SUSY PS chain (refresh)", "K4 TRIPLE-TESTED: the "
     "Hyper-K window survives the computed 210 spectrum (tree "
     "positivity 39/250, tau = 4.25e33 all percentiles) AND the "
     "computed 126bar spectrum (45/150 viable, same tau)",
     "output/audit9/dyn9b1b + dyn9b1c ledgers",
     f1b["PS"]["tree_positive"] == 39
     and f1b["PS"]["tau_16_50_84"] == ["4.25e+33"] * 3
     and k4c["viable"] == 45
     and all(abs(t / 4.25e33 - 1) < 0.01 for t in k4c["tau_percentiles"]),
     "three independent spectrum treatments agree")

s1d = L["9b1d"]["derivation_log"]
bank("K5r", "non-SUSY LR chain (refresh, CORRECTED by the ratio scan)",
     "the 210-only LR vacuum IS a tree-level local minimum in rare "
     "coupling regions at vev ratios near 0.6-0.7 (8/8800 sampled; "
     "He-Meljanac 1986 reproduced at claim level); the earlier "
     "fixed-ratio negative (0/249 at a/p = 0.8) was an artifact of "
     "the fixed ratio; the epsilon invariant neither rescues nor "
     "kills any sampled point; the adjoint route is an ALTERNATIVE, "
     "not a necessity",
     "output/audit9/dyn9b1d_lr_ratio_scan.json (corrects dyn9b1b/1c)",
     s1d["S4_conclusion"]["q1_positive"] is True
     and s1d["S3_scan"]["stats"]["pos_eps0"] == 8
     and s1d["S3_scan"]["stats"]["rescued_by_eps"] == 0
     and s1d["S3_scan"]["stats"]["killed_by_eps"] == 0
     and set(k for k, v in s1d["S3_scan"]["by_ratio"].items() if v)
     <= {"0.6", "0.65", "0.7"},
     "prior-art conflict RESOLVED as a fixed-ratio artifact")

s9b3 = L["9b3"]["derivation_log"]
r9b3 = s9b3["S2_prior_draw_regression"]
bank("K6r", "non-SUSY branch (refresh)", "K6 is replaced by a "
     "branch-tagged unflavored regression: ~6x harder (suppression 0.165, "
     "0/4000 target-band hits without boost, largest scanned hit fraction "
     "near x60); no likelihood or thermal kinetics make these viability "
     "probabilities, and absence of a gravitino ceiling is not a reheating "
     "prediction",
     "output/audit9/dyn9b3_nonsusy_leptogenesis.json",
     abs(s9b3["S1_central"]["suppression_vs_dyn7"] - 0.165) < 0.001
     and r9b3["target_band_hit_fraction"] == 0.0
     and r9b3["boost_target_band_hit_fractions"]["x60"] > 0.4
     and L["9b3"]["likelihood_applied"] is False
     and L["9b3"]["success_probability_publishable"] is False,
     f"median |eta_B| = {r9b3['median_abs_eta']:.2e}; regression only",
     promotable=False)

s9b2 = L["9b2"]["derivation_log"]
bank("K8r", "string-conditional (refresh)", "K8 remains an unpromoted "
     "conditional diagnostic: the fixed archival-kernel ceiling and "
     "optional top-like ansatz tension, together with the type-II estimate, "
     "disfavor the fixed-kernel renormalizable benchmark, while the "
     "leptogenesis comparison lacks a global flavor fit and flavored "
     "thermal kinetics; D3 itself remains unpromoted",
     "output/audit9/dyn9b2_nonsusy_flavor_refit.json",
     s9b2["S1_lock"]["so10_relations_force_top_like_Ynu"] is False
     and s9b2["S1_lock"]["chains"]["PS"]
     ["archival_kernel_suppression"] > 10
     and s9b2["S1_lock"]["chains"]["G_LR"]
     ["archival_kernel_suppression"] > 300
     and all(v["deficit_orders"] > 3
             for v in s9b2["S2_typeII"].values())
     and L["d3"]["promoted_to_paper"] is False,
     "conditional-dependency statement, not a promotion", promotable=False)

s4c = L["4c"]["derivation_log"]
bank("K2r", "conditional kernel-refit diagnostic", "K2 is refined by the "
     "kernel refit: the theta_23 tension is a frozen-anchor property "
     "(exact Majorana absorption at unchanged kernels; a 5% Y_e "
     "perturbation absorbs it alone) and contact essentiality is "
     "lifted to the perturbation level (cf never below 0.01 over "
     "1000 refits with eps up to 0.3)",
     "output/audit1/dyn4c_kernel_dirac_refit.json",
     L["4c"]["all_pass"]
     and min(v["cf_min"] for v in s4c["S2"].values()) > 0.01
     and s4c["S3"]["eps_absorb"] is not None
     and s4c["S3"]["eps_absorb"] <= 0.1,
     f"eps_e(theta23) = {s4c['S3']['eps_absorb']}", promotable=False)

zinv = s9b2["S3_invariance"]
bank("K11r", "benchmark bookkeeping (refresh)", "positive-real rescaling "
     "theorem: within the fixed-kernel y>0 family, normalized contact "
     "content is exactly invariant and only M_* moves; for complex y, "
     "zeta is phase-covariant rather than invariant",
     "output/audit9/dyn9b2_nonsusy_flavor_refit.json",
     abs(complex(zinv["zeta"][0], zinv["zeta"][1])
         - complex(0.1076472949, 0.0736514853)) < 1e-9
     and abs(zinv["contact_fraction"] - 0.1304275166688152) < 1e-9
     and zinv["domain_of_exact_zeta_invariance"]
     == "uniform positive-real y"
     and zinv["complex_rescaling"]["zeta_phase_invariant"] is False,
     "positive-real invariance and complex phase covariance both verified",
     promotable=False)

print("== DYN-8 section 5: papers carry the update ==")

tex_e = (ROUTE_E_ROOT / "tex"
         / "route_e_first_principles.tex").read_text()
tex_a = (ROOT / "paper" / "gut_framework.tex").read_text()
check("reconstruction paper contains the Falsifiability section with "
      "branch map and all eleven criteria K1-K11",
      "Dynamics Audits and Falsifiability" in tex_e
      and all(f"K{i}." in tex_e for i in range(1, 12))
      and "sec:branchmap" in tex_e)
check("reconstruction paper status ledger updated (mechanical evidence + "
      "still-open list) and manifest extended with the collector",
      "Mechanically executed dynamics evidence" in tex_e
      and "audit8_dyn8_falsifiability_collection.py" in tex_e)
check("main paper status ledger updated (deferred -> executed with "
      "findings) and Z_178 retirement note added",
      "executed July 2026" in tex_a and "Retirement note, July 2026"
      in tex_a)
check("self-containment discipline: the reconstruction paper uses "
      "working aliases only inside the terminology remark and literal "
      "file paths (spot check: no 'DYN-' outside texttt/remark context)",
      tex_e.count("DYN-0''--``DYN-9") == 1
      and "route A" not in tex_e.replace("``route~A,''", ""))
check("9b-4 REFRESH carried by both papers: the reconstruction paper's "
      "branch map records the executed re-derivation (vacuum / source / "
      "contact stages), K4 cites three spectrum treatments, K6 is "
      "branch-tagged, K8 is explicitly unpromoted, and the invariance theorem "
      "replaces the bookkeeping remark; the main paper's status note "
      "records the executed re-derivation",
      "three independent spectrum treatments" in " ".join(tex_e.split())
      and "not a promoted branch requirement" in " ".join(tex_e.split())
      and "scaling covariance, not an" in tex_e
      and "lock theorem" in tex_e
      and "\\zeta'=e^{2i\\arg y}\\zeta" in tex_e
      and "dyn9b*.json" in tex_a.replace("\\", ""))

# ---------------------------------------------------------------- ledger
n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "DYN-8 falsifiability collection",
    "dyn_item": "DYN-8",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "branch_map": {
        "kinematic_core": "SUSY-agnostic, machine-verified (44+24), "
                          "untouched by all dynamics results",
        "susy_minimal_slice": "EXCLUDED (slice family, not model class): "
                              "unification x proton decay + obstruction "
                              "map",
        "nonsusy_variant": "PRELIMINARY: PS 210-compatible marginal (4.3e33, "
                           "Hyper-K window, TRIPLE-TESTED against the "
                           "computed 210 and 126bar spectra); G_LR "
                           "alive (tau ~ 1e36); the 210-only LR vacuum "
                           "has rare sampled local minima near vev ratios "
                           "0.6--0.7, so the 45-adjoint route is an "
                           "alternative rather than a necessity; "
                           "among the two unpromoted D3 mechanisms, only "
                           "the scale-decoupled instanton benchmark survives "
                           "the price card; this is not an exhaustive source "
                           "classification; the "
                           "normalized contact content is "
                           "branch-independent by the invariance "
                           "theorem (M_* branch-local)",
        "string_conditional": "three pricing cards, unpromoted, must "
                              "not be cited as evidence",
    },
    "kill_bank": KILL_BANK,
    "upstream_checks_total": sum(v["checks_passed"] for v in L.values()),
    "collector_status": "preliminary_fail_closed",
    "physics_promotion_allowed": False,
    "blockers": [
        "DYN-5 has no valid interacting messenger action and permits X L H_u",
        "DYN-7 and DYN-9b-3 require branch-local flavored kinetics",
        "DYN-9b-2 has not performed a global non-SUSY Spin(10) flavor fit",
        "RE-SC3/4/5 remain unpromoted pricing cards; RE-SC4 inherits an "
        "invalid DYN-5 numerical ceiling",
    ],
    # negative-boundary flags
    "derives_nothing_new": True,
    "zeta_value_derived": False,
    "nonsusy_flavor_refit_still_open": True,
    "string_cards_promoted": False,
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "dyn8_falsifiability_collection.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-8: fail-closed falsifiability collection", "",
      f"{n_pass}/{len(CHECKS)} mechanical/disclosure checks pass.  The "
      "physics main line remains open; this collector derives nothing new; "
      "zeta is NOT derived.", "",
      "## Branch map", ""]
md += [f"- **{k}**: {v}" for k, v in payload["branch_map"].items()]
md += ["", "## Kill bank (branch-tagged, all re-verified)", ""]
md += [f"- **{e['id']}** [{e['branch']}] {e['statement']} "
       f"(`{e['source_ledger']}`; physics-promotable="
       f"`{str(e['physics_promotable']).lower()}`)" for e in KILL_BANK]
md += ["", "## Blocking physics work", ""]
md += [f"- {item}" for item in payload["blockers"]]
md += ["", f"Upstream: {payload['upstream_checks_total']} checks green "
       f"across {len(L)} source ledgers.", "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn8_falsifiability_collection.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-8: {n_pass}/{len(CHECKS)} mechanical/disclosure checks; "
      f"physics promotion BLOCKED; kill bank K1-K11 re-verified against "
      f"{len(L)} ledgers ({payload['upstream_checks_total']} upstream "
      f"checks); ledgers -> {OUT.relative_to(ROOT)}/")
if n_pass != len(CHECKS):
    raise SystemExit(1)
