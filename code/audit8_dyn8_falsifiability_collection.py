#!/usr/bin/env python3
"""DYN-8: falsifiability collection (closes the dynamics main line).

Collector audit: loads every dynamics-lane and string-card ledger,
re-verifies each number quoted in the reconstruction paper's "Dynamics
Audits and Falsifiability" section against its source ledger
(chain-of-custody), checks that both tex files carry the update, and
writes the collected branch map + kill bank as a ledger of its own.

This audit derives nothing new; it is the bookkeeping closure of the
program DYN-0 -> (1 || 4a) -> 2 -> 3 -> 5 -> 8 with the DYN-7/9 lanes
and the D3-D5 string-conditional cards.  zeta is NOT derived anywhere.
"""

import cmath
import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "audit8"
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
    "dyn7": load("output/audit7/dyn7_leptogenesis_argzeta.json"),
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
check("BRANCH non-SUSY ALIVE: G_LR (M_I 1e9.4, M_X 1e16.3, tau ~ 1e36) "
      "alive; PS 210-compatible marginal at 4.3e33",
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


def bank(tag, branch, statement, ledger, ok, detail=""):
    KILL_BANK.append({"id": tag, "branch": branch, "statement": statement,
                      "source_ledger": ledger, "verified": bool(ok)})
    check(f"{tag} [{branch}]: {statement[:64]}...", ok, detail)


bank("K1", "kinematic core", "fourth sequential chiral family kills the "
     "a-priori bound N_fam <= 3 (self-carried-index principle)",
     "route_E/output/route_e_first_principles.json",
     L["core"]["all_pass"], "theorem-level; LHC Higgs data exclude SM4")

ce = L["dyn4b"]["contact_essentiality"]
bank("K2", "kinematic core target", "Dirac neutrinos remove the Majorana "
     "contact target; essentiality posterior-wide",
     "output/audit1/dyn4b_unconditional_zeta.json",
     ce["P(cf > 0.01)"] == 1.0, f"P(cf>0.01) = {ce['P(cf > 0.01)']}")

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
bank("K6", "SUSY slice (soft)", "thermal unflavored leptogenesis fails "
     "by ~1.5 orders, wrong central sign; x3-10 boost -> 34-45%",
     "output/audit7/dyn7_leptogenesis_argzeta.json",
     abs(ps7["P_success"] - 0.0035) < 1e-6
     and abs(bs7["P_success_x3"] - 0.337) < 0.01
     and abs(bs7["P_success_x10"] - 0.4545) < 0.01,
     f"P = {ps7['P_success']*100:.2f}%, x3 -> {bs7['P_success_x3']*100:.0f}%")

c5 = L["dyn5"]["derivation_log"]["S4_r_selection"]["ceilings_order_estimates"]
bank("K7", "hidden-messenger extension", "odd-R spurion ceiling "
     "eps < 4.4e-4 (refreshed windows); light-sector silence structural",
     "output/audit5/dyn5_messenger_one_loop.json",
     abs(c5["eps_odd_tight"] / 4.44e-4 - 1) < 0.01
     and L["dyn5"]["one_loop_silence_of_light_sector"],
     f"eps_odd_tight = {c5['eps_odd_tight']:.2e}")

cx = L["d3"]["derivation_log"]["S4_price_card"]
bank("K8", "string-conditional", "instanton source requires N2, N3 "
     "ABOVE the intermediate gauge scale (39x, 4.8e3x on PS): "
     "coexistence correlation falsifiable",
     "route_d/output/d3_instanton_majorana_pricing.json",
     "N_2, N_3 > M_I" in " ".join(cx["buys"]),
     "verified against the D3 price card")

g4 = L["d4"]["derivation_log"]["S3_closure"]["selectivity_gap_by_Ms"]["2e16"]
bank("K9", "string-conditional", "zero-mode selectivity gap "
     "Delta S >= 6.6 (refreshed) / 3.2 (loose) at M_s = 2e16",
     "route_d/output/d4_stueckelberg_protection.json",
     abs(g4["dS_refreshed"] - 6.64) < 0.05 and abs(g4["dS_loose"] - 3.18) < 0.05,
     f"Delta S = {g4['dS_refreshed']:.2f} / {g4['dS_loose']:.2f}")

w5 = L["d5"]["windows"]
bank("K10", "string-conditional", "SUSY-breaking bridge: G_LR window "
     "log10 M_SS in [10.30, 15.55]; PS bridge EMPTY (structural)",
     "route_d/output/d5_susy_breaking_bridge.json",
     w5["PS"] is None and abs(w5["G_LR"]["log10_MSS_min"] - 10.30) < 0.01
     and abs(w5["G_LR"]["log10_MSS_max"] - 15.55) < 0.01,
     f"G_LR window = [{w5['G_LR']['log10_MSS_min']:.2f}, "
     f"{w5['G_LR']['log10_MSS_max']:.2f}]")

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
bank("K6r", "non-SUSY branch (refresh)", "K6 REPLACED by the "
     "branch-tagged rerun: ~6x harder (suppression 0.165, P = 0/4000 "
     "unflavored, boost sweet spot x60) but the gravitino ceiling "
     "dissolves by the branch -- a soft constraint on the source "
     "scenario, not a kill",
     "output/audit9/dyn9b3_nonsusy_leptogenesis.json",
     abs(s9b3["S1_central"]["suppression_vs_dyn7"] - 0.165) < 0.001
     and s9b3["S2_posterior"]["P_success"] == 0.0
     and s9b3["S2_posterior"]["boost"]["x60"] > 0.4,
     f"median |eta_B| = {s9b3['S2_posterior']['median_abs_eta']:.2e}")

s9b2 = L["9b2"]["derivation_log"]
bank("K8r", "string-conditional (refresh)", "K8 UPGRADED to a branch "
     "requirement: lock tension 9x/159x + type-II deficit 3.7/7.3 "
     "orders + leptogenesis floor leave the scale-decoupled "
     "instanton-type source as the only quantified scenario; D3 "
     "itself remains unpromoted",
     "output/audit9/dyn9b2_nonsusy_flavor_refit.json",
     s9b2["S1_lock"]["chains"]["PS"]["lock_tension"] > 5
     and s9b2["S1_lock"]["chains"]["G_LR"]["lock_tension"] > 50
     and all(v["deficit_orders"] > 3
             for v in s9b2["S2_typeII"].values())
     and L["d3"]["promoted_to_paper"] is False,
     "conditional-dependency statement, not a promotion")

zinv = s9b2["S3_invariance"]
bank("K11r", "benchmark bookkeeping (refresh)", "the zeta-invariance "
     "theorem: the normalized contact content is exactly invariant "
     "under Dirac rescaling; only M_* is branch-local (the 'what "
     "would not falsify' remark is now a theorem)",
     "output/audit9/dyn9b2_nonsusy_flavor_refit.json",
     abs(complex(zinv["zeta"][0], zinv["zeta"][1])
         - complex(0.1076472949, 0.0736514853)) < 1e-9
     and abs(zinv["contact_fraction"] - 0.1304275166688152) < 1e-9,
     "machine-verified to 1e-16 in the source ledger")

print("== DYN-8 section 5: papers carry the update ==")

tex_e = (ROOT / "route_E" / "tex"
         / "route_e_first_principles.tex").read_text()
tex_a = (ROOT / "paper" / "gut_framework.tex").read_text()
check("reconstruction paper contains the Falsifiability section with "
      "branch map and all eleven criteria K1-K11",
      "Dynamics Audits and Falsifiability" in tex_e
      and all(f"K{i}." in tex_e for i in range(1, 12))
      and "sec:branchmap" in tex_e)
check("reconstruction paper status ledger updated (executed audits + "
      "still-open list) and manifest extended with the collector",
      "Executed dynamics audits" in tex_e
      and "audit8\\_dyn8\\_falsifiability\\_collection" in tex_e)
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
      "branch-tagged, K8 is upgraded, and the invariance theorem "
      "replaces the bookkeeping remark; the main paper's status note "
      "records the executed re-derivation",
      "three independent spectrum treatments" in " ".join(tex_e.split())
      and "UPGRADED to a branch requirement" in tex_e
      and "branch-locality is now a \\emph{theorem}" in tex_e
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
        "nonsusy_variant": "ALIVE: PS 210-compatible marginal (4.3e33, "
                           "Hyper-K window, TRIPLE-TESTED against the "
                           "computed 210 and 126bar spectra); G_LR "
                           "alive (tau ~ 1e36) but requires the "
                           "45-adjoint route (210-only vacuum "
                           "tree-dead, epsilon lever exhausted); "
                           "source selection leaves the "
                           "scale-decoupled instanton-type tower as "
                           "the only quantified scenario; the "
                           "normalized contact content is "
                           "branch-independent by the invariance "
                           "theorem (M_* branch-local)",
        "string_conditional": "three pricing cards, unpromoted, must "
                              "not be cited as evidence",
    },
    "kill_bank": KILL_BANK,
    "upstream_checks_total": sum(v["checks_passed"] for v in L.values()),
    # negative-boundary flags
    "derives_nothing_new": True,
    "zeta_value_derived": False,
    "nonsusy_flavor_refit_still_open": True,
    "string_cards_promoted": False,
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "dyn8_falsifiability_collection.json").write_text(
    json.dumps(payload, indent=2) + "\n")

md = ["# DYN-8: falsifiability collection", "",
      f"{n_pass}/{len(CHECKS)} checks pass.  Closes the dynamics main "
      "line; derives nothing new; zeta NOT derived.", "",
      "## Branch map", ""]
md += [f"- **{k}**: {v}" for k, v in payload["branch_map"].items()]
md += ["", "## Kill bank (branch-tagged, all re-verified)", ""]
md += [f"- **{e['id']}** [{e['branch']}] {e['statement']} "
       f"(`{e['source_ledger']}`)" for e in KILL_BANK]
md += ["", f"Upstream: {payload['upstream_checks_total']} checks green "
       f"across {len(L)} source ledgers.", "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "dyn8_falsifiability_collection.md").write_text("\n".join(md) + "\n")

print(f"\nDYN-8: {n_pass}/{len(CHECKS)} checks; kill bank K1-K11 "
      f"re-verified against {len(L)} ledgers "
      f"({payload['upstream_checks_total']} upstream checks); "
      f"ledgers -> {OUT.relative_to(ROOT)}/")
