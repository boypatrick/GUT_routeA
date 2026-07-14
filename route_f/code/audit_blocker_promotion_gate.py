#!/usr/bin/env python3
"""Fail-closed promotion gate for the recovered Route-E dynamics blockers.

This audit deliberately passes when the repository *refuses* to promote
mechanically green but physically incomplete ledgers.  It also recomputes the
three numerical corrections used by the 2026-07-14 review:

* the SU(5)-normalised one-loop ``g1`` convention in DYN-9b-2;
* the SM Davidson--Ibarra bound with one, not two, square roots;
* the electroweak-scale lower edge hidden in the RE-SC5 toy window.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "route_f" / "output"
CHECKS: list[dict[str, object]] = []


def load(rel: str) -> dict:
    return json.loads((ROOT / rel).read_text())


def sha256(rel: str) -> str:
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def check(name: str, ok: bool, detail: str = "") -> None:
    record = {"name": name, "pass": bool(ok), "detail": detail}
    CHECKS.append(record)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}"
          + (f" ({detail})" if detail else ""))


REGISTRY_PATH = "route_E/code_dyn/dyn_claim_registry.json"
PATHS = {
    "core": "route_E/output/route_e_first_principles.json",
    "re_sc3": "route_d/output/d3_instanton_majorana_pricing.json",
    "re_sc4": "route_d/output/d4_stueckelberg_protection.json",
    "re_sc5": "route_d/output/d5_susy_breaking_bridge.json",
    "dyn5v": "output/audit5/dyn5_model_validity.json",
    "dyn7f": "output/audit7/dyn7_flavor_regime_gate.json",
    "dyn9b2": "output/audit9/dyn9b2_nonsusy_flavor_refit.json",
    "dyn9b3": "output/audit9/dyn9b3_nonsusy_leptogenesis.json",
    "dyn8": "output/audit8/dyn8_falsifiability_collection.json",
}

registry = load(REGISTRY_PATH)
L = {name: load(path) for name, path in PATHS.items()}

print("== promotion gate: Cartan/Killing logical boundary ==")
ad_e = 0
B_abelian = 1
killing_abelian = ad_e * ad_e
invariance_residual = 2 * ad_e * B_abelian
check(
    "one-dimensional abelian algebra defeats original H3 but fails H3+",
    B_abelian != 0 and invariance_residual == 0 and killing_abelian == 0
    and L["core"].get("original_h3_selects_genus") is False
    and L["core"].get("nfam_upper_bound_without_h3plus") == 3
    and L["core"].get("nfam_three_requires_killing_contact") is True,
    "B(e,e)=1 is invariant/nondegenerate; Killing=0; only N_fam<=3 is unconditional",
)

print("== promotion gate: recovered-card status ==")
check(
    "RE-SC3/4/5 exist and are mechanically green",
    all(L[name].get("all_pass") is True
        for name in ("re_sc3", "re_sc4", "re_sc5")),
    ", ".join(f"{name}={L[name]['checks_passed']}/"
              f"{L[name]['checks_total']}"
              for name in ("re_sc3", "re_sc4", "re_sc5")),
)
check(
    "the recovered string cards remain unpromoted pricing cards",
    all(L[name].get("promoted_to_paper") is False
        for name in ("re_sc3", "re_sc4", "re_sc5"))
    and all(L[name].get("physics_promotion_allowed") is False
            for name in ("re_sc3", "re_sc4", "re_sc5"))
    and all(card["physics_status_after_recovery"]
            == "unpromoted_pricing_only"
            for card in registry["external_string_cards"].values()),
)
check(
    "string-card boundary flags expose their non-derived inputs",
    L["re_sc3"].get("instanton_prefactors_A_i_fixed_to_one") is True
    and L["re_sc3"].get("S_ge_1_semiclassical_control_proved") is False
    and L["re_sc3"].get("three_cycles_requires_rank1_per_cycle_assumption")
    is True
    and L["re_sc4"].get("numerical_selectivity_gap_valid") is False
    and L["re_sc5"].get("MSS_below_Mstar_is_only_a_necessary_scale_condition")
    is True,
)

print("== promotion gate: exact additive-charge no-go ==")
# Allowed NLH and XN imply qL+qH=-qN and qX=-qN.  Therefore
# q(XLH)=qX+qL+qH=-2qN=2qX=q(XX), for any additive Abelian factor.
abelian_rows = []
for q_n, q_l in product(range(-5, 6), repeat=2):
    q_h = -q_n - q_l
    q_x = -q_n
    abelian_rows.append((q_x + q_l + q_h, 2 * q_x))
check(
    "allowed NLH and XN force q(XLH)=q(XX) for every sampled additive charge",
    all(q_xlh == q_xx for q_xlh, q_xx in abelian_rows),
    f"{len(abelian_rows)} exact integer assignments",
)
check(
    "RE-SC4 numerical selectivity gaps are blocked by invalid DYN-5 input",
    L["re_sc4"]["provenance"]["dyn5_ledger"]
    == "output/audit5/dyn5_messenger_one_loop.json"
    and L["dyn5v"].get("claim_status") == "invalid_pending_rederivation"
    and L["dyn5v"].get("physics_claim_valid") is False
    and L["re_sc4"].get("numerical_selectivity_gap_valid") is False
    and L["re_sc4"].get("physics_promotion_allowed") is False
    and any("not implied by equal charge" in item
            for item in L["re_sc4"]["price_card"]["cannot_buy"]),
    "the Abelian no-go survives; Delta-S pricing does not",
)

print("== promotion gate: DYN-9b-2 hypercharge normalization ==")


def run_yt(g1_initial: float, n: int = 4000) -> float:
    """Replay the explicit Euler integrator used by DYN-9b-2."""

    a1, a2, a3, y = g1_initial, 0.648, 1.166, 0.936
    h = math.log(2e16 / 173.0) / n
    for _ in range(n):
        b1 = (41 / 10) * a1**3 / (16 * math.pi**2)
        b2 = (-19 / 6) * a2**3 / (16 * math.pi**2)
        b3 = -7 * a3**3 / (16 * math.pi**2)
        by = y * (4.5 * y*y - 8 * a3*a3 - 2.25 * a2*a2
                  - 0.85 * a1*a1) / (16 * math.pi**2)
        a1, a2, a3, y = (a1 + h*b1, a2 + h*b2, a3 + h*b3,
                         y + h*by)
    return y


yt_correct = run_yt(0.462)
yt_double = run_yt(math.sqrt(5 / 3) * 0.462)
yt_ledger = L["dyn9b2"]["derivation_log"]["S1_lock"]["yt_MX"]
check(
    "DYN-9b-2 uses already-normalised g1=0.462 with b1=41/10",
    abs(yt_ledger - yt_correct) < 1e-12
    and abs(yt_correct - 0.4411648125938088) < 1e-12,
    f"correct={yt_correct:.8f}, double-normalised={yt_double:.8f}",
)
s1_lock = L["dyn9b2"]["derivation_log"]["S1_lock"]
s3_scale = L["dyn9b2"]["derivation_log"]["S3_invariance"]
check(
    "SO(10) matrix relations do not promote the optional top-like ansatz",
    s1_lock.get("so10_relations_force_top_like_Ynu") is False
    and s1_lock["counterexample"]
    == {"h_over_f": 3.0, "Ynu_over_f": 0.0, "Yu_over_f": 4.0}
    and L["dyn9b2"].get("top_like_Ynu_is_additional_ansatz") is True,
    "h=3f gives Y_nu=0 while Y_u=4f",
)
arch_ps = s1_lock["chains"]["PS"]["archival_kernel_suppression"]
arch_lr = s1_lock["chains"]["G_LR"]["archival_kernel_suppression"]
top_ps = s1_lock["chains"]["PS"]["top_like_ansatz_tension"]
top_lr = s1_lock["chains"]["G_LR"]["top_like_ansatz_tension"]
check(
    "fixed archival-kernel suppression is separated from top-like tension",
    19.0 < arch_ps < 20.0 and 342.0 < arch_lr < 343.0
    and 9.0 < top_ps < 10.0 and 169.0 < top_lr < 170.0,
    f"archival={arch_ps:.1f}/{arch_lr:.1f}; top-like={top_ps:.1f}/{top_lr:.1f}",
)
check(
    "zeta invariance is domain-limited and complex rescaling is phase-covariant",
    s3_scale["domain_of_exact_zeta_invariance"] == "uniform positive-real y"
    and s3_scale["complex_rescaling"]["zeta_phase_invariant"] is False
    and s3_scale["complex_rescaling"]["covariance_error"] < 1e-12
    and L["dyn9b2"].get("complex_rescaling_changes_zeta_phase") is True,
    f"complex covariance error={s3_scale['complex_rescaling']['covariance_error']:.2e}",
)

print("== promotion gate: Davidson--Ibarra and flavor regime ==")
repair = L["dyn7f"]["davidson_ibarra_repair"]
m1_heavy = L["dyn7f"]["input"]["M1_GeV"]
eps_di_sm = (3 * m1_heavy * repair["delta_m_eV"] * 1e-9
             / (16 * math.pi * 174.0**2))
eps_di_9b3 = L["dyn9b3"]["derivation_log"]["S1_central"]["chain"][
    "eps_DI_max"]
check(
    "SM Davidson--Ibarra bound uses physical m3-m1 after one square root",
    abs(eps_di_sm - 2.3077948550898078e-6) < 1e-18
    and abs(eps_di_9b3 - eps_di_sm) < 1e-16,
    f"epsilon_DI_SM={eps_di_sm:.12e}",
)
check(
    "DYN-9b-3 remains blocked in the tau-resolved regime",
    L["dyn9b3"].get("claim_status")
    == "blocked_missing_branch_thermal_inputs"
    and L["dyn9b3"].get("physics_posterior_valid") is False
    and L["dyn9b3"].get("likelihood_applied") is False
    and L["dyn9b3"].get("source_selected") is False
    and L["dyn9b3"].get("success_probability_publishable") is False
    and L["dyn9b3"].get("physics_promotion_allowed") is False
    and registry["claims"]["dyn-9b-3"]["status"]
    == "blocked_missing_branch_thermal_inputs",
)

print("== promotion gate: RE-SC5 lower-edge diagnostic ==")
alive_lr = [row for row in L["re_sc5"]["scan"]["G_LR"]
            if row.get("alive")]
edge = min(alive_lr, key=lambda row: row["log10_MSS"])
edge_10tev = min((row for row in alive_lr if row["log10_MI"] >= 4.0),
                 key=lambda row: row["log10_MSS"])
check(
    "published toy-window lower edge is only about 101 GeV",
    abs(edge["log10_MI"] - 2.0056405343009227) < 1e-12
    and 100.0 < 10**edge["log10_MI"] < 102.0,
    f"log10(MSS)={edge['log10_MSS']:.2f}, "
    f"MI={10**edge['log10_MI']:.2f} GeV",
)
check(
    "illustrative MI>=10 TeV diagnostic moves the grid edge to log10(MSS)=11.85",
    abs(edge_10tev["log10_MSS"] - 11.85) < 1e-12,
    "diagnostic only; an experimental branch-specific bound is still required",
)
check(
    "RE-SC5 toy window remains non-promotable despite its mechanical scan",
    L["re_sc5"].get("physics_status") == "preliminary_toy_scan"
    and L["re_sc5"].get("physics_promotion_allowed") is False
    and L["re_sc5"].get("branch_specific_experimental_MI_floor_applied")
    is False
    and L["re_sc5"].get("MSS_below_Mstar_is_only_a_necessary_scale_condition")
    is True,
)

print("== promotion gate: collector and registry ==")
check(
    "DYN-9b-2 is preliminary rather than a claimed global flavor fit",
    L["dyn9b2"].get("full_nonsusy_yukawa_fit_performed") is False
    and registry["claims"]["dyn-9b-2"]["status"] == "preliminary",
)
check(
    "DYN-8 passes disclosure checks while forbidding physics promotion",
    L["dyn8"].get("all_pass") is True
    and L["dyn8"].get("collector_status") == "preliminary_fail_closed"
    and L["dyn8"].get("physics_promotion_allowed") is False
    and registry["claims"]["dyn-8"]["status"] == "preliminary",
)
dyn8_text = json.dumps({"branch_map": L["dyn8"].get("branch_map"),
                        "kill_bank": L["dyn8"].get("kill_bank")})
check(
    "DYN-8 contains no stale promotion or exhaustive-source wording",
    "K8 UPGRADED" not in dyn8_text
    and "only quantified scenario" not in dyn8_text
    and "EMPTY (structural)" not in dyn8_text
    and all(item.get("physics_promotable") is False
            for item in L["dyn8"]["kill_bank"]
            if item["id"] in {"K2", "K2r", "K6", "K6r", "K7", "K8",
                              "K8r", "K9", "K10", "K11r"}),
)

passed = sum(item["pass"] for item in CHECKS)
payload = {
    "audit": "Route-F blocker promotion gate",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": passed == len(CHECKS),
    "checks_passed": passed,
    "checks_total": len(CHECKS),
    "checks": CHECKS,
    "mechanical_status": "pass" if passed == len(CHECKS) else "fail",
    "physics_promotion_allowed": False,
    "resolved_artifact_blockers": ["RE-SC3", "RE-SC4", "RE-SC5"],
    "resolved_logic_blockers": [
        "Original H3 demoted: invariant/nondegenerate contact does not select genus",
        "N_fam=3 is explicitly conditional on the H3+ Killing-contact axiom",
        "Top-like Y_nu is an explicit ansatz, not a consequence of Y_nu=h-3f and Y_u=h+f",
        "Exact complex-zeta invariance is restricted to uniform positive-real rescaling",
        "Fixed archival-tower scale ordering is not a generic instanton prediction",
    ],
    "open_physics_blockers": [
        "H3+: physical motivation or UV realization of the Killing-contact axiom",
        "DYN-5: explicit interacting messenger action and operator-complete symmetry",
        "RE-SC4: recompute numerical selectivity pricing after DYN-5 rederivation",
        "RE-SC5: branch-specific experimental MI floor and threshold/proton envelopes",
        "DYN-9b-2: global branch-local non-SUSY Spin(10) flavor fit",
        "DYN-9b-3: tau-resolved kinetic or density-matrix evolution with thermal inputs",
    ],
    "numerical_repairs": {
        "yt_MX_correct_g1": yt_correct,
        "yt_MX_double_normalised_historical": yt_double,
        "epsilon_DI_SM": eps_di_sm,
        "RE_SC5_lower_edge": edge,
        "RE_SC5_MI_ge_10TeV_diagnostic_edge": edge_10tev,
        "DYN9b2_archival_kernel_suppression": {"PS": arch_ps, "G_LR": arch_lr},
        "DYN9b2_top_like_ansatz_tension": {"PS": top_ps, "G_LR": top_lr},
        "DYN9b2_complex_rescaling_covariance_error":
            s3_scale["complex_rescaling"]["covariance_error"],
    },
    "source_sha256": {
        REGISTRY_PATH: sha256(REGISTRY_PATH),
        **{path: sha256(path) for path in PATHS.values()},
    },
}

OUT.mkdir(parents=True, exist_ok=True)
json_path = OUT / "blocker_promotion_gate.json"
md_path = OUT / "blocker_promotion_gate.md"
json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

lines = [
    "# Route-F blocker promotion gate",
    "",
    f"{passed}/{len(CHECKS)} fail-closed checks pass.  This is an expected "
    "mechanical pass with `physics_promotion_allowed=false`.",
    "",
    "## Resolved artifact blockers",
    "",
    "- RE-SC3, RE-SC4, and RE-SC5 now exist and are hash-tracked.",
    "- They remain `unpromoted_pricing_only`; existence is not theory closure.",
    "",
    "## Resolved logical blocker",
    "",
    "- Original H3 proves no genus selection: `B(e,e)=1` is a "
    "nondegenerate invariant contact on the one-dimensional abelian branch.",
    "- The theorem now states `N_fam<=3` unconditionally and "
    "`N_fam=3` only under the explicit H3+ Killing-contact axiom.",
    "- `Y_nu=h-3f`, `Y_u=h+f` does not force a top-like neutrino "
    "Dirac coupling (`h=3f` is a counterexample).",
    "- Exact zeta invariance holds for uniform positive-real rescaling; "
    "complex rescaling rotates its phase by `y^2/|y|^2`.",
    "",
    "## Open physics blockers",
    "",
]
lines.extend(f"- {item}" for item in payload["open_physics_blockers"])
lines.extend(["", "## Numerical repairs", "",
              f"- Correct `y_t(M_X)` = {yt_correct:.8f}; historical "
              f"double-normalised result = {yt_double:.8f}.",
              f"- Correct SM Davidson--Ibarra bound = {eps_di_sm:.12e}.",
              f"- Fixed archival-kernel suppression = {arch_ps:.1f}x (PS) / "
              f"{arch_lr:.1f}x (G_LR); optional top-like tension = "
              f"{top_ps:.1f}x / {top_lr:.1f}x.",
              f"- RE-SC5 toy lower edge has `M_I={10**edge['log10_MI']:.2f}` "
              f"GeV; the illustrative `M_I>=10 TeV` grid edge begins at "
              f"`log10(M_SS)={edge_10tev['log10_MSS']:.2f}`.",
              "", "## Checks", ""])
lines.extend(f"- [{'PASS' if item['pass'] else 'FAIL'}] {item['name']}"
             + (f" — {item['detail']}" if item["detail"] else "")
             for item in CHECKS)
md_path.write_text("\n".join(lines) + "\n")

print(f"blocker gate: {passed}/{len(CHECKS)}; "
      f"physics_promotion_allowed=false; outputs -> {OUT.relative_to(ROOT)}")
if passed != len(CHECKS):
    raise SystemExit(1)
