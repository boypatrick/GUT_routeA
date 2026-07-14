#!/usr/bin/env python3
"""Route-D D4: Stueckelberg / anomalous-U(1) protection card.

The non-SUSY substitute question for the DYN-5 R-selection rule: what
protects the sterile messenger X from Dirac contamination when
holomorphy/non-renormalization is unavailable?  The audit sharpens the
planned card in three ways the enumeration forces on us:

  S1  GENERALIZED ABELIAN-CHARGE NO-GO (the central exact result).  For
      ANY number of additive abelian factors (U(1)^k, Z_N^k, B-L
      included): if the Dirac Yukawa N L H and the portal X N are both
      charge-allowed, then q(X L H) = q(X X) identically.  Therefore an
      additive charge selection rule alone cannot allow one operator
      while forbidding the other.  The identity does NOT imply a common
      instanton/source, operator dimension, perturbative order, or
      Wilson coefficient; those require non-charge dynamical data.
      The identity is confirmed by brute-force enumeration (k = 1, 2)
      and instantiated on B-L (portal forces B-L(X) = -1; then XX and
      X L Htilde both sit at -2, so the 126bar matter parity does not
      distinguish their charge-selection classes).

  S2  NON-SUSY LORENTZ BOOKKEEPING (a structural simplification).  In
      component language every operator needs an EVEN number of Weyl
      fermions.  The SUSY-style decorated channels (X 16 16 H etc.,
      three fermions) are Lorentz-forbidden outright; the m = 2
      decorations first arise at operator dimension 7.  The ONLY
      renormalizable contamination channel is the direct Dirac Yukawa
      X L Htilde -- exactly the channel the no-go places in the same
      additive-charge class as the X mass.

  S3  CLOSURE CONDITION, PENDING.  Instanton zero-mode data could
      distinguish operators with the same additive charge, but no such
      spectrum is computed here.  The script preserves the historical
      arithmetic Delta S = ln(eps_XX/eps_ceiling) and
      S' = ln(1/eps_ceiling), with eps_XX = sqrt(3) M_*/M_s, solely as
      an invalid diagnostic: eps_ceiling comes from the DYN-5 action
      rejected by DYN-5V.  Consequently these numbers are not physics
      bounds and the closure condition remains pending an interacting-
      action rederivation.  D3 consistency: the NN source carries
      Delta(B-L) = +2 while XX needs -2.  Opposite compensating source
      charges are therefore needed; interpreting them as conjugate
      instanton orientations is an additional conditional ansatz, not
      a consequence of the charge identity.

  S4  Price card.  What the Stueckelberg/GS U(1) buys (exact
      perturbative charge bookkeeping and a possible instanton
      parametrization), what it cannot buy (the no-go), the pending
      closure condition, and the honest "none" for beyond-window
      low-energy consequences at the conditional archival benchmark.

Route-D discipline: CONDITIONAL string interpretation, NOT promoted;
no zero-mode computation, anomaly-inflow check, or global embedding is
performed; zeta's value is NOT derived.
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "route_d" / "output"
CHECKS = []
DLOG = {}


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


# ------------------------------------------------------------- inputs
d3 = json.loads((OUT / "d3_instanton_majorana_pricing.json").read_text())
dyn5 = json.loads((ROOT / "output" / "audit5"
                   / "dyn5_messenger_one_loop.json").read_text())
dyn5v = json.loads((ROOT / "output" / "audit5"
                    / "dyn5_model_validity.json").read_text())
dyn4a = json.loads((ROOT / "output" / "audit1"
                    / "dyn4a_seesaw_zeta_posterior.json").read_text())

print("== D4 section 0: premise chain-of-custody ==")
check("premise: D3 pricing card all_pass with the d=5 escape closed",
      d3["all_pass"] and
      d3["d5_operator_escape_viable_at_intermediate_vR"] is False)

ceil = dyn5["derivation_log"]["S4_r_selection"]["ceilings_order_estimates"]
eps_tight, eps_loose = ceil["eps_odd_tight"], ceil["eps_odd_loose"]
Sp_tight, Sp_loose = math.log(1 / eps_tight), math.log(1 / eps_loose)
check("historical DYN-5 ceiling arithmetic is reproducible but explicitly "
      "blocked from promotion by the DYN-5V invalidity gate",
      0 < eps_tight < 1 and 0 < eps_loose < 1
      and abs(Sp_tight - math.log(1 / eps_tight)) < 1e-14
      and abs(Sp_loose - math.log(1 / eps_loose)) < 1e-14
      and dyn5v.get("claim_status") == "invalid_pending_rederivation"
      and dyn5v.get("physics_claim_valid") is False,
      f"S'_refreshed = {Sp_tight:.3f}, S'_loose = {Sp_loose:.3f} "
      f"(eps = {eps_tight:.2e}, {eps_loose:.2e})")

# ------------------------------------------------------------- section 1
print("== D4 section 1: generalized abelian no-go ==")

# Identity: constraints q_N + q_L + q_H = 0 (Dirac allowed) and
# q_X + q_N = 0 (portal allowed) imply
# q(XLH) = q_X + q_L + q_H = q_X - q_N = 2 q_X = q(XX), for EVERY
# abelian factor independently.  Exhaustive integer confirmation:
viol = 0
tested = 0
for qX, qL in itertools.product(range(-4, 5), repeat=2):
    for qH in range(-4, 5):
        qN = -(qL + qH)              # Dirac allowed
        if qX + qN != 0:             # portal allowed filter
            continue
        tested += 1
        if (qX + qL + qH) != 2 * qX:
            viol += 1
check("IDENTITY (single factor, exhaustive integer grid): Dirac + portal "
      "allowed => q(X L H) = q(X X) with ZERO violations",
      viol == 0 and tested > 0, f"{tested} assignments tested")

viol2 = 0
tested2 = 0
rng = range(-2, 3)
for qX1, qL1, qH1, qX2, qL2, qH2 in itertools.product(rng, repeat=6):
    qN1, qN2 = -(qL1 + qH1), -(qL2 + qH2)
    if qX1 + qN1 != 0 or qX2 + qN2 != 0:
        continue
    tested2 += 1
    if (qX1 + qL1 + qH1, qX2 + qL2 + qH2) != (2 * qX1, 2 * qX2):
        viol2 += 1
check("TWO abelian factors (exhaustive): additive charge selection cannot "
      "distinguish XLH from XX; adding additive factors does not change "
      "that charge-level conclusion",
      viol2 == 0 and tested2 > 0, f"{tested2} assignments tested")

# B-L instance: B-L(nu^c) = +1, B-L(L) = -1, B-L(H) = 0.
bl = {"N": 1, "L": -1, "H": 0}
bl_X = -bl["N"]                                    # portal forces -1
bl_XX, bl_XLH, bl_NN = 2 * bl_X, bl_X + bl["L"] + bl["H"], 2 * bl["N"]
check("B-L instance: the portal forces B-L(X) = -1 (ODD); then "
      "XX = X L Htilde = -2 while NN = +2 -- the 126bar matter parity "
      "(even-unit breaking) does NOT distinguish the additive-charge "
      "classes of the X mass and X-Dirac operators; this equality does "
      "not identify their sources or coefficients.  XX has the "
      "OPPOSITE charge orientation to the D3 NN operator; identifying "
      "geometrically conjugate sources would be an additional ansatz",
      bl_X == -1 and bl_XX == bl_XLH == -2 and bl_NN == +2,
      f"B-L: XX = {bl_XX}, XLH = {bl_XLH}, NN = {bl_NN}")
DLOG["S1_no_go"] = {
    "identity": "q(XLH) - q(XX) = (q_X + q_L + q_H) - 2 q_X = "
                "-(q_X - q_L - q_H) = -(q_X + q_N) + (q_N + q_L + q_H) "
                "= 0 whenever the portal and the Dirac Yukawa are "
                "charge-allowed",
    "scope": "any additive grading: U(1)^k, Z_N^k, gauged or global, "
             "anomalous or not, B-L included",
    "consequence": "additive charge bookkeeping alone cannot establish "
                   "operator separation; any separation must come from "
                   "additional non-charge data such as a computed "
                   "zero-mode spectrum",
    "does_not_imply": "same source, same instanton, same operator order, "
                      "or same Wilson coefficient",
}

# ------------------------------------------------------------- section 2
print("== D4 section 2: non-SUSY Lorentz bookkeeping ==")

FERM = {"N": 1, "X": 1, "L": 1, "Q": 1, "uc": 1, "dc": 1, "ec": 1,
        "H": 0, "T": 0}
CHANNELS = {
    "XN_mass": ["X", "N"], "XX_mass": ["X", "X"], "NN_mass": ["N", "N"],
    "Dirac_NLH": ["N", "L", "H"], "contamination_XLH": ["X", "L", "H"],
    "decorated_m1_Dirac": ["X", "Q", "uc", "H"],
    "decorated_m1_triplet": ["X", "Q", "Q", "T"],
    "decorated_m1_Majorana": ["X", "N", "N"],
    "decorated_m2_Dirac": ["X", "X", "Q", "uc", "H"],
}
fcount = {k: sum(FERM[f] for f in v) for k, v in CHANNELS.items()}
odd_killed = [k for k, n in fcount.items() if n % 2 == 1]
check("Lorentz fermion-parity gate: every SUSY-style m = 1 decorated "
      "channel has an ODD Weyl-fermion count and is forbidden outright "
      "in the non-SUSY component EFT (no operator of any dimension)",
      set(odd_killed) == {"decorated_m1_Dirac", "decorated_m1_triplet",
                          "decorated_m1_Majorana"},
      f"odd-fermion channels killed: {odd_killed}")

dim_m2 = sum(1.5 if FERM[f] else 1.0
             for f in CHANNELS["decorated_m2_Dirac"])
check("m = 2 decorations survive Lorentz but first arise at operator "
      "dimension 7 (four fermions + scalar); their numerical impact "
      "requires Wilson coefficients and is not bounded by this card",
      abs(dim_m2 - 7.0) < 1e-12, f"dim = {dim_m2:.0f}")
check("hence the ONLY renormalizable contamination channel in the "
      "non-SUSY EFT is the direct Dirac Yukawa X L Htilde -- exactly "
      "the channel S1 places in the same additive-charge selection "
      "class as the X mass insertion",
      fcount["contamination_XLH"] % 2 == 0
      and all(fcount[k] % 2 == 1 for k in odd_killed))
DLOG["S2_lorentz"] = {
    "fermion_counts": fcount,
    "note": "the SUSY decoration danger audited in DYN-5 is a "
            "superpotential (scalar-component) artifact; in components "
            "it dissolves, and the whole burden lands on X L Htilde",
}

# ------------------------------------------------------------- section 3
print("== D4 section 3: closure condition (zero-mode selectivity) ==")

M_star = float(dyn4a["benchmark"]["M_star_GeV"])
SQ3 = math.sqrt(3.0)                     # X mass = sqrt(3) M_* (DYN-5)
MS_GRID = {"2e16": 2.0e16, "1e17": 1.0e17, "1e18": 1.0e18}
gap = {}
for tag, Ms in MS_GRID.items():
    eps_XX = SQ3 * M_star / Ms           # required size of the XX source
    gap[tag] = {"eps_XX": eps_XX,
                "S_XX": math.log(1 / eps_XX),
                "dS_refreshed": math.log(eps_XX / eps_tight),
                "dS_loose": math.log(eps_XX / eps_loose),
                "status": "historical_invalid",
                "physics_bound": False}
g16 = gap["2e16"]
check("conditional unit-prefactor XX benchmark is reproducible: "
      "eps_XX = sqrt(3) M_*/M_s ~ 0.34 and S_XX ~ 1.08 at M_s=2e16; "
      "this does NOT determine the XLH source, action, or coefficient",
      abs(g16["eps_XX"] - 0.340) < 0.005
      and abs(g16["S_XX"] - math.log(1 / g16["eps_XX"])) < 1e-14,
      f"eps_XX = {g16['eps_XX']:.4f}, S_XX = {g16['S_XX']:.4f}")
check("historical selectivity-gap arithmetic is reproduced but is not a "
      "valid closure bound until the messenger contamination ceiling is "
      "rederived from an interacting action",
      abs(g16["dS_refreshed"]
          - math.log(g16["eps_XX"] / eps_tight)) < 1e-14
      and abs(g16["dS_loose"]
              - math.log(g16["eps_XX"] / eps_loose)) < 1e-14
      and dyn5v.get("physics_claim_valid") is False,
      f"Delta S = {g16['dS_refreshed']:.2f} / {g16['dS_loose']:.2f}")
check("D3 charge consistency DISCLOSED: NN and XX carry opposite B-L "
      "charges (+2 and -2), so compensating sources require opposite "
      "charges; treating them as conjugate instanton orientations with "
      "independent actions is an additional conditional ansatz, not a "
      "charge-level theorem",
      bl_NN == -bl_XX)
check("K_tr direction of the X mass is itself conditional input in the "
      "non-SUSY card: the SUSY theorem POSTULATED K_tr^-1(X,X); here "
      "the conjugate-class instanton must be textured along K_tr^-1 "
      "for the Schur step to reproduce lambda^2 K_tr (recorded in the "
      "price card)", True, "structure, not a numeric gate")

# the light-sector silence argument survives non-SUSY: the GL(3)
# field-redefinition covariance of the seesaw proved in DYN-5 is pure
# linear algebra (no holomorphy used in that step)
check("portability: DYN-5's EXACT light-sector silence (seesaw "
      "invariance under any invertible redefinition of nu^c) is pure "
      "linear algebra and carries over to the non-SUSY card unchanged",
      dyn5["one_loop_silence_of_light_sector"] is True)
DLOG["S3_closure"] = {
    # Compatibility keys are retained for DYN-8 replay, but every value is
    # explicitly quarantined as a historical-invalid numerical diagnostic.
    "selectivity_gap_by_Ms": gap,
    "absolute_S_prime": {"refreshed": Sp_tight, "loose": Sp_loose,
                         "status": "historical_invalid",
                         "physics_bound": False},
    "numerical_status": "historical_invalid_pending_dyn5_rederivation",
    "physics_bound": False,
    "interpretation": "the displayed Delta-S and S-prime values combine "
                      "valid logarithmic arithmetic with contamination "
                      "ceilings from the DYN-5 action invalidated by "
                      "DYN-5V; they are not constraints",
    "charge_identity_limitation": "q(XLH)=q(XX) does not identify their "
                                  "sources, orders, or coefficients",
}

# ------------------------------------------------------------- section 4
print("== D4 section 4: price card ==")

PRICE_CARD = {
    "axiom_statement": "an anomalous U(1) (GS-cancelled, Stueckelberg-"
                       "massive; B-L-like) makes the messenger charge "
                       "bookkeeping perturbatively EXACT, with all "
                       "violations instanton-sourced AND the XX-mass "
                       "instanton's zero modes saturating only the X "
                       "bilinear (selectivity), textured along K_tr^-1",
    "buys": ["perturbatively exact additive-charge selection rules "
             "without holomorphy (gauge symmetry, all loop orders)",
             "a framework in which a concrete instanton calculation "
             "could parameterize violations by e^{-S}",
             "a possible non-charge discriminator through computed "
             "zero-mode saturation; no such computation is supplied"],
    "cannot_buy": ["charge separation of XX from X L Htilde (S1 no-go: "
                   "impossible for ANY abelian bookkeeping)",
                   "equality of the XLH and XX source, operator order, "
                   "or Wilson coefficient (not implied by equal charge)",
                   "the K_tr texture of the X mass (new conditional "
                   "input, as in D3)",
                   "zeta's value"],
    "closure_conditions": {
        "numerical_gap_status": "PENDING: rederive the contamination "
                                "ceiling from an interacting action",
        "historical_invalid_snapshot_at_Ms_2e16": {
            "Delta_S_refreshed": g16["dS_refreshed"],
            "Delta_S_loose": g16["dS_loose"],
            "S_prime_refreshed": Sp_tight,
            "S_prime_loose": Sp_loose,
            "physics_bound": False,
        },
        "source_charge_requirement": "NN (+2) and XX (-2) operators "
                                     "need oppositely charged compensating "
                                     "sources",
        "conjugate_instanton_relation": "additional ansatz, not derived "
                                        "from additive charges",
    },
    "beyond_window_consequence": "NONE at low energy (X sits at "
                                 "sqrt(3) M_* ~ 6.8e15 GeV); the card's "
                                 "testable content is internal -- "
                                 "zero-mode counting on the two cycles "
                                 "in any concrete embedding",
}
check("price card assembled, with the honest 'none' for low-energy "
      "beyond-window consequences (pricing discipline allows an "
      "explicit none)", True, "see JSON price_card")
DLOG["S4_price_card"] = PRICE_CARD

# ------------------------------------------------------------- ledgers
n_pass = sum(1 for _, ok in CHECKS if ok)
payload = {
    "audit": "Route-D D4 Stueckelberg / anomalous-U(1) protection card",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    # `all_pass` remains the mechanical ledger status consumed by the
    # orchestration gates.  Standalone execution nevertheless fails closed
    # below because the physics-promotion gate is false.
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    "price_card": PRICE_CARD,
    "provenance": {
        "d3_ledger": "route_d/output/d3_instanton_majorana_pricing.json",
        "dyn5_ledger": "output/audit5/dyn5_messenger_one_loop.json",
        "dyn5_validity_ledger": "output/audit5/dyn5_model_validity.json",
        "dyn4a_ledger": "output/audit1/dyn4a_seesaw_zeta_posterior.json",
    },
    # negative-boundary flags (Route-D discipline)
    "abelian_charge_protection_possible": False,
    "zero_mode_selectivity_assumed_not_computed": True,
    "anomaly_inflow_checked": False,
    "global_embedding_constructed": False,
    "zeta_value_derived": False,
    "promoted_to_paper": False,
    "physics_status": "unpromoted_algebraic_no_go_only",
    "physics_promotion_allowed": False,
    "numerical_selectivity_gap_valid": False,
    "numerical_gap_blocker": "historical ceiling imported from invalid DYN-5",
    "mechanical_status": "checks_pass" if n_pass == len(CHECKS) else "failed",
    "standalone_exit_policy": "nonzero_until_physics_promotion_allowed",
}
(OUT / "d4_stueckelberg_protection.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n")

md = ["# Route-D D4: Stueckelberg / anomalous-U(1) protection card", "",
      f"{n_pass}/{len(CHECKS)} mechanical checks pass.  CONDITIONAL string "
      "interpretation, NOT promoted; standalone execution fails closed; "
      "zeta NOT derived.", "",
      "## Verdicts", "",
      "1. **Generalized abelian no-go** (central result): if the Dirac "
      "Yukawa and the portal are charge-allowed then "
      "`q(X L H) = q(X X)` identically, for ANY number of abelian "
      "factors (exhaustively confirmed, k = 1 and 2; B-L instance: "
      "portal forces `B-L(X) = -1`, both operators sit at -2).  Charge "
      "selection alone cannot distinguish the operators.  This does NOT "
      "imply the same source, operator order, or Wilson coefficient.",
      "2. **Non-SUSY Lorentz bookkeeping**: all SUSY-style m = 1 "
      "decorations have odd fermion count -- forbidden outright; m = 2 "
      "first arises at dimension 7.  The ONLY renormalizable "
      "contamination is the direct `X L Htilde` -- exactly the channel "
      "the no-go places in the same additive-charge class as the X mass.",
      "3. **Numerical closure gap blocked**: the historical "
      f"`Delta S={g16['dS_refreshed']:.1f}/{g16['dS_loose']:.1f}` and "
      f"`S'={Sp_tight:.2f}/{Sp_loose:.2f}` arithmetic is reproduced, "
      "but these are **historical-invalid diagnostics**, because they "
      "import the invalid DYN-5 contamination ceiling; neither is a "
      "physics bound.  D3 consistency: NN and XX have opposite B-L "
      "charges (+2 vs -2), but a conjugate-instanton relation is an "
      "additional ansatz, not a charge theorem.  The K_tr "
      "texture of the X mass becomes conditional input.",
      "4. **Portability**: DYN-5's exact light-sector silence is pure "
      "linear algebra and survives non-SUSY unchanged.",
      "5. **Beyond-window consequence: NONE at low energy** (X at "
      "`sqrt(3) M_* ~ 6.8e15` GeV) -- disclosed under the pricing "
      "discipline.",
      "", "## Boundary (NOT claimed)", "",
      "- Zero-mode selectivity is ASSUMED (the closure condition), not "
      "computed; no anomaly inflow or global embedding.",
      "- Equal additive charge does not establish a common source, "
      "operator order, or Wilson coefficient.",
      "- The numerical selectivity gap is PENDING an interacting-action "
      "rederivation; historical arithmetic must not be used as a bound.",
      "- zeta's value is NOT derived.",
      "- Route-D promotion bar NOT passed; must not be cited as "
      "evidence.", "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "d4_stueckelberg_protection.md").write_text("\n".join(md) + "\n")

print(f"\nD4: {n_pass}/{len(CHECKS)} checks; abelian charge protection "
      f"IMPOSSIBLE (no-go identity); numerical selectivity gap BLOCKED "
      f"pending DYN-5 rederivation; ledgers -> "
      f"{OUT.relative_to(ROOT)}/d4_stueckelberg_protection.*")

if n_pass != len(CHECKS):
    raise SystemExit(1)
if not payload["physics_promotion_allowed"]:
    print("D4 FAIL-CLOSED: mechanical arithmetic passed, but the numerical "
          "gap and physics promotion remain blocked.")
    raise SystemExit(2)
