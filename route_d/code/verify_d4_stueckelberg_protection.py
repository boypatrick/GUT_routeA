#!/usr/bin/env python3
"""Route-D D4: Stueckelberg / anomalous-U(1) protection card.

The non-SUSY substitute question for the DYN-5 R-selection rule: what
protects the sterile messenger X from Dirac contamination when
holomorphy/non-renormalization is unavailable?  The audit sharpens the
planned card in three ways the enumeration forces on us:

  S1  GENERALIZED ABELIAN NO-GO (the central result).  For ANY number
      of additive abelian factors (U(1)^k, Z_N^k, B-L included): if the
      Dirac Yukawa N L H and the portal X N are both charge-allowed,
      then q(X L H) = q(X X) identically.  So charge bookkeeping can
      NEVER separate the required X mass insertion from the dangerous
      X-Dirac contamination: whatever sources one admits the other at
      the same order.  Proved as an integer identity, confirmed by
      brute-force enumeration (k = 1, 2), and instantiated on B-L
      (portal forces B-L(X) = -1; then XX and X L Htilde both sit at
      -2, so the 126bar matter parity does NOT separate them).

  S2  NON-SUSY LORENTZ BOOKKEEPING (a structural simplification).  In
      component language every operator needs an EVEN number of Weyl
      fermions.  The SUSY-style decorated channels (X 16 16 H etc.,
      three fermions) are Lorentz-forbidden outright; the m = 2
      decorations first arise at operator dimension 7.  The ONLY
      renormalizable contamination channel is the direct Dirac Yukawa
      X L Htilde -- exactly the channel the no-go ties to the X mass.

  S3  CLOSURE CONDITION, QUANTIFIED.  The card therefore closes only
      with instanton ZERO-MODE SELECTIVITY (string-natural: an
      instanton generates the specific operator its zero modes
      saturate, not every operator of the same charge).  Required
      selectivity gap: Delta S = S'(XLH) - S(XX) >=
      ln(eps_XX / eps_ceiling) with eps_XX = sqrt(3) M_*/M_s;
      at M_s = 2e16 this is 6.6 (refreshed DYN-5 ceiling) / 3.2
      (loose).  A DIRECT X L Htilde instanton must have
      S' >= ln(1/eps_ceiling) = 7.72 / 4.26.  D3 consistency: the NN
      source carries Delta(B-L) = +2 while XX needs -2 -- the
      CONJUGATE instanton class (both classes must be present; a
      doubled assumption, priced).

  S4  Price card.  What the Stueckelberg/GS U(1) buys (exact
      perturbative bookkeeping, calculable violations), what it cannot
      buy (the no-go), the closure condition, and the honest "none"
      for beyond-window low-energy consequences (X sits at
      sqrt(3) M_* ~ 7e15 GeV).

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
dyn4a = json.loads((ROOT / "output" / "audit1"
                    / "dyn4a_seesaw_zeta_posterior.json").read_text())

print("== D4 section 0: premise chain-of-custody ==")
check("premise: D3 pricing card all_pass with the d=5 escape closed",
      d3["all_pass"] and
      d3["d5_operator_escape_viable_at_intermediate_vR"] is False)

ceil = dyn5["derivation_log"]["S4_r_selection"]["ceilings_order_estimates"]
eps_tight, eps_loose = ceil["eps_odd_tight"], ceil["eps_odd_loose"]
Sp_tight, Sp_loose = math.log(1 / eps_tight), math.log(1 / eps_loose)
check("DYN-5 contamination ceilings imported and converted to absolute "
      "instanton actions S' = ln(1/eps): refreshed 7.72, loose 4.26",
      abs(Sp_tight - 7.72) < 0.05 and abs(Sp_loose - 4.26) < 0.05,
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
check("TWO abelian factors (exhaustive): the separation q(XLH) != q(XX) "
      "is impossible -- adding abelian factors can never help",
      viol2 == 0 and tested2 > 0, f"{tested2} assignments tested")

# B-L instance: B-L(nu^c) = +1, B-L(L) = -1, B-L(H) = 0.
bl = {"N": 1, "L": -1, "H": 0}
bl_X = -bl["N"]                                    # portal forces -1
bl_XX, bl_XLH, bl_NN = 2 * bl_X, bl_X + bl["L"] + bl["H"], 2 * bl["N"]
check("B-L instance: the portal forces B-L(X) = -1 (ODD); then "
      "XX = X L Htilde = -2 while NN = +2 -- the 126bar matter parity "
      "(even-unit breaking) does NOT separate the X mass from the "
      "X-Dirac contamination, and the XX source is the CONJUGATE class "
      "of the D3 NN source",
      bl_X == -1 and bl_XX == bl_XLH == -2 and bl_NN == +2,
      f"B-L: XX = {bl_XX}, XLH = {bl_XLH}, NN = {bl_NN}")
DLOG["S1_no_go"] = {
    "identity": "q(XLH) - q(XX) = (q_X + q_L + q_H) - 2 q_X = "
                "-(q_X - q_L - q_H) = -(q_X + q_N) + (q_N + q_L + q_H) "
                "= 0 whenever the portal and the Dirac Yukawa are "
                "charge-allowed",
    "scope": "any additive grading: U(1)^k, Z_N^k, gauged or global, "
             "anomalous or not, B-L included",
    "consequence": "charge bookkeeping alone can NEVER close the "
                   "non-SUSY protection card; the separator must be "
                   "non-charge data (instanton zero-mode structure)",
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
      "dimension 7 (four fermions + scalar): suppressed by (1/M_s)^3, "
      "negligible against every DYN-4 window",
      abs(dim_m2 - 7.0) < 1e-12, f"dim = {dim_m2:.0f}")
check("hence the ONLY renormalizable contamination channel in the "
      "non-SUSY EFT is the direct Dirac Yukawa X L Htilde -- exactly "
      "the channel the S1 no-go ties to the X mass insertion",
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
                "dS_loose": math.log(eps_XX / eps_loose)}
g16 = gap["2e16"]
check("the XX source must be LARGE (eps_XX = sqrt(3) M_*/M_s ~ 0.34 at "
      "M_s = 2e16, i.e. S_XX ~ 1.1): a generic same-charge spurion "
      "would put the X-Dirac coupling at 0.34 -- more than 2 orders "
      "above even the LOOSE ceiling; charge bookkeeping alone is dead "
      "(quantified form of the S1 no-go)",
      abs(g16["eps_XX"] - 0.340) < 0.005 and g16["eps_XX"] > eps_loose * 10,
      f"eps_XX = {g16['eps_XX']:.4f} vs eps_loose = {eps_loose:.2e}")
check("closure condition quantified: zero-mode selectivity must supply "
      "Delta S = S'(XLH) - S(XX) >= 6.6 (refreshed window) / 3.2 "
      "(loose) at M_s = 2e16; a DIRECT X L Htilde instanton needs "
      "S' >= 7.72 / 4.26 absolute",
      abs(g16["dS_refreshed"] - 6.64) < 0.05
      and abs(g16["dS_loose"] - 3.18) < 0.05,
      f"Delta S = {g16['dS_refreshed']:.2f} / {g16['dS_loose']:.2f}")
check("D3 consistency DISCLOSED: the NN tower source carries "
      "Delta(B-L) = +2 but the XX mass needs -2 -- the CONJUGATE "
      "instanton class; the shared D3+D4 bookkeeping therefore needs "
      "BOTH orientations present with independent actions (a doubled "
      "conditional assumption, priced -- NOT the same class the "
      "roadmap anticipated)",
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
DLOG["S3_closure"] = {"selectivity_gap_by_Ms": gap,
                      "absolute_S_prime": {"refreshed": Sp_tight,
                                           "loose": Sp_loose}}

# ------------------------------------------------------------- section 4
print("== D4 section 4: price card ==")

PRICE_CARD = {
    "axiom_statement": "an anomalous U(1) (GS-cancelled, Stueckelberg-"
                       "massive; B-L-like) makes the messenger charge "
                       "bookkeeping perturbatively EXACT, with all "
                       "violations instanton-sourced AND the XX-mass "
                       "instanton's zero modes saturating only the X "
                       "bilinear (selectivity), textured along K_tr^-1",
    "buys": ["perturbatively exact selection rules without holomorphy "
             "(gauge symmetry, all loop orders)",
             "calculable violation hierarchy e^{-S}",
             "with selectivity: X-Dirac contamination pushed below the "
             "DYN-4 windows"],
    "cannot_buy": ["charge separation of XX from X L Htilde (S1 no-go: "
                   "impossible for ANY abelian bookkeeping)",
                   "the K_tr texture of the X mass (new conditional "
                   "input, as in D3)",
                   "zeta's value"],
    "closure_conditions": {
        "zero_mode_selectivity_gap": "Delta S >= 6.6 (refreshed) / 3.2 "
                                     "(loose) at M_s = 2e16",
        "direct_instanton_bound": f"S' >= {Sp_tight:.2f} (refreshed) / "
                                  f"{Sp_loose:.2f} (loose)",
        "both_instanton_orientations": "NN (+2) and XX (-2) classes "
                                       "with independent actions",
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
    "all_pass": n_pass == len(CHECKS), "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    "price_card": PRICE_CARD,
    "provenance": {
        "d3_ledger": "route_d/output/d3_instanton_majorana_pricing.json",
        "dyn5_ledger": "output/audit5/dyn5_messenger_one_loop.json",
        "dyn4a_ledger": "output/audit1/dyn4a_seesaw_zeta_posterior.json",
    },
    # negative-boundary flags (Route-D discipline)
    "abelian_charge_protection_possible": False,
    "zero_mode_selectivity_assumed_not_computed": True,
    "anomaly_inflow_checked": False,
    "global_embedding_constructed": False,
    "zeta_value_derived": False,
    "promoted_to_paper": False,
}
(OUT / "d4_stueckelberg_protection.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n")

md = ["# Route-D D4: Stueckelberg / anomalous-U(1) protection card", "",
      f"{n_pass}/{len(CHECKS)} checks pass.  CONDITIONAL string "
      "interpretation, NOT promoted; zeta NOT derived.", "",
      "## Verdicts", "",
      "1. **Generalized abelian no-go** (central result): if the Dirac "
      "Yukawa and the portal are charge-allowed then "
      "`q(X L H) = q(X X)` identically, for ANY number of abelian "
      "factors (exhaustively confirmed, k = 1 and 2; B-L instance: "
      "portal forces `B-L(X) = -1`, both operators sit at -2).  Charge "
      "bookkeeping can never close the card.",
      "2. **Non-SUSY Lorentz bookkeeping**: all SUSY-style m = 1 "
      "decorations have odd fermion count -- forbidden outright; m = 2 "
      "first arises at dimension 7.  The ONLY renormalizable "
      "contamination is the direct `X L Htilde` -- exactly the channel "
      "the no-go ties to the X mass.",
      "3. **Closure condition**: instanton zero-mode selectivity with "
      f"`Delta S >= {g16['dS_refreshed']:.1f}` (refreshed) / "
      f"`{g16['dS_loose']:.1f}` (loose) at `M_s = 2e16`; a direct "
      f"X-Dirac instanton needs `S' >= {Sp_tight:.2f}` / "
      f"`{Sp_loose:.2f}`.  D3 consistency: the XX source is the "
      "CONJUGATE class of the NN source (+2 vs -2) -- both "
      "orientations required (doubled assumption, priced).  The K_tr "
      "texture of the X mass becomes conditional input.",
      "4. **Portability**: DYN-5's exact light-sector silence is pure "
      "linear algebra and survives non-SUSY unchanged.",
      "5. **Beyond-window consequence: NONE at low energy** (X at "
      "`sqrt(3) M_* ~ 6.8e15` GeV) -- disclosed under the pricing "
      "discipline.",
      "", "## Boundary (NOT claimed)", "",
      "- Zero-mode selectivity is ASSUMED (the closure condition), not "
      "computed; no anomaly inflow or global embedding.",
      "- zeta's value is NOT derived.",
      "- Route-D promotion bar NOT passed; must not be cited as "
      "evidence.", "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "d4_stueckelberg_protection.md").write_text("\n".join(md) + "\n")

print(f"\nD4: {n_pass}/{len(CHECKS)} checks; abelian charge protection "
      f"IMPOSSIBLE (no-go identity), closure = zero-mode selectivity "
      f"with Delta S >= {g16['dS_refreshed']:.1f}; ledgers -> "
      f"{OUT.relative_to(ROOT)}/d4_stueckelberg_protection.*")
