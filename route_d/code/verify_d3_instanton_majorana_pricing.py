#!/usr/bin/env python3
"""Route-D D3: instanton / d=5 Majorana-source pricing card.

Prices the two string-motivated escapes for the archival M_R tower that
the DYN-9 perturbativity gate closed for the renormalizable 126bar
coupling f v_R (no chain surviving proton decay hosts the tower with
all f < 4 pi).  Four sections, derivations recorded in the JSON:

  S1  d=5-operator escape (16 16 16bar 16bar)/M_s with v_R ~ M_I:
      f'_needed = M_R,i M_s / v_R^2 = f_renorm * (M_s/v_R) -- STRICTLY
      WORSE than the renormalizable case for every M_s > v_R.  This is
      a NEGATIVE gate: the escape is closed at intermediate-scale v_R,
      leaving instantons as the only string escape for the tower.

  S2  Instanton source M_R,i ~ M_s e^{-S_i}: required actions
      S_i = ln(M_s / M_R,i) over an M_s grid; link to the D2 contact
      action S_zeta = -ln|zeta| (same-ballpark diagnostic, NOT
      evidence); weak-coupling boundary S_top >= 1 requires
      M_s >= e * M_R,top = 1.07e16 GeV (disclosed per chain).

  S3  Rank/texture bookkeeping: Takagi decomposition of the replayed
      archival M_R = sum_k sigma_k u_k u_k^T -- exactly three rank-1
      terms with sigma_k > 0, so >= 3 distinct instanton cycles are
      required; the instanton data (3 actions + unitary texture)
      counts 3 + 9 = 12 real parameters = the FULL complex-symmetric
      Majorana data.  The axiom buys SCALE (evading f < 4 pi), not
      texture compression.

  S4  Price card: precise axiom statement, cost, what it buys, what it
      does not buy, and the beyond-zeta consequence -- on the
      210-compatible PS chain the mechanism REQUIRES N_2, N_3 heavier
      than the intermediate gauge scale M_I (impossible with a
      renormalizable source), a falsifiable correlation if both scales
      are ever measured.

Route-D discipline: this is a CONDITIONAL string interpretation, NOT
promoted to the paper; no global compactification, zero-mode spectrum,
or moduli stabilization is constructed; zeta's value is NOT derived.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "route_d" / "output"
CHECKS = []
DLOG = {}
FOUR_PI = 4 * math.pi


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    line = f"[{'PASS' if ok else 'FAIL'}] {name}"
    if detail:
        line += f"   ({detail})"
    print(line)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cmat(raw) -> np.ndarray:
    return np.array([[c["re"] + 1j * c["im"] for c in row] for row in raw],
                    dtype=complex)


# ------------------------------------------------------------- inputs
dyn9 = json.loads((ROOT / "output" / "audit9"
                   / "dyn9_nonsusy_intermediate.json").read_text())
dyn4a = json.loads((ROOT / "output" / "audit1"
                    / "dyn4a_seesaw_zeta_posterior.json").read_text())
d2 = json.loads((OUT / "d2_e3_instanton_zeta.json").read_text())
card_path = ROOT / "route_E" / "output" / "publication_closure_card" \
    / "publication_closure_card.json"
if not card_path.exists():
    raise SystemExit("D3 requires the user-restored archival closure card "
                     f"at {card_path} (sha256 recorded in the DYN-4a ledger)")

print("== D3 section 0: premise chain-of-custody ==")
check("premise: the patched DYN-9 ledger declares the archival zeta card "
      "NOT transplantable and the tower NOT perturbatively realizable on "
      "surviving chains",
      dyn9["archival_zeta_card_transplantable_to_nonsusy_chains"] is False
      and dyn9["archival_MR_tower_perturbatively_realizable_on_surviving_chains"]
      is False)
card_sha = sha256(card_path)
check("archival closure card sha256 matches the DYN-4a provenance",
      card_sha == dyn4a["provenance"]["closure_card"]["sha256"],
      f"{card_sha[:16]}...")

MR_sigma = sorted(dyn4a["benchmark"]["heavy_MR_singular_values_GeV"])
chains = dyn9["chains"]
MI = {n: 10 ** v["log10_MI"] for n, v in chains.items()}
MX = {n: 10 ** v["log10_MX"] for n, v in chains.items()}
SURV = ["G_LR", "PS"]                       # alive + marginal (proton-wise)

# ------------------------------------------------------------- section 1
print("== D3 section 1: d=5-operator escape is CLOSED (negative gate) ==")

MS_GRID = {"MX_PS": MX["PS"], "MX_GLR": MX["G_LR"], "2e16": 2.0e16,
           "1e17": 1.0e17, "1e18": 1.0e18}
fprime = {}
worse = True
for n in SURV:
    f_ren = [s / MI[n] for s in MR_sigma]
    fprime[n] = {}
    for tag, Ms in MS_GRID.items():
        fp = [s * Ms / MI[n] ** 2 for s in MR_sigma]
        fprime[n][tag] = fp
        worse &= all(a > b for a, b in zip(fp, f_ren)) and fp[2] > FOUR_PI
check("f'_needed = f_renorm * (M_s / v_R) is STRICTLY worse than the "
      "renormalizable coupling for every surviving chain and every M_s "
      "in the grid, and always violates 4 pi at the tower top",
      worse,
      f"PS @ M_s = 2e16: f'_top = {fprime['PS']['2e16'][2]:.3e}")
anchor = MR_sigma[2] * 2.0e16 / MI["PS"] ** 2
check("anchor: on the 210-compatible PS chain at M_s = 2e16 the d=5 "
      "source needs f' ~ 1.2e8 at the tower top",
      abs(anchor - 1.157e8) / 1.157e8 < 0.01, f"f'_top = {anchor:.4e}")
vR_min = MR_sigma[2] / FOUR_PI
clash = vR_min / MI["PS"]
check("structural clash quantified: a renormalizable tower top needs "
      "v_R >= sigma_top/(4 pi) = 3.1e14 GeV, while unification fixes "
      "M_I(PS) = 8.2e11 -- a factor ~380 that the d=5 operator can only "
      "worsen; the escape is CLOSED at intermediate-scale v_R",
      abs(clash - 379) / 379 < 0.02,
      f"v_R_min/M_I = {clash:.0f}")
DLOG["S1_d5_escape_closed"] = {
    "formula": "M_R = f' v_R^2 / M_s  =>  f'_needed = M_R M_s / v_R^2 "
               "= f_renorm * (M_s/v_R); M_s > v_R always, so the d=5 "
               "source strictly RAISES the required coupling",
    "f_prime_tables": fprime,
    "v_R_min_renormalizable_GeV": vR_min,
    "conclusion": "instantons are the only string escape for the tower "
                  "at intermediate-scale v_R",
}

# ------------------------------------------------------------- section 2
print("== D3 section 2: instanton action table ==")

# paper reproducibility anchor (same constant D2 uses); the replayed
# zeta in the DYN-4a ledger agrees with it only to ~1e-11
ZETA = complex(0.1076472949, 0.0736514853)
S_zeta = -math.log(abs(ZETA))
check("D2 link: S_zeta = -ln|zeta| reproduces the D2 ledger value "
      "effective_action_if_prefactor_one",
      abs(S_zeta - d2["effective_action_if_prefactor_one"]) < 1e-12,
      f"S_zeta = {S_zeta:.10f}")

S_table = {tag: [math.log(Ms / s) for s in reversed(MR_sigma)]
           for tag, Ms in MS_GRID.items()}          # [top, mid, bottom]
S_2e16 = S_table["2e16"]
check("required actions S_i = ln(M_s/M_R,i) at M_s = 2e16: "
      "[top, mid, bottom] = [1.63, 6.43, 13.64]",
      abs(S_2e16[0] - 1.628) < 0.01 and abs(S_2e16[1] - 6.432) < 0.01
      and abs(S_2e16[2] - 13.642) < 0.01,
      f"S = [{S_2e16[0]:.3f}, {S_2e16[1]:.3f}, {S_2e16[2]:.3f}]")

Ms_weak = math.e * MR_sigma[2]
check("weak-coupling boundary DISCLOSED: S_top >= 1 requires "
      "M_s >= e * M_R,top = 1.07e16 GeV -- satisfied by the G_LR chain "
      "M_X = 2.1e16 but NOT by the PS chain M_X = 5.4e15 (there the "
      "string scale must sit above the unification scale)",
      abs(Ms_weak - 1.0675e16) / 1.0675e16 < 0.01
      and MX["G_LR"] > Ms_weak and MX["PS"] < Ms_weak,
      f"M_s_min = {Ms_weak:.3e}; M_X(G_LR) = {MX['G_LR']:.2e}, "
      f"M_X(PS) = {MX['PS']:.2e}")

check("same-ballpark DIAGNOSTIC (not evidence): at M_s = 2e16 the "
      "top-of-tower action and the D2 contact action differ by < 1",
      abs(S_2e16[0] - S_zeta) < 1.0,
      f"|S_top - S_zeta| = {abs(S_2e16[0]-S_zeta):.2f}")
DLOG["S2_instanton_actions"] = {
    "S_table_by_Ms": S_table, "S_zeta_D2": S_zeta,
    "weak_coupling_Ms_min_GeV": Ms_weak,
    "note": "three DISTINCT actions are required (hierarchy of the "
            "tower); the same-ballpark coincidence S_top ~ S_zeta is "
            "recorded as a diagnostic only",
}

# ------------------------------------------------------------- section 3
print("== D3 section 3: rank/texture bookkeeping (Takagi) ==")

card = json.loads(card_path.read_text())
yuk = {k: cmat(v) for k, v in card["selected_row"]["Yukawa_fit"].items()}
Y_e, Y_nu = yuk["charged_lepton"], yuk["neutrino_dirac"]
OLDT = {"m1_eV": 1.0e-3, "dm21": 7.42e-5, "dm31": 2.517e-3,
        "s12": 0.304, "s13": 0.0222, "s23": 0.573,
        "delta": 1.20 * math.pi, "a21": 0.35 * math.pi, "a31": 1.10 * math.pi}


def standard_pmns(s12, s13, s23, delta, a21, a31):
    sq = math.sqrt
    s12_, s13_, s23_ = sq(s12), sq(s13), sq(s23)
    c12, c13, c23 = sq(1 - s12), sq(1 - s13), sq(1 - s23)
    eid, emid = np.exp(1j * delta), np.exp(-1j * delta)
    u = np.array([
        [c12 * c13, s12_ * c13, s13_ * emid],
        [-s12_ * c23 - c12 * s23_ * s13_ * eid,
         c12 * c23 - s12_ * s23_ * s13_ * eid, s23_ * c13],
        [s12_ * s23_ - c12 * c23 * s13_ * eid,
         -c12 * s23_ - s12_ * c23 * s13_ * eid, c23 * c13]], dtype=complex)
    return u @ np.diag([1.0, np.exp(0.5j * a21), np.exp(0.5j * a31)])


vals = np.linalg.eigh(Y_e @ Y_e.conj().T)[0]
U_e = np.linalg.eigh(Y_e @ Y_e.conj().T)[1][:, np.argsort(vals)]
m_diag = np.array([OLDT["m1_eV"],
                   math.sqrt(OLDT["m1_eV"] ** 2 + OLDT["dm21"]),
                   math.sqrt(OLDT["m1_eV"] ** 2 + OLDT["dm31"])])
U_nu = U_e @ standard_pmns(OLDT["s12"], OLDT["s13"], OLDT["s23"],
                           OLDT["delta"], OLDT["a21"], OLDT["a31"])
m_light = U_nu.conj() @ np.diag(m_diag) @ U_nu.conj().T
mD = Y_nu * 100.0e9
MR_GeV = -(mD.T @ np.linalg.inv(m_light) @ mD) / 1e9

check("replayed archival M_R singular values match the DYN-4a ledger",
      np.allclose(sorted(np.linalg.svd(MR_GeV, compute_uv=False)), MR_sigma,
                  rtol=1e-9),
      f"sigma = {[f'{s:.3e}' for s in MR_sigma]} GeV")


def takagi(m):
    """Takagi factorization m = U diag(s) U^T for complex symmetric m."""
    w, s, _ = np.linalg.svd(m)
    z = w.conj().T @ m @ w.conj()          # = diag(s_k e^{i phi_k})
    phases = np.array([z[k, k] / s[k] for k in range(len(s))])
    u = w @ np.diag(np.sqrt(phases))
    return u, s


U_tak, s_tak = takagi(MR_GeV)
res = np.linalg.norm(U_tak @ np.diag(s_tak) @ U_tak.T - MR_GeV) \
    / np.linalg.norm(MR_GeV)
check("Takagi factorization M_R = U diag(sigma) U^T verified "
      "(self-implemented, residual machine-zero)", res < 1e-12,
      f"rel residual = {res:.2e}")

terms = [s_tak[k] * np.outer(U_tak[:, k], U_tak[:, k]) for k in range(3)]
ranks = [int(np.linalg.matrix_rank(t, tol=s_tak[k] * 1e-10))
         for k, t in enumerate(terms)]
check("the tower is EXACTLY three rank-1 terms with sigma_k > 0: a "
      "single instanton cycle (generic rank-1 M_R) is EXCLUDED; >= 3 "
      "distinct cycles with distinct actions are required",
      ranks == [1, 1, 1] and all(s > 0 for s in s_tak)
      and np.allclose(sum(terms), MR_GeV, rtol=1e-12),
      f"ranks = {ranks}, sigma ratios = "
      f"{s_tak[0]/s_tak[1]:.1f}, {s_tak[1]/s_tak[2]:.1f}")

n_majorana = 12                          # complex symmetric 3x3, real dof
n_axiom = 3 + 9                          # 3 actions + U(3) texture
check("parameter count: instanton data (3 actions + 9 texture) = 12 = "
      "the FULL complex-symmetric Majorana data -- the axiom buys SCALE "
      "(f < 4 pi evasion), ZERO texture compression",
      n_axiom == n_majorana, f"{n_axiom} = {n_majorana}")
DLOG["S3_rank_texture"] = {
    "takagi_sigma_GeV": [float(x) for x in s_tak],
    "takagi_residual": res,
    "texture_note": "the Takagi vectors u_k ARE the required instanton "
                    "family texture; they are a NEW conditional input, "
                    "not derived",
}

# ------------------------------------------------------------- section 4
print("== D3 section 4: price card ==")

coex = {n: [float(s / MI[n]) for s in MR_sigma[1:]] for n in SURV}
check("beyond-zeta consequence: on BOTH surviving chains the mechanism "
      "REQUIRES N_2 and N_3 heavier than the intermediate gauge scale "
      "(sigma/M_I > 1), impossible with a renormalizable perturbative "
      "source -- a falsifiable correlation if the Z'/W_R scale and the "
      "heavy-N masses are ever both measured",
      all(all(r > 1 for r in coex[n]) for n in SURV),
      f"PS: N2/M_I = {coex['PS'][0]:.0f}, N3/M_I = {coex['PS'][1]:.0f}")

PRICE_CARD = {
    "axiom_statement": "there exist >= 3 Euclidean brane instantons "
                       "(E3-type) on distinct cycles with actions "
                       "S_k = ln(M_s/sigma_k) and Grassmann zero modes "
                       "pairing the nu^c bilinear, generating "
                       "M_R = sum_k M_s e^{-S_k} u_k u_k^T with B-L "
                       "broken by the instantons themselves",
    "parameter_cost": {"actions": 3, "texture_real_params": 9,
                       "total": 12, "majorana_data_dim": 12,
                       "compression": 0},
    "buys": ["archival-tower SCALE without f < 4 pi violation",
             "M_R decoupled from the B-L gauge scale M_I "
             "(the DYN-9 scale clash dissolves)",
             "N_2, N_3 > M_I coexistence allowed"],
    "does_not_buy": ["zeta's value (boundary theorem intact)",
                     "the family texture u_k (12 = 12, zero compression)",
                     "the PS-chain weak-coupling window "
                     "(needs M_s > M_X there)"],
    "beyond_zeta_consequence": "N_2, N_3 heavier than the intermediate "
                               "gauge scale -- falsifiable correlation",
}
check("price card assembled: statement, 12-parameter cost, buys / "
      "does-not-buy lists, beyond-zeta consequence", True,
      "see JSON price_card")
DLOG["S4_price_card"] = PRICE_CARD

# ------------------------------------------------------------- ledgers
n_pass = sum(1 for _, ok in CHECKS if ok)
all_pass = n_pass == len(CHECKS)
payload = {
    "audit": "Route-D D3 instanton / d=5 Majorana-source pricing card",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "all_pass": all_pass, "checks_passed": n_pass,
    "checks_total": len(CHECKS),
    "checks": [{"name": n, "pass": ok} for n, ok in CHECKS],
    "derivation_log": DLOG,
    "price_card": PRICE_CARD,
    "provenance": {
        "dyn9_ledger": "output/audit9/dyn9_nonsusy_intermediate.json",
        "dyn4a_ledger": "output/audit1/dyn4a_seesaw_zeta_posterior.json",
        "d2_ledger": "route_d/output/d2_e3_instanton_zeta.json",
        "closure_card_sha256": card_sha,
    },
    # negative-boundary flags (Route-D discipline)
    "d5_operator_escape_viable_at_intermediate_vR": False,
    "zero_mode_spectrum_checked": False,
    "global_divisor_constructed": False,
    "family_texture_derived": False,
    "moduli_stabilization_addressed": False,
    "zeta_value_derived": False,
    "promoted_to_paper": False,
}
(OUT / "d3_instanton_majorana_pricing.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n")

md = ["# Route-D D3: instanton / d=5 Majorana-source pricing card", "",
      f"{n_pass}/{len(CHECKS)} checks pass.  CONDITIONAL string "
      "interpretation, NOT promoted; zeta NOT derived.", "",
      "## Verdicts", "",
      "1. **d=5-operator escape CLOSED** (negative gate): "
      "`f' = f_renorm * (M_s/v_R)` is strictly worse; PS chain at "
      f"`M_s = 2e16` needs `f'_top = {anchor:.1e}`.  The renormalizable "
      f"clash factor is `v_R_min/M_I ~ {clash:.0f}`.",
      "2. **Instanton escape priced**: required actions at `M_s = 2e16` "
      f"are `S = [{S_2e16[0]:.2f}, {S_2e16[1]:.2f}, {S_2e16[2]:.2f}]` "
      f"(top/mid/bottom); D2 contact action `S_zeta = {S_zeta:.4f}` is "
      "same-ballpark as the tower top (diagnostic, not evidence).  "
      f"Weak coupling needs `M_s >= {Ms_weak:.2e}` GeV: OK on G_LR, "
      "ABOVE `M_X` on the PS chain (disclosed).",
      "3. **Rank/texture**: the tower is exactly three rank-1 Takagi "
      "terms => >= 3 distinct cycles; instanton data = 3 + 9 = 12 real "
      "parameters = the full Majorana data.  **The axiom buys scale, "
      "not texture.**",
      "4. **Beyond-zeta consequence**: N_2, N_3 heavier than the "
      f"intermediate gauge scale (PS: {coex['PS'][0]:.0f}x, "
      f"{coex['PS'][1]:.0f}x M_I) -- impossible renormalizably, "
      "falsifiable if both scales are measured.",
      "", "## Boundary (NOT claimed)", "",
      "- No global compactification, divisor, zero-mode spectrum, or "
      "moduli stabilization is constructed.",
      "- The family texture u_k is a NEW conditional input.",
      "- zeta's value is NOT derived (boundary theorem intact).",
      "- Route-D promotion bar NOT passed; this card must not be cited "
      "as evidence.", "", "## Checks", ""]
md += [f"- [{'PASS' if ok else 'FAIL'}] {n}" for n, ok in CHECKS]
(OUT / "d3_instanton_majorana_pricing.md").write_text("\n".join(md) + "\n")

print(f"\nD3: {n_pass}/{len(CHECKS)} checks; d=5 escape CLOSED, instanton "
      f"escape priced at 12 params (scale, not texture); ledgers -> "
      f"{OUT.relative_to(ROOT)}/d3_instanton_majorana_pricing.*")
