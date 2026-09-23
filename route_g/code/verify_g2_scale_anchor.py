#!/usr/bin/env python3
"""Synthetic scale-identifiability audit; no observed masses or physical unit choice.

All dimensionful detector and beam settings co-scale in the unit-rescaling
tests. A real externally calibrated energy reference breaks this degeneracy.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

import verify_g2_local_conversion as g2
import verify_g2_driven as driven
import verify_g2_resolution as resolution

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_scale_anchor"
SYNTHETIC_SCALES = [.1, 1.0, 10.0]
SIGMA_I, SIGMA_D = .1, .05


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scale_card(scale):
    return dict(g2.CARD, R=g2.CARD["R"]/scale,
                M5=g2.CARD["M5"]*scale, MD=g2.CARD["MD"]*scale,
                incoming_com_momentum=g2.CARD["incoming_com_momentum"]*scale)


def packet_at_scale(scale):
    card = scale_card(scale)
    labels, _, weights = driven.finite_packet()
    rows = []
    for l, weight in zip(labels, weights):
        enumeration = driven.enumerate_driven(card["incoming_phi_n"], int(l),
                                               card["incoming_com_momentum"], card=card,
                                               epsilon=driven.DRIVE["epsilon"],
                                               omega=driven.DRIVE["omega"]*scale)
        rows.extend(dict(row, weighted_rate=float(weight)*row["rate_coefficient"])
                    for row in enumeration["open_channels"])
    rate = sum(row["weighted_rate"] for row in rows)
    work_coefficient = sum(row["weighted_rate"]*row["pump_work"] for row in rows)
    for row in rows:
        row["conditional_share"] = row["weighted_rate"]/rate
    return dict(card=card, rows=rows, total_rate=rate, work_coefficient=work_coefficient,
                work_per_event=work_coefficient/rate)


def beam_components_at_scale(scale):
    card = scale_card(scale)
    p, iw, _ = resolution.incoming_quadrature(SIGMA_I)
    labels, _, weights = driven.finite_packet()
    out = []
    for l, sign in [(0, 1), (-2, -1)]:
        wl = weights[labels.tolist().index(l)]
        for q in [-1, 0, 1]:
            channels = [driven.driven_channel(1, l, 0, q, float(momentum)*scale,
                                              card=card, epsilon=driven.DRIVE["epsilon"],
                                              omega=driven.DRIVE["omega"]*scale) for momentum in p]
            out.append(dict(l=l, q=q, sign=sign,
                            means=np.array([row["outgoing_momentum"] for row in channels]),
                            rates=wl*iw*np.array([row["rate_coefficient"] for row in channels])))
    return out


def reconstruct_center_zero(y_minus, y_zero, y_plus):
    """Three externally anchored signed labels -1,0,+1, not sorted masses."""
    A = (y_plus-2*y_zero+y_minus)/2
    B = (y_plus-y_minus)/2
    C = y_zero
    if A <= 0:
        return dict(A=float(A), B=float(B), C=float(C), admissible=False,
                    reason="A must be positive for a real finite radius")
    alpha = B/(2*A)
    bulk_mass_squared = C-B*B/(4*A)
    return dict(A=float(A), B=float(B), C=float(C), R=float(1/np.sqrt(A)),
                Lambda=float(np.sqrt(A)), alpha=float(alpha),
                bulk_mass_squared=float(bulk_mass_squared),
                admissible=bool(bulk_mass_squared >= 0),
                positivity_scope="M5^2>=0 is the declared positive bulk-mass card; it is checked, not inferred from a three-point fit.")


def run():
    checks = []

    def check(name, passed, value=None, tolerance=None):
        row = dict(name=name, passed=bool(passed))
        if value is not None: row["measured"] = float(value)
        if tolerance is not None: row["tolerance"] = float(tolerance)
        checks.append(row)

    dependencies = dict(g2_source=Path(g2.__file__), driven_source=Path(driven.__file__),
                        resolution_source=Path(resolution.__file__),
                        g2_artifact=ROOT/"output"/"g2_local_conversion.json",
                        driven_artifact=ROOT/"output"/"g2_driven.json",
                        resolution_artifact=ROOT/"output"/"g2_resolution.json")
    hashes = {key: sha(path) for key, path in dependencies.items()}
    for stem in ["g2", "driven", "resolution"]:
        artifact = json.loads(dependencies[stem+"_artifact"].read_text())
        if stem == "g2":
            recorded = artifact["source_sha256"]
        elif stem == "driven":
            recorded = artifact["source_sha256"]["verify_g2_driven.py"]
        else:
            recorded = artifact["verification"]["source_sha256"]
        check("source_hash:"+stem, recorded == hashes[stem+"_source"])
    base_packet = packet_at_scale(1)
    base_keyed = {(r["l"], r["m"], r["j"], r["q"]): r for r in base_packet["rows"]}
    base_beam = beam_components_at_scale(1)
    base_classification = resolution.classify(base_beam, SIGMA_D)
    parent = json.loads(dependencies["resolution_artifact"].read_text())
    old = next(row for row in parent["scan"] if row["driven"] and row["sigma_i"] == SIGMA_I and row["sigma_d"] == SIGMA_D)
    check("baseline:resolution_artifact_recovery", abs(base_classification["bayes_error"]-old["bayes_error"]) < 1e-13)
    scaled = []
    spectral = []
    for scale in SYNTHETIC_SCALES:
        packet = packet_at_scale(scale)
        card = packet["card"]
        keyed = {(r["l"], r["m"], r["j"], r["q"]): r for r in packet["rows"]}
        check(f"scale={scale}:open_channel_support", set(keyed) == set(base_keyed))
        rate_errors, p_errors, share_errors = [], [], []
        for key, row in keyed.items():
            before = base_keyed[key]
            rate_error = abs(row["rate_coefficient"]*scale*scale/before["rate_coefficient"]-1)
            p_error = abs(row["outgoing_momentum"]/scale-before["outgoing_momentum"])
            share_error = abs(row["conditional_share"]-before["conditional_share"])
            rate_errors.append(rate_error); p_errors.append(p_error); share_errors.append(share_error)
            check(f"scale={scale}:channel={key}:K_inverse_square", rate_error < 1e-11, rate_error, 1e-11)
            check(f"scale={scale}:channel={key}:momentum_linear", p_error < 1e-11, p_error, 1e-11)
            check(f"scale={scale}:channel={key}:share_invariant", share_error < 1e-13)
            check(f"scale={scale}:channel={key}:work_linear", abs(row["pump_work"]/scale-before["pump_work"]) < 1e-14)
        check(f"scale={scale}:total_rate", abs(packet["total_rate"]*scale*scale/base_packet["total_rate"]-1) < 1e-13)
        check(f"scale={scale}:work_coefficient", abs(packet["work_coefficient"]*scale/base_packet["work_coefficient"]-1) < 1e-12)
        check(f"scale={scale}:mean_work", abs(packet["work_per_event"]/scale/base_packet["work_per_event"]-1) < 1e-12)
        kappa5 = 2*np.pi*card["R"]*card["g"]
        check(f"scale={scale}:kappa5", abs(kappa5*scale-2*np.pi*g2.CARD["R"]*g2.CARD["g"]) < 1e-14)
        beam = beam_components_at_scale(scale)
        classification = resolution.classify(beam, SIGMA_D*scale)
        check(f"scale={scale}:beam_error_invariant", abs(classification["bayes_error"]-base_classification["bayes_error"]) < 1e-12)
        check(f"scale={scale}:beam_priors_invariant", np.max(abs(np.array(classification["priors"])-np.array(base_classification["priors"]))) < 1e-13)
        check(f"scale={scale}:threshold_linear", len(classification["bayes_boundaries"]) == len(base_classification["bayes_boundaries"]) and np.max(abs(np.array(classification["bayes_boundaries"])/scale-np.array(base_classification["bayes_boundaries"]))) < 1e-11)
        check(f"scale={scale}:beam_rate_inverse_square", abs(classification["total_pair_rate"]*scale*scale/base_classification["total_pair_rate"]-1) < 1e-13)
        # Check an undriven explicit-card channel, independently of Floquet plumbing.
        undriven = g2.channel(1, 0, 0, .2*scale, card=card)
        old_undriven = g2.channel(1, 0, 0, .2)
        check(f"scale={scale}:stationary_K", abs(undriven["rate_coefficient"]*scale*scale/old_undriven["rate_coefficient"]-1) < 1e-12)
        scaled.append(dict(synthetic_scale=scale, scale_is_physical_unit_assignment=False,
                           card=card, omega=driven.DRIVE["omega"]*scale,
                           sigma_i=SIGMA_I*scale, sigma_d=SIGMA_D*scale, kappa5=kappa5,
                           open_channel_count=len(keyed), total_rate_coefficient=packet["total_rate"],
                           pump_work_rate_coefficient=packet["work_coefficient"],
                           mean_pump_work_per_scattering=packet["work_per_event"],
                           classification=classification,
                           max_channel_rate_scaling_error=max(rate_errors),
                           max_rescaled_momentum_error=max(p_errors),
                           max_conditional_share_error=max(share_errors)))
        # Synthetic spectrum, with genuinely signed labels (mass order differs).
        labels = [-1, 0, 1, 2]
        y = [g2.phi_mass(n, card)**2 for n in labels]
        fit = reconstruct_center_zero(*y[:3])
        third_difference = y[3]-3*y[2]+3*y[1]-y[0]
        predicted_fourth = 4*fit["A"]+2*fit["B"]+fit["C"]
        check(f"scale={scale}:three_anchor_admissible", fit["admissible"])
        for name, actual, expected in [("A", fit["A"], scale**2), ("alpha", fit["alpha"], .25),
                                       ("M5sq", fit["bulk_mass_squared"], .25*scale**2),
                                       ("Lambda", fit["Lambda"], scale)]:
            check(f"scale={scale}:three_anchor:{name}", abs(actual-expected) < 1e-12*max(1, abs(expected)))
        check(f"scale={scale}:fourth_anchor_closure", abs(third_difference)/scale**2 < 1e-12)
        distortion = .07*scale**2
        distorted_fourth = y[3]+distortion
        rejected_residual = distorted_fourth-predicted_fourth
        check(f"scale={scale}:distorted_fourth_reject", abs(rejected_residual)/scale**2 > .069)
        shift = .3*scale**2
        shifted = [value+shift for value in y]
        shifted_fit = reconstruct_center_zero(*shifted[:3])
        check(f"scale={scale}:common_shift_A", abs(shifted_fit["A"]-fit["A"])/scale**2 < 1e-12)
        check(f"scale={scale}:common_shift_B", abs(shifted_fit["B"]-fit["B"])/scale**2 < 1e-12)
        check(f"scale={scale}:common_shift_M5sq", abs(shifted_fit["bulk_mass_squared"]-fit["bulk_mass_squared"]-shift)/scale**2 < 1e-12)
        check(f"scale={scale}:common_shift_fourth_closure", abs(shifted[3]-3*shifted[2]+3*shifted[1]-shifted[0])/scale**2 < 1e-12)
        # One-mass anchoring is algebraically possible only after fixing a ratio.
        synthetic_anchor_mass = np.sqrt(y[1])
        inferred_scale = synthetic_anchor_mass/g2.phi_mass(0)
        check(f"scale={scale}:single_anchor_conditional_inference", abs(inferred_scale-scale) < 1e-12)
        for k in [-2, 1, 3]:
            shifted_card = dict(card, alpha=card["alpha"]+k)
            check(f"scale={scale}:holonomy_relabel:k={k}", max(abs(g2.phi_mass(n-k, shifted_card)-g2.phi_mass(n, card)) for n in labels)/scale < 1e-12)
        check(f"scale={scale}:neutral_detector_degeneracy", g2.detector_mass(1, card) == g2.detector_mass(-1, card))
        spectral.append(dict(synthetic_scale=scale, signed_labels=labels, squared_masses=y,
                             reconstruction=fit, predicted_fourth_squared_mass=predicted_fourth,
                             third_difference=third_difference,
                             deliberately_distorted_fourth_squared_mass=distorted_fourth,
                             distorted_fourth_residual=rejected_residual,
                             common_positive_squared_mass_shift=shift,
                             shifted_squared_masses=shifted, shifted_reconstruction=shifted_fit,
                             synthetic_single_mass_anchor=synthetic_anchor_mass,
                             inferred_scale_if_baseline_ratios_frozen=inferred_scale))
    check("admissibility:reject_nonpositive_A", not reconstruct_center_zero(1, 2, 1)["admissible"])
    # Positive observed y values need not imply the stipulated M5^2 >= 0 card.
    check("admissibility:reject_negative_bulk_M5sq", not reconstruct_center_zero(2.1, .1, .1)["admissible"])
    check("dependencies:unchanged", all(sha(path) == hashes[key] for key, path in dependencies.items()))
    failed = [row["name"] for row in checks if not row["passed"]]
    result = dict(status="Synthetic dimensional-consistency and identifiability audit, not a physical energy-scale determination.",
                  baseline_dimensionless_card=g2.CARD, baseline_drive=driven.DRIVE,
                  synthetic_scales=SYNTHETIC_SCALES,
                  scale_convention="Lambda=1/R_phys; the baseline R_hat=1. Natural units. 0.1,1,10 are synthetic unit-rescalings, not GeV values or recommended physical scales.",
                  equations=dict(spectrum="m_n^2=Lambda^2[0.5^2+(n+0.25)^2]; mu_j^2=Lambda^2[5^2+j^2] for the frozen baseline ratios.",
                                 scaling="R=1/Lambda; (M5,MD,p,Omega,sigma_i,sigma_d,t_classifier) scale as Lambda; g,alpha,epsilon and normalized compact weights remain fixed.",
                                 coupling="kappa5=2 pi R g scales as Lambda^(-1).",
                                 rate="Kq=|gq|^2 pf/(16 pi En El Eout) scales as Lambda^(-2).",
                                 work="q Omega scales as Lambda; sum w Kq q Omega scales as Lambda^(-1); its ratio to sum w Kq scales as Lambda.",
                                 probability="P_i=W_i/sum W is invariant. Under x=Lambda xhat, rho_Lambda(x)=rho_hat(x/Lambda)/Lambda; Bayes error is invariant when resolution and threshold co-scale.",
                                 nonidentifiability="Without a dimensional anchor, d P_i/d log Lambda=0 and dimensionless multinomial Fisher information for log Lambda is zero.",
                                 single_anchor="Lambda=m_n,calibrated/sqrt(M5_hat^2+(n+alpha)^2), only if dimensionless ratios and the mode identification are fixed externally.",
                                 three_anchor="y_n=A n^2+B n+C; for signed n=-1,0,+1: A=(y_1-2y_0+y_-1)/2, B=(y_1-y_-1)/2, C=y_0; R=A^(-1/2), alpha=B/(2A), M5^2=C-B^2/(4A).",
                                 fourth_anchor="y_2-3 y_1+3 y_0-y_-1=0, equivalently y_2=4A+2B+C.",
                                 universal_shift="y_n -> y_n+delta M^2 implies A,B unchanged and C,M5^2 shifted by delta M^2; second and third finite differences are unchanged."),
                  qualifications=["Dimensionful measured momenta, a calibrated drive frequency, a known mass, or an absolute rate with controlled couplings/luminosity can break scale degeneracy. None has been supplied here.",
                                  "Adopting a calibrated drive frequency together with the frozen Omega*R=.2 card defines a candidate R; the arbitrary drive-frequency card is not itself a derived radius or a measured mode spacing.",
                                  "Fixing a mass by hand calibrates a free scale; it is not a prediction of that mass.",
                                  "Three masses require assigned consecutive signed KK labels. Sorting observed masses does not provide those labels.",
                                  "Integer holonomy relabeling alpha->alpha+k,n->n-k leaves the spectrum invariant. Alpha is a class modulo integers, and orientation reversal can add ambiguity without a signed-mode reference.",
                                  "The assumed positive-bulk-mass card requires A>0 and M5^2>=0; three arbitrary observations need not satisfy it.",
                                  "A fourth signed-mode anchor supplies a real closure test. The displayed distorted example is synthetic, not an excluded experimental observation.",
                                  "The common positive mass-squared shift is a spectral algebra test only. No Higgs portal, detector interaction or threshold dynamics is added by this script.",
                                  "Neutral X masses still obey mu_j=mu_-j. A scale anchor does not solve the recoil-sign readout problem.",
                                  "Scale invariance of Bayes error presumes ALL dimensional beam, drive and detector settings co-scale. Holding a real detector resolution fixed does not preserve it.",
                                  "The inherited Gaussian beam, tree EFT, perfect mode labels and finite-time/absolute-rate limitations remain. Synthetic rescaling cannot establish EFT validity or calibrate hardware."],
                  scaled_predictions=scaled, synthetic_spectral_anchors=spectral,
                  checks=checks, summary=dict(checks=len(checks), passed=len(checks)-len(failed), failed=failed),
                  verification=dict(source_sha256=sha(__file__), dependencies=hashes))
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# G2 scale anchor: what would make the energy unit physical?", "",
             f"Verification: **{len(checks)-len(failed)}/{len(checks)} checks passed**. No measured masses, detector data or physical energy-unit assignment is used.", "",
             "## Scale is not determined by dimensionless event fractions", "",
             "Let Lambda=1/R_phys, with the existing dimensionless R=1. The frozen spectrum is m_n^2=Lambda^2[0.25+(n+0.25)^2], mu_j^2=Lambda^2[25+j^2]. M5, MD, p, Omega and all momentum resolutions scale with Lambda; R and kappa5=2 pi Rg scale with Lambda^-1. The dimensionless g, alpha and epsilon are unchanged.", "",
             "The local kernel Kq has dimension energy^-2 and scales with Lambda^-2. The work per sideband event q Omega scales with Lambda, while the work-weighted rate coefficient sum w Kq q Omega scales with Lambda^-1. Conditional fractions cancel the common rate scale. Under x=Lambda xhat, rho_Lambda(x)=rho_hat(x/Lambda)/Lambda, so integral min(rho_+,rho_-) dx is invariant if detector resolution and the classifier threshold also scale.", "",
             "Consequently dimensionless data alone have dP/dlog Lambda=0 and zero multinomial Fisher information for log Lambda. A calibrated drive frequency or momentum, a known mass, or an absolute rate with controlled coupling/luminosity is an external dimensional anchor; none is present in this test.", "",
             "| Synthetic rescaling | Open channels | Total K | Mean pump work/event | Conditional Bayes error |",
             "|---:|---:|---:|---:|---:|"]
    for row in scaled:
        lines.append(f"| {row['synthetic_scale']:g} | {row['open_channel_count']} | {row['total_rate_coefficient']:.8g} | {row['mean_pump_work_per_scattering']:.8g} | {100*row['classification']['bayes_error']:.5f}% |")
    lines += ["", "The numbers 0.1, 1 and 10 are synthetic unit-rescalings, not GeV values or recommended energy scales. Holding an actual detector resolution fixed would not leave the classifier error invariant.", "",
              "## Minimal and stronger anchor contracts", "",
              "A single independently measured and identified mass can set Lambda=m_calibrated/m_hat **only after** freezing the dimensionless ratios. This calibrates a free scale; it does not predict the anchor mass. Likewise, inserting a known drive frequency into the assumed Omega*R=.2 card only defines a candidate radius; that arbitrary card is not a derived radius or a measured spectral spacing.", "",
              "A stronger test uses three assigned consecutive SIGNED mode masses, not three masses sorted by size. For y_n=m_n^2 and n=-1,0,+1:", "",
              "A=(y_1-2y_0+y_-1)/2, B=(y_1-y_-1)/2, C=y_0; R=1/sqrt(A), alpha=B/(2A), M5^2=C-B^2/(4A).", "",
              "Require A>0 and M5^2>=0 for the declared positive-bulk-mass card. A fourth signed n=2 anchor must satisfy y_2-3y_1+3y_0-y_-1=0. Thus three points determine parameters, while a fourth provides a falsifiable closure test.", "",
              "The synthetic baseline reconstructs A=Lambda^2, alpha=.25 and M5^2=.25 Lambda^2. Deliberately increasing only the fourth squared mass by .07 Lambda^2 violates closure and is rejected. No experimental mass has been tested or excluded.", "",
              "## What a common mass shift can and cannot hide", "",
              "A uniform positive delta M^2 changes C and M5^2 but leaves A, B, alpha and the second/third finite differences unchanged. The synthetic shift .3 Lambda^2 verifies this algebra. A portal producing a truly mode-independent shift therefore cannot repair a failed spectral-spacing test merely by adjusting the common mass. This script does not derive a portal or its dynamics.", "",
              "Integer relabeling alpha->alpha+k,n->n-k leaves masses invariant, so alpha is only defined modulo integers until mode conventions are fixed. Uncalibrated orientation adds a sign ambiguity. The neutral detector still has mu_j=mu_-j; assigning an energy unit does not distinguish those recoil signs.", "",
              "Next physical input: specify an externally calibrated mass, drive frequency or momentum reference with a justified mode identification, or provide actual detector calibration data. Until then retain all energy requirements as dimensionless ratios. Any physical beam/EFT cutoff and detector coupling require independent justification, not a chosen numerical unit.", ""]
    OUT.with_suffix(".md").write_text("\n".join(lines))
    print(json.dumps(result["summary"]))
    print(json.dumps([dict(scale=row["synthetic_scale"], channels=row["open_channel_count"], K=row["total_rate_coefficient"], work=row["mean_pump_work_per_scattering"], error=row["classification"]["bayes_error"]) for row in scaled]))
    if failed:
        raise SystemExit("Failed checks: "+", ".join(failed))


if __name__ == "__main__":
    run()
