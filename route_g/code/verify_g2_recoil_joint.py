#!/usr/bin/env python3
"""Conditional joint event-rate predictions, importing the unchanged G2 action.

No absolute collision probability, new action, empirical fit or work reservoir.
The imported G2 main is never executed. Ordinary COM momentum and external
overlap follow exactly its declared common-overlap dilute-rate prescription.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

import verify_g2_local_conversion as g2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "g2_recoil_joint"


def conditionals(matrix):
    """P(column|row); null when the conditioning event has zero support."""
    return [(r / r.sum()).tolist() if r.sum() > 0 else None for r in matrix]


def information(matrix):
    independent = np.outer(matrix.sum(axis=1), matrix.sum(axis=0))
    nz = matrix > 0
    return float(np.sum(matrix[nz] * np.log2(matrix[nz] / independent[nz])))


def run():
    checks = []

    def check(name, passed, measured=None, tolerance=None):
        item = dict(name=name, passed=bool(passed))
        if measured is not None:
            item["measured"] = float(measured)
        if tolerance is not None:
            item["tolerance"] = float(tolerance)
        checks.append(item)

    def close(name, actual, expected, tolerance=2e-12):
        error = float(np.max(abs(np.asarray(actual) - np.asarray(expected))))
        check(name, error < tolerance, error, tolerance)

    card, packet = dict(g2.CARD), dict(g2.PACKET)
    n, pi = card["incoming_phi_n"], card["incoming_com_momentum"]
    ls = np.arange(packet["min_l"], packet["max_l"] + 1)
    coefficients = np.exp(-ls ** 2 / (4 * packet["width"] ** 2) - 1j * ls * packet["theta0"])
    coefficients /= np.linalg.norm(coefficients)
    weights = abs(coefficients) ** 2
    close("finite_packet:normalization", weights.sum(), 1, 2e-14)
    rows = []
    for l, weight in zip(ls, weights):
        for old in g2.enumerate_channels(n, int(l), pi)["open_channels"]:
            rows.append(dict(old, packet_weight=float(weight),
                             weighted_rate=float(weight * old["rate_coefficient"]),
                             r=abs(old["j"]), delta_j=old["j"] - int(l)))
    check("enumeration:unique_signed_pairs", len({(r["m"], r["j"]) for r in rows}) == len(rows))
    check("enumeration:one_elastic_row_per_l", sum(r["m"] == n for r in rows) == len(ls))
    ms, js, rs = [sorted({r[key] for r in rows}) for key in ("m", "j", "r")]
    mi, ji, ri = [{v: i for i, v in enumerate(labels)} for labels in (ms, js, rs)]

    # The contact vertex gives dOmega/(4 pi). Quadrature is a normalization
    # regression, not an angular acceptance or finite collision-wavepacket model.
    cs, ws = np.polynomial.legendre.leggauss(8)
    phis = 2 * np.pi * np.arange(16) / 16
    directions = np.array([[np.sqrt(1 - c*c)*np.cos(p), np.sqrt(1 - c*c)*np.sin(p), c]
                           for c in cs for p in phis])
    angular_weights = np.repeat(ws / 32, 16)
    close("angular:normalization", angular_weights.sum(), 1)
    close("angular:first_moment", angular_weights @ directions, np.zeros(3))
    close("angular:second_moment", directions.T @ (angular_weights[:, None] * directions), np.eye(3)/3)
    for row in rows:
        tag = f"event:l={row['l']}:m={row['m']}"
        check(tag + ":compact_transfer", row["delta_j"] == n - row["m"])
        close(tag + ":energy", row["final_phi_energy"] + row["final_detector_energy"], row["sqrt_s"])
        pvec = row["outgoing_momentum"] * directions
        close(tag + ":back_to_back", pvec + (-pvec), np.zeros_like(pvec))
        close(tag + ":phi_mass_shell", row["final_phi_energy"]**2 - row["outgoing_momentum"]**2, row["final_phi_mass"]**2)
        close(tag + ":X_mass_shell", row["final_detector_energy"]**2 - row["outgoing_momentum"]**2, row["final_detector_mass"]**2)
        row["angular_momentum_mean"] = [0.0, 0.0, 0.0]
        row["angular_phi_covariance"] = (np.eye(3)*row["outgoing_momentum"]**2/3).tolist()
        row["angular_phi_X_cross_covariance"] = (-np.eye(3)*row["outgoing_momentum"]**2/3).tolist()

    ensembles = {}
    for name, selected in (("all_scattering", rows), ("mode_conversion", [r for r in rows if r["m"] != n])):
        total = float(sum(r["weighted_rate"] for r in selected))
        ps = np.array([r["weighted_rate"] / total for r in selected])
        joint = np.zeros((len(ms), len(js)))
        coarse = np.zeros((len(ms), len(rs)))
        selected_l = np.zeros(len(ls))
        for row, p in zip(selected, ps):
            joint[mi[row["m"]], ji[row["j"]]] += p
            coarse[mi[row["m"]], ri[row["r"]]] += p
            selected_l[row["l"] - packet["min_l"]] += p
            row["probability_conditional_" + name] = float(p)
        close(name + ":joint_normalization", joint.sum(), 1)
        check(name + ":nonnegative", joint.min() >= 0, joint.min())
        close(name + ":coarse_normalization", coarse.sum(), 1)
        close(name + ":unchanged_mode_marginal", joint.sum(axis=1), coarse.sum(axis=1))
        close(name + ":selected_l_normalization", selected_l.sum(), 1)
        cond = {}
        for key, matrix in (("j_given_m", joint), ("m_given_j", joint.T),
                            ("r_given_m", coarse), ("m_given_r", coarse.T)):
            cond[key] = conditionals(matrix)
            for index, distribution in enumerate(cond[key]):
                if distribution is not None:
                    close(name + f":{key}:{index}:normalization", sum(distribution), 1)
                    close(name + f":{key}:{index}:Bayes", np.asarray(distribution)*matrix[index].sum(), matrix[index])
        values = np.array([[r["m"], r["j"], r["l"], r["delta_j"]] for r in selected])
        mean = ps @ values
        centered = values - mean
        cov = centered.T @ (ps[:, None]*centered)
        close(name + ":mean_j_identity", mean[1], mean[2] + n - mean[0])
        close(name + ":mean_transfer", mean[3], n - mean[0])
        close(name + ":variance_j_identity", cov[1, 1], cov[2, 2] + cov[0, 0] - 2*cov[0, 2])
        close(name + ":covariance_m_j_identity", cov[0, 1], cov[0, 2] - cov[0, 0])
        close(name + ":variance_transfer", cov[3, 3], cov[0, 0])
        close(name + ":covariance_mode_transfer", cov[0, 3], -cov[0, 0])
        correlation = cov[0, 3] / np.sqrt(cov[0, 0] * cov[3, 3])
        close(name + ":transfer_anticorrelation", correlation, -1)
        check(name + ":covariance_PSD", np.linalg.eigvalsh(cov).min() > -2e-12,
              np.linalg.eigvalsh(cov).min(), 2e-12)
        signed_mi, coarse_mi = information(joint), information(coarse)
        check(name + ":data_processing_abs_j", coarse_mi <= signed_mi + 2e-12, signed_mi-coarse_mi)
        posteriors = []
        for m in ms:
            for r in rs:
                events = sorted([e for e in selected if e["m"] == m and e["r"] == r], key=lambda e: e["j"])
                if not events:
                    continue
                denominator = sum(e["weighted_rate"] for e in events)
                posterior = [dict(l=e["l"], j=e["j"], probability=e["weighted_rate"]/denominator,
                                  outgoing_momentum=e["outgoing_momentum"]) for e in events]
                close(name + f":sign_posterior:m={m}:r={r}", sum(e["probability"] for e in posterior), 1)
                gap = max(e["outgoing_momentum"] for e in events)-min(e["outgoing_momentum"] for e in events)
                posteriors.append(dict(m=m, r=r, ambiguous_sign=len(events)>1,
                                       signed_recoil_posterior=posterior, momentum_separation=gap))
        ensembles[name] = dict(
            condition="scattering event" if name == "all_scattering" else "mode-changing scattering event",
            total_rate_coefficient=total, event_rows=len(selected),
            m_labels=ms, j_labels=js, r_labels=rs, l_labels=ls.tolist(),
            signed_joint_probability=joint.tolist(), abs_recoil_joint_probability=coarse.tolist(),
            marginal_m=joint.sum(axis=1).tolist(), marginal_j=joint.sum(axis=0).tolist(),
            marginal_r=coarse.sum(axis=0).tolist(), selected_initial_l=selected_l.tolist(),
            conditionals=cond, signed_recoil_posteriors_given_m_r=posteriors,
            moments=dict(labels=["m", "j", "l", "delta_j"], mean=mean.tolist(), covariance=cov.tolist(),
                         correlation_mode_transfer=float(correlation)),
            mutual_information_bits=dict(m_and_signed_j=signed_mi, m_and_abs_j=coarse_mi,
                                         loss_from_sign_coarse_graining=signed_mi-coarse_mi))

    prior_path = ROOT / "output" / "g2_local_conversion.json"
    prior = json.loads(prior_path.read_text())
    source_path = Path(g2.__file__).resolve()
    source_hash = hashlib.sha256(source_path.read_bytes()).hexdigest()
    check("frozen_G2:artifact_source_hash", prior["source_sha256"] == source_hash)
    close("frozen_G2:all_rate", ensembles["all_scattering"]["total_rate_coefficient"],
          prior["packet"]["total_rate_coefficient"], 1e-18)
    close("frozen_G2:conversion_rate", ensembles["mode_conversion"]["total_rate_coefficient"],
          prior["packet"]["all_mode_changing_rate_coefficient"], 1e-18)
    target = card["target_phi_m"]
    close("frozen_G2:target_rate", sum(r["weighted_rate"] for r in rows if r["m"] == target),
          prior["packet"]["target_rate_coefficient"], 1e-18)
    transformed = {}
    gauge_card = dict(card, alpha=card["alpha"]+1, incoming_phi_n=n-1, target_phi_m=target-1)
    for l, weight in zip(ls, weights):
        for r in g2.enumerate_channels(n-1, int(l), pi, gauge_card)["open_channels"]:
            transformed[(r["m"]+1, r["j"])] = float(weight*r["rate_coefficient"])
    original = {(r["m"], r["j"]): r["weighted_rate"] for r in rows}
    check("large_gauge:joint_support", transformed.keys() == original.keys())
    close("large_gauge:joint_rates", [transformed[k] for k in sorted(original)],
          [original[k] for k in sorted(original)], 1e-18)
    for name, phase in (("translated_center", -1.3*ls), ("arbitrary_phase", .31*ls**2+.22*ls**3)):
        new_weights = abs(coefficients*np.exp(1j*phase))**2
        rates = [new_weights[r["l"]-packet["min_l"]]*r["rate_coefficient"] for r in rows]
        close(name + ":joint_rates", rates, [r["weighted_rate"] for r in rows], 1e-18)
    rho, dephased = np.outer(coefficients, coefficients.conj()), np.diag(weights)
    pure_rates, mixed_rates = [], []
    for r in rows:
        effect = np.zeros((len(ls), len(ls)))
        i = r["l"] - packet["min_l"]
        effect[i, i] = r["rate_coefficient"]
        pure_rates.append(float(np.trace(effect@rho).real))
        mixed_rates.append(float(np.trace(effect@dephased).real))
    close("dephasing:joint_rates", pure_rates, mixed_rates, 1e-18)

    reconstructions = []
    for row in rows:
        # This is inference within the model, not an independent verification
        # of compact-momentum conservation, which is used in the reconstruction.
        m, r, pf = row["m"], row["r"], row["outgoing_momentum"]
        sqrt_s = np.hypot(g2.phi_mass(m), pf)+np.hypot(g2.detector_mass(r), pf)
        EX = sqrt_s-np.hypot(g2.phi_mass(n), pi)
        l2 = card["R"]**2*(EX**2-card["MD"]**2-pi**2)
        d = m-n
        item = dict(m=m, r=r, outgoing_momentum=pf, inferred_sqrt_s=float(sqrt_s),
                    inferred_incoming_X_energy=float(EX), inferred_l_squared=float(l2))
        close(f"inference:l={row['l']}:m={m}:l_squared", l2, row["l"]**2)
        if d:
            jhat = (l2-d**2-r**2)/(2*d)
            item.update(inferred_j=float(jhat), inferred_l=float(jhat+d), sign_status="unique conditional model inference")
            close(f"inference:l={row['l']}:m={m}:j", jhat, row["j"])
            close(f"inference:l={row['l']}:m={m}:l", jhat+d, row["l"])
        else:
            item.update(inferred_j=0 if r == 0 else None, inferred_l=0 if r == 0 else None,
                        sign_status="trivial j=l=0" if r == 0 else "elastic sign unresolved")
        reconstructions.append(item)
    for name, data in ensembles.items():
        for row in data["signed_recoil_posteriors_given_m_r"]:
            if row["ambiguous_sign"]:
                if row["m"] != n:
                    check(f"{name}:conversion_pair:m={row['m']}:r={row['r']}:distinct_pf",
                          row["momentum_separation"] > 1e-10, row["momentum_separation"])
                else:
                    close(f"{name}:elastic_pair:r={row['r']}:same_pf", row["momentum_separation"], 0)
    example = next(r for r in ensembles["mode_conversion"]["signed_recoil_posteriors_given_m_r"]
                   if r["m"] == 0 and r["r"] == 1)
    check("example:m0_r1:both_signs", {e["j"] for e in example["signed_recoil_posterior"]} == {-1, 1})
    ambiguous_conversion = [r for r in ensembles["mode_conversion"]["signed_recoil_posteriors_given_m_r"] if r["ambiguous_sign"]]
    fraction = ensembles["mode_conversion"]["total_rate_coefficient"]/ensembles["all_scattering"]["total_rate_coefficient"]
    result = dict(
        status="bounded tree-level conditional event-rate prediction, not absolute collision probabilities",
        card=card, packet=dict(packet, l_labels=ls.tolist(), incident_weights=weights.tolist()),
        assumptions=[
            "Unchanged G2 local action, exact finite prepared packet and tree contact amplitude.",
            "Common incoming COM momentum magnitude p_i and common dilute external-space overlap across l.",
            "K=sigma*v_Moller, weighted by |c_l|^2. All-scattering and mode-conversion ensembles normalized separately.",
            "Elastic m=n denotes scattered events only; unscattered identity probability is absent.",
            "Signed j is theoretical. Neutral-X mass determines only |j|; mass readout sums the orthogonal degenerate subspace.",
            "Exact outgoing COM momentum infers signed conversion labels only under the known shared spectrum and common p_i.",
            "Inference uses compact conservation and is not an independent test of that premise.",
            "No external-space collision packet, luminosity, duration, efficiency or absolute probability is specified."],
        equations=dict(
            weight="W_mj=|c_(m+j-n)|^2 K_(m+j-n,m) on open channels",
            joint="P_A=W/K_A; P_C=1_(m!=n) W/K_C",
            transfer="j-l=n-m; l=m+j-n",
            differential="dP_E(m,j,Omega)=P_E(m,j)dOmega/(4pi); p_Phi=+p_f u, p_X=-p_f u",
            coarse="P_E(m,r)=sum_(j:|j|=r)P_E(m,j), r=0 counted once",
            infer_s="S_f=sqrt(m_m^2+p_f^2)+sqrt(MD^2+r^2/R^2+p_f^2)",
            infer_l2="L^2=R^2[(S_f-sqrt(m_n^2+p_i^2))^2-MD^2-p_i^2]",
            infer_j="d=m-n !=0: j=(L^2-d^2-r^2)/(2d), l=j+d",
            means="<j>=<l>+n-<m> in the same selected-event ensemble",
            covariance="Var(j)=Var(l)+Var(m)-2Cov(l,m); Cov(m,j)=Cov(m,l)-Var(m)"),
        events=rows, ensembles=ensembles, conversion_fraction_of_scattering_rates=fraction,
        ideal_kinematic_reconstructions=reconstructions, example_m0_absj1=example,
        conversion_sign_pairs=dict(count=len(ambiguous_conversion),
                                   minimum_momentum_separation=min(r["momentum_separation"] for r in ambiguous_conversion)),
        coherence_limit="All inclusive joint rates are unchanged by packet phases, translation or dephasing; these observables do not certify localization or determine outgoing coherences.",
        verification=dict(angular_quadrature="8 Gauss-Legendre points in cos(theta), 16 uniform azimuths",
                          g2_source_sha256=source_hash,
                          g2_artifact_sha256=hashlib.sha256(prior_path.read_bytes()).hexdigest(),
                          source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        checks=checks, summary=dict(checks=len(checks), passed=sum(c["passed"] for c in checks),
                                   failed=[c["name"] for c in checks if not c["passed"]]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.with_suffix(".json").write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    write_report(result)
    print(json.dumps(result["summary"], indent=2))
    print(json.dumps(dict(conversion_fraction=fraction, example=example,
                          all_mode_marginal=ensembles["all_scattering"]["marginal_m"],
                          conversion_mode_marginal=ensembles["mode_conversion"]["marginal_m"],
                          conversion_sign_pairs=result["conversion_sign_pairs"]), indent=2))
    if result["summary"]["failed"]:
        raise SystemExit(1)


def table(lines, matrix, row_labels, column_labels, title, column_name):
    lines += ["", "### "+title, "", "| m / "+column_name+" | "+" | ".join(map(str, column_labels))+" | sum |",
              "|---:"+"|---:"*(len(column_labels)+1)+"|"]
    matrix = np.asarray(matrix)
    for label, row in zip(row_labels, matrix):
        lines.append("| "+str(label)+" | "+" | ".join(f"{v:.8g}" for v in row)+f" | {row.sum():.8g} |")
    lines.append("| sum | "+" | ".join(f"{v:.8g}" for v in matrix.sum(axis=0))+f" | {matrix.sum():.8g} |")


def write_report(result):
    A, C = (result["ensembles"][key] for key in ("all_scattering", "mode_conversion"))
    summary = result["summary"]
    lines = [
        "# Route G G2 — joint output mode and detector recoil", "",
        "A conditional event-rate prediction from the unchanged G2 tree action and preparation. No empirical masses, new parameters, time-dependent background or absolute collision probability enter.", "",
        "## Preparation and normalization", "",
        "R=1, M5=0.5, alpha=0.25, MD=5, g=0.05, initial Phi mode n=1 and common initial COM momentum p_i=0.2. The exact finite detector packet has l=-4,...,4, width 1.2, center 0.7. Values use engineering units, not experimental measurements.", "",
        r"$$W_{mj}=|c_{m+j-n}|^2 K_{m+j-n,m},\quad K_{lm}=\sigma_{lm}v_{lm},\quad j-l=n-m.$$", "",
        r"$$P_A(m,j)=W_{mj}/K_A,\quad K_A=\sum_{m,j}W_{mj};\qquad P_C(m,j)=\mathbf1_{m\ne n}W_{mj}/K_C,\quad K_C=\sum_{m\ne n,j}W_{mj}.$$", "",
        "A conditions on scattering; C conditions on mode-changing scattering. Both use the declared common-overlap dilute rate prescription. Other flux/overlap preparations can change these weights. Scattered elastic events m=n are included only in A; neither table includes the unscattered identity amplitude.", "",
        "| Quantity | Value |", "|---|---:|", f"| All open signed rows | {len(result['events'])} |",
        f"| Mode-changing rows | {C['event_rows']} |",
        f"| K_A (mass unit^-2) | {A['total_rate_coefficient']:.12g} |",
        f"| K_C (mass unit^-2) | {C['total_rate_coefficient']:.12g} |",
        f"| K_C/K_A | {result['conversion_fraction_of_scattering_rates']:.12g} |", "",
        "All open outputs are enumerated, not only the previous m=0 example. Signed j is a theoretical label; neutral detector mass alone identifies only r=|j|.",
    ]
    for name, data in result["ensembles"].items():
        lines += ["", "## "+name.replace("_", " ").title(), ""]
        table(lines, data["signed_joint_probability"], data["m_labels"], data["j_labels"], "Signed joint P(m,j)", "j")
        table(lines, data["abs_recoil_joint_probability"], data["m_labels"], data["r_labels"], "Mass-readout joint P(m,|j|)", "r")
        lines += ["", "| Event moment | Value |", "|---|---:|"]
        for label, value in zip(data["moments"]["labels"], data["moments"]["mean"]):
            lines.append(f"| mean {label} | {value:.12g} |")
        for key, value in data["mutual_information_bits"].items():
            lines.append(f"| {key}, bits | {value:.12g} |")
        lines += ["", "Mutual information quantifies correlation in this conditional event ensemble, not confidence in the theory. JSON provides all marginals, Bayes conditionals and covariance matrices.", "",
                  "| Initial l | Incident weight | Weight conditioned on this ensemble |", "|---:|---:|---:|"]
        for l, prior, selected in zip(data["l_labels"], result["packet"]["incident_weights"], data["selected_initial_l"]):
            lines.append(f"| {l} | {prior:.10g} | {selected:.10g} |")
    lines += [
        "", "## Conservation, event selection and the differential prediction", "",
        r"$$\langle j\rangle_E=\langle l\rangle_E+n-\langle m\rangle_E,\quad \mathrm{Var}_E(j)=\mathrm{Var}_E(l)+\mathrm{Var}_E(m)-2\mathrm{Cov}_E(l,m).$$", "",
        r"$$\mathrm{Cov}_E(m,j)=\mathrm{Cov}_E(m,l)-\mathrm{Var}_E(m),\quad j-l=n-m.$$", "",
        "Both sides must use the same selected-event ensemble E. The incident mean l is zero, but selection by channel rates biases l. Raw m-versus-j covariance is therefore not the transfer law. The transfer j-l and m have perfect anticorrelation whenever their variance is nonzero.", "",
        r"$$dP_E(m,j,\Omega)=P_E(m,j)\frac{d\Omega}{4\pi},\quad\mathbf p_\Phi=p_{f,lm}\hat u,\quad\mathbf p_X=-p_{f,lm}\hat u.$$", "",
        "Each JSON event records p_f, both outgoing energies, incoming sqrt(s_l), incident weight, channel rate and both ensemble fractions. The contact amplitude is isotropic. Mean three-momenta vanish, each covariance is p_f^2 times the unit matrix / 3, and the cross-covariance is its negative. The discrete and angular laws are normalized separately.", "",
        "## Mass-readout ambiguity and ideal kinematic sign inference", "",
        "Neutral X has mu_j=mu_-j. A mass measurement projects onto the full orthogonal degenerate subspace and sums probabilities; it does not implement a coherent sign-mixing readout. Thus measuring m and detector mass does not generally identify the sign. Additional ordinary COM momentum information can distinguish conversion branches.", "",
        r"$$S_f=\sqrt{m_m^2+p_f^2}+\sqrt{M_D^2+r^2/R^2+p_f^2},\quad L^2=R^2\left[\left(S_f-\sqrt{m_n^2+p_i^2}\right)^2-M_D^2-p_i^2\right].$$", "",
        r"$$d=m-n\ne0:\qquad j=\frac{L^2-d^2-r^2}{2d},\quad l=j+d.$$", "",
        "The final identity follows from l=j+d and j^2=r^2. It reconstructs each conversion row at ideal precision. This is model-dependent inference using the known common incoming COM momentum and shared spectrum. It uses compact conservation, so reconstruction is not an independent experimental verification of that premise. The predicted outgoing momentum lines, conditional fractions and isotropy are testable. Finite resolution or unknown incoming momentum can restore ambiguity; elastic d=0 remains sign-degenerate except r=0.", "",
        "### Concrete example: m=0 and |j|=1", "",
        "| Incoming l | Outgoing j | Conditional sign weight given (m=0,r=1) | Outgoing p_f |", "|---:|---:|---:|---:|",
    ]
    example = result["example_m0_absj1"]
    for row in example["signed_recoil_posterior"]:
        lines.append(f"| {row['l']} | {row['j']} | {row['probability']:.12g} | {row['outgoing_momentum']:.12g} |")
    lines += ["", f"The lines differ by **{example['momentum_separation']:.12g}** in the common mass unit. Sign weights agree when conditioned from A or C, because m=0 is already mode-changing. These are ideal lines, not detector-resolution or event-count forecasts.", "",
              "### All ambiguous mass-readout pairs", "",
              "| m | r | Possible j | p_f separation | Ideal-p_f result |", "|---:|---:|---|---:|---|"]
    for row in A["signed_recoil_posteriors_given_m_r"]:
        if row["ambiguous_sign"]:
            signs = ", ".join(str(p["j"]) for p in row["signed_recoil_posterior"])
            verdict = "resolved for declared conversion" if row["m"] != result["card"]["incoming_phi_n"] else "elastic sign unresolved"
            lines.append(f"| {row['m']} | {row['r']} | {signs} | {row['momentum_separation']:.12g} | {verdict} |")
    pairs = result["conversion_sign_pairs"]
    lines += [
        "", f"Among the **{pairs['count']} nonelastic sign pairs**, the smallest momentum separation is **{pairs['minimum_momentum_separation']:.12g}**. Elastic sign pairs have zero separation and are excluded from that minimum.", "",
        "## Limits and next physical choice", "",
        "Translation, arbitrary packet phases and dephasing preserve every inclusive joint rate here: each signed output fixes its incoming l. These observables test the spectrum/interaction but do not certify coherent localization or fix off-diagonal recoil coherence. A coherence-sensitive readout would be a separate experiment.", "",
        "Next choose finite-resolution discrimination or an absolute encounter probability. Discrimination needs a declared detector response and incoming momentum spread. Absolute probabilities additionally need ordinary-space collision packets or luminosity/overlap and duration. Neither is silently supplied. A powered time-dependent background remains unused.", "",
        "## Verification", "",
        f"**{summary['passed']}/{summary['checks']} checks passed.** This validates the implementation, not the empirical truth of this model.", "",
        "Checks cover the frozen G2 totals/hash, joint normalization and Bayes conditionals, coarse readout, event conservation, on-shell/back-to-back recoil, isotropic angular moments, selected mean/covariance identities, full-joint large-gauge relabelling, phase/translation/dephasing invariance, ideal model-based reconstruction and information loss under |j| coarse graining.", "",
        "Reproduce: python3 -B route_g/code/verify_g2_recoil_joint.py. JSON keeps every event, full-precision results, checks and source hashes. This imports frozen G2 functions without executing or rewriting their main.", "",
    ]
    OUT.with_suffix(".md").write_text("\n".join(lines))


if __name__ == "__main__":
    run()
