#!/usr/bin/env python3
"""Join two bounded LO screens into one figure; no new optimization."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/p54-physics-first-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/p54-physics-first-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RF = Path(__file__).resolve().parents[1]
OUT = RF / "output/p54_physics_first_screen.json"
FIG = RF / "output/figures/p54_physics_first_feasibility"
BLUE, ORANGE, INK = "#25628A", "#B86520", "#223340"


def decode(v):
    return np.asarray(v["real"])+1j*np.asarray(v["imag"])


def metrics(c, f):
    """Distance from interval center / half-width: <=1 is the actual gate."""
    p, n = f["discrepancy_policy"], f["neutrino_targets"]
    bands = np.asarray(n["bands_3sigma"])
    rb = np.asarray(n["ratio_3sigma_rectangle"])
    return [
        float(max(abs(np.asarray(c["charged_fractional_residuals"])) /
                  np.asarray(p["charged_fractional_widths"]))),
        float(max(abs(np.asarray(c["CKM_fractional_residuals"]))) / p["CKM_fractional_width"]),
        c["CKM_delta_error_rad"] / p["CKM_delta_rad_width"],
        float(max(abs(np.asarray(c["PMNS_sin2_theta12_theta23_theta13"])-bands.mean(axis=1)) /
                  ((bands[:, 1]-bands[:, 0])/2))),
        float(abs(c["mass_squared_ratio"]-rb.mean())/((rb[1]-rb[0])/2)),
    ]


def audit(s, f):
    checks = []
    def check(name, condition):
        checks.append(dict(name=name, passed=bool(condition)))
    def close(name, a, b):
        check(name, np.allclose(a, b, rtol=2e-9, atol=1e-14))
    check("scale_subreport_checks", s["passed"] == s["total"])
    check("flavor_subreport_checks", f["summary"]["all_pass"])
    for label, hashes in (("scale", s["source_hashes"]), ("flavor", f["source_sha256"])):
        for path, expected in hashes.items():
            check(label+"_source_hash_"+path,
                  hashlib.sha256((RF/path).read_bytes()).hexdigest() == expected)
    b, ip = s["baseline"], s["input_card"]
    a = np.array([1/ip["alpha_s"], ip["alpha_em_inverse"]*ip["sin2_theta_w"],
                  .6*ip["alpha_em_inverse"]*(1-ip["sin2_theta_w"])])
    t = 2*np.pi*(a[2]-.4*a[0]-.6*a[1])/(44/5)
    close("independent_D_parity_MI", s["MZ_GeV"]*np.exp(t), b["MI_GeV"])
    close("shared_comparison_scale", f["scales"]["comparison_mu_GeV"], b["MI_GeV"])
    close("shared_sigma_not_equal_MI", f["scales"]["baseline_sigma_GeV"], b["sigma_GeV"])
    for c in s["scenarios"]:
        close(c["id"]+"_same_sigma", c["sigma_GeV"], b["sigma_GeV"])
    for key in ("best", "fixed_sigma_refinement"):
        c = f[key]
        family = c["normalized_overlap_seesaw_family"]
        h, ff = decode(c["Hprime"]), decode(c["Fprime"])
        r, ss = c["r"], complex(decode(c["s"]))
        dmin = np.linalg.norm(ff, 2)/family["f_D_norm_cap"]
        room = 1-(1+r*r)*(np.linalg.norm(h, 2)/family["h_D_norm_cap"])**2
        dmax = np.sqrt(max(room, 0)/(1+abs(r*ss)**2))
        close(key+"_independent_overlap_interval", [dmin, dmax], family["d_interval_under_caps"])
        check(key+"_plot_decisions_equal_saved_gates",
              [v <= 1 for v in metrics(c, f)] == list(c["shape_screen_gates"].values()))
        close(key+"_sigma_interval_from_d", np.array([dmin, dmax])*family["sigma_over_d_GeV"],
              family["sigma_interval_GeV_central_dm31"])
    check("no_scale_only_physical_promotion", not any(c["physical_candidate"] for c in s["scenarios"]))
    check("frozen_benchmark_not_reused", not f["scalar_benchmark_reused_for_physical_fit"])
    check("no_matched_UV_claim", not f["joint_flavor_seesaw_proton_matched_UV_witness_found"])
    return checks


def figure(s, f):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelcolor": INK, "text.color": INK,
                         "svg.hashsalt": "p54-physics-first-v1"})
    fig = plt.figure(figsize=(15, 9), facecolor="#FCFCFA")
    gs = fig.add_gridspec(2, 2, left=.11, right=.97, bottom=.20, top=.81,
                          hspace=.75, wspace=.35, width_ratios=[1.1, 1])
    axa, axb, axc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])
    fig.text(.045, .95, "P54 | Physics-first feasibility", fontsize=23, weight="bold")
    fig.text(.045, .905, "Three flavor-shape starts + one fixed-scale refinement | nine correlated threshold cards", fontsize=12)
    headline = ("A conditional flavor/scale witness exists; proton and vacuum tests remain open."
                if f["joint_LO_flavor_and_scale_witness_found"] else
                "No joint witness in this bounded search. This is not a model-wide exclusion.")
    fig.text(.045, .865, headline, fontsize=12, color=ORANGE, weight="bold")
    names = ["Charged masses", "CKM angles", "CKM phase", "PMNS angles", r"$\Delta m^2$ ratio"]
    yy = np.arange(5)
    for key, offset, color, label in (("best", -.16, BLUE, "Shape-first"),
                                      ("fixed_sigma_refinement", .16, ORANGE, "Fixed-scale retry")):
        axa.barh(yy+offset, metrics(f[key], f), height=.29, color=color, label=label)
    axa.axvline(1, color=INK, linestyle="--", lw=1)
    axa.set(yticks=yy, yticklabels=names, xlabel="Largest interval-distance ratio in group (pass: ≤ 1)")
    axa.invert_yaxis()
    axa.set_title("A  |  Flavor discrepancies", loc="left", fontsize=13, weight="bold", pad=13)
    axa.legend(loc="lower right", frameon=False, fontsize=9)
    axa.grid(axis="x", alpha=.16)
    axa.set_axisbelow(True)
    for j, (key, color) in enumerate((("best", BLUE), ("fixed_sigma_refinement", ORANGE))):
        lo, hi = np.array(f[key]["normalized_overlap_seesaw_family"]["sigma_interval_GeV_central_dm31"])/1e13
        axb.plot([lo, hi], [j, j], color=color, lw=7, solid_capstyle="butt")
        axb.plot([lo, hi], [j, j], "|", color=color, ms=15, mew=2)
    sigma = s["baseline"]["sigma_GeV"]/1e13
    axb.axvline(sigma, color=INK, ls="--", lw=1.5)
    axb.text(sigma, -.52, "Gauge-required σ", ha="right", va="bottom", fontsize=10)
    axb.set(xscale="log", ylim=(-.8, 1.55), yticks=[0, 1],
            yticklabels=["Shape-first\nprofile", "Fixed-scale\nprofile"],
            xlabel=r"$\sigma$ [$10^{13}$ GeV] (log scale)")
    axb.invert_yaxis()
    axb.set_title("B  |  Same-matrix seesaw / overlap", loc="left", fontsize=13, weight="bold", pad=13)
    axb.grid(axis="x", alpha=.16)
    axb.text(0, -.35, "One texture per line, not the model's full allowed region.\n"
             "Central atmospheric scale; Dirac-unit caps ≤ 1. Not a vacuum test.",
             transform=axb.transAxes, fontsize=9, va="top")
    colors = {1/3: BLUE, 1.: INK, 3.: ORANGE}
    styles = {1.: "-", 2.: "--", 3.: ":"}
    angle = np.linspace(0, np.pi/2, 250)
    for row in s["scenarios"]:
        ratios = row["declared_parent_mass_ratios"]
        k20 = ratios["Phi54:(20prime,1,1)"]["ratio"]
        k15 = ratios["Sigma126:(15,2,2)"]["ratio"]
        fmax = row["allowed_effective_flavor_norm"]
        label = (r"$\kappa_{20}=1/3$" if k20 < 1 else r"$\kappa_{20}=$"+f"{k20:g}") if k15 == 1 else None
        axc.plot(fmax/row["gauge_only_A_L"]*np.cos(angle),
                 fmax/row["gauge_only_A_R"]*np.sin(angle),
                 color=colors[k20], ls=styles[k15], lw=1.6, label=label)
    axc.set(aspect="equal", xlabel=r"$|F_L|$ (UV coefficient convention)",
            ylabel=r"$|F_R|$ (UV coefficient convention)", xlim=(0, 1.30), ylim=(0, 1.42))
    axc.set_title("C  |  Gauge-proton amplitude conditions", loc="left", fontsize=13, weight="bold", pad=13)
    axc.grid(alpha=.15)
    axc.legend(loc="upper right", frameon=False, title=r"$\kappa_{33}=1/\kappa_{20}$", fontsize=10)
    axc.text(.05, .945, "Interior meets this channel only.\nNo amplitudes calculated.\nNo physical point certified.",
             transform=axc.transAxes, fontsize=9,
             va="top",
             bbox=dict(facecolor="white", alpha=.88, edgecolor="none"))
    axc.text(0, -.17, r"Solid / dashed / dotted: $\kappa_{15}=1,2,3$."+"\n"
             r"$M_X=(1.32-2.62)\times10^{15}$ GeV; $M_I=5.08\times10^{13}$ GeV."+"\n"
             "Whole-parent hypotheses; no scalar exchange or other channels.",
             transform=axc.transAxes, fontsize=9, va="top")
    fig.text(.045, .045, "LO only: two-matrix/type-I boundary at MI; SM target running below MI; one-loop gauge / parent leading logs.\n"
             "No PS Yukawa running, sequential seesaw or finite matching. A uses declared widths / frozen NuFIT intervals, not statistical χ².\n"
             "Sources: saved screening JSON; Mummidi–Patel; NuFIT 6.0; Super-K (2020); Yoo et al. lattice (2021). See report for conventions.",
             fontsize=9, va="bottom", color="#53626D")
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG.with_suffix(".png"), dpi=180, facecolor=fig.get_facecolor())
    fig.savefig(FIG.with_suffix(".svg"), facecolor=fig.get_facecolor(), metadata={"Date": None})
    plt.close(fig)


def main():
    paths = [RF/"output/p54_scales_lo_screen.json", RF/"output/p54_flavor_lo_screen.json"]
    s, f = [json.loads(p.read_text()) for p in paths]
    checks = audit(s, f)
    if not all(c["passed"] for c in checks):
        raise AssertionError([c for c in checks if not c["passed"]])
    figure(s, f)
    # Disclose inherited input-card mismatch; do not disguise it as precision
    # or silently change the established low-order target tolerances.
    b = s["baseline"]
    ps = b["alpha_PS_inverse_at_MI"]
    gauge123 = np.sqrt(4*np.pi/np.array([.4*ps[0]+.6*ps[2], ps[1], ps[0]]))
    target123 = np.asarray(f["target_record"]["gauge_at_high"])
    rows = [dict(source="p54_flavor_lo_screen.json", key=key, figure_A_metrics=metrics(f[key], f),
                 figure_B_sigma_interval_GeV=f[key]["normalized_overlap_seesaw_family"]["sigma_interval_GeV_central_dm31"])
            for key in ("best", "fixed_sigma_refinement")]
    result = dict(schema="p54-physics-first-join-v1", date="2026-09-22",
                  frozen_benchmark_regression_only=True, no_full_matching_started=True,
                  joint_flavor_scale_witness=f["joint_LO_flavor_and_scale_witness_found"],
                  common_scale_but_not_harmonized_precision_input_card=dict(
                      gauge_screen_g123_at_MI=gauge123.tolist(),
                      flavor_transport_g123_at_MI=target123.tolist(),
                      relative_discrepancy=(target123/gauge123-1).tolist(),
                      new_blocker=False),
                  matched_physical_witness=False, figure_rows=rows,
                  figure_C_scenarios=[dict(id=c["id"], source="p54_scales_lo_screen.json",
                                          A_L=c["gauge_only_A_L"], A_R=c["gauge_only_A_R"],
                                          norm_limit=c["allowed_effective_flavor_norm"]) for c in s["scenarios"]],
                  checks=checks, passed=sum(c["passed"] for c in checks), total=len(checks),
                  source_hashes={str(p.relative_to(RF)):hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in paths+[Path(__file__)]})
    OUT.write_text(json.dumps(result, indent=2)+"\n")
    print(f"Join: {result['passed']}/{result['total']} checks; joint LO witness={result['joint_flavor_scale_witness']}")
    print(FIG.with_suffix(".png"))


if __name__ == "__main__":
    main()
