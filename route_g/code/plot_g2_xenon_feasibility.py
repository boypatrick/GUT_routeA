#!/usr/bin/env python3
"""Plot the bounded G2-X feasibility test; never assume an incident flux.

The dense mesh only renders the verified analytic coupling/flux formulas.
It is not an additional Gaussian scan or a detector/source optimization.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/route_g_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np

from verify_g2_xenon_feasibility import CARD, tower_width_unit

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "output" / "g2_xenon_feasibility.json"
OUT = ROOT / "output" / "figures" / "g2_xenon_feasibility"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    d = json.loads(INPUT.read_text())
    assert not d["summary"]["failed"]
    assert d["card"] == CARD
    norm = d["normalization"]["counts_per_flux_at_lambda1_Lambda1"]
    gamma = d["Higgs_constraint"]["Gamma_tower_max_GeV"]
    threshold = CARD["mh_GeV"] / 10
    # Resolve every tower threshold, especially the final narrow approach.
    thresholds = CARD["mh_GeV"] / (2 * np.sqrt(25 + np.arange(64)**2))
    near = thresholds[:, None] * (1 + np.array([-1e-4, -1e-7, -1e-11, 0, 1e-7]))
    x = np.unique(np.r_[np.geomspace(1, 30, 700), near.ravel(), threshold])
    x = x[(x >= 1) & (x <= 30)]
    y = np.geomspace(CARD["lambda_X_scan_min"], CARD["lambda_X_scan_cap"], 240)
    vacuum = 50 * x*x / CARD["v_GeV"]**2
    widths = np.array([tower_width_unit(float(a))["unit_width_GeV"] for a in x])
    higgs = np.full_like(x, np.inf)
    np.divide(gamma, widths, out=higgs, where=widths > 0)
    higgs = np.sqrt(higgs)
    limit = np.minimum(np.minimum(vacuum, higgs), 1)
    required_flux = 3 * x[None, :]**2 / (norm * y[:, None]**2)
    allowed = y[:, None] <= limit[None, :]
    masked = np.ma.masked_where(~allowed, np.log10(required_flux))

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelsize": 11, "axes.titlesize": 12,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(8.5, 10))
    grid = fig.add_gridspec(2, 2, width_ratios=[1, .035], height_ratios=[1.13, 1],
                           left=.115, right=.89, top=.85, bottom=.20,
                           wspace=.065, hspace=.54)
    ax = fig.add_subplot(grid[0, 0])
    cax = fig.add_subplot(grid[0, 1])
    bx = fig.add_subplot(grid[1, 0])
    fig.suptitle("Route G2-X  |  Physical feasibility, not a realized detector",
                 x=.115, y=.971, ha="left", fontsize=13, weight="bold")
    fig.text(.115, .929,
             "X-only Higgs portal  /  public XENON1T S2-only NR response\n"
             r"Declared scale: $m_X=\sqrt{26}\,\Lambda$; no measured radius or incident flux",
             fontsize=10, linespacing=1.5, va="top")

    mesh = ax.pcolormesh(x, y, masked, cmap="viridis_r", norm=Normalize(6, 17),
                         shading="auto", rasterized=True)
    # Theory failure takes visual precedence where both restrictions apply.
    ax.fill_between(x, vacuum, 1, color="#d9d9dc", zorder=3)
    ax.fill_between(x, higgs, np.minimum(vacuum, 1),
                    where=higgs < np.minimum(vacuum, 1),
                    color="#fae4cf", hatch="///", edgecolor="#b6814c", linewidth=0, zorder=3)
    ax.plot(x, vacuum, color="#41434a", lw=1.35, zorder=4)
    finite = np.isfinite(higgs)
    ax.plot(x[finite], higgs[finite], color="#b56325", lw=1.4, zorder=4)
    ax.axvline(threshold, color="#606873", ls=":", lw=1.1, zorder=5)
    ax.text(2.1, .18, "Retained-action vacuum failure", fontsize=9, color="#41434a")
    ax.text(3.1, .004, "Conditional Higgs restriction", fontsize=9, color="#8c4515")
    ax.text(14.1, 1.7e-5, r"$m_h/10$", fontsize=9, color="white")
    ax.set(xscale="log", yscale="log", xlim=(1, 30), ylim=(1e-5, 1),
           xlabel=r"Scale $\Lambda=1/R$ [GeV]", ylabel=r"Portal coupling $\lambda_X$")
    ax.xaxis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 30]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_title("A  |  Required flux vs conditional restrictions", loc="left", pad=10)
    cb = fig.colorbar(mesh, cax=cax, ticks=[6, 8, 10, 12, 14, 16, 17])
    cb.set_label(r"$\log_{10}\,\Phi_3$  [flux in cm$^{-2}$ s$^{-1}$]", fontsize=10)
    ax.legend(handles=[Patch(facecolor="#d9d9dc", label=r"Theory: $M_{D,0}^2<0$"),
                       Patch(facecolor="#fae4cf", hatch="///", edgecolor="#b6814c",
                             label=r"ATLAS 2023: $B_{\rm inv}<0.107$ (95% CL)")],
              loc="upper left", bbox_to_anchor=(-.01, -.24), ncol=2,
              frameon=False, fontsize=8.4, handlelength=1.6, columnspacing=1.6)

    spectrum = d["S2_spectrum"]
    edges = np.r_[spectrum["start_PE"], spectrum["end_PE"][-1]]
    binwidth = np.array(spectrum["end_PE"]) - np.array(spectrum["start_PE"])
    profiles = np.array(spectrum["conditional_probability_by_true_sign"])
    assert np.allclose(profiles[0], profiles[1], rtol=0, atol=1e-14)
    assert np.allclose(profiles.sum(axis=1), 1, rtol=0, atol=1e-12)
    for shape, color, ls, width, label in zip(
            profiles, ["#293748", "#c1732c"], ["-", (0, (5, 3))], [2.5, 1.8],
            [r"$j=+1$", r"$j=-1$ (same shape)"]):
        bx.stairs(shape / binwidth, edges, color=color, ls=ls, lw=width, label=label)
    bx.set(xscale="log", xlim=(edges[0], edges[-1]), ylim=(0, .00061),
           xlabel="Observed S2 [photoelectrons, PE]",
           ylabel=r"Conditional density [PE$^{-1}$]")
    bx.set_title("B  |  Real response does not distinguish the two signs", loc="left", pad=11)
    bx.xaxis.set_major_locator(FixedLocator([150.027, 300, 600, 1000, 2000, 3000]))
    bx.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "150" if v < 151 else f"{v:g}"))
    bx.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    bx.legend(frameon=False, loc="lower left", fontsize=9)
    err = 100*d["factorization"]["S2_only_bayes_error"]
    bx.text(.98, .96, f"Best S2-only error: {err:.2f}%\n"
            "= always guess the majority sign\n"
            "No gain from finer S2 resolution",
            transform=bx.transAxes, va="top", ha="right", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="#e3e5e8", boxstyle="round,pad=.45"))
    bx.grid(axis="y", alpha=.18)
    fig.text(.115, .10,
             "Color = at-target flux for 3 expected accepted events, NOT a discovery threshold.\n"
             "No colored region is certified viable: production, transport and backgrounds are missing.\n"
             "NR window 0.7–50 keV; full S2 bins 150.027–3000 PE; exposure 356770 kg day.\n"
             "Sources: XENON1T 2019 / official 2020 response; ATLAS 2023; PDG; CIAAW.\n"
             "Scalar-current/Helm tree model; calibration-informed response, not new measured events.",
             fontsize=8, color="#49515d", linespacing=1.5, va="top")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT.with_suffix(".png"), dpi=210)
    fig.savefig(OUT.with_suffix(".pdf"))
    plt.close(fig)
    metadata = dict(input=str(INPUT.relative_to(ROOT)), input_sha256=sha(INPUT),
                    source_sha256=sha(__file__), kernel_sha256=sha(Path(__file__).with_name("verify_g2_xenon_feasibility.py")),
                    note="Dense analytic rendering only; no assumed beam flux, added Gaussian scan, or observed-event fit.",
                    display_scale_points=len(x), accepted_probabilities_normalized=True,
                    sign_shapes_identical=True)
    OUT.with_suffix(".json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(f"Rendered {OUT.with_suffix('.png')} and vector PDF; source hashes saved.")


if __name__ == "__main__":
    main()
