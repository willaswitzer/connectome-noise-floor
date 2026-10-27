"""Preprint figures, built only from committed result tables (no model runs).

Fig 1  Arm 1: cross-animal agreement of linear predictions per rung, vs the mismatched-target
       null; and linear vs LIF agreement on the case-study set.
Fig 2  Arm 2: prediction drift as a fraction of the left-right disagreement, by error type and
       multiple of its anchor; linear (rectified; NT flips unrectified) and LIF (R0, seed-specific).
Fig 4  Transmitter signs: classifier-rate errors (Fig. 2A confusion matrices) and glutamate made
       excitatory, linear (unrectified) vs LIF, R0.
Fig 3  Real proofreading (FlyWire v630 -> v783): drift per target vs the rung disagreement.

Usage: python src/make_figures.py  ->  results/figures/fig{1..4}.{png,pdf}
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = Path("results")
OUT = R / "figures"

# reference palette (light), validated: dataviz validate_palette.js
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
LINEAR, LIF = "#2a78d6", "#eb6834"
TYPE_COLOR = {"combined": "#2a78d6", "spur": "#eb6834", "flip": "#1baf7a",
              "merge": "#eda100", "split": "#e87ba4", "miss": "#008300"}
TYPE_MARK = {"combined": "o", "spur": "s", "flip": "D", "merge": "^", "split": "v", "miss": "P"}
TYPE_NAME = {"combined": "Combined", "spur": "Spurious synapses", "flip": "Transmitter sign flips",
             "merge": "False merges", "split": "False splits", "miss": "Missed synapses"}
RUNG_NAME = {"R0": "Left vs. right\n(one fly)", "R1": "FlyWire vs.\nmale CNS", "R2": "FlyWire vs.\nhemibrain"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "lines.linewidth": 2, "pdf.fonttype": 42})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def breakeven_line(ax, label=True, x=None):
    ax.axhline(1, color=INK, lw=1, ls=(0, (4, 3)))
    if label:
        ax.annotate("break-even (= left–right disagreement)", xy=(1, 1), xycoords=("axes fraction", "data"),
                    xytext=(0, 3), textcoords="offset points", ha="right", va="bottom", color=INK2, fontsize=8)


# ----------------------------------------------------------------------------- Fig 1
def fig1():
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.8), gridspec_kw={"width_ratios": [1.15, 1]})
    rng = np.random.default_rng(0)
    for i, rung in enumerate(("R0", "R1", "R2")):
        r = pd.read_csv(R / f"linear/a0.5_K5/rho_{rung}.csv")
        r = r[(r.set == "draw") & r.rho.notna()]
        h = json.loads((R / f"linear/a0.5_K5/h1_{rung}.json").read_text())
        lo, mid, hi = h["null_median_q025_q50_q975"]
        a.fill_between([i - 0.38, i + 0.38], lo, hi, color=GRID, lw=0, zorder=0)
        a.plot([i - 0.38, i + 0.38], [mid, mid], color=INK2, lw=1, zorder=1)
        x = i + rng.uniform(-0.28, 0.28, len(r))
        a.scatter(x, r.rho, s=9, color=LINEAR, alpha=0.55, lw=0, zorder=2)
        med = r.rho.median()
        a.plot([i - 0.36, i + 0.36], [med, med], color=INK, lw=2.5, zorder=3, solid_capstyle="round")
        a.annotate(f"{med:.2f}", (i + 0.38, med), xytext=(3, 0), textcoords="offset points",
                   va="center", fontsize=8.5, color=INK, fontweight="bold")
        a.annotate(f"n = {len(r)}", (i, -0.37), ha="center", fontsize=7.5, color=INK2)
    a.annotate("null: mismatched\ntargets (median,\n95% band)", (2.4, 0.1), xytext=(2.43, 0.1), fontsize=7.5,
               color=INK2, ha="left", va="center")
    a.set_xticks(range(3), [RUNG_NAME[k] for k in ("R0", "R1", "R2")])
    a.set_xlim(-0.6, 3.1)
    a.set_ylim(-0.42, 1.02)
    a.set_ylabel("Spearman ρ of predicted responses")
    a.set_title("A  Linear predictions replicate across animals")
    a.grid(axis="x", visible=False)

    c = pd.read_csv(R / "lif/primary/rho_case_study.csv").dropna(subset=["rho", "rho_linear"])
    for i, rung in enumerate(("R0", "R1", "R2")):
        s = c[c.rung == rung]
        for _, row in s.iterrows():
            b.plot([i - 0.17, i + 0.17], [row.rho_linear, row.rho], color=GRID, lw=1, zorder=1)
        b.scatter(np.full(len(s), i - 0.17), s.rho_linear, s=18, color=LINEAR, zorder=2, lw=0)
        b.scatter(np.full(len(s), i + 0.17), s.rho, s=18, color=LIF, zorder=2, lw=0, marker="s")
        for xx, v, col in ((i - 0.17, s.rho_linear.median(), LINEAR), (i + 0.17, s.rho.median(), LIF)):
            b.plot([xx - 0.11, xx + 0.11], [v, v], color=INK, lw=2.5, zorder=3, solid_capstyle="round")
    lower = (c.rho < c.rho_linear).sum()
    b.text(0.02, 0.03, f"LIF below linear in {lower} of {len(c)} comparisons", transform=b.transAxes,
           fontsize=8, color=INK2)
    b.scatter([], [], s=18, color=LINEAR, label="Linear")
    b.scatter([], [], s=18, color=LIF, marker="s", label="Spiking (LIF)")
    b.legend(loc="lower center", fontsize=8, handletextpad=0.2, bbox_to_anchor=(0.5, 0.1))
    b.set_xticks(range(3), [RUNG_NAME[k] for k in ("R0", "R1", "R2")])
    b.set_ylim(-0.55, 1.02)
    b.set_ylabel("Spearman ρ")
    b.set_title("B  Spiking predictions agree less (14 case studies)")
    b.grid(axis="x", visible=False)
    fig.tight_layout(w_pad=3)
    save(fig, "fig1_noise_ceiling")


# ----------------------------------------------------------------------------- Fig 2
def fig2():
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 4.0), sharey=True,
                               gridspec_kw={"width_ratios": [1.35, 1]})
    cv = pd.read_csv(R / "arm2/linear/curves.csv")
    for fam in ("combined", "spur", "flip", "merge", "split", "miss"):
        model = "norect" if fam == "flip" else "rect"
        s = cv[(cv.family == fam) & (cv.model == model) & (cv.rung == "R0")].sort_values("m")
        y = s.median_drift / s.D_R
        a.plot(s.m, y, color=TYPE_COLOR[fam], marker=TYPE_MARK[fam], ms=5, label=TYPE_NAME[fam])
    a.set_xscale("log")
    a.set_xticks([0.5, 1, 2, 5, 10, 25, 50], ["0.5", "1", "2", "5", "10", "25", "50"])
    a.set_xlabel("Error rate (multiple of each type's anchor)")
    a.set_ylabel("Prediction drift / left–right disagreement")
    a.set_title("A  Linear model")
    a.axvline(1, color=INK2, lw=0.8, ls=":")
    a.annotate("proofread\naccuracy", (1, 1.62), xytext=(3, 0), textcoords="offset points", fontsize=7.5,
               color=INK2, va="top")
    breakeven_line(a, label=False)
    a.annotate("break-even", xy=(0.52, 1), xytext=(0, 3), textcoords="offset points", fontsize=8, color=INK2)

    lc = pd.read_csv(R / "arm2/lif/curves.csv")
    for t in ("spur", "flip", "merge", "split", "miss"):
        s = lc[(lc.type == t) & (lc.rung == "R0") & (lc.row == "seed_specific")].sort_values("m")
        b.plot(s.m, s.median_drift / s.D_R, color=TYPE_COLOR[t], marker=TYPE_MARK[t], ms=5)
    b.set_xscale("log")
    b.set_xticks([1, 3, 10], ["1", "3", "10"])
    b.set_xlim(0.8, 40)
    b.set_xlabel("Error rate (multiple of each type's anchor)")
    for t, dy in (("flip", 0), ("merge", 0.05), ("split", -0.05), ("spur", 0), ("miss", 0)):
        s = lc[(lc.type == t) & (lc.rung == "R0") & (lc.row == "seed_specific")].sort_values("m")
        yv = (s.median_drift / s.D_R).iloc[-1]
        b.annotate(TYPE_NAME[t], (s.m.iloc[-1], yv), xytext=(6, dy * 100), textcoords="offset points",
                   va="center", fontsize=7.5, color=INK)
    b.set_title("B  Spiking model (LIF), 9 case studies")
    breakeven_line(b, label=False)
    for ax in (a, b):
        ax.set_ylim(0, 1.65)
        ax.grid(axis="x", visible=False)
    # direct labels at the right end of each linear curve; legend as well (>= 2 series)
    handles, labels = a.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=8, bbox_to_anchor=(0.5, -0.04),
               handletextpad=0.4, columnspacing=1.2)
    fig.text(0.01, -0.06, "Anchors (m = 1): missed synapses, false splits and false merges 0.8% of a neuron's synapses; "
             "spurious synapses 5% added; transmitter sign flips 6% of neurons. Linear: rectified model, R0, 107 targets; "
             "sign flips unrectified, 133 targets.\nLIF: 9 seed-specific case-study perturbations, R0; sign flips seed-exempt. "
             "Break-even = median drift equals the median left–right disagreement.", fontsize=7, color=INK2, va="top")
    fig.tight_layout(rect=(0, 0.06, 1, 1), w_pad=2)
    save(fig, "fig2_error_budget")


# ----------------------------------------------------------------------------- Fig 3
def fig4():
    ln = pd.read_csv(R / "arm2/linear/ntcm_summary.csv")
    lg = pd.read_csv(R / "arm2/linear/glu_summary.csv")
    fn = pd.read_csv(R / "arm2/lif/ntcm_summary.csv")
    fg = pd.read_csv(R / "arm2/lif/glu_summary.csv")
    lin = lambda d, c: float(d[(d.condition == c) & (d.model == "norect") & (d.rung == "R0")].drift_over_D.iloc[0])
    lif = lambda d, c: float(d[(d.condition == c) & (d.rung == "R0") & (d.row == "seed_specific")].drift_over_D.iloc[0])
    rows = [
        ("Classifier errors", "FAFB matrix, per neuron (1.1%)", lin(ln, "ntcm_fafb_neuron"), lif(fn, "ntcm_fafb_neuron")),
        ("Classifier errors", "FAFB matrix, per cell type (1.2%)", lin(ln, "ntcm_fafb_type"), lif(fn, "ntcm_fafb_type")),
        ("Classifier errors", "Hemibrain matrix, per neuron (3.5%)", lin(ln, "ntcm_hemi_neuron"), lif(fn, "ntcm_hemi_neuron")),
        ("Classifier errors", "Hemibrain matrix, per cell type (3.4%)", lin(ln, "ntcm_hemi_type"), lif(fn, "ntcm_hemi_type")),
        ("Glutamate excitatory", "10%, by receiving cell type", lin(lg, "glu_post_f10"), lif(fg, "glu_post_f10")),
        ("Glutamate excitatory", "25%, by sending cell type", lin(lg, "glu_pre_f25"), lif(fg, "glu_pre_f25")),
        ("Glutamate excitatory", "25%, by receiving cell type", lin(lg, "glu_post_f25"), lif(fg, "glu_post_f25")),
    ]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    y = np.arange(len(rows))[::-1].astype(float)
    y[4:] = y[4:] - 0.0
    gap = np.array([0, 0, 0, 0, 0.6, 0.6, 0.6])
    y = np.arange(len(rows))[::-1] + gap[::-1] * 0
    y = np.array([7.6, 6.6, 5.6, 4.6, 2.6, 1.6, 0.6])
    h = 0.34
    for yy, (_, lab, vl, vf) in zip(y, rows):
        ax.barh(yy + h / 2 + 0.02, vl, height=h, color=LINEAR, lw=0)
        ax.barh(yy - h / 2 - 0.02, vf, height=h, color=LIF, lw=0)
        ax.text(vl + 0.02, yy + h / 2 + 0.02, f"{vl:.2f}", va="center", fontsize=7.5, color=INK2)
        ax.text(vf + 0.02, yy - h / 2 - 0.02, f"{vf:.2f}", va="center", fontsize=7.5, color=INK)
    ax.set_yticks(y, [r[1] for r in rows])
    ax.text(-0.01, 8.35, "Transmitter classifier errors (share of signed neurons flipped)", transform=ax.get_yaxis_transform(),
            ha="right", fontsize=8.5, fontweight="bold", color=INK)
    ax.text(-0.01, 3.35, "Glutamate made excitatory", transform=ax.get_yaxis_transform(), ha="right",
            fontsize=8.5, fontweight="bold", color=INK)
    ax.axvline(1, color=INK, lw=1, ls=(0, (4, 3)))
    ax.text(1.01, 8.3, "break-even", fontsize=8, color=INK2, va="center")
    ax.set_xlim(0, 1.62)
    ax.set_ylim(0.1, 8.7)
    ax.set_xlabel("Prediction drift / left–right disagreement (R0)")
    ax.grid(axis="y", visible=False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=LINEAR, label="Linear (unrectified)"),
                       Patch(color=LIF, label="Spiking (LIF), 9 case studies")],
              loc="center right", fontsize=8, bbox_to_anchor=(1.0, 0.385), frameon=True,
              facecolor=SURFACE, edgecolor=SURFACE, framealpha=1)
    ax.set_title("Transmitter signs: classifier errors stay below break-even; glutamate's sign may not",
                 loc="left", x=-0.62)
    fig.tight_layout()
    save(fig, "fig4_transmitter_signs")


# ----------------------------------------------------------------------------- Fig 4
def fig3():
    v = pd.read_csv(R / "arm2/linear/v630_drift.csv")
    v = v[(v.condition == "v630") & (v.model == "rect") & v.perturbation.str.startswith("t:")]
    dis = pd.read_csv(R / "arm2/linear/rung_disagreement.csv")
    dis = dis[dis.model == "rect"]
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    rng = np.random.default_rng(0)
    for i, rung in enumerate(("R0", "R1", "R2")):
        s = v[(v.rung == rung) & v.drift.notna()]
        d = dis[(dis.rung == rung) & dis.perturbation.isin(s.perturbation)].dis.dropna()
        x = i + rng.uniform(-0.25, 0.25, len(s))
        ax.scatter(x, s.drift.clip(lower=1e-4), s=9, color=LINEAR, alpha=0.6, lw=0, zorder=2)
        ax.plot([i - 0.3, i + 0.3], [s.drift.median()] * 2, color=INK, lw=2.5, solid_capstyle="round", zorder=3)
        ax.plot([i - 0.3, i + 0.3], [d.median()] * 2, color=INK, lw=1, ls=(0, (4, 3)), zorder=3)
        ax.annotate(f"{s.drift.median():.3f}", (i + 0.32, s.drift.median()), xytext=(2, 0),
                    textcoords="offset points", va="center", fontsize=8, fontweight="bold")
        ax.annotate(f"{d.median():.2f}", (i + 0.32, d.median()), xytext=(2, 0),
                    textcoords="offset points", va="center", fontsize=8, color=INK2)
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1.5)
    ax.set_xticks(range(3), [RUNG_NAME[k] for k in ("R0", "R1", "R2")])
    ax.set_xlim(-0.55, 2.75)
    ax.set_ylabel("Drift (1 − ρ), log scale")
    ax.grid(axis="x", visible=False)
    ax.plot([], [], color=INK, lw=2.5, label="median drift, v630 → v783")
    ax.plot([], [], color=INK, lw=1, ls=(0, (4, 3)), label="median disagreement on that comparison")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Real proofreading moved predictions ~2% of the noise floor")
    fig.tight_layout()
    save(fig, "fig3_real_proofreading")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print("wrote", sorted(p.name for p in OUT.glob("*.png")))
