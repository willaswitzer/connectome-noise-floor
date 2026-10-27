"""Every headline number quoted in the preprint, recomputed from the committed result tables and
checked against the value in the text (rounded as printed). Run: python src/preprint_numbers.py"""
import json
from pathlib import Path

import pandas as pd

R = Path("results")


def lin_summary(f, cond, model, rung="R0"):
    d = pd.read_csv(R / f)
    return float(d[(d.condition == cond) & (d.model == model) & (d.rung == rung)].drift_over_D.iloc[0])


def lif_summary(f, cond, rung="R0", row="seed_specific"):
    d = pd.read_csv(R / f)
    return float(d[(d.condition == cond) & (d.rung == rung) & (d.row == row)].drift_over_D.iloc[0])


def numbers():
    h1 = {k: json.loads((R / f"linear/a0.5_K5/h1_{k}.json").read_text())["observed_median_rho"] for k in ("R0", "R1", "R2")}
    cs = pd.read_csv(R / "lif/primary/rho_case_study.csv")
    pred = json.loads((R / "arm2/linear/predictions.json").read_text())
    be = pd.read_csv(R / "arm2/linear/breakeven.csv")
    bec = lambda rung: float(be[(be.family == "combined") & (be.model == "rect") & (be.rung == rung)].m_star.iloc[0])
    v = pd.read_csv(R / "arm2/linear/v630_drift.csv")
    v = v[(v.condition == "v630") & (v.model == "rect") & v.perturbation.str.startswith("t:")]
    vm = lambda rung: float(v[v.rung == rung].drift.median())
    rows = [
        # (section, quantity, recomputed, value in the preprint, decimals)
        ("3.1", "Linear median rho, R0", h1["R0"], 0.52, 2),
        ("3.1", "Linear median rho, R1", h1["R1"], 0.65, 2),
        ("3.1", "Linear median rho, R2", h1["R2"], 0.59, 2),
        ("3.3", "LIF median rho, R0", cs[cs.rung == "R0"].rho.median(), 0.55, 2),
        ("3.3", "LIF median rho, R1", cs[cs.rung == "R1"].rho.median(), 0.32, 2),
        ("3.3", "LIF median rho, R2", cs[cs.rung == "R2"].rho.median(), 0.18, 2),
        ("4.1", "P5 combined drift (m = 1)", pred["P5"]["median_drift"], 0.29, 2),
        ("4.1", "D_R0 (like-for-like)", pred["P5"]["D_R0ipsi"], 0.48, 2),
        ("4.1", "P6 unproofread drift", pred["P6"]["median_drift"], 0.61, 2),
        ("4.2", "Combined break-even m*, R0", bec("R0"), 2.5, 1),
        ("4.2", "Combined break-even m*, R1", bec("R1"), 2.2, 1),
        ("4.3", "Real-proofreading drift, R0", vm("R0"), 0.008, 3),
        ("4.3", "Real-proofreading drift, R1", vm("R1"), 0.004, 3),
        ("4.3", "Real-proofreading drift, R2", vm("R2"), 0.002, 3),
        ("4.4", "Linear, FAFB-matrix flips per neuron / D_R0", lin_summary("arm2/linear/ntcm_summary.csv", "ntcm_fafb_neuron", "norect"), 0.12, 2),
        ("4.4", "LIF, FAFB-matrix flips per neuron / D_R0", lif_summary("arm2/lif/ntcm_summary.csv", "ntcm_fafb_neuron"), 0.46, 2),
        ("4.4", "LIF, FAFB-matrix flips per cell type / D_R0", lif_summary("arm2/lif/ntcm_summary.csv", "ntcm_fafb_type"), 0.64, 2),
        ("4.4", "LIF, hemibrain-matrix flips per neuron / D_R0", lif_summary("arm2/lif/ntcm_summary.csv", "ntcm_hemi_neuron"), 0.63, 2),
        ("4.4", "LIF, hemibrain-matrix flips per cell type / D_R0", lif_summary("arm2/lif/ntcm_summary.csv", "ntcm_hemi_type"), 0.75, 2),
        ("4.4", "Linear, 50% glutamate excitatory (receiver) / D_R0", lin_summary("arm2/linear/glu_summary.csv", "glu_post_f50", "norect"), 0.56, 2),
        ("4.4", "LIF, 10% glutamate excitatory (receiver) / D_R0", lif_summary("arm2/lif/glu_summary.csv", "glu_post_f10"), 0.54, 2),
        ("4.4", "LIF, 25% glutamate excitatory (sender) / D_R0", lif_summary("arm2/lif/glu_summary.csv", "glu_pre_f25"), 1.00, 2),
        ("4.4", "LIF, 25% glutamate excitatory (receiver) / D_R0", lif_summary("arm2/lif/glu_summary.csv", "glu_post_f25"), 1.40, 2),
    ]
    t = pd.DataFrame(rows, columns=["section", "quantity", "recomputed", "in_text", "decimals"])
    t["matches"] = [round(r, d) == round(x, d) for r, x, d in zip(t.recomputed, t.in_text, t.decimals)]
    return t


if __name__ == "__main__":
    t = numbers()
    with pd.option_context("display.width", 160, "display.max_colwidth", 60):
        print(t.round(4).to_string(index=False))
    print("all match" if t.matches.all() else f"MISMATCHES: {(~t.matches).sum()}")
