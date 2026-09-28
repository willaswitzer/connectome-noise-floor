"""H4: rung ordering R0 > R2 > R1 (PREREG A1.1, test defined in Amendment 5).

Like-for-like rung forms (A5.3): R0 = R0ipsi, R1 = R1ipsi_right, R2 = R2right.
Per drawn target: excess rho = rho - that rung's H1 null median (h1_<rung>.json), so each rung
is measured against its own mismatched-target baseline (the "same permutation baseline as
H1", A1.1). For each predicted pair (R0>R2, R2>R1, R0>R1): targets non-NA on both rungs,
one-sided Wilcoxon signed-rank on paired excess rho, Holm over the 3 pairs, median paired
difference with a 95% bootstrap CI over targets (1,000 resamples). Raw rho is reported too.

Usage: python src/h4_test.py [--cond a0.5_K5]"""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ap = argparse.ArgumentParser()
ap.add_argument("--cond", default="a0.5_K5")
args = ap.parse_args()
D = Path(f"results/linear/{args.cond}")
FORMS = {"R0": "R0ipsi", "R1": "R1ipsi_right", "R2": "R2right"}
PAIRS = [("R0", "R2"), ("R2", "R1"), ("R0", "R1")]
rng = np.random.default_rng(20260928)

rho, null = {}, {}
for r, f in FORMS.items():
    d = pd.read_csv(D / f"rho_{f}.csv")
    rho[r] = d[d["set"] == "draw"].set_index("perturbation")["rho"]
    null[r] = json.loads((D / f"h1_{f}.json").read_text())["null_median_q025_q50_q975"][1]

rows = []
for a, b in PAIRS:
    j = pd.concat([rho[a].rename("a"), rho[b].rename("b")], axis=1).dropna()
    ea, eb = j["a"] - null[a], j["b"] - null[b]
    diff = (ea - eb).to_numpy()
    p = float(wilcoxon(ea, eb, alternative="greater").pvalue) if len(j) else np.nan
    boot = [np.median(rng.choice(diff, len(diff))) for _ in range(1000)] if len(j) else [np.nan]
    rows.append({"pair": f"{a}>{b}", "n_targets": len(j), "null_median_a": null[a], "null_median_b": null[b],
                 "median_excess_diff": float(np.median(diff)) if len(j) else np.nan,
                 "ci95_low": float(np.quantile(boot, .025)), "ci95_high": float(np.quantile(boot, .975)),
                 "median_raw_rho_a": float(j["a"].median()), "median_raw_rho_b": float(j["b"].median()),
                 "p_one_sided": p})
out = pd.DataFrame(rows)
order = out["p_one_sided"].sort_values().index
m, running = len(out), 0.0
for rank, i in enumerate(order):
    running = max(running, min(1.0, (m - rank) * out.loc[i, "p_one_sided"]))
    out.loc[i, "p_holm"] = running
out.to_csv(D / "h4.csv", index=False)
print(out.to_string(index=False))
