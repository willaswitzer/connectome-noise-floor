"""alpha/K sensitivity (A6.3): collect H1 and H4 across conditions into one table (descriptive).

Usage: python src/summarize_alphak.py -> results/linear/sensitivity_alphaK.csv"""
import json
from pathlib import Path

import pandas as pd

R = Path("results/linear")
CONDS = ["a0.3_K5", "a0.7_K5", "a0.5_K3", "a0.5_K5", "a0.5_K10"]
RUNGS = ["R0", "R1", "R2", "R2right", "R0ipsi", "R1ipsi_right"]
rows = []
for c in CONDS:
    for r in RUNGS:
        h = json.loads((R / c / f"h1_{r}.json").read_text())
        rows.append({"condition": c, "test": "H1", "item": r, "n": h["n_targets"],
                     "median_rho": h["observed_median_rho"], "null_median": h["null_median_q025_q50_q975"][1],
                     "p": h["p_one_sided"]})
    h4 = pd.read_csv(R / c / "h4.csv")
    for _, x in h4.iterrows():
        rows.append({"condition": c, "test": "H4", "item": x["pair"], "n": x["n_targets"],
                     "median_excess_diff": x["median_excess_diff"], "ci95_low": x["ci95_low"],
                     "ci95_high": x["ci95_high"], "p": x["p_one_sided"], "p_holm": x["p_holm"]})
out = pd.DataFrame(rows)
out.to_csv(R / "sensitivity_alphaK.csv", index=False)
print(out.round(4).to_string(index=False))
