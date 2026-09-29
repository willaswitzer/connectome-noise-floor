"""EXPLORATORY E1 (not pre-registered; planned 2026-09-29 after H2/H3 were seen).

Renormalized edge filter: the H2 filter (A6.1) with input fractions recomputed on the filtered
graph (condition suffix _efrn). Separates "which edges" from "how much total drive".

For each model (rectified base, unrectified base) and primary rung (R0, R1, R2; A2.2 forms):
per drawn target, excess rho = rho - that condition's H1 null median for the rung (as in H2).
Comparisons, paired over drawn targets non-NA in both conditions:
    efrn - base   (does the renormalized filter still lower agreement?)
    efrn - ef     (how much of the H2 drop is lost drive rather than lost edges?)
Reported: n, median paired difference in excess rho with a 95% bootstrap CI (1,000 resamples,
default_rng(20260928)), raw rho medians, two-sided Wilcoxon p. No multiplicity correction and
no hypothesis test: descriptive only.

Usage: python src/e1_compare.py   -> results/linear/exploratory/e1.csv"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

R = Path("results/linear")
OUT = R / "exploratory"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20260928)


def excess(cond, rung):
    d = pd.read_csv(R / cond / f"rho_{rung}.csv")
    null = json.loads((R / cond / f"h1_{rung}.json").read_text())["null_median_q025_q50_q975"][1]
    s = d[d["set"] == "draw"].set_index("perturbation")["rho"]
    return s - null, s, null


rows = []
for base in ["a0.5_K5", "a0.5_K5_norect"]:
    for cmp_name, a, b in [("efrn-base", f"{base}_efrn", base), ("efrn-ef", f"{base}_efrn", f"{base}_ef")]:
        for rung in ["R0", "R1", "R2"]:
            ea, ra, na = excess(a, rung)
            eb, rb, nb = excess(b, rung)
            j = pd.concat([ea.rename("a"), eb.rename("b")], axis=1).dropna()
            diff = (j["a"] - j["b"]).to_numpy()
            boot = [np.median(rng.choice(diff, len(diff))) for _ in range(1000)]
            rows.append({"analysis": "E1 (exploratory)", "model": "unrectified" if "norect" in base else "rectified",
                         "comparison": cmp_name, "rung": rung, "n_targets": len(j),
                         "median_diff_excess": float(np.median(diff)),
                         "ci95_low": float(np.quantile(boot, .025)), "ci95_high": float(np.quantile(boot, .975)),
                         "null_median_a": na, "null_median_b": nb,
                         "median_raw_rho_a": float(ra.loc[j.index].median()),
                         "median_raw_rho_b": float(rb.loc[j.index].median()),
                         "p_two_sided_uncorrected": float(wilcoxon(j["a"], j["b"]).pvalue)})
out = pd.DataFrame(rows)
out.to_csv(OUT / "e1.csv", index=False)
print(out.to_string(index=False))
