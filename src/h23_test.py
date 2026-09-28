"""H2 (edge filter raises replicability) and H3 (sign ablation lowers it), Amendment 6.

For each primary rung (R0, R1, R2; A2.2 forms, as in H1): per drawn target, excess rho =
rho - that condition's H1 null median for the rung. Paired over targets non-NA in both
conditions. H2: one-sided Wilcoxon, filtered > base. H3: one-sided Wilcoxon, base > ablated.
Holm over the 3 rungs within each hypothesis and model. Median paired difference with a 95%
bootstrap CI (1,000 resamples, default_rng(20260928)).

Usage: python src/h23_test.py --base a0.5_K5   (reads <base>, <base>_ef, <base>_signabl)"""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ap = argparse.ArgumentParser()
ap.add_argument("--base", default="a0.5_K5")
args = ap.parse_args()
R = Path("results/linear")
rng = np.random.default_rng(20260928)


def excess(cond, rung):
    d = pd.read_csv(R / cond / f"rho_{rung}.csv")
    null = json.loads((R / cond / f"h1_{rung}.json").read_text())["null_median_q025_q50_q975"][1]
    s = d[d["set"] == "draw"].set_index("perturbation")["rho"]
    return s - null, s, null


def holm(p):
    p = np.asarray(p, float); order = np.argsort(p); out = np.empty_like(p); run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - rank) * p[i])); out[i] = run
    return out


rows = []
for hyp, cond, better in [("H2", f"{args.base}_ef", "cond"), ("H3", f"{args.base}_signabl", "base")]:
    block = []
    for rung in ["R0", "R1", "R2"]:
        eb, rb, nb = excess(args.base, rung)
        ec, rc, nc = excess(cond, rung)
        j = pd.concat([eb.rename("base"), ec.rename("cond")], axis=1).dropna()
        hi, lo = (j["cond"], j["base"]) if better == "cond" else (j["base"], j["cond"])
        diff = (hi - lo).to_numpy()
        p = float(wilcoxon(hi, lo, alternative="greater").pvalue) if len(j) else np.nan
        boot = [np.median(rng.choice(diff, len(diff))) for _ in range(1000)] if len(j) else [np.nan]
        block.append({"hypothesis": hyp, "rung": rung, "condition": cond, "n_targets": len(j),
                      "predicted": "cond>base" if better == "cond" else "base>cond",
                      "median_diff_excess": float(np.median(diff)) if len(j) else np.nan,
                      "ci95_low": float(np.quantile(boot, .025)), "ci95_high": float(np.quantile(boot, .975)),
                      "null_median_base": nb, "null_median_cond": nc,
                      "median_raw_rho_base": float(rb.loc[j.index].median()),
                      "median_raw_rho_cond": float(rc.loc[j.index].median()), "p_one_sided": p})
    ph = holm([b["p_one_sided"] for b in block])
    for b, q in zip(block, ph):
        b["p_holm"] = float(q)
    rows += block
out = pd.DataFrame(rows)
out.to_csv(R / args.base / "h2_h3.csv", index=False)
print(out.to_string(index=False))
