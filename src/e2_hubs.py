"""EXPLORATORY E2 (not pre-registered; planned 2026-09-29 after H2/H3 were seen).

Question: under sign ablation (A6.2), do responses mostly follow structural hub connectivity
or a seed-independent "generic" profile, and does that explain the H3 reversal (higher
cross-dataset agreement when every neuron is excitatory)?

Rung R1 (FlyWire vs male CNS, whole-unit readout, all 134 drawn targets). Rows per target are
the A4.1 rows for that target (seed unit removed, rows missing in either dataset removed, rows
zero in both removed). For each dataset and condition (base, sign-ablated; rectified and
unrectified), per drawn target:
  rho_indeg     Spearman(response, unit-mean number of presynaptic partners)
  rho_insyn     Spearman(response, unit-mean input synapse count)
  rho_eig       Spearman(response, unit-mean eigenvector centrality of the synapse-count matrix
                (post x pre, leading Perron vector by power iteration). The input-fraction
                operator itself is row-stochastic when every neuron is excitatory, so its
                leading eigenvector is flat and uninformative; the count matrix is used instead.)
  rho_generic   Spearman(response, leave-one-out mean of the other drawn targets' responses,
                each column scaled to unit L1 norm) in the same dataset
Cross-dataset:
  rho_cross     the A4.1 rho (FlyWire vs male CNS), recomputed here as a check
  rho_partial   partial Spearman: ranks of each dataset's response regressed on the ranks of
                that dataset's generic profile; Pearson correlation of the two residuals
Also: Spearman between the two datasets' full generic profiles (all drawn targets).
Descriptive only; no tests.

Usage: python src/e2_hubs.py -> results/linear/exploratory/e2_per_target.csv, e2_summary.csv"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from analyze_linear import spearman, MIN_ROWS
from linear_model import GRAPHS, load_graph, unit_means

D = Path("data/derived/linear")
OUT = Path("results/linear/exploratory")
OUT.mkdir(parents=True, exist_ok=True)
DATASETS = {"flywire": "R1_flywire", "male_cns": "R1_male_cns"}
CONDS = {"rectified": ("a0.5_K5", "a0.5_K5_signabl"), "unrectified": ("a0.5_K5_norect", "a0.5_K5_norect_signabl")}


def hub_scores(name, iters=300):
    g = load_graph(name, sign_ablate=True)
    e = pd.read_parquet(f"{GRAPHS}/{name}.parquet", columns=["bodyId_pre", "bodyId_post", "weight"])
    pre = np.searchsorted(g.ids, e["bodyId_pre"].to_numpy())
    post = np.searchsorted(g.ids, e["bodyId_post"].to_numpy())
    n = len(g.ids)
    import scipy.sparse as sp
    N = sp.csr_matrix((e["weight"].to_numpy(np.float64), (post, pre)), shape=(n, n))
    indeg = np.diff(N.indptr).astype(float)
    insyn = np.asarray(N.sum(axis=1)).ravel()
    x = np.ones(n) / np.sqrt(n)
    for i in range(iters):
        y = N @ x
        y /= np.linalg.norm(y)
        delta = np.abs(y - x).max()
        x = y
    print(f"{name}: eigenvector power iteration, last max change {delta:.2e}", flush=True)
    H = unit_means(np.column_stack([indeg, insyn, x]), g.unit)
    H.columns = ["indeg", "insyn", "eig"]
    H.index = H.index.astype(str)
    return H


def generic_profiles(df, keys):
    X = df[keys].to_numpy(float)
    X = X / np.nansum(np.abs(X), axis=0, keepdims=True)
    ok = ~np.isnan(X).all(axis=0) & (np.nansum(np.abs(X), axis=0) > 0)
    tot = np.nansum(X[:, ok], axis=1)
    cnt = ok.sum()
    gen = {}
    for j, k in enumerate(keys):
        gen[k] = (tot - (np.nan_to_num(X[:, j]) if ok[j] else 0)) / (cnt - (1 if ok[j] else 0))
    return pd.DataFrame(gen, index=df.index), pd.Series(tot / cnt, index=df.index)


def partial_rho(a, b, ga, gb):
    ra, rb, rga, rgb = (rankdata(v) for v in (a, b, ga, gb))
    res = []
    for r, g in ((ra, rga), (rb, rgb)):
        G = np.column_stack([np.ones_like(g), g])
        beta = np.linalg.lstsq(G, r, rcond=None)[0]
        res.append(r - G @ beta)
    d = np.sqrt((res[0] ** 2).sum() * (res[1] ** 2).sum())
    return float((res[0] * res[1]).sum() / d) if d > 0 else np.nan


targets = pd.read_csv("results/prereg/targets_DRAFT.csv")
keys = [f"t:{u}" for u in targets["unit_id"]]
hubs = {d: hub_scores(d) for d in DATASETS}

rows, gen_agree = [], []
for model, (cbase, cabl) in CONDS.items():
    for cond_label, cond in (("base", cbase), ("sign_ablated", cabl)):
        A = pd.read_parquet(D / cond / "R1_flywire.parquet")
        B = pd.read_parquet(D / cond / "R1_male_cns.parquet")
        idx = A.index.intersection(B.index)
        A, B = A.loc[idx], B.loc[idx]
        gA, gA_all = generic_profiles(A, keys)
        gB, gB_all = generic_profiles(B, keys)
        both = gA_all.notna() & gB_all.notna()
        gen_agree.append({"model": model, "condition": cond_label,
                          "rho_generic_flywire_vs_male_cns": spearman(gA_all[both].to_numpy(), gB_all[both].to_numpy())})
        HA, HB = hubs["flywire"].reindex(idx), hubs["male_cns"].reindex(idx)
        for k in keys:
            a, b = A[k].to_numpy(float), B[k].to_numpy(float)
            if np.isnan(a).all() or np.isnan(b).all():
                continue
            keep = ~idx.isin([k[2:]]) & ~np.isnan(a) & ~np.isnan(b) & ~((a == 0) & (b == 0))
            keep &= HA.notna().all(axis=1).to_numpy() & HB.notna().all(axis=1).to_numpy()
            if keep.sum() < MIN_ROWS:
                continue
            aa, bb = a[keep], b[keep]
            rec = {"model": model, "condition": cond_label, "perturbation": k, "n_rows": int(keep.sum()),
                   "rho_cross": spearman(aa, bb),
                   "rho_partial": partial_rho(aa, bb, gA[k].to_numpy()[keep], gB[k].to_numpy()[keep])}
            for dname, v, H, gen in (("flywire", aa, HA, gA), ("male_cns", bb, HB, gB)):
                for h in ("indeg", "insyn", "eig"):
                    rec[f"rho_{h}_{dname}"] = spearman(v, H[h].to_numpy()[keep])
                rec[f"rho_generic_{dname}"] = spearman(v, gen[k].to_numpy()[keep])
            rows.append(rec)

per = pd.DataFrame(rows)
per.to_csv(OUT / "e2_per_target.csv", index=False)
metrics = [c for c in per.columns if c.startswith("rho_")]
summ = per.groupby(["model", "condition"])[metrics].median().T
summ.columns = [f"{m}|{c}" for m, c in summ.columns]
# paired (ablated - base) on targets non-NA in both
for model in CONDS:
    p = per[per["model"] == model].pivot(index="perturbation", columns="condition", values=metrics)
    for m in metrics:
        j = p[m].dropna()
        summ.loc[m, f"{model}|median_paired_diff_ablated_minus_base"] = float((j["sign_ablated"] - j["base"]).median())
        summ.loc[m, f"{model}|n_paired"] = len(j)
summ.to_csv(OUT / "e2_summary.csv")
pd.DataFrame(gen_agree).to_csv(OUT / "e2_generic_agreement.csv", index=False)
pd.set_option("display.width", 250)
print(summ.round(3).to_string())
print(pd.DataFrame(gen_agree).to_string(index=False))
