"""Cross-dataset agreement of linear-model responses, and the H1 test (PREREG A4).

For each rung: per-perturbation Spearman rho between datasets (A4.1), then for the drawn
targets the mismatched-target permutation test (A4.2) and a cell-label-shuffle sanity null.

Usage: python src/analyze_linear.py [--cond a0.5_K5] [--perms 1000]
Reads data/derived/linear/<cond>/, writes results/linear/<cond>/."""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

MIN_ROWS = 10
RNG_SEED = 20260928


def spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d > 0 else np.nan


def pair_rho(a, b, exclude):
    """a, b: aligned 1-D arrays; exclude: boolean mask of seed rows. Returns (rho, n, reason)."""
    if np.isnan(a).all() or np.isnan(b).all():
        return np.nan, 0, "empty_seed"
    keep = ~exclude & ~np.isnan(a) & ~np.isnan(b) & ~((a == 0) & (b == 0))
    n = int(keep.sum())
    if n < MIN_ROWS:
        return np.nan, n, "lt10_rows"
    r = spearman(a[keep], b[keep])
    return r, n, ("constant_vector" if np.isnan(r) else "")


def derangement(n, rng):
    while True:
        p = rng.permutation(n)
        if not (p == np.arange(n)).any():
            return p


def r0_reindex(df, stim_side):
    """rows '<unit>|<left|right>' -> '<unit>|ipsi' / '<unit>|contra' relative to stim side."""
    u, s = df.index.str.split("|").str[0], df.index.str.split("|").str[1]
    rel = np.where(s == stim_side, "ipsi", "contra")
    out = df.copy()
    out.index = u + "|" + rel
    return out


def analyze(A, B, seed_rows, draw_keys, perms, rng):
    """A, B: DataFrames rows x perturbations (same columns); seed_rows: key -> set of row labels."""
    idx = A.index.intersection(B.index)
    A, B = A.loc[idx], B.loc[idx]
    Av, Bv = A.to_numpy(float), B.to_numpy(float)
    col = {k: j for j, k in enumerate(A.columns)}
    seedmask = {k: idx.isin(list(seed_rows[k])) for k in A.columns}

    rows = []
    for k in A.columns:
        r, n, why = pair_rho(Av[:, col[k]], Bv[:, col[k]], seedmask[k])
        rows.append({"perturbation": k, "set": "draw" if k.startswith("t:") else "case_study",
                     "rho": r, "n_rows": n, "na_reason": why})
    rho = pd.DataFrame(rows)

    keys = [k for k in draw_keys if k in col]
    obs = rho.set_index("perturbation").loc[keys, "rho"]
    observed = float(np.nanmedian(obs))
    null = np.empty(perms)
    for p in range(perms):
        pi = derangement(len(keys), rng)
        vals = [pair_rho(Av[:, col[a]], Bv[:, col[keys[j]]], seedmask[a] | seedmask[keys[j]])[0]
                for a, j in zip(keys, pi)]
        null[p] = np.nanmedian(vals)
    label_null = []
    for k in keys:
        a, b, ex = Av[:, col[k]], Bv[:, col[k]], seedmask[k]
        if np.isnan(a).all() or np.isnan(b).all():
            continue
        keep = ~ex & ~np.isnan(a) & ~np.isnan(b) & ~((a == 0) & (b == 0))
        if keep.sum() < MIN_ROWS:
            continue
        aa, bb = a[keep], b[keep]
        v = np.nanmean([spearman(aa, rng.permutation(bb)) for _ in range(perms)])
        if not np.isnan(v):
            label_null.append(v)
    h1 = {"n_targets": int(obs.notna().sum()), "n_na": int(obs.isna().sum()),
          "observed_median_rho": observed,
          "null_median_q025_q50_q975": [float(x) for x in np.quantile(null, [.025, .5, .975])],
          "p_one_sided": float((1 + (null >= observed).sum()) / (perms + 1)),
          "label_shuffle_null_median_rho": float(np.median(label_null)) if label_null else None,
          "perms": perms, "rng_seed": RNG_SEED}
    return rho, h1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", default="a0.5_K5")
    ap.add_argument("--perms", type=int, default=1000)
    args = ap.parse_args()
    D, OUT = Path(f"data/derived/linear/{args.cond}"), Path(f"results/linear/{args.cond}")
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    targets = pd.read_csv("results/prereg/targets_DRAFT.csv")
    fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
    case = pd.read_csv("results/prereg/case_study.csv")
    r2map = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]
    units = {f"t:{u}": {int(u)} for u in targets["unit_id"]}
    for _, row in case.iterrows():
        units[f"c:{row['name']}"] = set(fwu.loc[fwu["cell_type"].isin(row["flywire_types"].split("|")),
                                                "unit_id"].astype(int))
    draw_all = [f"t:{u}" for u in targets["unit_id"]]
    draw_B = [f"t:{u}" for u in targets.loc[targets["tierB"], "unit_id"]]

    rd = lambda f: pd.read_parquet(D / f"{f}.parquet")
    rungs = {
        "R1": (rd("R1_flywire"), rd("R1_male_cns"), {k: {str(u) for u in v} for k, v in units.items()}, draw_all),
        "R2": (rd("R2_flywire"), rd("R2_hemibrain"),
               {k: {str(int(r2map[u])) for u in v} for k, v in units.items()}, draw_B),
        "R2right": (rd("R2right_flywire"), rd("R2_hemibrain"),
                    {k: {str(int(r2map[u])) for u in v} for k, v in units.items()}, draw_B),
        "R0": (r0_reindex(rd("R0_flywire_left"), "left"), r0_reindex(rd("R0_flywire_right"), "right"),
               {k: {f"{u}|{s}" for u in v for s in ("ipsi", "contra")} for k, v in units.items()}, draw_all),
    }
    for name, (A, B, seeds, draw) in rungs.items():
        cols = [c for c in A.columns if c in B.columns]
        rho, h1 = analyze(A[cols], B[cols], seeds, draw, args.perms, rng)
        rho.to_csv(OUT / f"rho_{name}.csv", index=False)
        (OUT / f"h1_{name}.json").write_text(json.dumps(h1, indent=1))
        print(f"{name}: wrote rho_{name}.csv and h1_{name}.json", flush=True)


if __name__ == "__main__":
    main()
