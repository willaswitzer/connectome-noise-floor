"""Descriptive LIF readouts for the case-study set (Amendment 7). Not used for H1-H4.

Unit response = mean firing rate (Hz) over the unit's graph nodes and trials; silent nodes count
as 0. Units are those of the linear model: R1 units, R2 units, R0 "<unit>|ipsi/contra" (A2.2, A4.1).
Per perturbation and rung:
  rho, n_rows, na_reason   A4.1 (seed unit(s) removed, rows missing in either dataset removed,
                           double-zero rows removed, NA below 10 rows), all 30 trials
  rho_10                   the same from the first 10 trials only (A7.4)
  rho_linear               the linear-model rho for the same perturbation and rung (a0.5_K5)
  top{10,20,50}_overlap    |top-k units in A  intersect  top-k units in B| / k, over the A4.1 rows
  named readouts           rate in each dataset for the A3.1 named units, both > 0 or both = 0
Per run and dataset readout: rho between the 10-trial and 30-trial unit vectors (A7.4), active
neurons, spikes per trial.

Usage: python src/analyze_lif.py [--wsyn primary]   -> results/lif/<wsyn>/"""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_linear import pair_rho, spearman

ap = argparse.ArgumentParser()
ap.add_argument("--wsyn", default="primary", choices=["primary", "fixed"])
args = ap.parse_args()
L = Path("data/derived/lif")
OUT = Path(f"results/lif/{args.wsyn}")
OUT.mkdir(parents=True, exist_ok=True)
case = pd.read_csv("results/prereg/case_study.csv")
fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
r2map = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]
units_of = lambda types: set(fwu.loc[fwu["cell_type"].isin(types), "unit_id"].astype(int))
seed_units = {r["name"]: units_of(r["flywire_types"].split("|")) for _, r in case.iterrows()}
named = {r["name"]: (r["named_readouts"].split("|") if isinstance(r["named_readouts"], str) else [])
         for _, r in case.iterrows()}
RUNGS = {"R1": ("flywire", "male_cns"), "R2": ("flywire", "hemibrain"), "R0": ("flywire", "flywire")}

meta = {d: json.loads((L / d / "seeds.json").read_text()) for d in ("flywire", "male_cns", "hemibrain")}
labels = {d: pd.read_parquet(L / d / "labels.parquet") for d in meta}


def run_dir(d):
    return L / d / ("primary" if d == "flywire" else args.wsyn)


def node_rates(d, run, trials=None):
    r = pd.read_parquet(run_dir(d) / f"{run}.parquet")
    n = 30
    if trials is not None:
        r, n = r[r["trial"] < trials], trials
    return r.groupby("node_id")["n_spikes"].sum() / n


def unit_vector(d, run, scheme, trials=None, stim_side=None):
    lab = labels[d]
    rate = lab["node_id"].map(node_rates(d, run, trials)).fillna(0.0)
    key = lab[scheme]
    ok = key.notna()
    k = key[ok] if scheme == "r0" else key[ok].astype("int64").astype(str)  # stored as float
    v = rate[ok].groupby(k).mean()
    if scheme == "r0":
        u, s = v.index.str.split("|").str[0], v.index.str.split("|").str[1]
        v.index = u + "|" + np.where(s == stim_side, "ipsi", "contra")
    return v


def lookup(d, p, rung, stim):
    for x in meta[d]["index"]:
        if x["perturbation"] == p and x["rung"] == rung and x["stim"] == stim:
            return x["run"]
    return None


def linear_rho(rung):
    f = Path(f"results/linear/a0.5_K5/rho_{rung}.csv")
    return pd.read_csv(f).set_index("perturbation")["rho"] if f.exists() else pd.Series(dtype=float)


rows, named_rows, trial_rows = [], [], []
for rung, (da, db) in RUNGS.items():
    lin = linear_rho(rung)
    for p in seed_units:
        if rung == "R0":
            ka, kb, sa, sb, scheme = lookup(da, p, "R0", "left"), lookup(db, p, "R0", "right"), "left", "right", "r0"
            seeds = {f"{u}|{s}" for u in seed_units[p] for s in ("ipsi", "contra")}
        else:
            ka, kb, sa, sb = lookup(da, p, rung, "both"), lookup(db, p, rung, "both"), None, None
            scheme = "unit" if rung == "R1" else "r2"
            seeds = {str(u) for u in seed_units[p]} if rung == "R1" else {str(int(r2map[u])) for u in seed_units[p]}
        if (rung == "R2" and ka is None and kb is None):
            continue  # not on this rung (A3.1)
        rec = {"perturbation": p, "rung": rung, "run_a": ka, "run_b": kb,
               "rho_linear": float(lin.get(f"c:{p}", np.nan))}
        if ka is None or kb is None or not (run_dir(da) / f"{ka}.parquet").exists() \
                or not (run_dir(db) / f"{kb}.parquet").exists():
            rec.update({"rho": np.nan, "na_reason": "empty_seed" if (ka is None or kb is None) else "not_run"})
            rows.append(rec)
            continue
        vec = {}
        for T in (30, 10):
            A = unit_vector(da, ka, scheme, None if T == 30 else 10, sa)
            B = unit_vector(db, kb, scheme, None if T == 30 else 10, sb)
            idx = A.index.intersection(B.index)
            a, b = A.loc[idx].to_numpy(), B.loc[idx].to_numpy()
            ex = idx.isin(list(seeds))
            r, n, why = pair_rho(a, b, ex)
            vec[T] = (A, B)
            if T == 30:
                rec.update({"rho": r, "n_rows": n, "na_reason": why})
                keep = ~ex & ~((a == 0) & (b == 0))
                for k in (10, 20, 50):
                    ta = set(idx[keep][np.argsort(-a[keep])[:k]])
                    tb = set(idx[keep][np.argsort(-b[keep])[:k]])
                    rec[f"top{k}_overlap"] = len(ta & tb) / k
                for nm in named[p]:
                    for u in units_of([nm]):
                        key = str(u) if scheme == "unit" else (str(int(r2map[u])) if scheme == "r2" else f"{u}|ipsi")
                        ra, rb = float(A.get(key, np.nan)), float(B.get(key, np.nan))
                        named_rows.append({"perturbation": p, "rung": rung, "named": nm, "unit": key,
                                           "rate_a": ra, "rate_b": rb,
                                           "both_active_or_both_silent": (ra > 0) == (rb > 0),
                                           "ratio_b_over_a": rb / ra if ra > 0 else np.nan})
            else:
                rec["rho_10"] = r
        rows.append(rec)
        for side, d, k, T30, T10 in (("a", da, ka, vec[30][0], vec[10][0]), ("b", db, kb, vec[30][1], vec[10][1])):
            j = T30.index.intersection(T10.index)
            m = (T30.loc[j] != 0) | (T10.loc[j] != 0)
            info = json.loads((run_dir(d) / f"{k}.json").read_text())
            trial_rows.append({"perturbation": p, "rung": rung, "dataset": d, "run": k, "side": side,
                               "rho_10_vs_30": spearman(T30.loc[j][m].to_numpy(), T10.loc[j][m].to_numpy()),
                               "n_active_neurons": info["n_active_neurons"],
                               "frac_active": info["n_active_neurons"] / len(labels[d]),
                               "elapsed_s": info["elapsed_s"]})

rho = pd.DataFrame(rows)
rho["abs_diff_10_vs_30"] = (rho["rho"] - rho["rho_10"]).abs()
rho.to_csv(OUT / "rho_case_study.csv", index=False)
pd.DataFrame(named_rows).to_csv(OUT / "named_readouts.csv", index=False)
tr = pd.DataFrame(trial_rows).drop_duplicates(["dataset", "run", "rung"])
tr.to_csv(OUT / "trial_check_and_activity.csv", index=False)
summ = {"n_rho": int(rho["rho"].notna().sum()),
        "median_rho_by_rung": rho.groupby("rung")["rho"].median().round(4).to_dict(),
        "median_rho_linear_by_rung": rho.groupby("rung")["rho_linear"].median().round(4).to_dict(),
        "trial_check_median_abs_diff_rho": float(rho["abs_diff_10_vs_30"].median()),
        "trial_check_max_abs_diff_rho": float(rho["abs_diff_10_vs_30"].max()),
        "trial_check_median_rho_10_vs_30_unit_vectors": float(tr["rho_10_vs_30"].median())}
(OUT / "summary.json").write_text(json.dumps(summ, indent=1))
print(rho.round(3).to_string(index=False))
print(json.dumps(summ, indent=1))
