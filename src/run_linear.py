"""Run the linear model (A2.1) for every rung (A2.2) and save response vectors ONLY.

No correlations or summary statistics are computed here, so running this does not show
any results. Outputs go to data/derived/linear/a<alpha>_K<K>/ (gitignored):
    R1_flywire.parquet, R1_male_cns.parquet        rows: unit_id
    R2_flywire.parquet, R2_hemibrain.parquet       rows: r2_unit_id
    R2right_flywire.parquet                        rows: r2_unit_id, FlyWire right side only
    R0_flywire_left.parquet, R0_flywire_right.parquet
                                                   rows: "<unit_id>|<left|right>", R0-usable units
Columns are perturbations: "t:<unit_id>" for the 134-target draw, "c:<name>" for the
case-study set (A3.1). A seed with no neurons in a dataset gives an all-NaN column and is
listed in manifest.json. The seed unit is NOT removed here; the analysis step removes it (A2.1).

Usage: python src/run_linear.py [--alpha 0.5] [--K 5]"""
import argparse, json, subprocess, time
from pathlib import Path

import numpy as np
import pandas as pd

from linear_model import load_graph, propagate, unit_means

p = argparse.ArgumentParser()
p.add_argument("--alpha", type=float, default=0.5)
p.add_argument("--K", type=int, default=5)
p.add_argument("--batch", type=int, default=50)
args = p.parse_args()
OUT = Path(f"data/derived/linear/a{args.alpha}_K{args.K}")
OUT.mkdir(parents=True, exist_ok=True)

# ---------- perturbations ----------
targets = pd.read_csv("results/prereg/targets_DRAFT.csv")
fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
case = pd.read_csv("results/prereg/case_study.csv")
r2map = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]
r0units = set(pd.read_csv("data/derived/r0_units.csv").query("r0_usable")["unit_id"])

perts = {f"t:{u}": {int(u)} for u in targets["unit_id"]}
for _, row in case.iterrows():
    types = row["flywire_types"].split("|")
    us = set(fwu.loc[fwu["cell_type"].isin(types), "unit_id"].astype(int))
    assert us, f"case-study {row['name']}: no units for {types}"
    perts[f"c:{row['name']}"] = us
tierB = {f"t:{u}" for u in targets.loc[targets["tierB"], "unit_id"]}
case_r2 = {"c:" + n for n in case["name"]} - {"c:sugar_GRN", "c:bitter_GRN", "c:JO-A", "c:JO-CE"}
r2_perts = [k for k in perts if k in tierB or k in case_r2]
r2sets = {k: {int(r2map[x]) for x in perts[k]} for k in r2_perts}

side = pd.read_csv("data/flywire_783/annotations_v3_neurons.tsv", sep="\t",
                   usecols=["root_id", "side"], low_memory=False).set_index("root_id")["side"]
manifest = {"alpha": args.alpha, "K": args.K, "started": time.strftime("%Y-%m-%d %H:%M:%S"),
            "git": subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                                  text=True).stdout.strip(),
            "seed_neurons": {}, "empty_seeds": {}}


def run(g, name, keys, seed_sets, labels):
    t0 = time.time()
    seeds = [np.flatnonzero(m) for m in seed_sets]
    manifest["seed_neurons"][name] = {k: int(len(s)) for k, s in zip(keys, seeds)}
    manifest["empty_seeds"][name] = [k for k, s in zip(keys, seeds) if len(s) == 0]
    cols = []
    for i in range(0, len(keys), args.batch):
        R = propagate(g.A, seeds[i:i + args.batch], args.alpha, args.K)
        cols.append(unit_means(R, labels))
    df = pd.concat(cols, axis=1)
    df.columns = keys
    empty = manifest["empty_seeds"][name]
    df[empty] = np.nan
    df.index = df.index.astype(str)
    df.to_parquet(OUT / f"{name}.parquet")
    print(f"{name}: {len(keys)} perturbations x {len(df)} readout rows, "
          f"{len(empty)} empty seeds, {time.time() - t0:.0f}s", flush=True)


keys = list(perts)
num = lambda ser: ser.to_numpy(dtype="float64", na_value=np.nan)
isin = lambda arr, s: np.isin(arr, list(s))

# ---------- FlyWire: R1, R2, R2right, R0 ----------
g = load_graph("flywire")
u, r2 = num(g.unit), num(g.r2)
sd = pd.Series(g.ids).map(side).to_numpy()
run(g, "R1_flywire", keys, [isin(u, perts[k]) for k in keys], g.unit)

run(g, "R2_flywire", r2_perts, [isin(r2, r2sets[k]) for k in r2_perts], g.r2)
right = sd == "right"
run(g, "R2right_flywire", r2_perts, [isin(r2, r2sets[k]) & right for k in r2_perts],
    g.r2.where(pd.Series(right), pd.NA))

lab0 = pd.Series([f"{a}|{b}" if (pd.notna(a) and int(a) in r0units and b in ("left", "right")) else None
                  for a, b in zip(u, sd)], dtype=object)
for s in ("left", "right"):
    run(g, f"R0_flywire_{s}", keys, [isin(u, perts[k]) & (sd == s) for k in keys], lab0)
del g

# ---------- male CNS: R1 ----------
g = load_graph("male_cns")
run(g, "R1_male_cns", keys, [isin(num(g.unit), perts[k]) for k in keys], g.unit)
del g

# ---------- hemibrain: R2 ----------
g = load_graph("hemibrain")
run(g, "R2_hemibrain", r2_perts, [isin(num(g.r2), r2sets[k]) for k in r2_perts], g.r2)
del g

manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=1))
print("done ->", OUT)
