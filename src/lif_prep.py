"""Prepare LIF inputs for Amendment 7 (case-study set only).

For each dataset (flywire, male_cns, hemibrain) writes to data/derived/lif/<dataset>/:
    completeness.csv   node ids of the scoped graph (A2.7), in Shiu's `Completeness` format
    connectivity.parquet   Shiu's connectivity format: Presynaptic/Postsynaptic_ID and _Index,
                       Connectivity = synapse count, Excitatory = A3.4 sign (0 = no sign),
                       Excitatory x Connectivity = sign * count
    labels.parquet     per node: node_id, unit, r2, side (FlyWire: annotations v3 `side`;
                       male CNS: somaSide/rootSide), r0 label ("<unit>|<side>" for R0-usable units)
    seeds.json         run key -> {"perturbation", "rung", "stim", "nodes": [...]}
Seeds follow A2.2 exactly as run_linear.py builds them:
    R1  FlyWire and male CNS: all seed-unit neurons, both sides
    R0  FlyWire: left-side seed neurons (one run), right-side seed neurons (another run)
    R2  FlyWire and hemibrain: all neurons of the seed's R2 unit(s); case-study seeds without
        hemibrain neurons (sugar_GRN, bitter_GRN, JO-A, JO-CE) are excluded (A3.1)
Runs with identical node sets within a dataset are deduplicated (same simulation), and
seeds.json maps every (perturbation, rung, stim) to the run key that serves it.

Usage: python src/lif_prep.py [--datasets flywire,male_cns,hemibrain]"""
import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd

from linear_model import GRAPHS, SIGNS

ap = argparse.ArgumentParser()
ap.add_argument("--datasets", default="flywire,male_cns,hemibrain")
args = ap.parse_args()

case = pd.read_csv("results/prereg/case_study.csv")
fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
r2map = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]
r0units = set(pd.read_csv("data/derived/r0_units.csv").query("r0_usable")["unit_id"])
NO_HB = {"sugar_GRN", "bitter_GRN", "JO-A", "JO-CE"}
perts = {}
for _, row in case.iterrows():
    us = set(fwu.loc[fwu["cell_type"].isin(row["flywire_types"].split("|")), "unit_id"].astype(int))
    assert us, row["name"]
    perts[row["name"]] = us


def node_table(name):
    e = pd.read_parquet(f"{GRAPHS}/{name}.parquet")
    ids, inv = np.unique(np.concatenate([e["bodyId_pre"].to_numpy(), e["bodyId_post"].to_numpy()]),
                         return_inverse=True)
    pre, post = inv[:len(e)], inv[len(e):]
    lab = pd.DataFrame({"node_id": ids})
    for col, (cp, cq) in {"unit": ("unit_id_pre", "unit_id_post"), "r2": ("r2_unit_pre", "r2_unit_post")}.items():
        m = pd.concat([pd.Series(e[cp].to_numpy(), index=pre), pd.Series(e[cq].to_numpy(), index=post)])
        m = m[~m.index.duplicated()].sort_index()
        lab[col] = m.reindex(range(len(ids))).astype("Int64").to_numpy()
    f, key = SIGNS[name]
    s = pd.read_parquet(f"data/derived/{f}", columns=[key, "sign"]).set_index(key)["sign"]
    sign = pd.Series(ids).map(s).fillna(0).astype(int).to_numpy()
    con = pd.DataFrame({"Presynaptic_ID": ids[pre], "Postsynaptic_ID": ids[post],
                        "Presynaptic_Index": pre, "Postsynaptic_Index": post,
                        "Connectivity": e["weight"].to_numpy(np.int64)})
    con["Excitatory"] = sign[pre]
    con["Excitatory x Connectivity"] = con["Excitatory"] * con["Connectivity"]
    return lab, con


for name in args.datasets.split(","):
    out = Path(f"data/derived/lif/{name}")
    out.mkdir(parents=True, exist_ok=True)
    lab, con = node_table(name)
    if name == "flywire":
        side = pd.read_csv("data/flywire_783/annotations_v3_neurons.tsv", sep="\t",
                           usecols=["root_id", "side"], low_memory=False).set_index("root_id")["side"]
        lab["side"] = lab["node_id"].map(side)
    elif name == "male_cns":
        ma = pd.read_feather("data/male_cns/body-annotations-male-cns-v1.0-minconf-0.5.feather",
                             columns=["bodyId", "somaSide", "rootSide"]).set_index("bodyId")
        ms = ma["somaSide"].where(ma["somaSide"].isin(["L", "R"]), ma["rootSide"].where(ma["rootSide"].isin(["L", "R"])))
        lab["side"] = lab["node_id"].map(ms).map({"L": "left", "R": "right"})
    else:
        lab["side"] = None
    lab["r0"] = [f"{int(u)}|{s}" if (pd.notna(u) and int(u) in r0units and s in ("left", "right")) else None
                 for u, s in zip(lab["unit"], lab["side"])]
    # Shiu's model indexes neurons by row order of the completeness table; keep node order = index order
    pd.DataFrame({"Completed": True}, index=pd.Index(lab["node_id"], name=None)).to_csv(out / "completeness.csv")
    con.to_parquet(out / "connectivity.parquet")
    lab.to_parquet(out / "labels.parquet")

    u, r2, sd = lab["unit"], lab["r2"], lab["side"]
    wanted = []
    for p, us in perts.items():
        if name in ("flywire", "male_cns"):
            wanted.append((p, "R1", "both", u.isin(list(us))))
        if name == "flywire":
            for s in ("left", "right"):
                wanted.append((p, "R0", s, u.isin(list(us)) & (sd == s)))
        if name in ("flywire", "hemibrain") and p not in NO_HB:
            r2s = {int(r2map[x]) for x in us}
            wanted.append((p, "R2", "both", r2.isin(list(r2s))))
    runs, index = {}, []
    for p, rung, stim, mask in wanted:
        nodes = sorted(int(x) for x in lab.loc[mask.fillna(False).to_numpy(bool), "node_id"])
        key = hashlib.sha1(",".join(map(str, nodes)).encode()).hexdigest()[:12] if nodes else None
        if key:
            runs.setdefault(key, nodes)
        index.append({"perturbation": p, "rung": rung, "stim": stim, "run": key, "n_seed_neurons": len(nodes)})
    (out / "seeds.json").write_text(json.dumps({"runs": runs, "index": index}, indent=0))
    empty = [x for x in index if x["run"] is None]
    print(f"{name}: {len(lab)} nodes, {len(con)} edges, {len(runs)} distinct runs for {len(index)} "
          f"(perturbation, rung, stim); empty seeds: {[(x['perturbation'], x['rung'], x['stim']) for x in empty]}",
          flush=True)
