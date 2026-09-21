"""Draw perturbation targets in two tiers, fixed seeds.
Tier A: male CNS vs FlyWire, up to 20 primary units per class.
Tier B: all three datasets (hb_share > 0.5), up to 20 primary units per class."""
import pandas as pd, numpy as np, os

SEED = 20260921
d = pd.read_csv("data/derived/perturbation_candidates.csv")
cls = {"cb_sensory": "sensory", "ol_sensory": "sensory", "descending_neuron": "descending",
       "ascending_neuron": "ascending", "visual_projection": "visual_projection"}
d = d[d["primary"] & d["superclass"].isin(cls)].copy()
d["stratum"] = d["superclass"].map(cls)

def draw(pool, seed):
    rng = np.random.default_rng(seed)
    parts = []
    for st, g in pool.sort_values("unit_id").groupby("stratum"):
        k = min(20, len(g))
        idx = sorted(rng.choice(len(g), size=k, replace=False))
        parts.append(g.iloc[idx])
    return pd.concat(parts)

A = draw(d, SEED)
B = draw(d[d["hb_share"] > 0.5], SEED + 1)
ids = sorted(set(A["unit_id"]) | set(B["unit_id"]))
t = d[d["unit_id"].isin(ids)].copy()
t["tierA"] = t["unit_id"].isin(A["unit_id"])
t["tierB"] = t["unit_id"].isin(B["unit_id"])
cols = ["unit_id", "name", "stratum", "superclass", "n_male", "n_flywire", "ratio",
        "hb_share", "tierA", "tierB", "strict", "var"]
os.makedirs("results/prereg", exist_ok=True)
t.sort_values("unit_id")[cols].to_csv("results/prereg/targets_DRAFT.csv", index=False)

print("union targets", len(t), "| tier A", int(t["tierA"].sum()), "| tier B", int(t["tierB"].sum()),
      "| in both", int((t["tierA"] & t["tierB"]).sum()))
print(t.groupby("stratum").agg(n=("unit_id", "size"), A=("tierA", "sum"), B=("tierB", "sum"),
                               strict=("strict", "sum")))
runs = 2 * len(t) + int(t["tierB"].sum())
print("simulation runs (target x dataset):", runs)
for trials in (30, 10):
    print(f"LIF at {trials} trials x 1000 ms, 5 min per trial per thread: {runs*trials*5/60:.0f} CPU-hours")
