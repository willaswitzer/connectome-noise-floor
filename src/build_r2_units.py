"""R2 units (PREREG A2.8): connected components of the bipartite graph linking matched
units to the hemibrain bodyIds they map to (unit_hemibrain_map.csv). Units sharing any
hemibrain body fall in one R2 unit. r2_unit_id = smallest unit_id in the component, so ids
are deterministic. Units with no hemibrain body are singletons (r2_unit_id = unit_id).

Output: data/derived/r2_units.csv  (unit_id, r2_unit_id, r2_size)
        data/derived/r2_body_map.csv  (bodyId, r2_unit_id), one row per hemibrain body"""
import pandas as pd

m = pd.read_csv("data/derived/unit_hemibrain_map.csv")[["unit_id", "bodyId"]].drop_duplicates()
all_units = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")["unit_id"].unique()

parent = {u: u for u in all_units}
parent.update({u: u for u in m["unit_id"].unique()})


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


for _, units in m.groupby("bodyId")["unit_id"]:
    units = list(units)
    for u in units[1:]:
        a, b = find(units[0]), find(u)
        if a != b:
            parent[max(a, b)] = min(a, b)

r2 = pd.DataFrame({"unit_id": list(parent)})
r2["r2_unit_id"] = r2["unit_id"].map(find)
r2["r2_size"] = r2.groupby("r2_unit_id")["unit_id"].transform("size")
r2.sort_values("unit_id").to_csv("data/derived/r2_units.csv", index=False)

bm = m.assign(r2_unit_id=m["unit_id"].map(r2.set_index("unit_id")["r2_unit_id"]))
assert bm.groupby("bodyId")["r2_unit_id"].nunique().max() == 1
bm[["bodyId", "r2_unit_id"]].drop_duplicates().to_csv("data/derived/r2_body_map.csv", index=False)

t = pd.read_csv("results/prereg/targets_DRAFT.csv")
tb = t[t["tierB"]].merge(r2, on="unit_id", how="left")
shared = tb.groupby("r2_unit_id")["unit_id"].transform("size") > 1
print(f"units: {len(r2):,}; R2 units: {r2['r2_unit_id'].nunique():,}; "
      f"R2 units merging >1 unit: {int((r2.drop_duplicates('r2_unit_id')['r2_size'] > 1).sum()):,}")
print(f"Tier B targets in a merged R2 unit: {int((tb['r2_size'] > 1).sum())} / {len(tb)}; "
      f"sharing an R2 unit with another target: {int(shared.sum())}")
print(tb[tb["r2_size"] > 1][["name", "unit_id", "r2_unit_id", "r2_size"]].to_string(index=False))
