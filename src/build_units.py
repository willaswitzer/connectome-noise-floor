"""Build matched cell-type units between male CNS, FlyWire and hemibrain.
Rules (pre-registered): keep neurons with a superclass; drop every neuron with any
dimorphism label; merge many-to-many label pairs into composite units."""
import pandas as pd, networkx as nx

A = "data/male_cns/body-annotations-male-cns-v1.0-minconf-0.5.feather"
OUT = "data/derived/matched_units_malecns.csv"

a = pd.read_feather(A)
n = a[a["superclass"].notna() & a["dimorphism"].isna()]
b = n[n["hemibrainType"].notna() & n["flywireType"].notna()].copy()

G = nx.Graph()
for h, f in b[["hemibrainType", "flywireType"]].drop_duplicates().itertuples(index=False):
    G.add_edge(("h", h), ("f", f))
comps = list(nx.connected_components(G))
cid = {node: i for i, c in enumerate(comps) for node in c}

b["unit_id"] = b["flywireType"].map(lambda f: cid[("f", f)])
b["unit_size"] = b["unit_id"].map(lambda i: len(comps[i]))
b["is_1to1"] = b["unit_size"] == 2
b[["bodyId", "superclass", "hemibrainType", "flywireType",
   "unit_id", "unit_size", "is_1to1"]].to_csv(OUT, index=False)

print("units", len(comps), "| 1:1", int(sum(len(c) == 2 for c in comps)),
      "| neurons", len(b), "| neurons in 1:1 units", int(b["is_1to1"].sum()))
