"""Matched cell-type units: male CNS <-> FlyWire v783 (annotations v3.0.0).
Rules: male neurons need a superclass; any type carrying a dimorphism label in either
dataset is excluded; male names resolve to FlyWire cell types by exact match, then
hemibrain_type link, then suffix parent, then composite split; unresolved names are
dropped and counted; connected names and types merge into one unit."""
import re
import pandas as pd, networkx as nx

M = "data/male_cns/body-annotations-male-cns-v1.0-minconf-0.5.feather"
F = "data/flywire_783/annotations_v3_neurons.tsv"

a = pd.read_feather(M)
a = a[a["superclass"].notna()]
male_flag_names = set(a.loc[a["dimorphism"].notna(), "flywireType"].dropna())
m = a[a["dimorphism"].isna() & a["flywireType"].notna()].copy()
m = m[~m["flywireType"].isin(male_flag_names)]

f = pd.read_csv(F, sep="\t", low_memory=False)
f = f[f["cell_type"].notna()].copy()
fw_flag = set(f.loc[f["dimorphism"].notna() & (f["dimorphism"] != "isomorphic"), "cell_type"])
ct = set(f["cell_type"])
hb2ct = (f.dropna(subset=["hemibrain_type"])
           .groupby("hemibrain_type")["cell_type"].agg(set).to_dict())

def strip_suffix(x):
    g = re.match(r"^(.*?)_?[a-z]$", x)
    return g.group(1) if g else x

def resolve(x):
    if x in ct: return "exact", {x}
    if x in hb2ct: return "via_hemibrain_type", set(hb2ct[x])
    s = strip_suffix(x)
    if s != x and s in ct: return "suffix_parent", {s}
    core = re.sub(r"^\((.*)\)[a-z]?$", r"\1", x)
    if "," in core:
        parts = [p.strip() for p in core.split(",")]
        ok = [p if p in ct else (strip_suffix(p) if strip_suffix(p) in ct else None) for p in parts]
        if all(ok): return "composite", set(ok)
    return "unresolved", set()

names = sorted(set(m["flywireType"]))
res = {n: resolve(n) for n in names}
cnt_m = m["flywireType"].value_counts()

G = nx.Graph()
drop_unres = drop_flag = 0
for n, (rule, cts) in res.items():
    if rule == "unresolved":
        drop_unres += int(cnt_m[n]); continue
    if cts & fw_flag:
        drop_flag += int(cnt_m[n]); continue
    for c in sorted(cts):
        G.add_edge(("m", n), ("f", c))
comps = sorted((sorted(c) for c in nx.connected_components(G)), key=lambda c: c[0])
cid = {node: i for i, comp in enumerate(comps) for node in comp}

m["rule"] = m["flywireType"].map(lambda n: res[n][0])
m["unit_id"] = m["flywireType"].map(lambda n: cid.get(("m", n)))
mk = m[m["unit_id"].notna()].copy(); mk["unit_id"] = mk["unit_id"].astype(int)
f["unit_id"] = f["cell_type"].map(lambda c: cid.get(("f", c)))
fk = f[f["unit_id"].notna()].copy(); fk["unit_id"] = fk["unit_id"].astype(int)

nn = pd.Series({i: sum(1 for x in c if x[0] == "m") for i, c in enumerate(comps)}, name="n_male_names")
nt = pd.Series({i: sum(1 for x in c if x[0] == "f") for i, c in enumerate(comps)}, name="n_fw_types")
nm = mk.groupby("unit_id").size().rename("n_male")
nf = fk.groupby("unit_id").size().rename("n_flywire")
hb = fk.groupby("unit_id")["hemibrain_type"].apply(lambda x: x.notna().mean()).rename("hb_share")
s = pd.concat([nn, nt, nm, nf, hb], axis=1)
s[["n_male", "n_flywire"]] = s[["n_male", "n_flywire"]].fillna(0)
s["is_1to1"] = (s["n_male_names"] == 1) & (s["n_fw_types"] == 1)
s["ratio"] = s["n_male"] / s["n_flywire"].replace(0, float("nan"))

mk[["bodyId", "superclass", "flywireType", "rule", "unit_id"]].to_csv("data/derived/units_v2_male_neurons.csv", index=False)
fk[["root_id", "cell_type", "super_class", "hemibrain_type", "unit_id"]].to_csv("data/derived/units_v2_flywire_neurons.csv", index=False)
s.to_csv("data/derived/units_v2_summary.csv", index_label="unit_id")

print("male neurons in scope", len(m), "| kept", len(mk),
      "| dropped unresolved", drop_unres, "| dropped flagged", drop_flag)
print("units", len(s), "| 1:1", int(s["is_1to1"].sum()),
      "| largest unit (names + types)", int((s["n_male_names"] + s["n_fw_types"]).max()))
print("flywire neurons in units", len(fk))
print(mk["rule"].value_counts())
r = s["ratio"].dropna()
print("count ratio male/flywire per unit: median", round(r.median(), 2),
      "| units within 0.8-1.25:", round(r.between(0.8, 1.25).mean(), 3))
ok = s["ratio"].between(0.8, 1.25)
print("neuron-weighted share within 0.8-1.25:", round(s.loc[ok, "n_male"].sum() / s["n_male"].sum(), 3))
print("units with >50% of FlyWire neurons having a hemibrain type:", int((s["hb_share"] > 0.5).sum()))
