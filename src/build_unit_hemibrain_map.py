"""Map matched units to hemibrain bodyIds via FlyWire hemibrain_type strings.
A hemibrain_type string is split on commas; each component is either a neuPrint type
(all traced neurons of that type) or an hb<bodyId> placeholder (that single neuron).
Components matching nothing are written to a separate unresolved table. Does not touch
units_v2_* or the target draw. NT columns come from Eckstein et al. supplemental_data_3."""
import re
import pandas as pd

fw = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
hb = pd.read_parquet("data/hemibrain/neurons.parquet")[["bodyId", "type"]]
nt = pd.read_csv("data/hemibrain/neuron_nt.csv").rename(columns={"bodyid": "bodyId"})
nt = nt[["bodyId", "used_in_analysis", "conf_nt", "conf_nt_p", "top_nt"]]

by_type = hb.dropna(subset=["type"]).groupby("type")["bodyId"].apply(list).to_dict()
all_ids = set(hb["bodyId"])

pairs = (fw.dropna(subset=["hemibrain_type"])[["unit_id", "hemibrain_type"]]
         .drop_duplicates())
rows, unresolved = [], []
for unit_id, s in pairs.itertuples(index=False):
    for comp in (c.strip() for c in s.split(",")):
        if re.fullmatch(r"hb\d+", comp):
            b = int(comp[2:])
            ids, kind = ([b], "hb_id") if b in all_ids else ([], "hb_id")
        else:
            ids, kind = by_type.get(comp, []), "type"
        if not ids:
            unresolved.append((unit_id, s, comp, kind))
        rows += [(unit_id, s, comp, kind, b) for b in ids]

m = pd.DataFrame(rows, columns=["unit_id", "hemibrain_type", "component", "kind", "bodyId"])
m = m.drop_duplicates(["unit_id", "bodyId"]).merge(nt, on="bodyId", how="left")
m["used_in_analysis"] = m["used_in_analysis"].fillna(False).astype(bool)
m.to_csv("data/derived/unit_hemibrain_map.csv", index=False)
u = pd.DataFrame(unresolved, columns=["unit_id", "hemibrain_type", "component", "kind"])
u.to_csv("data/derived/unit_hemibrain_unresolved.csv", index=False)

print("units with a hemibrain string:", pairs["unit_id"].nunique())
print("units mapped to >=1 hemibrain bodyId:", m["unit_id"].nunique())
print("map rows (unit x bodyId):", len(m), "| distinct bodyIds:", m["bodyId"].nunique())
print("component kinds mapped:", m.groupby("kind")["component"].nunique().to_dict())
print("unresolved components:", len(u), "in", u["unit_id"].nunique(), "units;",
      u["kind"].value_counts().to_dict())
print("map rows with an NT row:", int(m["conf_nt"].notna().sum()),
      "| used_in_analysis:", int(m["used_in_analysis"].sum()))
