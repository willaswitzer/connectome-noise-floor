"""Perturbation-target eligibility for matched units.
primary: 1:1 unit, >=2 neurons per dataset, no technical or sex-linked morphology note.
strict:  primary AND count ratio in [0.8, 1.25] AND >=5 neurons per dataset (sensitivity)."""
import pandas as pd

s = pd.read_csv("data/derived/units_v2_summary.csv")
m = pd.read_csv("data/derived/units_v2_male_neurons.csv")
a = pd.read_feather("data/male_cns/body-annotations-male-cns-v1.0-minconf-0.5.feather")[["bodyId", "matchingNotes"]]
m = m.merge(a, on="bodyId", how="left")
n = m["matchingNotes"].fillna("").str.lower()
m["tech"] = n.str.contains("reconstruction state|too sparse|axon only|hemiorn")
m["morph"] = n.str.contains("difference in morphology|extra layer")
m["var"] = n.str.contains("variab|difference in number")

g = m.groupby("unit_id").agg(
    name=("flywireType", lambda x: "|".join(sorted(set(x)))[:60]),
    superclass=("superclass", lambda x: x.mode().iat[0]),
    tech=("tech", "sum"), morph=("morph", "sum"), var=("var", "sum"))
d = s.merge(g, on="unit_id")
lo = d[["n_male", "n_flywire"]].min(axis=1)
d["in_band"] = d["ratio"].between(0.8, 1.25)
d["base"] = d["is_1to1"] & (lo >= 2)
d["primary"] = d["base"] & (d["tech"] == 0) & (d["morph"] == 0)
d["strict"] = d["primary"] & d["in_band"] & (lo >= 5)
d["primary_hb"] = d["primary"] & (d["hb_share"] > 0.5)
d.to_csv("data/derived/perturbation_candidates.csv", index=False)

sc = ["cb_sensory", "ol_sensory", "descending_neuron", "ascending_neuron", "visual_projection"]
t = d[d["superclass"].isin(sc)].copy()
t["dropped_tech"] = t["base"] & (t["tech"] > 0)
t["dropped_morph"] = t["base"] & (t["morph"] > 0) & (t["tech"] == 0)
print(t.groupby("superclass")[["base", "dropped_tech", "dropped_morph", "primary", "strict", "primary_hb"]].sum())
print("all superclasses: base", int(d["base"].sum()), "| primary", int(d["primary"].sum()),
      "| strict", int(d["strict"].sum()))
print("units with morph note:", int((d["morph"] > 0).sum()), "| names:", list(d.loc[d["morph"] > 0, "name"].head(6)))
print("primary units with a variability note:", int((d["primary"] & (d["var"] > 0)).sum()))
