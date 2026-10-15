"""Exploratory (2026-10-15): neuron-level transmitter confusion of the v783 predictions
(annotations `top_nt`) against literature-known transmitters (`known_nt`), FlyWire v783.

Known transmitter per neuron: the classical transmitters (ACh, Glu, GABA, DA, 5-HT, OA) named
positively in `known_nt` (entries ending "-negative" ignored); neurons naming none or more than
one are excluded (as A3.4 tier 0). Histamine neurons are excluded (not a predicted class).
Reports row-normalized confusion, label and sign error (A3.4 signs), overall, with
top_nt_conf >= 0.5 (the A3.4 tier-1 floor), and by optic lobe vs. rest; plus the share of
cell types whose members' predictions are not all correct (clustering).
Usage: python src/nt_v783_known_check.py  ->  results/arm2/linear/nt_v783_known_check.txt
"""
import re
import pandas as pd

NT = {"acetylcholine": "ACh", "glutamate": "Glu", "gaba": "GABA", "dopamine": "DA",
      "serotonin": "5HT", "octopamine": "OA"}
SIGN = {"ACh": 1, "DA": 1, "5HT": 1, "OA": 1, "Glu": -1, "GABA": -1}
ORDER = ["ACh", "Glu", "GABA", "DA", "5HT", "OA"]


def known(s):
    if pd.isna(s):
        return None
    toks = {t.strip().lower() for t in re.split(r"[;,]", s)}
    pos = {NT[t] for t in toks if t in NT}
    return pos.pop() if len(pos) == 1 else None


a = pd.read_csv("data/flywire_783/annotations_v3_neurons.tsv", sep="\t", low_memory=False)
a["known"] = a["known_nt"].map(known)
a["pred"] = a["top_nt"].map(NT)
d = a.dropna(subset=["known", "pred"]).copy()
d["label_err"] = d["known"] != d["pred"]
d["sign_err"] = d["known"].map(SIGN) != d["pred"].map(SIGN)
d["ol"] = d["super_class"].isin(["optic", "visual_projection", "visual_centrifugal"])

out = []
def report(name, x):
    out.append(f"\n== {name}: {len(x)} neurons, {x['cell_type'].nunique()} types")
    cm = pd.crosstab(x["known"], x["pred"], normalize="index").reindex(index=ORDER, columns=ORDER).fillna(0)
    n = x["known"].value_counts().reindex(ORDER).fillna(0).astype(int)
    out.append(cm.round(3).assign(n=n).to_string())
    by = x.groupby("known")[["label_err", "sign_err"]].mean().reindex(ORDER)
    out.append("per known NT:\n" + by.round(3).to_string())
    out.append(f"overall label error {x.label_err.mean():.3f}, sign error {x.sign_err.mean():.3f}")
    t = x.dropna(subset=["cell_type"]).groupby("cell_type").agg(n=("sign_err", "size"),
                                                                  any_err=("sign_err", "any"),
                                                                  all_err=("sign_err", "all"))
    t = t[t.n >= 2]
    out.append(f"types with >=2 neurons: {len(t)}; with any sign error {t.any_err.mean():.3f}; "
               f"all members sign-wrong {t.all_err.mean():.3f}; "
               f"of types with any error, all wrong {t.all_err.sum() / max(t.any_err.sum(), 1):.3f}")

report("all known", d)
report("all known, top_nt_conf >= 0.5", d[d.top_nt_conf >= 0.5])
report("central (non-optic) known", d[~d.ol])
report("central known, top_nt_conf >= 0.5", d[~d.ol & (d.top_nt_conf >= 0.5)])
txt = "\n".join(out)
print(txt)
open("results/arm2/linear/nt_v783_known_check.txt", "w").write(txt + "\n")
