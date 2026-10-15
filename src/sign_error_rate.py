"""Label- and sign-error rates from the Eckstein Fig. 2A confusion matrices.

Signs follow PREREG A3.4: GABA, Glu negative; ACh, DA, 5-HT, OA positive.
Usage: python src/sign_error_rate.py data/anchors/eckstein_fig2a.csv [data/anchors/composition_v783_predicted.csv]
composition.csv (optional): columns nt,count  -> weights rows by a target
dataset's NT composition (e.g. FlyWire v783 neuron counts by predicted NT).
Without it, rows are weighted by the test-set n (Eckstein's own composition).
"""
import sys
import pandas as pd

NTS = ["ACh", "Glu", "GABA", "DA", "5HT", "OA"]
SIGN = {"ACh": 1, "DA": 1, "5HT": 1, "OA": 1, "Glu": -1, "GABA": -1}

m = pd.read_csv(sys.argv[1], comment="#")
m[NTS] = m[NTS].div(m[NTS].sum(axis=1), axis=0)  # renormalize rounding
comp = None
if len(sys.argv) > 2:
    comp = pd.read_csv(sys.argv[2]).set_index("nt")["count"]

rows = []
for (ds, lvl), g in m.groupby(["dataset", "level"], sort=False):
    g = g.set_index("true_nt")
    w = (comp.reindex(g.index).fillna(0) if comp is not None else g["n"]).astype(float)
    w = w / w.sum()
    label_err = 1 - pd.Series({t: g.loc[t, t] for t in g.index})
    sign_err = pd.Series({t: sum(g.loc[t, p] for p in NTS if SIGN[p] != SIGN[t]) for t in g.index})
    for t in g.index:
        rows.append(dict(dataset=ds, level=lvl, true_nt=t, n=int(g.loc[t, "n"]),
                         label_err=round(label_err[t], 3), sign_err=round(sign_err[t], 3)))
    rows.append(dict(dataset=ds, level=lvl, true_nt="WEIGHTED", n=int(g["n"].sum()),
                     label_err=round((label_err * w).sum(), 3), sign_err=round((sign_err * w).sum(), 3)))
out = pd.DataFrame(rows)
print(out.to_string(index=False))
