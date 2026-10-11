"""Arm-2 real-error condition analysis (A8.6; implementation note 2026-10-10).
drift_real per drawn target and rung (same definition as arm2_analyze.py), m_eq on each family
curve, and measured per-neuron loss/gain. Outputs: results/arm2/linear/v630_*.csv/json"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import arm2_analyze as a

OUT = a.OUT
PHOTO = {"t:6298", "t:6299"}   # R1-6 and R7: seed sets < 80% mapped


def m_eq(ms, med, x):
    ms, med = np.asarray(ms, float), np.asarray(med, float)
    if x < med[0]:
        return "< %g" % ms[0]
    hit = np.flatnonzero(med >= x)
    if len(hit) == 0:
        return "> %g" % ms[-1]
    i = hit[0]
    if i == 0:
        return float(ms[0])
    f = (x - med[i - 1]) / (med[i] - med[i - 1])
    return float(np.exp(np.log(ms[i - 1]) + f * (np.log(ms[i]) - np.log(ms[i - 1]))))


def main():
    seeds, draw = a.seed_rows()
    rows = []
    for model in ("rect", "norect"):
        other, base = a.other_side(model)
        for cond in ("v630", "v630_shiusign"):
            t = a.drift_table(a.D / model / cond / "k0", base, other, seeds)
            t["model"], t["condition"] = model, cond
            rows.append(t)
    dr = pd.concat(rows, ignore_index=True)
    dr.to_csv(OUT / "v630_drift.csv", index=False)

    curves = pd.read_csv(OUT / "curves.csv")
    res = []
    for model in ("rect", "norect"):
        for cond in ("v630", "v630_shiusign"):
            for rung in ("R0", "R0ipsi", "R1", "R2"):
                s = dr[(dr.model == model) & (dr.condition == cond) & (dr.rung == rung)]
                s = s[s.perturbation.isin(draw[rung])]
                for excl in (False, True):
                    v = s[~s.perturbation.isin(PHOTO)] if excl else s
                    med = float(np.nanmedian(v.drift))
                    r = {"model": model, "condition": cond, "rung": rung, "exclude_R1-6_R7": excl,
                         "n": int(v.drift.notna().sum()), "median_drift_real": med}
                    fams = ["flip"] if model == "norect" else ["combined", "lopsided", "miss", "split", "merge", "spur"]
                    for fam in fams:
                        fm = "norect" if fam == "flip" else "rect"
                        c = curves[(curves.family == fam) & (curves.model == fm) & (curves.rung == rung)].sort_values("m")
                        if len(c):
                            r[f"m_eq_{fam}"] = m_eq(c.m, c.median_drift, med)
                    res.append(r)
    pd.DataFrame(res).to_csv(OUT / "v630_meq.csv", index=False)

    lg = pd.read_csv("data/derived/arm2/v630_lossgain.csv")
    b = lg[lg.in_both & (lg.syn_v783 > 0) & (lg.syn_v630 > 0)]
    big = b[b.syn_v783 >= 100]
    out = {"n_neurons": int(len(b)),
           "loss_pooled": float(1 - b.matched.sum() / b.syn_v783.sum()),
           "gain_pooled": float(1 - b.matched.sum() / b.syn_v630.sum()),
           "loss_median_neuron_ge100": float(np.median(1 - big.matched / big.syn_v783)),
           "gain_median_neuron_ge100": float(np.median(1 - big.matched / big.syn_v630))}
    (OUT / "v630_lossgain.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
