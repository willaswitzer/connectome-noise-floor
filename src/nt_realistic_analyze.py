"""Exploratory (2026-10-15): drift for the confusion-matrix NT conditions (run_nt_realistic.py).

Same drift, rows, seeds and D_R as arm2_analyze.py (A8.4). For each condition x model x rung:
  - median drift over drawn targets (drift averaged over the 5 draws first), D_R, drift / D_R
  - D_R - median drift with a 95% bootstrap CI (1,000 resamples of targets, default_rng(20261015))
  - m_eq: the median drift placed on the existing flip curve of the same model and rung
    (results/arm2/linear/curves.csv), linear in log m; "< 0.5" below the first grid point
  - mean realized flip fraction (of signed nodes, the A8.1 denominator)
Outputs: results/arm2/linear/ntcm_drift.csv, ntcm_summary.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from arm2_analyze import D, OUT, seed_rows, other_side, drift_table, rung_disagreement
from run_nt_realistic import CONDITIONS, MODELS

RNG = 20261015


def m_eq(curve, y):
    c = curve.sort_values("m")
    ms, md = c["m"].to_numpy(float), c["median_drift"].to_numpy(float)
    if y < md[0]:
        return "< %g" % ms[0]
    if y > md[-1]:
        return "> %g" % ms[-1]
    i = int(np.flatnonzero(md >= y)[0])
    if i == 0:
        return "%.2f" % ms[0]
    f = (y - md[i - 1]) / (md[i] - md[i - 1])
    return "%.2f" % np.exp(np.log(ms[i - 1]) + f * (np.log(ms[i]) - np.log(ms[i - 1])))


def main():
    seeds, draw = seed_rows()
    tabs = []
    for model in MODELS:
        other, base = other_side(model)
        for cond in CONDITIONS:
            for kdir in sorted((D / model / cond).glob("k*")):
                if not (kdir / "meta.json").exists():
                    continue
                t = drift_table(kdir, base, other, seeds)
                meta = json.loads((kdir / "meta.json").read_text())
                t["model"], t["condition"], t["k"] = model, cond, int(kdir.name[1:])
                t["flip_frac_signed"] = meta["flip_frac_signed"]
                tabs.append(t)
    dr = pd.concat(tabs, ignore_index=True)
    g = dr.groupby(["model", "condition", "rung", "perturbation"])
    drift = g["drift"].agg(lambda v: np.nanmean(v) if v.notna().any() else np.nan).reset_index()
    drift["n_draws"] = g["k"].nunique().values
    drift.to_csv(OUT / "ntcm_drift.csv", index=False)

    curves = pd.read_csv(OUT / "curves.csv")
    rng = np.random.default_rng(RNG)
    rows = []
    for model in MODELS:
        dis = rung_disagreement(model)
        for cond in CONDITIONS:
            ff = dr[(dr.model == model) & (dr.condition == cond)].groupby("k")["flip_frac_signed"].first()
            for rung in ("R0", "R1", "R2"):
                s = drift[(drift.model == model) & (drift.condition == cond) & (drift.rung == rung)]
                s = s.set_index("perturbation")["drift"]
                dd = dis[dis.rung == rung].set_index("perturbation")["dis"]
                ok = [k for k in draw[rung] if k in s.index and pd.notna(s[k]) and pd.notna(dd.get(k))]
                x, y = dd.loc[ok].to_numpy(float), s.loc[ok].to_numpy(float)
                n = len(ok)
                if n == 0:
                    continue
                bd = np.array([np.median(x[i]) - np.median(y[i])
                               for i in (rng.integers(0, n, n) for _ in range(1000))])
                cv = curves[(curves.family == "flip") & (curves.model == model) & (curves.rung == rung)]
                rows.append({"condition": cond, "model": model, "rung": rung, "n_targets": n,
                             "flip_frac_signed": float(ff.mean()),
                             "median_drift": float(np.median(y)), "p90_drift": float(np.quantile(y, 0.9)),
                             "D_R": float(np.median(x)), "drift_over_D": float(np.median(y) / np.median(x)),
                             "D_minus_drift": float(np.median(x) - np.median(y)),
                             "ci_lo": float(np.quantile(bd, 0.025)), "ci_hi": float(np.quantile(bd, 0.975)),
                             "m_eq_flip_curve": m_eq(cv, float(np.median(y))) if len(cv) else ""})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "ntcm_summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(out.round(4).to_string(index=False))
    combined_ntfix(seeds, draw)


def combined_ntfix(seeds, draw):
    """run_combined_ntfix.py: break-even of the combined curve and the P5/P6 statistics with the
    flip component at the confusion-matrix rate (rectified model, as A8.4-A8.5)."""
    from arm2_analyze import breakeven
    from run_combined_ntfix import conditions as ntfix_conditions
    conds = ntfix_conditions()
    other, base = other_side("rect")
    tabs = []
    for c in conds:
        for kdir in sorted((D / "rect" / c).glob("k*")):
            if (kdir / "meta.json").exists():
                t = drift_table(kdir, base, other, seeds)
                t["condition"], t["k"] = c, int(kdir.name[1:])
                tabs.append(t)
    if not tabs:
        return
    dr = pd.concat(tabs, ignore_index=True)
    drift = dr.groupby(["condition", "rung", "perturbation"])["drift"].agg(
        lambda v: np.nanmean(v) if v.notna().any() else np.nan).reset_index()
    drift.to_csv(OUT / "ntfix_drift.csv", index=False)
    dis = rung_disagreement("rect")
    curve = [c for c in conds if c.startswith("combinedntfix")]
    ms = [conds[c][1] for c in curve]
    rng = np.random.default_rng(RNG)
    rows = []
    for rung in ("R0", "R1", "R2"):
        piv = drift[drift.rung == rung].pivot(index="perturbation", columns="condition", values="drift")
        dd = dis[dis.rung == rung].set_index("perturbation")["dis"]
        ok = [k for k in draw[rung] if k in piv.index and pd.notna(dd.get(k))
              and piv.loc[k, curve].notna().any()]
        res = breakeven(piv.loc[ok, curve].to_numpy(float), dd.loc[ok].to_numpy(float), ms, rng)
        rows.append({"rung": rung, "D_R": res["D_R"], "m_star": res["m_star"], "flag": res["flag"],
                     "ci_lo": res["ci_lo"], "ci_hi": res["ci_hi"], "n_targets": res["n_targets"],
                     "abs_rate": res["m_star"] * 0.008,
                     **{f"drift_m{m:g}": v for m, v in zip(ms, res["median_drift"])}})
    be = pd.DataFrame(rows)
    be.to_csv(OUT / "ntfix_breakeven.csv", index=False)
    # P5 / P6 statistics on R0ipsi, as A8.5 (no Holm: exploratory)
    dd = dis[dis.rung == "R0ipsi"].set_index("perturbation")["dis"]
    pred = {}
    for name, c in (("P5_ntfix", "combinedntfix_m1"), ("P6_ntfix", "unproofreadntfix")):
        s = drift[(drift.rung == "R0ipsi") & (drift.condition == c)].set_index("perturbation")["drift"]
        ok = [k for k in draw["R0ipsi"] if k in s.index and pd.notna(s[k]) and pd.notna(dd.get(k))]
        if not ok:
            continue
        x, y = dd.loc[ok].to_numpy(float), s.loc[ok].to_numpy(float)
        bd = np.array([np.median(x[i]) - np.median(y[i])
                       for i in (np.random.default_rng(RNG).integers(0, len(ok), (1000, len(ok))))])
        pred[name] = {"condition": c, "n_targets": len(ok), "D_R0ipsi": float(np.median(x)),
                      "median_drift": float(np.median(y)), "diff": float(np.median(x) - np.median(y)),
                      "ci95": [float(v) for v in np.quantile(bd, [0.025, 0.975])]}
    (OUT / "ntfix_predictions.json").write_text(json.dumps(pred, indent=1))
    print(be.round(4).to_string(index=False))
    print(json.dumps(pred, indent=1))


if __name__ == "__main__":
    main()
