"""Exploratory (2026-10-15): LIF drift for the confusion-matrix NT conditions
(run_arm2_lif.py --ntcm ntcm_hemi_neuron,ntcm_hemi_type; k = 0-2; 10 trials; seed-exempt flips).

Drift, rows, reference (first 10 trials of the Amendment-7 FlyWire runs) and D_R (1 - rho_10)
exactly as arm2_lif_analyze.py. For each condition x rung (R0, R1; R2 descriptive) x row
(seed_specific: 9 perturbations; pn_shared: 5 PN seeds): median drift (drift averaged over the
3 draws first), D_R, drift / D_R, D_R - drift with a 95% bootstrap CI (1,000 resamples of
perturbations, default_rng(20261015)), and m_eq on the A8.1 LIF flip curve
(results/arm2/lif/curves.csv).
Output: results/arm2/lif/ntcm_summary.csv, ntcm_drift_per_rep.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import arm2_lif_analyze as la
from nt_realistic_analyze import m_eq

RNG = 20261015


def main():
    rows = []
    for cdir in sorted(la.RUNS.glob("ntcm_*")):
        for rdir in sorted(cdir.iterdir()):
            if rdir.is_dir() and any(rdir.glob("*.parquet")):
                t = la.drift_for_dir(rdir)
                t["condition"], t["rep"] = cdir.name, rdir.name
                t["n_runs"] = len(list(rdir.glob("*.parquet")))
                rows.append(t)
    dr = pd.concat(rows, ignore_index=True)
    dr.to_csv(la.OUT / "ntcm_drift_per_rep.csv", index=False)
    d = dr.groupby(["condition", "rung", "perturbation"])["drift"].mean().reset_index()
    nrep = dr.groupby("condition")["rep"].nunique()
    rho = pd.read_csv("results/lif/primary/rho_case_study.csv").set_index(["perturbation", "rung"])["rho_10"]
    curves = pd.read_csv(la.OUT / "curves.csv")
    rng = np.random.default_rng(RNG)
    out = []
    for cond in sorted(d.condition.unique()):
        for rung in ("R0", "R1", "R2"):
            for row, ps in (("seed_specific", la.SPECIFIC), ("pn_shared", la.PN)):
                s = d[(d.condition == cond) & (d.rung == rung)].set_index("perturbation")["drift"]
                ok = [p for p in ps if p in s.index and pd.notna(s[p]) and (p, rung) in rho.index
                      and pd.notna(rho[(p, rung)])]
                if not ok:
                    continue
                y = s.loc[ok].to_numpy(float)
                x = (1 - rho.loc[[(p, rung) for p in ok]]).to_numpy(float)
                n = len(ok)
                bd = np.array([np.median(x[i]) - np.median(y[i])
                               for i in (rng.integers(0, n, n) for _ in range(1000))])
                cv = curves[(curves.type == "flip") & (curves.rung == rung) & (curves.row == row)]
                out.append({"condition": cond, "rung": rung, "row": row, "n": n,
                            "n_reps": int(nrep[cond]), "median_drift": float(np.median(y)),
                            "D_R": float(np.median(x)), "drift_over_D": float(np.median(y) / np.median(x)),
                            "D_minus_drift": float(np.median(x) - np.median(y)),
                            "ci_lo": float(np.quantile(bd, 0.025)), "ci_hi": float(np.quantile(bd, 0.975)),
                            "m_eq_lif_flip_curve": m_eq(cv, float(np.median(y))) if len(cv) else "",
                            "target": "descriptive" if rung == "R2" else "break-even"})
    res = pd.DataFrame(out)
    res.to_csv(la.OUT / "ntcm_summary.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
