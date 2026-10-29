"""Exploratory (2026-10-27): the left-right noise floor under the glutamate-modified networks.

Glutamate's sign is the same unknown in every animal, so if part of glutamate is excitatory it is
excitatory in both hemispheres. This recomputes the left-right disagreement (R0) inside each
glutamate-modified network, from the left- and right-side stimulation runs it already has, and
compares drift (glutamate network vs. the original) with that disagreement as well as with the
original D_R0.

LIF: R0 rho exactly as analyze_lif.py rho_10 (10 trials, r0 scheme, seed rows removed), per draw;
     check: on the Amendment-7 runs it reproduces rho_10 in results/lif/primary/rho_case_study.csv.
Linear: R0 and R0ipsi rho exactly as Arm 1 (analyze_linear.pair_rho on the half_* views), per draw;
     check: on the base run it reproduces the Arm-1 D.
D under a condition = 1 - rho, averaged over the 3 draws per perturbation; summaries are medians
over the same perturbations as the drift, with a 95% bootstrap CI of D_glu - drift
(1,000 resamples, default_rng(20261027)).
Outputs: results/arm2/{lif,linear}/glu_noise_floor.csv
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze_linear import pair_rho
import arm2_lif_analyze as la
from arm2_analyze import D as DLIN, views, seed_rows, rung_disagreement

RNG = 20261027
LIF_CONDS = ["glu_post_f10", "glu_post_f25", "glu_pre_f25"]
LIN_CONDS = [f"glu_{s}_f{f}" for s in ("pre", "post") for f in (10, 25, 50)]


def lif_r0_rho(d, p):
    kl, kr = la.lookup(p, "R0", "left"), la.lookup(p, "R0", "right")
    if not (kl and kr and (d / f"{kl}.parquet").exists() and (d / f"{kr}.parquet").exists()):
        return np.nan
    A, B = la.vec(d / f"{kl}.parquet", "r0", "left"), la.vec(d / f"{kr}.parquet", "r0", "right")
    idx = A.index.intersection(B.index)
    seeds = {f"{u}|{s}" for u in la.seed_units[p] for s in ("ipsi", "contra")}
    return pair_rho(A.loc[idx].to_numpy(), B.loc[idx].to_numpy(), idx.isin(list(seeds)))[0]


def boot(x, y, rng):
    n = len(x)
    bd = np.array([np.median(x[i]) - np.median(y[i]) for i in (rng.integers(0, n, n) for _ in range(1000))])
    return float(np.quantile(bd, .025)), float(np.quantile(bd, .975))


def lif():
    ref = pd.read_csv("results/lif/primary/rho_case_study.csv").set_index(["perturbation", "rung"])["rho_10"]
    chk = {p: lif_r0_rho(la.A7, p) for p in la.seed_units}
    diff = max(abs(chk[p] - ref[(p, "R0")]) for p in chk if pd.notna(chk[p]) and pd.notna(ref.get((p, "R0"))))
    print(f"LIF check: max |rho - rho_10| on Amendment-7 runs = {diff:.2e}")
    drift = pd.read_csv(la.OUT / "glu_drift_per_rep.csv")
    drift = drift[drift.rung == "R0"].groupby(["condition", "perturbation"]).drift.mean()
    rng = np.random.default_rng(RNG)
    rows = []
    for c in LIF_CONDS:
        Dg = {p: np.nanmean([1 - lif_r0_rho(la.RUNS / c / f"k{k}", p) for k in range(3)]) for p in la.seed_units}
        for row, ps in (("seed_specific", la.SPECIFIC), ("pn_shared", la.PN)):
            ok = [p for p in ps if pd.notna(Dg[p]) and pd.notna(drift.get((c, p))) and pd.notna(ref.get((p, "R0")))]
            dg = np.array([Dg[p] for p in ok]); d0 = np.array([1 - ref[(p, "R0")] for p in ok])
            y = np.array([drift[(c, p)] for p in ok])
            lo, hi = boot(dg, y, rng)
            rows.append({"condition": c, "row": row, "n": len(ok), "median_drift": float(np.median(y)),
                         "D_R0_original": float(np.median(d0)), "D_R0_glu_network": float(np.median(dg)),
                         "drift_over_D_original": float(np.median(y) / np.median(d0)),
                         "drift_over_D_glu_network": float(np.median(y) / np.median(dg)),
                         "ci_lo_Dglu_minus_drift": lo, "ci_hi_Dglu_minus_drift": hi})
    out = pd.DataFrame(rows)
    out.to_csv(la.OUT / "glu_noise_floor.csv", index=False)
    return out


def lin_r0(dirpath, seeds, draw):
    v = views(dirpath)
    out = {}
    for rung in ("R0", "R0ipsi"):
        A, B = v[rung]
        idx = A.index.intersection(B.index)
        for p in draw[rung]:
            if p in A.columns:
                out[(rung, p)] = pair_rho(A.loc[idx, p].to_numpy(), B.loc[idx, p].to_numpy(),
                                          idx.isin(list(seeds[rung][p])))[0]
    return out


def linear():
    seeds, draw = seed_rows()
    gd = pd.read_csv("results/arm2/linear/glu_drift.csv")
    rng = np.random.default_rng(RNG)
    rows = []
    for model in ("norect", "rect"):
        base = lin_r0(DLIN / model / "base" / "k0", seeds, draw)
        arm1 = rung_disagreement(model).set_index(["rung", "perturbation"])["dis"]
        diff = max(abs((1 - base[k]) - arm1[k]) for k in base if pd.notna(base[k]) and pd.notna(arm1.get(k)))
        print(f"linear {model} check: max |D - Arm-1 D| on base run = {diff:.1e}")
        for c in LIN_CONDS:
            per = [lin_r0(DLIN / model / c / f"k{k}", seeds, draw) for k in range(3)]
            for rung in ("R0", "R0ipsi"):
                dr = gd[(gd.model == model) & (gd.condition == c) & (gd.rung == rung)].set_index("perturbation").drift
                ok = [p for p in draw[rung] if pd.notna(dr.get(p)) and pd.notna(base.get((rung, p)))]
                dg = np.array([np.nanmean([1 - q[(rung, p)] for q in per]) for p in ok])
                d0 = np.array([1 - base[(rung, p)] for p in ok])
                y = dr.loc[ok].to_numpy(float)
                m = ~np.isnan(dg)
                lo, hi = boot(dg[m], y[m], rng)
                rows.append({"model": model, "condition": c, "rung": rung, "n": int(m.sum()),
                             "median_drift": float(np.median(y[m])), "D_original": float(np.median(d0[m])),
                             "D_glu_network": float(np.median(dg[m])),
                             "drift_over_D_original": float(np.median(y[m]) / np.median(d0[m])),
                             "drift_over_D_glu_network": float(np.median(y[m]) / np.median(dg[m])),
                             "ci_lo_Dglu_minus_drift": lo, "ci_hi_Dglu_minus_drift": hi})
    out = pd.DataFrame(rows)
    out.to_csv("results/arm2/linear/glu_noise_floor.csv", index=False)
    return out


if __name__ == "__main__":
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(lif().round(3).to_string(index=False))
        print(linear().round(3).to_string(index=False))
