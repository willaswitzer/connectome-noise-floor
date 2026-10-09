"""Arm-2 linear drift, break-even and P5/P6 (PREREG A8.4-A8.5).

Inputs: data/derived/arm2/linear/<model>/<condition>/k<k>/ (run_arm2_linear.py), the unperturbed
"base" run of the same model, and the Arm-1 other-dataset vectors (data/derived/linear/<cond>/).

Drift for target t on rung R (A8.4): 1 - rho(perturbed, unperturbed FlyWire) on R's own A4.1 rows:
  rows present in both datasets of R and non-missing in the other dataset, seed rows removed,
  rows zero in both the perturbed and unperturbed vector removed, NA if < 10 rows (A2.3).
  - A perturbed vector that is identically zero while the unperturbed one is not gets drift 1
    (rho taken as 0): the error abolished the response.
  - R0 and R0ipsi compare two FlyWire vectors (left- and right-side stimulation). Both are
    perturbed; drift is the mean of the left-stim and right-stim drifts.
  Drift is averaged over injection seeds (mean of non-NA values).
Rung disagreement D_R = median over the same targets of 1 - raw rho_R (Arm-1 A4.1 tables).
Break-even m*_R = smallest m where median drift >= D_R, linear in log m between grid points.
95% CI: 1,000 bootstrap resamples of targets (drift and D_R jointly), default_rng(20260930).

Models: rect (primary for all but flip), norect (flip primary, compared with norect D_R),
rect_ex (flip sensitivity, rect D_R). Targets: drawn targets only (R0, R0ipsi, R1: 134;
R2: 58 Tier B); case-study perturbations are written to the drift table but not used here.

Outputs: results/arm2/linear/drift.csv (per target x condition x model x rung),
         results/arm2/linear/curves.csv, breakeven.csv, predictions.json"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze_linear import pair_rho, r0_reindex, side_rows
import inject_errors as ie

D = Path("data/derived/arm2/linear")
OUT = Path("results/arm2/linear")
ARM1 = {"rect": "a0.5_K5", "rect_ex": "a0.5_K5", "norect": "a0.5_K5_norect"}
RUNGS = ["R0", "R0ipsi", "R1", "R2"]
B_BOOT = 1000
RNG_SEED = 20260930


# ----------------------------------------------------------------------------- rows and seeds

def seed_rows():
    targets = pd.read_csv("results/prereg/targets_DRAFT.csv")
    fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
    case = pd.read_csv("results/prereg/case_study.csv")
    r2map = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]
    units = {f"t:{u}": {int(u)} for u in targets["unit_id"]}
    for _, row in case.iterrows():
        units[f"c:{row['name']}"] = set(fwu.loc[fwu["cell_type"].isin(row["flywire_types"].split("|")),
                                                "unit_id"].astype(int))
    seeds = {
        "R1": {k: {str(u) for u in v} for k, v in units.items()},
        "R2": {k: {str(int(r2map[u])) for u in v} for k, v in units.items()},
        "R0": {k: {f"{u}|{s}" for u in v for s in ("ipsi", "contra")} for k, v in units.items()},
        "R0ipsi": {k: {str(u) for u in v} for k, v in units.items()},
    }
    draw = {"R0": [f"t:{u}" for u in targets["unit_id"]],
            "R0ipsi": [f"t:{u}" for u in targets["unit_id"]],
            "R1": [f"t:{u}" for u in targets["unit_id"]],
            "R2": [f"t:{u}" for u in targets.loc[targets["tierB"], "unit_id"]]}
    return seeds, draw


R0_UNITS = None


def r0_rows(df):
    """half_* rows of R0-usable units (A2.4), as run_linear.py's R0 readout."""
    global R0_UNITS
    if R0_UNITS is None:
        R0_UNITS = {str(u) for u in pd.read_csv("data/derived/r0_units.csv").query("r0_usable")["unit_id"]}
    keep = df.index.str.split("|").str[0].isin(R0_UNITS)
    return df[keep]


def views(d):
    """Load one run directory into the rung views: {rung: [(vector table, side tag), ...]}.
    R1/R2 have one view; R0/R0ipsi have a left-stim and a right-stim view."""
    rd = lambda f: pd.read_parquet(d / f"{f}.parquet").astype(np.float64)
    hl, hr = rd("half_left"), rd("half_right")
    return {"R1": [rd("R1")], "R2": [rd("R2")],
            "R0": [r0_reindex(r0_rows(hl), "left"), r0_reindex(r0_rows(hr), "right")],
            "R0ipsi": [side_rows(hl, "left"), side_rows(hr, "right")]}


def other_side(model):
    """For each rung, the table whose non-missing rows define R's A4.1 row set on the other side."""
    c = ARM1[model]
    base = views(D / ("rect" if model != "norect" else "norect") / "base" / "k0")
    return {"R1": [pd.read_parquet(f"data/derived/linear/{c}/R1_male_cns.parquet")],
            "R2": [pd.read_parquet(f"data/derived/linear/{c}/R2_hemibrain.parquet")],
            # R0 / R0ipsi: the other side is the opposite-side stimulation of FlyWire itself
            "R0": [base["R0"][1], base["R0"][0]],
            "R0ipsi": [base["R0ipsi"][1], base["R0ipsi"][0]]}, base


# ----------------------------------------------------------------------------- drift

def drift_one(p, u, other, seedmask):
    """p, u, other: aligned 1-D arrays; returns drift or NaN."""
    ex = seedmask | np.isnan(other)
    if np.isnan(p).all() or np.isnan(u).all():
        return np.nan
    keep = ~ex & ~np.isnan(p) & ~np.isnan(u)
    if keep.sum() >= 10 and np.all(p[keep] == 0) and np.any(u[keep] != 0):
        return 1.0
    r, n, why = pair_rho(p, u, ex)
    return 1.0 - r if not np.isnan(r) else np.nan


def drift_table(run_dir, base, other, seeds):
    P = views(run_dir)
    rows = []
    for rung in RUNGS:
        per = []
        for Pv, Uv, Ov in zip(P[rung], base[rung], other[rung]):
            idx = Uv.index.intersection(Ov.index)
            cols = [c for c in Uv.columns if c in Ov.columns]
            Pa = Pv.reindex(index=idx, columns=cols).to_numpy(float)
            Ua = Uv.loc[idx, cols].to_numpy(float)
            Oa = Ov.loc[idx, cols].to_numpy(float)
            per.append({c: drift_one(Pa[:, j], Ua[:, j], Oa[:, j], idx.isin(list(seeds[rung][c])))
                        for j, c in enumerate(cols)})
        for c in per[0]:
            vals = [d[c] for d in per if c in d]
            v = np.nanmean(vals) if not np.all(np.isnan(vals)) else np.nan
            rows.append({"rung": rung, "perturbation": c, "drift": v})
    return pd.DataFrame(rows)


def rung_disagreement(model, seeds=None, base=None, other=None):
    """1 - raw rho_R per perturbation, read from the Arm-1 A4.1 tables (results/linear/<cond>/
    rho_<rung>.csv). These were computed from the float64 Arm-1 vectors; recomputing from the
    float32 Arm-2 base run reproduces them to <= 1e-3 (ties among tiny R0 values)."""
    rows = []
    for rung in RUNGS:
        r = pd.read_csv(f"results/linear/{ARM1[model]}/rho_{rung}.csv")
        rows.append(pd.DataFrame({"rung": rung, "perturbation": r["perturbation"],
                                  "rho_R": r["rho"], "dis": 1 - r["rho"]}))
    return pd.concat(rows, ignore_index=True)


# ----------------------------------------------------------------------------- break-even

def crossing(ms, med, D):
    """Smallest m with med >= D, linear in log m. Returns (m*, flag) with flag in
    {'', 'le_min', 'gt_max'}."""
    ms, med = np.asarray(ms, float), np.asarray(med, float)
    hit = np.flatnonzero(med >= D)
    if len(hit) == 0:
        return np.inf, "gt_max"
    i = hit[0]
    if i == 0:
        return ms[0], "le_min"
    x0, x1 = np.log(ms[i - 1]), np.log(ms[i])
    y0, y1 = med[i - 1], med[i]
    f = (D - y0) / (y1 - y0) if y1 != y0 else 0.0
    return float(np.exp(x0 + f * (x1 - x0))), ""


def breakeven(curve_drift, dis, ms, rng):
    """curve_drift: targets x len(ms) array (NaN allowed); dis: per-target 1 - rho_R."""
    med = np.nanmedian(curve_drift, axis=0)
    D = np.nanmedian(dis)
    m_star, flag = crossing(ms, med, D)
    n = len(dis)
    boots = np.empty(B_BOOT)
    for b in range(B_BOOT):
        i = rng.integers(0, n, n)
        boots[b] = crossing(ms, np.nanmedian(curve_drift[i], axis=0), np.nanmedian(dis[i]))[0]
    # "> max m" is inf; quantiles are taken with inf as a sentinel and mapped back
    big = 1e12
    lo, hi = np.quantile(np.where(np.isinf(boots), big, boots), [0.025, 0.975])
    lo, hi = (np.inf if lo >= big else lo), (np.inf if hi >= big else hi)
    return {"D_R": float(D), "m_star": m_star, "flag": flag, "ci_lo": float(lo), "ci_hi": float(hi),
            "n_targets": int(n), "median_drift": [float(x) for x in med]}


FAMILIES = {  # family -> (condition prefix, anchor used for the absolute rate)
    "miss": ("miss", 0.008), "split": ("split", 0.008), "merge": ("merge", 0.008),
    "spur": ("spur", 0.05), "flip": ("flip", 0.06),
    "combined": ("combined", 0.008), "lopsided": ("lopsided", 0.008),
}


def compute_drift(seeds, grid):
    tables = []
    for model in ("rect", "norect", "rect_ex"):
        other, base = other_side(model)
        for cdir in sorted((D / model).iterdir()):
            if cdir.name == "base" or cdir.name not in grid:
                continue
            for kdir in sorted(cdir.glob("k*")):
                if not (kdir / "meta.json").exists():
                    continue
                t = drift_table(kdir, base, other, seeds)
                t["model"], t["condition"], t["k"] = model, cdir.name, int(kdir.name[1:])
                tables.append(t)
        print(f"{model}: drift computed", flush=True)
    dr = pd.concat(tables, ignore_index=True)
    dr.to_csv(OUT / "drift_per_seed.csv", index=False)
    g = dr.groupby(["model", "condition", "rung", "perturbation"])
    drift = g["drift"].agg(lambda v: np.nanmean(v) if v.notna().any() else np.nan).reset_index()
    drift["n_seeds"] = g["k"].nunique().values
    drift.to_csv(OUT / "drift.csv", index=False)
    return drift


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    reuse = "--reuse-drift" in sys.argv
    seeds, draw = seed_rows()
    grid = {lab: (kind, m, t) for lab, kind, m, t, _, _ in ie.grid()}
    dis_tabs = {m: rung_disagreement(m) for m in ("rect", "norect", "rect_ex")}
    if reuse:
        drift = pd.read_csv(OUT / "drift.csv")
    else:
        drift = compute_drift(seeds, grid)
    pd.concat([d.assign(model=m) for m, d in dis_tabs.items()]).to_csv(OUT / "rung_disagreement.csv",
                                                                       index=False)

    rng = np.random.default_rng(RNG_SEED)
    rows, curves = [], []
    for fam, (prefix, anchor) in FAMILIES.items():
        models = ["norect", "rect_ex"] if fam == "flip" else ["rect"]
        conds = sorted([l for l, (kind, m, t) in grid.items()
                        if (kind == "single" and t == fam) or kind == fam], key=lambda l: grid[l][1])
        for model in models:
            dis = dis_tabs[model]
            for rung in ("R0", "R1", "R2"):
                tg = [k for k in draw[rung]]
                dsub = drift[(drift.model == model) & (drift.rung == rung)]
                piv = dsub.pivot(index="perturbation", columns="condition", values="drift")
                have = [c for c in conds if c in piv.columns]
                if not have:
                    continue
                ms = [grid[c][1] for c in have]
                dd = dis[dis.rung == rung].set_index("perturbation")["dis"]
                ok = [k for k in tg if k in piv.index and k in dd.index and pd.notna(dd[k])
                      and piv.loc[k, have].notna().any()]
                cd = piv.loc[ok, have].to_numpy(float)
                res = breakeven(cd, dd.loc[ok].to_numpy(float), ms, rng)
                res.update(family=fam, model=model, rung=rung, m_grid=ms, conditions=have,
                           complete=len(have) == len(conds),
                           abs_rate=res["m_star"] * anchor if np.isfinite(res["m_star"]) else np.inf)
                rows.append(res)
                for m, c, v in zip(ms, have, res["median_drift"]):
                    curves.append({"family": fam, "model": model, "rung": rung, "m": m,
                                   "condition": c, "median_drift": v, "D_R": res["D_R"]})
    pd.DataFrame(rows).to_csv(OUT / "breakeven.csv", index=False)
    pd.DataFrame(curves).to_csv(OUT / "curves.csv", index=False)

    # P5 / P6 (A8.5): rectified model, R0 in its A5.3 form (R0ipsi)
    rng = np.random.default_rng(RNG_SEED)
    dis = dis_tabs["rect"]
    dd = dis[dis.rung == "R0ipsi"].set_index("perturbation")["dis"]
    pred = {}
    for name, cond in (("P5", "combined_m1"), ("P6", "unproofread")):
        s = drift[(drift.model == "rect") & (drift.rung == "R0ipsi") & (drift.condition == cond)]
        s = s.set_index("perturbation")["drift"]
        ok = [k for k in draw["R0ipsi"] if k in s.index and pd.notna(s[k]) and pd.notna(dd.get(k))]
        if not ok:
            continue
        x, y = dd.loc[ok].to_numpy(float), s.loc[ok].to_numpy(float)
        diff = np.median(x) - np.median(y)
        n = len(ok)
        bd = np.empty(B_BOOT)
        for b in range(B_BOOT):
            i = rng.integers(0, n, n)
            bd[b] = np.median(x[i]) - np.median(y[i])
        pred[name] = {"condition": cond, "n_targets": n, "D_R0ipsi": float(np.median(x)),
                      "median_drift": float(np.median(y)), "diff": float(diff),
                      "ci95": [float(v) for v in np.quantile(bd, [0.025, 0.975])],
                      "lower_bound_one_sided_95": float(np.quantile(bd, 0.05)),
                      "p_boot_one_sided": float((1 + (bd <= 0).sum()) / (B_BOOT + 1))}
    ps = sorted(pred, key=lambda k: pred[k]["p_boot_one_sided"])
    for j, k in enumerate(ps):   # Holm over the two predictions
        adj = max(min(1.0, (len(ps) - i) * pred[ps[i]]["p_boot_one_sided"]) for i in range(j + 1))
        pred[k]["p_holm"] = adj
        # supported: one-sided bootstrap test of D_R0 - drift > 0 at Holm-adjusted 0.05 (the
        # unadjusted form of this is the one-sided 95% lower bound > 0)
        pred[k]["supported"] = bool(adj < 0.05)
    (OUT / "predictions.json").write_text(json.dumps(pred, indent=1))
    print("wrote", OUT, flush=True)


if __name__ == "__main__":
    main()
