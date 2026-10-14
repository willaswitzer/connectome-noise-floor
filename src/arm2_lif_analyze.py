"""Arm-2 LIF drift and break-even (PREREG A8.3-A8.4), case-study set, FlyWire.

Unit response = mean rate over the unit's graph nodes and the 10 trials (Amendment 7 units;
fragment nodes have no unit). Drift for perturbation p on rung R = 1 - rho(perturbed,
unperturbed FlyWire) on R's A4.1 rows, where "unperturbed" is the first 10 trials of the
Amendment-7 FlyWire run (A8.4, like for like) and the row set is limited to rows present in R's
other dataset (R1: male CNS units; R2: hemibrain R2 units; R0: the opposite-side FlyWire vector).
Seed rows removed; double-zero rows removed; NA below 10 rows; an abolished response (perturbed
all zero, unperturbed not) gets drift 1. R0 drift = mean of the left- and right-stim drifts.
Drift is averaged over injection seeds. Noise floor = the same drift for the 3 unperturbed reruns.
D_R = median over the same perturbations of 1 - rho_10 (Amendment-7 first-10-trial rho).
Break-even (A8.4, LIF): rungs R0 and R1; primary = the 9 seed-specific perturbations; the 5 PN
seeds are a separate shared-state row; R2 descriptive only.

Outputs: results/arm2/lif/{drift.csv, noise_floor.csv, curves.csv, breakeven.csv}"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from analyze_linear import pair_rho
import arm2_analyze as a

L = Path("data/derived/lif/flywire")
A7 = Path("data/derived/lif/flywire/primary")          # Amendment-7 FlyWire runs (30 trials)
RUNS = Path("data/derived/arm2/lif")
OUT = Path("results/arm2/lif")
PN = ["DA1_uPN", "DM1_uPN", "DC1_uPN", "VA1v_uPN", "VP2_hot_PN"]
TRIALS = 10

case = pd.read_csv("results/prereg/case_study.csv")
fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
r2map = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]
seed_units = {r["name"]: set(fwu.loc[fwu["cell_type"].isin(r["flywire_types"].split("|")), "unit_id"].astype(int))
              for _, r in case.iterrows()}
SPECIFIC = [p for p in seed_units if p not in PN]
meta = json.loads((L / "seeds.json").read_text())
lab = pd.read_parquet(L / "labels.parquet")
male_units = set(pd.read_parquet("data/derived/linear/a0.5_K5/R1_male_cns.parquet").index.astype(str))
hb_r2 = set(pd.read_parquet("data/derived/linear/a0.5_K5/R2_hemibrain.parquet").index.astype(str))


def lookup(p, rung, stim):
    for x in meta["index"]:
        if x["perturbation"] == p and x["rung"] == rung and x["stim"] == stim:
            return x["run"]
    return None


def vec(path, scheme, stim_side=None):
    r = pd.read_parquet(path)
    r = r[r["trial"] < TRIALS]
    rate = lab["node_id"].map(r.groupby("node_id")["n_spikes"].sum() / TRIALS).fillna(0.0)
    key = lab[scheme]
    ok = key.notna()
    k = key[ok] if scheme == "r0" else key[ok].astype("int64").astype(str)
    v = rate[ok].groupby(k.to_numpy()).mean()
    if scheme == "r0":
        u, s = v.index.str.split("|").str[0], v.index.str.split("|").str[1]
        v.index = u + "|" + np.where(s == stim_side, "ipsi", "contra")
    return v


def drift(p_vec, u_vec, rows, seeds):
    idx = u_vec.index.intersection(rows)
    pp = p_vec.reindex(idx).fillna(0.0).to_numpy()
    uu = u_vec.loc[idx].to_numpy()
    ex = idx.isin(list(seeds))
    keep = ~ex & ~((pp == 0) & (uu == 0))
    if keep.sum() >= 10 and np.all(pp[keep] == 0) and np.any(uu[keep] != 0):
        return 1.0
    r, n, why = pair_rho(pp, uu, ex)
    return 1 - r if not np.isnan(r) else np.nan


def drift_for_dir(d):
    """Drift per (perturbation, rung) for one run directory (condition/rep)."""
    out = []
    for p in seed_units:
        # R1 and R2 use the both-sides runs
        for rung, scheme, rows, seeds in (
                ("R1", "unit", male_units, {str(u) for u in seed_units[p]}),
                ("R2", "r2", hb_r2, {str(int(r2map[u])) for u in seed_units[p]})):
            k = lookup(p, rung, "both")
            if k is None or not (d / f"{k}.parquet").exists():
                continue
            out.append({"perturbation": p, "rung": rung,
                        "drift": drift(vec(d / f"{k}.parquet", scheme), vec(A7 / f"{k}.parquet", scheme), rows, seeds)})
        kl, kr = lookup(p, "R0", "left"), lookup(p, "R0", "right")
        if kl and kr and (d / f"{kl}.parquet").exists() and (d / f"{kr}.parquet").exists():
            ul, ur = vec(A7 / f"{kl}.parquet", "r0", "left"), vec(A7 / f"{kr}.parquet", "r0", "right")
            seeds = {f"{u}|{s}" for u in seed_units[p] for s in ("ipsi", "contra")}
            dl = drift(vec(d / f"{kl}.parquet", "r0", "left"), ul, ur.index, seeds)
            dr = drift(vec(d / f"{kr}.parquet", "r0", "right"), ur, ul.index, seeds)
            out.append({"perturbation": p, "rung": "R0", "drift": np.nanmean([dl, dr])})
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for cdir in sorted(RUNS.iterdir()):
        if not cdir.is_dir():
            continue
        for rdir in sorted(cdir.iterdir()):
            if rdir.is_dir() and any(rdir.glob("*.parquet")):
                t = drift_for_dir(rdir)
                t["condition"], t["rep"] = cdir.name, rdir.name
                rows.append(t)
    dr = pd.concat(rows, ignore_index=True)
    dr.to_csv(OUT / "drift_per_rep.csv", index=False)
    g = dr.groupby(["condition", "rung", "perturbation"])["drift"]
    drift_t = g.mean().reset_index().assign(n_reps=g.count().values)
    drift_t.to_csv(OUT / "drift.csv", index=False)
    drift_t[drift_t.condition == "base"].to_csv(OUT / "noise_floor.csv", index=False)

    rho = pd.read_csv("results/lif/primary/rho_case_study.csv").set_index(["perturbation", "rung"])["rho_10"]
    rng = np.random.default_rng(a.RNG_SEED)
    res, curves = [], []
    import inject_errors as ie
    for t in ie.TYPES:
        conds = [(m, f"{t}_m{m:g}") for m in (1, 3, 10) if f"{t}_m{m:g}" in set(drift_t.condition)]
        if not conds:
            continue
        for rung in ("R0", "R1", "R2"):
            for row, ps in (("seed_specific", SPECIFIC), ("pn_shared", PN)):
                piv = drift_t[(drift_t.rung == rung)].pivot(index="perturbation", columns="condition", values="drift")
                ok = [p for p in ps if p in piv.index and (p, rung) in rho.index and pd.notna(rho[(p, rung)])]
                if not ok:
                    continue
                cd = piv.reindex(index=ok, columns=[c for _, c in conds]).to_numpy(float)
                dis = (1 - rho.loc[[(p, rung) for p in ok]]).to_numpy(float)
                floor = drift_t[(drift_t.condition == "base") & (drift_t.rung == rung) & drift_t.perturbation.isin(ok)].drift
                r = a.breakeven(cd, dis, [m for m, _ in conds], rng)
                r.update(type=t, rung=rung, row=row, n=len(ok), m_grid=[m for m, _ in conds],
                         noise_floor_median=float(np.nanmedian(floor)) if len(floor) else np.nan,
                         target=("descriptive" if rung == "R2" else "break-even"),
                         abs_rate=r["m_star"] * ie.ANCHOR[t] if np.isfinite(r["m_star"]) else np.inf)
                res.append(r)
                for (m, c), v in zip(conds, r["median_drift"]):
                    curves.append({"type": t, "rung": rung, "row": row, "m": m, "median_drift": v,
                                   "D_R": r["D_R"]})
    pd.DataFrame(res).to_csv(OUT / "breakeven.csv", index=False)
    pd.DataFrame(curves).to_csv(OUT / "curves.csv", index=False)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
