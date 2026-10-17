"""Exploratory (2026-10-16): how much of the left-right noise floor (R0) is within-type sign
disagreement, i.e. transmitter-prediction noise rather than biology?

Every signed v783 neuron whose A3.4 sign differs from its unit's majority sign (ties -> +1) is set
to the majority. The graph is otherwise unchanged. The linear model (A2.1, alpha 0.5, K 5) is run
rectified and unrectified on the Arm-1 seed sets; R0 and R0ipsi agreement is then computed exactly
as analyze_linear.py (pair_rho, seed rows removed) and compared, target by target, with the same
quantity on the unharmonized base run.

Outputs: results/linear/exploratory/sign_harmonized_r0.csv (per target) and
         sign_harmonized_r0_summary.json
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import inject_errors as ie
from linear_model import propagate, unit_means
from run_arm2_linear import readout_setup, pad, ALPHA, K, BATCH
from analyze_linear import pair_rho, r0_reindex, side_rows
from arm2_analyze import seed_rows, r0_rows

OUT = Path("results/linear/exploratory")


def minority_neurons(base):
    s = pd.read_parquet("data/derived/signs_flywire.parquet")[["root_id", "sign"]]
    u = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")[["root_id", "unit_id"]]
    d = u.merge(s, on="root_id")
    d = d[d.sign != 0]
    maj = d.groupby("unit_id")["sign"].agg(lambda x: 1 if (x > 0).sum() >= (x < 0).sum() else -1)
    ids = d.loc[d.sign != d.unit_id.map(maj), "root_id"].to_numpy()
    pos = pd.Series(np.arange(base.n), index=base.ids)
    return np.sort(pos.reindex(ids).dropna().astype(np.int64).to_numpy())


def run(base, flipped, jobs, rect):
    sign = base.sign.copy()
    sign[flipped] *= -1
    parts = [tuple(x.astype(np.int32) for x in (base.pre, base.post, base.w))]
    pert = ie.Perturbed(parts, sign, base.n, np.zeros(0, np.int64), flipped, {}, [], 0)
    A = ie.linear_operator(pert, consume=True)
    out = {}
    for name in ("half_left", "half_right"):
        keys, seeds, labs = jobs[name]
        lab = pad(labs, pert.n_total)
        cols = [unit_means(propagate(A, seeds[i:i + BATCH], ALPHA, K, rectify=rect), lab)
                for i in range(0, len(keys), BATCH)]
        df = pd.concat(cols, axis=1)
        df.columns = keys
        df[[k for k, s in zip(keys, seeds) if len(s) == 0]] = np.nan
        df.index = df.index.astype(str)
        out[name] = df.astype(np.float64)
    return out


def rhos(v, seeds, draw):
    hl, hr = v["half_left"], v["half_right"]
    forms = {"R0": (r0_reindex(r0_rows(hl), "left"), r0_reindex(r0_rows(hr), "right")),
             "R0ipsi": (side_rows(hl, "left"), side_rows(hr, "right"))}
    rows = []
    for rung, (A, B) in forms.items():
        idx = A.index.intersection(B.index)
        for k in draw[rung]:
            if k not in A.columns:
                continue
            r, n, why = pair_rho(A.loc[idx, k].to_numpy(), B.loc[idx, k].to_numpy(),
                                 idx.isin(list(seeds[rung][k])))
            rows.append({"rung": rung, "perturbation": k, "rho": r})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = ie.load_flywire_base()
    ids, jobs = readout_setup()
    flipped = minority_neurons(base)
    seeds, draw = seed_rows()
    res, summ = [], {"n_harmonized": int(len(flipped)), "n_signed": int((base.flip_weight > 0).sum())}
    for model, rect in (("rect", True), ("norect", False)):
        r0 = rhos(run(base, np.zeros(0, np.int64), jobs, rect), seeds, draw).rename(columns={"rho": "rho_base"})
        r1 = rhos(run(base, flipped, jobs, rect), seeds, draw).rename(columns={"rho": "rho_harmonized"})
        m = r0.merge(r1, on=["rung", "perturbation"]).assign(model=model)
        res.append(m)
        for rung in ("R0", "R0ipsi"):
            x = m[m.rung == rung].dropna()
            summ[f"{model}_{rung}"] = {
                "n": int(len(x)), "D_base": float(np.median(1 - x.rho_base)),
                "D_harmonized": float(np.median(1 - x.rho_harmonized)),
                "median_change_in_rho": float(np.median(x.rho_harmonized - x.rho_base)),
                "max_abs_change_in_rho": float(np.max(np.abs(x.rho_harmonized - x.rho_base)))}
    pd.concat(res).to_csv(OUT / "sign_harmonized_r0.csv", index=False)
    (OUT / "sign_harmonized_r0_summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
