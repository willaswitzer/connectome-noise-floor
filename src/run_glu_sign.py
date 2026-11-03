"""Exploratory (2026-10-16): sensitivity of the linear model to the sign of glutamate.

A3.4 makes glutamate inhibitory (-1). In the fly, glutamate inhibits through GluClalpha and
excites through ionotropic receptors, and its sign depends on the receptors of the receiving
neuron (Eckstein et al. 2024, Discussion). Two ways of making a fraction f of glutamatergic
transmission excitatory, f in {0.10, 0.25, 0.50}, 3 draws each (default_rng(20261016 + k)):

  glu_pre_f   sender side: a fraction f of glutamatergic presynaptic units (cell types; untyped
              neurons individually) become +1 onto all targets.
  glu_post_f  receiver side: a fraction f of postsynaptic units (all units, glutamatergic or not)
              receive ALL their glutamatergic input as excitatory (edge-level sign).

Glutamatergic = A3.4 transmitter 'glutamate' (any tier). Models: unrectified (norect, compared
with the unrectified D_R) and rectified (rect, rectified D_R). No seed exemption: a glutamatergic
seed gives no rectified response in the base run (A5.1), so its target is NA in both runs.
Outputs: data/derived/arm2/linear/<model>/<condition>/k<k>/ (as run_arm2_linear.py); with
--analyze: results/arm2/linear/glu_summary.csv, glu_drift.csv
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import inject_errors as ie
from linear_model import propagate, unit_means
from run_arm2_linear import readout_setup, pad, OUT, ALPHA, K, BATCH

FS = [0.10, 0.25, 0.50]
CONDS = [f"glu_{side}_f{int(f * 100)}" for side in ("pre", "post") for f in FS]
# exploratory 2026-10-30: receiver-side f = 0.15, 0.20 to locate the LIF crossing. Appended, so the
# spawn index (and every draw) of the six conditions above is unchanged.
CONDS += ["glu_post_f15", "glu_post_f20"]
SEED = 20261016
MODELS = ["norect", "rect"]


def node_tables(base):
    s = pd.read_parquet("data/derived/signs_flywire.parquet").set_index("root_id").reindex(base.ids)
    glu = (s["nt"].to_numpy(dtype=object) == "glutamate") & (base.sign < 0)
    u = pd.read_csv("data/derived/units_v2_flywire_neurons.csv").set_index("root_id")["unit_id"]
    unit = pd.Series(base.ids).map(u).to_numpy(dtype=float)
    # untyped nodes become their own group
    grp = np.where(np.isnan(unit), -(np.arange(base.n) + 1), unit)
    return glu, grp


def draw(base, cond, k, glu, grp):
    """The draw for (cond, k): (node signs after sender-side changes, receiver mask or None, info).
    Shared by the linear runner and the LIF runner (run_arm2_lif.py --glu), so both use the same draws."""
    side, f = cond.split("_")[1], int(cond.split("_f")[1]) / 100
    rng = np.random.default_rng(SEED + k).spawn(len(CONDS))[CONDS.index(cond)]
    sign = base.sign.copy()
    post_sel = None
    if side == "pre":
        g = np.unique(grp[glu])
        chosen = g[rng.random(len(g)) < f]
        flip = glu & np.isin(grp, chosen)
        sign[flip] = 1
        info = {"n_glu_neurons_excitatory": int(flip.sum())}
    else:
        g = np.unique(grp)
        chosen = g[rng.random(len(g)) < f]
        post_sel = np.isin(grp, chosen)
        info = {"n_post_nodes_glu_excitatory": int(post_sel.sum())}
    info["n_glu_neurons"] = int(glu.sum())
    return sign, post_sel, info


def operator(base, cond, k, glu, grp):
    sign, post_sel, info = draw(base, cond, k, glu, grp)
    parts = [tuple(x.astype(np.int32) for x in (base.pre, base.post, base.w))]
    pert = ie.Perturbed(parts, sign, base.n, np.zeros(0, np.int64), np.zeros(0, np.int64), {}, [], k)
    A = ie.linear_operator(pert, consume=True)
    if post_sel is not None:
        rows = np.repeat(np.arange(A.shape[0]), np.diff(A.indptr))
        m = post_sel[rows] & glu[A.indices]
        A.data[m] *= -1
        info["n_edges_flipped"] = int(m.sum())
    return A, info


def run(conds):
    base = ie.load_flywire_base()
    ids, jobs = readout_setup()
    glu, grp = node_tables(base)
    for cond in conds:
        for k in range(3):
            todo = [mo for mo in MODELS if not (OUT / mo / cond / f"k{k}" / "meta.json").exists()]
            if not todo:
                continue
            A, info = operator(base, cond, k, glu, grp)
            for mo in todo:
                d = OUT / mo / cond / f"k{k}"
                d.mkdir(parents=True, exist_ok=True)
                t = time.time()
                for name, (keys, seeds, labs) in jobs.items():
                    lab = pad(labs, base.n)
                    cols = [unit_means(propagate(A, seeds[i:i + BATCH], ALPHA, K, rectify=mo == "rect"), lab)
                            for i in range(0, len(keys), BATCH)]
                    df = pd.concat(cols, axis=1)
                    df.columns = keys
                    df[[kk for kk, s in zip(keys, seeds) if len(s) == 0]] = np.nan
                    df.index = df.index.astype(str)
                    df.astype(np.float32).to_parquet(d / f"{name}.parquet")
                (d / "meta.json").write_text(json.dumps({"condition": cond, "k": k, "model": mo, **info}, indent=1))
                print(f"{cond} k{k} {mo}: {info} {time.time() - t:.0f}s", flush=True)


def analyze():
    from arm2_analyze import seed_rows, other_side, drift_table, rung_disagreement, D
    seeds, draw = seed_rows()
    tabs = []
    for mo in MODELS:
        other, b = other_side(mo)
        for c in CONDS:
            for kd in sorted((D / mo / c).glob("k*")):
                if (kd / "meta.json").exists():
                    t = drift_table(kd, b, other, seeds)
                    t["model"], t["condition"], t["k"] = mo, c, int(kd.name[1:])
                    tabs.append(t)
    dr = pd.concat(tabs, ignore_index=True)
    drift = dr.groupby(["model", "condition", "rung", "perturbation"])["drift"].agg(
        lambda v: np.nanmean(v) if v.notna().any() else np.nan).reset_index()
    drift.to_csv("results/arm2/linear/glu_drift.csv", index=False)
    rng = np.random.default_rng(SEED)
    rows = []
    for mo in MODELS:
        dis = rung_disagreement(mo)
        for c in CONDS:
            for rung in ("R0", "R1", "R2"):
                s = drift[(drift.model == mo) & (drift.condition == c) & (drift.rung == rung)].set_index("perturbation")["drift"]
                dd = dis[dis.rung == rung].set_index("perturbation")["dis"]
                ok = [p for p in draw[rung] if p in s.index and pd.notna(s[p]) and pd.notna(dd.get(p))]
                if not ok:
                    continue
                x, y = dd.loc[ok].to_numpy(float), s.loc[ok].to_numpy(float)
                n = len(ok)
                bd = np.array([np.median(x[i]) - np.median(y[i]) for i in (rng.integers(0, n, n) for _ in range(1000))])
                rows.append({"condition": c, "model": mo, "rung": rung, "n": n, "median_drift": float(np.median(y)),
                             "D_R": float(np.median(x)), "drift_over_D": float(np.median(y) / np.median(x)),
                             "ci_lo": float(np.quantile(bd, .025)), "ci_hi": float(np.quantile(bd, .975))})
    out = pd.DataFrame(rows)
    out.to_csv("results/arm2/linear/glu_summary.csv", index=False)
    with pd.option_context("display.width", 200):
        print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyze", action="store_true")
    a = ap.parse_args()
    analyze() if a.analyze else run(CONDS)
