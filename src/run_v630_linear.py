"""Arm-2 real-error condition (PREREG A8.6; implementation note 2026-10-10): FlyWire v630 as a
"perturbed v783". Writes FlyWire-side vectors only, in the same layout as run_arm2_linear.py:
  data/derived/arm2/linear/<rect|norect>/<v630|v630_shiusign>/k0/{R1,R2,half_left,half_right}.parquet
Also writes data/derived/arm2/v630_graph_summary.json (structure only) and the per-neuron
loss/gain table data/derived/arm2/v630_lossgain.csv.

Usage: python src/run_v630_linear.py"""
import json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import inject_errors as ie
from linear_model import build_operator, propagate, unit_means
from run_arm2_linear import readout_setup, pad, ALPHA, K, BATCH, OUT

V630 = "Drosophila_brain_model/2023_03_23_connectivity_630_final.parquet"
V783 = "data/flywire_783_buhmann/Connectivity_783.parquet"


def read_shiu(path):
    """Shiu-format connectivity -> (pre_id, post_id, weight), one row per pair, as A2.7 builds
    flywire.parquet (Connectivity_783 is already one row per pair; summed defensively)."""
    e = pd.read_parquet(path, columns=["Presynaptic_ID", "Postsynaptic_ID", "Connectivity"])
    e = e.groupby(["Presynaptic_ID", "Postsynaptic_ID"], as_index=False)["Connectivity"].sum()
    return e["Presynaptic_ID"].to_numpy(), e["Postsynaptic_ID"].to_numpy(), e["Connectivity"].to_numpy(np.int64)


def main():
    t0 = time.time()
    # ---- check: the same reader reproduces the Amendment-7 v783 totals
    p, q, w = read_shiu(V783)
    check = {"v783_edges": int(len(w)), "v783_synapses": int(w.sum())}
    assert check == {"v783_edges": 15091983, "v783_synapses": 54492922}, check
    del p, q, w

    base = ie.load_flywire_base()
    ids, jobs = readout_setup()
    assert np.array_equal(ids, base.ids)
    n = base.n

    # ---- mapping
    m = pd.read_csv("data/derived/v630_map.csv")
    fan = m["root_id_630"].value_counts()
    merged = set(fan[fan > 1].index)
    one = m[~m["root_id_630"].isin(merged)]
    one = one[one["root_id_783"].isin(set(ids))]
    r630_to_idx = pd.Series(np.searchsorted(ids, one["root_id_783"].to_numpy()),
                            index=one["root_id_630"].to_numpy())

    pre630, post630, w630 = read_shiu(V630)
    exc = (pd.read_parquet(V630, columns=["Presynaptic_ID", "Excitatory"])
           .drop_duplicates("Presynaptic_ID").set_index("Presynaptic_ID")["Excitatory"])
    all630 = np.unique(np.concatenate([pre630, post630]))
    extra = all630[~np.isin(all630, r630_to_idx.index.to_numpy())]
    extra_idx = pd.Series(n + np.arange(len(extra)), index=extra)
    idx = pd.concat([r630_to_idx, extra_idx])
    pre = idx.reindex(pre630).to_numpy(np.int64)
    post = idx.reindex(post630).to_numpy(np.int64)
    n_total = n + len(extra)

    # ---- signs
    sign = np.concatenate([base.sign, np.zeros(len(extra))])
    m_sign = pd.Series(base.sign, index=ids)
    mm = m[m["root_id_630"].isin(merged) & m["root_id_783"].isin(set(ids))].copy()
    mm["s"] = mm["root_id_783"].map(m_sign)
    agree = mm.groupby("root_id_630")["s"].agg(lambda s: s.iloc[0] if s.nunique() == 1 else 0.0)
    is_merge = np.isin(extra, list(merged))
    sign_extra = np.where(is_merge, pd.Series(extra).map(agree).fillna(0).to_numpy(), 0.0)
    sign[n:] = sign_extra
    sign_shiu = sign.copy()
    unm = ~is_merge
    sign_shiu[n:][unm] = pd.Series(extra[unm]).map(exc).fillna(0).to_numpy()

    summ = {**check, "v630_edges": int(len(w630)), "v630_synapses": int(w630.sum()),
            "v630_nodes": int(len(all630)), "mapped_one_to_one_nodes": int(len(r630_to_idx)),
            "mapped_in_v630_graph": int(np.isin(r630_to_idx.index, all630).sum()),
            "extra_nodes": int(len(extra)), "extra_merges": int(is_merge.sum()),
            "extra_unmapped": int(unm.sum()),
            "v783_nodes_absent_in_v630": int(n - np.isin(r630_to_idx.index, all630).sum())}
    Path("data/derived/arm2").mkdir(parents=True, exist_ok=True)
    Path("data/derived/arm2/v630_graph_summary.json").write_text(json.dumps(summ, indent=1))
    print(summ, flush=True)

    # ---- measured per-neuron loss / gain (one-to-one neurons present in both graphs)
    both = np.zeros(n, bool)
    both[r630_to_idx[np.isin(r630_to_idx.index, all630)].to_numpy()] = True
    k0 = base.pre * n_total + base.post
    ok1 = (pre < n) & (post < n)
    k1 = pre[ok1] * n_total + post[ok1]
    s1 = pd.Series(w630[ok1], index=k1).groupby(level=0).sum()
    cur = s1.reindex(k0).fillna(0).to_numpy()
    kept = np.minimum(base.w, cur)
    ends = lambda a, b, ww, nn: np.bincount(a, weights=ww, minlength=nn) + np.bincount(b, weights=ww, minlength=nn)
    S0 = ends(base.pre, base.post, base.w.astype(float), n)
    Kp = ends(base.pre, base.post, kept, n)
    S1 = ends(pre[ok1], post[ok1], w630[ok1].astype(float), n) + \
        np.bincount(pre[(pre < n) & (post >= n)], weights=w630[(pre < n) & (post >= n)], minlength=n) + \
        np.bincount(post[(post < n) & (pre >= n)], weights=w630[(post < n) & (pre >= n)], minlength=n)
    lg = pd.DataFrame({"root_id": ids, "syn_v783": S0, "syn_v630": S1, "matched": Kp, "in_both": both})
    lg.to_csv("data/derived/arm2/v630_lossgain.csv", index=False)

    # ---- propagate
    for cond, sg in (("v630", sign), ("v630_shiusign", sign_shiu)):
        A = build_operator(pre, post, w630.astype(float), sg, n_total)
        for mo in ("rect", "norect"):
            d = OUT / mo / cond / "k0"
            d.mkdir(parents=True, exist_ok=True)
            for name, (keys, seeds, labels) in jobs.items():
                lab = pad(labels, n_total)
                cols = [unit_means(propagate(A, seeds[i:i + BATCH], ALPHA, K, rectify=mo == "rect"), lab)
                        for i in range(0, len(keys), BATCH)]
                df = pd.concat(cols, axis=1)
                df.columns = keys
                df[[kk for kk, s in zip(keys, seeds) if len(s) == 0]] = np.nan
                df.index = df.index.astype(str)
                df.astype(np.float32).to_parquet(d / f"{name}.parquet")
            (d / "meta.json").write_text(json.dumps({"condition": cond, "model": mo, "rectify": mo == "rect",
                                                     "summary": summ}, indent=1))
            print(f"{cond} {mo}: done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
