"""Exploratory (not pre-registered; 2026-10-15, after all Arm-2 results were seen):
transmitter-sign errors at rates taken from the Eckstein et al. 2024 Fig. 2A confusion matrices.

Why. The A8.1 NT anchor (6% of signed neurons flipped E<->I) came from the 94% neuron-level
label accuracy. Label errors are not sign errors: the classifier's main confusion (Glu -> GABA)
keeps the sign. Tier-0 neurons (A3.4: type-level known transmitter, 63% of signed v783 nodes)
do not take the classifier's prediction at all. This script injects sign flips at the rates the
confusion matrices imply, only where a prediction is used.

Error model.
  Eligible: v783 graph nodes with sign != 0 and A3.4 tier 1-3 (sign from a prediction).
  p_flip(neuron) = cross-sign mass of the confusion-matrix row of its assigned transmitter
  (A3.4 signs). Neuron-level matrices; their 5-HT and OA rows have n = 1 neuron, so those two
  rows are taken from the synapse-level matrix of the same dataset. The row of the PREDICTED
  transmitter is used as the flip probability (no Bayes inversion; the prior is unknown).
  Matrices (data/anchors/eckstein_fig2a.csv):
    fafb   FAFB-Catmaid (same imagery as FlyWire): primary
    hemi   hemibrain: pessimistic sensitivity
  Draws:
    neuron  independent Bernoulli per eligible neuron
    type    one Bernoulli per cell type (signs_flywire `type`), with the rate of the type's
            majority eligible transmitter; if it fires, every eligible member flips. Neurons
            with no type are drawn individually. Rationale: Eckstein et al. 2024 Fig. 4D-F and
            S2: errors are shared within types and across datasets.
  5 draws per condition, default_rng(20261015 + k); one child stream per condition.

Models: norect (NT primary, compared with norect D_R) and rect_ex (seed-exempt flips,
compared with rect D_R), as A8.2. Readouts, seeds and output layout as run_arm2_linear.py:
data/derived/arm2/linear/<model>/<condition>/k<k>/.

Usage: python src/run_nt_realistic.py [--only ntcm_fafb_neuron,...]
"""
import argparse, json, os, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import inject_errors as ie
from linear_model import propagate, unit_means
from run_arm2_linear import readout_setup, pad, OUT, ALPHA, K, BATCH

ANCHORS = Path("data/anchors/eckstein_fig2a.csv")
SIGNS = Path("data/derived/signs_flywire.parquet")
NT_MAP = {"acetylcholine": "ACh", "glutamate": "Glu", "gaba": "GABA", "dopamine": "DA",
          "serotonin": "5HT", "octopamine": "OA"}
NTS = ["ACh", "Glu", "GABA", "DA", "5HT", "OA"]
SIGN = {"ACh": 1, "DA": 1, "5HT": 1, "OA": 1, "Glu": -1, "GABA": -1}
DATASET = {"fafb": "FAFB-catmaid", "hemi": "hemibrain"}
CONDITIONS = [f"ntcm_{d}_{draw}" for d in ("fafb", "hemi") for draw in ("neuron", "type")]
SEED = 20261015
KS = range(5)
MODELS = ["norect", "rect_ex"]


def cross_sign_rates(dataset):
    """Per assigned transmitter: P(predicted sign != true sign), neuron level, with the 5-HT and
    OA rows from the synapse level (neuron-level n = 1)."""
    m = pd.read_csv(ANCHORS, comment="#")
    m = m[m["dataset"] == DATASET[dataset]]
    out = {}
    for nt in NTS:
        lvl = "synapse" if nt in ("5HT", "OA") else "neuron"
        row = m[(m["level"] == lvl) & (m["true_nt"] == nt)].iloc[0]
        tot = sum(row[p] for p in NTS)
        out[nt] = float(sum(row[p] for p in NTS if SIGN[p] != SIGN[nt]) / tot)
    return out


def eligible(base):
    s = pd.read_parquet(SIGNS).set_index("root_id")
    s = s.reindex(base.ids)
    ok = (s["sign"].fillna(0).to_numpy() != 0) & s["tier"].isin([1, 2, 3]).to_numpy()
    nt = s["nt"].map(NT_MAP).to_numpy(dtype=object)
    ok &= pd.notna(nt)
    # the graph's sign must agree with the table's (same A3.4 assignment)
    assert np.all(np.sign(base.sign[ok]) == s["sign"].to_numpy()[ok]), "sign table mismatch"
    return np.flatnonzero(ok), nt, s["type"].to_numpy(dtype=object)


def draw_flips(label, base, k):
    _, dataset, how = label.split("_")
    rates = cross_sign_rates(dataset)
    idx, nt, typ = eligible(base)
    rng = np.random.default_rng(SEED + k).spawn(len(CONDITIONS))[CONDITIONS.index(label)]
    p = np.array([rates[nt[i]] for i in idx])
    if how == "neuron":
        flip = idx[rng.random(len(idx)) < p]
    else:
        df = pd.DataFrame({"i": idx, "type": typ[idx], "nt": nt[idx], "p": p})
        typed = df[pd.notna(df["type"])]
        untyped = df[pd.isna(df["type"])]
        grp = typed.groupby("type", sort=True)
        maj = grp["nt"].agg(lambda x: x.value_counts().sort_index().idxmax())
        hit = rng.random(len(maj)) < maj.map(rates).to_numpy()
        fired = set(maj.index[hit])
        flip_t = typed.loc[typed["type"].isin(fired), "i"].to_numpy()
        flip_u = untyped["i"].to_numpy()[rng.random(len(untyped)) < untyped["p"].to_numpy()]
        flip = np.concatenate([flip_t, flip_u])
    return np.sort(flip.astype(np.int64)), rates, len(idx)


def perturbed(base, flipped, k):
    sign = base.sign.copy()
    sign[flipped] *= -1
    parts = [tuple(x.astype(np.int32) for x in (base.pre, base.post, base.w))]
    # no structural change: every original synapse end retained, none added (for ie.audit)
    retained = ie._ends(base.pre, base.post, base.w, base.n)
    return ie.Perturbed(parts, sign, base.n, np.zeros(0, np.int64), flipped, {}, [], k,
                        retained, np.zeros(base.n))


def run(labels):
    t0 = time.time()
    base = ie.load_flywire_base()
    ids, jobs = readout_setup()
    assert np.array_equal(ids, base.ids)
    n_signed = int((base.flip_weight > 0).sum())
    print(f"setup {time.time() - t0:.0f}s", flush=True)
    for label in labels:
        for k in KS:
            todo = [mo for mo in MODELS if not (OUT / mo / label / f"k{k}" / "meta.json").exists()]
            if not todo:
                continue
            flipped, rates, n_elig = draw_flips(label, base, k)
            pert = perturbed(base, flipped, k)
            fixes = {name: [ie.seed_fix_columns(pert, s) for s in seeds]
                     for name, (_, seeds, _) in jobs.items()}
            A = ie.linear_operator(pert, consume=True)
            for mo in todo:
                d = OUT / mo / label / f"k{k}"
                d.mkdir(parents=True, exist_ok=True)
                rect, exempt = mo != "norect", mo == "rect_ex"
                t2 = time.time()
                for name, (keys, seeds, labs) in jobs.items():
                    lab = pad(labs, pert.n_total)
                    cols = []
                    for i in range(0, len(keys), BATCH):
                        fc = fixes[name][i:i + BATCH] if exempt else None
                        R = propagate(A, seeds[i:i + BATCH], ALPHA, K, rectify=rect, fix_cols=fc)
                        cols.append(unit_means(R, lab))
                    df = pd.concat(cols, axis=1)
                    df.columns = keys
                    df[[kk for kk, s in zip(keys, seeds) if len(s) == 0]] = np.nan
                    df.index = df.index.astype(str)
                    df.astype(np.float32).to_parquet(d / f"{name}.parquet")
                meta = {"condition": label, "k": k, "model": mo, "rectify": rect,
                        "seed_exempt": exempt, "alpha": ALPHA, "K": K,
                        "cross_sign_rates": rates, "n_eligible": n_elig,
                        "n_flipped": int(len(flipped)),
                        "flip_frac_eligible": len(flipped) / n_elig,
                        "flip_frac_signed": len(flipped) / n_signed,
                        "seconds_propagate": round(time.time() - t2, 1),
                        "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
                        "git": os.environ.get("ARM2_GIT", "")}
                (d / "meta.json").write_text(json.dumps(meta, indent=1))
                print(f"{label} k{k} {mo}: {len(flipped)} flips "
                      f"({100 * len(flipped) / n_signed:.2f}% of signed), "
                      f"{time.time() - t2:.0f}s", flush=True)
            del A, pert


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    run([l for l in CONDITIONS if not a.only or l in a.only.split(",")])
