"""Arm-2 linear sweep (PREREG Amendment 8, A8.2-A8.3): FlyWire-side response vectors ONLY.

For every grid condition (inject_errors.grid()) and injection seed k = 0..4, the injected v783
graph is propagated with the A2.1 model (alpha = 0.5, K = 5) on the Arm-1 seed sets, and the
FlyWire-side readouts are saved. No drift, correlation or summary is computed here.

Models (A8.2):
  rect      rectified primary, for missed / split / merge / spurious, and for combined,
            lopsided and unproofread (their flip component seed-exempt; implementation note
            2026-10-07)
  norect    unrectified, for the flip type (primary for NT)
  rect_ex   rectified with seed-exempt flips, for the flip type (sensitivity)
The unperturbed graph is run once per model as condition "base" (k = 0), through the same code,
and checked against the Arm-1 vectors.

Readouts (same seed sets and row labels as run_linear.py):
  R1.parquet           rows unit_id, columns all perturbations ("t:", "c:")
  R2.parquet           rows r2_unit_id, columns R2 perturbations
  half_left.parquet    rows "<unit>|<side>", left-side stimulation (A5.3 forms)
  half_right.parquet   rows "<unit>|<side>", right-side stimulation
R0 (A2.2) uses the same one-sided seed sets as the half forms; its rows are the half rows of
R0-usable units, so it is derived from half_* in the analysis rather than stored twice.
Fragment nodes (false splits) belong to no unit and are never read out.
Values are stored as float32.

Output: data/derived/arm2/linear/<model>/<condition>/k<k>/  (+ meta.json written last)
Resumable: a (model, condition, k) with meta.json is skipped.

Usage:
  python src/run_arm2_linear.py                 # whole grid, one subprocess per condition
  python src/run_arm2_linear.py --only combined_m1,unproofread
  python src/run_arm2_linear.py --one <condition>   # internal: run one condition, all k"""
import argparse, json, os, subprocess, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import inject_errors as ie
from linear_model import load_graph, propagate, unit_means

OUT = Path("data/derived/arm2/linear")
KS = range(5)
ALPHA, K = 0.5, 5
BATCH = 50

# priority order: P5/P6 conditions and the baseline first
FIRST = ["base", "combined_m1", "unproofread"]


def models_for(label, kind, t):
    if label == "base":
        return ["rect", "norect"]
    if kind == "single" and t == "flip":
        return ["norect", "rect_ex"]
    return ["rect"]


def conditions():
    g = {lab: (kind, m, t, rates, capped) for lab, kind, m, t, rates, capped in ie.grid()}
    g["base"] = ("base", 0, None, {}, [])
    order = FIRST + [l for l in g if l not in FIRST]
    return {l: g[l] for l in order}


def perturbations():
    targets = pd.read_csv("results/prereg/targets_DRAFT.csv")
    fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
    case = pd.read_csv("results/prereg/case_study.csv")
    r2map = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]
    perts = {f"t:{u}": {int(u)} for u in targets["unit_id"]}
    for _, row in case.iterrows():
        us = set(fwu.loc[fwu["cell_type"].isin(row["flywire_types"].split("|")), "unit_id"].astype(int))
        perts[f"c:{row['name']}"] = us
    tierB = {f"t:{u}" for u in targets.loc[targets["tierB"], "unit_id"]}
    case_r2 = {"c:" + n for n in case["name"]} - {"c:sugar_GRN", "c:bitter_GRN", "c:JO-A", "c:JO-CE"}
    r2_perts = [k for k in perts if k in tierB or k in case_r2]
    r2sets = {k: {int(r2map[x]) for x in perts[k]} for k in r2_perts}
    return perts, r2_perts, r2sets


def readout_setup():
    """Node labels and seed sets on the v783 graph, exactly as run_linear.py builds them."""
    g = load_graph("flywire")
    perts, r2_perts, r2sets = perturbations()
    side = pd.read_csv("data/flywire_783/annotations_v3_neurons.tsv", sep="\t",
                       usecols=["root_id", "side"], low_memory=False).set_index("root_id")["side"]
    u = g.unit.to_numpy(dtype="float64", na_value=np.nan)
    r2 = g.r2.to_numpy(dtype="float64", na_value=np.nan)
    sd = pd.Series(g.ids).map(side).to_numpy()
    isin = lambda arr, s: np.isin(arr, list(s))
    keys = list(perts)
    labh = pd.Series([f"{int(a)}|{b}" if (pd.notna(a) and b in ("left", "right")) else None
                      for a, b in zip(u, sd)], dtype=object)
    jobs = {
        "R1": (keys, [np.flatnonzero(isin(u, perts[k])) for k in keys], g.unit),
        "R2": (r2_perts, [np.flatnonzero(isin(r2, r2sets[k])) for k in r2_perts], g.r2),
    }
    for s in ("left", "right"):
        jobs[f"half_{s}"] = (keys, [np.flatnonzero(isin(u, perts[k]) & (sd == s)) for k in keys], labh)
    return g.ids, jobs


def pad(labels, n_total):
    """Extend a per-node label series with NA for fragment nodes."""
    lab = pd.Series(labels).reset_index(drop=True)
    if n_total > len(lab):
        lab = pd.concat([lab, pd.Series([pd.NA] * (n_total - len(lab)), dtype=lab.dtype)],
                        ignore_index=True)
    return lab


def run_condition(label):
    kind, m, t, rates, capped = conditions()[label]
    models = models_for(label, kind, t)
    ks = [0] if label == "base" else list(KS)
    todo = [(mo, k) for mo in models for k in ks
            if not (OUT / mo / label / f"k{k}" / "meta.json").exists()]
    if not todo:
        print(f"{label}: done already", flush=True)
        return
    t0 = time.time()
    base = ie.load_flywire_base()
    ids, jobs = readout_setup()
    assert np.array_equal(ids, base.ids), "graph node order differs from injector base"
    print(f"{label}: setup {time.time() - t0:.0f}s", flush=True)
    for k in sorted({k for _, k in todo}):
        t1 = time.time()
        pert = ie.inject(base, rates, k, capped)
        aud = ie.audit(base, pert)
        nt, sign, flipped = pert.n_total, pert.sign, pert.flipped
        fixes = {name: [ie.seed_fix_columns(pert, s) for s in seeds]
                 for name, (_, seeds, _) in jobs.items()} if len(flipped) else None
        A = ie.linear_operator(pert, consume=True)
        t_inj = time.time() - t1
        for mo in [mo for mo, kk in todo if kk == k]:
            d = OUT / mo / label / f"k{k}"
            d.mkdir(parents=True, exist_ok=True)
            rect = mo != "norect"
            exempt = mo in ("rect", "rect_ex") and fixes is not None
            t2 = time.time()
            empty = {}
            for name, (keys, seeds, labels) in jobs.items():
                lab = pad(labels, nt)
                cols = []
                for i in range(0, len(keys), BATCH):
                    fc = fixes[name][i:i + BATCH] if exempt else None
                    R = propagate(A, seeds[i:i + BATCH], ALPHA, K, rectify=rect, fix_cols=fc)
                    cols.append(unit_means(R, lab))
                df = pd.concat(cols, axis=1)
                df.columns = keys
                e = [kk for kk, s in zip(keys, seeds) if len(s) == 0]
                df[e] = np.nan
                empty[name] = e
                df.index = df.index.astype(str)
                df.astype(np.float32).to_parquet(d / f"{name}.parquet")
            meta = {"condition": label, "kind": kind, "m": m, "type": t, "k": k, "model": mo,
                    "rectify": rect, "seed_exempt": bool(exempt), "rates": rates,
                    "capped": capped, "alpha": ALPHA, "K": K, "audit": aud,
                    "n_fragments": int(len(pert.frag_parent)), "n_flipped": int(len(flipped)),
                    "nnz": int(A.nnz), "empty_seeds": empty,
                    "seconds_inject": round(t_inj, 1), "seconds_propagate": round(time.time() - t2, 1),
                    "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "git": os.environ.get("ARM2_GIT") or subprocess.run(
                        ["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                        text=True).stdout.strip()}
            (d / "meta.json").write_text(json.dumps(meta, indent=1))
            print(f"{label} k{k} {mo}: inject {t_inj:.0f}s, propagate {time.time() - t2:.0f}s",
                  flush=True)
        del A, pert


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--one", default=None)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    if a.one:
        run_condition(a.one)
        return
    labels = list(conditions())
    if a.only:
        want = a.only.split(",")
        labels = [l for l in labels if l in want]
    for l in labels:
        r = subprocess.run([sys.executable, __file__, "--one", l])
        if r.returncode != 0:
            print(f"{l}: FAILED (exit {r.returncode})", flush=True)


if __name__ == "__main__":
    main()
