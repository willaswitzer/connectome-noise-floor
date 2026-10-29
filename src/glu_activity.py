"""Exploratory (2026-10-27): network activity of the glutamate and activity-matched LIF runs.

For every run in a condition directory, compare with the run of the same key in the Amendment-7
reference (data/derived/lif/flywire/primary, first 10 trials): total spikes per trial and number
of neurons that spike at least once (summed over trials, / trials). Reports per condition the
median over runs of the ratio (condition / reference).
Output: results/arm2/lif/glu_activity_matched.csv (glu_activity.csv, 2026-10-17, is left as written)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import arm2_lif_analyze as la

CONDS = ["glu_post_f10", "glu_post_f25", "glu_pre_f25", "scaleonly_ws0.85", "glu_post_f25_ws0.85"]


def stats(path, ntr=10):
    d = pd.read_parquet(path, columns=["node_id", "trial", "n_spikes"])
    d = d[d.trial < ntr]
    nt = d.trial.nunique()
    return d.n_spikes.sum() / nt, (d.groupby("node_id").n_spikes.sum() > 0).sum()


def main():
    rows = []
    for c in CONDS:
        for rd in sorted((la.RUNS / c).iterdir()):
            for f in sorted(rd.glob("*.parquet")):
                ref = la.A7 / f.name
                if not ref.exists():
                    continue
                s, a = stats(f)
                s0, a0 = stats(ref)
                rows.append({"condition": c, "rep": rd.name, "run": f.stem,
                             "spike_ratio": s / s0 if s0 else np.nan, "active_ratio": a / a0 if a0 else np.nan})
    per = pd.DataFrame(rows)
    out = per.groupby("condition").agg(n_runs=("run", "size"), spike_ratio_median=("spike_ratio", "median"),
                                       spike_ratio_max=("spike_ratio", "max"),
                                       active_ratio_median=("active_ratio", "median"),
                                       active_ratio_max=("active_ratio", "max")).reindex(CONDS).reset_index()
    out.to_csv(la.OUT / "glu_activity_matched.csv", index=False)
    print(out.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
