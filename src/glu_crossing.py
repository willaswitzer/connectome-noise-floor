"""Exploratory (2026-10-30): where glutamate's sign reaches the LIF noise floor.

The receiver-side conditions pick a fraction f of postsynaptic units (cell types; untyped neurons
singly), and unit sizes vary widely, so the share of glutamatergic connections actually made
excitatory varies between draws (e.g. 11-18% at f = 0.15). This regresses each draw's R0 drift / D_R0
(median drift over the row's perturbations / median D_R0, D_R0 = 1 - rho_10) on that share
(n_edges_flipped in the linear meta.json / all glutamatergic connections of the v783 graph),
over the 12 receiver-side LIF draws (f = 0.10, 0.15, 0.20, 0.25; k = 0-2), and reports where the
fitted line reaches 1, with a 95% bootstrap interval over draws (2,000 resamples,
default_rng(20261030)).
Outputs: results/arm2/lif/glu_crossing_per_draw.csv, glu_crossing.csv
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import arm2_lif_analyze as la
import inject_errors as ie
import run_glu_sign as g

FILES = ["glu_drift_per_rep.csv", "glu_post_f15_drift_per_rep.csv", "glu_post_f20_drift_per_rep.csv"]
LIN = Path("data/derived/arm2/linear/norect")


def main():
    base = ie.load_flywire_base()
    glu, _ = g.node_tables(base)
    n_glu = int(glu[base.pre].sum())
    rho = pd.read_csv("results/lif/primary/rho_case_study.csv").set_index(["perturbation", "rung"])["rho_10"]
    d = pd.concat([pd.read_csv(la.OUT / f) for f in FILES])
    d = d[(d.rung == "R0") & d.condition.str.startswith("glu_post")]
    rows = []
    for (c, r), grp in d.groupby(["condition", "rep"]):
        m = json.loads((LIN / c / r / "meta.json").read_text())
        for row, ps in (("seed_specific", la.SPECIFIC), ("pn_shared", la.PN)):
            g2 = grp[grp.perturbation.isin(ps)]
            D = np.median([1 - rho[(p, "R0")] for p in g2.perturbation])
            rows.append({"condition": c, "rep": r, "row": row, "glu_edges_excitatory": m["n_edges_flipped"],
                         "glu_edge_frac": m["n_edges_flipped"] / n_glu,
                         "drift_over_D": float(np.median(g2.drift) / D)})
    per = pd.DataFrame(rows)
    per.to_csv(la.OUT / "glu_crossing_per_draw.csv", index=False)
    rng = np.random.default_rng(20261030)
    out = []
    for row in ("seed_specific", "pn_shared"):
        s = per[per.row == row]
        x, y = s.glu_edge_frac.to_numpy(), s.drift_over_D.to_numpy()
        b, a = np.polyfit(x, y, 1)
        bs = []
        for _ in range(2000):
            i = rng.integers(0, len(x), len(x))
            bb, aa = np.polyfit(x[i], y[i], 1)
            bs.append((1 - aa) / bb)
        out.append({"row": row, "n_draws": len(s), "n_glu_edges": n_glu, "slope": b, "intercept": a,
                    "r": float(np.corrcoef(x, y)[0, 1]), "crossing_glu_edge_frac": (1 - a) / b,
                    "ci_lo": float(np.quantile(bs, .025)), "ci_hi": float(np.quantile(bs, .975))})
    res = pd.DataFrame(out)
    res.to_csv(la.OUT / "glu_crossing.csv", index=False)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
