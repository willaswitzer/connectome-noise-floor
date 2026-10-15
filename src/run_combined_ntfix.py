"""Exploratory (2026-10-15, after Arm-2 results were seen): the combined condition with the
NT-flip component at the confusion-matrix sign-error rate instead of the A8.1 6% anchor.

Identical to the A8.1 combined condition (segmentation 0.8%*m, spurious 5%*m, same injector,
seeds, rectified model with seed-exempt flips) except flip = FLIP_ANCHOR * m, where FLIP_ANCHOR
is the FAFB neuron-level cross-sign rate weighted by v783 tier 1-3 composition, expressed as a
fraction of all signed nodes (src/sign_error_rate.py; results/arm2/linear/ntcm_summary.csv).
Flips are still drawn by the A8.1 rule (weight 1 - top_nt_conf over all signed nodes).
Conditions: combinedntfix_m{0.5,1,2,3,5}, and unproofreadntfix (the A8.1 unproofread condition
with the same flip change). Output layout as run_arm2_linear.py.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_arm2_linear as r

FLIP_ANCHOR = 0.011
MS = [0.5, 1, 2, 3, 5]


def conditions():
    out = {}
    for m in MS:
        rates = {"miss": 0.004 * m, "split": 0.004 * m, "merge": 0.008 * m,
                 "spur": 0.05 * m, "flip": FLIP_ANCHOR * m}
        out[f"combinedntfix_m{m:g}"] = ("combined", m, None, rates, [])
    mu = r.ie.M_UNPROOFREAD
    out["unproofreadntfix"] = ("unproofread", mu, None,
                               {"miss": 0.004 * mu, "split": 0.004 * mu, "merge": 0.008 * mu,
                                "spur": 0.05, "flip": FLIP_ANCHOR}, [])
    return out


if __name__ == "__main__":
    r.conditions = conditions
    for label in conditions():
        r.run_condition(label)
