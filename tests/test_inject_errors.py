import sys
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, "src")
import inject_errors as ie
from linear_model import build_operator, propagate


def toy_base(n=400, deg=60, seed=0, n_np=4):
    """Random graph with heavy-tailed weights and n_np neuropils; neurons assigned to a home
    neuropil, edges placed in the presynaptic neuron's home neuropil."""
    rng = np.random.default_rng(seed)
    pre = np.repeat(np.arange(n), deg)
    post = rng.integers(0, n, size=len(pre))
    k = pre != post
    pre, post = pre[k], post[k]
    w = rng.geometric(0.25, size=len(pre)).astype(np.int64)
    pre, post, w = ie.aggregate(pre, post, w, n)
    ids = np.arange(n, dtype=np.int64) + 10**6
    home = np.arange(n) % n_np
    sign = np.where(rng.random(n) < 0.35, -1.0, 1.0)
    sign[:5] = 0
    conf = rng.uniform(0.3, 1.0, n)
    conf[5:8] = np.nan
    princeton = pd.DataFrame({"pre_root_id": ids[pre], "post_root_id": ids[post],
                              "neuropil": np.array([f"NP{h}" for h in home])[pre], "syn_count": w})
    return ie.build_base(ids, pre, post, w, sign, conf, princeton)


@pytest.fixture(scope="module")
def base():
    return toy_base()


def test_grid_caps_and_unproofread():
    g = {lab: (rates, capped) for lab, _, _, _, rates, capped in ie.grid()}
    assert "spur_m25" not in g and "spur_m50" not in g and "spur_m10" in g   # 5% x m <= 1.0
    assert "flip_m5" in g and "flip_m10" not in g                           # 6% x m <= 0.5
    assert "miss_m50" in g and abs(g["miss_m50"][0]["miss"] - 0.4) < 1e-12
    r, c = g["combined_m50"]
    assert set(c) == {"spur", "flip"} and r["spur"] == 1.0 and r["flip"] == 0.5
    r, _ = g["unproofread"]
    assert abs(r["miss"] + r["split"] - (1 - 0.799)) < 1e-12 and abs(r["merge"] - 0.201) < 1e-12
    assert r["spur"] == 0.05 and r["flip"] == 0.06
    r, _ = g["lopsided_m10"]
    assert abs((r["miss"] + r["split"]) / r["merge"] - 3) < 1e-12
    assert abs(r["miss"] + r["split"] + r["merge"] - 0.16) < 1e-12


def test_reproducible_and_seed_dependent(base):
    a = ie.inject(base, {"miss": 0.1, "split": 0.1}, k=0)
    b = ie.inject(base, {"miss": 0.1, "split": 0.1}, k=0)
    c = ie.inject(base, {"miss": 0.1, "split": 0.1}, k=1)
    assert np.array_equal(a.w, b.w) and np.array_equal(a.pre, b.pre)
    assert not (len(a.w) == len(c.w) and np.array_equal(a.w, c.w))


def test_common_random_numbers_across_conditions(base):
    # the flip component draws the same neurons whether or not other components are present
    a = ie.inject(base, {"flip": 0.1}, k=2)
    b = ie.inject(base, {"flip": 0.1, "miss": 0.05, "spur": 0.05}, k=2)
    assert np.array_equal(a.flipped, b.flipped)


@pytest.mark.parametrize("p", [0.008, 0.08, 0.4])
def test_miss_rate(base, p):
    out = ie.inject(base, {"miss": p}, k=0)
    a = ie.audit(base, out)
    assert abs(a["recall_loss"] - p) < max(0.15 * p, 0.002)
    assert a["foreign_frac"] == 0 and a["n_fragments"] == 0
    assert out.w.sum() < base.w.sum()


@pytest.mark.parametrize("p", [0.04, 0.2])
def test_split_rate_and_fragments(base, p):
    out = ie.inject(base, {"split": p}, k=0)
    a = ie.audit(base, out)
    assert abs(a["recall_loss"] - p) < 0.2 * p
    assert out.w.sum() == base.w.sum()                    # synapses moved, not removed
    assert out.n_total == base.n + len(out.frag_parent)
    assert (out.pre >= base.n).any() and (out.post >= base.n).any()
    # fragment sign = parent sign
    assert np.array_equal(out.sign[base.n:], base.sign[out.frag_parent])
    # each neuron has at most one fragment
    assert len(np.unique(out.frag_parent)) == len(out.frag_parent)


@pytest.mark.parametrize("p", [0.04, 0.2, 0.4])
def test_merge_foreign_fraction(base, p):
    out = ie.inject(base, {"merge": p}, k=0)
    n = base.n
    # receivers' own synapses are untouched; added synapses ~ S*p/(1-p) per receiver
    added = out.w.sum() - base.w.sum()
    S = base.syn_per_node()
    expect = (S * p / (1 - p)).sum()
    # every copied edge is attached to its receiver, so added synapses = summed need
    assert abs(added - expect) < 0.1 * expect
    a = ie.audit(base, out)
    assert a["recall_loss"] == 0
    # receiver side ~p, partner side ~p as well -> all-ends foreign fraction ~ 2p/(1+p)
    assert abs(a["foreign_frac"] - 2 * p / (1 + p)) < 0.15 * p
    assert a["n_fragments"] == 0 and np.all(out.pre != out.post)


def test_merge_donor_shares_top_neuropil(base):
    rng = np.random.default_rng(5)
    recv = np.arange(base.n)
    d = ie._pick_donors(recv, base.top_np, rng, base.n)
    assert np.all(d != recv)
    assert np.all(base.top_np[d] == base.top_np[recv])


@pytest.mark.parametrize("p", [0.05, 0.5, 1.0])
def test_spurious_count_and_neuropil(base, p):
    out = ie.inject(base, {"spur": p}, k=0)
    added = out.w.sum() - base.w.sum()
    assert added == round(p * base.w.sum())
    assert np.all(out.pre != out.post)


def test_spurious_respects_neuropil_profiles():
    # two neuropils, disjoint neuron sets: spurious edges never cross between them
    n = 40
    ids = np.arange(n) + 1000
    pre = np.r_[np.arange(0, 19), np.arange(20, 39)]; post = pre + 1
    w = np.full(len(pre), 5, np.int64)
    p = pd.DataFrame({"pre_root_id": ids[pre], "post_root_id": ids[post],
                      "neuropil": np.where(pre < 20, "A", "B"), "syn_count": w})
    b = ie.build_base(ids, pre, post, w, np.ones(n), np.full(n, .9), p)
    out = ie.inject(b, {"spur": 1.0}, k=3)
    assert np.all((out.pre < 20) == (out.post < 20))


@pytest.mark.parametrize("p", [0.06, 0.3])
def test_flip_fraction_eligibility_and_weighting(base, p):
    out = ie.inject(base, {"flip": p}, k=0)
    elig = base.flip_weight > 0
    assert len(out.flipped) == round(p * elig.sum())
    assert np.all(base.sign[out.flipped] != 0)            # unsigned neurons never flipped
    assert np.array_equal(out.sign[out.flipped], -base.sign[out.flipped])
    other = np.setdiff1d(np.arange(base.n), out.flipped)
    assert np.array_equal(out.sign[other], base.sign[other])
    assert np.array_equal(out.w, base.w)
    # weighting: flipped neurons have lower confidence (higher 1 - conf) on average
    assert base.flip_weight[out.flipped].mean() > base.flip_weight[elig].mean()


def test_missing_conf_gets_median_weight(base):
    fw = base.flip_weight
    elig = base.sign != 0
    assert np.all(fw[5:8] > 0)          # nodes 5-7: signed, NaN conf -> eligible


def test_combined_hits_f1_target(base):
    m = 10
    rates, _ = ie.condition_rates("combined", m)
    rates = {k: v for k, v in rates.items() if k in ("miss", "split", "merge")}
    out = ie.inject(base, rates, k=0)
    a = ie.audit(base, out)
    assert abs(a["recall_loss"] - 0.008 * m) < 0.15 * 0.008 * m


def test_seed_fix_columns_restore_seed_sign(base):
    out = ie.inject(base, {"flip": 0.5, "split": 0.3}, k=0)
    A = build_operator(out.pre, out.post, out.w.astype(float), out.sign, out.n_total)
    A0_sign = out.sign.copy()
    seed = out.flipped[:3]
    fix = ie.seed_fix_columns(out, seed)
    assert set(seed) <= set(fix)
    A0_sign[fix] *= -1
    A0 = build_operator(out.pre, out.post, out.w.astype(float), A0_sign, out.n_total)
    for rect in (True, False):
        R_fix = propagate(A, [seed, out.flipped[3:4]], 0.5, 5, rectify=rect,
                          fix_cols=[fix, ie.seed_fix_columns(out, out.flipped[3:4])])
        R_ref = propagate(A0, [seed], 0.5, 5, rectify=rect)
        assert np.allclose(R_fix[:, 0], R_ref[:, 0])
    # the first seed's restored sign matches the original graph's sign
    assert np.all(A0_sign[seed] == base.sign[seed])


def test_linear_operator_matches_build_operator(base):
    out = ie.inject(base, {"miss": 0.05, "split": 0.1, "merge": 0.1, "spur": 0.2, "flip": 0.1}, k=4)
    A1 = ie.linear_operator(out)
    A2 = build_operator(out.pre, out.post, out.w.astype(float), out.sign, out.n_total)
    assert abs(A1 - A2).max() < 1e-12
