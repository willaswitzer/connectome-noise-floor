import sys, numpy as np, pandas as pd
sys.path.insert(0, "src")
from analyze_linear import pair_rho, derangement, r0_reindex, analyze


def test_pair_rho_excludes_seed_and_double_zeros():
    a = np.array([9., 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0])
    b = np.array([0., 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 5])
    ex = np.zeros(13, bool); ex[0] = True
    r, n, why = pair_rho(a, b, ex)
    assert n == 11 and why == "" and r < 1  # row 1 (0,0) dropped, row 12 (0,5) kept


def test_pair_rho_min_rows():
    r, n, why = pair_rho(np.arange(5.), np.arange(5.), np.zeros(5, bool))
    assert np.isnan(r) and why == "lt10_rows"


def test_derangement():
    rng = np.random.default_rng(0)
    for _ in range(50):
        p = derangement(7, rng)
        assert sorted(p) == list(range(7)) and not (p == np.arange(7)).any()


def test_r0_reindex():
    df = pd.DataFrame({"x": [1, 2]}, index=["5|left", "5|right"])
    assert list(r0_reindex(df, "left").index) == ["5|ipsi", "5|contra"]
    assert list(r0_reindex(df, "right").index) == ["5|contra", "5|ipsi"]


def test_matched_beats_mismatched_on_synthetic():
    rng = np.random.default_rng(1)
    base = rng.gamma(1.0, size=(200, 12))
    A = pd.DataFrame(base, index=[str(i) for i in range(200)], columns=[f"t:{i}" for i in range(12)])
    B = A + rng.normal(0, 0.1, A.shape)
    rho, h1 = analyze(A, B, {k: set() for k in A.columns}, list(A.columns), 200, rng)
    assert h1["observed_median_rho"] > 0.9 and h1["p_one_sided"] < 0.01
    assert abs(h1["label_shuffle_null_median_rho"]) < 0.1


def test_constant_vector_reason():
    a = np.zeros(20); b = np.arange(20.)
    r, n, why = pair_rho(a, b, np.zeros(20, bool))
    assert np.isnan(r) and why == "constant_vector"


def test_side_rows():
    from analyze_linear import side_rows
    df = pd.DataFrame({"x": [1, 2, 3]}, index=["5|left", "5|right", "7|left"])
    assert list(side_rows(df, "left").index) == ["5", "7"]
