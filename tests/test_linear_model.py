import sys, numpy as np, pandas as pd
sys.path.insert(0, "src")
from linear_model import build_operator, propagate, unit_means


def test_hand_computed_chain():
    # 0 -> 1 (4 syn, excitatory), 2 -> 1 (4 syn, inhibitory), 1 -> 3 (2 syn), 0 -> 3 (6 syn)
    pre = np.array([0, 2, 1, 0]); post = np.array([1, 1, 3, 3]); w = np.array([4., 4., 2., 6.])
    sign = np.array([1., 1., -1., 1.])  # node 2 inhibitory
    A = build_operator(pre, post, w, sign, 4)
    # input fractions: onto 1: 0->1 = 0.5, 2->1 = -0.5; onto 3: 1->3 = 0.25, 0->3 = 0.75
    R = propagate(A, [np.array([0])], alpha=0.5, K=2)[:, 0]
    # hop1: r1 = [0, .25, 0, .375]; hop2: node3 gets .5*.25*.25 = .03125, node1 gets 0
    assert np.allclose(R, [0, 0.25, 0, 0.375 + 0.03125])


def test_rectification_and_inhibition():
    pre = np.array([0, 2]); post = np.array([1, 1]); w = np.array([1., 3.])
    sign = np.array([1., 1., -1.])
    A = build_operator(pre, post, w, sign, 3)
    R = propagate(A, [np.array([0, 2])], alpha=1.0, K=1)[:, 0]
    assert R[1] == 0.0  # net input .25 - .75 < 0 -> rectified


def test_sign_zero_blocks_output():
    A = build_operator(np.array([0, 1]), np.array([1, 2]), np.array([1., 1.]), np.array([1., 0., 1.]), 3)
    R = propagate(A, [np.array([0])], alpha=1.0, K=3)[:, 0]
    assert R[1] > 0 and R[2] == 0


def test_unit_means():
    R = np.array([[1., 0.], [3., 2.], [5., 4.]])
    out = unit_means(R, pd.array([7, 7, pd.NA], dtype="Int64"))
    assert out.loc[7].tolist() == [2.0, 1.0]
