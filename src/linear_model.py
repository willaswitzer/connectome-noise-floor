"""Linear propagation model (PREREG A2.1).

    W_ij = s_i * n_ij / sum_i' n_i'j      (input fraction onto j, over all its presynaptic
                                            partners in the dataset's scoped graph)
    r_0  = seed vector (1 on every seed neuron)
    r_k  = max(0, alpha * W^T r_{k-1})     (rectified each hop)
    R    = sum_{k=1..K} r_k                (primary alpha = 0.5, K = 5)

s_i comes from data/derived/signs_<dataset>.parquet (A3.4); neurons with no sign have
s_i = 0 (they receive input but pass none on). Many seeds are propagated at once as the
columns of a dense matrix. Unit responses are the mean of R over the unit's neurons that
are nodes of the graph (A2.1); the caller decides which units to read out and excludes the
seed unit (A2.1, A2.3)."""
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp

GRAPHS = "data/derived/graphs"
SIGNS = {"flywire": ("signs_flywire.parquet", "root_id"),
         "flywire_princeton": ("signs_flywire.parquet", "root_id"),
         "male_cns": ("signs_male_cns.parquet", "bodyId"),
         "male_cns_cns": ("signs_male_cns.parquet", "bodyId"),
         "hemibrain": ("signs_hemibrain.parquet", "bodyId"),
         "hemibrain_all": ("signs_hemibrain.parquet", "bodyId")}


@dataclass
class Graph:
    name: str
    ids: np.ndarray            # node id per index
    A: sp.csr_matrix           # post x pre, entries s_pre * n / in_total_post
    unit: pd.Series            # node index -> unit_id (Int64, NA if unmatched)
    r2: pd.Series              # node index -> r2_unit_id


def build_operator(pre, post, weight, sign_of_pre, n):
    """pre/post: node indices; weight: synapse counts; sign_of_pre: sign per node index."""
    in_total = np.bincount(post, weights=weight, minlength=n)
    vals = sign_of_pre[pre] * weight / in_total[post]
    return sp.csr_matrix((vals.astype(np.float64), (post, pre)), shape=(n, n))


def load_graph(name, sign_col="sign"):
    e = pd.read_parquet(f"{GRAPHS}/{name}.parquet")
    ids, inv = np.unique(np.concatenate([e["bodyId_pre"].to_numpy(), e["bodyId_post"].to_numpy()]),
                         return_inverse=True)
    pre, post = inv[:len(e)], inv[len(e):]
    f, key = SIGNS[name]
    s = pd.read_parquet(f"data/derived/{f}", columns=[key, sign_col]).set_index(key)[sign_col]
    sign = pd.Series(ids).map(s).fillna(0).to_numpy(dtype=np.float64)
    A = build_operator(pre, post, e["weight"].to_numpy(dtype=np.float64), sign, len(ids))

    def node_map(col_pre, col_post):
        m = pd.concat([pd.Series(e[col_pre].to_numpy(), index=pre),
                       pd.Series(e[col_post].to_numpy(), index=post)])
        m = m[~m.index.duplicated()].sort_index()
        return m.reindex(range(len(ids))).astype("Int64")

    return Graph(name, ids, A, node_map("unit_id_pre", "unit_id_post"),
                 node_map("r2_unit_pre", "r2_unit_post"))


def propagate(A, seeds, alpha=0.5, K=5, rectify=True):
    """seeds: list of arrays of node indices (one per perturbation). Returns N x B array R.
    rectify=False is the unrectified sensitivity model (Amendment 5): negative values are
    kept as suppression."""
    n, B = A.shape[0], len(seeds)
    r = np.zeros((n, B))
    for b, idx in enumerate(seeds):
        r[idx, b] = 1.0
    R = np.zeros_like(r)
    AT = A  # A is post x pre, so A @ r gives input to each post neuron
    for _ in range(K):
        r = alpha * (AT @ r)
        if rectify:
            r = np.maximum(0.0, r)
        R += r
    return R


def unit_means(R, labels):
    """labels: per-node group label (NA = not read out). Returns DataFrame groups x B."""
    lab = pd.Series(labels).reset_index(drop=True)
    ok = lab.notna().to_numpy()
    codes, groups = pd.factorize(lab[ok])
    M = sp.csr_matrix((np.ones(ok.sum()), (codes, np.flatnonzero(ok))), shape=(len(groups), len(lab)))
    counts = np.asarray(M.sum(axis=1)).ravel()
    return pd.DataFrame((M @ R) / counts[:, None], index=pd.Index(groups, name="unit"))
