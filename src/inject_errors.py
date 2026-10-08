"""Arm-2 reconstruction-error injector (PREREG Amendment 8, A8.1-A8.3).

Acts on the FlyWire v783 scoped graph (A2.7; Buhmann synapses, no threshold). Structural only:
this module produces perturbed edge lists and signs; it computes no model output.

Error types (rates are fractions; anchor x m gives the A8.1 rate):
  miss   missed synapses: binomial thinning of every edge count, p = rate.
  split  false splits: each edge END (the pre side for the presynaptic neuron, the post side for
         the postsynaptic neuron) is detached independently with probability rate and moved to
         that neuron's single new fragment node. Expected detached synapses per neuron = rate x
         its synapses. Fragments stay in the graph, belong to no unit, and carry the parent's sign.
  merge  false merges: neuron i receives F_i = S_i * rate / (1 - rate) foreign synapses, so the
         foreign fraction of its synapses is `rate`. Foreign synapses are copies of a donor's
         edge ends (out-edge d->k becomes i->k, in-edge k->d becomes k->i); each donor edge end is
         copied with probability min(1, F_i / S_d). Donor: uniform among neurons sharing i's top
         neuropil (Princeton table). If one donor is too small, further donors are drawn for the
         remainder. Donors keep their edges. Copies that would be i->i are dropped.
  spur   spurious synapses: rate x total synapses are added one at a time; neuropil drawn
         proportional to its synapse count, pre proportional to output synapses there, post
         proportional to input synapses there. Self-pairs are redrawn.
  flip   NT mislabel: round(rate x eligible) neurons with sign != 0 have their sign flipped, drawn
         without replacement with probability proportional to (1 - top_nt_conf).

Conditions (A8.1): single type at m; combined; lopsided; unproofread. Caps: 0.5 for miss, split,
merge, flip; 1.0 for spur. Single-type points above a cap are dropped; in combined conditions a
capped component is held at its cap and flagged.

Randomness: injection seed k uses default_rng(20260930 + k); each component draws from its own
child stream (rng.spawn), so a given k gives the same draws for a component across conditions.
Order: miss -> split (on surviving original synapses, at split/(1 - miss) so recall loss is
miss + split exactly in expectation) -> merge (copies of ORIGINAL donor edges) -> spur -> flip.
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

ANCHOR = {"miss": 0.008, "split": 0.008, "merge": 0.008, "spur": 0.05, "flip": 0.06}
CAP = {"miss": 0.5, "split": 0.5, "merge": 0.5, "spur": 1.0, "flip": 0.5}
TYPES = list(ANCHOR)
M_GRID = [0.5, 1, 2, 3, 5, 10, 25, 50]
M_UNPROOFREAD = (1 - 0.799) / (1 - 0.992)   # 25.125 (A8.1, P6)
SEED_BASE = 20260930


# ----------------------------------------------------------------------------- conditions

def condition_rates(kind, m=None, t=None):
    """Return (rates dict, capped list) or (None, reason) if the point is dropped.
    kind: 'single' (needs t), 'combined', 'lopsided', 'unproofread'."""
    if kind == "single":
        r = ANCHOR[t] * m
        if r > CAP[t]:
            return None, f"{t} rate {r:.3f} > cap {CAP[t]}"
        return {t: r}, []
    if kind == "combined":
        rates = {"miss": 0.004 * m, "split": 0.004 * m, "merge": 0.008 * m,
                 "spur": 0.05 * m, "flip": 0.06 * m}
    elif kind == "lopsided":
        rates = {"miss": 0.006 * m, "split": 0.006 * m, "merge": 0.004 * m,
                 "spur": 0.05 * m, "flip": 0.06 * m}
    elif kind == "unproofread":
        mu = M_UNPROOFREAD
        rates = {"miss": 0.004 * mu, "split": 0.004 * mu, "merge": 0.008 * mu,
                 "spur": 0.05, "flip": 0.06}
    else:
        raise ValueError(kind)
    capped = [k for k, v in rates.items() if v > CAP[k]]
    for k in capped:
        rates[k] = CAP[k]
    return rates, capped


def grid():
    """All linear Arm-2 grid points (A8.3) as (label, kind, m, t, rates, capped)."""
    out = []
    for t in TYPES:
        for m in M_GRID:
            rates, capped = condition_rates("single", m, t)
            if rates is not None:
                out.append((f"{t}_m{m:g}", "single", m, t, rates, capped))
    for kind in ("combined", "lopsided"):
        for m in M_GRID:
            rates, capped = condition_rates(kind, m)
            out.append((f"{kind}_m{m:g}", kind, m, None, rates, capped))
    rates, capped = condition_rates("unproofread")
    out.append(("unproofread", "unproofread", M_UNPROOFREAD, None, rates, capped))
    return out


# ----------------------------------------------------------------------------- base graph

@dataclass
class Base:
    ids: np.ndarray          # node id per index (sorted)
    pre: np.ndarray          # edge pre index (int64)
    post: np.ndarray         # edge post index
    w: np.ndarray            # synapse count (int64)
    sign: np.ndarray         # per node, float (+1, -1, 0)
    flip_weight: np.ndarray  # per node, 1 - top_nt_conf (0 where not eligible)
    top_np: np.ndarray       # per node, top-neuropil code (-1 = not in Princeton table)
    np_out: list = field(default_factory=list)  # per neuropil: (node idx array, out-synapse weights)
    np_in: list = field(default_factory=list)   # per neuropil: (node idx array, in-synapse weights)
    np_total: np.ndarray = None                 # per neuropil synapse count
    np_names: list = field(default_factory=list)

    @property
    def n(self):
        return len(self.ids)

    def syn_per_node(self):
        return (np.bincount(self.pre, weights=self.w, minlength=self.n)
                + np.bincount(self.post, weights=self.w, minlength=self.n))


def build_base(ids, pre, post, w, sign, conf, princeton):
    """princeton: DataFrame with pre_root_id, post_root_id, neuropil, syn_count (any rows; rows
    with an endpoint outside `ids` are ignored). conf: per node top_nt_conf (NaN allowed)."""
    n = len(ids)
    sign = np.asarray(sign, dtype=np.float64)
    conf = np.asarray(conf, dtype=np.float64)
    elig = sign != 0
    fw = 1.0 - conf
    # signed neurons without a confidence get the median flip weight of the eligible set
    med = np.nanmedian(fw[elig]) if elig.any() else 0.0
    fw = np.where(np.isnan(fw), med, fw)
    fw = np.where(elig, fw, 0.0)

    p = princeton
    pi = np.searchsorted(ids, p["pre_root_id"].to_numpy())
    qi = np.searchsorted(ids, p["post_root_id"].to_numpy())
    pi_ok = (pi < n) & (ids[np.minimum(pi, n - 1)] == p["pre_root_id"].to_numpy())
    qi_ok = (qi < n) & (ids[np.minimum(qi, n - 1)] == p["post_root_id"].to_numpy())
    ok = pi_ok & qi_ok
    pi, qi = pi[ok], qi[ok]
    npc, np_names = pd.factorize(p["neuropil"].to_numpy()[ok])
    sc = p["syn_count"].to_numpy()[ok].astype(np.float64)
    K = len(np_names)
    # per-node, per-neuropil synapse totals via a sparse table
    import scipy.sparse as sp
    O = sp.csr_matrix((sc, (npc, pi)), shape=(K, n))
    I = sp.csr_matrix((sc, (npc, qi)), shape=(K, n))
    tot = (O + I).tocsc()
    top = np.full(n, -1, dtype=np.int64)
    has = np.diff(tot.indptr) > 0
    T = tot.toarray() if K * n < 5e7 else None
    if T is not None:
        top[has] = T[:, has].argmax(axis=0)
    else:  # column-wise argmax on CSC
        for j in np.flatnonzero(has):
            s, e = tot.indptr[j], tot.indptr[j + 1]
            top[j] = tot.indices[s + np.argmax(tot.data[s:e])]
    np_out, np_in = [], []
    for k in range(K):
        r = O.getrow(k); np_out.append((r.indices.astype(np.int64), r.data))
        r = I.getrow(k); np_in.append((r.indices.astype(np.int64), r.data))
    np_total = np.asarray(O.sum(axis=1)).ravel()
    return Base(ids=np.asarray(ids), pre=np.asarray(pre, np.int64), post=np.asarray(post, np.int64),
                w=np.asarray(w, np.int64), sign=sign, flip_weight=fw, top_np=top,
                np_out=np_out, np_in=np_in, np_total=np_total, np_names=list(np_names))


def load_flywire_base(graph="data/derived/graphs/flywire.parquet",
                      signs="data/derived/signs_flywire.parquet",
                      ann="data/flywire_783/annotations_v3_neurons.tsv",
                      princeton="data/flywire_783/connections_princeton_no_threshold.csv.gz"):
    e = pd.read_parquet(graph, columns=["bodyId_pre", "bodyId_post", "weight"])
    ids, inv = np.unique(np.concatenate([e["bodyId_pre"].to_numpy(), e["bodyId_post"].to_numpy()]),
                         return_inverse=True)
    pre, post = inv[:len(e)], inv[len(e):]
    s = pd.read_parquet(signs, columns=["root_id", "sign"]).set_index("root_id")["sign"]
    sign = pd.Series(ids).map(s).fillna(0).to_numpy(np.float64)
    c = (pd.read_csv(ann, sep="\t", usecols=["root_id", "top_nt_conf"], low_memory=False)
         .drop_duplicates("root_id").set_index("root_id")["top_nt_conf"])
    conf = pd.Series(ids).map(c).to_numpy(np.float64)
    p = pd.read_csv(princeton, usecols=["pre_root_id", "post_root_id", "neuropil", "syn_count"])
    return build_base(ids, pre, post, e["weight"].to_numpy(np.int64), sign, conf, p)


# ----------------------------------------------------------------------------- components

def _miss(pre, post, w, p, rng):
    w2 = rng.binomial(w, 1.0 - p)
    k = w2 > 0
    return pre[k], post[k], w2[k]


def _split(pre, post, w, q, n, rng):
    """Detach each edge end with prob q to the owning neuron's fragment. Returns new arrays, the
    parent index of each fragment (fragments are numbered n, n+1, ...)."""
    sp_ = rng.random(len(pre)) < q
    sq_ = rng.random(len(pre)) < q
    owners = np.unique(np.concatenate([pre[sp_], post[sq_]]))
    frag_of = np.full(n, -1, dtype=np.int64)
    frag_of[owners] = n + np.arange(len(owners))
    pre2 = np.where(sp_, frag_of[pre], pre)
    post2 = np.where(sq_, frag_of[post], post)
    return pre2, post2, w.copy(), owners


def _csr(idx, n):
    order = np.argsort(idx, kind="stable")
    ptr = np.concatenate([[0], np.cumsum(np.bincount(idx, minlength=n))])
    return order, ptr


def _pick_donors(recv, top_np, rng, n):
    """Uniform donor sharing the receiver's top neuropil, != receiver. Fallback: any node."""
    order = np.argsort(top_np, kind="stable")
    codes = top_np[order]
    starts = np.searchsorted(codes, top_np[recv], side="left")
    ends = np.searchsorted(codes, top_np[recv], side="right")
    size = ends - starts
    d = np.empty(len(recv), dtype=np.int64)
    pool_ok = (top_np[recv] >= 0) & (size > 1)
    # draw among the size-1 others: pick r in [0, size-1), skip own position
    r = (rng.random(len(recv)) * np.maximum(size - 1, 1)).astype(np.int64)
    cand = order[np.minimum(starts + r, len(order) - 1)]
    # if the candidate is the receiver itself, take the last element of the group instead
    last = order[np.maximum(ends - 1, 0)]
    cand = np.where(cand == recv, last, cand)
    d[pool_ok] = cand[pool_ok]
    g = ~pool_ok
    if g.any():
        x = (rng.random(g.sum()) * (n - 1)).astype(np.int64)
        x = x + (x >= recv[g])
        d[g] = x
    return d


def _merge(pre0, post0, w0, rate, n, top_np, rng, max_rounds=20, chunk=20000):
    """Return foreign edges (pre, post, w) to ADD, built from the original edges."""
    S = (np.bincount(pre0, weights=w0, minlength=n) + np.bincount(post0, weights=w0, minlength=n))
    need = S * rate / (1.0 - rate)
    o_ord, o_ptr = _csr(pre0, n)
    i_ord, i_ptr = _csr(post0, n)
    add_pre, add_post, add_w = [], [], []
    recv = np.flatnonzero(need > 0)
    rem = need[recv]
    for _ in range(max_rounds):
        if len(recv) == 0:
            break
        d = _pick_donors(recv, top_np, rng, n)
        Sd = S[d]
        q = np.where(Sd > 0, np.minimum(1.0, rem / np.maximum(Sd, 1e-12)), 0.0)
        for c0 in range(0, len(recv), chunk):
            rr, dd, qq = recv[c0:c0 + chunk], d[c0:c0 + chunk], q[c0:c0 + chunk]
            for ordr, ptr, is_out in ((o_ord, o_ptr, True), (i_ord, i_ptr, False)):
                deg = ptr[dd + 1] - ptr[dd]
                tot = int(deg.sum())
                if tot == 0:
                    continue
                rep = np.repeat(np.arange(len(dd)), deg)
                off = np.arange(tot) - np.repeat(np.cumsum(deg) - deg, deg)
                eidx = ordr[ptr[dd][rep] + off]
                keep = rng.random(tot) < qq[rep]
                eidx, rcv = eidx[keep], rr[rep[keep]]
                if is_out:   # d -> k  becomes  i -> k
                    a, b = rcv, post0[eidx]
                else:        # k -> d  becomes  k -> i
                    a, b = pre0[eidx], rcv
                ok = a != b
                add_pre.append(a[ok]); add_post.append(b[ok]); add_w.append(w0[eidx][ok])
        rem = rem - q * Sd
        more = rem > 0.5          # stop once less than half a synapse is outstanding
        recv, rem = recv[more], rem[more]
    if not add_pre:
        z = np.zeros(0, np.int64)
        return z, z, z
    return np.concatenate(add_pre), np.concatenate(add_post), np.concatenate(add_w)


def _spur(base, total, rate, rng, chunk=4_000_000):
    """Spurious synapses, generated in chunks and aggregated to (pre, post, count) within each
    chunk so memory stays bounded at rate 1.0 (~54M synapses on v783)."""
    n_add = int(round(rate * total))
    if n_add == 0:
        z = np.zeros(0, np.int64)
        return z, z, z
    n = base.n
    per = rng.multinomial(n_add, base.np_total / base.np_total.sum())
    keys, cnts = [], []
    for k, cnt in enumerate(per):
        if cnt == 0:
            continue
        oi, ow = base.np_out[k]
        ii, iw = base.np_in[k]
        oc, ic = np.cumsum(ow) / ow.sum(), np.cumsum(iw) / iw.sum()
        left = int(cnt)
        while left > 0:
            m = min(left, chunk)
            x = oi[np.minimum(np.searchsorted(oc, rng.random(m), side="right"), len(oi) - 1)]
            y = ii[np.minimum(np.searchsorted(ic, rng.random(m), side="right"), len(ii) - 1)]
            ok = x != y                        # self-pairs are redrawn on the next pass
            u, c = np.unique(x[ok] * n + y[ok], return_counts=True)
            keys.append(((u // n).astype(np.int32), (u % n).astype(np.int32)))
            cnts.append(c.astype(np.int32))
            left -= int(ok.sum())
            del x, y, ok, u, c
    # pairs are aggregated within each chunk only; a pair drawn in two chunks appears twice,
    # which weight_matrix() and edges() sum
    pre = np.concatenate([a for a, _ in keys]); post = np.concatenate([b for _, b in keys])
    return pre, post, np.concatenate(cnts)


def _flip(base, rate, rng):
    elig = np.flatnonzero(base.flip_weight > 0)
    k = int(round(rate * len(elig)))
    if k == 0:
        return np.zeros(0, np.int64)
    p = base.flip_weight[elig] / base.flip_weight[elig].sum()
    return np.sort(rng.choice(elig, size=k, replace=False, p=p))


def aggregate(pre, post, w, n_total):
    """Sum duplicate (pre, post) rows."""
    return _aggregate_keys(pre.astype(np.int64) * n_total + post.astype(np.int64), w, n_total)


def _aggregate_keys(key, w, n_total):
    order = np.argsort(key)
    key = key[order]; w = w[order]
    del order
    starts = np.concatenate([[0], np.flatnonzero(key[1:] != key[:-1]) + 1])
    ww = np.add.reduceat(w, starts).astype(np.int64)
    u = key[starts]
    return u // n_total, u % n_total, ww


# ----------------------------------------------------------------------------- injection

@dataclass
class Perturbed:
    parts: list               # [(pre, post, w), ...]: kept original edges, then each addition
    sign: np.ndarray          # length n_total (nodes + fragments), after flips
    n_nodes: int              # original node count; fragments are n_nodes .. n_total-1
    frag_parent: np.ndarray   # parent node index of each fragment
    flipped: np.ndarray       # node indices whose sign was flipped
    rates: dict
    capped: list
    k: int
    retained_ends: np.ndarray = None   # per original node: original synapse ends still attached
    added_ends: np.ndarray = None      # per original node: ends of merge/spurious additions
    _agg: tuple = None

    @property
    def n_total(self):
        return self.n_nodes + len(self.frag_parent)

    def edges(self):
        """Aggregated (pre, post, w), one row per pair (cached). Memory-heavy at the largest
        merge/spurious rates; the linear model uses weight_matrix() instead."""
        if self._agg is None:
            nt = self.n_total
            key = np.concatenate([a.astype(np.int64) * nt + b for a, b, _ in self.parts])
            w = np.concatenate([c for _, _, c in self.parts])
            self._agg = _aggregate_keys(key, w, nt)
        return self._agg

    pre = property(lambda self: self.edges()[0])
    post = property(lambda self: self.edges()[1])
    w = property(lambda self: self.edges()[2])

    def weight_matrix(self, consume=False):
        """Synapse-count matrix W (post x pre, CSR, float64) with duplicate pairs summed,
        built from one int32 COO. consume=True releases the parts as soon as they are copied
        (the Perturbed then has no edges left), which keeps peak memory bounded at the
        largest grid points."""
        import scipy.sparse as sp
        nt = self.n_total
        rows = np.concatenate([b for _, b, _ in self.parts])
        cols = np.concatenate([a for a, _, _ in self.parts])
        data = np.concatenate([c for _, _, c in self.parts]).astype(np.float64)
        if consume:
            self.parts.clear(); self._agg = None
        W = sp.coo_matrix((data, (rows, cols)), shape=(nt, nt))
        del rows, cols, data
        W = W.tocsr()
        W.sum_duplicates()
        return W


def linear_operator(pert, sign=None, consume=False, chunk=2_000_000):
    """A2.1 operator from a Perturbed graph: A_ij = s_j * n_ij / sum_j' n_ij' (rows = post).
    Same result as linear_model.build_operator on the aggregated edges. Scaled in place, in
    chunks of stored entries."""
    W = pert.weight_matrix(consume=consume)
    s = pert.sign if sign is None else sign
    in_total = np.asarray(W.sum(axis=1)).ravel()
    inv = np.divide(1.0, in_total, out=np.zeros_like(in_total), where=in_total > 0)
    ptr = W.indptr
    for r0 in range(0, W.shape[0], 20000):
        r1 = min(r0 + 20000, W.shape[0])
        a, b = ptr[r0], ptr[r1]
        W.data[a:b] *= np.repeat(inv[r0:r1], np.diff(ptr[r0:r1 + 1]))
        W.data[a:b] *= s[W.indices[a:b]]
    return W


def inject(base, rates, k, capped=()):
    """Apply the rate dict (keys from TYPES) with injection seed k."""
    rng = np.random.default_rng(SEED_BASE + k)
    streams = dict(zip(TYPES, rng.spawn(len(TYPES))))
    n = base.n
    pre, post, w = base.pre, base.post, base.w
    miss, split = rates.get("miss", 0.0), rates.get("split", 0.0)
    if miss > 0:
        pre, post, w = _miss(pre, post, w, miss, streams["miss"])
    frag_parent = np.zeros(0, np.int64)
    if split > 0:
        pre, post, w, frag_parent = _split(pre, post, w, split / (1.0 - miss), n, streams["split"])
    retained = _ends(pre, post, w, n)
    parts = [(pre, post, w)]
    if rates.get("merge", 0) > 0:
        parts.append(_merge(base.pre, base.post, base.w, rates["merge"], n, base.top_np,
                            streams["merge"]))
    if rates.get("spur", 0) > 0:
        parts.append(_spur(base, int(base.w.sum()), rates["spur"], streams["spur"]))
    added = sum((_ends(p[0], p[1], p[2], n) for p in parts[1:]), np.zeros(n))
    parts = [tuple(x.astype(np.int32) for x in p) for p in parts]   # n_total < 2**31; w fits
    flipped = (_flip(base, rates["flip"], streams["flip"]) if rates.get("flip", 0) > 0
               else np.zeros(0, np.int64))
    sign = base.sign.copy()
    sign[flipped] *= -1
    sign = np.concatenate([sign, sign[frag_parent]])   # fragments carry the parent's (final) sign
    return Perturbed(parts, sign, n, frag_parent, flipped, dict(rates), list(capped), k,
                     retained, added)


def seed_fix_columns(pert, seed_idx):
    """Seed-exempt NT flips (A8.2): node indices (seed neurons and their fragments) whose sign was
    flipped and must be restored for a perturbation seeded at seed_idx."""
    if len(pert.flipped) == 0:
        return np.zeros(0, np.int64)
    s = np.intersect1d(seed_idx, pert.flipped)
    if len(s) == 0:
        return s
    frags = pert.n_nodes + np.flatnonzero(np.isin(pert.frag_parent, s))
    return np.concatenate([s, frags])


# ----------------------------------------------------------------------------- audit

def _ends(pre, post, w, n):
    """Synapse ends per original node (fragments excluded)."""
    a, b = pre < n, post < n
    return (np.bincount(pre[a], weights=w[a], minlength=n)
            + np.bincount(post[b], weights=w[b], minlength=n))


def audit(base, pert):
    """Realized rates (structure only). Per neuron, on its own synapse ends:
    recall loss = 1 - (original ends still on the neuron) / (original ends);
    foreign     = added ends / (retained + added ends), counting both ends of every added
                  synapse (a merge copy is foreign for the receiver AND for its partner).
    Pooled over neurons (synapse-weighted) and as the median over neurons with >= 100 synapses."""
    S0 = base.syn_per_node()
    R, A = pert.retained_ends, pert.added_ends
    has, big = S0 > 0, S0 >= 100
    return {"recall_loss": float(1 - R[has].sum() / S0[has].sum()),
            "foreign_frac": float(A.sum() / (R.sum() + A.sum())),
            "recall_loss_median_neuron": float(np.median(1 - R[big] / S0[big])),
            "foreign_frac_median_neuron": float(np.median(A[big] / np.maximum(R[big] + A[big], 1))),
            "syn_total_ratio": float(sum(int(c.sum(dtype=np.int64)) for _, _, c in pert.parts) / base.w.sum()),
            "n_fragments": int(len(pert.frag_parent)), "n_flipped": int(len(pert.flipped)),
            "flip_frac": float(len(pert.flipped) / max((base.flip_weight > 0).sum(), 1))}
