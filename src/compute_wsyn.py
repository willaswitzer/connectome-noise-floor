"""Per-dataset LIF Wsyn by structural rescaling (PREREG A1.4; edge definition in Amendment 2).

    Wsyn_k = 0.275 mV * median_FW / median_k

median_X = median, over matched edges, of the unit-pair edge strength in dataset X.
Unit-pair edge strength = total synapses from unit u to unit v divided by the number
of connected neuron pairs (u_i -> v_j with >=1 synapse), i.e. mean synapses per
connected pair, so it does not scale with how many neurons each dataset puts in a unit.
Matched edge = unit pair (u != v) present in both datasets and passing the conserved-edge
filter in both: strength > 10 synapses OR total u->v synapses >= 1% of v's total input
from matched units.
Structural only: reads no simulation output. Run after build_matched_graphs.py.

Output: results/wsyn.csv"""
import pandas as pd

G = "data/derived/graphs"
W_REF = 0.275  # mV, Shiu et al. 2024


def unit_edges(path):
    e = pd.read_parquet(path, columns=["weight", "unit_id_pre", "unit_id_post"]).dropna()
    e = e[e["unit_id_pre"] != e["unit_id_post"]]
    g = e.groupby(["unit_id_pre", "unit_id_post"])["weight"].agg(total="sum", n_pairs="size")
    g["strength"] = g["total"] / g["n_pairs"]
    g["frac_input"] = g["total"] / g.groupby(level="unit_id_post")["total"].transform("sum")
    g["conserved"] = (g["strength"] > 10) | (g["frac_input"] >= 0.01)
    return g


def rescale(ref, other):
    j = ref.join(other, how="inner", lsuffix="_ref", rsuffix="_k")
    j = j[j["conserved_ref"] & j["conserved_k"]]
    m_ref, m_k = j["strength_ref"].median(), j["strength_k"].median()
    return len(j), m_ref, m_k, W_REF * m_ref / m_k


fw = unit_edges(f"{G}/flywire.parquet")
rows = [("flywire", len(fw), None, None, W_REF)]
for name in ["male_cns", "hemibrain", "flywire_princeton"]:
    n, m_ref, m_k, w = rescale(fw, unit_edges(f"{G}/{name}.parquet"))
    rows.append((name, n, m_ref, m_k, w))
    print(f"{name}: matched conserved edges={n:,}  median FW={m_ref:.2f}  median {name}={m_k:.2f}  "
          f"Wsyn={w:.4f} mV")

pd.DataFrame(rows, columns=["dataset", "n_matched_edges", "median_flywire", "median_dataset",
                            "wsyn_mV"]).to_csv("results/wsyn.csv", index=False)
