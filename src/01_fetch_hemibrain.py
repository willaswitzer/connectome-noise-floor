"""Fetch hemibrain v1.2.1 neurons + connections from neuPrint into parquet.

Setup (once):
    pip install neuprint-python pandas pyarrow tqdm
    # Get a token: log in at https://neuprint.janelia.org with Google,
    # click your account (top right) -> "Auth Token" -> copy.
    export NEUPRINT_TOKEN="<paste token>"

Run:
    python src/01_fetch_hemibrain.py

Outputs (data/hemibrain/):
    neurons.parquet   one row per traced neuron (bodyId, type, instance, roiInfo, NT columns if present)
    edges.parquet     one row per neuron pair (bodyId_pre, bodyId_post, weight, roi info)
"""
import os
from pathlib import Path

from dotenv import load_dotenv
load_dotenv()

from neuprint import Client, NeuronCriteria as NC, fetch_neurons, fetch_adjacencies

DATASET = "hemibrain:v1.2.1"
OUT = Path("data/hemibrain")
MIN_WEIGHT = 1  # keep everything raw; apply the >10 synapse / >1% filters downstream

token = os.environ.get("NEUPRINT_TOKEN")
if not token:
    raise SystemExit("Set NEUPRINT_TOKEN first (see docstring).")

OUT.mkdir(parents=True, exist_ok=True)
client = Client("neuprint.janelia.org", dataset=DATASET, token=token)
print("Connected:", client)

# 1. Neurons: traced only (the proofread set used in cross-connectome typing)
neuron_df, roi_counts = fetch_neurons(NC(status="Traced"))
print(f"{len(neuron_df)} traced neurons")
print("Columns:", list(neuron_df.columns))  # check for predictedNt / consensusNt here
neuron_df.to_parquet(OUT / "neurons.parquet")
roi_counts.to_parquet(OUT / "neuron_roi_counts.parquet")

# 2. Connections: fetch_adjacencies batches internally and is the documented
# path for large neuron sets (fetch_simple_connections explicitly is not).
body_ids = neuron_df["bodyId"].tolist()
_, roi_conn_df = fetch_adjacencies(
    sources=NC(bodyId=body_ids),
    targets=NC(bodyId=body_ids),
    min_total_weight=MIN_WEIGHT,
    export_dir=str(OUT / "adjacency_export"),  # also dumps CSVs as a checkpoint
)

# roi_conn_df is per-(bodyId_pre, bodyId_post, roi); collapse to one row per
# neuron pair, matching what downstream analysis (edge filtering) expects.
edges = (
    roi_conn_df.groupby(["bodyId_pre", "bodyId_post"], as_index=False)["weight"]
    .sum()
)
print(f"{len(edges)} edges, total synapses {edges['weight'].sum():,}")
edges.to_parquet(OUT / "edges.parquet")
roi_conn_df.to_parquet(OUT / "edges_by_roi.parquet")

# 3. Sanity checks
print("Typed neurons:", neuron_df["type"].notna().sum(), "of", len(neuron_df))
print("Unique types:", neuron_df["type"].nunique())
print("Edges with weight >10:", (edges["weight"] > 10).sum())
