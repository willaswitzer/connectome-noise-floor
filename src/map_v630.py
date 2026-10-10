"""Map FlyWire v783 neurons to their v630 root IDs (PREREG A8.6). Run on the laptop, which has
CAVE access; the compute workspace cannot reach the FlyWire servers.

Each v783 neuron's annotation `supervoxel_id` is looked up in the chunkedgraph at the
materialization-630 timestamp. Structural only: reads no model output.

Output: data/derived/v630_map.csv  (root_id_783, supervoxel_id, root_id_630)
        plus a short summary printed at the end (one-to-one count, v630 merges, coverage).
Resumable: finished batches are kept in data/derived/v630_map.partial.csv.

Usage: .venv/bin/python src/map_v630.py [--datastack flywire_fafb_public] [--batch 10000]
If CAVE asks for a token: in Python, `from caveclient import CAVEclient;
CAVEclient().auth.setup_token(make_new=False)` and follow the link, then
`CAVEclient().auth.save_token(token="<token>")`."""
import argparse, time
from pathlib import Path

import numpy as np
import pandas as pd
from caveclient import CAVEclient

ap = argparse.ArgumentParser()
ap.add_argument("--datastack", default="flywire_fafb_public")
ap.add_argument("--version", type=int, default=630)
ap.add_argument("--batch", type=int, default=10000)
a = ap.parse_args()

OUT = Path("data/derived/v630_map.csv")
PART = Path("data/derived/v630_map.partial.csv")

ann = pd.read_csv("data/flywire_783/annotations_v3_neurons.tsv", sep="\t",
                  usecols=["root_id", "supervoxel_id"], low_memory=False).drop_duplicates("root_id")
ann = ann.rename(columns={"root_id": "root_id_783"})
no_sv = ann["supervoxel_id"].isna() | (ann["supervoxel_id"] == 0)
print(f"v783 neurons: {len(ann):,}; without a supervoxel_id: {int(no_sv.sum()):,}", flush=True)
todo = ann[~no_sv].copy()
todo["supervoxel_id"] = todo["supervoxel_id"].astype(np.uint64)

done = pd.read_csv(PART) if PART.exists() else pd.DataFrame(columns=["root_id_783", "supervoxel_id", "root_id_630"])
todo = todo[~todo["root_id_783"].isin(done["root_id_783"])]
print(f"to look up: {len(todo):,} (already done: {len(done):,})", flush=True)

client = CAVEclient(a.datastack)
ts = client.materialize.get_timestamp(a.version)
print(f"datastack {a.datastack}, version {a.version}, timestamp {ts}", flush=True)

t0 = time.time()
for i in range(0, len(todo), a.batch):
    chunk = todo.iloc[i:i + a.batch]
    for attempt in range(5):
        try:
            roots = client.chunkedgraph.get_roots(chunk["supervoxel_id"].to_numpy(), timestamp=ts)
            break
        except Exception as e:  # network hiccups: back off and retry
            print(f"  batch {i}: {type(e).__name__}: {e}; retry {attempt + 1}/5", flush=True)
            time.sleep(10 * (attempt + 1))
    else:
        raise SystemExit("giving up after 5 retries; rerun to resume")
    out = chunk.assign(root_id_630=np.asarray(roots, dtype=np.uint64))
    out.to_csv(PART, mode="a", header=not PART.exists(), index=False)
    n = i + len(chunk)
    print(f"  {n:,}/{len(todo):,}  {time.time() - t0:.0f}s", flush=True)

m = pd.read_csv(PART)
m.to_csv(OUT, index=False)
m = m[m["root_id_630"] != 0]
fan = m["root_id_630"].value_counts()
merged = set(fan[fan > 1].index)
one = m[~m["root_id_630"].isin(merged)]
print("\nsummary")
print(f"  mapped to a v630 root: {len(m):,} of {len(ann):,}")
print(f"  one-to-one: {len(one):,}")
print(f"  v630 roots receiving >1 v783 neuron (v630 merges): {len(merged):,} "
      f"covering {int(m['root_id_630'].isin(merged).sum()):,} v783 neurons")
c630 = Path("Drosophila_brain_model/2023_03_23_connectivity_630_final.parquet")
if c630.exists():
    e = pd.read_parquet(c630)
    cols = [c for c in ("Presynaptic_ID", "Postsynaptic_ID") if c in e.columns] or list(e.columns[:2])
    ids = pd.unique(pd.concat([e[c] for c in cols]))
    print(f"  one-to-one roots present in the v630 connectivity table ({cols}): "
          f"{int(one['root_id_630'].isin(ids).sum()):,}")
print(f"wrote {OUT}")
