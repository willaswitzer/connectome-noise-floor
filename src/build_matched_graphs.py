"""Build common-schema matched connectivity graphs, one parquet per dataset.

Schema: bodyId_pre, bodyId_post, weight, unit_id_pre, unit_id_post,
r2_unit_pre, r2_unit_post, dataset.

Graph scope per dataset follows PREREG A2.7; r2_unit_* follows A2.8 (run
build_r2_units.py first). Primary files: flywire, male_cns (brain only),
hemibrain (Traced/Roughly traced + typed Leaves). Sensitivity files:
flywire_princeton, male_cns_cns (whole CNS), hemibrain_all (all bodies).
weight is the raw synapse count between the pair. FlyWire is written twice
(PREREG A1.3): flywire.parquet from Buhmann et al. synapses (primary, Shiu's
Connectivity_783.parquet, already one row per pair) and flywire_princeton.parquet
from the Princeton table (sensitivity; per-neuropil syn_count summed to one row
per pair). Male CNS and hemibrain are already one row per pair.

unit_id comes from the v2 matching tables (build_units_v2.py /
build_unit_hemibrain_map.py), which share one unit_id space across all three
datasets -- matched_units_malecns.csv (build_units.py, v1) is NOT used here
since it is built from a different graph and its unit_id values do not
correspond to units_v2_flywire_neurons.csv's.

unit_id is null where a neuron isn't in a matched unit at all, and also
where a hemibrain bodyId maps to more than one unit_id (comma-joined
hemibrain_type fan-out; see PREREG.md limitations) -- both cases are treated
as unmatched rather than guessed, since unmatched/untyped neurons still
carry signal in simulation and should be kept, not dropped.

Each dataset's edges are restricted to its own traced/annotated neuron set:
FlyWire's connections file and hemibrain's edges.parquet (built from
status="Traced") are already scoped that way; the male CNS weights table is
not (it's the full autosegmentation, ~88M bodyIds) and is filtered here to
both endpoints being in the annotation table."""
import pandas as pd

OUT_DIR = "data/derived/graphs"


def report(name, edges):
    n_neurons = pd.unique(pd.concat([edges["bodyId_pre"], edges["bodyId_post"]])).shape[0]
    both = edges["unit_id_pre"].notna() & edges["unit_id_post"].notna()
    neither = edges["unit_id_pre"].isna() & edges["unit_id_post"].isna()
    print(f"{name}: neurons={n_neurons:,} edges={len(edges):,} "
          f"both_matched={int(both.sum()):,} ({both.mean():.1%}) "
          f"neither_matched={int(neither.sum()):,} ({neither.mean():.1%})")


# --- FlyWire ---
# PREREG A1.3: Buhmann et al. synapses (the table Shiu's Wsyn was calibrated on)
# are PRIMARY -> flywire.parquet. Princeton (2025) synapses are the SENSITIVITY
# table -> flywire_princeton.parquet. Both are unthresholded, one row per pair.
fw_id2unit = pd.read_csv("data/derived/units_v2_flywire_neurons.csv").set_index("root_id")["unit_id"]
unit2r2 = pd.read_csv("data/derived/r2_units.csv").set_index("unit_id")["r2_unit_id"]


def add_r2_from_units(edges):
    edges["r2_unit_pre"] = edges["unit_id_pre"].map(unit2r2).astype("Int64")
    edges["r2_unit_post"] = edges["unit_id_post"].map(unit2r2).astype("Int64")


def finish_flywire(edges, name, path):
    edges["unit_id_pre"] = edges["bodyId_pre"].map(fw_id2unit).astype("Int64")
    edges["unit_id_post"] = edges["bodyId_post"].map(fw_id2unit).astype("Int64")
    add_r2_from_units(edges)
    edges["dataset"] = name
    report(name, edges)
    edges.to_parquet(path, index=False)


fw_edges = pd.read_parquet("data/flywire_783_buhmann/Connectivity_783.parquet",
                           columns=["Presynaptic_ID", "Postsynaptic_ID", "Connectivity"])
fw_edges = fw_edges.rename(columns={"Presynaptic_ID": "bodyId_pre", "Postsynaptic_ID": "bodyId_post",
                                    "Connectivity": "weight"})
finish_flywire(fw_edges, "flywire", f"{OUT_DIR}/flywire.parquet")
del fw_edges

fw_edges = pd.read_csv("data/flywire_783/connections_princeton_no_threshold.csv.gz")
fw_edges = (fw_edges.groupby(["pre_root_id", "post_root_id"], as_index=False)["syn_count"].sum()
            .rename(columns={"pre_root_id": "bodyId_pre", "post_root_id": "bodyId_post",
                              "syn_count": "weight"}))
finish_flywire(fw_edges, "flywire_princeton", f"{OUT_DIR}/flywire_princeton.parquet")
del fw_edges

# --- male CNS (A2.7) ---
# The raw weights table covers the full autosegmentation (88.4M bodyIds). Nodes are
# restricted to annotated bodies with a superclass. Primary drops VNC superclasses so
# activity cannot loop through the VNC (FlyWire has none); sensitivity keeps the whole CNS.
VNC = {"vnc_intrinsic", "vnc_sensory", "vnc_motor", "vnc_efferent"}
mc_ann = pd.read_feather("data/male_cns/body-annotations-male-cns-v1.0-minconf-0.5.feather",
                         columns=["bodyId", "superclass"]).dropna(subset=["superclass"])
cns_ids = set(mc_ann["bodyId"])
brain_ids = set(mc_ann.loc[~mc_ann["superclass"].isin(VNC), "bodyId"])

mc_edges = pd.read_feather("data/male_cns/connectome-weights-male-cns-v1.0-minconf-0.5.feather")
mc_edges = mc_edges[mc_edges["body_pre"].isin(cns_ids) & mc_edges["body_post"].isin(cns_ids)]
mc_edges = mc_edges.rename(columns={"body_pre": "bodyId_pre", "body_post": "bodyId_post"})

mc_id2unit = pd.read_csv("data/derived/units_v2_male_neurons.csv").set_index("bodyId")["unit_id"]
mc_edges["unit_id_pre"] = mc_edges["bodyId_pre"].map(mc_id2unit).astype("Int64")
mc_edges["unit_id_post"] = mc_edges["bodyId_post"].map(mc_id2unit).astype("Int64")
add_r2_from_units(mc_edges)

mc_edges["dataset"] = "male_cns_cns"
report("male_cns_cns", mc_edges)
mc_edges.to_parquet(f"{OUT_DIR}/male_cns_cns.parquet", index=False)

mc_edges = mc_edges[mc_edges["bodyId_pre"].isin(brain_ids) & mc_edges["bodyId_post"].isin(brain_ids)].copy()
mc_edges["dataset"] = "male_cns"
report("male_cns", mc_edges)
mc_edges.to_parquet(f"{OUT_DIR}/male_cns.parquet", index=False)
del mc_edges

# --- hemibrain (A2.7, A2.8) ---
# status="Traced" returns 97,900 bodies, 76,161 of them cropped "Leaves" fragments.
# Primary nodes: statusLabel Traced / Roughly traced, plus any body mapped to a unit.
# unit_id_*: ambiguous bodies (mapping to >1 unit) stay null (sensitivity behaviour).
# r2_unit_*: every mapped body gets its R2 unit (primary for R2).
hb_edges = pd.read_parquet("data/hemibrain/edges.parquet")
hb_neu = pd.read_parquet("data/hemibrain/neurons.parquet", columns=["bodyId", "statusLabel"])

hb_map = pd.read_csv("data/derived/unit_hemibrain_map.csv")[["bodyId", "unit_id"]].drop_duplicates()
dup_counts = hb_map["bodyId"].value_counts()
ambiguous = set(dup_counts[dup_counts > 1].index)
hb_id2unit = hb_map[~hb_map["bodyId"].isin(ambiguous)].set_index("bodyId")["unit_id"]
hb_id2r2 = pd.read_csv("data/derived/r2_body_map.csv").set_index("bodyId")["r2_unit_id"]

hb_edges["unit_id_pre"] = hb_edges["bodyId_pre"].map(hb_id2unit).astype("Int64")
hb_edges["unit_id_post"] = hb_edges["bodyId_post"].map(hb_id2unit).astype("Int64")
hb_edges["r2_unit_pre"] = hb_edges["bodyId_pre"].map(hb_id2r2).astype("Int64")
hb_edges["r2_unit_post"] = hb_edges["bodyId_post"].map(hb_id2r2).astype("Int64")
print(f"hemibrain: {len(ambiguous):,} bodyIds map to >1 unit_id -- null unit_id, "
      "but assigned an R2 unit")

hb_edges["dataset"] = "hemibrain_all"
report("hemibrain_all", hb_edges)
hb_edges.to_parquet(f"{OUT_DIR}/hemibrain_all.parquet", index=False)

keep = set(hb_neu.loc[hb_neu["statusLabel"].isin(["Traced", "Roughly traced"]), "bodyId"]) | set(hb_id2r2.index)
hb_edges = hb_edges[hb_edges["bodyId_pre"].isin(keep) & hb_edges["bodyId_post"].isin(keep)].copy()
hb_edges["dataset"] = "hemibrain"
report("hemibrain", hb_edges)
hb_edges.to_parquet(f"{OUT_DIR}/hemibrain.parquet", index=False)
