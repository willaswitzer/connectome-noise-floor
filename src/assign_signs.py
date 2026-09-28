"""Per-neuron transmitter sign for each dataset (PREREG A3.4; A2.1 s_i).

Sign map: acetylcholine, dopamine, serotonin, octopamine -> +1;
          gaba, glutamate, histamine -> -1. Anything else is not a call.
Tiers (first that yields a transmitter wins):
  0  known NT for the neuron's matched unit, shared by all three datasets:
     strict majority (>50%) of parsed FlyWire `known_nt` over the unit's FlyWire neurons,
     else strict majority of male CNS `ground_truth` over its male neurons.
     If both exist and disagree, no tier 0. A `known_nt` entry counts only if it names
     exactly one mapped transmitter (neuropeptides, NO, tyramine, glycine and '-negative'
     tokens ignored). Hemibrain bodies mapping to >1 unit get tier 0 only if all agree.
  1  neuron-level prediction passing the dataset's QC:
     hemibrain top_nt if used_in_analysis; FlyWire top_nt if top_nt_conf >= 0.5;
     male CNS predicted_nt if predicted_nt_confidence >= 0.5.
  2  strict majority of tier-1 calls within the dataset's own cell type.
  3  strict majority of any named prediction within the type (hemibrain conf_nt;
     FlyWire top_nt; male CNS predicted_nt), flagged low-confidence.
  4  no call -> sign 0.
Also writes the predictions-only sensitivity (tiers 1-4, no tier 0).
"Majority" is implemented as strict (>50% of calls); ties fall through to the next tier.
Structural only: reads no model output.

Outputs: data/derived/signs_{flywire,male_cns,hemibrain}.parquet
         results/sign_summary.csv"""
import re
import pandas as pd

SIGN = {"acetylcholine": 1, "dopamine": 1, "serotonin": 1, "octopamine": 1,
        "gaba": -1, "glutamate": -1, "histamine": -1}


def clean(x):
    return x if isinstance(x, str) and x in SIGN else None


def parse_known(s):
    if not isinstance(s, str):
        return None
    toks = {t.strip().lower() for t in re.split(r"[;,]", s)}
    hits = {t for t in toks if t in SIGN}
    return hits.pop() if len(hits) == 1 else None


def strict_majority(s):
    s = s.dropna()
    if s.empty:
        return None
    vc = s.value_counts()
    return vc.index[0] if vc.iloc[0] > len(s) / 2 else None


def tiers(df, key, t1, t3, t0=None):
    """df has columns: id, type, t1 (tier-1 call or None), t3 (any named call or None)."""
    grp2 = df.dropna(subset=["type"]).groupby("type")["t1"].agg(strict_majority)
    grp3 = df.dropna(subset=["type"]).groupby("type")["t3"].agg(strict_majority)
    out = {}
    for use0 in ([True, False] if t0 is not None else [False]):
        nt = pd.Series(None, index=df.index, dtype=object)
        tier = pd.Series(4, index=df.index)
        cands = ([(0, df["t0"])] if use0 else []) + [
            (1, df["t1"]), (2, df["type"].map(grp2)), (3, df["type"].map(grp3))]
        for k, c in cands:
            fill = nt.isna() & c.notna()
            nt[fill], tier[fill] = c[fill], k
        out[use0] = (nt, tier)
    res = pd.DataFrame({key: df["id"], "type": df["type"]})
    nt, tier = out[t0 is not None]
    res["nt"], res["tier"], res["sign"] = nt, tier, nt.map(SIGN).fillna(0).astype(int)
    nt, tier = out[False]
    res["nt_predonly"], res["tier_predonly"] = nt, tier
    res["sign_predonly"] = nt.map(SIGN).fillna(0).astype(int)
    return res


# ---------- load ----------
fwa = pd.read_csv("data/flywire_783/annotations_v3_neurons.tsv", sep="\t", low_memory=False,
                  usecols=["root_id", "cell_type", "top_nt", "top_nt_conf", "known_nt"])
fwu = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")[["root_id", "unit_id"]]
fwa = fwa.merge(fwu, on="root_id", how="left")

mca = pd.read_feather("data/male_cns/body-annotations-male-cns-v1.0-minconf-0.5.feather",
                      columns=["bodyId", "type", "superclass"]).dropna(subset=["superclass"])
mcn = pd.read_feather("data/male_cns/body-neurotransmitters-male-cns-v1.0.feather",
                      columns=["body", "predicted_nt", "predicted_nt_confidence", "ground_truth"])
mca = mca.merge(mcn.rename(columns={"body": "bodyId"}), on="bodyId", how="left")
mcu = pd.read_csv("data/derived/units_v2_male_neurons.csv")[["bodyId", "unit_id"]]
mca = mca.merge(mcu, on="bodyId", how="left")

hbn = pd.read_parquet("data/hemibrain/neurons.parquet", columns=["bodyId", "type"])
hbt = pd.read_csv("data/hemibrain/neuron_nt.csv", usecols=["bodyid", "top_nt", "conf_nt", "used_in_analysis"])
hbn = hbn.merge(hbt.rename(columns={"bodyid": "bodyId"}), on="bodyId", how="left")
hbm = pd.read_csv("data/derived/unit_hemibrain_map.csv")[["bodyId", "unit_id"]].drop_duplicates()

# ---------- tier 0 per unit ----------
fw_known = fwa.assign(k=fwa["known_nt"].map(parse_known)).dropna(subset=["unit_id"]) \
              .groupby("unit_id")["k"].agg(strict_majority)
mc_known = mca.assign(k=mca["ground_truth"].map(clean)).dropna(subset=["unit_id"]) \
              .groupby("unit_id")["k"].agg(strict_majority)
k = pd.concat([fw_known.rename("fw"), mc_known.rename("mc")], axis=1)
k["t0"] = k["fw"].where(k["fw"].notna(), k["mc"])
k.loc[k["fw"].notna() & k["mc"].notna() & (k["fw"] != k["mc"]), "t0"] = None
unit_t0 = k["t0"]
n_conflict = int((k["fw"].notna() & k["mc"].notna() & (k["fw"] != k["mc"])).sum())

hb_t0 = (hbm.assign(t0=hbm["unit_id"].map(unit_t0))
         .groupby("bodyId")["t0"].agg(lambda s: s.iloc[0] if s.notna().all() and s.nunique() == 1 else None))

# ---------- per dataset ----------
fw = pd.DataFrame({"id": fwa["root_id"], "type": fwa["cell_type"],
                   "t0": fwa["unit_id"].map(unit_t0),
                   "t1": fwa["top_nt"].where(fwa["top_nt_conf"] >= 0.5).map(clean),
                   "t3": fwa["top_nt"].map(clean)})
mc = pd.DataFrame({"id": mca["bodyId"], "type": mca["type"],
                   "t0": mca["unit_id"].map(unit_t0),
                   "t1": mca["predicted_nt"].where(mca["predicted_nt_confidence"] >= 0.5).map(clean),
                   "t3": mca["predicted_nt"].map(clean)})
hb = pd.DataFrame({"id": hbn["bodyId"], "type": hbn["type"],
                   "t0": hbn["bodyId"].map(hb_t0),
                   "t1": hbn["top_nt"].where(hbn["used_in_analysis"].eq(True)).map(clean),
                   "t3": hbn["conf_nt"].map(clean)})

res = {"flywire": tiers(fw, "root_id", "t1", "t3", t0=True),
       "male_cns": tiers(mc, "bodyId", "t1", "t3", t0=True),
       "hemibrain": tiers(hb, "bodyId", "t1", "t3", t0=True)}
for name, r in res.items():
    r.to_parquet(f"data/derived/signs_{name}.parquet", index=False)

# ---------- summary ----------
rows = []
for name, r in res.items():
    for col in ["tier", "tier_predonly"]:
        vc = r[col].value_counts()
        rows.append({"dataset": name, "rule": "primary" if col == "tier" else "predonly",
                     **{f"tier{t}": int(vc.get(t, 0)) for t in range(5)},
                     "n": len(r),
                     "frac_neg": round((r["sign" if col == "tier" else "sign_predonly"] == -1).mean(), 3)})
summ = pd.DataFrame(rows)

# unit-level sign agreement, FlyWire vs male CNS (matched units)
def unit_sign(r, idcol, units, col):
    x = r.merge(units, left_on=idcol, right_on=units.columns[0])
    return x.groupby("unit_id")[col].agg(lambda s: strict_majority(s[s != 0]))

agree = []
for col in ["sign", "sign_predonly"]:
    a = unit_sign(res["flywire"], "root_id", fwu, col)
    b = unit_sign(res["male_cns"], "bodyId", mcu, col)
    j = pd.concat([a.rename("fw"), b.rename("mc")], axis=1).dropna()
    agree.append((col, len(j), round((j["fw"] != j["mc"]).mean(), 4)))

summ.to_csv("results/sign_summary.csv", index=False)
print(summ.to_string(index=False))
print(f"\ntier-0 units: {int(unit_t0.notna().sum()):,} of {len(k):,} with any known NT; "
      f"FlyWire/male CNS known-NT conflicts skipped: {n_conflict}")
for col, n, d in agree:
    print(f"FlyWire vs male CNS unit-sign disagreement ({col}): {d:.2%} of {n:,} units")
tyr = fwa["known_nt"].fillna("").str.contains(r"\btyramine\b|\bglycine\b") & fwa["known_nt"].map(parse_known).isna()
print(f"FlyWire neurons whose known_nt names only tyramine/glycine (no tier 0 from them): {int(tyr.sum())}")
