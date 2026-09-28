"""R0 (FlyWire left vs. right) eligibility, per PREREG A1.2 / Amendment 2.

A matched unit is R0-usable if it has >=1 FlyWire neuron with side 'left' and
>=1 with side 'right' (annotations v3.0.0 `side`). A target is R0-eligible if its
seed unit is R0-usable. The R0 readout vector uses only R0-usable units.
Structural only: no simulation output is read.

Outputs: results/prereg/r0_flags.csv (one row per draft target)
         data/derived/r0_units.csv (every unit with its L/R neuron counts)"""
import pandas as pd

units = pd.read_csv("data/derived/units_v2_flywire_neurons.csv")
ann = pd.read_csv("data/flywire_783/annotations_v3_neurons.tsv", sep="\t",
                  usecols=["root_id", "side"], low_memory=False)
u = units.merge(ann, on="root_id", how="left")

counts = (u.assign(side=u["side"].fillna("missing"))
          .pivot_table(index="unit_id", columns="side", values="root_id",
                       aggfunc="count", fill_value=0))
for c in ["left", "right", "center", "missing"]:
    if c not in counts:
        counts[c] = 0
counts = counts[["left", "right", "center", "missing"]].rename(columns=lambda c: f"n_{c}")
counts["r0_usable"] = (counts["n_left"] > 0) & (counts["n_right"] > 0)
counts.reset_index().to_csv("data/derived/r0_units.csv", index=False)

t = pd.read_csv("results/prereg/targets_DRAFT.csv")
t = t.merge(counts.reset_index(), on="unit_id", how="left")
t["r0_usable"] = t["r0_usable"].fillna(False).astype(bool)
t[["unit_id", "name", "stratum", "tierA", "tierB", "n_left", "n_right", "n_center",
   "n_missing", "r0_usable"]].rename(columns={"r0_usable": "r0"}).to_csv(
    "results/prereg/r0_flags.csv", index=False)

print(f"units: {len(counts):,}; R0-usable: {int(counts['r0_usable'].sum()):,} "
      f"({counts['r0_usable'].mean():.1%})")
print(f"neurons with side missing: {int(counts['n_missing'].sum()):,}")
print(f"targets R0-eligible: {int(t['r0_usable'].sum())} / {len(t)}")
print(t.groupby("stratum")["r0_usable"].agg(["sum", "count"]))
print(t.groupby("tierB")["r0_usable"].agg(["sum", "count"]))
print("not eligible:", t.loc[~t["r0_usable"], ["name", "n_left", "n_right", "n_center"]].to_string(index=False))
