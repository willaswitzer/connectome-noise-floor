"""Resolve male CNS flywireType names to Codex (FlyWire v783) type names.
Rules in order: exact primary; exact additional; composite split; suffix strip."""
import pandas as pd, re

u = pd.read_csv("data/derived/matched_units_malecns.csv")
c = pd.read_csv("data/flywire_783/consolidated_cell_types.csv.gz")
prim = set(c["primary_type"].dropna())
add = set()
for s in c["additional_type(s)"].dropna():
    for p in re.split(r"[;,]", str(s)):
        add.add(p.strip())

def strip_suffix(x):
    m = re.match(r"^(.*?)_?[a-z]$", x)
    return m.group(1) if m else x

def resolve_one(x):
    if x in prim: return "exact_primary", [x]
    if x in add: return "exact_additional", [x]
    core = re.sub(r"^\((.*)\)[a-z]?$", r"\1", x)
    if "," in core:
        parts = [p.strip() for p in core.split(",")]
        out = []
        for p in parts:
            if p in prim or p in add: out.append(p)
            elif strip_suffix(p) in prim or strip_suffix(p) in add: out.append(strip_suffix(p))
        if len(out) == len(parts): return "composite", out
        return "composite_partial", out
    s = strip_suffix(x)
    if s != x and (s in prim or s in add): return "suffix_strip", [s]
    return "unresolved", []

names = sorted(set(u["flywireType"]))
rows = []
for n in names:
    r, t = resolve_one(n)
    rows.append((n, r, ";".join(t)))
al = pd.DataFrame(rows, columns=["male_name", "rule", "codex_types"])
al.to_csv("data/derived/type_alias_table.csv", index=False)

cnt = u["flywireType"].value_counts()
al["neurons"] = al["male_name"].map(cnt)
print(al.groupby("rule").agg(names=("male_name", "count"), neurons=("neurons", "sum")))
print("total neurons", int(cnt.sum()))
for r in ["suffix_strip", "composite", "composite_partial", "unresolved"]:
    s = al[al["rule"] == r].sample(min(8, (al["rule"] == r).sum()), random_state=1)
    print("\n==", r); print(s[["male_name", "codex_types"]].to_string(index=False))
