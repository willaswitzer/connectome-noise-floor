# Pre-registration: connectome noise ceiling

Draft target list: `targets_DRAFT.csv` (134 targets).

## Variability-noted targets

Nine targets have `var > 0`: at least one male-CNS neuron in the unit carries a
`matchingNotes` annotation containing "variab" or "difference in number"
(see `src/eligibility.py`). The `var` column is not part of the `primary` or
`strict` eligibility rules.

**Decision:** these 9 targets are kept in the primary analysis.

**Sensitivity check:** the primary metric will be reported again restricted to
targets with `var == 0` (125 targets).

## Metric

For each target unit and each dataset pair (male CNS–FlyWire for all 134 targets; all three pairwise + three-way comparisons among male CNS/FlyWire/hemibrain for the 58 Tier-B targets), activate the target and record the steady-state response (linear model) or time-averaged spike rate (LIF model) of every matched unit shared across the datasets being compared. Compute Spearman's rank correlation between the two response vectors. Report median + IQR of these correlations per stratum and pooled, separately for: linear vs. LIF, unrestricted vs. conserved-edge-filtered (>10 synapses or >1% of target input) connectomes, and correctly-signed vs. sign-ablated neurotransmitters.

## Primary hypotheses (pre-registered — anything else is exploratory)

- **H1:** Median cross-dataset Spearman rho is significantly above a label-shuffled permutation baseline, for both model classes.
- **H2:** Conserved-edge filtering increases median replicability vs. the unfiltered connectome.
- **H3:** Sign-ablating neurotransmitter identity decreases median replicability vs. the correctly-signed model.
- **H4:** Three-way replicability (Tier B, all three datasets) is ≤ two-way replicability (male CNS–FlyWire) on the same 58 units.

**Exploratory, not pre-registered:** LIF vs. linear differences and per-stratum breakdowns.

## Limitations

- Tier B ascending stratum has only 4 targets — treat three-way ascending results as suggestive, not confirmatory.
- `hb_share` is a type-matching estimate, not yet checked against real hemibrain connectivity — re-verify once hemibrain data is fetched.
- Reported run count (326) covers dataset×target combinations only, not the additional multiples from model class, edge-filtering, and sign-ablation conditions.
- Confirm whether Shiu's sugar-neuron→MN9 checkpoint pair is included in the 134 targets; if not, it must be validated as a separate step before Week 5.
