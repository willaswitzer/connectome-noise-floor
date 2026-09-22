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
- Resolved: the Shiu checkpoint pair is not among the 134 targets (see Week 5 checkpoint section) — it runs as a standalone single-dataset validation.
- `build_units_v2.py`'s hb2ct matching does not split comma-joined hemibrain_type strings; measured impact is 3 of 7,579 candidate units changing membership (2 of 6,903 primary) and 5 new units, affecting 0 of the 134 draft targets. Left unfixed since the cost of touching a pre-registered artifact outweighs the near-zero effect on the actual target set.

## Week 5 checkpoint (Shiu reproduction)

The sugar-GRN → MN9 validation is a single-dataset FlyWire reproduction check, run separately from the 134 cross-connectome targets — it is not part of the pre-registered target draw and does not require hemibrain or male CNS.

- **Stimulus:** right-hemisphere labellar sugar GRNs, FlyWire type LB3 (unit 4935 in `units_v2_flywire_neurons.csv`), root IDs listed in `Drosophila_brain_model/example.ipynb` cell 3 / `figures.ipynb`. 20 of 21 root IDs matched our annotation release; one (720575940620900446) was not found and needs checking against the annotation version used.
- **Readout:** MN9, FlyWire type CB0701 (unit 1436), root ID 720575940660219265 (used in `example.ipynb` cell 17, target: 67.0 Hz at 100 Hz sugar stimulation). A second candidate ID (720575940645521262) appears once in `figures.ipynb` cell 17 as "right" but is inconsistent with other notebook comments — use ...219265 as primary and treat the second as a robustness check if time allows.
- **Pass criterion:** LIF model reproduces ~67 Hz MN9 firing rate at 100 Hz sugar-GRN stimulation before proceeding to Week 6–8 cross-connectome runs.
- **Result (Sep 21, 2026):** reproduced. MN9 fired at 66.43 Hz at 100 Hz sugar-GRN stimulation (target: 67.0 Hz, per example.ipynb cell 17), within expected tolerance given brian2 2.10.1/numpy 2.5 vs. the pinned 2.5.1/1.24 environment. Individual silencing of the 3 most active sugar GRNs gave 64.70, 68.73, 70.50 Hz — non-monotonic, which is worth a second look before relying on single-neuron ablation results, but doesn't affect the baseline pass criterion.

## Neurotransmitter sign assignment

Every neuron in the simulation needs a sign, assigned in this priority order:

1. **Primary:** `top_nt` where `used_in_analysis == True` (Eckstein et al.'s own QC'd calls).
2. **Fallback:** for a neuron without that, the majority `top_nt` among other `used_in_analysis == True` neurons of the same type.
3. **Deeper fallback:** if no `used_in_analysis` neuron exists anywhere in that type, the majority of any named `conf_nt` (excluding "unknown") across all matched neurons of that type — tag the unit's sign as low-confidence in output.
4. If a type has zero named calls at all, it has no assignable sign and is excluded — never defaulted to a guessed transmitter.

## Amendment 1 — 2026-09-22 (before any results)

**Status at time of amendment:** No cross-connectome or error-injection results have been computed or viewed. The only model run so far is the Shiu `example.ipynb` environment check on FlyWire v630, which does not touch the analysis targets. Everything below is pre-registered. The original decisions 1–5 still stand unless a section below changes them.

**Reason:** The project scope expanded (rev. 2, Sep 22). The noise-ceiling study becomes Arm 1, and a reconstruction-error arm (Arm 2) is added. Arm 2 is compared against Arm 1's cross-animal disagreement.

---

### A1.1 Comparison hierarchy (replaces the H4 framing)

All comparisons are made at the **matched cell-type level**, using the same units and readouts in both arms.

| Rung | Pair | Varies | Units |
|---|---|---|---|
| R0 | FlyWire left vs. FlyWire right | hemisphere only (same animal, same pipeline) | FlyWire types present on both sides |
| R1 | FlyWire vs. male CNS | individual, pipeline, sex | `build_units_v2.py` (7,579 units) |
| R2 | FlyWire vs. hemibrain | individual, pipeline (both female) | `build_unit_hemibrain_map.py` |
| R3 (exploratory) | hemibrain vs. male CNS | individual, pipeline, sex | via FlyWire-type bridge |

**H4 (revised):** Replicability ranks R0 > R2 > R1, per model and per readout. The test is a paired comparison over the targets shared by each pair of rungs, with the same permutation baseline as H1.

**Interpretation limits (pre-stated):**
- R0 contains hemisphere asymmetry *and* left/right type-matching noise. It is not a pure noise floor.
- The R1 and R2 gaps mix biology, reconstruction pipeline, and matching error. We do not attribute them to biology alone.
- R1 additionally includes sex. Units flagged as sexually dimorphic or sex-specific stay excluded, as in the original decision. Known dimorphic circuits are reported separately.

### A1.2 Coverage per pair

Each target in `targets_DRAFT.csv` gets a coverage flag for R0, R1, and R2. A target enters a rung's analysis only if **both** its seed type and ≥90% of its readout units lie inside both datasets' reconstructed volumes. <!-- CHECK: is 90% the right threshold? Pick before running. -->

- **Hemibrain:** the volume is taken from its ROI list. If the SEZ is outside it (to verify), sugar GRN → MN9 is excluded from R2.
- **Reporting:** the paper reports the target count for each rung. The draw itself is **not** redone.

### A1.3 Global synaptic weight (Wsyn)

Both conditions are pre-registered, and both are reported:
- **(a) Fixed:** Wsyn is held at the FlyWire-calibrated value on every connectome. This is the **primary** condition.
- **(b) Recalibrated:** Wsyn is re-tuned per connectome, using the same calibration target (sugar → MN9 rate). It is applied only on connectomes that contain the circuit. Otherwise it is matched to FlyWire's median firing rate over a fixed stimulus set. <!-- CHECK: fallback calibration target acceptable? -->

If (a) and (b) disagree in the ordering of the rungs, the disagreement is reported as a main-text caveat.

### A1.4 Edge weights

- **Primary:** normalized input fractions, meaning synapses from i to j divided by the total input synapses onto j, within each dataset.
- **Sensitivity check:** raw synapse counts.

**Resolved against the existing weight spec:** the pre-amendment spec uses raw synapse counts as the connection weight — `model.py` sets `syn.w = df_con['Excitatory x Connectivity'] * w_syn`, i.e. weight scales linearly with raw synapse count per edge, and the Metric section's conserved-edge filter also thresholds on raw counts (`>10 synapses`). **This amendment changes the primary weight metric**: normalized input fractions become primary, and raw synapse counts move from the implicit default to the sensitivity check. Everything already run (the Shiu checkpoint) used raw-count weighting and is unaffected, since it predates this amendment and is a single-dataset reproduction check, not a cross-connectome target.

### A1.5 Model assumptions (stated, not tested)

The LIF model has zero basal firing and no gap junctions. Its E/I sign comes from predicted neurotransmitter via the four-tier rule. The linear model shares the sign rule. These are named in the paper as the main structure-to-function limitation.

### A1.6 Readouts (shared by both arms)

1. **Primary (unchanged):** Spearman correlation of response vectors per stratum.
2. **Secondary:** sign and magnitude agreement for named downstream targets, including sugar → MN9 where covered.
3. **Secondary:** rank agreement of the top-k responding cell types. k = 20 is primary, and k = 10 and 50 are sensitivity checks. <!-- CHECK: k choice -->

### A1.7 H3 (sign ablation), scope clarified

H3 stays in Arm 1 as a within-model sanity check. The formal comparison of sign errors against cross-animal disagreement is done in Arm 2 (neurotransmitter-mislabel condition) and is not double-counted.

### A1.8 Arm 2 (stub; must be completed and committed before any Arm 2 output is viewed)

**Reference connectome:** FlyWire v783. The sweep measures sensitivity to errors *added on top of* FlyWire's residual errors.

| Error type | Unit of rate | Injection rule |
|---|---|---|
| False merge | merges per 1,000 neurons | join the synapse lists of same-neuropil or spatially adjacent pairs |
| False split | splits per 1,000 neurons | divide synapses by a spatial cut |
| Missed synapse | fraction of synapses removed | biased toward weak edges |
| Spurious synapse | fraction of synapses added | added between spatially co-located pairs |
| NT mislabel | fraction of neurons with a flipped sign | weighted by low prediction confidence |
| Section loss (optional) | slab thickness | split all neurites crossing the slab |

- Each type is swept separately, then combined. There are ≥5 seeds per rate.
- The "realistic" anchor for each type is fixed from a published source before any results. <!-- TODO: fill in citations and values -->
- **Break-even rate:** the rate at which the median drift (1 − readout agreement against unperturbed FlyWire) equals the R0, R1, and R2 disagreement. It is reported per error type, model, and readout, with bootstrap CIs over targets.

**Real-error condition:** FlyWire v630 vs. v783, run on the **same synapse table**. If the synapse tables can't be matched, this condition is dropped and reported as a limitation.
