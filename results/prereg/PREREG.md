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

Revised same day (Sep 22) before any results; supersedes c5315b7.

**Status at time of amendment.** No cross-connectome or error-injection results have been computed or viewed. The only model output so far is the Week-5 Shiu checkpoint (MN9 at 66.43 Hz vs. the published 67.0 Hz). That is a single-dataset validation outside the target draw. Everything below is pre-registered. Decisions 1–5 and the Tier A/B target draw stand unless changed below.

**Reason.** The project scope expanded (rev. 2, Sep 22). The noise-ceiling study becomes Arm 1. A reconstruction-error arm (Arm 2) is added and compared against Arm 1's cross-animal disagreement. Two public projects surfaced during the pre-submission literature check (brainlab, Sep 15, 2026; flybench). They show that absolute-weight LIF output depends on each dataset's synapse-detection sensitivity and on which FlyWire synapse table is used. A1.3 and A1.4 respond to that.

---

### A1.1 Comparison hierarchy (replaces the H4 framing)

All comparisons are made at the **matched cell-type level**, with the same units and readouts in both arms.

| Rung | Pair | What varies | Units |
|---|---|---|---|
| R0 | FlyWire left vs. FlyWire right | hemisphere only (same animal, same pipeline) | FlyWire types present on both sides (annotations v3.0.0 `side`) |
| R1 | FlyWire vs. male CNS | individual, pipeline, sex | `build_units_v2.py` (7,579 units) |
| R2 | FlyWire vs. hemibrain | individual, pipeline (both female) | `build_unit_hemibrain_map.py` |
| R3 (exploratory) | hemibrain vs. male CNS | individual, pipeline, sex | bridged through FlyWire types |

**H4 (revised).** Replicability ranks R0 > R2 > R1, per model and per readout. The test is a paired comparison over the targets each pair of rungs shares, against the same permutation baseline as H1.

**Interpretation limits (stated in advance):**
- R0 contains hemisphere asymmetry and left/right matching noise. It is not a pure noise floor.
- The R1 and R2 gaps mix biology, reconstruction pipeline, and matching error. They are not attributed to biology alone.
- R1 also varies sex. Dimorphic and sex-specific units stay excluded, as in the original decisions. Known dimorphic circuits are reported separately.

### A1.2 Coverage per pair

The existing eligibility rules and draw are unchanged.
- **Tier A (80 targets):** male CNS vs. FlyWire. These are the R1 targets.
- **Tier B (58 targets):** all three datasets, `hb_share > 0.5`. These are the R2 targets.
- **R0:** each target gets a flag for whether its seed type and readout units exist in both FlyWire hemispheres. R0 uses every flagged target.

The hemibrain volume excludes the optic lobe, the periesophageal neuropils, and the gnathal ganglia (Scheffer et al. 2020). Sugar GRN → MN9 is therefore excluded from R2 and run on R0 and R1 only. Per-rung target counts are reported.

### A1.3 Synapse table (FlyWire)

- **Primary:** FlyWire v783 with Buhmann et al. synapses. This is the table used to calibrate Shiu's Wsyn and to run the Week-5 checkpoint. It is also required for the real-error condition (A1.9) to compare like with like.
- **Sensitivity:** FlyWire v783 with Princeton (2025) synapses.

The hemibrain and male CNS use their own native synapse tables.

### A1.4 Edge weights and global weight (Wsyn)

Each model uses its own weighting:
- **Linear model:** normalized input fractions, i.e. synapses i→j divided by total input synapses onto j, within each dataset. Raw counts are a sensitivity check.
- **LIF model:** raw synapse counts × Wsyn.

**Primary Wsyn rule: structural rescaling.** For connectome k:

Wsyn_k = 0.275 mV × (median synapse count on matched edges in FlyWire ÷ median synapse count on the same matched edges in k)

- "Matched edges" means edges between matched units, present in both datasets, and passing the conserved-edge filter.
- The calibration uses no simulation output. That way it corrects for detector sensitivity without tuning for agreement.
- The same rule applies to the Princeton sensitivity table relative to Buhmann.

**Sensitivity:** Wsyn is fixed at 0.275 mV on every connectome. If the rung ordering differs between the primary and sensitivity conditions, that difference is reported as a main-text caveat.

**Resolved against the existing weight spec:** the pre-amendment spec (and `model.py`'s `syn.w = df_con['Excitatory x Connectivity'] * w_syn`) used raw synapse counts as the weight, with no separate rule for the linear vs. LIF model. This amendment splits the two: the **linear model's** primary weighting changes from raw counts to normalized input fractions (raw counts move to its sensitivity check); the **LIF model** keeps raw synapse counts as its edge weight, now scaled by the structural-rescaling Wsyn rule above instead of a fixed constant. The Week-5 Shiu checkpoint predates this amendment, used a fixed Wsyn and raw-count LIF weighting, and is unaffected.

### A1.5 Model assumptions (stated, not tested)

- **LIF:** zero basal firing and no gap junctions.
- **Both models:** E/I sign comes from the predicted neurotransmitter, using the four-tier rule.

These are named in the paper as the main structure-to-function limitation.

### A1.6 Readouts (shared by both arms)

1. **Primary (unchanged):** Spearman correlation of response vectors per stratum.
2. **Secondary:** sign and magnitude agreement for named downstream targets, including sugar → MN9 where the rung covers it.
3. **Secondary:** rank agreement of the top-20 responding cell types. k = 10 and k = 50 are sensitivity checks.

### A1.7 H3 (sign ablation), scope clarified

H3 stays in Arm 1 as a within-model sanity check. The formal comparison of sign errors against cross-animal disagreement belongs to Arm 2 (neurotransmitter-mislabel condition). The two are not double-counted.

### A1.8 Arm 2 error model (stub; the protocol is finalized before any Arm 2 output is viewed)

**Reference connectome:** FlyWire v783 (Buhmann synapses). The sweep measures sensitivity to errors added *on top of* FlyWire's residual errors.

| Error type | Unit of rate | Injection rule | Realistic anchor (fixed now) |
|---|---|---|---|
| Missed synapse | fraction of a neuron's synapses removed | biased toward weak edges and distal twigs | Unproofread vs. proofread segments: median assigned-synapse F1 is 0.799 vs. 0.992 (Dorkenwald et al. 2022, Nat Methods). Twig attachment is the main accuracy limit (Dorkenwald et al. 2024). The anchor is a ~1% loss for proofread data and ~20% for unproofread data. |
| Spurious synapse | fraction of synapses added | added between spatially co-located pairs | Buhmann et al. 2021 report ~95% connection accuracy at the ≥5-synapse threshold. The anchor is ~5% of edges. |
| NT mislabel | fraction of neurons with a flipped sign | weighted toward low prediction confidence (`conf_nt_p`) | Accuracy is 94% at the neuron level (Eckstein et al. 2024). The anchor is 6% of neurons. |
| False merge | merges per 1,000 neurons | join the synapse lists of same-neuropil or spatially adjacent pairs | The anchor is derived from the proofreading-round F1 precision component (Dorkenwald et al. 2022). <!-- extract the precision value before finalizing --> |
| False split | splits per 1,000 neurons | divide synapses with a spatial cut | Same source as false merges, using the recall component. <!-- extract before finalizing --> |
| Section loss (optional) | slab thickness | split every neurite crossing the slab | Descriptive only; there is no anchor. |

- **Sweep:** each type separately, then combined. Rates run from 0 to at least 3× the anchor, with ≥5 seeds per rate.
- **Break-even rate:** the rate at which median drift (1 − readout agreement with unperturbed FlyWire) equals the R0, R1, and R2 disagreement. It is reported per error type, model, and readout, with bootstrap CIs over targets.

### A1.9 Real-error condition

The comparison is FlyWire v630 vs. v783, both with Buhmann synapses.

- **Primary:** neurons present in both releases, mapped between root IDs. This isolates proofreading edits to existing neurons, the real-world counterpart of Arm 2's injected errors.
- **Secondary (if time):** all neurons. This adds completeness differences (127,978 neurons in v630 vs. 139,255 in v783).

If root IDs can't be mapped reliably, the condition is dropped and reported as a limitation. It is already on the cut list.

### A1.10 Out of scope

BANC (v888) is not added. It would extend R2 to the SEZ, but its ~0.35× synapse deficit (brainlab, 2026) would add a new confound. It is noted as future work.
