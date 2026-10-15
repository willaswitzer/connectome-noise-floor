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

## Amendment 2 — 2026-09-28 (before any results)

**Status at time of amendment.** No cross-connectome, linear-model or error-injection output
has been computed or viewed. The only model output remains the Week-5 Shiu checkpoint.
The R0 counts below come from annotations only (`src/r0_flags.py`). Everything in Amendment 1
stands unless changed here.

**Reason.** Amendment 1 names a linear model but never specifies it, and three
implementation details were left open: the R1 target count, the R0 eligibility rule,
and what "matched edges" means in the Wsyn rule. These are fixed here before the
linear model is run.

### A2.1 Linear model (specified)

- **Weights:** W_ij = s_i · n_ij / Σ_i' n_i'j. Here n_ij = synapses i→j, and the denominator is
  j's total input from all presynaptic neurons in that dataset's graph, matched or not (A1.4).
- **Sign s_i:** comes from the four-tier NT rule. GABA and glutamate give −1; every other
  assigned transmitter gives +1, the same binary convention as Shiu's `Excitatory` column.
  Neurons with no assignable sign (tier 4) have s_i = 0: they receive input but pass none on.
- **Propagation:** truncated rectified cascade. r_0 = seed vector (1 on every seed neuron,
  0 elsewhere); r_k[j] = max(0, α · Σ_i W_ij · r_{k−1}[i]); response R = Σ_{k=1..K} r_k.
  **Primary: α = 0.5, K = 5.** Sensitivity: α ∈ {0.3, 0.7}, K ∈ {3, 10}.
- This replaces the phrase "steady-state response" in the Metric section. A fixed point
  of a rectified signed network need not exist or be unique. The truncated cascade is
  deterministic and cheap enough to run on every target and every Arm-2 sweep point.
- **Unit response:** the mean of R over the unit's neurons in that dataset. Unmatched
  neurons propagate activity but are not read out. The seed unit is excluded from its
  own response vector.

### A2.2 Stimulation and readout per rung

- **R1 (FlyWire vs. male CNS):** stimulate all seed-unit neurons in each dataset,
  both sides. Read out over all matched units present in both datasets.
- **R2 (FlyWire vs. hemibrain):** same as R1, on units present in both.
  **Sensitivity:** FlyWire stimulated and read out on right-side neurons only, to
  mirror hemibrain's mostly right-hemisphere volume. If the R2 rung position differs
  between primary and sensitivity, it is reported as a main-text caveat.
- **R0 (FlyWire left vs. right):** stimulate the seed unit's left-side neurons in one
  run and its right-side neurons in another. The response vector is indexed by
  (unit, ipsilateral | contralateral), where ipsilateral means neurons on the stimulated side. It is
  computed over R0-usable units only; midline (`center`) neurons are excluded.

### A2.3 Correlation set

The Spearman correlation is computed over the rung's shared units, excluding units whose
response is zero in both members of the pair. If fewer than 10 units remain, that
target×rung is reported as NA and counted.
**Sensitivity:** include double-zero units.

### A2.4 Target counts per rung (clarifies A1.2)

Tier A (80) and Tier B (58) do not overlap, and together they make up the 134 targets. All 134 have
matched units in FlyWire and male CNS.
- **R1:** 134 targets. A1.2's "Tier A (80): these are the R1 targets" is corrected.
  The Metric section already says 134.
- **R2:** the 58 Tier B targets.
- **R0:** 134 targets. A unit is R0-usable if it has ≥1 FlyWire neuron with `side` = left
  and ≥1 with `side` = right (annotations v3.0.0); 7,476 of 7,579 units qualify. A target is
  R0-eligible if its seed unit is R0-usable; all 134 are. See `results/prereg/r0_flags.csv`.

### A2.5 Matched edges for the Wsyn rule (A1.4)

- **Unit-pair edge strength:** total u→v synapses ÷ number of connected neuron pairs (u_i→v_j,
  ≥1 synapse), i.e. mean synapses per connected pair. This keeps the strength independent of
  how many neurons each dataset assigns to a unit.
- **Matched edge:** a unit pair (u ≠ v) present in both datasets and passing the
  conserved-edge filter in both. The filter is strength > 10 OR u→v synapses ≥ 1% of v's
  total input from matched units.
- Wsyn_k = 0.275 mV × median_FlyWire ÷ median_k over those edges. Computed by
  `src/compute_wsyn.py`, which reads no simulation output.

### A2.6 Implementation record

- The FlyWire graph (`data/derived/graphs/flywire.parquet`) was first built on 2026-09-22
  from the Princeton table, contrary to A1.3. No model was run on it. The graph is
  rebuilt from Buhmann synapses (primary), with Princeton written separately as
  `flywire_princeton.parquet` (sensitivity).
- The project planning doc "rev. 3" (2026-09-28) proposed a ~15-target canonical
  sensory set, a ≥5-synapse edge threshold everywhere, and MN9-tuned Wsyn as a
  sensitivity. None of these is adopted; the 134-target draw and Amendment 1 stand.

### A2.7 Graph scope (which bodies are nodes)

The first graph build (2026-09-28) used different scopes per dataset. Hemibrain
`status="Traced"` returns 97,900 bodies, of which 76,161 are cropped "Leaves" fragments.
The male CNS graph kept all 188,778 annotated bodies, including ~20k VNC neurons and bodies
with no superclass. FlyWire v783 is brain-only proofread neurons. Scopes are fixed as:

- **FlyWire v783:** all neurons in the release (unchanged).
- **Male CNS, primary:** brain neurons only, i.e. bodies with a superclass that is not
  `vnc_intrinsic`, `vnc_sensory`, `vnc_motor` or `vnc_efferent`. Ascending and descending
  neurons are kept, but their VNC partners are not, so activity cannot loop through the VNC
  (FlyWire has no VNC). Bodies with no superclass are excluded, consistent with the
  166,700-neuron definition used in `build_units_v2.py`.
  **Sensitivity:** whole CNS, all bodies with a superclass.
- **Hemibrain, primary:** bodies with statusLabel `Traced` or `Roughly traced`
  (21,739), plus any `Leaves` body that maps to a matched unit (2,880). Untyped Leaves
  fragments are excluded. **Sensitivity:** all 97,900 bodies.
- Input-fraction denominators (A2.1) are computed within the scoped graph.
- These are structural choices, made after seeing only node and edge counts
  from `build_matched_graphs.py`. No model output has been computed.

### A2.8 R2 units: merge units that share hemibrain bodies

Hemibrain types are often coarser than FlyWire types. In `unit_hemibrain_map.csv`,
3,092 of 17,868 mapped hemibrain bodies map to more than one unit (up to 10), touching 816
units. Dropping those bodies, as the first graph build did, leaves 5 of the 58 Tier B seeds
with no hemibrain neurons (LC10d, LC26, OCG01c, VST1, DNp41) and shrinks a sixth (DNpe009).

- **Primary:** for R2 (and R3), the comparison unit is an **R2 unit**: a connected
  component of the bipartite graph linking units to the hemibrain bodies they map to
  (the grouping of Schlegel et al. 2024). In both datasets, a target's R2 seed is the
  R2 unit containing its seed unit, and FlyWire responses are aggregated to R2 units.
  R0 and R1 keep the original units.
- **Sensitivity:** drop ambiguous bodies (the first build's behavior). Targets left with
  no hemibrain seed neurons are reported as NA.
- R2 reports how many targets share an R2 unit with another target, and how many
  R2 units merge more than one original unit.
- This is distinct from the `build_units_v2.py` hb2ct comma-join limitation already
  in PREREG (3 of 7,579 units). That limitation stands.

## Amendment 3 — 2026-09-28 (before any results)

**Status at time of amendment.** Unchanged from Amendment 2. No linear-model, cross-connectome
or error-injection output exists. Counts below come from annotation and matching tables only.

**Reason.** (1) Add a named case-study set alongside the 134-target draw, to support the
named-target readout (A1.6 item 2) and interpretable figures. (2) State how the case-study set
overlaps the draw and how that overlap is counted. (3) State that readouts are brain-only in
every cross-dataset comparison. (4) Fix a sign-rule gap: A2.1 treats histamine as excitatory,
and the four-tier rule's tier-1 column (`used_in_analysis`) exists only in the hemibrain table.

### A3.1 Case-study set (secondary; fixed names)

Seeds are fixed FlyWire v783 cell types (annotations v3.0.0). A seed with several types
stimulates all of them together. Male CNS and hemibrain seeds are the matched units
(R1, R0) or R2 units (R2) containing those types.

| Name | Seed types | Named readout |
|---|---|---|
| sugar_GRN | LB3 | MN9 (CB0701) |
| bitter_GRN | LB1a,LB1d; LB1b; LB1c; LB1e | MN9 (CB0701) |
| DA1_uPN | DA1_lPN | none |
| DM1_uPN | DM1_lPN | none |
| DC1_uPN | DC1_adPN | none |
| VA1v_uPN | VA1v_adPN | none |
| LC4 | LC4 | DNp01, DNp02, DNp11 |
| LPLC2 | LPLC2 | DNp01, DNp02, DNp11 |
| LC6 | LC6 | none |
| LC10a | LC10a | DNa02 |
| JO-A | JO-A1–A4 | none |
| JO-CE | all JO-C* and JO-E* types (14) | none |
| VP2_hot_PN | VP2_adPN | none |
| VP4_dry_PN | VP4_vPN | none |

- A separate water GRN is not included, because FlyWire v783 types sugar and water GRNs together (LB3,
  sub-class "sugar/water").
- Named readouts are reported as the response of that unit, with sign agreement across
  datasets (A1.6 item 2). Perturbations without a named readout use items 1 and 3 only.
- **Coverage per rung.** R1: 14 of 14. R0: 14 of 14 (every seed has neurons on both sides).
  R2: 10 of 14. sugar_GRN, bitter_GRN, JO-A and JO-CE have no hemibrain neurons and are
  excluded. MN9 has no hemibrain neurons, so the gustatory named readout exists only
  on R0 and R1.
- **Role.** The case-study set is descriptive. It does not enter H1–H4, their permutation
  baselines, or the primary break-even rates. Break-even rates for the case-study set are
  reported as secondary. It is also the fixed short list for LIF runs in Arm 2.

### A3.2 Overlap with the 134-target draw

- **bitter_GRN** contains two drawn targets, `LB1a,LB1d` and `LB1c`. The case-study run
  stimulates all four LB1 types together, so it is a different run from either drawn target.
- **LC10a** is not drawn, but its R2 unit (4937, which merges LC10a–f) also contains the drawn
  targets LC10b and LC10d. On R2, those three stimulate the same R2 unit, so it is one run,
  reported once and labelled with all three names.
- No other case-study seed shares a unit or an R2 unit with a drawn target.
- **Counting.** H1–H4 are computed on the 134 drawn targets only. Case-study results are never
  pooled with them. For LC10b and LC10d on R2, the shared run counts once per drawn target, as
  the draw already specified.

### A3.3 Brain-only readout in every cross-dataset comparison

- The readout vector contains only matched units (A2.1). No matched male CNS neuron has a VNC
  superclass (0 of 139,980), and FlyWire and hemibrain have no VNC. So the readout is
  brain cell types in every comparison. This is stated in the methods.
- In the primary male CNS graph (A2.7), VNC neurons are not nodes. Descending-neuron outputs and
  ascending-neuron inputs are cut at the neck, as they are in FlyWire and hemibrain.
- In the whole-CNS sensitivity graph, VNC neurons propagate activity but are never read out.
  Any primary-vs-sensitivity difference is reported separately for the descending (38) and
  ascending (24) strata, because those strata are where VNC loops can act.

### A3.4 Transmitter sign, per dataset (replaces the sign sentence in A2.1)

**Transmitter → sign.** Acetylcholine, dopamine, serotonin and octopamine give +1. GABA,
glutamate and histamine give −1. Histamine is inhibitory in Drosophila (HisCl channels), and
A2.1 wrongly assigned it +1. The male CNS predictor has a histamine class and FlyWire's does
not, so without a shared rule the same types would get opposite signs across datasets.
Drawn targets affected: R1-6, R7, HBeyelet.

**Tiers** (applied per neuron; the first tier that gives a transmitter wins):
0. **Known transmitter (type-level, shared across datasets).** If the unit has a known
   transmitter, it is used for that unit's neurons in all three datasets. The known
   transmitter comes from the majority of FlyWire `known_nt` among member neurons, or, if that
   is absent, from the majority of male CNS `ground_truth`. Entries naming two different
   classical transmitters, and "-negative" entries, are ignored. If FlyWire and male CNS
   disagree, tier 0 is skipped.
1. **Neuron-level prediction passing the dataset's own QC.** Hemibrain: `top_nt` where
   `used_in_analysis` (unchanged). FlyWire: `top_nt` where `top_nt_conf` ≥ 0.5 (82.5% of
   neurons). Male CNS: `predicted_nt` where `predicted_nt_confidence` ≥ 0.5. The two confidence
   scales are not calibrated to each other; 0.5 is a floor, not an equivalence.
2. **Type majority** of tier-1 calls in the same dataset (unchanged).
3. **Type majority** of any named prediction in the same dataset, flagged low-confidence (unchanged).
4. No sign: excluded from propagation (s = 0), as in A2.1.

**Sensitivity:** predictions only (tiers 1–4, no tier 0). For FlyWire this is close to
Shiu et al.'s per-neuron convention. The Week-5 checkpoint used Shiu's `Excitatory` column
and is unaffected; A3.4 applies to all new runs of both models.

## Implementation notes (append-only; not amendments)

- **2026-09-28, first linear run (commit 687afaf), recorded from `manifest.json` seed counts only;
  no correlations computed or viewed.** Target HBeyelet (unit 4656, Tier B) has 4 left and 2 right
  FlyWire neurons in annotations v3.0.0, but the 2 right neurons have no synapses in the Buhmann
  table, so they are not graph nodes. R0 eligibility (A2.4) was computed from annotations, not graph
  nodes. HBeyelet therefore has an empty right-side seed. Per A2.3, R0 and R2-right are NA for
  this target and counted as such: R0 covers 133 of 134 targets, and the R2 right-side sensitivity
  covers 57 of 58 Tier B targets. R1 and R2 primary are unaffected.
- Unit 4383 (DNge153, not a target) has no FlyWire graph node, so the FlyWire R1 readout has
  7,578 of 7,579 units. It drops out of the shared-unit set for any comparison involving FlyWire.

## Amendment 4 — 2026-09-28 (before any correlation is computed)

**Status at time of amendment.** Linear-model response vectors for the primary condition
(α = 0.5, K = 5) have been written to `data/derived/linear/a0.5_K5/` (commit 687afaf).
Only `manifest.json` seed counts have been read (see Implementation notes). No response
values, correlations or summary statistics have been computed or viewed.

**Reason.** H1 names a "label-shuffled permutation baseline" without defining the shuffle,
and the analysis steps between response vectors and H1 are not fully specified.

### A4.1 Comparison per target and rung

For target i on a rung with datasets (A, B):
1. Take i's response vectors in A and B (R1: unit rows; R2: R2-unit rows; R0: rows
   indexed (unit, ipsi | contra), where ipsi means the stimulated side, left-stim vs. right-stim).
2. Keep rows present and non-missing in both.
3. Remove the seed unit(s) of every perturbation being compared. For R0, both sides are
   removed. For R2, the seed R2 unit is removed.
4. Remove rows that are zero in both (A2.3). If fewer than 10 rows remain, the result is NA.
5. ρ = Spearman correlation, with ties given average ranks.

### A4.2 H1 test: mismatched-target null (primary)

- **Statistic:** the median of ρ(i_A, i_B) over the rung's targets. Only drawn targets
  count (R0 133, R1 134, R2 58). Case-study perturbations never enter.
- **Null:** 1,000 random derangements π of the targets. Each gives the median of
  ρ(i_A, π(i)_B), computed with A4.1. Step 3 removes the seed units of both i and π(i).
- **p** = (1 + #{null median ≥ observed median}) / 1001, one-sided.
- The null asks whether a target's response in one dataset matches the same target in the
  other dataset better than it matches a different target. A cell-label shuffle only asks
  whether the vectors share any structure at all, which hub cell types guarantee.
- **Secondary (sanity check):** a cell-label shuffle of B's vector within the retained rows,
  1,000 times per target. Reported, not tested.
- **Multiplicity:** H1 is tested per rung (R0, R1, R2) and per model (linear, LIF), 6 tests,
  with a Holm correction at α = 0.05.
- **Seed:** permutations use numpy default_rng(20260928).

### A4.3 Outputs

- `results/linear/<condition>/rho_<rung>.csv`: one row per perturbation (draw and case study),
  with ρ, the number of rows used, and an NA reason.
- `results/linear/<condition>/h1_<rung>.json`: observed median, the null quantiles
  (2.5, 50, 97.5), p, and the median of the cell-label-shuffle null.
- H2–H4 comparisons are specified before their runs in a later amendment. H4 already
  names "the same permutation baseline as H1" (A1.1), which now means A4.2.

## Amendment 5 — 2026-09-28 (after H1 medians were seen; before any H4 test was run)

**Status at time of amendment.** The primary linear condition (α = 0.5, K = 5, rectified)
has been analyzed under Amendment 4. Results seen: per-target ρ and H1 summaries for R0, R1,
R2 and R2right. Medians: R0 0.52, R1 0.65, R2 0.59, R2right 0.55; p = 0.001 on every rung.
Also seen: the case-study R1 ρ values. **Everything below was specified after those H1
medians were seen and before any H4 test was run.** No R0ipsi, R1ipsi, R1half or unrectified
output exists yet.

### A5.1 H1 scope: inhibitory seeds produce no response (a model property)

- The pre-registered H1 population is all drawn targets. Every target whose seed unit is
  inhibitory returned zero response everywhere: 27 of 134 on R0 and R1, and 10 of 58 on R2.
  With zero baseline activity and rectification at every hop (A2.1), an inhibitory
  perturbation has nothing to suppress. These targets are NA under A2.3.
- H1 as computed therefore covers excitatory-seed targets: R0 107, R1 106, R2 41. It is
  reported that way, with the n and the reason stated. The draw is **not** re-drawn or
  replaced; its value is that it was random.
- This is reported as a finding about zero-baseline models. The Shiu LIF model shares the
  property (zero basal firing, A1.5).
- Other NAs are recorded, not changed. `constant_vector` means a response that is identically
  zero in one dataset, mostly descending-neuron targets whose outputs lie outside the
  hemibrain volume. `lt10_rows` is defined in A2.3.

### A5.2 Sensitivity model: unrectified linear cascade

- This is the A2.1 model with rectification removed. Negative values are kept as suppression.
  α, K, weights and signs are unchanged. With |W| row sums ≤ 1 and α = 0.5 the cascade is stable.
- It runs on **all** targets and all rungs, including the A5.3 forms. H1 and H4 are
  repeated on it and reported next to the primary model. It does not replace the primary
  model: inhibitory targets become measurable here, and the excitatory set is compared
  between the two models.
- Spearman ρ, the double-zero rule and the nulls are as in A4.

### A5.3 Like-for-like rung forms for H4, a method control, and the H4 test

**Why.** R0 as defined in A2.2 stimulates one side and splits the readout into ipsi and contra
rows. R1 and R2 stimulate both sides and read out whole units. The observed order
(R1 > R2 > R0) could be a construction effect, so it is not interpreted.

**Forms used for H4** (every rung: one-sided stimulation, readout on the stimulated side,
one row per unit):
- **R0 → R0ipsi:** FlyWire left-stim / left-side readout vs. FlyWire right-stim / right-side
  readout, per unit.
- **R1 → R1ipsi_right:** FlyWire right-stim / right-side readout vs. male CNS right-stim /
  right-side readout. The right side is chosen to match R2. Male CNS side = `somaSide`, else
  `rootSide` (side known for 99.6% of matched neurons). Sensitivity: **R1ipsi_left**.
- **R2 → R2right** (already pre-registered, A2.2): FlyWire right-stim / right-side readout
  vs. hemibrain. Hemibrain side labels are unreliable (28% of mapped bodies have no side
  suffix), so hemibrain uses its whole volume, which is mostly right hemisphere.

**Method control (reported, not tested): R1half.** FlyWire vs. male CNS, both built exactly
like the original R0: left-stim, readout split into ipsi and contra. If R1half ≈ R0, the
original R0 median reflects the construction. If R1half ≈ R1, the within-animal agreement is
genuinely lower. An R2half is not run, because hemibrain side labels are unreliable (above).

**H4 test** (`src/h4_test.py`):
- Per drawn target, excess ρ = ρ − that rung form's H1 null median. Each rung is scored against
  its own mismatched-target baseline, which is how "the same permutation baseline as H1"
  (A1.1) is implemented.
- Pairs R0 > R2, R2 > R1, R0 > R1. For each pair, use targets non-NA on both, a one-sided
  Wilcoxon signed-rank test on paired excess ρ, and Holm correction over the 3 pairs within
  each model.
- Report the median paired difference with a 95% bootstrap CI (1,000 resamples,
  `default_rng(20260928)`) and the raw ρ medians.
- "Ordering supported" requires both R0 > R2 and R2 > R1 to be significant after Holm.
  Any other outcome, including a reversed order, is reported as found.

**R2 coverage.** Rows are already restricted to units present in both datasets (A4.1). What
remains is truncation of arbors at the hemibrain volume edge, which lowers R2 for a coverage
reason. It is stated as a limitation, not corrected.

### A5.4 Case-study correction: VP4_dry_PN

`VP4_vPN` is GABAergic (tier 0), a selection error in A3.1. It was found because its R1 result
was NA (zero response). It is replaced by `VP4+VL1_l2PN`, the only cholinergic PN type
innervating VP4. That type also innervates VL1, so it is not a pure dry-air channel; this is
stated. No other case-study result was used to make this choice.

### A5.5 Open item for the Arm-2 protocol (A1.8)

In the rectified model, flipping a seed neuron's sign to inhibitory zeroes the whole response
(drift = 1), and flipping a downstream inhibitory neuron to excitatory reveals activity that was
invisible. The NT-mislabel curve would therefore partly measure the rectifier. A1.8 must choose,
before any injector exists, between (a) excluding seed-unit flips from drift, and (b) running the
NT-mislabel condition on the unrectified model (A5.2).
- **2026-09-28, Amendment-5 run (commit 5f91290).** New empty seeds, read from `manifest.json`:
  DNpe020 (unit 4469) has no male CNS neurons with a side label (`somaSide`/`rootSide` both
  missing), so it is NA on R1ipsi_left, R1ipsi_right and R1half. Case-study VP2_hot_PN has a
  single male CNS neuron (left side), so it is NA on R1ipsi_right. R0ipsi inherits the HBeyelet
  NA. The RuntimeWarning "Mean of empty slice" printed during the run comes from targets whose
  label-shuffle null is all NaN (constant vectors). It is harmless and does not affect any ρ or p.

## Amendment 6 — 2026-09-28 (before any H2 or H3 run)

**Status at time of amendment.** Linear results for H1 (all rung forms) and H4, rectified and
unrectified, have been computed and seen (commit 148f0dd). No edge-filtered or sign-ablated
model has been run. The only numbers from those conditions are the edge-survival fractions
below, which come from graph structure alone.

**Reason.** H2 and H3 are named in the Metric section, but the exact procedures are not fixed.

### A6.1 H2 — edge filter

- **Filter:** in each dataset's scoped graph (A2.7), keep a neuron-to-neuron edge i→j if
  n_ij > 10 synapses or n_ij ≥ 1% of j's total input. The filter is applied to each dataset
  independently.
- The Metric section calls this "conserved-edge filtered". The filter is the stated threshold
  only. It does not additionally require the edge to exist in the other dataset, because that
  would make each dataset's graph depend on which rung it is being compared on.
- **Weights:** surviving edges keep their unfiltered input fraction (n_ij ÷ j's unfiltered total
  input). Filtering only removes edges; it does not renormalize the rest.
- **Survival (structure only):** FlyWire keeps 22.3% of edges and 58.6% of synapses; male CNS
  keeps 19.2% and 60.3%; hemibrain keeps 18.9% and 59.6%.
- **Test:** for each primary rung (R0, R1, R2, in the A2.2 forms used for H1), compute
  excess ρ = ρ − that condition's H1 null median for the rung. Pair drawn targets that are
  non-NA in both conditions. Run a one-sided Wilcoxon signed-rank test with the filtered
  condition predicted higher. Apply Holm over the 3 rungs, within each model.
- **Relation to Currier & Clandinin (2025):** their result predicts a positive H2. That
  prediction is noted here, before the run.

### A6.2 H3 — sign ablation

- **Ablation:** every neuron is set excitatory (s_i = +1), including neurons with no assigned
  sign (tier 4). Everything else is as in the base model.
- In the rectified model, inhibitory-seed targets respond under ablation but not in the base
  model. So H3 is tested only on targets non-NA in both conditions, which means excitatory
  seeds for the rectified model and all targets for the unrectified one.
- **Test:** same as A6.1, with the base condition predicted higher. One-sided Wilcoxon on paired
  excess ρ, Holm over the 3 rungs within each model.
- As A1.7 states, H3 is a within-model sanity check. The Arm-2 NT-mislabel condition is
  separate.

### A6.3 Models, outputs, and α/K sensitivity

- H2 and H3 are run for the rectified primary model and, as a sensitivity check, the
  unrectified model (A5.2). LIF follows once its Arm-1 scope is fixed (open item below).
- For each hypothesis, report the median paired difference in excess ρ with a 95% bootstrap CI
  (1,000 resamples, `default_rng(20260928)`), the raw ρ medians, the one-sided p, and the Holm p.
  Output: `results/linear/<base>/h2_h3.csv` (`src/h23_test.py`).
- The α/K sensitivity conditions from A2.1 (α ∈ {0.3, 0.7} at K = 5; K ∈ {3, 10} at α = 0.5) are run
  for H1 and H4 and reported descriptively. The primary conclusions rest on α = 0.5, K = 5.

### A6.4 Open item: LIF Arm-1 scope

- The Metric section runs LIF on every target, and the H1 Holm correction (A4.2) counts LIF as
  the second model. A3.1 makes the case-study set the fixed LIF list for Arm 2.
- Whether LIF Arm-1 runs on all 134 drawn targets or on a smaller fixed subset must be decided,
  with a runtime estimate from the re-run Shiu checkpoint, before any LIF cross-dataset run.
- **2026-09-29, H2/H3 results (Amendment 6; commits e331490 and next).** Both hypotheses not
  supported, in both models; both reversed. Median paired difference in excess ρ (95% CI):
  - H2 (edge filter), rectified: R0 −0.012 (−0.041, 0.010); R1 −0.075 (−0.099, −0.066);
    R2 −0.093 (−0.144, −0.056). Unrectified: R0 −0.099 (−0.109, −0.081); R1 −0.100 (−0.119,
    −0.088); R2 −0.102 (−0.135, −0.073).
  - H3 (sign ablation), rectified: R0 −0.253; R1 −0.229; R2 −0.152 (ablated higher; all CIs
    exclude 0). Unrectified: R0 −0.240; R1 −0.269; R2 −0.313. Raw median ρ under ablation
    0.91–0.96 on every rung in both models.
  - The Currier & Clandinin–based prediction for H2 (A6.1) did not hold.
- **Environment.** On 2026-09-28 at about 16:05, installing the unused `connectome-drift`
  skeleton into the shared venv downgraded pandas 3.0.6 → 2.3.3. Every result from commit
  05034b3 onward used pandas 2.3.3. The matched units and target draw (2026-09-21/22) were
  built with 3.0.6. The skeleton package was uninstalled on 2026-09-29 and exact package
  versions are pinned in `requirements-lock.txt`. No numerical effect is expected; this is
  recorded for reproducibility.

## Exploratory analyses (not pre-registered; planned 2026-09-29 after H2/H3 were seen)

- **E1, renormalized edge filter.** H2 filtering removes about 40% of each neuron's input
  without renormalizing (A6.1). E1 repeats H2 with input fractions recomputed on the filtered
  graph, to separate "which edges" from "how much total drive".
- **E2, hub tracking under sign ablation.** Tests whether all-excitatory responses mainly follow
  structural hub connectivity (response vs. in-degree and vs. the leading eigenvector of the
  weight matrix), as an explanation for the H3 reversal.
- Both are labelled exploratory in all outputs and in the paper.
- **2026-09-29, α/K sensitivity (A6.3; descriptive).** Rectified model, α ∈ {0.3, 0.7} at K = 5 and
  K ∈ {3, 10} at α = 0.5, rung forms R0, R1, R2, R2right, R0ipsi, R1ipsi_right. Table:
  `results/linear/sensitivity_alphaK.csv` (`src/summarize_alphak.py`).
  - H1: p = 0.001 on every rung form at every grid point. Median ρ rises with K (R1: 0.59 at K = 3,
    0.65 at K = 5, 0.66 at K = 10) and moves little with α (R1: 0.64 at α = 0.3, 0.66 at α = 0.7).
  - H4: the same pattern as the primary condition at all four grid points. R0 > R2 and R0 > R1 are
    significant after Holm. R2 > R1 is reversed, with median paired excess difference −0.014 to
    −0.049 and every CI below 0. The observed order R0 > R1 > R2 does not depend on α or K.
  - At K = 3 the mismatched-target null medians are near or below 0 on the one-sided forms
    (R0ipsi −0.071, R1ipsi_right −0.025), compared with 0.07–0.08 at K = 5.
- **2026-09-29, compute environment.** The runs from this date onward were computed in a Linux
  x86_64 environment with the `requirements-lock.txt` versions of numpy, pandas, scipy and pyarrow
  (Python 3.12.3), not on the laptop. As a check, the primary condition (α = 0.5, K = 5) was re-run
  there first. R1/R2 response vectors match commit 148f0dd's to ≤ 6e-17, and `rho_R1.csv` and
  `h1_R1.json` are identical.
- **2026-09-29, E1 results (exploratory).** Code: `run_linear.py --edge-filter-renorm` (condition
  suffix `_efrn`) and `src/e1_compare.py`. Output: `results/linear/exploratory/e1.csv`. Median
  paired difference in excess ρ (95% bootstrap CI), with two-sided Wilcoxon, uncorrected:
  - Renormalized filter vs. base, rectified: R0 −0.029 (−0.052, 0.018); R1 −0.087 (−0.098, −0.074);
    R2 −0.099 (−0.139, −0.055). Unrectified: R0 −0.082; R1 −0.092; R2 −0.084, all CIs below 0.
  - Renormalized vs. unrenormalized filter, rectified: R0 −0.002; R1 −0.006; R2 −0.002.
    Unrectified: +0.019, +0.012, +0.006.
  - Reading: renormalizing restores total drive but does not recover agreement (≤ 0.02 of the
    ~0.08–0.10 H2 drop). So the H2 reversal comes from which edges the filter removes, not from lost
    drive. The weak edges removed by the filter carry signal that is consistent across datasets.
- **2026-09-29, E2 results (exploratory).** Code: `src/e2_hubs.py`. Outputs:
  `results/linear/exploratory/e2_*.csv`. Specification details fixed when the code was written:
  - It runs on rung R1 only. Rows are the A4.1 rows for each target.
  - The "leading eigenvector" is taken from the synapse-count matrix. The input-fraction operator is
    row-stochastic when every neuron is excitatory, so its leading eigenvector is flat.
  - A seed-independent "generic" profile is added: the leave-one-out mean response over the other
    drawn targets. So is a partial ρ that controls for it.
  - Results, median over targets, rectified base → sign-ablated:
    - Response vs. in-degree: FlyWire 0.08 → 0.08, male CNS 0.09 → 0.11.
    - Response vs. eigenvector centrality: −0.03 → −0.11 (FlyWire), −0.05 → −0.10 (male CNS).
    - Response vs. generic profile: 0.13 → 0.18 (FlyWire), 0.20 → 0.20 (male CNS).
    - Cross-dataset ρ 0.65 → 0.95. Partial ρ (generic profile removed) 0.64 → 0.95.
    - The unrectified base gives 0.46 → 0.95. The sign-ablated responses are identical in the two
      models, as expected: there is nothing negative to rectify.
    - Row counts are nearly equal across conditions (median 7,393 vs. 7,577), so the rise is not a
      row-set effect.
  - Reading: the hub explanation is not supported. All-excitatory responses do not follow in-degree,
    eigenvector centrality or a shared generic profile. They stay target-specific, and their
    cross-dataset agreement is about as high after the generic profile is removed. The H3 reversal
    therefore comes from the signs. The next exploratory question is whether it comes from
    cancellation under mixed signs, from sign disagreements between datasets for matched types, or
    both. That is not run yet.
- **2026-09-29, Week-5 checkpoint rerun on v783 (A3.4 signs) and LIF runtime (input to A6.4).**
  Code: `src/shiu_checkpoint_783.py`. Output: `results/checkpoints/shiu_783/<condition>/`.
  - **Setup.** The Sep 21 checkpoint used Shiu's v630 files (`2023_03_23_*_630_final`, see
    `example_executed.ipynb` cell 0), not v783. This rerun keeps the Shiu model unchanged
    (w_syn 0.275 mV, 100 Hz Poisson input, 30 trials × 1 s). Stimulus: Shiu's sugar-GRN list. It has
    20 IDs in v783; 720575940620900446 is not a v783 root. Readout: MN9, 720575940660219265.
  - **MN9 rate** (mean ± SEM over 30 trials; target 67.0):
    - v630, Shiu signs, 21 GRNs (control, same environment): 67.9 ± 1.1 Hz
    - v783, Shiu `Excitatory` column: 60.4 ± 0.8 Hz
    - v783, A3.4 signs: 60.8 ± 1.0 Hz
    - v783, A3.4 predictions-only signs: 62.1 ± 0.8 Hz
  - **Reading.**
    - The environment reproduces the Sep 21 result (67.9 vs. 66.4 Hz).
    - A3.4 signs do not move MN9 relative to Shiu's signs on the same graph. A3.4 flips the sign of
      3.0% of synapses (14,658 presynaptic neurons) and zeroes 0.5% (769 neurons with no sign).
    - The drop of about 7 Hz (about 11%) comes from v630 → v783. That change includes one fewer
      stimulated GRN, which is not separated here.
    - The v783 baseline for all LIF work is therefore about 61 Hz, not 67. This is a replicability
      finding in its own right: the same model on two proofreading snapshots of the same animal.
      It is also a data point for the Arm-2 real-error condition (A1.9).
  - **Runtime.**
    - Single core, after the first compile: about 10 s per 1-s trial on v783. A 30-trial run takes
      about 5 min, including network build.
    - Peak memory is about 4 GB per worker. Two parallel workers exceeded 8 GB, and trials in one
      process leak memory, so the script runs 10-trial chunks in fresh subprocesses.
    - The Sep 21 laptop run took 311 s for 30 trials with all cores. Taking about 5 min per
      (perturbation × dataset or side) run on the laptop:
      - All drawn targets on R0, R1 and R2, which is the H1 LIF family in A4.2: about 650 runs, about
        55 h. Adding the H4 one-sided forms brings it to about 790 runs, about 65 h.
      - R1 only, all 134 targets: about 270 runs, about 22 h.
      - Case-study set only, on R0/R1/R2: about 76 runs, about 6–7 h.
    - The male CNS graph has about 1.45× the connections of FlyWire (1.9× the synapses), so its runs will be slower.
    - Ten trials instead of 30 divides every figure by about 3.
    - The A6.4 decision is still open. These are the inputs for it.

## Amendment 7 — 2026-09-29 (before any LIF cross-dataset run; resolves A6.4)

**Status at time of amendment.** No LIF cross-dataset output exists. The only LIF output is the
v783 checkpoint rerun, which is single-dataset (Implementation notes, 2026-09-29). Linear results
for H1–H4, α/K and E1/E2 have been seen.

**Reason.** A6.4 left the LIF Arm-1 scope open. At about 5 min per 30-trial run, the full H1 LIF
family would take about 55 h of laptop compute.

### A7.1 Scope

- LIF runs on the case-study set only (A3.1). Perturbations per rung: R0 14, R1 14, R2 10.
- That is 76 runs at primary Wsyn:
  - R0: 28 (left and right stimulation)
  - R1: 28 (FlyWire and male CNS)
  - R2: 20 (FlyWire and hemibrain)
- The Wsyn sensitivity condition adds 24 runs (male CNS 14, hemibrain 10). FlyWire's Wsyn is
  0.275 mV in both conditions, so its runs are shared.
- There are no LIF runs on the 134-target draw in Arm 1.
- Graphs, seeds and readout units are those of the linear model: the scoped graphs (A2.7), the A2.2
  seed and readout definitions, and the brain-only readout (A3.3).

### A7.2 Model

- The Shiu LIF model (`model.py`) is unchanged except for weights and signs.
- Signs: A3.4. Neurons with no sign have w = 0.
- Edge weight: raw synapse count × Wsyn.
- Primary Wsyn: structural rescaling (A1.4), taken from `results/wsyn.csv`: FlyWire 0.275 mV,
  male CNS 0.189 mV, hemibrain 0.262 mV.
- Sensitivity: 0.275 mV on every connectome (A1.4).
- Stimulation: Poisson input at 150 Hz, Shiu's default, on every seed neuron. 30 trials of 1 s.
  Trials run in 10-trial subprocesses (Implementation notes, 2026-09-29).
- Unit response: the mean firing rate over the unit's neurons and trials. Per-trial rates are stored
  for every neuron.

### A7.3 Readouts and status

- Readouts are A1.6 items 1–3:
  - Spearman ρ of unit response vectors, per A4.1, with the seed unit(s) removed and the
    double-zero rule applied
  - named-target sign and magnitude agreement
  - top-20 rank agreement, with k = 10 and k = 50 as sensitivity checks
- LIF results are descriptive. They are not used for H1–H4.
- **Deviation from the pre-registration, stated as such.** The Metric section and H1 name both
  model classes, and A4.2 applies Holm over 6 tests (3 rungs × 2 models). H1 is now tested for the
  linear model only, so its Holm family is the 3 linear rungs. The raw linear H1 p-values are 0.001
  on every rung, so Holm-adjusted p = 0.003 each. The conclusion does not change.
- Case-study perturbations have no mismatched-target null of their own. The A4.2 null uses drawn
  targets only. LIF ρ values are therefore reported raw, next to the linear ρ for the same
  perturbation and rung.

### A7.4 Trial-count check (input to A1.8)

- For every run, the Spearman ρ between the unit response vector from the first 10 trials and the
  vector from all 30 is reported.
- For every perturbation and rung, the cross-dataset ρ computed from the first 10 trials is
  compared with the ρ from all 30. Report the median and maximum absolute difference.
- This informs the Arm-2 trial count. No threshold is fixed here.
- **2026-09-29, LIF implementation (Amendment 7), recorded before any LIF cross-dataset ρ was
  computed.**
  - **Code.** `src/lif_prep.py` converts each scoped graph (A2.7) into Shiu's input format, with
    A3.4 signs, and builds the A2.2 seed sets. `src/run_lif.py` runs the simulations and is
    resumable. `src/analyze_lif.py` computes the readouts.
  - **Inputs.** The FlyWire graph matches Shiu's `Connectivity_783.parquet` in edge count
    (15,091,983), synapse total (54,492,922) and node count (138,639). None of the three graphs
    applies a minimum-synapse threshold; PREREG sets none.
  - **Runs.** Identical seed sets are simulated once (9 FlyWire R2 seeds equal their R1 seeds). That
    gives 67 distinct runs at primary Wsyn (FlyWire 43, male CNS 14, hemibrain 10) plus 24 at fixed
    Wsyn.
  - **Trials per process.** With large seed sets, brian2 slows down about 10× after the first
    trial in a process (sugar_GRN, 122 seeds: 23 s for the first trial, then about 265 s per trial)
    and keeps using more memory. Each trial therefore runs in its own subprocess, with 2 in
    parallel on FlyWire and hemibrain. The first hemibrain run (DA1_uPN R2) was made with 10-trial
    processes. It is kept, because the model and trials are the same.
  - **Runtime.** About 5–6 min per 30-trial run with 2 parallel processes. The whole set is roughly
    12–14 h, longer than the checkpoint-based estimate, because whole-type seeds at 150 Hz recruit
    far more activity than 20 GRNs at 100 Hz.
  - **Stimulus weight.** Shiu sets it as w_syn × f_poi. It is held at 0.275 mV × 250 = 68.75 mV in
    every condition by setting f_poi = 250 × 0.275 / Wsyn. Otherwise the Wsyn rescaling (A1.4)
    would also weaken how hard the stimulus drives male CNS seeds. This is an implementation choice
    and was not in A7.2.
  - **Activity level.** Seen from single-dataset activity only, before any cross-dataset
    comparison. Some seeds recruit a large part of the network. Examples at 150 Hz on DA1 PNs
    (both sides):
    - hemibrain: 7,481 of 24,616 nodes active, mean 86 Hz (the protocol run, 30 trials)
    - FlyWire: 5,152 of 138,639 active (1-trial diagnostic)
    - male CNS: 7,898 of 146,176 active (1-trial diagnostic)
    - Sugar GRNs in FlyWire recruit 611 (1-trial diagnostic).
  - **Structural correlate.** The median net excitatory synapse count onto a neuron is 10 in
    FlyWire, 25 in male CNS and 36 in hemibrain; median total input is 200, 336 and 312. Multiplied
    by Wsyn, that is about 2.8, 4.7 and 9.4 mV of median net excitatory drive per unit of presynaptic
    activity.
  - **Handling.** Wsyn is not re-tuned, because A1.4 excludes calibrating on simulation output.
    Activity (fraction of nodes active, spikes per trial) is reported for every run as a
    descriptive LIF result. The diagnostics were written to a scratch folder and are not used.
- **2026-09-30, LIF case-study results (Amendment 7; descriptive, not used for H1–H4).**
  - **Code.** `src/run_lif.py`, `src/analyze_lif.py`. Outputs: `results/lif/<primary|fixed>/`, with
    per-trial spike counts for every run in `results/lif/runs/`.
  - **Changes after the 2026-09-29 note** (implementation only; model and protocol unchanged):
    - **Stimulus input.** From run 31 of 43 on FlyWire onward, every run uses one PoissonGroup plus
      Synapses for all seeds, instead of Shiu's one PoissonInput per seed. The stimulus is the same:
      an independent 150 Hz Poisson train per seed, adding 68.75 mV to v, with no refractory
      period. The reason is speed. brian2 compiles each per-seed PoissonInput separately, which cost
      tens of minutes per new seed set of several hundred neurons.
    - **Equivalence check.** Two FlyWire runs (sugar_GRN R1; a 15-seed run with 5.8k neurons
      active) were repeated with the new input. Unit vectors from the two implementations correlate
      at ρ = 0.999, the same as the split-half ρ within the old implementation (0.999). Total rates
      differ by < 0.1%, and MN9 is 100.5 vs. 100.3 Hz.
    - Run metadata records `poisson` = `group` for new runs. It is absent for runs made with Shiu's
      PoissonInput.
    - FlyWire and male CNS trials run one at a time; two parallel FlyWire trials ran the
      workspace out of memory.
  - **Cross-dataset ρ (A4.1 rows), median [min, max], primary Wsyn:**
    - R0 0.55 [−0.10, 0.78], n = 14
    - R1 0.32 [0.04, 0.54], n = 14
    - R2 0.18 [−0.46, 0.48], n = 10
    - The linear medians for the same perturbations are 0.70, 0.79 and 0.73.
    - LIF agreement is lower than linear for 36 of 38 perturbation × rung pairs. The exceptions are
      sugar_GRN and bitter_GRN on R0.
    - The LIF order is R0 > R1 > R2, as in the linear H4 result.
  - **Fixed Wsyn (0.275 mV everywhere):** R0 is unchanged (FlyWire only), R1 0.22, R2 0.18. Male CNS
    at 0.275 mV is more active: a median of 10.6% of nodes active, vs. 5.8% at 0.189 mV.
  - **Activity (fraction of graph nodes active), median [min, max], primary:**
    - FlyWire 0.8% [0.0, 4.2]
    - male CNS 5.8% [0.1, 7.2]
    - hemibrain 30.3% [0.8, 32.9]
  - **Seed independence.** Responses to the five PN seeds (DA1, DM1, DC1, VA1v uPNs and VP2 hot
    PN) are nearly the same within each dataset. Between different seeds, ρ is 0.96 (FlyWire), 0.97
    (male CNS) and 0.99 (hemibrain), and the Jaccard overlap of active units is 0.90–0.96. Their
    cross-dataset ρ therefore measures agreement of a shared high-activity state, not how each
    seed's own signal propagates. That explains why their ρ is almost identical across seeds
    (R1 0.46–0.48, R2 0.16–0.19).
  - **Named readouts** (A1.6 item 2; R1/R2 need the target active in both datasets or silent in
    both):
    - sugar → MN9: active in FlyWire (85 Hz) and male CNS (41 Hz).
    - LC4 → DNp01, DNp02, DNp11: active in all three datasets.
    - LPLC2 → DNp01: active in all three; DNp02 and DNp11 silent in all three.
    - LC10a → DNa02: active everywhere; hemibrain 3.8 Hz vs. FlyWire 100 Hz on R2.
    - 15 of 16 R1/R2 checks agree. The disagreement is DNp11 under LPLC2 on R2, at 0.1 vs. 0 Hz.
  - **Trial-count check (A7.4).**
    - The unit vector from the first 10 trials correlates with the 30-trial vector at a median
      ρ = 0.999 (minimum 0.89).
    - Cross-dataset ρ from 10 trials differs from the 30-trial value by a median 0.004. The maximum
      difference is 0.021 on R1/R2 and 0.138 on R0 (bitter_GRN, which has few active units).

## Amendment 8 — 2026-09-30 (before any Arm-2 code or output; finalizes A1.8, resolves A5.5)

**Status at time of amendment.** No error injector, perturbed graph, v630 graph or Arm-2 output
exists. Seen so far: linear H1–H4, α/K, E1, E2, the v783 checkpoint, and the LIF case-study
results (Amendment 7). The rung disagreements used as break-even targets are therefore already
known. The error model, grid, estimator and predictions below are fixed before any drift is
computed.

### A8.1 Error types and injection rules (replaces the A1.8 table)

- Injections act on the FlyWire v783 scoped graph (A2.7; Buhmann synapses; no threshold).
- Rates are multiples *m* of a per-type anchor.
- The synapse table has no coordinates. "Co-located" therefore means sharing a neuropil in
  `connections_princeton_no_threshold`.

| Type | Rate unit | Rule | Anchor (m = 1) |
|---|---|---|---|
| Missed synapse | fraction of a neuron's synapses removed | Each synapse is removed independently (binomial thinning of every edge count). | 0.8% |
| False split | fraction of a neuron's synapses detached | Whole partner edges are chosen uniformly and moved to a new untyped fragment node. The fragment stays in the graph and belongs to no unit. | 0.8% |
| False merge | fraction of a neuron's synapses that are foreign | Edges of a random donor that shares the neuron's top neuropil are copied onto it. The donor keeps its edges. | 0.8% |
| Spurious synapse | synapses added, as a fraction of the total | Neuropil drawn ∝ its synapse count; pre ∝ output synapses there; post ∝ input synapses there; +1 synapse. Repeated. | 5% (Buhmann et al. 2021) |
| NT mislabel | fraction of neurons with E↔I flipped | Drawn with probability ∝ (1 − `top_nt_conf`). Neurons with no sign are not eligible. | 6% (Eckstein et al. 2024) |

- **Merge and split anchors.** Precision and recall are set equal to the published proofread
  assigned-synapse F1: median 0.992 (Dorkenwald et al. 2022, Nat Methods, Ext. Data Fig. 6b;
  n = 60 neurons). Each synapse-level anchor is therefore 0.8%. This replaces A1.8's "~1%".
- **Combined.** The segmentation errors are scaled together so that a neuron's assigned-synapse
  F1 against the original graph is 1 − 0.008m:
  - recall loss 0.8%·m, split equally between missed synapses and false splits
  - precision loss 0.8%·m (false merge)
  - spurious synapses at 5%·m and NT flips at 6%·m
- **Lopsided (descriptive).** The same total as the combined condition, but with recall loss to
  precision loss 3:1 (1.2%·m vs. 0.4%·m). This covers the observation that residual errors are
  mostly orphan twigs (Dorkenwald et al. 2024).
- **Unproofread condition (for P6).** Segmentation errors at m_u = (1 − 0.799)/(1 − 0.992) ≈ 25,
  with spurious synapses and NT flips at m = 1. Proofreading changes segmentation only, not
  synapse detection or transmitter prediction.
- **Caps.** Removal, detachment, merge and flip fractions are capped at 0.5; spurious synapses at
  1.0. Single-type grid points above a cap are dropped. In the combined condition, a capped
  component is held at its cap and flagged.
- **Limitations** (recorded here):
  - The neuron sample is small and hand-picked (n = 60).
  - Uniform sampling understates rare, large errors. There is no hub-concentrated condition.
- **Section loss** is dropped from Arm 2. It is named as future work, relevant to serial-section
  EM and expansion microscopy.

### A8.2 NT mislabel and the rectifier (resolves A5.5)

- **Linear primary:** the unrectified model (A5.2), compared against the unrectified rung
  disagreements.
- **Sensitivity:** the rectified model, with seed neurons exempt from flips.
- **LIF:** seed-exempt only.
  - The LIF NT row is therefore not like-for-like with the linear NT row. It is footnoted as such.
  - NT break-even rates are reported against their own disagreement targets and are not rescaled.
- The other types use the rectified primary model.

### A8.3 Grid, perturbations, trials

- **Linear:**
  - m ∈ {0.5, 1, 2, 3, 5, 10, 25, 50}, plus the unproofread condition.
  - 5 injection seeds per (type, m), using `default_rng(20260930 + k)`.
  - Perturbations: the 134-target draw on its A2.4 rung sets. The case-study set is reported
    separately.
- **LIF:**
  - Case-study set, FlyWire only.
  - m ∈ {1, 3, 10}, 3 injection seeds, per type. The combined condition only if time allows; it is
    cut first.
  - 10 trials per run (A7.4; 2026-09-30 LIF results note: median |Δρ| from 30 → 10 trials is
    0.004).
  - Noise floor: each perturbation is re-run unperturbed with 3 fresh Poisson seeds.

### A8.4 Drift and break-even

- **Drift** for target *t* on rung R: 1 − ρ(perturbed, unperturbed FlyWire), on the FlyWire side of
  R's own A4.1 rows (same seeds, readout units, seed removal and double-zero rule). It is averaged
  over injection seeds.
- **Rung disagreement:** D_R = median over the same targets of 1 − raw ρ_R.
- **LIF, like for like:** both sides use 10-trial vectors. D_R comes from the first 10 trials of the
  Amendment-7 runs.
- **Break-even** m*_R: the smallest m at which median drift ≥ D_R, interpolated linearly in log m.
  If it is never reached, it is reported as "> max m".
  - 95% CI from 1,000 bootstrap resamples of targets (drift and D_R resampled jointly),
    `default_rng(20260930)`.
  - Reported per type × model × rung, both as m* and as an absolute rate.
  - Secondary readouts: A1.6 items 2–3, and per-target break-even.
- **LIF targets (post hoc):** R2 is not a LIF break-even target.
  - Hemibrain LIF has a median of 30% of nodes active vs. 0.8% in FlyWire, and A1.4 forbids
    retuning Wsyn.
  - R2 LIF drift is reported descriptively.
  - LIF break-even on R0 and R1 uses the 9 seed-specific case-study perturbations as primary. The
    5 PN seeds are reported as a separate shared-state row.

### A8.5 Predictions (confirmatory; rectified linear model; R0 in its A5.3 like-for-like form)

- **P5:** combined-condition drift at m = 1 is below D_R0.
- **P6:** unproofread-condition drift is below D_R0. Failure to support P6 is treated as the
  substantive result.
- **Test:** one-sided. A prediction is supported if the 95% bootstrap lower bound of
  D_R0 − drift is > 0. Holm correction is applied over the two predictions.

### A8.6 Real-error condition (specifies A1.9)

- **v630 graph:** `2023_03_23_connectivity_630_final.parquet`, processed exactly as v783 (A2.7,
  A3.4 signs by mapped type, no threshold).
  - Check: rebuilding v783 through the same code reproduces the Amendment-7 edge and synapse
    totals.
- **Mapping:** each v783 neuron's annotation `supervoxel_id` is looked up at `mat_630` (CAVE).
  - Accepted when the map is one-to-one.
  - A v630 root that receives more than one v783 neuron is counted as a v630 merge, then excluded.
- **Drop rule:** the condition is dropped if fewer than 80% of the neurons in the readout units and
  seed sets map.
- **Outputs:**
  - drift_real on each rung's rows (linear: drawn targets; LIF: case-study set)
  - the equivalent multiple m_eq on the combined and per-type curves
  - the measured per-neuron synapse loss and gain, compared with the injected rates at m_eq
- **2026-10-07, Arm-2 injector (Amendment 8), recorded before any Arm-2 model output.**
  - **Status.** No perturbed-graph model output exists. The only numbers below are structural
    audits of the injected graphs (synapse counts and fractions) and runtimes. No drift has been
    computed.
  - **Code.** `src/inject_errors.py` (injector, grid, audit) and `tests/test_inject_errors.py`
    (36 tests: rates, caps, reproducibility, eligibility, seed exemption, operator equivalence).
    `linear_model.propagate` gains an optional `fix_cols` argument for seed-exempt flips. Existing
    tests pass unchanged.
  - **Choices not fixed by A8.1, fixed here:**
    - **Neuropils.** Per-neuron, per-neuropil synapse counts come from the Princeton table,
      restricted to pairs whose two ends are nodes of the v783 graph (79 neuropils). The top
      neuropil is the one with the most input + output synapses. 2 graph nodes do not appear
      in the table; their merge donors are drawn from all nodes.
    - **False split.** Each edge end is detached independently with probability equal to the
      rate: the pre end for the presynaptic neuron, the post end for the postsynaptic one. The
      expected detached fraction of each neuron's synapses then equals the rate. All of a
      neuron's detached ends go to one fragment node. Fragments carry the parent's sign
      (after any flip) and are never stimulated.
    - **False merge.** Receiver i gets F_i = S_i · r / (1 − r) foreign synapses, so the foreign
      fraction of its own synapses is r. Each of the donor's edge ends is copied with probability
      min(1, F_i / S_donor). If one donor is too small, more donors are drawn until less than half
      a synapse is outstanding. Copies that would be i → i are dropped.
    - **Spurious synapses.** round(rate × 54,492,922) synapses are added. Self-pairs are redrawn.
    - **NT flips.** Eligible: sign ≠ 0 (137,752 nodes), tier-0 neurons included, as A8.1 says.
      The 204 signed nodes with no `top_nt_conf` get the median weight of the eligible set.
      Count = round(rate × eligible), drawn without replacement.
    - **Combined order.** Missed → split → merge → spurious → flip. The split rate is applied to
      the surviving synapses as split / (1 − missed), so the expected recall loss is exactly
      missed + split. Merge copies use the original donor edges, so foreign synapses are never
      thinned or split.
    - **Seeds.** Injection seed k uses `default_rng(20260930 + k)`, k = 0–4. Each component draws
      from its own child stream (`spawn`, order missed, split, merge, spurious, flip). A given k
      therefore draws the same flips, merges and so on whether that component is injected alone
      or within a combined condition.
    - **Seed exemption (A8.2).** Flips are drawn once per (condition, k) for the whole graph. For
      each perturbation, the flipped signs of its own seed neurons (and their fragments) are
      restored for that perturbation only. This is exact and is implemented as a column sign
      change during propagation. It applies to the rectified sensitivity row and, because they
      run on the rectified model, to the flip component of the combined, lopsided and unproofread
      conditions (P5, P6).
  - **Grid (A8.3, after caps).** Single-type points dropped: spurious m = 25, 50; flip m = 10, 25,
    50. Capped and flagged: combined and lopsided m = 10 (flip held at 0.5) and m = 25, 50 (spurious
    at 1.0, flip at 0.5). That leaves 52 conditions × 5 seeds.
  - **Structural audit (k = 0).** Realized rates match the targets:
    - Recall loss at m = 1 is 0.80% (missed), 0.79% (split) and 0.80% (combined). Unproofread
      recall loss is 20.1%.
    - Spurious synapses at m = 1 are +5.0%. Flips at m = 1 are 8,265 neurons (6.0%).
    - Two properties to note:
      - Split errors are many small fragments. At m = 1, 89,956 neurons have a fragment, and the
        median neuron loses 0.42% against a mean of 0.79%, because whole heavy-tailed edges are
        moved. This is the uniform-sampling limitation already recorded in A8.1.
      - Merge copies add foreign synapses to the donor's partners as well as to the receiver.
        Counting both ends of every synapse, the foreign fraction is ≈ 2r/(1 + r): 1.6% at m = 1
        rather than 0.8%. The receiver-side rate, which A8.1 defines, is r. Both rates are
        reported.
    - The largest point, combined m = 50, has 3.1× the original synapses (75.7M stored edges).
  - **Runtime.** On 2 cores, a 50-perturbation batch takes about 11 s at m = 1 and about 55 s at
    combined m = 50, plus about 30 s to load the base graph. The FlyWire side needs about 510
    perturbation columns per injected graph. The R0 and A5.3 one-sided forms share the same seed
    sets, so they are propagated once.
- **2026-10-08, Arm-2 linear analysis (A8.4–A8.5), fixed before any drift was computed or viewed.**
  - **Code.** `src/run_arm2_linear.py` writes FlyWire-side vectors only (float32).
    `src/arm2_analyze.py` computes drift, curves, break-even and P5/P6. Outputs go to
    `results/arm2/linear/`.
  - **Check.** The unperturbed graph, run through the Arm-2 runner, reproduces the Arm-1 FlyWire
    vectors (max relative difference 3.7e-08). Running drift on the base run against itself gives
    0 on every rung.
  - **Choices not fixed by A8.4:**
    - **R0 and R0ipsi.** Both sides of these rungs are FlyWire, so both are perturbed. Drift is
      computed separately for the left-stimulation and right-stimulation vectors on the rung's
      A4.1 rows, and the two are averaged.
    - **Rows.** A4.1 step 2 is applied against the other dataset of the rung: male CNS for R1,
      hemibrain for R2, and the opposite-side FlyWire vector for R0 and R0ipsi. Seed rows are
      removed, and so are rows that are zero in both the perturbed and unperturbed vector.
    - **Abolished responses.** A perturbed vector that is identically zero while the unperturbed
      one is not gets drift = 1 (ρ taken as 0), not NA.
    - **D_R** is read from the Arm-1 A4.1 tables (`results/linear/<cond>/rho_<rung>.csv`), using
      the model-matched condition: unrectified for the NT primary, rectified otherwise. It is
      taken over the same targets as the drift. Recomputing it from the float32 Arm-2 base run
      reproduces the tables to ≤ 1e-3; the differences come from ties among tiny R0 values.
    - **Break-even.** m* is the first crossing of median drift over D_R. If the first grid point
      already crosses, it is reported as "≤ min m"; if no point crosses, as "> max m". In the
      bootstrap, "> max m" is treated as ∞ when computing CIs. The absolute rate is m* × the
      family's anchor (combined and lopsided: 0.8%, the recall/precision anchor).
    - **P5/P6 test.** The statistic is median(1 − ρ_R0ipsi) − median(drift) over the same drawn
      targets. One-sided bootstrap p = (1 + #{resampled statistic ≤ 0}) / 1001; Holm is applied
      over the two predictions, and a prediction is supported if Holm p < 0.05. The two-sided
      95% CI and the one-sided 95% lower bound are also reported. With one test, this is the
      same as A8.5's "lower bound > 0"; Holm makes the threshold stricter for the larger p.
  - **Order of events.** This note was approved at 2026-10-08 14:17 ET, before any drift had been
    computed. The computer holding this file was offline at that moment, so the note was appended
    after P5/P6 and the partial break-even table had been computed (about 14:30 ET). The text above
    is unchanged from the approved draft.
- **2026-10-08, Arm-2 linear results (Amendment 8; A8.4–A8.5).** Code: `src/run_arm2_linear.py`,
  `src/arm2_analyze.py`. Outputs: `results/arm2/linear/`. Full grid: 287 runs, which is 52 conditions ×
  5 injection seeds plus the two unperturbed base runs, with no failures. Computed in the Linux
  lock-version environment.
  - **Predictions** (rectified model, R0 in its A5.3 form R0ipsi, n = 107 drawn targets;
    D_R0ipsi = 0.476):
    - **P5 supported.** Combined condition, m = 1: median drift 0.294. D_R0 − drift = 0.181 (95% CI
      0.147–0.213; one-sided 95% lower bound 0.154). Holm p = 0.002.
    - **P6 not supported.** Unproofread condition: median drift 0.612. D_R0 − drift = −0.136 (95% CI
      −0.156 to −0.100). Holm p = 1.0. Per A8.5, this is the substantive result. Segmentation at
      unproofread accuracy, with spurious synapses and NT flips at their anchors, changes the model's
      predictions more than the left–right difference within one fly does.
  - **Break-even m\* [95% CI]** (absolute rate in parentheses; "> max" means no crossing up to the
    largest grid point that survived the caps):
    - Combined: R0 2.50 [2.03, 3.27] (2.0% synapse-level error); R1 2.18 [1.79, 2.85] (1.7%); R2 7.80
      [5.11, 16.1] (6.2%).
    - Lopsided (descriptive): R0 2.63, R1 2.22, R2 8.11, close to the combined values.
    - NT flips, unrectified (NT primary, compared with the unrectified D_R): R0 3.24 [3.04, 3.46]
      (19% of neurons); R1 3.50 [3.22, 3.83] (21%); R2 > 5.
    - NT flips, rectified with seeds exempt (sensitivity): > 5 on every rung.
    - Spurious synapses: R0 6.16 [4.59, 9.00] (31%); R1 4.49 [3.61, 6.41] (22%); R2 > 10.
    - False merges: R0 18.4 [16.1, 25.1] (15% foreign on the receiver side); R1 16.1 [13.5, 21.4]
      (13%); R2 47.1 [35.5, > 50].
    - Missed synapses and false splits: > 50 on every rung, that is, still below disagreement at 40%
      loss.
  - **Median drift on R0 at m = 1, by type:** missed 0.004, split 0.025, merge 0.073, spurious 0.194,
    NT unrectified 0.151, NT rectified with seeds exempt 0.103, combined 0.294.
    - The model is nearly insensitive to losing synapses and sensitive to gaining them.
    - In the combined condition, spurious synapses and NT flips dominate.
    - R2 break-evens are the largest because D_R2 is large relative to drift on its 41 targets. R2
      comparisons are also limited by hemibrain truncation (A5.3).
  - **Merge rates.** Counting both ends of every synapse, the foreign fraction is about 2r/(1 + r);
    see the 2026-10-07 implementation note. The receiver-side rates are given above.
- **2026-10-08, exploratory check of the spurious-synapse effect (not pre-registered; seen after
  the results above).** The question was whether spurious-synapse drift is a Spearman artifact:
  units that were silent before and now respond weakly would enter the rank correlation.
  - Rung R1, injection seed 0. Drift was recomputed only on rows that are nonzero in both the
    perturbed and the unperturbed vector.
  - spur m = 1: 0.161 on all A4.1 rows vs. 0.117 on rows nonzero in both. On median, 4% of rows
    become newly nonzero, and they carry 0.3% of the response mass.
  - Combined m = 1: 0.226 vs. 0.179. Merge m = 1: 0.056 vs. 0.041.
  - Reading: about a quarter of the spurious drift comes from newly active rows. The rest is
    reordering among units already active. It is a real sensitivity of rank agreement to added
    weak inputs, consistent with E1 (weak edges carry shared signal), and is not mainly a row-set
    artifact.
- **2026-10-10, Arm-2 real-error condition (A8.6), choices fixed before any v630 model output.**
  - **Mapping.** `src/map_v630.py` was run on the laptop, which has CAVE access, on 2026-10-09.
    Every v783 annotation `supervoxel_id` was looked up in `flywire_fafb_public` at the
    materialization-630 timestamp (2023-03-21 08:10 UTC).
    - All 139,248 neurons returned a root. 136,682 map one-to-one.
    - 1,055 v630 roots each receive more than one v783 neuron (v630 merges), covering 2,566 v783
      neurons.
    - 127,023 of the one-to-one roots are present in the v630 connectivity table.
  - **Drop rule (A8.6).** The rule is read as pooled over the neurons in the readout units and the
    seed sets. "Mapped" means one-to-one and present in the v630 connectivity table.
    - Readout-unit neurons that are v783 graph nodes: 91.6% of 132,877 map, so the condition is kept.
    - Per drawn target, 132 of 134 seed sets are ≥ 80% mapped (median 100%). The exceptions are R1-6
      (unit 6298; 30% of 8,452 neurons) and R7 (unit 6299; 36% of 1,342). These photoreceptors were
      largely reconstructed after v630. drift_real is reported with and without these two targets.
    - Case-study seed sets are 94–100% mapped.
  - **v630 graph.** Shiu's `2023_03_23_connectivity_630_final.parquet` (Buhmann synapses, no
    threshold) is read through the same code as v783.
    - One-to-one roots are relabelled to their v783 IDs.
    - v630 merges and v630 roots with no v783 counterpart are kept as unlabeled nodes. They are
      never seeded or read out, like split fragments in the injector.
    - v783 neurons absent from the v630 table are isolated nodes: seeding one has no effect.
  - **Signs.**
    - Mapped neurons take their v783 A3.4 sign.
    - A v630 merge takes the sign of its v783 neurons if they all agree, otherwise 0.
    - Unmapped v630 roots take 0, which is A3.4's "no call" (primary). The sensitivity condition
      uses Shiu's v630 `Excitatory` column for these roots.
  - **Check.** The same reader applied to `Connectivity_783.parquet` must reproduce the Amendment-7
    totals: 15,091,983 edges and 54,492,922 synapses.
  - **Models and readouts.** As in the sweep: rectified for drift_real against the combined and
    segmentation curves, unrectified for comparison with the unrectified NT curve. Seeds, readout
    units, rows and the drift definition are those of the 2026-10-08 analysis note. There is a
    single "injection", so no averaging over seeds.
  - **m_eq.** The median drift_real over drawn targets on each rung is placed on the median-drift
    curve of each family by linear interpolation in log m. If it is below the smallest grid point,
    it is reported as "< 0.5"; if above the largest, as "> max".
  - **Measured loss and gain.** For neurons that are one-to-one and present in both graphs, synapse
    ends are matched by mapped partner and pair, the same rule as the injector audit:
    - loss = 1 − matched / v783 ends
    - gain = 1 − matched / v630 ends
    Both are reported pooled and as a per-neuron median, and are compared with the injected rates
    at m_eq.
- **2026-10-10, Arm-2 real-error results (A8.6).** Code: `src/run_v630_linear.py`,
  `src/arm2_v630_analyze.py`. Outputs: `results/arm2/linear/v630_*`.
  - **Check passed.** The reader reproduces the v783 totals: 15,091,983 edges and 54,492,922
    synapses.
  - **v630 graph.** 14,687,178 edges, 52,793,639 synapses and 127,400 nodes.
    - 127,015 of its nodes map one-to-one to v783 graph nodes.
    - 385 are unlabeled: 239 v630 merges and 146 unmapped roots.
    - 11,624 v783 graph nodes have no v630 counterpart.
  - **Measured change, v783 → v630,** on the 127,015 one-to-one neurons, synapse ends matched by
    mapped partner:
    - loss 2.85% pooled (median neuron 1.5%)
    - gain 0.82% pooled (median neuron 0.24%)
    - The edit profile is mostly losses, the opposite of the combined condition, which adds 5%
      spurious synapses and 6% NT flips at m = 1. Neither of those differs between the two
      releases here, since detection and signs are shared.
  - **drift_real,** median over drawn targets [90th percentile, max]:
    - Rectified: R0 0.008 [0.016, 0.032]; R1 0.004 [0.009, 0.023]; R2 0.002 [0.008, 0.019].
    - Unrectified: R0 0.007; R1 0.006; R2 0.004 (max 0.12).
    - For comparison, D_R is 0.35–0.48 (rectified) and 0.54–0.70 (unrectified).
    - Excluding R1-6 and R7 does not change any median. The photoreceptors are inhibitory, so
      they are NA in the rectified model. In the unrectified model, n drops from 133 to 131 on R0.
    - The sensitivity signs (Shiu's `Excitatory` for unmapped roots) change medians by ≤ 0.001.
  - **m_eq.**
    - Below the smallest grid point (< 0.5) on the combined, lopsided, split, merge and spurious
      curves, and on the unrectified NT curve.
    - On the missed-synapse curve: R0 1.9, R1 1.6, R2 2.1. That corresponds to 1.3–1.7% injected
      loss, against 2.85% measured loss.
  - **Reading.** The proofreading between v630 and v783 changed the linear
    model's predictions by about 2% of the left–right disagreement (0.008 vs. 0.48). The synthetic error
    model's main drivers, spurious synapses and NT flips, are absent from a release-to-release
    comparison by construction. This condition therefore bounds segmentation edits on already
    proofread neurons. It does not bound the full error budget.
- **2026-10-14, Arm-2 LIF sweep (A8.3–A8.4): implementation and results.** The implementation
  choices below were made before the LIF drift was computed; the results follow.
  - **Code.** `src/run_arm2_lif.py` (runs), `src/arm2_lif_analyze.py` (drift, break-even). Outputs:
    `results/arm2/lif/`. Computed in the Linux lock-version environment, with Brian2 2.10.1.
  - **Runs.** FlyWire only. Each graph runs the 43 distinct case-study seed sets of the Amendment-7
    FlyWire runs, with 10 trials per run.
    - Noise floor: the unperturbed v783 graph, re-run 3 times.
    - Each error type at m ∈ {1, 3, 10} with injection seeds k = 0–2. NT flips at m = 10 are dropped
      by the cap. The combined condition was not run (A8.3: cut first).
    - Total: 45 graphs, 1,935 runs.
    - Model as in Amendment 7. The setup check reproduced MN9 at 100 Hz for sugar GRNs.
  - **Implementation choices:**
    - **Drift reference.** The unperturbed vector is the first 10 trials of the Amendment-7 FlyWire
      run (A8.4, like for like).
    - **Rows.** Rows are limited to those present in the rung's other dataset (R1: male CNS units;
      R2: hemibrain R2 units). R0 uses the opposite-side FlyWire vector, and its drift is the mean
      of the left- and right-stimulation drifts. The rest is as in the 2026-10-08 linear analysis
      note. Fragment nodes have no unit and are never read out.
    - **Seed exemption (A8.2).** For each run, the flipped signs of its own seed neurons and their
      fragments are restored, by negating their outgoing weights in that run only.
    - **Trials.** One trial per subprocess, as in Amendment 7. Two trials run in parallel; on out of
      memory a trial is re-run alone, and a graph switches to one at a time once half its trials
      fail in parallel. From spurious m = 10 onward, workers read only the three connectivity
      columns `model.py` uses (`Presynaptic_Index`, `Postsynaptic_Index`,
      `Excitatory x Connectivity`), as int32, because that graph did not fit in memory otherwise.
      The values are identical.
    - The sweep was interrupted 5 times (3 workspace restarts, 2 out-of-memory crashes). Finished runs were kept, and
      interrupted runs were re-run from scratch.
  - **Noise floor** (unperturbed reruns vs. Amendment 7), median drift: R0 0.005, R1 0.003, R2 0.003
    (max 0.10, R0).
  - **Median drift, 9 seed-specific perturbations, m = 1 / 3 / 10** (D_R: R0 0.54, R1 0.77, R2 0.75):
    - NT flips: R0 0.27 / 0.70; R1 0.25 / 0.63; R2 0.22 / 0.57
    - Merges: R0 0.03 / 0.12 / 0.35; R1 0.03 / 0.09 / 0.27
    - Splits: R0 0.04 / 0.11 / 0.32; R1 0.05 / 0.08 / 0.24
    - Spurious: R0 0.02 / 0.04 / 0.11; R1 0.01 / 0.04 / 0.10
    - Missed: R0 0.010 / 0.017 / 0.046; R1 0.007 / 0.011 / 0.035
  - **Break-even m\* [95% CI]** (R0 and R1 are targets; R2 is descriptive, A8.4):
    - NT flips, seed-specific: R0 2.0 [1.5, > 3] (12% of neurons); R1 > 3.
    - NT flips, PN shared-state row: R0 1.2 [1.2, 1.2] (7%); R1 1.9 [1.8, 2.2] (12%).
    - Every other type: > 10 on every rung and row (> 3 for NT flips).
  - **Reading.**
    - In LIF only transmitter-sign errors reach the left–right disagreement within the grid.
    - Segmentation errors and spurious synapses stay below a third of D_R even at 10× their anchors.
    - This is the opposite of the linear model on spurious synapses: linear crossed at m ≈ 4.5–6,
      LIF drift is 0.10–0.11 at m = 10.
    - LIF break-even rates rest on 9 perturbations; the CIs are wide and open-ended.
- **2026-10-15, exploratory: transmitter-sign error rates from the Eckstein et al. 2024
  confusion matrices (not pre-registered; designed and run after all Arm-2 results above were
  seen).** Source: Eckstein et al. 2024, Fig. 2A, which gives per-transmitter confusion
  matrices; the published accuracies are empirical on annotated neurons and depend on the
  transmitter. A8.1 and the results above are unchanged; this
  entry adds conditions and does not replace any confirmatory result.
  - **Why.** The A8.1 NT anchor (6% of signed neurons, E↔I) was the 94% neuron-level *label*
    accuracy read as a *sign* error rate. Two things make the sign rate lower:
    - The classifier's main confusion is glutamate → GABA, which keeps the sign under A3.4.
    - Tier-0 neurons (A3.4: type-level known transmitter) do not use the prediction. They are
      86,811 of 138,335 signed v783 neurons (63%).
  - **Anchor data.** `data/anchors/eckstein_fig2a.csv`: the four Fig. 2A matrices (FAFB-Catmaid
    and hemibrain; per synapse and per neuron; neuron split, test set), transcribed from the
    published figure and checked cell by cell. `src/sign_error_rate.py` gives label and sign
    (cross-sign) error per transmitter.
  - **Rates.** Weighted by the transmitter mix of v783 tier 1–3 neurons
    (`data/anchors/composition_v783_predicted.csv`; ACh 30,714, Glu 12,244, GABA 7,642, 5-HT 712,
    DA 180, OA 32):

    | Matrix | Label error | Sign error, tier 1–3 | Sign error, all signed |
    |---|---|---|---|
    | FAFB, per neuron (n = 135) | 3.6% | 2.9% | 1.1% |
    | hemibrain, per neuron (n = 207) | 14.9% | 9.3% | 3.5% |
    | FAFB, per synapse | 10.1% | 7.8% | — |
    | hemibrain, per synapse | 22.7% | 13.0% | — |

    Under FAFB, most sign error sits on glutamate neurons (Glu → 5-HT 0.11). Under hemibrain, it
    sits mostly on the many cholinergic neurons (ACh → Glu/GABA 0.10). The neuron-level 5-HT and
    OA rows have n = 1, so those two rows come from the synapse level. The FAFB neuron-level ACh
    row (n = 68) has no errors, which is optimistic.
  - **Conditions** (`src/run_nt_realistic.py`). Only tier 1–3 neurons are eligible. Each flips
    with the cross-sign probability of its assigned transmitter's row (the row of the predicted
    transmitter, no Bayes inversion).
    - Matrices: FAFB (primary; same imagery as FlyWire) and hemibrain (pessimistic).
    - Draws: per neuron, independent; or per cell type, with one draw per `type` at the rate of
      the type's majority transmitter, applied to every eligible member. Untyped neurons are
      drawn individually. The per-type draw follows Eckstein et al. 2024 Fig. 4D–F and S2:
      errors are shared within types and across datasets.
    - 5 draws each, `default_rng(20261015 + k)`.
    - Models, seeds, rows and drift as A8.2/A8.4 and the 2026-10-08 analysis note: unrectified
      (NT primary) and rectified with seeds exempt.
  - **Check.** In this environment, `flip_m1` re-run through `run_arm2_linear.py` reproduces
    `drift.csv` to 1e-16 on every row.
  - **Results** (`results/arm2/linear/ntcm_summary.csv`). Median drift over drawn targets, as a
    fraction of D_R, for R0 / R1 / R2:

    | Condition | Flipped (of signed) | Unrectified | Rectified, seeds exempt |
    |---|---|---|---|
    | FAFB, per neuron | 1.1% | 0.12 / 0.10 / 0.06 | 0.12 / 0.10 / 0.05 |
    | FAFB, per type | 1.2% | 0.14 / 0.12 / 0.07 | 0.13 / 0.11 / 0.06 |
    | hemibrain, per neuron | 3.5% | 0.42 / 0.34 / 0.21 | 0.31 / 0.22 / 0.11 |
    | hemibrain, per type | 3.4% | 0.45 / 0.40 / 0.25 | 0.34 / 0.27 / 0.18 |

    - Every 95% bootstrap CI of D_R − drift lies above 0. At classifier-implied rates, transmitter
      errors stay below every rung's disagreement in both linear models.
    - Per-type draws add up to 0.03 of drift at the same flip fraction.
    - These flips do more per neuron than the A8.1 rule (weighted by 1 − `top_nt_conf`). On the
      A8.1 flip curve, the FAFB conditions sit at m_eq ≤ 0.7 (below the 0.5 grid point on three
      rows), the drift of up to about 4% flips under A8.1, from about 1.1% actual flips. The
      hemibrain conditions sit at m_eq ≈ 1.3–2.3. Break-even expressed as "% of neurons flipped" therefore depends on which neurons flip.
  - **Combined condition with the NT component corrected** (`src/run_combined_ntfix.py`;
    rectified model, seed-exempt flips). The flip component is 1.1%·m of signed nodes instead of
    6%·m. It is still drawn by the A8.1 rule. Everything else is as A8.1.
    - Break-even m\* [95% CI]: R0 3.77 [2.88, 4.97] (3.0% synapse-level error); R1 2.97 [2.35,
      4.02] (2.4%); R2 > 5, the top of this grid. A8.1 gave 2.50 (2.0%), 2.18 (1.7%) and 7.8.
    - P5 statistic: drift 0.249 vs. D_R0ipsi 0.476. The difference is 0.226 (95% CI 0.195–0.260),
      against 0.181 in A8.5.
    - P6 statistic: drift 0.605. The difference is −0.130 (CI −0.147 to −0.092), against −0.136.
      The P6 result does not depend on the NT anchor: unproofread segmentation alone exceeds
      D_R0.
  - **Not rerun: LIF.** The LIF NT break-even (R0 m\* 2.0, 12% of neurons; PN shared-state row
    1.2, 7%) was computed under the A8.1 rule. On the linear curves, the FAFB conditions match
    A8.1 m ≈ 0.5–0.7 and the hemibrain conditions m ≈ 1.3–2.3. If the LIF model responds to the
    two rules the same way, FAFB-rate errors would sit below LIF break-even on both LIF rows. The
    hemibrain-rate errors would sit near the R0 seed-specific break-even and above the PN row.
    This is an untested mapping, not a LIF result.
  - **Limitations.**
    - The matrices are small (135 and 207 test neurons) and come from cell types with known
      transmitters. Types that have never been characterized may be predicted worse.
    - FAFB ground truth used manually placed synapses, while FlyWire predictions use automatic
      synapses, so FAFB rates may be optimistic.
    - Tier-0 known transmitters are treated as error-free.
    - R0–R2 disagreement cannot detect transmitter errors shared across datasets. Under
      predictions-only signs, a shared misprediction raises cross-dataset agreement. Tier 0
      removes this only for types with a known transmitter.
- **2026-10-15, correction to the entry above.** Its source sentence was reworded to cite only
  Eckstein et al. 2024, Fig. 2A. No analysis, number or condition changed.
- **2026-10-15, exploratory: v783 predictions against known transmitters (check on the
  Fig. 2A rates).** `src/nt_v783_known_check.py`, output `results/arm2/linear/nt_v783_known_check.txt`.
  Compares `top_nt` with `known_nt` (one positively named classical transmitter; histamine
  excluded) on 68,690 v783 neurons of 1,013 known types. Fig. 2A has 135 test neurons.
  - **Sign error, by subset.** Each row is the per-transmitter rate, then that rate weighted by
    the tier 1–3 mix (as in the entry above):

    | Subset | Neurons | Sign error | Weighted, tier 1–3 | Weighted, all signed |
    |---|---|---|---|---|
    | All known | 68,690 | 2.1% | 2.5% | 0.9% |
    | All known, `top_nt_conf` ≥ 0.5 | 59,608 | 0.7% | 1.1% | 0.4% |
    | Central (non-optic) known | 13,410 | 4.5% | 4.0% | 1.5% |
    | Central known, conf ≥ 0.5 | 10,693 | 1.3% | 1.4% | 0.5% |

  - The FAFB-matrix condition (2.9% of tier 1–3; 1.1% of signed) lies inside this range. The
    hemibrain-matrix condition (9.3%; 3.5%) lies above all of it, so it remains a pessimistic
    bound.
  - Large label errors that keep the sign: ACh → DA (Kenyon cells, a known misprediction) and
    Glu → GABA. Small classes are poor: 5-HT 75–80% label error, OA 25–45% (n = 34–79).
  - **Clustering.** Among types with at least 2 neurons and any sign error, 21–50% have every
    member wrong. Real errors lie between the per-neuron and per-type conditions, which bracket
    them.
  - **Caveat.** Known types are not a random sample (80% optic lobe; well-studied central
    types). Tier 1–3 neurons are by definition of unknown type and may be predicted worse.
    The central, unfiltered row (4.0% of tier 1–3) is the closest proxy.
