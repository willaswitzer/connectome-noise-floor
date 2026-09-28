# Related work notes

## Currier & Clandinin, Cell 2025 — "Infrequent strong connections constrain connectomic predictions of neuronal function"
bioRxiv 10.1101/2025.03.06.641774; Cell PMID 40460825; data: Dryad 10.5061/dryad.pg4f4qs1j.

**What they did**
- Two-photon imaging (SPARC-L sparse labeling) of visual responses in 43 medulla cell types
  (510 neurons, 54 flies; 2D ternary noise, blue and UV).
- Predicted each type's response from its **direct presynaptic inputs only** (one synaptic step),
  weighting measured input responses by connectome input fraction. No multi-hop simulation,
  no LIF, no whole-brain model.
- Connectome (STAR Methods, "Connectome analyses"): input fractions and the strong-connection
  result (Fig. 6, S8) use the **male optic lobe** (Nern et al. 2024, FAMB v1.0, via the Reiser lab
  cell type explorer). Type-to-type Jaccard distances (Fig. 7) use **FlyWire v783** (Matsliah et
  al. 2024). Imaged flies were **female**. So the core predictions pair female physiology with a
  male wiring diagram; the paper does not treat this cross-animal, cross-sex mismatch as a variable.
- No cross-animal comparison of predictions, no reconstruction-error analysis.
- TO CHECK: whether the Nern optic lobe is the same animal as MaleCNS v1.0 (believed yes).

**Findings**
- Connectivity predicts orientation, direction and color selectivity; it does not predict
  receptive-field size (arbor size vs. RF area R = 0.19).
- "Strong" input = >5% of a cell's input (~15–30 synapses), about the top 10% of connections.
  Strong connections show higher pre/post response correlation.
- A "strong-input-dominant" model beats the connection-weighted sum for 28/43 types (4/43 the
  reverse). Removing the single strongest input reduces prediction quality (28/43 types).
- Weak inputs (<5%) do not shape key computations, individually or collectively.
- Asymmetry: similar physiology predicts similar connectivity, not the reverse.

**Relation to this project**
- Hook for the intro: they pair female physiology with male connectivity without measuring the
  cost; Arm 1 (R1) measures exactly that cost at the level of model predictions.
- Different question: they ask which connections *within one connectome* determine function,
  scored against measured activity. We ask whether a whole-brain model's predictions *replicate
  across animals* and under reconstruction error, with no claim about measured activity.
- Arm 1: H2 (conserved-edge filter, >10 synapses or >1% of input, raises replicability) is in
  the same spirit; cite them there. Their threshold is 5%.
- Arm 2 prediction (candidate for the A1.8 protocol, to be fixed before any Arm-2 output):
  errors that mainly touch weak edges (missed and spurious synapses) should cause little drift;
  errors that alter strong edges (false merges/splits, NT sign flips on strong inputs) should
  dominate. Operationalize by also reporting drift stratified by the input fraction of the
  edges each injected error touches.

**Draft paragraph (preprint, related work)**

Currier and Clandinin (2025) tested connectomic predictions of neuronal function against
measured activity: they imaged the visual responses of 43 medulla cell types in female flies and
compared them with predictions built from each type's direct presynaptic inputs in a male
optic-lobe connectome (Nern et al., 2024). Connectivity predicted
orientation, direction and color selectivity but not receptive-field size, and the predictive
signal came from rare strong connections, inputs above 5% of a cell's total input, while weaker
inputs contributed little individually or collectively. Their question is which connections
within a single connectome determine function. We ask a different one: whether the multi-hop
predictions of the same whole-brain model replicate across connectomes from different animals,
and how much reconstruction error they tolerate relative to that cross-animal variability; we
make no claim about agreement with measured activity. Their design also illustrates why that
question matters: connectome-to-function studies routinely borrow wiring from a different
animal, here of the other sex, without measuring what that substitution costs. Their result nonetheless makes a specific
prediction for our error budget: errors that mainly remove or add weak connections, such as
missed and spurious synapses, should shift model predictions little, whereas errors that alter
strong connections, such as false merges and splits or sign errors on strong inputs, should
dominate the drift.
