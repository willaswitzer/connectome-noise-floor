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
