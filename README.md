# Connectome noise ceiling

**Do connectome reconstruction errors matter more or less for brain emulation than natural
differences between individual brains?**

This repository runs the same perturbations through linear and spiking (Shiu et al. 2024 LIF)
whole-brain models built from three *Drosophila* connectomes: FlyWire v783, hemibrain v1.2.1 and
MaleCNS v1.0. It also compares the left and right hemispheres within FlyWire. Two arms:

- **Arm 1 (noise ceiling).** Cross-animal agreement of predicted responses, on matched cell-type
  units (7,579). Rungs: R0 FlyWire left vs. right, R1 FlyWire vs. MaleCNS, R2 FlyWire vs. hemibrain.
- **Arm 2 (error budget).** Realistic reconstruction errors injected into FlyWire v783 at increasing
  rates, to find the break-even rate at which prediction drift equals the R0/R1/R2 disagreement.
  Real-error condition: FlyWire v630 vs. v783.

This is a test of the Sandberg & Bostrom "noise floor" criterion. The results are about
replicability and sensitivity of model predictions. They are not a validation against neural
activity.

The analysis was pre-registered: [`results/prereg/PREREG.md`](results/prereg/PREREG.md) is the
binding spec, with dated, append-only amendments. Anything not in H1–H4 is labelled exploratory.

## Reproduce the figures and every number in the preprint

Needs Python 3.11+ and the committed result tables only (no connectome downloads, a few seconds).

```bash
git clone https://github.com/willaswitzer/connectome-noise-ceiling.git
cd connectome-noise-ceiling
python -m venv .venv && source .venv/bin/activate
pip install -r requirements-repro.txt

python src/preprint_numbers.py           # recomputes the 28 headline numbers, checks them against the text
jupyter nbconvert --to notebook --execute notebooks/reproduce_figures.ipynb   # Figs 1-4
pytest -q tests                          # unit tests for the linear model, analysis and error injector
```

Figures are written to `results/figures/`.

## Rebuild from the raw connectomes

This step needs the full environment (`requirements-lock.txt`, macOS / Python 3.12) and the
upstream data below, placed under `data/` (git-ignored). The LIF runs also need a checkout of
[philshiu/Drosophila_brain_model](https://github.com/philshiu/Drosophila_brain_model) at
`./Drosophila_brain_model`. That checkout supplies the model code and the FlyWire v630/v783
connectivity files.

| Dataset | Source | Goes to |
|---|---|---|
| FlyWire v783 annotations, connections | [flywire.ai](https://codex.flywire.ai/api/download), [flyconnectome/flywire_annotations](https://github.com/flyconnectome/flywire_annotations) v3.0.0 | `data/flywire_783/` |
| FlyWire v783 / v630 Buhmann connectivity | Shiu et al. repo above | `Drosophila_brain_model/` |
| hemibrain v1.2.1 | neuPrint, via `src/01_fetch_hemibrain.py` (needs `NEUPRINT_TOKEN`) | `data/hemibrain/` |
| MaleCNS v1.0 | `gs://flyem-male-cns/v1.0/connectome-data/flat-connectome/` | `data/male_cns/` |
| v630 ↔ v783 root-ID map | `src/map_v630.py` (needs a CAVE token) | `data/derived/v630_map.csv` |

The pipeline runs in the order given by the PREREG amendments, roughly:

1. units and graphs: `build_units_v2.py`, `build_r2_units.py`, `build_matched_graphs.py`,
   `assign_signs.py`, `compute_wsyn.py`
2. Arm 1: `run_linear.py`, `analyze_linear.py`, `lif_prep.py`, `run_lif.py`, `analyze_lif.py`
3. Arm 2: `inject_errors.py`, `run_arm2_linear.py`, `run_arm2_lif.py`, `arm2_*analyze.py`
4. exploratory: transmitter error model and glutamate sign (`run_nt_realistic.py`,
   `run_glu_sign.py`, `glu_crossing.py`)

Each script's docstring states its inputs, outputs and the PREREG section it implements. The LIF
sweeps are resumable and skip finished runs.

## Layout

```
src/                 analysis code (one script per step)
results/prereg/      PREREG.md (binding spec), target draw, case-study set
results/linear/      Arm 1 linear-model outputs per (alpha, K, condition)
results/lif/         Arm 1 LIF outputs
results/arm2/        Arm 2 drift, break-even and exploratory outputs
results/figures/     preprint figures
data/anchors/        published statistics used by the error model (transcribed, with source)
notebooks/           reproduce_figures.ipynb
tests/               pytest suite
```

## Data and credit

The derived tables in `results/` and `data/anchors/` are released under CC BY 4.0. They are derived
from FlyWire (Dorkenwald et al. 2024; Schlegel et al. 2024), hemibrain (Scheffer et al. 2020),
MaleCNS v1.0 (Berg et al. 2025) and the FlyWire neurotransmitter predictions (Eckstein et al.
2024). Please cite those sources too. `data/anchors/eckstein_fig2a.csv` holds confusion-matrix
values transcribed from Eckstein et al. 2024, Fig. 2A. The LIF model is from Shiu et al. 2024
(MIT licence) and is used unchanged except where `src/run_lif.py` says otherwise.

## Licence and citation

Code: MIT (see `LICENSE`). Cite via `CITATION.cff` (GitHub's "Cite this repository") and the
preprint.
