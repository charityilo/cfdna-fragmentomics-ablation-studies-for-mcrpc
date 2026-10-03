# cfdna-fragmentomics-ablation-studies-for-mcrpc
Contains Bioinformatics Projects
# BIOINFORMATICS
# Fragmentomics versus ichorCNA: analysis notebooks

Code for the technical report *Fragmentomic and Copy-Number-Derived Tumour Fraction Contributions to ctDNA Burden Prediction in Metastatic Castration-Resistant Prostate Cancer: An Ablation Study*.

## Contents

| File | Purpose | Report sections |
|---|---|---|
| `data/elife-71569-supp1-v2.xlsx` | Supplementary File 1 (version 2) of Renaud et al. (2022), *eLife* 11:e71569, redistributed under the article's Creative Commons Attribution licence. SHA-256 `f0693eee0b6b2c0ddcdcfd609ba9ab087828db73313a135b05afb4859c4698e6` |
| `common.py` | Data loading, checksum and the 14 fragmentomic features |
| `engine.py` | Models, cross-validation, nested recalibration, metrics and the corrected t-test |
| `Data verification and features Notebook.ipynb` 01 | File verification, reproduction of published values, features, pairing test  | Tables 3.1–3.2, Section 4.1, Figure 4.1 |
| `Ablation, calibration and leakage Notebook.ipynb`02 | Analyses 1–2 and the leakage sensitivity analysis (about 30 minutes on one core) | Tables 4.1–4.4, Figures 4.2–4.4 |
| `Cross-platform transfer and cohort-size stability Notebook.ipynb` 03 | Analyses 3–4 (about 25 minutes on one core) | Tables 4.5–4.6, Figures 4.5–4.6 |
| `robustness.py` | XGBoost and support vector regression (tuned by nested cross-validation, and with default settings), and grouped permutation importance | 3.6, 3.7.5 |
| `Robustness and explainability Notebook.ipynb` o4 | Robustness to the choice of algorithm and to tuning, and feature importance (about 20 minutes on one core; run after 02) | Table 4.7, Figure 4.7 |
| `Feature Engineering Guide.md` | Feature specification written before the analysis | 
| `requirements.txt` | Pinned package versions |

## How to run

```
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter lab
```

Open the notebooks in order (01, then 02, then 03, then 04). In each, choose *Kernel → Restart Kernel and Run All Cells*, so that the saved outputs come from a clean run. Tables are written to `results/` and figures to `figures/`.

The slow steps are cached (`results/ablation.pkl`, `leakage.csv`, `transfer.csv` and `stability.csv`). Re-running a notebook reloads these files instead of repeating the computation; delete a file to force it to be recomputed. If these notebooks replace an earlier copy, results already cached in `results/` are reused, so re-running them only regenerates the tables and figures under the current file names.

## Output files

Files are named after their numbers in the report.

## Random seeds

Cross-validation partitions: 2026. Model fitting and inner cross-validation: the split index (the recalibration's inner folds use the split index + 1). Leakage analysis, random exclusion: 7, re-seeded for each k. Pairing permutation test: 1. Stability subsets: 100 × n + draw (subset) and 1,000 × n + draw (partition).

