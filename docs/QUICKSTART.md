# Quickstart Guide

This guide summarises the minimum steps required to reproduce the ALS perturbation and embedding analysis workflows delivered for the Helical Coding Challenge.

---

## 1. Environment Preparation

- Python 3.10 (3.9 compatible if Scanpy installed via conda).
- GPU with CUDA 11.8+ strongly recommended for GeneFormer embeddings.
- Create an isolated environment:

```bash
python -m venv .venv
# On Windows: .venv\Scripts\activate
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

- Optional: if you have access to Helical's SDK, install the GeneFormer package.

```bash
pip install helical
```

---

## 2. Required Data Assets

1. Download the pre-filtered ALS dataset (`counts_combined_filtered_BA4_sALS_PN.h5ad`) from the Helical S3 bucket referenced in the challenge brief.
2. Place the file under `localtools/` or adjust `INPUT_H5AD` in `main.py` / the notebooks to point to your copy.
3. Verify integrity using `tools/check_geneformer.py` and `tools/check_task3_prerequisites.py` if you need automated sanity checks.

---

## 3. Running the Workflows

### Preferred order

Detailed notes for each notebook are available in `docs/TASK1.md`, `docs/TASK2.md`, and `docs/TASK3.md`.

1. `third attempt/task1_perturbation_workflow.ipynb` - build and validate perturbation helpers.
2. `third attempt/task2_als_perturbations.ipynb` - apply ALS perturbations and embed with GeneFormer.
3. `third attempt/task3_embedding_analysis.ipynb` - analyse latent-space shifts, clustering quality, and rescue metrics.

Each notebook exposes configuration cells at the top for data paths, sampling, and runtime shortcuts.

### Command-line alternative

For automated perturbation generation without notebooks:

```bash
python main.py \
  --input localtools/counts_combined_filtered_BA4_sALS_PN.h5ad \
  --targets SOD1 TARDBP FUS C9orf72 \
  --subset 'Condition == "ALS" and CellClass == "Neuron" and Region == "BA4"' \
  --output als_simulated.h5ad
```

Set `MAKE_PLOTS = False` in `main.py` when running headless or on servers without graphical backends.

---

## 4. Outputs and Interpretation

- Notebook artefacts are stored in `third attempt/outputs/`. Expect figures (UMAPs, clustering comparisons, volcano plots) and metrics tables (`embedding_metrics.csv`).
- Intermediate embeddings (`ALS_embeddings.npy` plus `metadata.csv`) are saved for reuse in external analyses.
- `main.py` writes `als_simulated.h5ad`, containing baseline, knock-up, and knock-down AnnData slices labelled via `obs['perturbation']`.

---

## 5. Troubleshooting

- **GeneFormer import issues:** run `tools/check_geneformer.py` to see which dependency path fails. Use `force_pca=True` with `generate_embeddings` for a deterministic fallback.
- **Large memory usage:** down-sample within notebooks using provided configuration toggles or adjust `sample_size` arguments in helper functions.
- **Missing outputs:** use `tools/check_task3_prerequisites.py` to confirm Task 2 artefacts before executing Task 3.
- **Scanpy plotting on servers:** disable plotting in scripts (`MAKE_PLOTS = False`) or use `sc.settings.set_figure_params(dpi=100, figsize=(6,5))` to reduce resource requirements.

---

## 6. Suggested Next Steps

- Capture final figures for presentation from `third attempt/outputs/images/`.
- Export summary statistics (`embedding_metrics.csv`) alongside notebook HTML exports for review.
- Version-tag working checkpoints once reruns complete to track reproducibility.

