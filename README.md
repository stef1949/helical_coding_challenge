# Helical Coding Challenge - Computational Biology

Author: Steph Ritchie (MSc Genetic Manipulation & Molecular Biosciences)  
Project: In-silico gene perturbation and embedding analysis for ALS  
Date: October 2025

---

## Overview

This repository contains the Helical Computational Biology Coding Challenge submission focused on amyotrophic lateral sclerosis (ALS). The work simulates knock-up and knock-down perturbations for ALS-associated genes, embeds the resulting expression profiles with GeneFormer_V2 (gf-12L-95M-i4096), and interprets how those perturbations reorganise cell states in latent space.

Deliverables include:

- Modular perturbation utilities for AnnData objects (dense or sparse).
- GeneFormer-ready preprocessing and embedding workflows.
- Quantitative metrics for clustering quality, disease-to-healthy trajectory shifts, and neighbourhood composition.
- Reusable notebooks that document the full analysis pipeline end-to-end.

---

## ALS Fold-Change Reference Data

The following table summarises **representative ALS fold changes** across multi-omic layers (gene, protein, and metabolite level) from peer-reviewed human studies. These values can be used for benchmarking or validating simulated perturbation magnitudes.

| **Biomarker / Gene / Protein** | **Sample Type** | **log₂ Fold Change (ALS vs Control)** | **Approx. Fold Change** | **Direction** | **Pathway / Relevance** | **Source / DOI** |
|--------------------------------|-----------------|----------------------------------------|--------------------------|----------------|---------------------------|------------------|
| **Neurofilament light (NfL)** | CSF | +2.3 to +3.5 | ↑ ~5–10× | Upregulated | Neuronal damage marker; diagnostic biomarker | [Benatar et al., *Nat. Med.* 2023, doi:10.1038/s41591-021-01295-3](https://doi.org/10.1038/s41591-021-01295-3) |
| **Chitotriosidase (CHIT1)** | CSF | +3.3 | ↑ ~10× | Upregulated | Microglial activation marker | [Varghese et al., *JAMA Neurol.* 2021, doi:10.1001/jamaneurol.2020.3965](https://doi.org/10.1001/jamaneurol.2020.3965) |
| **GFAP** | CSF / Plasma | +0.58 | ↑ ~1.5× | Upregulated | Astrocytosis; glial reactivity | [Huss et al., *Neurology* 2023, doi:10.1212/WNL.0000000000201447](https://doi.org/10.1212/WNL.0000000000201447) |
| **GAD2 (mRNA)** | Spinal cord | −1.78 | ↓ ~0.29× | Downregulated | GABAergic neuron loss | [Zhang et al., *Genes* 2020, 11(4):448, doi:10.3390/genes11040448](https://doi.org/10.3390/genes11040448) |
| **CALB1 (mRNA)** | Spinal cord | −1.96 | ↓ ~0.25× | Downregulated | Motor neuron marker; selective vulnerability | [Zhang et al., *Genes* 2020, 11(4):448, doi:10.3390/genes11040448](https://doi.org/10.3390/genes11040448) |
| **GABRE (mRNA)** | Spinal cord | −1.39 | ↓ ~0.38× | Downregulated | Inhibitory signalling | [Zhang et al., *Genes* 2020, 11(4):448, doi:10.3390/genes11040448](https://doi.org/10.3390/genes11040448) |
| **miR-206** | Serum / Muscle | +1.5 | ↑ ~2.8× | Upregulated | Muscle regeneration and denervation response | [Waller et al., *Neurology* 2017, doi:10.1212/WNL.0b013e318272f45d](https://doi.org/10.1212/WNL.0b013e318272f45d) |
| **miR-199a-5p** | Serum / Plasma | −0.9 | ↓ ~0.53× | Downregulated | Mitochondrial regulation; neuroinflammation | [Freischmidt et al., *Cell Death Discov.* 2020, doi:10.1038/s41420-020-00397-6](https://doi.org/10.1038/s41420-020-00397-6) |
| **Lactate** | CSF / Plasma | +0.9 | ↑ ~1.9× | Upregulated | Mitochondrial dysfunction marker | [Blasco et al., *Neurobiol. Aging* 2018, doi:10.1016/j.neurobiolaging.2018.08.003](https://doi.org/10.1016/j.neurobiolaging.2018.08.003) |

**Summary:**
- **RNA-seq data** (post-mortem or single-cell) shows moderate fold changes (~0.25–0.4× down for neuronal genes; ~2× up for glial genes).
- **Protein biomarkers** (NfL, CHIT1) display the largest fold changes (up to 10×), making them the most robust clinical indicators.
- **Metabolites and miRNAs** exhibit intermediate changes (1.5–3×), reflecting metabolic and compensatory dysregulation.

---

## Repository Layout

| Path | Description |
| --- | --- |
| `analysis.py` | Helper class for analysing latent embeddings (distances, neighbourhoods, clustering metrics). |
| `data_utils.py` | Utilities for ALS gene lists, dataset filtering, and GeneFormer preprocessing. |
| `main.py` | Scriptable entry point that generates simulated perturbations and optional UMAP plots. |
| `perturbation.py` | Higher-level perturbation helpers (batch simulations, validation, combined matrices). |
| `third attempt/task1_perturbation_workflow.ipynb` | Notebook that builds the perturbation pipeline and validates per-gene edits. |
| `third attempt/task2_als_perturbations.ipynb` | Notebook that applies perturbations to ALS genes and embeds them with GeneFormer_V2. |
| `third attempt/task3_embedding_analysis.ipynb` | Notebook that interprets embedding shifts (dimensionality reduction, clustering, rescue scores). |
| `third attempt/outputs/` | Generated figures, metrics tables, and intermediate artefacts saved by the notebooks. |
| `requirements.txt` | Python dependencies known to work for the challenge environment. |

---

## Task Documentation

- [Task 1 - In-Silico Perturbation Workflow](docs/TASK1.md)
- [Task 2 - ALS Gene Perturbations](docs/TASK2.md)
- [Task 3 - Embedding Analysis and Biological Interpretation](docs/TASK3.md)

---

## Getting Started

### Prerequisites

- Python 3.10 (3.9 works if Scanpy is installed from conda-forge).
- GPU with CUDA 11.8 or newer is recommended for GeneFormer embeddings.
- Python packages: PyTorch (CUDA build if available), transformers, anndata, scanpy, numpy, pandas, seaborn, matplotlib, scikit-learn, tqdm.
- Access to GeneFormer weights via the `helical` package or a local checkpoint.

### Environment Setup

```bash
git clone https://github.com/stef1949/helical-coding-challenge.git
cd helical-coding-challenge
python -m venv .venv
# On Windows use: .venv\Scripts\activate
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

Optional extras:

- Install `helical` when GeneFormer weights are available locally: `pip install helical`.
- Adjust the `torch` wheel URL in `requirements.txt` to match your CUDA toolkit.

### Dataset Access

The ALS dataset (GSE174332, BA4 region) is provided as a pre-filtered `.h5ad` file via the Helical S3 bucket referenced inside the notebooks. Place the file in `localtools/` (default expected path in `main.py`) or update the `INPUT_H5AD` constant to point to your copy.

---

## Running the Workflows

### Notebook-driven analysis

1. `third attempt/task1_perturbation_workflow.ipynb`  
   Builds and validates the perturbation utilities, illustrating per-gene knock-up and knock-down edits on AnnData objects.

2. `third attempt/task2_als_perturbations.ipynb`  
   Applies the workflow to canonical ALS genes (SOD1, TARDBP, FUS, C9orf72), embeds the perturbations with GeneFormer_V2, and records latent representations to disk.

3. `third attempt/task3_embedding_analysis.ipynb`  
   Quantifies latent-space consequences with dimensionality reduction, neighbourhood composition, silhouette scores, and disease-to-healthy rescue metrics.

Each notebook contains configuration cells at the top for setting data paths, sampling options, and output locations.

### Command-line batch generation

Use `main.py` when you want to generate perturbations without launching notebooks:

```bash
python main.py \
  --input localtools/counts_combined_filtered_BA4_sALS_PN.h5ad \
  --targets SOD1 TARDBP FUS C9orf72 \
  --subset 'Condition == "ALS" and CellClass == "Neuron" and Region == "BA4"' \
  --output als_simulated.h5ad
```

The script resolves gene symbols against `var_names`, creates knock-up and knock-down variants, concatenates them with the baseline AnnData object, and writes the combined data. Disable plotting by setting `MAKE_PLOTS = False` near the top of the file when running on headless infrastructure.

### Reusable utilities

- `GenePerturbation` (`perturbation.py`) batch generates perturbations, assembles combined AnnData matrices, and validates per-gene edits.
- `EmbeddingAnalyzer` (`analysis.py`) computes rescue scores, neighbour composition, clustering metrics, and provides helper plots for embedding interpretation.
- `data_utils.py` hosts quality-of-life helpers for ALS gene lists, data filtering, and GeneFormer preprocessing (including preservation of raw counts in `.layers['counts']`).

---

## Outputs and Reporting

Generated artefacts are written to `third attempt/outputs/`:

- `images/` folder with UMAPs, clustering comparisons, volcano plots, and scenario dashboards.
- `embedding_metrics.csv` containing clustering quality and rescue statistics.
- `ALS_embeddings.npy` and `metadata.csv` for downstream reuse in external tooling.

Key visualisations referenced in presentations:

- `task3_optimal_clusters.png` summarises the k-means versus Leiden comparison.
- `task3_als_control_distances.png` highlights cell types with the largest disease-control separation.
- Scenario comparison dashboards illustrate disease modelling, rescue, and control perturbations.

---

## Reproducibility Notes

- Random seeds are fixed inside notebooks and helper modules for deterministic sampling.
- GeneFormer runs require consistent model checkpoints; pin `geneformer` or `helical` versions when sharing results.
- Ensure the CUDA toolkit matches the installed PyTorch wheel when using GPU acceleration.
- Perturbation functions copy the input AnnData to avoid modifying raw datasets in place.

---

## References

- GeneFormer documentation: https://helical.readthedocs.io/en/latest/models/geneformer/
- ALS dataset (GSE174332): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332
- Scanpy: https://scanpy.readthedocs.io/

---

## Citation

If you reuse or extend this work, please cite:

```
Ritchie, S. (2025). In-silico Gene Perturbation and Embedding Analysis for ALS. Helical Computational Biology Challenge Submission.
```
