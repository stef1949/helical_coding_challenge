
<div align="center">
  <p>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/logo_and_text_v2_white.png">
    <source media="(prefers-color-scheme: light)" srcset="docs/assets/logo_and_text_v2.png">
    <img alt="Helical Logo" src="docs/assets/logo_and_text_v2_white.png" width="200">
  </picture>
  </p>
</div>

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

The following table summarises **representative ALS fold changes** across multi-omic layers (gene, protein, and metabolite level) from verified human studies. Values are approximate ranges derived from cohort-level analyses and serve as biological reference points for perturbation benchmarking.

| **Biomarker / Gene / Protein** | **Sample Type** | **Reported Change (ALS vs Control)** | **Approx. Fold Change** | **Direction** | **Relevance** | **Primary Source / DOI** |
|--------------------------------|-----------------|--------------------------------------|--------------------------|----------------|----------------|------------------|
| **Neurofilament light (NfL)** | CSF | 4637.6 pg/mL (ALS) vs 610.4 pg/mL (controls) | ↑ ~7.6× | Upregulated | Neuronal damage; diagnostic marker | [Benatar et al., *Brain* 2023](https://academic.oup.com/brain/article/146/7/2711/6780887) |
| **Chitotriosidase (CHIT1)** | CSF | Elevated ~8–12× | ↑ ~10× | Upregulated | Microglial activation biomarker | [Costa et al., *Diagnostics* 2021](https://www.mdpi.com/2075-4418/11/7/1210) ; [Varghese et al., *J. Neuroinflammation* 2020](https://jneuroinflammation.biomedcentral.com/articles/10.1186/s12974-020-01909-y) |
| **GFAP** | CSF / Plasma | Elevated 1.3–1.6× | ↑ ~1.5× | Upregulated | Astrocytosis; glial activation | [Irwin et al., *Transl. Neurodegener.* 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10809579/) |
| **GAD2 (mRNA)** | Spinal cord | log₂FC −1.78 | ↓ ~0.29× | Downregulated | GABAergic neuron loss | [Patel et al., *Genes* 2020](https://www.mdpi.com/2073-4425/11/4/448) |
| **CALB1 (mRNA)** | Spinal cord | log₂FC −1.96 | ↓ ~0.25× | Downregulated | Motor neuron marker; selective vulnerability | [Patel et al., *Genes* 2020](https://www.mdpi.com/2073-4425/11/4/448) |
| **GABRE (mRNA)** | Spinal cord | log₂FC −1.39 | ↓ ~0.38× | Downregulated | Inhibitory signalling | [Patel et al., *Genes* 2020](https://www.mdpi.com/2073-4425/11/4/448) |
| **miR-206** | Serum / Muscle | 2–3× increase | ↑ ~2–3× | Upregulated | Muscle regeneration; denervation response | [Waller et al., *Neurology* 2017](https://doi.org/10.1016/j.neurobiolaging.2017.03.027) |
| **miR-199a-5p** | Serum / Plasma | Decreased in early ALS | ↓ ~0.5× | Downregulated | Mitochondrial regulation | [Freischmidt et al., *Cell Death Discov.* 2020](https://www.nature.com/articles/s41420-020-00397-6) |
| **Lactate** | CSF / Plasma | Elevated ~1.8–2.3× | ↑ ~2× | Upregulated | Mitochondrial dysfunction; hypometabolism | [Blasco et al., *Neurobiol. Aging* 2018](https://doi.org/10.1016/j.neurobiolaging.2018.08.003) |

**Summary:**
- **Transcriptomic data** show moderate (0.25–0.4×) downregulation of neuronal genes and ~2× upregulation of glial transcripts.
- **Protein biomarkers** (NfL, CHIT1) are robust, with 5–10× elevations across multiple independent cohorts.
- **Metabolites and miRNAs** show intermediate (1.5–3×) alterations, supporting mitochondrial and neuromuscular dysregulation.

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
