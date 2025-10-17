# 🧬 Helical Coding Challenge — Computational Biology

### Author: **Steph Ritchie (MSc Genetic Manipulation & Molecular Biosciences)**

**Project:** *In-Silico Gene Perturbation and Embedding Analysis for ALS*
**Date:** October 2025

---

## 📖 Overview

This repository contains the implementation of the **Helical Computational Biology Coding Challenge**, focused on **simulation-based gene perturbation** and **latent-space embedding analysis** using **GeneFormer_V2 (gf-12L-95M-i4096)**.
The goal is to simulate *in-silico* gene up-/down-regulation (“knock-up” / “knock-down”) events, apply them to disease-specific targets in **amyotrophic lateral sclerosis (ALS)**, and interpret their downstream effects in an embedding space derived from foundation models.

---

## 🧩 Repository Structure

| File                                                 | Description                                                                                                                                                  |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `task1_perturbation_workflow.ipynb`                  | Defines a scalable *in-silico perturbation pipeline* for gene knock-up and knock-down simulation.                                                            |
| `task2_als_perturbations.ipynb`                      | Applies the workflow to **ALS-related genes**, generates perturbations on control vs. diseased expression profiles, and embeds them using **GeneFormer_V2**. |
| `task3_embedding_analysis.ipynb`                     | Performs dimensionality reduction, clustering, and neighborhood analysis to interpret latent-space shifts induced by perturbations.                          |
| `Helical Coding Challenge Computational Biology.pdf` | Official challenge brief (Helical AI).                                                                                                                       |
| `slides/` *(optional)*                               | One-slide summaries per task for final submission.                                                                                                           |

---

## ⚙️ Environment Setup

### Requirements

* Python ≥ 3.10
* `torch --index-url https://download.pytorch.org/whl/cu126`, `transformers`, `anndata`, `scanpy`, `numpy`, `pandas`, `seaborn`, `matplotlib`, `scikit-learn`
* `helical` or `geneformer` Python package (local install for foundation model access)
* GPU recommended (CUDA ≥ 11.8)

### Installation

```bash
git clone https://github.com/stef1949/helical-coding-challenge.git
cd helical-coding-challenge
pip install -r requirements.txt
```

If using the **Helical local SDK**:

```bash
pip install helical
```

---

## 🧠 Tasks

### **Task 1 — In-Silico Perturbation Workflow**

Develops a flexible simulation framework to mimic *knock-up* and *knock-down* effects:

* Parameterized for multiple gene targets
* Includes scaling options for magnitude of perturbation
* Outputs standardized expression matrices ready for model embedding

**Notebook:** `task1_perturbation_workflow.ipynb`
**Outputs:** `perturbed_expression.h5ad`

---

### **Task 2 — ALS Gene Perturbations**

Applies the perturbation workflow to **ALS-specific genes** (e.g. *SOD1, TARDBP, FUS, C9orf72*), integrating both **control** and **diseased** cell populations.
Embeds perturbation effects using **GeneFormer_V2 (gf-12L-95M-i4096)** into a biologically informed latent space.

**Notebook:** `task2_als_perturbations.ipynb`
**Dataset:**

* `GSE174332` (Motor cortex, BA4)
* Provided via:
  [Helical Candidate Dataset (S3)](https://s3.eu-west-2.amazonaws.com/helical-candidate-datasets/counts_combined_filtered_BA4...)
  **Outputs:**
* `ALS_embeddings.npy`
* `metadata.csv`

---

### **Task 3 — Embedding Interpretation**

Analyzes the embedding space for biological insights:

* **Dimensionality reduction:** UMAP / PCA
* **Cluster metrics:** Silhouette score, Davies–Bouldin index
* **Neighborhood shifts:** Quantifies trajectory movement between perturbed and baseline embeddings
* **Biological interpretation:** Highlights potential compensatory or dysregulated gene networks in ALS

**Notebook:** `task3_embedding_analysis.ipynb`
**Outputs:**

* UMAP plots
* Perturbation shift heatmaps
* Summary metrics table (`embedding_metrics.csv`)

---

## 📊 Evaluation Summary

| Criterion                         | Implementation Focus                                                                      |
| --------------------------------- | ----------------------------------------------------------------------------------------- |
| **Scientific rigor & creativity** | Designed modular, biologically grounded perturbation pipeline with parameterized control. |
| **Proper dataset use**            | Integrated ALS dataset (GSE174332) via preprocessed count matrices.                       |
| **Interpretation depth**          | Quantified embedding displacement and gene-level contributions.                           |
| **Reproducibility**               | Fully documented Jupyter notebooks with fixed random seeds and version metadata.          |

---

## 🧪 Example Output (Task 3)

![UMAP of ALS Perturbations](third attempt\outputs\task3_scenario1_disease_rescue_disease_rescue_embeddings.png)
**Figure:** Visualization of *in-silico* knock-up/down embeddings showing distinct clustering and trajectory shifts between ALS-perturbed and healthy control states.

---

## 🧭 References

* GeneFormer: [https://helical.readthedocs.io/en/latest/models/geneformer/](https://helical.readthedocs.io/en/latest/models/geneformer/)
* Dataset: [GSE174332 — ALS motor cortex single-cell RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174332)

---

## 📦 Citation

If reproducing or extending this work:

```
Ritchie, S. (2025). In-Silico Gene Perturbation and Embedding Analysis for ALS. Helical Computational Biology Challenge Submission.
```