# Task 2 - ALS Gene Perturbations

Documentation for `third attempt/task2_als_perturbations.ipynb`.

---

## Objectives

- Apply the Task 1 perturbation workflow to ALS-associated genes across disease and control populations.
- Simulate therapeutic scenarios (disease rescue, disease modelling, control perturbations).
- Generate GeneFormer embeddings that capture transcriptional consequences of each scenario.

---

## Inputs and Resources

- Perturbation helpers and QC utilities authored in Task 1 (`perturbation.py`, `geneformer_helper.py`).
- ALS dataset: `localtools/counts_combined_filtered_BA4_sALS_PN.h5ad`.
- Literature-derived ALS gene panel (e.g. `SOD1`, `TARDBP`, `FUS`, `C9orf72`, `OPTN`, `TBK1`, `SQSTM1`, `VCP`, `ANG`, `VAPB`).
- GeneFormer model weights (`gf-12L-95M-i4096`) available via the `helical` SDK or Hugging Face transformers.

---

## Workflow Outline

1. **Load dataset and perturbation helpers**  
   - Reuse the validated Task 1 functions; ensure baseline AnnData objects are loaded into memory.

2. **Identify ALS gene targets**  
   - Cross-reference literature panel with dataset `var_names` (report missing genes).

3. **Subset cells by condition**  
   - Prepare ALS and control cohorts; keep track of tissue region and cell class metadata for later stratification.

4. **Define perturbation scenarios**  
   - Disease-to-healthy rescue (knock down overexpressed genes in ALS neurons).  
   - Disease modelling (knock up pathogenic genes in control cells).  
   - Control scenario (mild perturbations to evaluate baseline robustness).

5. **Apply perturbations**  
   - Use batch utilities to create per-gene and multi-gene simulations.  
   - Attach descriptive metadata (scenario labels, fold changes, targeted cell types).

6. **Generate GeneFormer embeddings**  
   - Initialise the model via `geneformer_helper.py`, falling back to PCA if weights are unavailable.  
   - Store embeddings in `.obsm['X_geneformer']` and optionally in `.X` for compatibility.

7. **Persist results**  
   - Save scenario `.h5ad` files:  
     - `third attempt/perturbed_scenario1_disease_rescue.h5ad`  
     - `third attempt/perturbed_scenario2_disease_rescue.h5ad`  
     - `third attempt/perturbed_scenario3_disease_rescue.h5ad`  
   - Export NumPy arrays (`ALS_embeddings.npy`) and metadata tables (`metadata.csv`) for rapid reuse.

8. **Quick visual diagnostics**  
   - Produce UMAP/t-SNE previews to verify separation of conditions and perturbation labels.

---

## Key Outputs

- Scenario-specific AnnData files with perturbation metadata and embeddings.
- Embedding arrays (`ALS_embeddings.npy`) aligned with companion metadata (`metadata.csv`).
- Figures saved under `third attempt/outputs/images/` highlighting condition shifts.

---

## Validation Checklist

- Inspect log output from `geneformer_helper.py` to confirm which backend produced embeddings.
- Ensure each scenario file contains `sim_kind`, `perturbation`, and `Condition` columns.
- Confirm embedding dimensionality (256 for GeneFormer; document if PCA fallback is used).
- Verify saved arrays and metadata share the same ordering and counts.

---

## Suggested Enhancements

- Add perturbation magnitude sweeps to evaluate non-linear responses.
- Record random seeds and configuration JSON alongside outputs for reproducibility.
- Integrate automated checks (e.g. `tools/check_geneformer.py`) before long embedding runs.
