# Task 1 - In-Silico Perturbation Workflow

Documentation for `third attempt/task1_perturbation_workflow.ipynb`.

---

## Objectives

- Construct a reusable pipeline for gene knock-up and knock-down simulations on single-cell RNA-seq data.
- Support single and multi-gene perturbations, cell-type specific targeting, and magnitude control.
- Validate the perturbation outputs before passing datasets to downstream embedding models.

---

## Inputs and Dependencies

- Base dataset: `localtools/counts_combined_filtered_BA4_sALS_PN.h5ad`.
- Helper modules: `perturbation.py` (class `GenePerturbation`), `data_utils.py` for filtering, and Scanpy utilities.
- Optional configuration toggles inside the notebook for sampling, noise injection, and output paths.

---

## Workflow Outline

1. **Load and explore the dataset**  
   - Inspect ALS versus control composition and confirm availability of target genes.

2. **Design the perturbation API**  
   - Implement helper functions that work with dense and sparse `AnnData` matrices.
   - Ensure metadata columns (`perturbed_gene`, `perturbation_type`, `fold_change`) describe each simulated edit.

3. **Single-gene demonstrations**  
   - Apply knock-down and knock-up to a representative gene and visualise before/after distributions.

4. **Multi-gene and cell-type specific perturbations**  
   - Showcase batched edits across ALS-associated genes.
   - Restrict perturbations to selected cell classes to mimic targeted interventions.

5. **Validation and QC**  
   - Compare expression histograms pre/post perturbation.
   - Calculate summary statistics (mean shifts, fraction of altered cells) for each scenario.

6. **Persist outputs**  
   - Write example results (e.g. `perturbed_expression.h5ad`) for integration with Task 2.

---

## Key Outputs

- Annotated `AnnData` objects containing original and perturbed expression matrices.
- QC figures saved under `third attempt/outputs/images/` (cell-type distribution, before/after comparisons).
- Utility functions ready for reuse in scripts such as `main.py`.

---

## Validation Checklist

- Confirm all requested genes are present before perturbation (`filter_genes_present` helper).
- Ensure `AnnData.obs` gains the expected perturbation metadata columns.
- Verify sparsity handling: conversions between CSR/LIL occur only when needed.
- Review QC plots to confirm biological realism (no negative counts, expected shift magnitude).

---

## Suggested Enhancements

- Parameterise noise models (e.g. log-normal scaling) for more realistic knock-up simulations.
- Add unit tests covering sparse and dense matrices to guard against regression.
- Extend the notebook to export perturbation recipes (JSON) for batch automation.
