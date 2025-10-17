# Task 3 - Embedding Analysis and Biological Interpretation

Documentation for `third attempt/task3_embedding_analysis.ipynb`.

---

## Objectives

- Analyse GeneFormer (or PCA fallback) embeddings produced in Task 2.
- Characterise clustering structure, disease versus control separation, and perturbation impact.
- Derive biological insight through differential expression and pathway enrichment hooks.

---

## Inputs and Dependencies

- Scenario `.h5ad` files generated in Task 2 (expects embeddings in `.obsm['X_geneformer']` or `.X`).
- Metrics and plotting libraries: Scanpy, scikit-learn, seaborn, matplotlib, statsmodels.
- Optional pathway resources for enrichment analysis (external services such as g:Profiler or Enrichr).

---

## Workflow Outline

1. **Data loading and sanity checks**  
   - Detect embedding location (auto-switch between `.obsm` and `.X`).  
   - Confirm required metadata columns (`Condition`, `CellType`, `Scenario`, etc.).

2. **Exploratory analysis**  
   - Summaries of cell counts per condition/cell type.  
   - Diagnostic plots for magnitude and variance of embeddings.

3. **Dimensionality reduction**  
   - Compute PCA (if not already present) and UMAP for 2D visualisation.  
   - Save annotated plots for each scenario and colouring scheme.

4. **Clustering and validation**  
   - Run Leiden and k-means clustering over embeddings.  
   - Evaluate cluster quality via silhouette scores, Davies-Bouldin index, and within-cluster dispersion.

5. **Cluster characterisation**  
   - Map clusters to known cell types and compute enrichment tables.  
   - Highlight shifts in cluster composition between ALS and control samples.

6. **Disease versus control comparisons**  
   - Quantify centroid distances and nearest-neighbour mixing between phenotypes.  
   - Measure rescue scores for perturbation scenarios.

7. **Differential expression and pathways**  
   - Perform per-cell-type differential expression on perturbed versus baseline groups.  
   - Prepare ranked gene lists for downstream pathway tools (enrichment performed externally).

8. **Scenario comparison dashboard**  
   - Aggregate metrics and plots to compare disease modelling, rescue, and control interventions.

9. **Summary and recommendations**  
   - Capture actionable insights and possible follow-up analyses.

---

## Key Outputs

- Figures in `third attempt/outputs/images/` (UMAPs, clustering comparison, volcano plots, ALS-control distance charts).
- Aggregated metrics table: `third attempt/outputs/embedding_metrics.csv`.
- Notebook narrative capturing interpretation highlights and suggested next steps.

---

## Validation Checklist

- Confirm embeddings are normalised or scaled consistently before clustering.  
- Review silhouette and Davies-Bouldin scores for unexpected degradations.  
- Ensure differential expression uses appropriate covariates or filters to avoid artefacts.  
- Cross-check saved figures against notebook widgets to verify correct legend labels.

---

## Suggested Enhancements

- Automate pathway enrichment using APIs (if network access is permitted) and embed results into the notebook.  
- Incorporate permutation testing for rescue scores to quantify statistical significance.  
- Export summary JSON files with metric values for integration into reporting dashboards.
