# simulate_perturbations.py
# Works with dense or sparse AnnData. No np.clip on sparse.
import os
from typing import Iterable, List, Optional

import numpy as np
from scipy.sparse.csgraph import connected_components
import anndata as ad
from scipy import sparse as sp

# Optional plotting (set MAKE_PLOTS=False if you don't want figures)
MAKE_PLOTS = True
try:
    import scanpy as sc
except Exception:
    MAKE_PLOTS = False


# ----------------------------
# Config: edit these as needed
# ----------------------------
INPUT_H5AD = "localtools\counts_combined_filtered_BA4_sALS_PN.h5ad"   # path to your dataset
OUTPUT_BASENAME = "als_simulated"                  # base for outputs
# The biological subset (pandas query string on .obs)
SUBSET_QUERY = 'Condition == "ALS" and CellClass == "Neuron" and Region == "BA4"'

# Targets to perturb. We'll resolve case-insensitively against var_names.
TARGET_GENES_SYMBOLS = ["SOD1", "TARDBP", "FUS", "C9orf72"]  # will be filtered to those present

# If your raw counts live in a layer, put its name here (e.g., "counts"); else None uses .X
OPERATE_ON_LAYER: Optional[str] = None   # e.g., "counts" if present


# ------------------------------------------------
# Utilities
# ------------------------------------------------
def resolve_genes(
    adata: ad.AnnData,
    requested: Iterable[str],
    prefer_symbols: bool = True,
    ensid_col: str = "ENSID",
) -> List[str]:
    """
    Return the list of gene column names (matching adata.var_names) that correspond to 'requested'.
    - Case-insensitive match on symbols in var_names.
    - If not found and prefer_symbols is False, try matching ENSG IDs via adata.var[ensid_col].
    """
    requested = list(requested)
    # quick exits
    if len(requested) == 0:
        return []

    # Build case-insensitive map for symbols
    var_lower_to_true = {g.lower(): g for g in adata.var_names}
    resolved = []

    # First pass: try symbols (case-insensitive)
    for r in requested:
        key = r.lower()
        if key in var_lower_to_true:
            resolved.append(var_lower_to_true[key])

    # Second pass (optional): try ENSG IDs in var column
    if not prefer_symbols and ensid_col in adata.var.columns:
        ensid_series = adata.var[ensid_col].astype(str).str.upper()
        ensid_to_symbol = {ensg: sym for ensg, sym in zip(ensid_series, adata.var_names)}
        for r in requested:
            r_up = str(r).upper()
            if r_up in ensid_to_symbol and ensid_to_symbol[r_up] not in resolved:
                resolved.append(ensid_to_symbol[r_up])

    return resolved


def get_row_indices_from_query(adata: ad.AnnData, subset_query: Optional[str]) -> np.ndarray:
    """Return integer row indices given a pandas query on .obs (or all rows if None)."""
    if subset_query is None or subset_query.strip() == "":
        return np.arange(adata.n_obs)
    idx = adata.obs.query(subset_query).index
    # convert index labels to integer positions
    return adata.obs_names.get_indexer(idx)


def simulate_perturbations(
    adata: ad.AnnData,
    targets: List[str],
    direction: str = "up",            # "up" or "down"
    fold_change: float = 2.0,         # 2.0 => 2x up; "down" uses reciprocal
    subset_query: Optional[str] = None,
    layer: Optional[str] = None       # e.g., "counts"; if None uses X
) -> ad.AnnData:
    """
    Multiply the expression of `targets` by a fold-change in a selected subset.
    Works for both dense and sparse matrices. Does not clip (not needed for positive FC).
    """
    if not targets:
        raise ValueError("No valid target genes to perturb (after resolving against var_names).")

    fc = float(fold_change) if direction.lower() == "up" else 1.0 / float(fold_change)

    rows = get_row_indices_from_query(adata, subset_query)
    cols = adata.var_names.get_indexer(targets)
    if np.any(cols < 0):
        raise ValueError("Some targets could not be located in var_names after resolution.")

    rows = np.asarray(rows, dtype=np.intp)
    cols = np.asarray(cols, dtype=np.intp)
    # choose source matrix
    Xsrc = adata.layers[layer] if layer else adata.X
    is_sparse = sp.issparse(Xsrc)

    if is_sparse:
        X_mod = Xsrc.tocsr(copy=True)
        tgt_rows = np.unique(rows)
        col_mask = np.zeros(X_mod.shape[1], dtype=bool)
        col_mask[cols] = True
        for i in tgt_rows:
            row_start = X_mod.indptr[i]
            row_end = X_mod.indptr[i + 1]
            if row_start == row_end:
                continue
            row_cols = X_mod.indices[row_start:row_end]
            mask = col_mask[row_cols]
            if mask.any():
                X_mod.data[row_start:row_end][mask] *= fc
    else:
        X_mod = np.array(Xsrc, copy=True)
        X_mod[np.ix_(rows, cols)] *= fc

    # Build new AnnData:
    # - If we operated on X: put X_mod into X
    # - If we operated on a layer: keep X intact, copy layers, and write X_mod to that layer
    out = ad.AnnData(
        X_mod if layer is None else (adata.X.copy() if not sp.issparse(adata.X) else adata.X.copy()),
        obs=adata.obs.copy(),
        var=adata.var.copy(),
        obsm=adata.obsm.copy(),
        varm=adata.varm.copy() if hasattr(adata, "varm") else None,
        obsp=adata.obsp.copy() if hasattr(adata, "obsp") else None,
        layers=None,
    )
    if layer is not None:
        out.layers = dict(adata.layers)
        out.layers[layer] = X_mod

    out.obs["perturbation"] = f"{'+'.join(targets)}_{direction}_{fold_change}x"
    out.obs["perturbed_subset"] = subset_query or "all"
    return out


def main():
    # ----------------------------
    # Load data
    # ----------------------------
    if not os.path.exists(INPUT_H5AD):
        raise FileNotFoundError(f"Cannot find {INPUT_H5AD}")

    print(f"Loading: {INPUT_H5AD}")
    adata = ad.read_h5ad(INPUT_H5AD)
    print("Shape:", adata.shape)
    print("obs columns:", list(adata.obs.columns)[:10], "...")
    print("var columns:", list(adata.var.columns)[:10], "...")

    # ----------------------------
    # Resolve targets present
    # ----------------------------
    targets = resolve_genes(adata, TARGET_GENES_SYMBOLS, prefer_symbols=True, ensid_col="ENSID")
    if not targets:
        raise RuntimeError("None of the target gene symbols were found in var_names.")
    print("Resolved targets:", targets)

    # ----------------------------
    # Simulate up/down
    # ----------------------------
    print("Simulating knock_up...")
    pert_up = simulate_perturbations(
        adata,
        targets,
        direction="up",
        fold_change=2.0,
        subset_query=SUBSET_QUERY,
        layer=OPERATE_ON_LAYER
    )

    print("Simulating knock_down...")
    pert_dn = simulate_perturbations(
        adata,
        targets,
        direction="down",
        fold_change=2.0,
        subset_query=SUBSET_QUERY,
        layer=OPERATE_ON_LAYER
    )

    # ----------------------------
    # Concatenate baseline + perturbed
    # ----------------------------
    print("Concatenating...")
    combined = ad.concat(
        {"baseline": adata, "knock_up": pert_up, "knock_down": pert_dn},
        label="sim_kind",
        axis=0,
        join="outer",
    )
    combined.obs_names_make_unique()

    # ----------------------------
    # Save results
    # ----------------------------
    out_path = f"{OUTPUT_BASENAME}.h5ad"
    print(f"Writing: {out_path}")
    # Use lzf (fast) or gzip (smaller). lzf is fine for iterative work.
    combined.write_h5ad(out_path, compression="lzf")
    print("Done.")

    # ----------------------------
    # Optional: quick UMAP using existing representation if available
    # If you plan to embed with Geneformer later, set use_rep to that obsm key.
    # ----------------------------
    if MAKE_PLOTS:
        print("Computing neighbors/UMAP on PCA (for a quick sanity plot)...")
        # If you already have a preferred layer/embedding, adjust here:
        # e.g., sc.pp.neighbors(combined, use_rep="X_geneformer")
        sc.pp.normalize_total(combined, target_sum=1e4, inplace=True)
        sc.pp.log1p(combined)
        sc.pp.highly_variable_genes(combined, flavor="seurat", n_top_genes=3000, subset=True)
        sc.pp.scale(combined, max_value=10, zero_center=False)
        sc.tl.pca(combined, n_comps=50)
        sc.pp.neighbors(combined, n_neighbors=15, use_rep="X_pca")
        sc.tl.umap(combined, min_dist=0.3)
        # color keys must exist in obs
        color_keys = [k for k in ["Condition", "CellType", "CellClass", "sim_kind", "perturbation"] if k in combined.obs.columns]
        sc.pl.umap(combined, color=color_keys, wspace=0.4, frameon=False, ncols=3, show=False, save=f"_{OUTPUT_BASENAME}.png")
        print("Saved UMAP to `figures/umap_als_simulated.png` (Scanpy default folder).")


if __name__ == "__main__":
    main()
