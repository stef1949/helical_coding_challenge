"""
Data Loading and Preprocessing Utilities
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from typing import List, Optional


def identify_als_genes() -> List[str]:
    """
    Return list of known ALS-associated genes.
    """
    als_genes = [
        'SOD1',      # Superoxide dismutase 1
        'C9orf72',   # Most common genetic cause
        'TARDBP',    # TDP-43
        'FUS',       # Fused in sarcoma
        'OPTN',      # Optineurin
        'TBK1',      # TANK binding kinase 1
        'SQSTM1',    # Sequestosome 1
        'VCP',       # Valosin containing protein
        'ANG',       # Angiogenin
        'VAPB',      # VAMP associated protein B/C
        'CHCHD10',   # Coiled-coil-helix domain 10
        'NEK1',      # NIMA related kinase 1
    ]
    return als_genes


def filter_genes_present(adata: ad.AnnData, gene_list: List[str]) -> List[str]:
    """
    Filter gene list to only those present in dataset.
    """
    present_genes = [g for g in gene_list if g in adata.var_names]
    missing_genes = [g for g in gene_list if g not in adata.var_names]
    
    if missing_genes:
        print(f"Warning: {len(missing_genes)} genes not found")
        print(f"  {', '.join(missing_genes)}")
    
    print(f"Found {len(present_genes)} / {len(gene_list)} genes")
    return present_genes


def preprocess_for_geneformer(
    adata: ad.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_genes: Optional[int] = None
) -> ad.AnnData:
    """
    Preprocess single-cell data for GeneFormer.
    """
    adata = adata.copy()
    
    print(f"Starting: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Filter
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes:
        sc.pp.filter_cells(adata, max_genes=max_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"After filtering: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Store raw counts (materialize to in-memory arrays/sparse matrices; handle h5/backed datasets)
    x = adata.X
    if x is None:
        adata.layers['counts'] = None
    else:
        try:
            if sp.issparse(x):
                # materialize as an in-memory CSR sparse matrix
                try:
                    adata.layers['counts'] = x.tocsr()
                except Exception:
                    adata.layers['counts'] = sp.csr_matrix(x)
            else:
                # materialize array-like / h5-backed datasets to a NumPy array
                adata.layers['counts'] = np.asarray(x)
        except Exception:
            # fallback: store the original object to avoid raising on unsupported types
            adata.layers['counts'] = x

    return adata