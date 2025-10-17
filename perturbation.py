"""
Gene Perturbation Utilities
Simulates in-silico knock-up and knock-down experiments on single-cell data
"""

import numpy as np
import pandas as pd
import anndata as ad
from typing import List, Dict, Union, Tuple
import scanpy as sc
from copy import deepcopy


class GenePerturbation:
    """
    A class to simulate gene perturbations (knock-up/knock-down) in single-cell data.
    """
    
    def __init__(self, adata: ad.AnnData):
        """
        Initialize with an AnnData object.
        
        Parameters:
        -----------
        adata : ad.AnnData
            Annotated data matrix containing gene expression counts
        """
        self.adata = adata.copy()
        self.original_adata = adata.copy()
        
    def perturb_gene(
     self, 
     gene_name: str, 
     perturbation_type: str = 'knockdown',
     fold_change: float = None,
     cell_indices: np.ndarray = None
 ) -> ad.AnnData:
     """
     Simulate perturbation of a single gene.
     """
     adata_perturbed = self.adata.copy()
     
     # Set default fold changes
     if fold_change is None:
         fold_change = 0.1 if perturbation_type == 'knockdown' else 5.0
         
     # Validate
     if perturbation_type not in ['knockdown', 'knockup']:
         raise ValueError("perturbation_type must be 'knockdown' or 'knockup'")
         
     if gene_name not in adata_perturbed.var_names:
         raise ValueError(f"Gene {gene_name} not found in dataset")
         
     # Get gene index
     gene_idx = np.where(adata_perturbed.var_names == gene_name)[0][0]
     
     # Select cells to perturb
     if cell_indices is None:
         cell_indices = np.arange(adata_perturbed.n_obs)
     else:
         cell_indices = np.asarray(cell_indices)
         if cell_indices.dtype == bool:
             cell_indices = np.where(cell_indices)[0]
         else:
             cell_indices = cell_indices.astype(int)
     
     # FIXED: Properly handle both sparse and dense matrices
     import scipy.sparse as sp
     
     if sp.issparse(adata_perturbed.X):
         # Convert to LIL format for efficient modification
         X_mod = adata_perturbed.X.tolil(copy=True)
         for cell_idx in np.atleast_1d(cell_indices):
             X_mod[cell_idx, gene_idx] = X_mod[cell_idx, gene_idx] * fold_change
         adata_perturbed.X = X_mod.tocsr()  # Convert back to CSR
     else:
         # Dense array
         if adata_perturbed.X is not None:
             adata_perturbed.X[cell_indices, gene_idx] = \
                 adata_perturbed.X[cell_indices, gene_idx] * fold_change
     
     # Add perturbation metadata
     adata_perturbed.obs['perturbed_gene'] = 'none'
     adata_perturbed.obs['perturbation_type'] = 'none'
     adata_perturbed.obs['fold_change'] = 1.0
     
     # Only set for perturbed cells
     cell_names = adata_perturbed.obs.index[cell_indices]
     adata_perturbed.obs.loc[cell_names, 'perturbed_gene'] = gene_name
     adata_perturbed.obs.loc[cell_names, 'perturbation_type'] = perturbation_type
     adata_perturbed.obs.loc[cell_names, 'fold_change'] = fold_change
     
     # Store perturbation info
     adata_perturbed.uns['perturbation'] = {
         'gene': gene_name,
         'type': perturbation_type,
         'fold_change': fold_change,
         'n_cells_perturbed': len(cell_indices)
     }
     
     return adata_perturbed
    
    def batch_perturbation(
        self,
        gene_list: List[str],
        perturbation_type: str = 'knockdown',
        fold_change: float = None,
        cells_per_gene: int = None
    ) -> Dict[str, ad.AnnData]:
        """
        Apply perturbations to multiple genes.
        
        Parameters:
        -----------
        gene_list : List[str]
            List of genes to perturb
        perturbation_type : str
            'knockdown' or 'knockup'
        fold_change : float
            Multiplicative factor for expression change
        cells_per_gene : int
            Number of cells to perturb per gene. If None, use all cells
            
        Returns:
        --------
        Dict[str, ad.AnnData]
            Dictionary mapping gene names to perturbed AnnData objects
        """
        perturbed_data = {}
        
        for gene in gene_list:
            if gene not in self.adata.var_names:
                print(f"Warning: Gene {gene} not found in dataset, skipping...")
                continue
                
            # Select random subset of cells if specified
            if cells_per_gene is not None:
                n_cells = min(cells_per_gene, self.adata.n_obs)
                cell_indices = np.random.choice(
                    self.adata.n_obs, 
                    size=n_cells, 
                    replace=False
                )
            else:
                cell_indices = None
                
            # Perturb gene
            adata_perturbed = self.perturb_gene(
                gene_name=gene,
                perturbation_type=perturbation_type,
                fold_change=fold_change,
                cell_indices=cell_indices
            )
            
            perturbed_data[gene] = adata_perturbed
            
        return perturbed_data
    
    def create_perturbation_matrix(
        self,
        gene_list: List[str],
        perturbation_types: List[str] = ['knockdown', 'knockup'],
        fold_changes: Dict[str, float] = None
    ) -> ad.AnnData:
        """
        Create a combined dataset with multiple perturbations for comparison.
        
        Parameters:
        -----------
        gene_list : List[str]
            Genes to perturb
        perturbation_types : List[str]
            Types of perturbations to apply
        fold_changes : Dict[str, float]
            Custom fold changes for each perturbation type
            
        Returns:
        --------
        ad.AnnData
            Combined dataset with all perturbations
        """
        if fold_changes is None:
            fold_changes = {'knockdown': 0.1, 'knockup': 5.0}
            
        # Start with original data
        combined_data = [self.original_adata.copy()]
        combined_data[0].obs['condition'] = 'control'
        combined_data[0].obs['perturbed_gene'] = 'none'
        combined_data[0].obs['perturbation_type'] = 'none'
        
        # Add perturbations
        for gene in gene_list:
            if gene not in self.adata.var_names:
                continue
                
            for pert_type in perturbation_types:
                adata_pert = self.perturb_gene(
                    gene_name=gene,
                    perturbation_type=pert_type,
                    fold_change=fold_changes.get(pert_type)
                )
                adata_pert.obs['condition'] = f'{gene}_{pert_type}'
                combined_data.append(adata_pert)
        
        # Concatenate all datasets
        adata_combined = ad.concat(
            combined_data,
            join='outer',
            label='perturbation_id',
            keys=[f'pert_{i}' for i in range(len(combined_data))]
        )
        
        return adata_combined


def validate_perturbation(
    adata_original: ad.AnnData,
    adata_perturbed: ad.AnnData,
    gene_name: str
) -> Dict:
    """
    Validate that perturbation was applied correctly.
    
    Parameters:
    -----------
    adata_original : ad.AnnData
        Original data
    adata_perturbed : ad.AnnData
        Perturbed data
    gene_name : str
        Name of perturbed gene
        
    Returns:
    --------
    Dict
        Validation statistics
    """
    gene_idx = np.where(adata_original.var_names == gene_name)[0][0]
    
    # Get expression values
    if hasattr(adata_original.X, 'toarray'):
        expr_orig = adata_original.X.toarray()[:, gene_idx]
        expr_pert = adata_perturbed.X.toarray()[:, gene_idx]
    else:
        expr_orig = adata_original.X[:, gene_idx]
        expr_pert = adata_perturbed.X[:, gene_idx]
    
    # Calculate statistics
    stats = {
        'gene': gene_name,
        'mean_expression_original': np.mean(expr_orig),
        'mean_expression_perturbed': np.mean(expr_pert),
        'fold_change_achieved': np.mean(expr_pert) / (np.mean(expr_orig) + 1e-10),
        'cells_affected': np.sum(expr_orig != expr_pert),
        'zero_expression_cells_original': np.sum(expr_orig == 0),
        'zero_expression_cells_perturbed': np.sum(expr_pert == 0)
    }
    
    return stats
