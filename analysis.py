"""
Embedding Analysis Utilities
Functions for analyzing perturbation effects in embedding space
"""

import numpy as np
import pandas as pd
import anndata as ad
from typing import List, Dict, Tuple, Optional
import scanpy as sc
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.neighbors import NearestNeighbors


class EmbeddingAnalyzer:
    """
    Analyze perturbation effects in embedding space.
    """
    
    def __init__(self, adata: ad.AnnData, embedding_key: str = 'X_geneformer'):
        """
        Initialize with embedded data.
        
        Parameters:
        -----------
        adata : ad.AnnData
            Data with embeddings in .obsm
        embedding_key : str
            Key in adata.obsm containing embeddings
        """
        self.adata = adata
        self.embedding_key = embedding_key
        
        if embedding_key not in adata.obsm:
            raise ValueError(f"{embedding_key} not found in adata.obsm")
            
        self.embeddings = adata.obsm[embedding_key]
        
    def compute_disease_to_healthy_shift(
        self,
        disease_condition: str,
        healthy_condition: str,
        perturbed_gene: str,
        condition_col: str = 'condition'
    ) -> Dict:
        """
        Compute how perturbation shifts disease state toward healthy state.
        
        Parameters:
        -----------
        disease_condition : str
            Label for disease cells
        healthy_condition : str
            Label for healthy cells
        perturbed_gene : str
            Gene that was perturbed
        condition_col : str
            Column in adata.obs containing condition labels
            
        Returns:
        --------
        Dict
            Statistics about disease->healthy shift
        """
        # Get masks
        disease_mask = self.adata.obs[condition_col] == disease_condition
        healthy_mask = self.adata.obs[condition_col] == healthy_condition
        perturbed_mask = self.adata.obs['perturbed_gene'] == perturbed_gene
        
        # Get embeddings
        disease_emb = self.embeddings[disease_mask]
        healthy_emb = self.embeddings[healthy_mask]
        perturbed_emb = self.embeddings[perturbed_mask]
        
        # Compute centroids
        disease_centroid = np.mean(disease_emb, axis=0)
        healthy_centroid = np.mean(healthy_emb, axis=0)
        perturbed_centroid = np.mean(perturbed_emb, axis=0)
        
        # Compute distances
        disease_to_healthy_dist = np.linalg.norm(healthy_centroid - disease_centroid)
        disease_to_perturbed_dist = np.linalg.norm(perturbed_centroid - disease_centroid)
        perturbed_to_healthy_dist = np.linalg.norm(healthy_centroid - perturbed_centroid)
        
        # Compute rescue score (how much closer to healthy state)
        rescue_score = (disease_to_healthy_dist - perturbed_to_healthy_dist) / disease_to_healthy_dist
        
        # Compute direction alignment
        disease_to_healthy_vec = healthy_centroid - disease_centroid
        disease_to_perturbed_vec = perturbed_centroid - disease_centroid
        
        # Cosine similarity between vectors
        alignment = np.dot(disease_to_healthy_vec, disease_to_perturbed_vec) / (
            np.linalg.norm(disease_to_healthy_vec) * np.linalg.norm(disease_to_perturbed_vec)
        )
        
        return {
            'perturbed_gene': perturbed_gene,
            'disease_to_healthy_distance': disease_to_healthy_dist,
            'disease_to_perturbed_distance': disease_to_perturbed_dist,
            'perturbed_to_healthy_distance': perturbed_to_healthy_dist,
            'rescue_score': rescue_score,
            'direction_alignment': alignment,
            'disease_centroid': disease_centroid,
            'healthy_centroid': healthy_centroid,
            'perturbed_centroid': perturbed_centroid
        }
    
    def neighborhood_analysis(
        self,
        cell_indices: np.ndarray,
        reference_condition: str,
        condition_col: str = 'condition',
        n_neighbors: int = 50
    ) -> Dict:
        """
        Analyze neighborhood composition around perturbed cells.
        
        Parameters:
        -----------
        cell_indices : np.ndarray
            Indices of cells to analyze
        reference_condition : str
            Condition to compare against (e.g., 'healthy')
        condition_col : str
            Column containing condition labels
        n_neighbors : int
            Number of neighbors to consider
            
        Returns:
        --------
        Dict
            Neighborhood statistics
        """
        # Fit nearest neighbors
        nn = NearestNeighbors(n_neighbors=n_neighbors + 1, metric='euclidean')
        nn.fit(self.embeddings)
        
        # Find neighbors for each cell
        distances, indices = nn.kneighbors(self.embeddings[cell_indices])
        
        # Analyze neighborhood composition
        ref_mask = self.adata.obs[condition_col] == reference_condition
        neighbor_is_reference = np.array(ref_mask.values[indices[:, 1:]])  # Exclude self
        
        reference_fraction = np.mean(neighbor_is_reference, axis=1)
        
        return {
            'mean_reference_fraction': np.mean(reference_fraction),
            'std_reference_fraction': np.std(reference_fraction),
            'median_reference_fraction': np.median(reference_fraction),
            'reference_fractions': reference_fraction,
            'neighbor_distances': distances[:, 1:]
        }
    
    def compute_silhouette_scores(
        self,
        labels: np.ndarray,
        sample_size: Optional[int] = None
    ) -> Dict:
        """
        Compute silhouette scores for clustering quality.
        """
        from sklearn.metrics import silhouette_score, silhouette_samples
        
        # Ensure labels are a numpy array for reliable boolean indexing and typing
        labels = np.asarray(labels)
        
        if sample_size and len(labels) > sample_size:
            indices = np.random.choice(len(labels), sample_size, replace=False)
            embeddings_sample = self.embeddings[indices]
            labels_sample = labels[indices]
        else:
            embeddings_sample = self.embeddings
            labels_sample = labels
            
        overall_score = silhouette_score(embeddings_sample, labels_sample)
        sample_scores = np.asarray(silhouette_samples(embeddings_sample, labels_sample))
        
        # Compute mean score per label using explicit boolean masks and handle empty groups
        mean_score_per_label = {}
        for label in np.unique(labels_sample):
            mask = labels_sample == label
            if np.any(mask):
                mean_score_per_label[label] = float(np.mean(sample_scores[mask]))
            else:
                mean_score_per_label[label] = float('nan')
        
        return {
            'overall_score': overall_score,
            'sample_scores': sample_scores,
            'mean_score_per_label': mean_score_per_label
        }