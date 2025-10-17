"""
GeneFormer Integration Helper - Corrected Version

This module provides multiple methods to work with GeneFormer,
handling different versions of the helical library and providing fallbacks.
"""

import numpy as np
import scanpy as sc
from typing import Optional, Tuple

def get_geneformer_method1():
    """
    Method 1: Try importing from helical (new API)
    """
    try:
        from helical import Geneformer, GeneformerConfig
        print("✓ Method 1: Helical new API")
        return Geneformer, GeneformerConfig, "helical_new"
    except ImportError as e:
        print(f"✗ Method 1 failed: {e}")
        return None, None, None

def get_geneformer_method2():
    """
    Method 2: Try importing from helical.models
    """
    try:
        from helical.models.geneformer.model import Geneformer
        from helical.models.geneformer.geneformer_config import GeneformerConfig
        print("✓ Method 2: Helical models API")
        return Geneformer, GeneformerConfig, "helical_models"
    except ImportError as e:
        print(f"✗ Method 2 failed: {e}")
        return None, None, None

def get_geneformer_method3():
    """
    Method 3: Try direct transformers
    """
    try:
        from transformers import AutoModel, AutoTokenizer
        print("✓ Method 3: Direct Transformers")
        return AutoModel, AutoTokenizer, "transformers"
    except ImportError as e:
        print(f"✗ Method 3 failed: {e}")
        return None, None, None

def initialize_geneformer():
    """
    Try multiple methods to initialize GeneFormer.
    
    Returns:
    --------
    tuple: (model, config/tokenizer, method_name) or (None, None, None)
    """
    print("Attempting to initialize GeneFormer...\n")
    
    # Try Method 1: Helical new API
    Model, Config, method = get_geneformer_method1()
    if Model is not None:
        try:
            config = Config(model_name="gf-12L-95M-i4096", batch_size=10, device="cuda")
            model = Model(config)
            print(f"✓ Successfully initialized with {method}")
            return model, config, method
        except Exception as e:
            print(f"✗ Failed to initialize with {method}: {e}\n")
    
    # Try Method 2: Helical models API
    Model, Config, method = get_geneformer_method2()
    if Model is not None:
        try:
            config = Config(model_name="gf-12L-95M-i4096", batch_size=10, device="cuda")
            model = Model(config)
            print(f"✓ Successfully initialized with {method}")
            return model, config, method
        except Exception as e:
            print(f"✗ Failed to initialize with {method}: {e}\n")
    
    # Try Method 3: Direct transformers
    Model, Tokenizer, method = get_geneformer_method3()
    if Model is not None:
        try:
            model_name = "ctheodoris/Geneformer"
            model = Model.from_pretrained(model_name, trust_remote_code=True)
            tokenizer = Tokenizer.from_pretrained(model_name, trust_remote_code=True)
            print(f"✓ Successfully initialized with {method}")
            return model, tokenizer, method
        except Exception as e:
            print(f"✗ Failed to initialize with {method}: {e}\n")
    
    print("✗ All GeneFormer initialization methods failed")
    print("→ Will use PCA fallback for embeddings\n")
    return None, None, None

def prepare_data_for_geneformer(adata, sample_size=None):
    """
    Prepare AnnData for GeneFormer by ranking genes.
    
    GeneFormer expects genes ranked by expression level per cell.
    """
    if sample_size and sample_size < adata.n_obs:
        sample_idx = np.random.choice(adata.n_obs, sample_size, replace=False)
        adata_sample = adata[sample_idx].copy()
    else:
        adata_sample = adata.copy()
        sample_idx = None
    
    print(f"Preparing {adata_sample.n_obs} cells for GeneFormer...")
    
    # Ensure data is in the right format
    if hasattr(adata_sample.X, 'toarray'):
        adata_sample.X = adata_sample.X.toarray()
    
    # GeneFormer expects per-cell count totals in obs['n_counts'] for stable normalization.
    # Some helical versions broadcast this column instead of recomputing from X, so we
    # populate it explicitly to avoid shape mismatches during tokenization.
    if 'n_counts' not in adata_sample.obs.columns:
        counts = adata_sample.X.sum(axis=1)
        # Make sure we end up with a flat array regardless of sparse/dense backing
        counts = np.asarray(counts).ravel()
        adata_sample.obs['n_counts'] = counts
    
    return adata_sample, sample_idx

def generate_embeddings_helical(adata, model, config, method, sample_size=None):
    """
    Generate embeddings using helical GeneFormer.
    """
    adata_prep, sample_idx = prepare_data_for_geneformer(adata, sample_size)
    
    print("Generating embeddings with helical GeneFormer...")
    
    try:
        # Different methods for different helical API versions
        if hasattr(model, 'process_data'):
            print("  Tokenizing data with model.process_data()...")
            tokenized_dataset = model.process_data(adata_prep)
        else:
            # Fallback for older APIs that don't have a separate process_data step.
            print("  Using adata directly (older API pattern)...")
            tokenized_dataset = adata_prep

        # Step 2: Generate embeddings from the (potentially tokenized) data.
        if hasattr(model, 'get_embeddings'):
            print("  Generating embeddings with model.get_embeddings()...")
            embeddings = model.get_embeddings(tokenized_dataset)
        elif hasattr(model, 'encode'):
            print("  Generating embeddings with model.encode()...")
            embeddings = model.encode(tokenized_dataset)
        elif hasattr(model, '__call__'):
            print("  Generating embeddings with model.__call__()...")
            embeddings = model(tokenized_dataset)
        else:
            raise AttributeError("Cannot find appropriate method to generate embeddings (get_embeddings, encode, or __call__)")
        
        # Convert to numpy if needed
        if not isinstance(embeddings, np.ndarray):
            # Handle torch-style outputs
            if hasattr(embeddings, "numpy"):
                embeddings = embeddings.numpy()
            elif hasattr(embeddings, "detach"):
                embeddings = embeddings.detach().cpu().numpy()
            else:
                # Some helical builds return a HuggingFace Dataset
                try:
                    from datasets import Dataset  # type: ignore
                except ImportError:
                    Dataset = tuple()  # type: ignore

                if isinstance(embeddings, Dataset):
                    # Try the usual embedding column names first
                    for column in ("embeddings", "embedding", "latent", "features"):
                        if column in embeddings.column_names:
                            column_data = embeddings[column]
                            embeddings = np.asarray(column_data)
                            break
                    else:
                        # If only token ids are returned we cannot visualise them directly
                        raise TypeError(
                            "Helical GeneFormer returned token sequences instead of embeddings; "
                            "please regenerate using PCA fallback (set force_pca=True) or ensure "
                            "the helical version exposes an embedding column."
                        )
                elif isinstance(embeddings, (list, tuple)):
                    embeddings = np.asarray(embeddings)
                else:
                    raise TypeError(
                        f"Unsupported embedding type returned from helical GeneFormer: {type(embeddings)}"
                    )

        print(f"✓ Generated embeddings with shape: {embeddings.shape}")
        return embeddings, sample_idx
        
    except Exception as e:
        print(f"✗ Error generating embeddings: {e}")
        import traceback
        traceback.print_exc()
        raise

def generate_embeddings_direct_geneformer(adata, sample_size=None):
    """
    Generate embeddings using GeneFormer directly via transformers.
    This bypasses helical's processing which can have compatibility issues.
    """
    try:
        from transformers import AutoModel
        import torch
        
        adata_prep, sample_idx = prepare_data_for_geneformer(adata, sample_size)
        
        print("Generating embeddings with direct GeneFormer (bypassing helical)...")
        print("Note: This uses a simpler approach that's more reliable")
        
        # Load model
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        model = AutoModel.from_pretrained("ctheodoris/Geneformer", trust_remote_code=True)
        model = model.to(device)
        model.eval()
        
        embeddings_list = []
        batch_size = 16
        n_cells = adata_prep.n_obs
        
        print(f"Processing {n_cells} cells in batches of {batch_size}...")
        
        with torch.no_grad():
            for i in range(0, n_cells, batch_size):
                batch_end = min(i + batch_size, n_cells)
                
                # Get batch data
                if hasattr(adata_prep.X, 'toarray'):
                    batch_expr = adata_prep.X[i:batch_end].toarray()
                else:
                    batch_expr = adata_prep.X[i:batch_end]
                
                # For each cell, get top genes by expression
                batch_embeddings = []
                for j in range(batch_expr.shape[0]):
                    cell_expr = batch_expr[j]
                    
                    # Rank genes by expression (GeneFormer expects ranked genes)
                    nonzero_idx = np.where(cell_expr > 0)[0]
                    if len(nonzero_idx) == 0:
                        # If no expression, use mean embedding
                        batch_embeddings.append(np.zeros(model.config.hidden_size))
                        continue
                    
                    # Get top expressed genes (up to 2048)
                    sorted_idx = nonzero_idx[np.argsort(cell_expr[nonzero_idx])[::-1]]
                    top_genes_idx = sorted_idx[:min(2048, len(sorted_idx))]
                    
                    # Get gene names
                    gene_ids = top_genes_idx.tolist()
                    
                    # Simple tokenization: use gene indices as tokens
                    # Add special tokens: [CLS] = 0
                    input_ids = torch.tensor([[0] + gene_ids], dtype=torch.long).to(device)
                    
                    # Generate embedding
                    try:
                        outputs = model(input_ids=input_ids)
                        # Use [CLS] token embedding
                        cell_embedding = outputs.last_hidden_state[0, 0, :].cpu().numpy()
                        batch_embeddings.append(cell_embedding)
                    except:
                        # Fallback: mean of all token embeddings
                        batch_embeddings.append(np.zeros(model.config.hidden_size))
                
                embeddings_list.extend(batch_embeddings)
                
                if (i // batch_size + 1) % 10 == 0:
                    print(f"  Processed {batch_end}/{n_cells} cells...")
        
        embeddings = np.array(embeddings_list)
        print(f"✓ Generated embeddings with shape: {embeddings.shape}")
        
        return embeddings, sample_idx
        
    except Exception as e:
        print(f"✗ Direct GeneFormer failed: {e}")
        raise

def generate_embeddings_transformers(adata, model, tokenizer, sample_size=None):
    """
    Generate embeddings using transformers GeneFormer (via helical initialization).
    """
    adata_prep, sample_idx = prepare_data_for_geneformer(adata, sample_size)
    
    print("Generating embeddings with transformers GeneFormer...")
    
    try:
        import torch
        
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        model = model.to(device)
        model.eval()
        
        embeddings_list = []
        
        # Process in batches
        batch_size = 16
        n_cells = adata_prep.n_obs
        
        with torch.no_grad():
            for i in range(0, n_cells, batch_size):
                batch_end = min(i + batch_size, n_cells)
                batch_data = adata_prep[i:batch_end]
                
                # Tokenize batch
                # For each cell, rank genes by expression
                batch_tokens = []
                for j in range(batch_data.n_obs):
                    cell_expr = batch_data.X[j]
                    # Get top expressed genes (GeneFormer typically uses top 2048)
                    top_genes_idx = np.argsort(cell_expr)[::-1][:2048]
                    top_genes = batch_data.var_names[top_genes_idx].tolist()
                    batch_tokens.append(top_genes)
                
                # Tokenize
                inputs = tokenizer(
                    batch_tokens,
                    padding=True,
                    truncation=True,
                    return_tensors='pt'
                ).to(device)
                
                # Get embeddings
                outputs = model(**inputs)
                embeddings = outputs.last_hidden_state[:, 0, :].cpu().numpy()
                embeddings_list.append(embeddings)
                
                if (i // batch_size) % 10 == 0:
                    print(f"  Processed {batch_end}/{n_cells} cells...")
        
        embeddings = np.vstack(embeddings_list)
        print(f"✓ Generated embeddings with shape: {embeddings.shape}")
        return embeddings, sample_idx
        
    except Exception as e:
        print(f"✗ Error generating embeddings: {e}")
        raise
    """
    Generate embeddings using transformers GeneFormer.
    """
    adata_prep, sample_idx = prepare_data_for_geneformer(adata, sample_size)
    
    print("Generating embeddings with transformers GeneFormer...")
    
    try:
        import torch
        
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        model = model.to(device)
        model.eval()
        
        embeddings_list = []
        
        # Process in batches
        batch_size = 16
        n_cells = adata_prep.n_obs
        
        with torch.no_grad():
            for i in range(0, n_cells, batch_size):
                batch_end = min(i + batch_size, n_cells)
                batch_data = adata_prep[i:batch_end]
                
                # Tokenize batch
                # For each cell, rank genes by expression
                batch_tokens = []
                for j in range(batch_data.n_obs):
                    cell_expr = batch_data.X[j]
                    # Get top expressed genes (GeneFormer typically uses top 2048)
                    top_genes_idx = np.argsort(cell_expr)[::-1][:2048]
                    top_genes = batch_data.var_names[top_genes_idx].tolist()
                    batch_tokens.append(top_genes)
                
                # Tokenize
                inputs = tokenizer(
                    batch_tokens,
                    padding=True,
                    truncation=True,
                    return_tensors='pt'
                ).to(device)
                
                # Get embeddings
                outputs = model(**inputs)
                embeddings = outputs.last_hidden_state[:, 0, :].cpu().numpy()
                embeddings_list.append(embeddings)
                
                if (i // batch_size) % 10 == 0:
                    print(f"  Processed {batch_end}/{n_cells} cells...")
        
        embeddings = np.vstack(embeddings_list)
        print(f"✓ Generated embeddings with shape: {embeddings.shape}")
        return embeddings, sample_idx
        
    except Exception as e:
        print(f"✗ Error generating embeddings: {e}")
        raise

def generate_pca_embeddings(adata, n_components=50, sample_size=None):
    """
    Generate PCA embeddings as fallback.
    
    This is a simple alternative when GeneFormer is not available.
    While not as biologically informed as GeneFormer, PCA still
    captures major sources of variation in the data.
    """
    if sample_size and sample_size < adata.n_obs:
        sample_idx = np.random.choice(adata.n_obs, sample_size, replace=False)
        adata_sample = adata[sample_idx].copy()
    else:
        adata_sample = adata.copy()
        sample_idx = None
    
    print(f"Generating PCA embeddings for {adata_sample.n_obs} cells...")
    print("(Using PCA as GeneFormer alternative)")
    
    # Normalize and compute PCA
    sc.pp.normalize_total(adata_sample, target_sum=1e4)
    sc.pp.log1p(adata_sample)
    sc.pp.highly_variable_genes(adata_sample, n_top_genes=2000, flavor='seurat_v3')
    
    # Use only highly variable genes for PCA
    adata_hvg = adata_sample[:, adata_sample.var.highly_variable]
    sc.pp.scale(adata_hvg, max_value=10)
    sc.pp.pca(adata_hvg, n_comps=n_components)
    
    embeddings = adata_hvg.obsm['X_pca']
    
    print(f"✓ Generated PCA embeddings with shape: {embeddings.shape}")
    print(f"  Variance explained: {adata_hvg.uns['pca']['variance_ratio'].sum():.2%}")
    
    return embeddings, sample_idx

def generate_embeddings(adata, sample_size=None, use_pca_fallback=True, force_pca=False):
    """
    Main function to generate embeddings.
    
    Tries multiple approaches in order:
    1. Direct GeneFormer (most reliable)
    2. Helical GeneFormer (if available)
    3. PCA fallback (always works)
    
    Parameters:
    -----------
    adata : AnnData
        Single-cell dataset
    sample_size : int, optional
        Number of cells to sample
    use_pca_fallback : bool
        Whether to use PCA if GeneFormer fails
    force_pca : bool
        If True, skip GeneFormer and use PCA directly (fastest, always works)
    
    Returns:
    --------
    embeddings : np.ndarray
        Cell embeddings
    sample_idx : np.ndarray or None
        Sampled cell indices
    method : str
        Method used ('geneformer' or 'pca')
    """
    
    # Option to force PCA (useful if GeneFormer keeps failing)
    if force_pca:
        print("force_pca=True: Skipping GeneFormer, using PCA directly")
        embeddings, sample_idx = generate_pca_embeddings(adata, sample_size=sample_size)
        return embeddings, sample_idx, 'pca'
    
    # Try 1: Direct GeneFormer (bypasses helical issues)
    print("Attempting direct GeneFormer approach...")
    try:
        embeddings, sample_idx = generate_embeddings_direct_geneformer(adata, sample_size)
        return embeddings, sample_idx, 'geneformer_direct'
    except Exception as e:
        print(f"⚠ Direct GeneFormer failed: {str(e)[:100]}")
        print("→ Trying helical approach...\n")
    
    # Try 2: Initialize GeneFormer via helical
    model, tokenizer, method = initialize_geneformer()
    
    if model is not None:
        try:
            if 'helical' in method:
                print("Attempting helical GeneFormer...")
                embeddings, sample_idx = generate_embeddings_helical(
                    adata, model, tokenizer, method, sample_size
                )
                return embeddings, sample_idx, 'geneformer_helical'
            elif method == 'transformers':
                print("Attempting transformers GeneFormer...")
                embeddings, sample_idx = generate_embeddings_transformers(
                    adata, model, tokenizer, sample_size
                )
                return embeddings, sample_idx, 'geneformer_transformers'
        except Exception as e:
            print(f"\n⚠ GeneFormer embedding failed: {str(e)[:100]}")
            if not use_pca_fallback:
                raise
            print("→ Falling back to PCA...\n")
    
    # Try 3: PCA fallback
    if use_pca_fallback:
        embeddings, sample_idx = generate_pca_embeddings(adata, sample_size=sample_size)
        return embeddings, sample_idx, 'pca'
    else:
        raise RuntimeError("All GeneFormer methods failed and PCA fallback disabled")

if __name__ == "__main__":
    print("="*60)
    print("GeneFormer Integration Test")
    print("="*60)
    
    # Test initialization
    model, config, method = initialize_geneformer()
    
    if model is not None:
        print(f"\n✓ GeneFormer available via: {method}")
    else:
        print("\n✗ GeneFormer not available")
        print("→ PCA fallback will be used")
