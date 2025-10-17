#!/usr/bin/env python3
"""
Task 3 Prerequisites Checker
Verifies that all required files from Task 2 are present and correctly formatted
"""

import sys
from pathlib import Path

def check_task3_prerequisites():
    """Check if all Task 3 prerequisites are met"""
    
    print("=" * 70)
    print("                TASK 3 PREREQUISITES CHECKER")
    print("=" * 70)
    print()
    
    all_good = True
    
    # Check 1: Python packages
    print("📦 Checking Python packages...")
    required_packages = [
        'scanpy', 'anndata', 'numpy', 'pandas', 
        'matplotlib', 'seaborn', 'scipy', 
        'sklearn', 'umap', 'statsmodels'
    ]
    
    missing_packages = []
    for package in required_packages:
        try:
            if package == 'sklearn':
                __import__('sklearn')
            else:
                __import__(package)
            print(f"   ✅ {package}")
        except ImportError:
            print(f"   ❌ {package} - NOT INSTALLED")
            missing_packages.append(package)
            all_good = False
    
    if missing_packages:
        print(f"\n⚠️  Install missing packages with:")
        print(f"   pip install {' '.join(missing_packages)}")
    print()
    
    # Check 2: Original data file
    print("📁 Checking original data file...")
    data_file = Path('localtools/counts_combined_filtered_BA4_sALS_PN.h5ad')
    if data_file.exists():
        print(f"   ✅ {data_file} found")
        print(f"      Size: {data_file.stat().st_size / (1024**2):.1f} MB")
    else:
        print(f"   ❌ {data_file} NOT FOUND")
        print(f"      This should be your original ALS dataset")
        all_good = False
    print()
    
    # Check 3: Task 2 output files
    print("📂 Checking Task 2 output files...")
    scenario_files = [
        'third attempt/perturbed_scenario1_disease_rescue.h5ad',
        'third attempt/perturbed_scenario2_disease_rescue.h5ad',
        'third attempt/perturbed_scenario3_disease_rescue.h5ad'
    ]
    
    found_files = []
    for file_path in scenario_files:
        p = Path(file_path)
        if p.exists():
            print(f"   ✅ {p.name} found")
            print(f"      Size: {p.stat().st_size / (1024**2):.1f} MB")
            found_files.append(file_path)
        else:
            print(f"   ❌ {p.name} NOT FOUND")
            all_good = False
    print()
    
    # Check 4: Verify file contents (if scanpy available)
    if found_files and 'scanpy' not in missing_packages:
        print("🔍 Verifying file contents...")
        import scanpy as sc
        import numpy as np
        
        for file_path in found_files:
            print(f"\n   Checking: {Path(file_path).name}")
            try:
                adata = sc.read_h5ad(file_path)
                print(f"      ✅ Loads successfully")
                print(f"      📊 Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} features")
                
                # Check for embeddings
                if adata.X.shape[1] == 256:
                    print(f"      ✅ Contains GeneFormer embeddings (256-dim) in X")
                elif 'X_geneformer' in adata.obsm:
                    emb_shape = adata.obsm['X_geneformer'].shape
                    print(f"      ✅ Contains embeddings in obsm: {emb_shape}")
                elif any('embed' in k.lower() or 'geneformer' in k.lower() for k in adata.obsm.keys()):
                    print(f"      ✅ Contains embeddings in obsm")
                    print(f"         Keys: {list(adata.obsm.keys())}")
                else:
                    print(f"      ⚠️  No obvious embeddings found")
                    print(f"         X shape: {adata.X.shape}")
                    print(f"         obsm keys: {list(adata.obsm.keys())}")
                    print(f"      This might still work - Task 3 will attempt to find them")
                
                # Check metadata
                if 'CellType' in adata.obs.columns:
                    print(f"      ✅ Contains CellType metadata ({adata.obs['CellType'].nunique()} types)")
                if 'Condition' in adata.obs.columns:
                    conditions = adata.obs['Condition'].value_counts().to_dict()
                    print(f"      ✅ Contains Condition metadata: {conditions}")
                    
            except Exception as e:
                print(f"      ❌ Error loading file: {e}")
                all_good = False
        print()
    
    # Check 5: System resources
    print("💻 Checking system resources...")
    try:
        import psutil
        
        mem = psutil.virtual_memory()
        mem_gb = mem.total / (1024**3)
        mem_avail_gb = mem.available / (1024**3)
        
        print(f"   Total RAM: {mem_gb:.1f} GB")
        print(f"   Available RAM: {mem_avail_gb:.1f} GB")
        
        if mem_gb < 16:
            print(f"   ⚠️  Less than 16 GB RAM - consider subsampling")
        elif mem_avail_gb < 8:
            print(f"   ⚠️  Less than 8 GB available - close other programs")
        else:
            print(f"   ✅ Sufficient memory available")
        
        # Check disk space
        disk = psutil.disk_usage('.')
        disk_free_gb = disk.free / (1024**3)
        print(f"   Free disk space: {disk_free_gb:.1f} GB")
        
        if disk_free_gb < 2:
            print(f"   ⚠️  Less than 2 GB free - may need more space for outputs")
        else:
            print(f"   ✅ Sufficient disk space")
            
    except ImportError:
        print(f"   ℹ️  Install psutil to check system resources: pip install psutil")
    except Exception as e:
        print(f"   ⚠️  Could not check resources: {e}")
    print()
    
    # Final summary
    print("=" * 70)
    if all_good:
        print("✅ ALL CHECKS PASSED - Ready to run Task 3!")
        print()
        print("Next steps:")
        print("   1. jupyter notebook task3_embedding_analysis.ipynb")
        print("   2. Cell → Run All")
        print("   3. Wait ~45-60 minutes")
        print("   4. Check outputs/ folder for results")
    else:
        print("❌ SOME CHECKS FAILED - Please address issues above")
        print()
        print("Common solutions:")
        print("   • Install missing packages: pip install scanpy anndata ...")
        print("   • Ensure Task 2 completed successfully")
        print("   • Check file paths are correct")
        print("   • Verify you're in the correct directory")
    print("=" * 70)
    
    return all_good

if __name__ == "__main__":
    success = check_task3_prerequisites()
    sys.exit(0 if success else 1)
