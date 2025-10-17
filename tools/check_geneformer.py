"""
GeneFormer Setup and Testing Script

This script helps verify GeneFormer installation and tests basic functionality.
"""

import sys
import subprocess

def check_package(package_name):
    """Check if a package is installed."""
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False

def install_package(package_name):
    """Install a package using pip."""
    print(f"Installing {package_name}...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", package_name, "--break-system-packages"])

def check_geneformer_setup():
    """
    Check GeneFormer setup and dependencies.
    """
    print("="*60)
    print("GeneFormer Setup Check")
    print("="*60)
    
    required_packages = {
        'numpy': 'numpy',
        'pandas': 'pandas',
        'torch': 'torch',
        'scanpy': 'scanpy',
    }
    
    optional_packages = {
        'helical': 'helical',
        'transformers': 'transformers',
    }
    
    print("\n1. Checking Required Packages:")
    print("-" * 60)
    all_required = True
    for display_name, package_name in required_packages.items():
        if check_package(package_name):
            print(f"✓ {display_name:15} - Installed")
        else:
            print(f"✗ {display_name:15} - NOT INSTALLED")
            all_required = False
    
    print("\n2. Checking Optional Packages (GeneFormer):")
    print("-" * 60)
    helical_available = check_package('helical')
    transformers_available = check_package('transformers')
    
    if helical_available:
        print(f"✓ {'helical':15} - Installed")
    else:
        print(f"✗ {'helical':15} - NOT INSTALLED")
    
    if transformers_available:
        print(f"✓ {'transformers':15} - Installed")
    else:
        print(f"✗ {'transformers':15} - NOT INSTALLED")
    
    print("\n3. Testing GeneFormer Import Methods:")
    print("-" * 60)
    
    # Test Method 1: Helical new API
    try:
        from helical import Geneformer, GeneformerConfig
        print("✓ Method 1: from helical import Geneformer")
        method1_works = True
    except ImportError as e:
        print(f"✗ Method 1 failed: {str(e)[:50]}")
        method1_works = False
    
    # Test Method 2: Helical models API
    try:
        from helical.models.geneformer.model import Geneformer
        from helical.models.geneformer.geneformer_config import GeneformerConfig
        print("✓ Method 2: from helical.models.geneformer")
        method2_works = True
    except ImportError as e:
        print(f"✗ Method 2 failed: {str(e)[:50]}")
        method2_works = False
    
    # Test Method 3: Direct transformers
    try:
        from transformers import AutoModel, AutoTokenizer
        print("✓ Method 3: from transformers (direct)")
        method3_works = True
    except ImportError as e:
        print(f"✗ Method 3 failed: {str(e)[:50]}")
        method3_works = False
    
    geneformer_available = method1_works or method2_works or method3_works
    
    print("\n4. PyTorch Configuration:")
    print("-" * 60)
    if check_package('torch'):
        import torch
        print(f"PyTorch version: {torch.__version__}")
        print(f"CUDA available: {torch.cuda.is_available()}")
        if torch.cuda.is_available():
            print(f"CUDA version: {torch.version.cuda}")
            print(f"GPU devices: {torch.cuda.device_count()}")
            for i in range(torch.cuda.device_count()):
                print(f"  - Device {i}: {torch.cuda.get_device_name(i)}")
        else:
            print("Note: Running on CPU. GPU recommended for faster processing.")
    
    print("\n5. GeneFormer Status:")
    print("-" * 60)
    if geneformer_available:
        print("✓ GeneFormer can be imported")
        if method1_works:
            print("  Primary method: helical (new API)")
        elif method2_works:
            print("  Primary method: helical.models")
        elif method3_works:
            print("  Primary method: transformers (direct)")
    else:
        print("✗ GeneFormer cannot be imported")
        print("\nTo install GeneFormer:")
        if not helical_available:
            print("  pip install helical")
        if not transformers_available:
            print("  pip install transformers")
    
    print("\n6. Recommendations:")
    print("-" * 60)
    if not all_required:
        print("⚠ Install missing required packages:")
        for display_name, package_name in required_packages.items():
            if not check_package(package_name):
                print(f"  pip install {package_name}")
    
    if not geneformer_available:
        print("\n⚠ For GeneFormer functionality, install ONE of:")
        print("  Option 1 (Recommended): pip install helical")
        print("  Option 2 (Alternative): pip install transformers")
        print("\nNote: The notebook includes automatic PCA fallback")
        print("      if GeneFormer is unavailable")
    
    if all_required and geneformer_available:
        print("✓ All systems ready! You can proceed with Task 2.")
    elif all_required:
        print("✓ Required packages ready.")
        print("  GeneFormer unavailable but PCA fallback will work.")
    
    print("\n" + "="*60)
    
    return geneformer_available

def test_geneformer_loading():
    """
    Test loading GeneFormer model using helper module.
    """
    print("\n" + "="*60)
    print("Testing GeneFormer Model Loading")
    print("="*60)
    
    try:
        # Try to import and use the helper module
        import sys
        sys.path.insert(0, '.')
        from geneformer_helper import initialize_geneformer
        
        print("\nAttempting to load GeneFormer model...")
        print("Note: First load may take time to download model (~400MB)\n")
        
        model, tokenizer, method = initialize_geneformer()
        
        if model is not None:
            print(f"\n✓ GeneFormer model loaded successfully via {method}!")
            return True
        else:
            print("\n✗ Could not load GeneFormer with any method")
            return False
        
    except ImportError:
        print("✗ Cannot import geneformer_helper module")
        print("  Make sure geneformer_helper.py is in the same directory")
        return False
    except Exception as e:
        print(f"✗ Error loading GeneFormer: {e}")
        print("\nTroubleshooting:")
        print("  1. Check internet connection (model downloads from HuggingFace)")
        print("  2. Ensure sufficient disk space (~1GB)")
        print("  3. Try: pip install --upgrade helical")
        print("  4. Alternative: pip install transformers")
        return False

def create_test_data():
    """
    Create minimal test data for GeneFormer.
    """
    print("\n" + "="*60)
    print("Creating Test Data")
    print("="*60)
    
    try:
        import scanpy as sc
        import numpy as np
        
        # Create minimal test dataset
        n_cells = 100
        n_genes = 1000
        
        X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
        
        adata = sc.AnnData(X)
        adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
        adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]
        
        print(f"✓ Created test dataset: {n_cells} cells × {n_genes} genes")
        
        return adata
        
    except Exception as e:
        print(f"✗ Error creating test data: {e}")
        return None

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Check GeneFormer setup')
    parser.add_argument('--test-load', action='store_true', 
                       help='Test loading GeneFormer model (requires download)')
    parser.add_argument('--install', action='store_true',
                       help='Attempt to install missing packages')
    
    args = parser.parse_args()
    
    # Run basic checks
    geneformer_available = check_geneformer_setup()
    
    # Optionally test loading
    if args.test_load:
        success = test_geneformer_loading()
        if success:
            print("\n✓ GeneFormer is fully operational!")
        else:
            print("\n⚠ GeneFormer loading failed - will use PCA fallback")
    
    # Optionally install packages
    if args.install:
        print("\n" + "="*60)
        print("Installing Missing Packages")
        print("="*60)
        
        packages_to_install = ['scanpy', 'torch', 'helical']
        for package in packages_to_install:
            if not check_package(package):
                try:
                    install_package(package)
                    print(f"✓ {package} installed")
                except Exception as e:
                    print(f"✗ Failed to install {package}: {e}")
    
    # Final summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    if geneformer_available:
        print("✓ GeneFormer is available and ready to use")
    else:
        print("ℹ GeneFormer not available")
        print("  The notebook will automatically use PCA fallback")
        print("  This still works well for the challenge!")
