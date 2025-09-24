#!/usr/bin/env python3
"""
Demonstration script for Alabama River Fish Population Genomics Python tools.

This script demonstrates the usage of all three main analysis scripts with
example data, showing how to perform the complete workflow.

Usage:
    python demo_scripts.py
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

def create_example_ne_data():
    """Create example NeEstimator data for demonstration."""
    print("Creating example NeEstimator data...")
    
    # Example data similar to what NeEstimator would output
    data = {
        'Population': ['Pop1', 'Pop2', 'Pop3', 'Pop4', 'Pop5'],
        'r2': [0.025, 0.030, 0.028, 0.032, 0.027],
        'exp_r2': [0.010, 0.012, 0.011, 0.013, 0.010],
        'EffDF': [50, 45, 52, 48, 55],
        'Ne_original': [150.5, 120.3, 135.7, 110.2, 145.8]
    }
    
    df = pd.DataFrame(data)
    df.to_csv('demo_neestimator_data.csv', index=False)
    print("✓ Created demo_neestimator_data.csv")
    return 'demo_neestimator_data.csv'

def create_example_genepop_data():
    """Create example GENEPOP data for demonstration."""
    print("Creating example GENEPOP data...")
    
    # Create a simple GENEPOP format file
    genepop_content = """Demo genetic data for Alabama River Fish
Locus1 Locus2 Locus3 Locus4 Locus5
POP
Ind1, 001002 003004 002001 001003 002002
Ind2, 001001 003003 002002 001001 002001
Ind3, 002002 004004 001001 002002 001001
Ind4, 001002 003003 002001 001002 002002
POP
Ind5, 002002 004003 001002 002001 001003
Ind6, 001001 003004 002001 001001 002002
Ind7, 002001 004004 001001 002002 002001
Ind8, 001002 003003 002002 001003 001001
POP
Ind9, 001001 003004 002001 001001 002001
Ind10, 002002 004003 001002 002002 001003
"""
    
    with open('demo_genetic_data.gen', 'w') as f:
        f.write(genepop_content)
    
    print("✓ Created demo_genetic_data.gen")
    return 'demo_genetic_data.gen'

def demo_ne_correction():
    """Demonstrate Ne correction functionality."""
    print("\n" + "="*60)
    print("DEMONSTRATING Ne CORRECTION")
    print("="*60)
    
    # Create example data
    data_file = create_example_ne_data()
    
    # Import and run Ne correction
    from correct_ne_chromosome_number import NeCorrector
    
    print(f"\nRunning Ne correction on {data_file}...")
    corrector = NeCorrector(data_file, chromosomes=24)
    
    # Process all steps
    results = corrector.process_all_steps()
    
    if results is not None:
        # Print summary
        corrector.print_summary()
        
        # Save results
        output_file = 'demo_ne_corrected.csv'
        corrector.save_results(output_file)
        
        print(f"\n✓ Demo complete! Results saved to {output_file}")
        print("First few rows of corrected data:")
        display_cols = ['Ne_new', 'NeCorrected', 'lowCIcorrected', 'highCIcorrected']
        available_cols = [col for col in display_cols if col in results.columns]
        if available_cols:
            print(results[available_cols].head(3))
    
    return results

def demo_amova_analysis():
    """Demonstrate AMOVA analysis functionality."""
    print("\n" + "="*60)
    print("DEMONSTRATING AMOVA ANALYSIS")
    print("="*60)
    
    # Create example genetic data
    genepop_file = create_example_genepop_data()
    
    # Import AMOVA classes
    from amova_al_river_fish import GenePopReader, AMOVAAnalysis
    
    print(f"\nReading genetic data from {genepop_file}...")
    
    # Read genetic data
    reader = GenePopReader(genepop_file)
    genetic_data = reader.read_genepop()
    
    print(f"✓ Loaded genetic data:")
    print(f"  - {len(genetic_data['individuals'])} individuals")
    print(f"  - {len(set(genetic_data['populations']))} populations")
    print(f"  - {len(genetic_data['loci'])} loci")
    
    # Perform AMOVA
    print("\nPerforming AMOVA analysis...")
    amova = AMOVAAnalysis(genetic_data['data'], genetic_data['populations'])
    amova_results = amova.perform_amova()
    
    print("✓ AMOVA Results:")
    print(f"  - Among population variance: {amova_results['among_pop_variance']:.6f}")
    print(f"  - Within population variance: {amova_results['within_pop_variance']:.6f}")
    print(f"  - F_ST: {amova_results['fst']:.6f}")
    print(f"  - Number of populations: {amova_results['n_populations']}")
    print(f"  - Number of individuals: {amova_results['n_individuals']}")
    
    # Perform randomization test (smaller number for demo)
    print("\nPerforming randomization test (100 permutations for demo)...")
    rand_results = amova.randomization_test(n_permutations=100)
    print(f"✓ P-value: {rand_results['p_value']:.4f}")
    
    print("\n✓ AMOVA demo complete!")
    return amova_results, rand_results

def demo_dapc_analysis():
    """Demonstrate DAPC analysis functionality."""
    print("\n" + "="*60)
    print("DEMONSTRATING DAPC ANALYSIS")
    print("="*60)
    
    # Use the same genetic data file
    genepop_file = 'demo_genetic_data.gen'
    
    # Import DAPC classes
    from dapc_al_river_fish import GenePopReader, DAPCAnalysis
    
    print(f"\nReading genetic data from {genepop_file}...")
    
    # Read genetic data
    reader = GenePopReader(genepop_file)
    genetic_data = reader.read_genepop()
    
    print(f"✓ Loaded genetic data for DAPC:")
    print(f"  - {len(genetic_data['individuals'])} individuals")
    print(f"  - {len(set(genetic_data['populations']))} populations")
    
    # Perform DAPC analysis
    print("\nPerforming DAPC analysis...")
    dapc = DAPCAnalysis(genetic_data['genetic_matrix'], 
                       genetic_data['populations'],
                       n_pca_components=5)  # Small number for demo
    
    # Find clusters
    print("Finding optimal clusters...")
    cluster_results = dapc.find_clusters(max_n_clust=3)  # Small number for demo
    print(f"✓ Optimal clusters: {cluster_results['optimal_clusters']}")
    print(f"✓ Silhouette score: {max(cluster_results['silhouette_scores']):.4f}")
    
    # Perform DAPC
    print("Performing DAPC...")
    dapc_results = dapc.perform_dapc(n_pca_components=3)
    
    print(f"✓ DAPC Results:")
    print(f"  - PCA components: {dapc_results['n_pca_components']}")
    print(f"  - LDA components: {dapc_results['n_lda_components']}")
    print(f"  - PCA variance explained: {dapc_results['explained_variance_pca'].sum():.4f}")
    
    print("\n✓ DAPC demo complete!")
    return dapc_results

def cleanup_demo_files():
    """Clean up demonstration files."""
    demo_files = [
        'demo_neestimator_data.csv',
        'demo_genetic_data.gen',
        'demo_ne_corrected.csv'
    ]
    
    print(f"\nCleaning up demo files...")
    for file in demo_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"✓ Removed {file}")

def main():
    """Run all demonstrations."""
    print("Alabama River Fish Population Genomics - Python Tools Demo")
    print("="*70)
    print("This demo shows the basic functionality of all three analysis scripts.")
    print("It uses small example datasets to demonstrate the workflow.\n")
    
    try:
        # Demo 1: Ne correction
        ne_results = demo_ne_correction()
        
        # Demo 2: AMOVA analysis
        amova_results, rand_results = demo_amova_analysis()
        
        # Demo 3: DAPC analysis  
        dapc_results = demo_dapc_analysis()
        
        print("\n" + "="*70)
        print("ALL DEMOS COMPLETED SUCCESSFULLY!")
        print("="*70)
        print("The Python scripts are working correctly and ready to use")
        print("with your actual genetic data files.")
        print("\nFor more information, see the README.md file or run:")
        print("  python amova_al_river_fish.py --help")
        print("  python correct_ne_chromosome_number.py --help")
        print("  python dapc_al_river_fish.py --help")
        
    except Exception as e:
        print(f"\nError during demo: {e}")
        print("Please check that all required packages are installed:")
        print("  pip install -r requirements.txt")
        
    finally:
        # Clean up demo files
        cleanup_demo_files()

if __name__ == "__main__":
    main()