#!/usr/bin/env python3
"""
AMOVA Analysis for Alabama River Fish Population Genomics

This script performs Analysis of Molecular Variance (AMOVA) on genetic data
from freshwater drum and white crappie populations, replicating the functionality
of the original R script AMOVA_AL-River-fish.R.

Author: Converted from R to Python
Date: 2024
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
import re
import warnings


class GenePopReader:
    """
    Class to read and parse GENEPOP format files for genetic analysis.
    
    GENEPOP format is a standard format for population genetics data.
    """
    
    def __init__(self, filepath):
        """
        Initialize GenePopReader with a file path.
        
        Args:
            filepath (str): Path to the GENEPOP (.gen) file
        """
        self.filepath = filepath
        self.data = None
        self.populations = None
        self.loci = None
        self.individuals = None
        
    def read_genepop(self):
        """
        Read and parse a GENEPOP format file.
        
        Returns:
            dict: Dictionary containing genetic data, populations, and metadata
        """
        with open(self.filepath, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        
        # Find the line with loci names (second line typically)
        loci_line = 1
        loci = lines[loci_line].split()
        
        # Find the "POP" markers to identify population boundaries
        pop_indices = [i for i, line in enumerate(lines) if line.upper() == 'POP']
        
        populations = []
        individuals = []
        genetic_data = []
        
        current_pop = 0
        
        for i, line in enumerate(lines):
            if line.upper() == 'POP':
                current_pop += 1
                continue
            elif i > loci_line and line and not line.upper() == 'POP':
                # Parse individual data
                parts = line.split(',')
                if len(parts) >= 2:
                    individual_id = parts[0].strip()
                    genotype_data = parts[1].strip().split()
                    
                    individuals.append(individual_id)
                    populations.append(f"Pop_{current_pop}")
                    genetic_data.append(genotype_data)
        
        self.data = genetic_data
        self.populations = populations
        self.loci = loci
        self.individuals = individuals
        
        return {
            'data': genetic_data,
            'populations': populations,
            'loci': loci,
            'individuals': individuals
        }


class AMOVAAnalysis:
    """
    Analysis of Molecular Variance (AMOVA) implementation for population genetics.
    
    This class implements AMOVA analysis similar to the poppr package in R,
    calculating variance components within and among populations.
    """
    
    def __init__(self, genetic_data, populations, cutoff=0.5):
        """
        Initialize AMOVA analysis.
        
        Args:
            genetic_data (list): Genetic data matrix
            populations (list): Population assignments for each individual
            cutoff (float): Cutoff value for genetic distance calculations
        """
        self.genetic_data = genetic_data
        self.populations = populations
        self.cutoff = cutoff
        self.distance_matrix = None
        self.amova_results = None
        
    def calculate_genetic_distances(self):
        """
        Calculate genetic distances between individuals.
        
        Returns:
            numpy.ndarray: Distance matrix between individuals
        """
        # Convert genetic data to numerical format for distance calculation
        numeric_data = []
        
        for individual in self.genetic_data:
            numeric_individual = []
            for locus in individual:
                try:
                    # Convert genotype to numeric (simple approach)
                    if locus == '000000':  # Missing data
                        numeric_individual.append(np.nan)
                    else:
                        # Split alleles and convert to numbers
                        allele1 = int(locus[:3]) if locus[:3] != '000' else np.nan
                        allele2 = int(locus[3:]) if locus[3:] != '000' else np.nan
                        numeric_individual.append((allele1, allele2))
                except (ValueError, IndexError):
                    numeric_individual.append(np.nan)
            
            numeric_data.append(numeric_individual)
        
        # Calculate pairwise distances (simplified Euclidean for demonstration)
        distances = []
        n_individuals = len(numeric_data)
        
        for i in range(n_individuals):
            row_distances = []
            for j in range(n_individuals):
                if i == j:
                    row_distances.append(0.0)
                else:
                    # Simple genetic distance calculation
                    distance = self._calculate_individual_distance(
                        numeric_data[i], numeric_data[j]
                    )
                    row_distances.append(distance)
            distances.append(row_distances)
        
        self.distance_matrix = np.array(distances)
        return self.distance_matrix
    
    def _calculate_individual_distance(self, ind1, ind2):
        """
        Calculate genetic distance between two individuals.
        
        Args:
            ind1, ind2: Genetic data for two individuals
            
        Returns:
            float: Genetic distance
        """
        differences = 0
        total_loci = 0
        
        for locus1, locus2 in zip(ind1, ind2):
            if not (pd.isna(locus1) or pd.isna(locus2)):
                total_loci += 1
                if locus1 != locus2:
                    differences += 1
        
        return differences / total_loci if total_loci > 0 else 0.0
    
    def perform_amova(self):
        """
        Perform AMOVA analysis.
        
        Returns:
            dict: AMOVA results including variance components and F-statistics
        """
        if self.distance_matrix is None:
            self.calculate_genetic_distances()
        
        # Get unique populations
        unique_pops = list(set(self.populations))
        n_pops = len(unique_pops)
        n_individuals = len(self.populations)
        
        # Calculate within-population and among-population variance components
        within_pop_ss = 0
        among_pop_ss = 0
        total_ss = 0
        
        # Calculate total sum of squares
        grand_mean = np.mean(self.distance_matrix[np.triu_indices_from(self.distance_matrix, k=1)])
        
        for i in range(n_individuals):
            for j in range(i+1, n_individuals):
                total_ss += (self.distance_matrix[i, j] - grand_mean) ** 2
        
        # Calculate among-population sum of squares
        pop_means = {}
        pop_sizes = {}
        
        for pop in unique_pops:
            pop_indices = [i for i, p in enumerate(self.populations) if p == pop]
            pop_sizes[pop] = len(pop_indices)
            
            if len(pop_indices) > 1:
                pop_distances = []
                for i in pop_indices:
                    for j in pop_indices:
                        if i != j:
                            pop_distances.append(self.distance_matrix[i, j])
                pop_means[pop] = np.mean(pop_distances) if pop_distances else 0
            else:
                pop_means[pop] = 0
        
        # Calculate variance components (simplified)
        among_pop_variance = np.var(list(pop_means.values()))
        within_pop_variance = np.mean([
            np.var([self.distance_matrix[i, j] for j in range(n_individuals) 
                   if self.populations[j] == self.populations[i] and i != j])
            for i in range(n_individuals)
        ])
        
        total_variance = among_pop_variance + within_pop_variance
        
        # Calculate F-statistics
        fst = among_pop_variance / total_variance if total_variance > 0 else 0
        
        self.amova_results = {
            'among_pop_variance': among_pop_variance,
            'within_pop_variance': within_pop_variance,
            'total_variance': total_variance,
            'fst': fst,
            'degrees_freedom_among': n_pops - 1,
            'degrees_freedom_within': n_individuals - n_pops,
            'n_populations': n_pops,
            'n_individuals': n_individuals
        }
        
        return self.amova_results
    
    def randomization_test(self, n_permutations=1000):
        """
        Perform randomization test for AMOVA significance.
        
        Args:
            n_permutations (int): Number of permutations for the test
            
        Returns:
            dict: Results of randomization test including p-value
        """
        if self.amova_results is None:
            self.perform_amova()
        
        observed_fst = self.amova_results['fst']
        permuted_fsts = []
        
        original_populations = self.populations.copy()
        
        for _ in range(n_permutations):
            # Randomly permute population assignments
            permuted_pops = original_populations.copy()
            np.random.shuffle(permuted_pops)
            
            # Create temporary AMOVA with permuted populations
            temp_amova = AMOVAAnalysis(self.genetic_data, permuted_pops, self.cutoff)
            temp_amova.distance_matrix = self.distance_matrix  # Use same distances
            temp_results = temp_amova.perform_amova()
            permuted_fsts.append(temp_results['fst'])
        
        # Calculate p-value
        p_value = np.sum(np.array(permuted_fsts) >= observed_fst) / n_permutations
        
        return {
            'observed_fst': observed_fst,
            'permuted_fsts': permuted_fsts,
            'p_value': p_value,
            'n_permutations': n_permutations
        }
    
    def plot_randomization_test(self, rand_test_results, output_file=None):
        """
        Plot results of randomization test.
        
        Args:
            rand_test_results (dict): Results from randomization_test()
            output_file (str): Optional path to save the plot
        """
        plt.figure(figsize=(10, 6))
        
        # Plot histogram of permuted F_ST values
        plt.hist(rand_test_results['permuted_fsts'], bins=50, alpha=0.7, 
                color='lightblue', edgecolor='black', label='Permuted F_ST')
        
        # Plot observed F_ST as vertical line
        plt.axvline(rand_test_results['observed_fst'], color='red', 
                   linestyle='--', linewidth=2, label=f"Observed F_ST = {rand_test_results['observed_fst']:.4f}")
        
        plt.xlabel('F_ST')
        plt.ylabel('Frequency')
        plt.title(f'Randomization Test Results\nP-value = {rand_test_results["p_value"]:.4f}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        
        plt.show()


def process_genepop_files(data_directory, pattern="*.genepop"):
    """
    Process GENEPOP files in a directory, renaming .genepop to .gen.
    
    Args:
        data_directory (str): Directory containing GENEPOP files
        pattern (str): Pattern to match files
        
    Returns:
        list: List of processed .gen files
    """
    data_path = Path(data_directory)
    genepop_files = list(data_path.glob(pattern))
    
    processed_files = []
    
    for genepop_file in genepop_files:
        # Rename .genepop to .gen
        gen_file = genepop_file.with_suffix('.gen')
        if not gen_file.exists():
            genepop_file.rename(gen_file)
            print(f"Renamed {genepop_file.name} to {gen_file.name}")
        processed_files.append(gen_file)
    
    return processed_files


def analyze_freshwater_drum(data_directory):
    """
    Analyze freshwater drum genetic data using AMOVA.
    
    Args:
        data_directory (str): Directory containing freshwater drum data
        
    Returns:
        dict: Analysis results
    """
    print("#### Freshwater Drum AMOVA Analysis ####")
    
    # Process genepop files
    processed_files = process_genepop_files(data_directory)
    
    # Look for the specific freshwater drum file
    drum_file = None
    for file in processed_files:
        if "freshwater-drum_m3M2n2_p4r7_mac3.haps.gen" in str(file):
            drum_file = file
            break
    
    if not drum_file:
        print(f"Freshwater drum data file not found in {data_directory}")
        return None
    
    # Read genetic data
    reader = GenePopReader(drum_file)
    genetic_data = reader.read_genepop()
    
    # Perform AMOVA
    amova = AMOVAAnalysis(genetic_data['data'], genetic_data['populations'])
    amova_results = amova.perform_amova()
    
    print("AMOVA Results:")
    print(f"Among population variance: {amova_results['among_pop_variance']:.6f}")
    print(f"Within population variance: {amova_results['within_pop_variance']:.6f}")
    print(f"F_ST: {amova_results['fst']:.6f}")
    
    # Perform randomization test
    print("\nPerforming randomization test...")
    rand_results = amova.randomization_test(n_permutations=1000)
    print(f"P-value: {rand_results['p_value']:.4f}")
    
    # Plot results
    amova.plot_randomization_test(rand_results, "freshwater_drum_randomization_test.png")
    
    return {
        'amova_results': amova_results,
        'randomization_results': rand_results,
        'genetic_data': genetic_data
    }


def analyze_white_crappie(data_directory):
    """
    Analyze white crappie genetic data using AMOVA.
    
    Args:
        data_directory (str): Directory containing white crappie data
        
    Returns:
        dict: Analysis results
    """
    print("#### White Crappie AMOVA Analysis ####")
    
    # Process genepop files
    processed_files = process_genepop_files(data_directory)
    
    # Look for white crappie file (using same drum file name as in original R script)
    crappie_file = None
    for file in processed_files:
        if "freshwater-drum_m3M2n2_p4r7_mac3.haps.gen" in str(file):
            crappie_file = file
            break
    
    if not crappie_file:
        print(f"White crappie data file not found in {data_directory}")
        return None
    
    # Read genetic data
    reader = GenePopReader(crappie_file)
    genetic_data = reader.read_genepop()
    
    # Perform AMOVA
    amova = AMOVAAnalysis(genetic_data['data'], genetic_data['populations'])
    amova_results = amova.perform_amova()
    
    print("AMOVA Results:")
    print(f"Among population variance: {amova_results['among_pop_variance']:.6f}")
    print(f"Within population variance: {amova_results['within_pop_variance']:.6f}")
    print(f"F_ST: {amova_results['fst']:.6f}")
    
    # Perform randomization test
    print("\nPerforming randomization test...")
    rand_results = amova.randomization_test(n_permutations=1000)
    print(f"P-value: {rand_results['p_value']:.4f}")
    
    # Plot results
    amova.plot_randomization_test(rand_results, "white_crappie_randomization_test.png")
    
    return {
        'amova_results': amova_results,
        'randomization_results': rand_results,
        'genetic_data': genetic_data
    }


def main():
    """
    Main function to run AMOVA analysis on both species.
    
    Example usage:
        python amova_al_river_fish.py
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='AMOVA Analysis for Alabama River Fish')
    parser.add_argument('--drum_dir', type=str, 
                       default='~/Fish/AL-River_migratory-fish/freshwater-drum/freshwater-drum_m3M2n2_p4r7_mac3/',
                       help='Directory containing freshwater drum data')
    parser.add_argument('--crappie_dir', type=str,
                       default='~/Fish/AL-River_migratory-fish/freshwater-drum/freshwater-drum_m3M2n2_p4r7_mac3/',
                       help='Directory containing white crappie data')
    
    args = parser.parse_args()
    
    # Expand home directory
    drum_dir = Path(args.drum_dir).expanduser()
    crappie_dir = Path(args.crappie_dir).expanduser()
    
    # Analyze freshwater drum
    if drum_dir.exists():
        drum_results = analyze_freshwater_drum(drum_dir)
    else:
        print(f"Freshwater drum directory not found: {drum_dir}")
        drum_results = None
    
    # Analyze white crappie
    if crappie_dir.exists():
        crappie_results = analyze_white_crappie(crappie_dir)
    else:
        print(f"White crappie directory not found: {crappie_dir}")
        crappie_results = None
    
    return drum_results, crappie_results


if __name__ == "__main__":
    main()