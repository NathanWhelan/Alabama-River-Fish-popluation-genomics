#!/usr/bin/env python3
"""
Discriminant Analysis of Principal Components (DAPC) for Alabama River Fish

This script performs DAPC analysis on genetic data from white crappie and 
freshwater drum populations, replicating the functionality of the original 
R script DAPC_AL-River-Fish.R.

DAPC combines PCA dimensionality reduction with Linear Discriminant Analysis
to identify genetic clusters and visualize population structure.

Author: Converted from R to Python
Date: 2024
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from pathlib import Path
import glob
import re
import warnings


class GenePopReader:
    """
    Class to read and parse GENEPOP format files for DAPC analysis.
    """
    
    def __init__(self, filepath):
        """
        Initialize GenePopReader.
        
        Args:
            filepath (str): Path to the GENEPOP (.gen) file
        """
        self.filepath = filepath
        self.data = None
        self.populations = None
        self.individuals = None
        self.genetic_matrix = None
        
    def read_genepop(self):
        """
        Read and parse GENEPOP file into a format suitable for DAPC.
        
        Returns:
            dict: Genetic data and metadata
        """
        with open(self.filepath, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        
        # Parse header and loci information
        title = lines[0] if lines else "Genetic Data"
        loci_line = 1
        loci = lines[loci_line].split() if len(lines) > 1 else []
        
        # Find population boundaries
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
                # Parse individual genetic data
                parts = line.split(',')
                if len(parts) >= 2:
                    individual_id = parts[0].strip()
                    genotype_str = parts[1].strip()
                    
                    # Convert genotypes to numerical format
                    genotype_numeric = self._parse_genotypes(genotype_str)
                    
                    individuals.append(individual_id)
                    populations.append(f"Pop_{current_pop}")
                    genetic_data.append(genotype_numeric)
        
        # Convert to matrix format
        self.genetic_matrix = np.array(genetic_data)
        self.populations = populations
        self.individuals = individuals
        
        return {
            'genetic_matrix': self.genetic_matrix,
            'populations': populations,
            'individuals': individuals,
            'loci': loci,
            'title': title
        }
    
    def _parse_genotypes(self, genotype_str):
        """
        Parse genotype string into numerical format for analysis.
        
        Args:
            genotype_str (str): Raw genotype string from GENEPOP file
            
        Returns:
            list: Numerical representation of genotypes
        """
        genotypes = genotype_str.split()
        numeric_genotypes = []
        
        for genotype in genotypes:
            if genotype == '000000' or genotype == '0000':  # Missing data
                numeric_genotypes.extend([0, 0])
            else:
                try:
                    if len(genotype) == 6:  # Format: AAABBB
                        allele1 = int(genotype[:3])
                        allele2 = int(genotype[3:])
                    elif len(genotype) == 4:  # Format: AABB
                        allele1 = int(genotype[:2])
                        allele2 = int(genotype[2:])
                    else:
                        allele1 = allele2 = 0
                    
                    numeric_genotypes.extend([allele1, allele2])
                except (ValueError, IndexError):
                    numeric_genotypes.extend([0, 0])
        
        return numeric_genotypes


class DAPCAnalysis:
    """
    Discriminant Analysis of Principal Components implementation.
    
    This class performs DAPC analysis by first applying PCA for dimensionality
    reduction, then using Linear Discriminant Analysis for classification.
    """
    
    def __init__(self, genetic_matrix, populations, n_pca_components=None, random_state=42):
        """
        Initialize DAPC analysis.
        
        Args:
            genetic_matrix (numpy.ndarray): Genetic data matrix
            populations (list): Population assignments
            n_pca_components (int): Number of PCA components to retain
            random_state (int): Random state for reproducibility
        """
        self.genetic_matrix = genetic_matrix
        self.populations = populations
        self.n_pca_components = n_pca_components
        self.random_state = random_state
        
        # Analysis components
        self.scaler = StandardScaler()
        self.pca = None
        self.lda = None
        self.clusters = None
        
        # Results
        self.pca_data = None
        self.dapc_data = None
        self.explained_variance_pca = None
        self.explained_variance_lda = None
        
    def find_clusters(self, max_n_clust=5, n_pca_components=250):
        """
        Find optimal number of clusters using K-means clustering.
        
        Args:
            max_n_clust (int): Maximum number of clusters to test
            n_pca_components (int): Number of PCA components for clustering
            
        Returns:
            dict: Clustering results
        """
        print(f"Finding optimal clusters (max: {max_n_clust})...")
        
        # Prepare data for clustering
        scaled_data = self.scaler.fit_transform(self.genetic_matrix)
        
        # Apply PCA for dimensionality reduction
        pca_for_clustering = PCA(n_components=min(n_pca_components, scaled_data.shape[1]))
        pca_data = pca_for_clustering.fit_transform(scaled_data)
        
        # Test different numbers of clusters
        silhouette_scores = []
        inertias = []
        cluster_range = range(2, max_n_clust + 1)
        
        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, random_state=self.random_state)
            cluster_labels = kmeans.fit_predict(pca_data)
            
            silhouette_avg = silhouette_score(pca_data, cluster_labels)
            silhouette_scores.append(silhouette_avg)
            inertias.append(kmeans.inertia_)
        
        # Find optimal number of clusters (highest silhouette score)
        optimal_idx = np.argmax(silhouette_scores)
        optimal_clusters = cluster_range[optimal_idx]
        
        print(f"Optimal number of clusters: {optimal_clusters}")
        print(f"Silhouette score: {silhouette_scores[optimal_idx]:.4f}")
        
        # Perform final clustering with optimal number
        kmeans = KMeans(n_clusters=optimal_clusters, random_state=self.random_state)
        self.clusters = kmeans.fit_predict(pca_data)
        
        return {
            'optimal_clusters': optimal_clusters,
            'cluster_labels': self.clusters,
            'silhouette_scores': silhouette_scores,
            'inertias': inertias,
            'cluster_range': list(cluster_range)
        }
    
    def perform_dapc(self, n_pca_components=None):
        """
        Perform DAPC analysis.
        
        Args:
            n_pca_components (int): Number of PCA components to retain
            
        Returns:
            dict: DAPC results
        """
        if self.clusters is None:
            print("Warning: No clusters found. Running find_clusters first...")
            self.find_clusters()
        
        print("Performing DAPC analysis...")
        
        # Scale the data
        scaled_data = self.scaler.fit_transform(self.genetic_matrix)
        
        # Step 1: PCA
        if n_pca_components is None:
            n_pca_components = min(self.n_pca_components or 50, scaled_data.shape[1])
        
        self.pca = PCA(n_components=n_pca_components)
        pca_data = self.pca.fit_transform(scaled_data)
        self.pca_data = pca_data
        self.explained_variance_pca = self.pca.explained_variance_ratio_
        
        print(f"PCA: {n_pca_components} components explain {self.explained_variance_pca.sum():.4f} of variance")
        
        # Step 2: Linear Discriminant Analysis
        n_classes = len(np.unique(self.clusters))
        n_lda_components = min(n_classes - 1, pca_data.shape[1])
        
        self.lda = LinearDiscriminantAnalysis(n_components=n_lda_components)
        dapc_data = self.lda.fit_transform(pca_data, self.clusters)
        self.dapc_data = dapc_data
        
        # Calculate explained variance for LDA components
        if hasattr(self.lda, 'explained_variance_ratio_'):
            self.explained_variance_lda = self.lda.explained_variance_ratio_
        else:
            # Calculate manually if not available
            total_var = np.sum(self.lda.coef_**2, axis=1)
            self.explained_variance_lda = total_var / total_var.sum()
        
        print(f"LDA: {n_lda_components} discriminant axes")
        
        return {
            'pca_data': pca_data,
            'dapc_data': dapc_data,
            'explained_variance_pca': self.explained_variance_pca,
            'explained_variance_lda': self.explained_variance_lda,
            'n_pca_components': n_pca_components,
            'n_lda_components': n_lda_components
        }
    
    def plot_dapc_results(self, title="DAPC Analysis", output_file=None, 
                         color_by_population=True, figsize=(15, 10)):
        """
        Create comprehensive DAPC visualization plots.
        
        Args:
            title (str): Plot title
            output_file (str): Optional path to save the plot
            color_by_population (bool): Whether to color by population or cluster
            figsize (tuple): Figure size
        """
        if self.dapc_data is None:
            print("Error: DAPC analysis not performed. Run perform_dapc() first.")
            return
        
        # Color scheme
        colors = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', 
                 '#0072B2', '#D55E00', '#CC79A7', 'grey', 'darkblue', 'green']
        
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # 1. Main DAPC scatter plot
        ax1 = fig.add_subplot(gs[0:2, 0:2])
        
        if color_by_population:
            # Color by original populations
            unique_pops = list(set(self.populations))
            pop_colors = {pop: colors[i % len(colors)] for i, pop in enumerate(unique_pops)}
            
            for pop in unique_pops:
                mask = np.array(self.populations) == pop
                if self.dapc_data.shape[1] >= 2:
                    ax1.scatter(self.dapc_data[mask, 0], self.dapc_data[mask, 1], 
                              c=pop_colors[pop], label=pop, alpha=0.7, s=50)
                else:
                    ax1.scatter(range(sum(mask)), self.dapc_data[mask, 0], 
                              c=pop_colors[pop], label=pop, alpha=0.7, s=50)
        else:
            # Color by clusters
            unique_clusters = np.unique(self.clusters)
            cluster_colors = {cluster: colors[i % len(colors)] for i, cluster in enumerate(unique_clusters)}
            
            for cluster in unique_clusters:
                mask = self.clusters == cluster
                if self.dapc_data.shape[1] >= 2:
                    ax1.scatter(self.dapc_data[mask, 0], self.dapc_data[mask, 1], 
                              c=cluster_colors[cluster], label=f'Cluster {cluster+1}', alpha=0.7, s=50)
                else:
                    ax1.scatter(range(sum(mask)), self.dapc_data[mask, 0], 
                              c=cluster_colors[cluster], label=f'Cluster {cluster+1}', alpha=0.7, s=50)
        
        if self.dapc_data.shape[1] >= 2:
            ax1.set_xlabel(f'LD1 ({self.explained_variance_lda[0]:.1%} of variance)')
            ax1.set_ylabel(f'LD2 ({self.explained_variance_lda[1]:.1%} of variance)')
        else:
            ax1.set_xlabel('Sample Index')
            ax1.set_ylabel(f'LD1 ({self.explained_variance_lda[0]:.1%} of variance)')
        
        ax1.set_title(title)
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        
        # 2. PCA scree plot
        ax2 = fig.add_subplot(gs[0, 2])
        n_components_to_show = min(20, len(self.explained_variance_pca))
        ax2.bar(range(1, n_components_to_show + 1), 
                self.explained_variance_pca[:n_components_to_show])
        ax2.set_xlabel('PC')
        ax2.set_ylabel('Variance Explained')
        ax2.set_title('PCA Scree Plot')
        ax2.grid(True, alpha=0.3)
        
        # 3. LDA scree plot
        ax3 = fig.add_subplot(gs[1, 2])
        if len(self.explained_variance_lda) > 1:
            ax3.bar(range(1, len(self.explained_variance_lda) + 1), 
                    self.explained_variance_lda)
        else:
            ax3.bar([1], self.explained_variance_lda)
        ax3.set_xlabel('LD')
        ax3.set_ylabel('Variance Explained')
        ax3.set_title('LDA Scree Plot')
        ax3.grid(True, alpha=0.3)
        
        # 4. Cumulative variance explained
        ax4 = fig.add_subplot(gs[2, :])
        cumvar_pca = np.cumsum(self.explained_variance_pca)
        ax4.plot(range(1, len(cumvar_pca) + 1), cumvar_pca, 'o-', label='PCA')
        
        if len(self.explained_variance_lda) > 1:
            cumvar_lda = np.cumsum(self.explained_variance_lda)
            ax4.plot(range(1, len(cumvar_lda) + 1), cumvar_lda, 's-', label='LDA')
        
        ax4.set_xlabel('Component')
        ax4.set_ylabel('Cumulative Variance Explained')
        ax4.set_title('Cumulative Variance Explained')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        ax4.set_ylim(0, 1)
        
        plt.suptitle(f'{title} - Comprehensive Results', fontsize=16)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        
        plt.show()
    
    def save_results(self, output_prefix):
        """
        Save DAPC results to files.
        
        Args:
            output_prefix (str): Prefix for output files
        """
        if self.dapc_data is None:
            print("Error: No results to save. Run perform_dapc() first.")
            return
        
        # Save DAPC coordinates
        dapc_df = pd.DataFrame(self.dapc_data, 
                              columns=[f'LD{i+1}' for i in range(self.dapc_data.shape[1])])
        dapc_df['Individual'] = self.individuals
        dapc_df['Population'] = self.populations
        dapc_df['Cluster'] = self.clusters + 1  # Convert to 1-based indexing
        
        dapc_file = f"{output_prefix}_dapc_coordinates.csv"
        dapc_df.to_csv(dapc_file, index=False)
        print(f"DAPC coordinates saved to {dapc_file}")
        
        # Save explained variance
        variance_df = pd.DataFrame({
            'Component': [f'PC{i+1}' for i in range(len(self.explained_variance_pca))],
            'PCA_Variance': self.explained_variance_pca
        })
        
        if len(self.explained_variance_lda) > 0:
            lda_variance = np.zeros(len(self.explained_variance_pca))
            lda_variance[:len(self.explained_variance_lda)] = self.explained_variance_lda
            variance_df['LDA_Variance'] = lda_variance
        
        variance_file = f"{output_prefix}_variance_explained.csv"
        variance_df.to_csv(variance_file, index=False)
        print(f"Variance explained saved to {variance_file}")


def process_genepop_files(data_directory, pattern="*.genepop"):
    """
    Process GENEPOP files in directory, renaming .genepop to .gen.
    
    Args:
        data_directory (str): Directory containing files
        pattern (str): File pattern to match
        
    Returns:
        list: List of processed .gen files
    """
    data_path = Path(data_directory)
    genepop_files = list(data_path.glob(pattern))
    
    processed_files = []
    for genepop_file in genepop_files:
        gen_file = genepop_file.with_suffix('.gen')
        if not gen_file.exists():
            genepop_file.rename(gen_file)
            print(f"Renamed {genepop_file.name} to {gen_file.name}")
        processed_files.append(gen_file)
    
    return processed_files


def analyze_white_crappie(data_directory):
    """
    Perform DAPC analysis on white crappie data.
    
    Args:
        data_directory (str): Directory containing white crappie data
        
    Returns:
        DAPCAnalysis: Analysis object with results
    """
    print("#### White Crappie DAPC Analysis ####")
    
    data_path = Path(data_directory)
    
    # Process genepop files
    process_genepop_files(data_directory)
    
    # Look for populations.haps.gen file
    target_file = data_path / "populations.haps.gen"
    if not target_file.exists():
        print(f"populations.haps.gen not found in {data_directory}")
        return None
    
    # Read genetic data
    reader = GenePopReader(target_file)
    genetic_data = reader.read_genepop()
    
    # Perform DAPC analysis
    dapc = DAPCAnalysis(genetic_data['genetic_matrix'], 
                       genetic_data['populations'],
                       n_pca_components=250)
    
    # Find clusters
    cluster_results = dapc.find_clusters(max_n_clust=5)
    
    # Perform DAPC
    dapc_results = dapc.perform_dapc()
    
    # Create visualizations
    dapc.plot_dapc_results(title="White Crappie DAPC - Clustered by Population",
                          output_file="white_crappie_DAPC_clustered_by_population.png")
    
    dapc.plot_dapc_results(title="White Crappie DAPC - Plotted by Cluster",
                          output_file="white_crappie_DAPC_plotted_by_cluster.png",
                          color_by_population=False)
    
    # Save results
    dapc.save_results("white_crappie")
    
    return dapc


def analyze_freshwater_drum(data_directory):
    """
    Perform DAPC analysis on freshwater drum data.
    
    Args:
        data_directory (str): Directory containing freshwater drum data
        
    Returns:
        DAPCAnalysis: Analysis object with results
    """
    print("#### Freshwater Drum DAPC Analysis ####")
    
    data_path = Path(data_directory)
    
    # Look for freshwater drum file
    target_file = data_path / "freshwater-drum_m3M2n2_p4r7_mac3.haps.gen"
    if not target_file.exists():
        print(f"freshwater-drum_m3M2n2_p4r7_mac3.haps.gen not found in {data_directory}")
        return None
    
    # Read genetic data
    reader = GenePopReader(target_file)
    genetic_data = reader.read_genepop()
    
    # Perform DAPC analysis
    dapc = DAPCAnalysis(genetic_data['genetic_matrix'], 
                       genetic_data['populations'],
                       n_pca_components=250)
    
    # Find clusters
    cluster_results = dapc.find_clusters(max_n_clust=5)
    
    # Perform DAPC
    dapc_results = dapc.perform_dapc()
    
    # Create visualizations
    dapc.plot_dapc_results(title="Freshwater Drum DAPC - Clustered by Population",
                          output_file="freshwater_drum_DAPC_clustered_by_population.png")
    
    dapc.plot_dapc_results(title="Freshwater Drum DAPC - Plotted by Cluster",
                          output_file="freshwater_drum_DAPC_plotted_by_cluster.png",
                          color_by_population=False)
    
    # Save results
    dapc.save_results("freshwater_drum")
    
    return dapc


def main():
    """
    Main function to run DAPC analysis on both species.
    
    Example usage:
        python dapc_al_river_fish.py --crappie_dir /path/to/crappie/data --drum_dir /path/to/drum/data
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='DAPC Analysis for Alabama River Fish')
    parser.add_argument('--crappie_dir', type=str,
                       default='~/Fish/AL-River_migratory-fish/white-crappie/white_crappie_m3M1n1_p4r7_mac3_multi/',
                       help='Directory containing white crappie data')
    parser.add_argument('--drum_dir', type=str,
                       default='~/Fish/AL-River_migratory-fish/freshwater-drum/freshwater-drum_m3M2n2_p4r7_mac3/',
                       help='Directory containing freshwater drum data')
    
    args = parser.parse_args()
    
    # Expand home directories
    crappie_dir = Path(args.crappie_dir).expanduser()
    drum_dir = Path(args.drum_dir).expanduser()
    
    results = {}
    
    # Analyze white crappie
    if crappie_dir.exists():
        results['white_crappie'] = analyze_white_crappie(crappie_dir)
    else:
        print(f"White crappie directory not found: {crappie_dir}")
    
    # Analyze freshwater drum
    if drum_dir.exists():
        results['freshwater_drum'] = analyze_freshwater_drum(drum_dir)
    else:
        print(f"Freshwater drum directory not found: {drum_dir}")
    
    print("\nDAPC Analysis Complete!")
    return results


if __name__ == "__main__":
    main()