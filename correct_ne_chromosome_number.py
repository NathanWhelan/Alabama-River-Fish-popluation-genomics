#!/usr/bin/env python3
"""
Correct Ne for Chromosome Number - Alabama River Fish Population Genomics

This script processes NeEstimator results and corrects effective population size (Ne)
estimates for chromosome number, replicating the functionality of the original R script
Correct-Ne-for-chromosome-number_AL-River-Fish.R.

The correction follows the method from Waples, Larson, Waples (2016) Heredity:
y = 0.09775 + 0.21888*ln(Chr)

Author: Converted from R to Python
Date: 2024
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings


class NeCorrector:
    """
    Class to correct effective population size (Ne) estimates for chromosome number
    and calculate confidence intervals using jackknife resampling.
    """
    
    def __init__(self, data_file, chromosomes=24):
        """
        Initialize Ne corrector.
        
        Args:
            data_file (str): Path to CSV file with NeEstimator results
            chromosomes (int): Number of chromosome pairs (default: 24 for fish species)
        """
        self.data_file = data_file
        self.chromosomes = chromosomes
        self.data = None
        self.correction_factor = None
        
    def load_data(self):
        """
        Load NeEstimator results from CSV file.
        
        Returns:
            pandas.DataFrame: Loaded data
        """
        try:
            self.data = pd.read_csv(self.data_file)
            print(f"Loaded data with {len(self.data)} rows and {len(self.data.columns)} columns")
            print("Columns:", list(self.data.columns))
            return self.data
        except FileNotFoundError:
            print(f"Error: File {self.data_file} not found")
            return None
        except Exception as e:
            print(f"Error loading data: {e}")
            return None
    
    def calculate_correction_factor(self):
        """
        Calculate Ne correction factor using Waples et al. (2016) formula.
        
        Returns:
            float: Correction factor
        """
        # Correction factor from Waples, Larson, Waples. 2016 Heredity
        # y = 0.09775 + 0.21888*ln(Chr)
        # where Chr = # of chromosome pairs
        
        self.correction_factor = 0.09775 + (0.21888 * np.log(self.chromosomes))
        print(f"Correction factor for {self.chromosomes} chromosome pairs: {self.correction_factor:.6f}")
        return self.correction_factor
    
    def step1_recalculate_ne(self):
        """
        STEP 1: Recalculate r^2' and Ne with higher precision.
        
        This step recalculates values that NeEstimator may not output with
        sufficient decimal places.
        """
        if self.data is None:
            print("Error: Data not loaded. Call load_data() first.")
            return
        
        print("STEP 1: Recalculating r^2' and Ne...")
        
        # Calculate r^2' (r2p)
        self.data['r2p'] = self.data['r2'] - self.data['exp_r2']
        
        # Recalculate Ne with higher precision
        # Handle cases where square root argument is negative
        sqrt_arg = 0.308**2 - 2.08 * self.data['r2p']
        
        # Set negative square root arguments to zero
        sqrt_term = np.where(sqrt_arg > 0, np.sqrt(sqrt_arg), 0)
        
        self.data['Ne_new'] = (0.308 + sqrt_term) / (2 * self.data['r2p'])
        
        # Handle division by zero or very small values
        self.data['Ne_new'] = self.data['Ne_new'].replace([np.inf, -np.inf], np.nan)
        
        print(f"Recalculated Ne for {len(self.data)} populations")
        
    def step2_apply_chromosome_correction(self):
        """
        STEP 2: Apply chromosome number correction to Ne estimates.
        """
        if self.correction_factor is None:
            self.calculate_correction_factor()
        
        print("STEP 2: Applying chromosome number correction...")
        
        # Correct Ne by dividing by correction factor
        self.data['NeCorrected'] = self.data['Ne_new'] / self.correction_factor
        
        print(f"Applied correction factor {self.correction_factor:.6f} to Ne estimates")
        
    def step3_adjust_confidence_intervals(self):
        """
        STEP 3: Adjust jackknife confidence intervals using effective degrees of freedom.
        """
        print("STEP 3: Adjusting confidence intervals...")
        
        # (1) Calculate chi-square quantiles at 0.025 and 0.975
        self.data['df025'] = self.data['EffDF'].apply(
            lambda x: stats.chi2.ppf(0.025, df=x, loc=0, scale=1) if pd.notna(x) and x > 0 else np.nan
        )
        
        self.data['df975'] = self.data['EffDF'].apply(
            lambda x: stats.chi2.ppf(0.975, df=x, loc=0, scale=1) if pd.notna(x) and x > 0 else np.nan
        )
        
        # (2) Calculate adjusted r2'
        self.data['adjR2p'] = self.data['NeCorrected'].apply(
            lambda x: (0.308/x) - (0.52/(x**2)) if pd.notna(x) and x != 0 else np.nan
        )
        
        # (3) Calculate adjusted overall r2
        self.data['adjobR2'] = self.data['adjR2p'] + self.data['exp_r2']
        
        # (4) Recalculate CI r2 values
        self.data['r2025'] = (self.data['EffDF'] * self.data['adjobR2']) / self.data['df025']
        self.data['r2975'] = (self.data['EffDF'] * self.data['adjobR2']) / self.data['df975']
        
        # (5) Calculate r2' for confidence intervals
        self.data['r2p025'] = self.data['r2025'] - self.data['exp_r2']
        self.data['r2p975'] = self.data['r2975'] - self.data['exp_r2']
        
        # (6) Calculate corrected confidence intervals
        # Low CI (using r2p975 for conservative estimate)
        sqrt_arg_low = 0.308**2 - 2.08 * self.data['r2p975']
        sqrt_term_low = np.where(sqrt_arg_low > 0, np.sqrt(sqrt_arg_low), 0)
        self.data['lowCIcorrected'] = (0.308 + sqrt_term_low) / (2 * self.data['r2p975'])
        
        # High CI (using r2p025)
        sqrt_arg_high = 0.308**2 - 2.08 * self.data['r2p025']
        sqrt_term_high = np.sqrt(sqrt_arg_high)  # Should always be positive for high CI
        self.data['highCIcorrected'] = (0.308 + sqrt_term_high) / (2 * self.data['r2p025'])
        
        # Handle infinite values (negative high CIs represent "infinity" estimates)
        self.data['lowCIcorrected'] = self.data['lowCIcorrected'].replace([np.inf, -np.inf], np.nan)
        self.data['highCIcorrected'] = self.data['highCIcorrected'].replace([np.inf, -np.inf], np.nan)
        
        print("Confidence intervals adjusted")
        
    def process_all_steps(self):
        """
        Execute all three steps of Ne correction process.
        
        Returns:
            pandas.DataFrame: Processed data with corrected Ne values
        """
        if self.data is None:
            self.load_data()
            
        if self.data is None:
            return None
        
        # Execute all steps
        self.step1_recalculate_ne()
        self.step2_apply_chromosome_correction()
        self.step3_adjust_confidence_intervals()
        
        return self.data
    
    def save_results(self, output_file):
        """
        Save corrected results to CSV file.
        
        Args:
            output_file (str): Path for output CSV file
        """
        if self.data is None:
            print("Error: No data to save. Run process_all_steps() first.")
            return
        
        try:
            self.data.to_csv(output_file, index=False)
            print(f"Results saved to {output_file}")
        except Exception as e:
            print(f"Error saving results: {e}")
    
    def print_summary(self):
        """
        Print summary of correction results.
        """
        if self.data is None:
            print("Error: No data to summarize.")
            return
        
        print("\n" + "="*60)
        print("Ne CORRECTION SUMMARY")
        print("="*60)
        print(f"Number of populations: {len(self.data)}")
        print(f"Chromosome pairs used: {self.chromosomes}")
        print(f"Correction factor: {self.correction_factor:.6f}")
        
        if 'Ne_new' in self.data.columns and 'NeCorrected' in self.data.columns:
            print(f"\nOriginal Ne range: {self.data['Ne_new'].min():.2f} - {self.data['Ne_new'].max():.2f}")
            print(f"Corrected Ne range: {self.data['NeCorrected'].min():.2f} - {self.data['NeCorrected'].max():.2f}")
            print(f"Mean correction: {(self.data['Ne_new'] / self.data['NeCorrected']).mean():.4f}")
        
        print("\nKey columns in output:")
        key_columns = ['Ne_new', 'NeCorrected', 'lowCIcorrected', 'highCIcorrected']
        for col in key_columns:
            if col in self.data.columns:
                print(f"  {col}: {self.data[col].notna().sum()} valid values")
        
        print("="*60)
    
    def plot_ne_comparison(self, output_file=None):
        """
        Create plots comparing original and corrected Ne values.
        
        Args:
            output_file (str): Optional path to save the plot
        """
        if self.data is None or 'Ne_new' not in self.data.columns:
            print("Error: Data not available for plotting.")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Ne Correction Analysis', fontsize=16)
        
        # 1. Before vs After correction
        axes[0, 0].scatter(self.data['Ne_new'], self.data['NeCorrected'], alpha=0.7)
        axes[0, 0].plot([self.data['Ne_new'].min(), self.data['Ne_new'].max()], 
                       [self.data['Ne_new'].min(), self.data['Ne_new'].max()], 
                       'r--', label='No correction')
        axes[0, 0].set_xlabel('Original Ne')
        axes[0, 0].set_ylabel('Corrected Ne')
        axes[0, 0].set_title('Original vs Corrected Ne')
        axes[0, 0].legend()
        
        # 2. Correction factor effect
        correction_ratio = self.data['Ne_new'] / self.data['NeCorrected']
        axes[0, 1].hist(correction_ratio, bins=20, alpha=0.7, edgecolor='black')
        axes[0, 1].axvline(correction_ratio.mean(), color='red', linestyle='--', 
                          label=f'Mean: {correction_ratio.mean():.4f}')
        axes[0, 1].set_xlabel('Correction Ratio (Original/Corrected)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Distribution of Correction Ratios')
        axes[0, 1].legend()
        
        # 3. Confidence intervals
        if 'lowCIcorrected' in self.data.columns and 'highCIcorrected' in self.data.columns:
            valid_data = self.data.dropna(subset=['NeCorrected', 'lowCIcorrected', 'highCIcorrected'])
            if len(valid_data) > 0:
                x_pos = range(len(valid_data))
                # Calculate error bars, ensuring no negative values
                lower_err = np.maximum(0, valid_data['NeCorrected'] - valid_data['lowCIcorrected'])
                upper_err = np.maximum(0, valid_data['highCIcorrected'] - valid_data['NeCorrected'])
                
                axes[1, 0].errorbar(x_pos, valid_data['NeCorrected'],
                                   yerr=[lower_err, upper_err],
                                   fmt='o', alpha=0.7, capsize=3)
                axes[1, 0].set_xlabel('Population Index')
                axes[1, 0].set_ylabel('Corrected Ne')
                axes[1, 0].set_title('Ne with Confidence Intervals')
        
        # 4. Log-scale comparison if values span large range
        if self.data['NeCorrected'].max() / self.data['NeCorrected'].min() > 10:
            axes[1, 1].scatter(np.log10(self.data['Ne_new']), 
                              np.log10(self.data['NeCorrected']), alpha=0.7)
            axes[1, 1].plot([np.log10(self.data['Ne_new']).min(), np.log10(self.data['Ne_new']).max()], 
                           [np.log10(self.data['Ne_new']).min(), np.log10(self.data['Ne_new']).max()], 
                           'r--', label='No correction')
            axes[1, 1].set_xlabel('Log10(Original Ne)')
            axes[1, 1].set_ylabel('Log10(Corrected Ne)')
            axes[1, 1].set_title('Log-scale Ne Comparison')
            axes[1, 1].legend()
        else:
            axes[1, 1].axis('off')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {output_file}")
        
        plt.show()


def main():
    """
    Main function to process NeEstimator results and correct for chromosome number.
    
    Example usage:
        python correct_ne_chromosome_number.py --input NeEstimator_white-crappie_all.csv --chromosomes 24
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='Correct Ne estimates for chromosome number')
    parser.add_argument('--input', type=str, 
                       default='NeEstimator_white-crappie_all.csv',
                       help='Input CSV file with NeEstimator results')
    parser.add_argument('--output', type=str,
                       default='NeEst_corrected.csv',
                       help='Output CSV file for corrected results')
    parser.add_argument('--chromosomes', type=int, default=24,
                       help='Number of chromosome pairs (default: 24)')
    parser.add_argument('--plot', action='store_true',
                       help='Generate comparison plots')
    
    args = parser.parse_args()
    
    # Create Ne corrector
    corrector = NeCorrector(args.input, chromosomes=args.chromosomes)
    
    # Process all correction steps
    results = corrector.process_all_steps()
    
    if results is not None:
        # Print summary
        corrector.print_summary()
        
        # Save results
        corrector.save_results(args.output)
        
        # Generate plots if requested
        if args.plot:
            plot_file = args.output.replace('.csv', '_comparison.png')
            corrector.plot_ne_comparison(plot_file)
        
        print(f"\nProcessing complete. Results saved to {args.output}")
        
        # Display first few rows of results
        print("\nFirst 5 rows of corrected data:")
        display_cols = ['Ne_new', 'NeCorrected', 'lowCIcorrected', 'highCIcorrected']
        available_cols = [col for col in display_cols if col in results.columns]
        if available_cols:
            print(results[available_cols].head())
        
    else:
        print("Error: Could not process data")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())