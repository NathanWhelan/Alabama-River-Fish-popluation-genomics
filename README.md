# Alabama-River-Fish-popluation-genomics

Scripts used for population genomics in the paper "Evaluating the Effects of Three Alabama River Dams on Movements and Population Connectivity of Freshwater Drum and White Crappie Using Otolith Microchemistry and Genetics"

## Overview

This repository contains both the original R scripts and new Python implementations for analyzing genetic data from Alabama River fish populations. The Python scripts maintain the same statistical rigor while leveraging modern Python libraries for genomics and data analysis.

## Files

### Original R Scripts
- `AMOVA_AL-River-fish.R` - Analysis of Molecular Variance (AMOVA) for population genetics
- `Correct-Ne-for-chromosome-number_AL-River-Fish.R` - Effective population size (Ne) correction for chromosome number
- `DAPC_AL-River-Fish.R` - Discriminant Analysis of Principal Components for population structure

### Python Scripts
- `amova_al_river_fish.py` - Python implementation of AMOVA analysis
- `correct_ne_chromosome_number.py` - Python implementation of Ne correction
- `dapc_al_river_fish.py` - Python implementation of DAPC analysis
- `requirements.txt` - Python dependencies

## Installation

### Python Environment Setup

1. **Clone the repository:**
```bash
git clone https://github.com/NathanWhelan/Alabama-River-Fish-popluation-genomics.git
cd Alabama-River-Fish-popluation-genomics
```

2. **Create a virtual environment (recommended):**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies:**
```bash
pip install -r requirements.txt
```

### Required Python Libraries
- pandas >= 1.5.0 - Data manipulation and analysis
- numpy >= 1.24.0 - Numerical computing
- scipy >= 1.10.0 - Scientific computing and statistics
- scikit-learn >= 1.2.0 - Machine learning algorithms
- scikit-bio >= 0.5.8 - Biological data analysis
- matplotlib >= 3.6.0 - Plotting and visualization
- seaborn >= 0.12.0 - Statistical data visualization

## Usage

### 1. AMOVA Analysis (`amova_al_river_fish.py`)

Performs Analysis of Molecular Variance to assess genetic differentiation within and among populations.

**Basic usage:**
```bash
python amova_al_river_fish.py
```

**With custom directories:**
```bash
python amova_al_river_fish.py --drum_dir /path/to/drum/data --crappie_dir /path/to/crappie/data
```

**Features:**
- Reads GENEPOP format files (.gen)
- Calculates F-statistics and variance components
- Performs randomization tests (1000 permutations)
- Generates publication-ready plots
- Processes both freshwater drum and white crappie data

**Expected input files:**
- `freshwater-drum_m3M2n2_p4r7_mac3.haps.gen`
- Population genetic data in GENEPOP format

**Output:**
- AMOVA results with F-statistics
- Randomization test p-values
- Visualization plots (PNG format)

### 2. Ne Correction (`correct_ne_chromosome_number.py`)

Corrects effective population size estimates for chromosome number using the Waples et al. (2016) method.

**Basic usage:**
```bash
python correct_ne_chromosome_number.py --input NeEstimator_results.csv
```

**With options:**
```bash
python correct_ne_chromosome_number.py \
    --input NeEstimator_white-crappie_all.csv \
    --output NeEst_corrected.csv \
    --chromosomes 24 \
    --plot
```

**Parameters:**
- `--input`: Input CSV file with NeEstimator results
- `--output`: Output CSV file for corrected results
- `--chromosomes`: Number of chromosome pairs (default: 24)
- `--plot`: Generate comparison plots

**Features:**
- Recalculates Ne with higher precision
- Applies chromosome number correction factor
- Adjusts confidence intervals using chi-square distribution
- Generates before/after comparison plots
- Handles edge cases (negative square roots, infinite values)

**Expected input columns:**
- `r2`: Observed r²
- `exp_r2`: Expected r²
- `EffDF`: Effective degrees of freedom
- `Ne_new`: Original Ne estimate

### 3. DAPC Analysis (`dapc_al_river_fish.py`)

Performs Discriminant Analysis of Principal Components to identify genetic clusters and visualize population structure.

**Basic usage:**
```bash
python dapc_al_river_fish.py
```

**With custom directories:**
```bash
python dapc_al_river_fish.py \
    --crappie_dir /path/to/white-crappie/data \
    --drum_dir /path/to/freshwater-drum/data
```

**Features:**
- Automatic optimal cluster detection using silhouette analysis
- PCA for dimensionality reduction
- Linear Discriminant Analysis for classification
- Comprehensive visualization suite
- Exports results in CSV format

**Expected input files:**
- `populations.haps.gen` (for white crappie)
- `freshwater-drum_m3M2n2_p4r7_mac3.haps.gen` (for freshwater drum)

**Output:**
- DAPC coordinate files (CSV)
- Variance explained summaries
- Multiple visualization plots:
  - Main DAPC scatter plots
  - PCA and LDA scree plots
  - Cumulative variance plots

## Data Format Requirements

### GENEPOP Format (.gen files)
The Python scripts expect genetic data in GENEPOP format:
```
Title line
Locus1 Locus2 Locus3 ...
POP
Individual1, 001002 003004 002001 ...
Individual2, 001001 003003 002002 ...
POP
Individual3, 002002 004004 001001 ...
...
```

### NeEstimator Results Format
CSV file with columns:
- `r2`: Observed linkage disequilibrium
- `exp_r2`: Expected r² under null hypothesis
- `EffDF`: Effective degrees of freedom
- `Ne_new`: Estimated effective population size

## Statistical Methods

### AMOVA
- Calculates within and among population variance components
- Computes F-statistics (FST)
- Performs permutation tests for significance testing
- Uses genetic distance matrices

### Ne Correction
- Applies Waples et al. (2016) chromosome correction: CF = 0.09775 + 0.21888×ln(Chr)
- Recalculates confidence intervals using chi-square distribution
- Handles statistical edge cases appropriately

### DAPC
- Combines PCA dimensionality reduction with Linear Discriminant Analysis
- Automated cluster detection using K-means and silhouette analysis
- Maximizes between-group vs within-group genetic variance

## Example Workflow

```bash
# 1. Set up environment
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# 2. Run AMOVA analysis
python amova_al_river_fish.py --drum_dir ./data/drum/ --crappie_dir ./data/crappie/

# 3. Correct Ne estimates
python correct_ne_chromosome_number.py \
    --input NeEstimator_results.csv \
    --chromosomes 24 \
    --plot

# 4. Perform DAPC analysis
python dapc_al_river_fish.py --crappie_dir ./data/crappie/ --drum_dir ./data/drum/
```

## Troubleshooting

### Common Issues

1. **Missing input files**: Ensure GENEPOP files are in the correct directories
2. **File format errors**: Check that genetic data follows GENEPOP format specifications
3. **Memory issues**: For large datasets, consider reducing PCA components
4. **Missing dependencies**: Install all required packages using `pip install -r requirements.txt`

### Performance Notes

- DAPC analysis may take several minutes for large datasets
- Randomization tests (AMOVA) use 1000 permutations by default
- PCA components are automatically selected based on data dimensions

## Citation

If you use these scripts, please cite the original paper:
"Evaluating the Effects of Three Alabama River Dams on Movements and Population Connectivity of Freshwater Drum and White Crappie Using Otolith Microchemistry and Genetics"

## License

MIT License - see LICENSE file for details.

## Contact

For questions about the Python implementations, please open an issue on GitHub.
