# SILAC Peptide Decay Modeling

This repository contains R scripts and supporting files used for analyzing proteome turnover using SILAC-labeled mass spectrometry data. The core analysis estimates peptide degradation rates over time using nonlinear curve fitting, enabling comparisons across biological conditions (e.g., acute vs chronic T-cell exhaustion).

> 📝 **This analysis was used in the manuscript:  
> _Comprehensive analysis of proteome turnover dynamics during T cell exhaustion._**

## Overview

The workflow includes:

1. **Data Preprocessing**  
   - Cleaning peptide-level quantification data  
   - Removing post-translational modifications (e.g., +57 for carbamidomethylation)  
   - Restructuring isotope-labeled intensities

2. **Model Fitting**  
   - Nonlinear least squares fitting of exponential decay models per peptide  
   - Filtering peptides by model fit (adjusted R² ≥ 0.9)

3. **Summary & Export**  
   - Calculating half-lives per protein condition/replicate  
   - Generating visualizations of percent heavy/light across time  
   - Exporting cleaned and summarized results

## Directory Structure

```
├── scripts/
│   └── silac_decay_modeling.R        # Full decay modeling pipeline
│   └── preprocess_ratios.R           # Peptide quantification and cleaning
├── output/                           # Processed results (not committed)
├── data/                             # Raw inputs (not committed)
├── README.md
└── .gitignore
```

## Requirements

Install dependencies using:

```r
install.packages(c("tidyverse", "broom", "cowplot", "qs", "furrr", "modelr"))
```

Also install `nlfitr` and `nplyr` if they are custom or from GitHub:
```r
# Example (if hosted on GitHub)
# devtools::install_github("yourusername/nlfitr")
```

## Input

- Input `.csv` must include columns like `Peptide`, `Protein`, `Area`, `Isotope_Label_Type`, `Replicate`, etc.
- Quantitative filtering is done using a `Quantitative` logical column.

## Output

- Ratio plots (boxplots) of heavy/light isotope per replicate
- Half-life estimates per protein condition
- Cleaned `.tsv` tables for downstream analysis

## Citation

If you use this code, please cite:

**Comprehensive analysis of proteome turnover dynamics during T cell exhaustion**  
(Authors list, Journal, Year, DOI or link if available)

## License

MIT License (or specify your lab's preferred license)

## Contact

Maintained by: Dr. Michael Bauer  
Department of Biomedical Informatics  
University of Arkansas for Medical Sciences
