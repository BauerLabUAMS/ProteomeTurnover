# SILAC Peptide Decay Modeling

This repository contains R scripts and supporting files used for analyzing SILAC-labeled proteomics data to estimate peptide degradation rates using nonlinear modeling. The analysis workflow is designed for paired acute and chronic T-cell conditions and includes model fitting, filtering, and summary statistics for protein turnover.

## Overview

The main script performs the following steps:

1. **Preprocessing**: Calculates percent light labeling from heavy ratios and filters for valid timepoints.
2. **Model Fitting**: Uses nonlinear least squares (NLS) to fit an exponential decay model to each peptide.
3. **Model Evaluation**: Calculates pseudo and adjusted R² to assess model quality.
4. **Averaging and Refitting**: Averages signal across replicates, re-fits the decay model, and recalculates half-lives.
5. **Export**: Outputs half-life estimates for downstream analysis and visualization.

## Directory Structure
