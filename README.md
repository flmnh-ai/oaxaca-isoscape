# Oaxaca Isoscape Project

This repository contains the analysis code for the Oaxaca Strontium Isoscape project, which models and maps strontium isotope ratios (87Sr/86Sr) in Mexico with a focus on the Oaxaca region. The project uses Bayesian Additive Regression Trees (BART) to create high-resolution isoscapes for archaeological origin assignment studies.

## Project Structure

```
oaxaca-isoscape/
├── R/
│   ├── 01_preprocess_data.R      # Data preprocessing (run once)
│   ├── 02_fit_model.R            # BART model fitting (run once)
│   ├── helpers.R                 # Utility functions
│   ├── aerosol_deposition.R      # MERRA-2 data preprocessing
│   └── archive/                  # Experimental code (reference only)
├── notebooks/
│   ├── analysis.qmd              # Analysis and figure generation
│   └── archive/                  # Old notebooks (reference only)
├── data/
│   ├── raw/                      # Original data files
│   └── derived/                  # Processed data files
├── outputs/
│   ├── figures/                  # Publication figures
│   └── models/                   # Saved BART model
├── DATA.md                       # Input data documentation
├── CLAUDE.md                     # Developer guidance
└── README.md                     # This file
```

## Running the Analysis

The analysis is organized as **two R scripts** for data generation and **one Quarto notebook** for visualization.

### Step 1: Data Preprocessing

Process global isotope data and environmental predictors:

```bash
Rscript R/01_preprocess_data.R
```

**Note:** Run `Rscript R/aerosol_deposition.R` first if aerosol data hasn't been preprocessed.

**Creates:**
- `data/derived/mexico_predictors.rds` (predictor rasters for Mexico)
- `data/derived/pts_combined.rds` (isotope observations)
- `data/derived/oaxaca_plants.rds` (Oaxaca plant samples)
- `data/derived/dat.rds` (full modeling dataset)

**Runtime:** ~30-60 minutes

### Step 2: Model Fitting

Fit BART model with hyperparameter tuning:

```bash
Rscript R/02_fit_model.R
```

**Creates:**
- `outputs/models/bart_final.rds` (final BART model)
- `data/derived/train.rds` (training dataset)
- `data/derived/test.rds` (test dataset)

**Runtime:** Several hours (uses Bayesian optimization + 8 MCMC chains)

### Step 3: Visualization

Generate all publication figures:

```bash
quarto render notebooks/analysis.qmd
```

**Creates:**
- All figures in `outputs/figures/`
- HTML report with analysis

**Runtime:** ~10-20 minutes (excluding prediction generation, which is cached)

## Required Packages

### Core packages
- **tidyverse** - data manipulation
- **tidymodels** - modeling framework
- **dbarts** - Bayesian Additive Regression Trees

### Spatial data
- **sf**, **terra**, **stars** - spatial data structures
- **rnaturalearth** - base maps
- **dggridR** - discrete global grids

### Visualization
- **patchwork** - plot composition
- **scico** - scientific color palettes
- **ggrepel** - label positioning
- **ggnewscale** - multiple fill scales
- **RColorBrewer** - color palettes

### Utilities
- **here** - path management
- **readxl** - Excel file reading
- **bundle**, **butcher** - model serialization
- **coda** - MCMC diagnostics

### Optional (for satellite imagery)
- **rgee** - Google Earth Engine interface
- **reticulate** - Python integration

## Data Sources

This analysis integrates multiple global and regional datasets:

### Isotope Data
- Global strontium isotope databases (Bataille et al. 2024; Wang et al. 2021; Scaffidi & Knudson 2020)
- CAMBIO human isotope database
- Oaxaca plant sampling campaign
- Monte Alban archaeological individuals

### Environmental Predictors
- **Bedrock geology** - Bataille et al. (2020) global isoscape predictors
- **Soil properties** - SoilGrids 1km aggregated data
- **Climate** - CHELSA V2.1 bioclimatic variables (1981-2010)
- **Elevation** - SRTM 90m data
- **Aerosol deposition** - MERRA-2 reanalysis

**For complete data documentation with citations and URLs, see [DATA.md](DATA.md)**

## Main Outputs

1. **BART prediction model** for Mexico strontium isotope ratios
2. **High-resolution isoscape** (1km resolution) with uncertainty estimates
3. **Publication figures:**
   - Study area geological/satellite maps
   - Training data distributions
   - Isoscape predictions for Mexico and Oaxaca
   - Monte Alban origin assignment analysis
4. **Model diagnostics:** convergence checks, variable importance, performance metrics

## Reproducibility

To reproduce the entire analysis from scratch:

```bash
# 1. Preprocess aerosol data (if not already done)
Rscript R/aerosol_deposition.R

# 2. Preprocess all data (~30-60 min)
Rscript R/01_preprocess_data.R

# 3. Fit BART model (several hours)
Rscript R/02_fit_model.R

# 4. Generate figures (~10-20 min)
quarto render notebooks/analysis.qmd
```

**System Requirements:**
- R version 4.0 or higher
- At least 16GB RAM (for global data processing)
- Multi-core CPU recommended (model fitting uses 8 threads)

## Documentation

- **[DATA.md](DATA.md)** - Complete documentation of all input data sources with citations
- **[CLAUDE.md](CLAUDE.md)** - Developer guidance for working with this codebase
- **[notebooks/archive/](notebooks/archive/)** - Previous notebook-based workflow (reference only)

## Key Technical Details

- **Data transformation:** Strontium isotope ratios are logit-transformed before modeling to handle bounded [0.703, 0.78] range
- **Spatial projections:** Eckert IV for global analysis, Albers Equal Area for Mexico predictions
- **Dimensionality reduction:** Principal Component Analysis reduces environmental predictors (bedrock, soil, climate, aerosols) to key components
- **Cross-validation:** Spatial holdout - train on Americas, test on Mexico/Oaxaca
- **Model specification:** 400 trees, 8 MCMC chains, 500 post-burnin samples with thinning

## Citation

If you use this code or methodology, please cite:

> Gauthier, N. (2025). Oaxaca Strontium Isoscape. [Paper in review]

## Author

**Nick Gauthier**
University of Florida

## License

See [LICENSE](LICENSE) file for details.
