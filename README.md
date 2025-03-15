# Oaxaca Isoscape Project

This repository contains the analysis code for the Oaxaca Strontium Isoscape project, which models and maps strontium isotope ratios in Mexico with a focus on the Oaxaca region.

## Project Structure

```
oaxaca-isoscape/
├── data/
│   ├── raw/              # Original data files
│   └── derived/          # Processed data files
├── R/
│   └── helpers.R         # Helper functions used across analyses
├── analysis/
│   ├── 01-data-prep.qmd  # Data preparation and exploration
│   ├── 02-modeling.qmd   # Model building and evaluation
│   ├── 03-visualization.qmd  # Final visualizations and maps
├── outputs/
│   ├── figures/          # Exported figures
│   └── models/           # Saved model objects
└── README.md             # This file
```

## Running the Analysis

The analysis is structured as a sequence of Quarto documents that should be run in order:

1. First, run `analysis/01-data-prep.qmd` to process all data sources
2. Then, run `analysis/02-modeling.qmd` to build and evaluate models
3. Finally, run `analysis/03-visualization.qmd` to create maps and visualizations

## Required Packages

This analysis requires the following R packages:

- tidyverse (data manipulation)
- sf, terra, stars, rnaturalearth (spatial data)
- tidymodels, dbarts (modeling)
- patchwork (plot composition)

## Data Sources

- Global isotope data from Bataille et al. (2020) and others
- Oaxaca isotope samples (OAX Sr Result Table)
- CHELSA climate data
- Global soil data from SoilGrids
- Elevation data

## Output

The main outputs of this analysis are:

1. A strontium isoscape prediction model for Mexico
2. High-resolution maps of predicted strontium isotope ratios
3. Uncertainty estimations for predictions
4. Analysis of model performance and variable importance

## Author

Nick Gauthier
