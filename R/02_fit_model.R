# ==============================================================================
# Oaxaca Isoscape Model Fitting
# ==============================================================================
# This script performs BART model hyperparameter tuning and training
# Generates:
#   - outputs/models/bart_final.rds (Final BART model for predictions)
#   - data/derived/train.rds (Training dataset)
#   - data/derived/test.rds (Test dataset)
#
# This is a computationally expensive script that should only be run once
# Hyperparameter tuning uses Bayesian optimization with spatial cross-validation
# Final model training uses 8 MCMC chains with 500 post-burnin samples
# ==============================================================================

# Load packages ----------------------------------------------------------------
library(tidymodels)
library(here)
library(sf)
library(coda)
library(butcher)
library(bundle)
library(future)
source(here('R/helpers.R'))

# Define spatial regions -------------------------------------------------------
mexico_bbox <- st_bbox(c(xmin = -106, xmax = -87, ymin = 13, ymax = 23),
                       crs = 4326) |>
  st_transform(eck4)

train_bbox <- st_bbox(c(xmin = -115, xmax = -50, ymin = -90, ymax = 90),
                       crs = 4326) |>
  st_transform(eck4)

# Load and prepare data --------------------------------------------------------
message("Loading and preparing data...")

dat <- readRDS(here('data/derived/dat.rds')) |>
  filter(!if_any(c(bedrock_age_sr_gradient:distance), is.na)) |>
  mutate(Sr = logit_transform(Sr)) |>
  mutate(mexico = lengths(st_intersects(geometry, st_as_sfc(mexico_bbox))) > 0) |>
  st_crop(train_bbox) |>
  select(-sediment) # Sediment variable not used in final model

# Create train/test splits
train <- dat |>
  filter(!mexico) |>
  select(-mexico)

test <- dat |>
  filter(mexico) |>
  select(-mexico) |>
  mutate(type2 = replace_na(type2, ''),
         material = replace_na(material, ''))

# Save splits for use in analysis notebook
saveRDS(train, here('data/derived/train.rds'))
saveRDS(test, here('data/derived/test.rds'))

message("Train: ", nrow(train), " observations")
message("Test: ", nrow(test), " observations")

# Create data splits for cross-validation -------------------------------------
# Strategy: Train on North + South America, test on Mesoamerica
splits <- make_splits(train, test)
splits_oaxaca <- make_splits(bind_rows(train, filter(test, !oaxaca)), filter(test, oaxaca))
splits_mexico <- make_splits(train, filter(test, !oaxaca))

# Create 5-fold spatial cross-validation
mexico_folds <- manual_rset(list(splits_mexico,
                                 splits_mexico,
                                 splits_mexico,
                                 splits_mexico,
                                 splits_mexico), ids = as.character(1:5))

# ==============================================================================
# MODEL SPECIFICATION
# ==============================================================================

# Define recipe ----------------------------------------------------------------
# Logit-transformed Sr ~ environmental predictors
# Normalize predictors for better BART convergence
rec <- recipe(Sr ~ ., train) |>
  update_role(c(geometry, oaxaca, type, type2, material, seqnum),
              new_role = 'other') |>
  update_role_requirements('other', bake = FALSE) |>
  step_impute_knn(starts_with('soil')) |> # Impute coastal/urban NAs
  step_normalize(all_numeric_predictors())

# Define BART model with tunable hyperparameters ------------------------------
bart_mod <- parsnip::bart(mode = 'regression',
                 trees = tune(),
                 prior_terminal_node_coef = tune(),
                 prior_terminal_node_expo = tune(),
                 prior_outcome_range = tune()
              ) %>%
  set_engine('dbarts', nskip = 800)

bart_wflw <- workflow(rec, bart_mod)

# ==============================================================================
# HYPERPARAMETER TUNING
# ==============================================================================

# Initial grid search ----------------------------------------------------------
message("\n========================================")
message("Starting hyperparameter grid search...")
message("========================================\n")

bart_params <- crossing(
  trees = c(50, 100, 200, 300, 400),
  prior_outcome_range = c(8, 12, 16, 18, 20, 24, 32),
  prior_terminal_node_coef = c(0.85, 0.9, .95, .99),
  prior_terminal_node_expo = c(1, 1.5, 2)
)

set.seed(1111)
plan(multisession)

bart_grid <- tune_grid(
   bart_wflw,
   resamples = mexico_folds,
   grid = bart_params,
   control = control_grid(pkgs = 'sf')
 )

plan(sequential)

message("\nBest grid results:")
print(show_best(bart_grid, n = 5))

# Bayesian optimization --------------------------------------------------------
message("\n========================================")
message("Starting Bayesian optimization...")
message("========================================\n")

control <- control_bayes(
  pkgs = 'sf',
  no_improve = 30,
  uncertain = 5,
  verbose_iter = TRUE,
  parallel_over = 'everything'
)

set.seed(1111)
plan(multisession)

bart_bayes <- tune_bayes(
   bart_wflw,
   resamples = mexico_folds,
   iter = 60,
   initial = bart_grid,
   control = control
 )

plan(sequential)

message("\nBest Bayesian optimization results:")
print(show_best(bart_bayes, n = 5))
print(select_by_one_std_err(bart_bayes, trees, 'rmse'))

# ==============================================================================
# FINAL MODEL TRAINING
# ==============================================================================

# Define final workflow with optimized hyperparameters ------------------------
# Based on hyperparameter tuning results:
# trees = 400, prior_outcome_range = 15,
# prior_terminal_node_expo = 1.5, prior_terminal_node_coef = 0.99
message("\n========================================")
message("Training final model...")
message("========================================\n")

bart_wflw_final <- parsnip::bart(mode = 'regression',
                        trees = 400,
                        prior_outcome_range = 15,
                        prior_terminal_node_expo = 1.5,
                        prior_terminal_node_coef = 0.99
) %>%
  set_engine('dbarts',
             nskip = 800,    # Burn-in samples
             ndpost = 500,   # Post-burn-in samples
             keepevery = 4,  # Thinning interval
             nchain = 8,     # Number of MCMC chains
             nthread = 8,    # Parallel threads
             combinechains = FALSE
  ) %>%
  workflow(rec, .)

# Evaluation models (for validation metrics) ----------------------------------
message("Fitting evaluation models for cross-validation...")

set.seed(42)
results_mexico <- last_fit(bart_wflw_final, splits_mexico)
results_oaxaca <- last_fit(bart_wflw_final, splits_oaxaca)

# Check convergence
message("\nChecking convergence for Mexico model...")
convergence(results_mexico)

message("\nChecking convergence for Oaxaca model...")
convergence(results_oaxaca)

# Print evaluation metrics
message("\n=== Model Performance ===")

message("\nTrain on Americas, test on Mesoamerica (logit scale):")
set.seed(42)
metrics_mex_logit <- extract_workflow(results_mexico) |>
  augment(filter(test, !oaxaca)) |>
  summarize(rmse = rmse_vec(.pred, Sr),
           rsq = rsq_vec(.pred, Sr))
print(metrics_mex_logit)

message("\nTrain on Americas, test on Mesoamerica (original scale):")
set.seed(42)
metrics_mex_orig <- extract_workflow(results_mexico) |>
  augment(filter(test, !oaxaca)) |>
  mutate(across(.pred:Sr, ~inverse_logit_transform(.x))) |>
  summarize(rmse = rmse_vec(.pred, Sr),
           rsq = rsq_vec(.pred, Sr))
print(metrics_mex_orig)

message("\nTrain on Americas + Mesoamerica, test on Oaxaca (logit scale):")
set.seed(42)
metrics_oax_logit <- extract_workflow(results_oaxaca) |>
  augment(filter(test, oaxaca)) |>
  summarize(rmse = rmse_vec(.pred, Sr),
           rsq = rsq_vec(.pred, Sr))
print(metrics_oax_logit)

message("\nTrain on Americas + Mesoamerica, test on Oaxaca (original scale):")
set.seed(42)
metrics_oax_orig <- extract_workflow(results_oaxaca) |>
  augment(filter(test, oaxaca)) |>
  mutate(across(.pred:Sr, ~inverse_logit_transform(.x))) |>
  summarize(rmse = rmse_vec(.pred, Sr),
           rsq = rsq_vec(.pred, Sr))
print(metrics_oax_orig)

# Fit final model on all data --------------------------------------------------
message("\n========================================")
message("Fitting final model on full dataset...")
message("========================================\n")

set.seed(42)
bart_final <- fit(bart_wflw_final, dat)

# Check convergence
message("\nChecking convergence for final model...")
convergence(bart_final)

# Variable importance
message("\nCalculating variable importance...")
vars <- varimp(bart_final)
print(vars)

# Save final model -------------------------------------------------------------
message("\nSaving final model...")

bart_final |>
  butcher() |>
  bundle() |>
  saveRDS(here('outputs/models/bart_final.rds'))

message("\n========================================")
message("Model fitting complete!")
message("========================================")
message("\nCreated files:")
message("  - outputs/models/bart_final.rds")
message("  - data/derived/train.rds")
message("  - data/derived/test.rds")
message("\nUse the saved model in the analysis notebook for predictions.")
