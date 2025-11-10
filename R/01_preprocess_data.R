# ==============================================================================
# Oaxaca Isoscape Data Preprocessing
# ==============================================================================
# This script processes global isotope data and environmental predictors
# to create the modeling dataset. Run this script once to generate:
#   - data/derived/mexico_predictors.rds (Mexico-cropped predictor rasters)
#   - data/derived/pts_combined.rds (All isotope observations)
#   - data/derived/dat.rds (Full modeling dataset with predictors)
#   - data/derived/oaxaca_plants.rds (Oaxaca plant samples)
#
# Prerequisites: Run R/aerosol_deposition.R first to create aerosols.rds
# ==============================================================================

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(readxl)
library(stars)
library(terra)
library(sf)
library(rnaturalearth)
library(patchwork)
library(here)
library(dggridR)
source(here('R/helpers.R'))

# ==============================================================================
# GLOBAL RASTER PREDICTORS
# ==============================================================================

# Bedrock isotope predictors ---------------------------------------------------
# Get predictor files from Bataille et al. (2020)
pred_files <- here('data/raw/isotope_predictors') |>
  list.files(pattern = '\\.tif$', full.names = TRUE)

pred_names <- pred_files |>
  str_split_i('/', 11) |>
  str_remove('\\.tif$') |>
  str_remove('_reproj') |>
  str_remove('r\\.')

sediment_tibble <- tibble(
  id = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
  category = c('none', 'colluvial', 'eolian', 'glacial', 'lacustrine',
               'marine', 'organic', 'evaporite', 'pyroclastics', 'coastal', 'ice')
)

# Process variables (different origins and extents require step-by-step processing)
vars1 <- rast(pred_files[c(4:7)])
names(vars1) <- pred_names[c(4:7)]
levels(vars1[[2]]) <- sediment_tibble
names(vars1[[2]]) <- 'sediment'
names(vars1[[4]]) <- 'sr_median'

vars2 <- rast(pred_files[c(1:3, 8:9)])
crs(vars2) <- crs(vars1)
names(vars2) <- pred_names[c(1:3, 8:9)]

predictors <- c(vars1, vars2)
bedrock <- predictors[[4:9]]
other_predictors <- predictors[[1:3]]

# Calculate principal components of bedrock predictors
set.seed(42)
bedrock_pcs <- prcomp(bedrock, scale. = TRUE, maxcell = 10000000, retx = FALSE)
bedrock_pc_names <- c("bedrock_age_sr_gradient", "bedrock_young_highSrVar")

# Soil predictors --------------------------------------------------------------
# SoilGrids 250m aggregated to 1km
soils <- here('data/raw/soilgrids') |>
  list.files(full.names = TRUE) |>
  rast()

names(soils) <- names(soils) |>
  str_replace('-', '_')

soil_pcs <- prcomp(soils, scale. = TRUE, maxcell = 10000000, retx = FALSE)
soil_pc_names <- c(
  "soil_moisture_and_bulk_density",
  "soil_clay_cec_gradient",
  "soil_clay_cec_texture"
)

# Elevation --------------------------------------------------------------------
# SRTM elevation data from geodata package
elevation <- geodata::elevation_global(0.5, path = here('data/raw'))
names(elevation) <- 'elevation'

# Climate ----------------------------------------------------------------------
# CHELSA V2.1 climate data
climate <- list.files(here('data/raw/CHELSA_V2_bio_clim'), full.names = TRUE) |>
  rast()

names(climate) <- names(climate) |>
  str_remove('CHELSA_') |>
  str_remove('_1981-2010_V.2.1') |>
  str_remove('1981-2010_') |>
  str_remove('_V.2.1')

climate_pcs <- prcomp(climate, scale. = TRUE, maxcell = 10000000, retx = FALSE)
climate_pc_names <- c(
  "climate_temp_dryness_gradient",
  "climate_moisture_precipitation",
  "climate_temperature_extremes",
  "climate_seasonal_precipitation",
  "climate_precipitation_intensity_timing",
  "climate_radiation_vs_wind_dryness",
  "climate_subtle_moisture_regime"
)

# Aerosol Deposition -----------------------------------------------------------
# Dry and wet aerosol deposition from MERRA-2
# Note: Must run R/aerosol_deposition.R first
aerosols_stars <- readRDS(here('data/derived/aerosols.rds'))

aerosols <- rast(aerosols_stars)
names(aerosols) <- names(aerosols_stars)

aerosol_pcs <- prcomp(aerosols, scale. = TRUE, maxcell = 10000000, retx = FALSE)
aerosol_pc_names <- c(
  "aerosol_carbonaceous_sulfate_abundance",
  "aerosol_sea_salt_signature",
  "aerosol_dust_concentration",
  "aerosol_dry_wet_sulfate_organic_contrast",
  "aerosol_dry_carbon_dominance"
)

# ==============================================================================
# MEXICO PREDICTOR DATA
# ==============================================================================
# Crop to Mexico bounding box to reduce memory footprint for spatial predictions

# Function to progressively fill gaps in rasters
fillNA <- function(x, method = 'mean') {
  focal(x, w = 3, fun = method, na.policy = 'only', na.rm = TRUE) |>
    focal(w = 3, fun = method, na.policy = 'only', na.rm = TRUE) |>
    focal(w = 3, fun = method, na.policy = 'only', na.rm = TRUE) |>
    focal(w = 3, fun = method, na.policy = 'only', na.rm = TRUE) |>
    focal(w = 3, fun = method, na.policy = 'only', na.rm = TRUE) |>
    focal(w = 3, fun = method, na.policy = 'only', na.rm = TRUE) |>
    focal(w = 3, fun = method, na.policy = 'only', na.rm = TRUE) |>
    focal(w = 3, fun = method, na.policy = 'only', na.rm = TRUE)
}

# Create Mexico predictor rasters
albers <- '+proj=aea +lon_0=-96 +lat_1=16.5 +lat_2=21.5 +lat_0=19 +datum=WGS84 +units=m +no_defs'
mexico_bbox_eck <- st_bbox(c(xmin = -106, xmax = -87, ymin = 15, ymax = 23), crs = 4326) |>
  st_as_sfc() |>
  st_buffer(500000) |>
  st_transform(eck4) |>
  st_bbox()

mexico_bbox <- st_bbox(c(xmin = -106, xmax = -87, ymin = 15, ymax = 23), crs = 4326) |>
  st_transform(albers)

mexico_bedrock <- crop(bedrock, mexico_bbox_eck) |>
  fillNA(method = 'modal') |>
  predict(bedrock_pcs, index = 1:2, cores = 4) |>
  project(albers) |>
  crop(mexico_bbox)
names(mexico_bedrock) <- bedrock_pc_names

mexico_other <- crop(other_predictors, mexico_bbox_eck) |>
  fillNA(method = 'modal') |>
  project(albers) |>
  crop(mexico_bbox)

mexico_climate <- crop(climate, st_transform(mexico_bbox, 4326)) |>
  predict(climate_pcs, index = 1:7, cores = 4) |>
  project(mexico_bedrock)
names(mexico_climate) <- climate_pc_names

mexico_soils <- crop(soils, st_transform(mexico_bbox, crs(soils))) |>
  predict(soil_pcs, index = 1:3, cores = 4) |>
  project(mexico_bedrock) |>
  fillNA('mean')
names(mexico_soils) <- soil_pc_names

mexico_aerosols <- crop(aerosols, st_transform(mexico_bbox, 4326), snap = 'out') |>
  predict(aerosol_pcs, index = 1:5, cores = 4) |>
  project(mexico_bedrock, method = 'cubicspline')
names(mexico_aerosols) <- aerosol_pc_names

mexico_elevation <- crop(elevation, st_transform(mexico_bbox, 4326)) |>
  project(mexico_bedrock)

# Combine all predictors and mask out ocean pixels
mexico_predictors <- c(mexico_bedrock, mexico_other, mexico_aerosols, mexico_climate, mexico_soils, mexico_elevation) %>%
  mask(., .[[14]])

# Check rasters (optional)
plot(mexico_bedrock)
plot(mexico_aerosols)
plot(mexico_climate)
plot(mexico_soils)
plot(mexico_other)
plot(mexico_predictors[[1:12]])
plot(mexico_predictors[[13:21]])

# Save Mexico predictors
saveRDS(mexico_predictors, here('data/derived/mexico_predictors.rds'))

# ==============================================================================
# GLOBAL ISOTOPE OBSERVATIONS
# ==============================================================================

# Dataset 1: FACETS database
pts1 <- here('data/raw/facets-2024-0180supplb.xlsx') |>
  read_excel(sheet = 2, guess_max = 1000000) |>
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) |>
  st_transform(eck4) |>
  select(Sr = `87Sr/86Sr`, type = Type_1, type2 = Type2, material = Material) |>
  mutate(type = tolower(type), type2 = tolower(type2), material = tolower(material)) |>
  filter(type != 'food',
         type != 'insect') |>
  mutate(type = if_else(type %in% c('soil', 'minerals'), 'soil', type),
         type = if_else(type %in% c('mammals', 'mammal', 'reptile', 'invertebrate', 'amphibian', 'biomineral'), 'animal', type))

# Dataset 2: Wang et al. database (remove points near pts1)
not_near <- function(x, y, dist) !st_is_within_distance(x, y, dist)
pts2 <- here('data/raw/Wang et al./data/Wang_et_al_Dataset_S1.xlsx') |>
  read_excel(skip = 2, guess_max = 1000000) |>
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) |>
  st_transform(eck4) |>
  st_filter(st_union(pts1), .predicate = not_near, dist = 500) |>
  select(Sr = `87Sr/86Sr`, material = `Material Type`, type2 = `Taxa of plants and animals`) |>
  mutate(material = tolower(material), type2 = tolower(type2),
         type = if_else(material == 'plant', 'plant', if_else(material == 'Soil leachates', 'soil', 'animal'))) |>
  filter(material != 'insect')

# Dataset 3: Additional database
pts3 <- here('data/raw/1-s2.0-S0031018220302947-mmc1.xlsx') |>
  read_excel(sheet = 2, guess_max = 1000000) |>
  filter(Source.of.data.reference.ID != 144, # location centroids not correct
         Source.of.data.reference.ID != 45) |>  # extreme variation in seasonal water at one site
  filter(`Geolocation Uncertainty` != 'Very High' | is.na(`Geolocation Uncertainty`)) |>
  mutate(Latitude = if_else(Country %in% c('Peru', 'Chile') & Latitude > 0,
                            Latitude * -1, Latitude)) |> # some human samples have wrong hemisphere
  select(Sr = `87Sr/86Sr`, type = `Type 1`, type2 = `Type 2`, type3 = `Type 3`,
         material = Material, Longitude, Latitude) |>
  filter(!is.na(Sr)) |>
  mutate(across(type:material, tolower),
         Sr = str_replace(Sr, ',', '.'),
         Sr = parse_number(Sr)) |>
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) |>
  st_transform(eck4) |>
  filter(!(type %in% c('bird', 'insect', 'dust'))) |>
  mutate(type2 = if_else(is.na(type2) & type3 %in% c('snail', 'snail shell'), 'snail', type2),
         type2 = if_else(is.na(type2), type, type2),
         type2 = if_else(type2 == "–", type, type2),
         type = if_else(type2 == 'human', 'human', type)) |>
  select(-type3)

# Extract water samples
water <- pts3 |>
  filter(type == 'water',
         !(type2 %in% c('snow', 'tap', 'rain', 'throughfall', 'mineral', 'minewater', 'sewage')))

# Human samples from CAMBIO database
cambio <- read_excel(here('data/raw/CAMBIO Humans Database 2024 - Local Sr Only.xlsx'),
                     sheet = 1,  guess_max = 10000) |>
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) |>
  st_transform(eck4) |>
  select(Sr = `87Sr/86Sr`, type2 = Culture) |>
  mutate(type = 'human')

humans <- pts3 |>
  filter(type == 'human') |>
  st_filter(st_union(cambio), .predicate = not_near, dist = 2) |>
  bind_rows(cambio)

# Combine all global isotope observations
pts <- bind_rows(pts1, pts2, water, humans) |>
  mutate(oaxaca = FALSE)

# Oaxaca plant samples ---------------------------------------------------------
oaxaca_plants <- read_excel(here('data/raw/OAX Sr Result Table_Updated.xlsx'), sheet = 2) %>%
  slice_head(n = -5) |> # remove empty rows from end of spreadsheet
  st_as_sf(coords = c('UTM-E', 'UTM-N'), crs = 32614) |>
  st_transform(eck4) |>
  rename(Locality = Site,
         elevation = `Elevation (masl)`,
         rooting_depth = `Plant Rooting Depth`) |>
  filter(Locality != 'San Martin Tilcajete') # remove contaminated site

# Visualize Oaxaca samples
ggplot(oaxaca_plants, aes(Locality, Sr)) +
  geom_point()

# Save oaxaca_plants for visualization
saveRDS(oaxaca_plants, here('data/derived/oaxaca_plants.rds'))

# Combine all isotope observations (global + Oaxaca)
pts_combined <- oaxaca_plants |>
  select(Sr) |>
  mutate(type = 'plant',
         oaxaca = TRUE) |>
  bind_rows(pts) |>
  filter(between(Sr, .703, .78))

# Save combined observations
saveRDS(pts_combined, here('data/derived/pts_combined.rds'))

# Optional: Interactive map view
# mapview::mapview(pts_combined, zcol = 'type')

# ==============================================================================
# EXTRACT PREDICTORS AT OBSERVATION LOCATIONS
# ==============================================================================

# Extract elevation
elevation_extract <- terra::extract(elevation,
                                    st_transform(pts_combined, 4326),
                                    ID = FALSE,
                                    search_radius = 5000,
                                    na.rm = TRUE) |>
  select(elevation)

# Helper function for extraction
extract_fun <- function(x, rast, crs = 4326){
  terra::extract(rast[[x]], st_transform(pts_combined, crs = crs),
                 ID = FALSE,
                 search_radius = 5000,
                 na.rm = TRUE)[,1, drop = FALSE]
}

# Extract and transform climate data
climate_extract <- map_dfc(1:nlyr(climate), extract_fun, climate) %>%
  predict(climate_pcs, .) %>%
  .[,1:7] %>%
  as_tibble() |>
  setNames(climate_pc_names)

# Extract and transform soil data
soils_extract <- map_dfc(1:nlyr(soils), extract_fun, soils, crs(soils)) %>%
  predict(soil_pcs, .) %>%
  .[,1:3] |>
  as_tibble() |>
  setNames(soil_pc_names)

# Extract and transform aerosol data
aerosol_extract <- terra::extract(aerosols,
                                  st_transform(pts_combined, 4326),
                                  ID = FALSE,
                                  method = 'bicubic',
                                  na.rm = TRUE) %>%
  predict(aerosol_pcs, .) %>%
  .[,1:5] |>
  as_tibble() |>
  setNames(aerosol_pc_names)

# Extract and transform bedrock data
bedrock_extract <- map_dfc(1:nlyr(bedrock), extract_fun, bedrock, eck4) %>%
  predict(bedrock_pcs, .) %>%
  .[,1:2] |>
  as_tibble() |>
  setNames(bedrock_pc_names)

# Extract other predictors
other_predictors[['sediment']]  <- as.int(other_predictors[['sediment']])
other_extract <-  map_dfc(1:nlyr(other_predictors), extract_fun, other_predictors, eck4)

# Combine all predictors with observations
dat <- pts_combined |>
  bind_cols(bedrock_extract) |>
  bind_cols(climate_extract) |>
  bind_cols(aerosol_extract) |>
  bind_cols(elevation_extract) |>
  bind_cols(other_extract) |>
  bind_cols(soils_extract) |>
  mutate(type = factor(type)) |>
  mutate(elevation = as.numeric(elevation)) |>
  mutate(sediment = factor(sediment, labels = sediment_tibble$category[1:10])) |>
  select(-geometry) # hack to move geometry to last column

# ==============================================================================
# SPATIAL BINNING
# ==============================================================================
# Use discrete global grid for equal-area spatial bins

dggs <- dgconstruct(res = 4)
grid <- dgearthgrid(dggs) |>
  st_wrap_dateline() |>
  st_transform(eck4)

coords <- dat |>
  st_transform(4326) |>
  st_coordinates()

dat$seqnum <- dgGEO_to_SEQNUM(dggs, coords[,1], coords[,2])$seqnum

# Visualize spatial binning
plot(grid)
ggplot() +
  geom_sf(data = grid) +
  geom_sf(data = dat, aes(color = seqnum))

# Save final modeling dataset
saveRDS(dat, here('data/derived/dat.rds'))

message("Data preprocessing complete!")
message("Created files:")
message("  - data/derived/mexico_predictors.rds")
message("  - data/derived/pts_combined.rds")
message("  - data/derived/oaxaca_plants.rds")
message("  - data/derived/dat.rds")
