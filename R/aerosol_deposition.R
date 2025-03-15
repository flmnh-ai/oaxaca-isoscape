library(stars)
library(here)
library(dplyr)

process_raster_1 <- function(raster){
  read_ncdf(raster, eps = 0.001) |>
    transmute(dust_wet = DUWT001 + DUWT002 + DUWT003 + DUWT004 + DUWT005,
              dust_dry = DUDP001 + DUDP002 + DUDP003 + DUDP004 + DUDP005,
              seasalt_wet = SSWT001 + SSWT002 + SSWT003 + SSWT004 + SSWT005,
              seasalt_dry = SSDP001 + SSDP002 + SSDP003 + SSDP004 + SSDP005)
}


process_raster_2 <- function(raster){
  read_ncdf(raster, eps = 0.001) |>
    transmute(black_carbon_wet = BCWT001 + BCWT002,
              black_carbon_dry = BCDP001 + BCDP002 ,
              sulfate_dry = SUDP001 + SUDP002 + SUDP003 + SUDP004,
              sulfate_wet = SUWT001 + SUWT002 + SUWT003 + SUWT004,
              organic_carbon_dry = OCDP001 + OCDP002,
              organic_carbon_wet = OCWT001 + OCWT002)
}


files1 <- list.files(here('data/raw/MERRA2'), pattern = '.nc4', full.names = TRUE) |>
  str_subset('DUD') |>
  map(process_raster_1)

combined_stars1 <- do.call(c, c(files1, list(along = "time"))) |>
  st_apply(1:2, mean, na.rm = TRUE)

saveRDS(combined_stars, here('data/derived/deposition.rds'))

files2 <- list.files(here('data/raw/MERRA2'), pattern = '.nc4', full.names = TRUE) |>
  str_subset('BCD') %>%
  map(process_raster_2)

combined_stars2 <- do.call(c, c(files2, list(along = "time"))) |>
  st_apply(1:2, mean, na.rm = TRUE)

combined_stars <- c(combined_stars1, combined_stars2)


saveRDS(combined_stars, here('data/derived/aerosols.rds'))
