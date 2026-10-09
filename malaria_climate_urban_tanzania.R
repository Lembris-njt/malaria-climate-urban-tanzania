################################################################################
#
#  Climate, vegetation and urbanisation as drivers of malaria in Tanzania. A spatio-temporal Bayesian model with distributed lag non-linear terms
#
#  Author : Lembris Njotto
#  Data   : monthly malaria cases for 184 councils, 2014-2025
#
#  The script runs the whole analysis from the start to the end:
#
#    Part 1  Packages, folders and settings
#    Part 2  Council boundaries
#    Part 3  Climate covariates from TerraClimate (population weighted)
#    Part 4  Urbanisation, population density and ecozones
#    Part 5  Vegetation (MODIS EVI) and El Nino (ONI)
#    Part 6  Putting the analysis dataset together
#    Part 7  Descriptive analysis
#    Part 8  Lagged exposures and cross-basis matrices
#    Part 9  Model set-up
#    Part 10 Step 1: which climate drivers?
#    Part 11 Step 2: which interactions?
#    Part 12 Final model and sensitivity model
#    Part 13 Results: random effects
#    Part 14 Results: exposure-lag-response curves
#    Part 15 Results: cases attributable to EVI and rainfall
#    Part 16 Results: relative risk and exceedance maps
#    Part 17 Validation: calibration, holdout and cross-validation
#
#  No data are included. The climate, vegetation and population data are open and the code below shows where to get them.
#
################################################################################


# ==============================================================================
# Part 1. Packages, folders and settings
# ==============================================================================

# INLA is not on CRAN. Install it once with:
# install.packages("INLA",
#   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)

library(INLA)
library(data.table)
library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(exactextractr)
library(ncdf4)
library(geodata)
library(spdep)
library(dlnm)
library(tsModel)
library(ggpubr)
library(RColorBrewer)
library(corrplot)
library(car)

# Folders. Put the raw files in data_dir; everything the script makes goes
# into out_dir.
data_dir <- "Data"
out_dir <- "outputs"
model_dir <- file.path(out_dir, "models")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
for (d in c(data_dir, model_dir, fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Switches. Building the dataset needs the raw files and a Google Earth
# Engine account; once Data/analysis_data.csv exists it can be skipped.
# The leave-one-month-out cross-validation is 132 model fits, so it is off
# by default.
build_dataset <- TRUE
run_cross_validation <- FALSE

years <- 2014:2025


# ==============================================================================
# Part 2. Council boundaries
# ==============================================================================

# Shapefile of the 184 councils (columns: Region, name)
councils <- st_read(file.path(data_dir, "your_area_geometry.shp")) |>
  rename(Council = name) |>
  st_make_valid()


# ==============================================================================
# Part 3. Climate covariates from TerraClimate
# ==============================================================================
#
# Six monthly variables: maximum and minimum temperature (tmax, tmin),
# vapour pressure (vap), Palmer drought severity index (PDSI), precipitation
# (ppt) and wind speed (ws). One NetCDF file per variable and year.
#
# Instead of a simple area mean over each council we use a population
# weighted mean, so that the climate value reflects where people live and
# not empty land such as game reserves. Weights are WorldPop 1 km rasters:
# the 2015 raster for 2014-2017 and the 2020 raster for 2018-2025.

tc_vars <- c("tmax", "tmin", "vap", "PDSI", "ppt", "ws")
tc_dir <- file.path(data_dir, "terraclimate_nc")
dir.create(tc_dir, showWarnings = FALSE)

# Download the files (about 100 MB each, only the ones that are missing)
download_terraclimate <- function(var, yr) {
  f <- file.path(tc_dir, paste0("TerraClimate_", var, "_", yr, ".nc"))
  if (!file.exists(f)) {
    url <- paste0(
      "http://thredds.northwestknowledge.net:8080/thredds/fileServer/",
      "TERRACLIMATE_ALL/data/TerraClimate_", var, "_", yr, ".nc"
    )
    download.file(url, f, mode = "wb")
  }
  f
}

# Population rasters (GPW via the geodata package), cropped to Tanzania
pop_dir <- file.path(data_dir, "worldpop")
dir.create(pop_dir, showWarnings = FALSE)

get_population_raster <- function(yr) {
  f <- file.path(pop_dir, paste0("tza_ppp_", yr, "_1km.tif"))
  if (!file.exists(f)) {
    r <- geodata::population(year = yr, res = 0.5, path = tempdir())
    r <- crop(r, vect(councils))
    writeRaster(r, f, overwrite = TRUE)
  }
  rast(f)
}

# Read one TerraClimate file and turn it into a raster with one layer per
# month. The values are read as raw integers so that missing values can be
# masked before the scale factor is applied.
read_terraclimate <- function(nc_path, var) {
  nc <- nc_open(nc_path)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  origin <- as.Date(str_extract(nc$dim$time$units, "\\d{4}-\\d{2}-\\d{2}"))
  dates <- origin + ncvar_get(nc, "time")

  # Only read the box around Tanzania
  bb <- st_bbox(councils)
  ix <- which(lon >= bb["xmin"] - 0.1 & lon <= bb["xmax"] + 0.1)
  iy <- which(lat >= bb["ymin"] - 0.1 & lat <= bb["ymax"] + 0.1)

  scale <- ncatt_get(nc, var, "scale_factor")$value
  offset <- ncatt_get(nc, var, "add_offset")$value
  fill <- ncatt_get(nc, var, "_FillValue")$value

  vals <- ncvar_get(nc, var,
    start = c(min(ix), min(iy), 1),
    count = c(length(ix), length(iy), length(dates)),
    raw_datavals = TRUE
  )
  nc_close(nc)

  vals[vals == fill] <- NA
  vals <- vals * scale + offset

  lon <- lon[ix]
  lat <- lat[iy]
  half_cell <- 1 / 48
  r <- rast(
    nrows = length(lat), ncols = length(lon), nlyr = length(dates),
    xmin = min(lon) - half_cell, xmax = max(lon) + half_cell,
    ymin = min(lat) - half_cell, ymax = max(lat) + half_cell,
    crs = "EPSG:4326"
  )
  for (i in seq_along(dates)) {
    m <- t(vals[, , i])
    if (lat[1] < lat[length(lat)]) m <- m[nrow(m):1, ]
    values(r[[i]]) <- as.vector(m)
  }
  names(r) <- format(dates, "%Y-%m")
  r
}

# Population weighted mean of every layer for every council.
# weighted mean = sum(value * pop * coverage) / sum(pop * coverage)
# A few tiny boundary cells give NaN, so they are dropped (na.rm).
pop_weighted_mean <- function(r, pop) {
  cells <- exact_extract(r, st_transform(councils, crs(r)),
    fun = NULL, weights = pop,
    include_cols = c("Region", "Council"), progress = FALSE
  )
  layers <- setdiff(names(cells[[1]]), c("Region", "Council", "weight", "coverage_fraction"))

  one_council <- function(d) {
    w <- d$weight * d$coverage_fraction
    means <- sapply(layers, function(l) {
      ok <- !is.na(d[[l]]) & !is.na(w)
      sum(d[[l]][ok] * w[ok]) / sum(w[ok])
    })
    data.frame(Region = d$Region[1], Council = d$Council[1], period = layers, value = means)
  }
  bind_rows(lapply(cells, one_council))
}

if (build_dataset) {
  pop_2015 <- get_population_raster(2015)
  pop_2020 <- get_population_raster(2020)
  pop_2015[is.na(pop_2015)] <- 0
  pop_2020[is.na(pop_2020)] <- 0

  climate_list <- list()
  for (var in tc_vars) {
    for (yr in years) {
      cat("TerraClimate", var, yr, "\n")
      r <- read_terraclimate(download_terraclimate(var, yr), var)
      pop <- if (yr <= 2017) pop_2015 else pop_2020
      climate_list[[paste(var, yr)]] <- pop_weighted_mean(r, pop) |> mutate(variable = var)
    }
  }

  climate <- bind_rows(climate_list) |>
    pivot_wider(names_from = variable, values_from = value) |>
    mutate(
      date = ym(period),
      year = year(date),
      month = month(date)
    ) |>
    rename(pdsi = PDSI) |>
    select(Region, Council, year, month, tmax, tmin, vap, pdsi, ppt, ws)

  # Relative humidity from vapour pressure and mean temperature (Tetens).
  # vap is in kPa, the saturation pressure is in kPa as well.
  climate <- climate |>
    mutate(
      tmean = (tmax + tmin) / 2,
      vap_sat = 0.6108 * exp(17.27 * tmean / (tmean + 237.3)),
      rh_pct = pmin(pmax(vap / vap_sat * 100, 0), 100)
    ) |>
    select(-tmean, -vap_sat)

  cat("Climate rows:", nrow(climate), "(expected", 184 * 144, ")\n")
  print(summary(climate[, c("tmax", "tmin", "pdsi", "ppt", "ws", "rh_pct")]))
}


# ==============================================================================
# Part 4. Urbanisation, population density and ecozones
# ==============================================================================
#
# Urbanisation = % of the council population living in urban cells, using
# the GHSL settlement model (SMOD, 1 km). SMOD classes >= 21 are urban,
# class 10 is water. As for the climate weights, the 2015 epoch is used for
# 2014-2017 and the 2020 epoch for 2018-2025. Download the two SMOD tiles
# for the globe (Mollweide, 1 km) from https://human-settlement.emergency.copernicus.eu

smod_dir <- file.path(data_dir, "ghsl_smod")

urban_mask <- function(epoch) {
  r <- rast(file.path(smod_dir, paste0("GHS_SMOD_E", epoch, "_GLOBE_R2023A_54009_1000_V1_0.tif")))
  r <- crop(r, project(vect(councils), crs(r)))
  ifel(r >= 21, 1, ifel(r == 10, NA, 0))
}

urban_share <- function(pop, mask) {
  mask <- resample(project(mask, crs(pop)), pop, method = "near")
  cnc <- st_transform(councils, crs(pop))
  total <- exact_extract(pop, cnc, fun = "sum", progress = FALSE)
  urban <- exact_extract(pop * mask, cnc, fun = "sum", progress = FALSE)
  data.frame(
    Council = councils$Council,
    urban_pop = urban,
    urban = pmin(pmax(urban / total * 100, 0), 100)
  )
}

# Ecozones used for the seasonal random effect (Here we use Region 1 ... as example of regions names)
ecozones <- tribble(
  ~Region, ~ecozone_name,
  "Region 1", "Northern",
  "Region 2", "Northern",
  "Region 3", "Central",
  "Region 4", "Central",
  "Region 5", "Eastern",
  "Region 6", "Eastern",
  "Region 7", "Southern_Highlands",
  "Region 8", "Southern_Highlands",
  "Region 9", "Southern",
  "Region 10", "Southern",
  "Region 11", "Lake",
  "Region 12", "Lake",
  "Region 13", "Western",
  "Region 14", "Western",
)
zone_order <- c("Northern", "Central", "Eastern", "Southern_Highlands", "Southern", "Lake", "Western")
ecozones$ecozone_code <- match(ecozones$ecozone_name, zone_order)

if (build_dataset) {
  mask_2015 <- urban_mask(2015)
  mask_2020 <- urban_mask(2020)
  urban <- bind_rows(
    urban_share(pop_2015, mask_2015) |> crossing(year = 2014:2017),
    urban_share(pop_2020, mask_2020) |> crossing(year = 2018:2025)
  )

  # Council area in km2 (UTM zone 36S), used for population density
  area <- councils |>
    st_transform(32736) |>
    mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6) |>
    st_drop_geometry() |>
    select(Council, area_km2)
}


# ==============================================================================
# Part 5. Vegetation (MODIS EVI) and El Nino (ONI)
# ==============================================================================
#
# EVI: MODIS MOD13A3 v061, monthly at 1 km, council means computed in Google
# Earth Engine through the rgee package. You need an Earth Engine account.
# EVI is kept as an area mean: it describes the landscape (breeding habitat)
# rather than the place where people live.
#
# ONI: Oceanic Nino Index from NOAA CPC, one national value per month.

if (build_dataset) {
  library(rgee)
  ee_Initialize() # run ee_Authenticate() once on a new computer

  council_ee <- sf_as_ee(councils |> st_transform(4326) |> select(Council))

  evi_images <- ee$ImageCollection("MODIS/061/MOD13A3")$
    select("EVI")$
    filterDate("2014-01-01", "2026-01-01")

  council_evi <- function(img) {
    date <- ee$Date(img$get("system:time_start"))
    img$multiply(0.0001)$ # MODIS scale factor
      reduceRegions(collection = council_ee, reducer = ee$Reducer$mean(), scale = 1000)$
      map(function(f) f$set(list(year = date$get("year"), month = date$get("month"))))
  }

  evi <- ee_as_sf(evi_images$map(council_evi)$flatten(), via = "drive") |>
    st_drop_geometry() |>
    transmute(Council, year = as.integer(year), month = as.integer(month), evi = round(mean, 5))

  # ONI. The table gives 3-month seasons; each season is given to its
  # middle month (DJF -> January, JFM -> February, and so on).
  oni_raw <- read.table("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt", header = TRUE)
  seasons <- c("DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ")
  oni <- oni_raw |>
    transmute(year = as.integer(YR), month = match(SEAS, seasons), oni = as.numeric(ANOM)) |>
    filter(year %in% years)
}


# ==============================================================================
# Part 6. Putting the analysis dataset together
# ==============================================================================
#
# The case file is one row per council and month with the columns
#   Region, Council, year, month, Mal_cases, attendance, population
# where Mal_cases are confirmed malaria cases seen at outpatient departments
# and population is the council population for that year.

analysis_file <- file.path(data_dir, "analysis_data.csv")

if (build_dataset) {
  cases <- read.csv(file.path(data_dir, "your_malaria_data_monthly.csv"))

  data <- cases |>
    left_join(climate, by = c("Region", "Council", "year", "month")) |>
    left_join(urban, by = c("Council", "year")) |>
    left_join(area, by = "Council") |>
    left_join(evi, by = c("Council", "year", "month")) |>
    left_join(oni, by = c("year", "month")) |>
    left_join(ecozones, by = "Region") |>
    left_join(st_drop_geometry(councils)[, c("Council", "council_code")], by = "Council") |>
    mutate(
      pop_density = population / area_km2,
      time = (year - min(years)) * 12 + month
    ) |>
    arrange(time, council_code)

  # Every council must have every month and no missing covariate
  stopifnot(nrow(data) == 184 * 144)
  stopifnot(!anyNA(data[, c("ecozone_code", "urban", "evi", "oni", "tmax", "ppt")]))

  write.csv(data, analysis_file, row.names = FALSE)
}

data <- fread(analysis_file)


# ==============================================================================
# Part 7. Descriptive analysis
# ==============================================================================

# Map objects used by all the figures. The council_index (row number in the
# shapefile) is the index for the spatial random effect, so the adjacency
# graph below and the data must use the same order.
council_map <- councils |>
  st_transform(32736) |>
  st_make_valid() |>
  st_buffer(0) |>
  st_transform(4326) |>
  left_join(distinct(data, Council, ecozone_name, ecozone_code), by = "Council")

ecozone_map <- council_map |>
  group_by(ecozone_name) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

model_data <- data[year > 2014] # 2014 is only used for the lags

# Malaria incidence rate (MIR) per 100,000 people
model_data[, mir := Mal_cases / population * 1e5]

# 7a. National monthly MIR with the main climate variables
national <- model_data[, .(
  MIR = sum(Mal_cases) / sum(population) * 1e5,
  Rainfall = weighted.mean(ppt, population),
  Tmax = weighted.mean(tmax, population),
  RH = weighted.mean(rh_pct, population),
  PDSI = weighted.mean(pdsi, population),
  EVI = mean(evi),
  ONI = mean(oni)
), by = .(year, month)]
national[, date := as.Date(paste(year, month, 1, sep = "-"))]

fig_national <- national |>
  pivot_longer(c(MIR, Rainfall, Tmax, RH, PDSI, EVI, ONI)) |>
  mutate(name = factor(name, levels = c("MIR", "Rainfall", "Tmax", "RH", "PDSI", "EVI", "ONI"))) |>
  ggplot(aes(date, value)) +
  geom_line(colour = "grey20") +
  facet_wrap(~name, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.placement = "outside", strip.background = element_blank())
ggsave(file.path(fig_dir, "fig01_national_time_series.pdf"), fig_national, width = 18, height = 22, units = "cm")

# 7b. Annual MIR per council
annual <- model_data[, .(mir = sum(Mal_cases) / mean(population) * 1e5), by = .(council_code, year)]
fig_annual <- council_map |>
  left_join(annual, by = "council_code") |>
  ggplot() +
  geom_sf(aes(fill = mir), colour = NA) +
  geom_sf(data = ecozone_map, fill = NA, colour = "grey20", linewidth = 0.3) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, trans = "sqrt", name = "Annual MIR\nper 100,000") +
  facet_wrap(~year, ncol = 4) +
  theme_void()
ggsave(file.path(fig_dir, "fig02_annual_mir_maps.pdf"), fig_annual, width = 30, height = 25, units = "cm")

# 7c. Urbanisation map (January 2015)
fig_urban <- council_map |>
  left_join(model_data[time == 13, .(council_code, urban)], by = "council_code") |>
  ggplot() +
  geom_sf(aes(fill = urban), colour = "grey40", linewidth = 0.1) +
  geom_sf(data = ecozone_map, fill = NA, colour = "grey20", linewidth = 0.5) +
  scale_fill_distiller(palette = "Blues", direction = 1, name = "% urban") +
  theme_void()
ggsave(file.path(fig_dir, "fig03_urban_map.pdf"), fig_urban, width = 14, height = 14, units = "cm")

# 7d. Correlation between covariates and variance inflation factors
covariates <- c("tmin", "tmax", "ppt", "rh_pct", "pdsi", "ws", "urban", "evi", "oni")
cor_mat <- cor(model_data[, ..covariates])
write.csv(round(cor_mat, 2), file.path(tab_dir, "covariate_correlation.csv"))

pdf(file.path(fig_dir, "fig04_covariate_correlation.pdf"), width = 7, height = 7)
corrplot(cor_mat,
  method = "color", type = "upper", addCoef.col = "black", number.cex = 0.7,
  col = colorRampPalette(c("#2166AC", "white", "#D73027"))(200)
)
dev.off()

# VIF for the covariates that enter the final model (tmax, not tmin)
vif_fit <- lm(Mal_cases ~ tmax + ppt + rh_pct + pdsi + ws + urban + evi + oni, data = model_data)
print(vif(vif_fit))


# ==============================================================================
# Part 8. Lagged exposures and cross-basis matrices
# ==============================================================================
#
# Each climate variable gets a distributed lag non-linear term over lags 0-6
# months. The lags are made on the full data (from January 2014) so that
# January 2015 already has six months of history; then 2014 is dropped.
#
# Exposure dimension: natural splines with knots at the 10th and 90th
# percentiles (PDSI: fixed knots at -2 and 0, the drought boundaries).
# Lag dimension: natural spline with one knot at lag 3 (PDSI: two equally
# spaced knots).

nlag <- 6

make_lags <- function(x) {
  lags <- tsModel::Lag(x, group = data$council_code, k = 0:nlag)
  lags[data$year > 2014, ]
}

lag_tmin <- make_lags(data$tmin)
lag_tmax <- make_lags(data$tmax)
lag_ppt <- make_lags(data$ppt)
lag_rh <- make_lags(data$rh_pct)
lag_pdsi <- make_lags(data$pdsi)
lag_ws <- make_lags(data$ws)
lag_evi <- make_lags(data$evi)

data <- model_data # from here on only 2015-2025

knots_10_90 <- function(x) quantile(x, c(0.10, 0.90))

cross_basis <- function(lagged, x, name) {
  cb <- crossbasis(lagged,
    argvar = list(fun = "ns", knots = knots_10_90(x)),
    arglag = list(fun = "ns", knots = nlag / 2)
  )
  colnames(cb) <- paste0(name, ".", colnames(cb))
  cb
}

basis_tmin <- cross_basis(lag_tmin, data$tmin, "basis_tmin")
basis_tmax <- cross_basis(lag_tmax, data$tmax, "basis_tmax")
basis_ppt <- cross_basis(lag_ppt, data$ppt, "basis_ppt")
basis_rh <- cross_basis(lag_rh, data$rh_pct, "basis_rh")
basis_ws <- cross_basis(lag_ws, data$ws, "basis_ws")
basis_evi <- cross_basis(lag_evi, data$evi, "basis_evi")

basis_pdsi <- crossbasis(lag_pdsi,
  argvar = list(fun = "ns", knots = c(-2, 0)),
  arglag = list(fun = "ns", knots = equalknots(0:nlag, 2))
)
colnames(basis_pdsi) <- paste0("basis_pdsi.", colnames(basis_pdsi))

# Interaction terms: a cross-basis multiplied by a centred modifier.
# The centring value is where the main cross-basis describes the effect,
# e.g. urban_basis1_* gives the effect at the 75th percentile of urbanisation.
# 1 = 75th, 2 = 50th, 3 = 25th percentile.
centre <- function(x, p) x - quantile(x, p)
pcts <- c(0.75, 0.50, 0.25)

interaction_basis <- function(basis, modifier, name) {
  out <- basis * modifier
  colnames(out) <- paste0(name, ".", colnames(basis))
  out
}

for (k in 1:3) {
  u <- centre(data$urban, pcts[k])
  e <- centre(data$evi, pcts[k])
  assign(paste0("urban_basis", k, "_pdsi"), interaction_basis(basis_pdsi, u, paste0("urban_basis", k, "_pdsi")))
  assign(paste0("urban_basis", k, "_ppt"), interaction_basis(basis_ppt, u, paste0("urban_basis", k, "_ppt")))
  assign(paste0("urban_basis", k, "_evi"), interaction_basis(basis_evi, u, paste0("urban_basis", k, "_evi")))
  assign(paste0("evi_basis", k, "_ppt"), interaction_basis(basis_ppt, e, paste0("evi_basis", k, "_ppt")))
}


# ==============================================================================
# Part 9. Model set-up
# ==============================================================================

# Adjacency graph for the BYM2 spatial effect. Islands (councils with no
# neighbour, e.g. Mafia) are linked to their nearest council.
nb <- poly2nb(as_Spatial(council_map$geometry), snap = 0.01)
islands <- which(card(nb) == 0)
if (length(islands) > 0) {
  centroids <- st_coordinates(st_centroid(st_geometry(council_map)))
  nearest <- knn2nb(knearneigh(centroids, k = 1))
  for (i in islands) {
    j <- nearest[[i]]
    nb[[i]] <- as.integer(j)
    nb[[j]] <- sort(as.integer(c(nb[[j]], i)))
  }
}
graph_file <- file.path(out_dir, "council.graph")
nb2INLA(graph_file, nb)

# Model variables
df <- data.frame(
  Y = data$Mal_cases,
  E = data$population / 1e5, # offset: rate per 100,000
  T1 = data$month, # month of the year (seasonality)
  T2 = data$year - 2014, # year (replicate for the spatial effect)
  T3 = data$time - 12, # month 1..132 (long-term trend)
  S1 = match(data$Council, council_map$Council), # council index (same order as the graph)
  S2 = data$ecozone_code, # ecozone (replicate for seasonality)
  Vu = data$urban,
  Voni = data$oni,
  Vpost2025 = as.integer(data$year >= 2025)
)
E <- df$E

# Penalised complexity prior for all random effect precisions
pc_prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

# Negative binomial model fitted with INLA. full = FALSE is used for the
# screening (faster); full = TRUE keeps what is needed to draw posterior
# samples later.
fit_model <- function(formula, data = df, full = FALSE) {
  inla(formula,
    data = data,
    family = "nbinomial",
    offset = log(E),
    control.inla = list(strategy = "simplified.laplace"),
    control.compute = list(dic = TRUE, config = full, cpo = FALSE, return.marginals = FALSE),
    control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
    control.predictor = list(link = 1, compute = full),
    verbose = FALSE
  )
}

# Fit a model once and keep it on disk, so the script can be re-run
fit_or_load <- function(name, formula, full = FALSE) {
  f <- file.path(model_dir, paste0(name, ".rds"))
  if (file.exists(f)) {
    return(readRDS(f))
  }
  cat("Fitting", name, "...\n")
  m <- fit_model(formula, full = full)
  saveRDS(m, f)
  m
}

# Base model:
#   seasonal curve per ecozone (cyclic RW1)
#   BYM2 spatial effect for each year
#   national long-term trend (RW2)
base_formula <- Y ~ 1 +
  f(T1, replicate = S2, model = "rw1", cyclic = TRUE, constr = TRUE, scale.model = TRUE, hyper = pc_prior) +
  f(S1, model = "bym2", replicate = T2, graph = graph_file, scale.model = TRUE, hyper = pc_prior) +
  f(T3, model = "rw2", constr = TRUE, scale.model = TRUE, hyper = pc_prior)

add_terms <- function(formula, terms) {
  if (length(terms) == 0 || terms == "") return(formula)
  update(formula, as.formula(paste("~ . +", terms)))
}


# ==============================================================================
# Part 10. Step 1: which climate drivers?
# ==============================================================================

step1 <- c(
  "base" = "",
  "tmin" = "basis_tmin",
  "tmax" = "basis_tmax",
  "ppt" = "basis_ppt",
  "rh" = "basis_rh",
  "pdsi" = "basis_pdsi",
  "ws" = "basis_ws",
  "tmin + ppt" = "basis_tmin + basis_ppt",
  "tmax + ppt" = "basis_tmax + basis_ppt",
  "tmin + rh" = "basis_tmin + basis_rh",
  "tmin + pdsi" = "basis_tmin + basis_pdsi",
  "tmin + ppt + rh" = "basis_tmin + basis_ppt + basis_rh",
  "tmin + ppt + pdsi" = "basis_tmin + basis_ppt + basis_pdsi",
  "tmin + rh + pdsi" = "basis_tmin + basis_rh + basis_pdsi",
  "tmin + ppt + rh + pdsi" = "basis_tmin + basis_ppt + basis_rh + basis_pdsi",
  "tmin + ppt + rh + pdsi + ws" = "basis_tmin + basis_ppt + basis_rh + basis_pdsi + basis_ws",
  "tmax + ppt + rh" = "basis_tmax + basis_ppt + basis_rh",
  "tmax + ppt + pdsi" = "basis_tmax + basis_ppt + basis_pdsi",
  "tmax + rh + pdsi" = "basis_tmax + basis_rh + basis_pdsi",
  "tmax + ppt + rh + pdsi" = "basis_tmax + basis_ppt + basis_rh + basis_pdsi",
  "tmax + ppt + rh + pdsi + ws" = "basis_tmax + basis_ppt + basis_rh + basis_pdsi + basis_ws",
  "evi" = "basis_evi",
  "oni" = "Voni",
  "tmin + ppt + rh + pdsi + ws + evi" = "basis_tmin + basis_ppt + basis_rh + basis_pdsi + basis_ws + basis_evi",
  "tmin + ppt + rh + pdsi + ws + oni" = "basis_tmin + basis_ppt + basis_rh + basis_pdsi + basis_ws + Voni",
  "tmin + ppt + rh + pdsi + ws + evi + oni" = "basis_tmin + basis_ppt + basis_rh + basis_pdsi + basis_ws + basis_evi + Voni",
  "tmax + ppt + rh + pdsi + ws + evi" = "basis_tmax + basis_ppt + basis_rh + basis_pdsi + basis_ws + basis_evi",
  "tmax + ppt + rh + pdsi + ws + oni" = "basis_tmax + basis_ppt + basis_rh + basis_pdsi + basis_ws + Voni",
  "tmax + ppt + rh + pdsi + ws + evi + oni" = "basis_tmax + basis_ppt + basis_rh + basis_pdsi + basis_ws + basis_evi + Voni"
)

dic_step1 <- data.frame(model = names(step1), DIC = NA)
for (i in seq_along(step1)) {
  m <- fit_or_load(paste0("step1_", i - 1), add_terms(base_formula, step1[i]))
  dic_step1$DIC[i] <- round(m$dic$dic)
}
dic_step1$dDIC_vs_base <- dic_step1$DIC - dic_step1$DIC[1]
write.csv(dic_step1, file.path(tab_dir, "step1_dic.csv"), row.names = FALSE)
print(dic_step1[order(dic_step1$DIC), ])

# Keep tmax or tmin, whichever gives the lower DIC in the full model.
# In our data tmax was better.
dic_full_tmax <- dic_step1$DIC[dic_step1$model == "tmax + ppt + rh + pdsi + ws + evi + oni"]
dic_full_tmin <- dic_step1$DIC[dic_step1$model == "tmin + ppt + rh + pdsi + ws + evi + oni"]
temp_term <- if (dic_full_tmax < dic_full_tmin) "basis_tmax" else "basis_tmin"

climate_terms <- paste(temp_term, "+ basis_ppt + basis_rh + basis_pdsi + basis_ws + basis_evi + Voni")
climate_formula <- add_terms(base_formula, climate_terms)


# ==============================================================================
# Part 11. Step 2: which interactions?
# ==============================================================================
#
# Does urbanisation change the effect of drought, rainfall or vegetation,
# and does vegetation change the effect of rainfall? Each interaction is
# tried with the modifier centred at its 75th, 50th and 25th percentile,
# then the best ones are combined.

step2 <- c(
  "PDSI x urban (Q75)" = "urban_basis1_pdsi + Vu",
  "PDSI x urban (Q50)" = "urban_basis2_pdsi + Vu",
  "PDSI x urban (Q25)" = "urban_basis3_pdsi + Vu",
  "PPT x urban (Q75)" = "urban_basis1_ppt + Vu",
  "PPT x urban (Q50)" = "urban_basis2_ppt + Vu",
  "PPT x urban (Q25)" = "urban_basis3_ppt + Vu",
  "EVI x urban (Q75)" = "urban_basis1_evi + Vu",
  "EVI x urban (Q50)" = "urban_basis2_evi + Vu",
  "EVI x urban (Q25)" = "urban_basis3_evi + Vu",
  "EVI x PPT (EVI Q75)" = "evi_basis1_ppt",
  "EVI x PPT (EVI Q50)" = "evi_basis2_ppt",
  "EVI x PPT (EVI Q25)" = "evi_basis3_ppt",
  "PPT x urban + EVI x PPT" = "urban_basis3_ppt + evi_basis3_ppt + Vu",
  "EVI x urban + PDSI x urban" = "urban_basis2_evi + urban_basis2_pdsi + Vu",
  "EVI x urban + PPT x urban + EVI x PPT" = "urban_basis1_evi + urban_basis3_ppt + evi_basis3_ppt + Vu"
)

dic_climate <- min(dic_full_tmax, dic_full_tmin)
dic_step2 <- data.frame(model = names(step2), DIC = NA)
for (i in seq_along(step2)) {
  m <- fit_or_load(paste0("step2_", i), add_terms(climate_formula, step2[i]))
  dic_step2$DIC[i] <- round(m$dic$dic)
}
dic_step2$dDIC_vs_climate_model <- dic_step2$DIC - dic_climate
write.csv(dic_step2, file.path(tab_dir, "step2_dic.csv"), row.names = FALSE)
print(dic_step2)

best <- which.min(dic_step2$DIC)
cat("Selected model:", dic_step2$model[best], "\n")


# ==============================================================================
# Part 12. Final model and sensitivity model
# ==============================================================================
#
# The final model has all three interactions:
#   EVI x urbanisation  (centred at urban Q75)
#   PPT x urbanisation  (centred at urban Q25)
#   PPT x EVI           (centred at EVI Q25)

final_terms <- step2[best]
final_formula <- add_terms(climate_formula, final_terms)
model_final <- fit_or_load("final_model", final_formula, full = TRUE)

summary(model_final)

fixed <- model_final$summary.fixed
write.csv(round(fixed, 4), file.path(tab_dir, "final_fixed_effects.csv"))
hyper <- data.frame(parameter = rownames(model_final$summary.hyperpar), model_final$summary.hyperpar, check.names = FALSE)
write.csv(hyper, file.path(tab_dir, "final_hyperparameters.csv"), row.names = FALSE)

cat(sprintf(
  "Urbanisation main effect: RR %.3f (%.3f - %.3f) per 1%% urban\n",
  exp(fixed["Vu", "mean"]), exp(fixed["Vu", "0.025quant"]), exp(fixed["Vu", "0.975quant"])
))

# Sensitivity: is there an extra jump in 2025 not explained by the model?
model_post2025 <- fit_or_load("final_model_post2025", add_terms(final_formula, "Vpost2025"))
p25 <- model_post2025$summary.fixed["Vpost2025", ]
cat(sprintf(
  "Post-2025 indicator: RR %.3f (%.3f - %.3f), dDIC = %+.0f\n",
  exp(p25$mean), exp(p25$`0.025quant`), exp(p25$`0.975quant`),
  model_post2025$dic$dic - model_final$dic$dic
))


# ==============================================================================
# Part 13. Results: random effects
# ==============================================================================

nyear <- length(unique(df$T2))
ncouncil <- nrow(council_map)

# Seasonal effect for each ecozone
season <- model_final$summary.random$T1 |>
  mutate(
    month = rep(1:12, times = 7),
    ecozone_name = rep(zone_order, each = 12)
  )

fig_season <- ggplot(season, aes(month)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "#2166AC", alpha = 0.3) +
  geom_line(aes(y = mean), colour = "#2166AC") +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  facet_wrap(~ecozone_name, ncol = 4) +
  labs(x = "Month", y = "Contribution to log(MIR)") +
  theme_bw()
ggsave(file.path(fig_dir, "fig05_seasonal_effect.pdf"), fig_season, width = 22, height = 10, units = "cm")

# Spatial effect for each year. For BYM2 INLA returns, per year, first the
# combined effect for every council and then the structured part only;
# we map the combined effect.
spatial <- model_final$summary.random$S1 |>
  mutate(
    year = rep(2015:2025, each = 2 * ncouncil),
    part = rep(rep(c("combined", "structured"), each = ncouncil), times = nyear),
    council_index = rep(rep(1:ncouncil, 2), times = nyear)
  ) |>
  filter(part == "combined")
spatial$council_code <- council_map$council_code[spatial$council_index]

lim <- max(abs(spatial$mean))
fig_spatial <- council_map |>
  left_join(spatial, by = "council_code") |>
  ggplot() +
  geom_sf(aes(fill = mean), colour = NA) +
  geom_sf(data = ecozone_map, fill = NA, colour = "grey40", linewidth = 0.3) +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu")), limits = c(-lim, lim), name = "Contribution\nto log(MIR)") +
  facet_wrap(~year, ncol = 4) +
  theme_void()
ggsave(file.path(fig_dir, "fig06_spatial_effect_by_year.pdf"), fig_spatial, width = 30, height = 25, units = "cm")

# National trend
trend <- model_final$summary.random$T3 |>
  mutate(date = as.Date("2015-01-01") %m+% months(ID - 1))
fig_trend <- ggplot(trend, aes(date)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "grey50", alpha = 0.3) +
  geom_line(aes(y = mean)) +
  labs(x = NULL, y = "Contribution to log(MIR)") +
  theme_bw()
ggsave(file.path(fig_dir, "fig07_national_trend.pdf"), fig_trend, width = 18, height = 9, units = "cm")


# ==============================================================================
# Part 14. Results: exposure-lag-response curves
# ==============================================================================
#
# We draw 500 samples from the joint posterior of the final model and use
# the mean and covariance of the cross-basis coefficients in crosspred().
# Curves show the relative risk (RR) summed over lags 0-6 compared with a
# reference value: the median for most drivers, the 10th percentile for
# rainfall and 0 for PDSI.

set.seed(20260101)
post <- inla.posterior.sample(500, model_final)
latent_names <- rownames(post[[1]]$latent)

coef_samples <- function(prefix) {
  idx <- grep(paste0("^", prefix, "[.]"), latent_names)
  sapply(post, function(s) s$latent[idx, 1])
}

predict_basis <- function(basis, samples, at, cen, bylag = NULL) {
  args <- list(basis,
    coef = rowMeans(samples), vcov = cov(t(samples)),
    model.link = "log", at = at, cen = cen
  )
  if (!is.null(bylag)) args$bylag <- bylag
  do.call(crosspred, args)
}

temp <- sub("basis_", "", temp_term)
temp_label <- if (temp == "tmax") "Maximum temperature (°C)" else "Minimum temperature (°C)"

drivers <- list(
  list(name = temp, x = data[[temp]], ref = median(data[[temp]]), top = 0.99, label = temp_label, col = "#D73027"),
  list(name = "ppt", x = data$ppt, ref = quantile(data$ppt, 0.10), top = 0.999, label = "Precipitation (mm/month)", col = "#2166AC"),
  list(name = "rh", x = data$rh_pct, ref = median(data$rh_pct), top = 0.99, label = "Relative humidity (%)", col = "#1A9641"),
  list(name = "pdsi", x = data$pdsi, ref = 0, top = NA, label = "PDSI (drought index)", col = "#FF7F00"),
  list(name = "ws", x = data$ws, ref = median(data$ws), top = 0.99, label = "Wind speed (m/s)", col = "#762A83"),
  list(name = "evi", x = data$evi, ref = median(data$evi), top = 0.99, label = "EVI (vegetation greenness)", col = "#8B4513")
)

exposure_range <- function(d, n) {
  if (d$name == "pdsi") return(seq(-10, 10, length.out = n))
  seq(quantile(d$x, 0.01), quantile(d$x, d$top), length.out = n)
}

plot_cumulative <- function(cp, d) {
  ggplot(data.frame(x = cp$predvar, rr = cp$allRRfit, lo = cp$allRRlow, hi = cp$allRRhigh), aes(x)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = d$col, alpha = 0.2) +
    geom_line(aes(y = rr), colour = d$col, linewidth = 0.8) +
    geom_rug(data = data.frame(x = sample(d$x, 600)), sides = "b", alpha = 0.1) +
    scale_y_log10() +
    labs(x = d$label, y = "RR (lags 0-6)") +
    theme_bw()
}

plot_surface <- function(cp, ylab) {
  s <- expand.grid(lag = seq(0, nlag, 0.25), exposure = cp$predvar)
  s$rr <- as.vector(t(cp$matRRfit))
  lim <- max(abs(log(range(s$rr))))
  ggplot(s, aes(lag, exposure)) +
    geom_tile(aes(fill = log(rr))) +
    geom_contour(aes(z = rr), breaks = 1, colour = "grey30", linetype = "dashed") +
    scale_fill_gradientn(
      colours = rev(brewer.pal(11, "RdBu")), limits = c(-lim, lim),
      name = "RR", labels = function(b) round(exp(b), 2)
    ) +
    scale_x_continuous(breaks = 0:nlag) +
    labs(x = "Lag (months)", y = ylab) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

cumulative_plots <- list()
surface_plots <- list()
for (d in drivers) {
  prefix <- paste0("basis_", d$name)
  b <- coef_samples(prefix)
  cb <- get(prefix)
  cp <- predict_basis(cb, b, exposure_range(d, 100), d$ref)
  cp_lag <- predict_basis(cb, b, exposure_range(d, 80), d$ref, bylag = 0.25)
  cumulative_plots[[d$name]] <- plot_cumulative(cp, d)
  surface_plots[[d$name]] <- plot_surface(cp_lag, d$label)
}

ggsave(file.path(fig_dir, "fig08_cumulative_rr.pdf"),
  ggarrange(plotlist = cumulative_plots, ncol = 3, nrow = 2, labels = LETTERS[1:6]),
  width = 32, height = 18, units = "cm"
)
ggsave(file.path(fig_dir, "fig09_lag_surfaces.pdf"),
  ggarrange(plotlist = surface_plots, ncol = 3, nrow = 2, labels = LETTERS[1:6]),
  width = 32, height = 22, units = "cm"
)

# Interactions. With an interaction the coefficients for a given level of
# the modifier are   beta_main + (level - centring value) * beta_interaction
urban_q <- quantile(data$urban, c(0.25, 0.50, 0.75))
evi_q <- quantile(data$evi, c(0.25, 0.50, 0.75))
urban_labels <- paste0(c("Rural", "Intermediate", "Highly urban"), " (", round(urban_q, 1), "% urban)")
evi_labels <- paste0(c("Sparse", "Moderate", "Dense"), " vegetation (EVI ", round(evi_q, 2), ")")

b_evi <- coef_samples("basis_evi")
b_evi_urban <- coef_samples("urban_basis1_evi") # centred at urban Q75
b_ppt <- coef_samples("basis_ppt")
b_ppt_urban <- coef_samples("urban_basis3_ppt") # centred at urban Q25
b_ppt_evi <- coef_samples("evi_basis3_ppt") # centred at EVI Q25

curves_by_level <- function(basis, samples_list, labels, at, ref) {
  bind_rows(lapply(seq_along(labels), function(j) {
    cp <- predict_basis(basis, samples_list[[j]], at, ref)
    data.frame(x = cp$predvar, rr = cp$allRRfit, lo = cp$allRRlow, hi = cp$allRRhigh, level = labels[j])
  })) |>
    mutate(level = factor(level, levels = labels))
}

plot_levels <- function(curves, xlab, cols, legend) {
  ggplot(curves, aes(x, colour = level, fill = level)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(aes(y = rr), linewidth = 0.9) +
    scale_y_log10() +
    scale_colour_manual(values = cols, name = legend) +
    scale_fill_manual(values = cols, name = legend) +
    labs(x = xlab, y = "RR (lags 0-6)") +
    theme_bw() +
    theme(legend.position = "bottom", legend.direction = "vertical")
}

at_evi <- seq(quantile(data$evi, 0.01), quantile(data$evi, 0.99), length.out = 100)
at_ppt <- seq(quantile(data$ppt, 0.01), quantile(data$ppt, 0.999), length.out = 100)
ref_evi <- median(data$evi)
ref_ppt <- quantile(data$ppt, 0.10)

# EVI effect by urbanisation
evi_by_urban <- curves_by_level(
  basis_evi,
  lapply(urban_q, function(u) b_evi + (u - urban_q[3]) * b_evi_urban),
  urban_labels, at_evi, ref_evi
)

# Rainfall effect by urbanisation (at EVI Q25) and by vegetation (at urban Q25)
ppt_by_urban <- curves_by_level(
  basis_ppt,
  lapply(urban_q, function(u) b_ppt + (u - urban_q[1]) * b_ppt_urban),
  urban_labels, at_ppt, ref_ppt
)
ppt_by_evi <- curves_by_level(
  basis_ppt,
  lapply(evi_q, function(e) b_ppt + (e - evi_q[1]) * b_ppt_evi),
  evi_labels, at_ppt, ref_ppt
)

urban_cols <- c("#1A9641", "#FD8D3C", "#D73027")
evi_cols <- c("#D4B483", "#5DAD5C", "#1A4D00")
fig_interactions <- ggarrange(
  plot_levels(evi_by_urban, "EVI (vegetation greenness)", urban_cols, "Urbanisation"),
  plot_levels(ppt_by_urban, "Precipitation (mm/month)", urban_cols, "Urbanisation"),
  plot_levels(ppt_by_evi, "Precipitation (mm/month)", evi_cols, "Vegetation"),
  ncol = 3, labels = c("A", "B", "C")
)
ggsave(file.path(fig_dir, "fig10_interactions.pdf"), fig_interactions, width = 36, height = 14, units = "cm")

# EVI lag surfaces at the three urbanisation levels
evi_surfaces <- lapply(urban_q, function(u) {
  cp <- predict_basis(basis_evi, b_evi + (u - urban_q[3]) * b_evi_urban,
    seq(quantile(data$evi, 0.01), quantile(data$evi, 0.99), length.out = 80), ref_evi,
    bylag = 0.25
  )
  plot_surface(cp, "EVI (vegetation greenness)")
})
ggsave(file.path(fig_dir, "fig11_evi_surfaces_by_urban.pdf"),
  ggarrange(plotlist = evi_surfaces, ncol = 3, labels = urban_labels, font.label = list(size = 9)),
  width = 36, height = 12, units = "cm"
)


# ==============================================================================
# Part 15. Results: cases attributable to EVI and rainfall
# ==============================================================================
#
# Backward attributable number (Gasparrini and Leone 2014). For each
# council-month, the contribution of the exposure history over lags 0-6
# compared with the reference is
#     c = (w - w_ref) %*% beta
# where w is the cross-basis row of that council-month and w_ref the row for
# a constant exposure at the reference. beta includes the interactions, so it
# is different for each council-month. Attributable cases = y * (1 - exp(-c))
# summed over council-months with c > 0. Uncertainty from the 500 samples.

reference_row <- function(basis, ref) {
  var_basis <- do.call(onebasis, c(list(x = ref), attr(basis, "argvar")))
  lag_basis <- do.call(onebasis, c(list(x = 0:nlag), attr(basis, "arglag")))
  as.vector(kronecker(var_basis, t(colSums(lag_basis))))
}

attributable <- function(basis, ref, contribution) {
  W <- unclass(basis) - matrix(reference_row(basis, ref), nrow(basis), ncol(basis), byrow = TRUE)
  y <- data$Mal_cases
  an <- sapply(1:500, function(s) {
    c_it <- contribution(W, s)
    pos <- c_it > 0
    sum(y[pos] * (1 - exp(-c_it[pos])))
  })
  c(
    cases = mean(an),
    lower = unname(quantile(an, 0.025)),
    upper = unname(quantile(an, 0.975)),
    percent = mean(an) / sum(y) * 100
  )
}

u75 <- centre(data$urban, 0.75)
u25 <- centre(data$urban, 0.25)
e25 <- centre(data$evi, 0.25)

attr_evi <- attributable(basis_evi, ref_evi, function(W, s) {
  as.vector(W %*% b_evi[, s]) + u75 * as.vector(W %*% b_evi_urban[, s])
})
attr_ppt <- attributable(basis_ppt, ref_ppt, function(W, s) {
  as.vector(W %*% b_ppt[, s]) + u25 * as.vector(W %*% b_ppt_urban[, s]) + e25 * as.vector(W %*% b_ppt_evi[, s])
})

attr_table <- data.frame(
  driver = c("EVI above the median", "Rainfall above the 10th percentile"),
  round(rbind(attr_evi, attr_ppt), 1)
)
write.csv(attr_table, file.path(tab_dir, "attributable_cases.csv"), row.names = FALSE)
print(attr_table)

rm(post)
gc()


# ==============================================================================
# Part 16. Results: relative risk and exceedance maps
# ==============================================================================
#
# 500 joint posterior samples of the expected count for every council-month
# (the linear predictor includes the offset, so exp() gives counts).
# Relative risk per council and year uses internal standardisation:
#   expected = population x national rate of that year
#   RR = sum of fitted cases / sum of expected cases
# Exceedance = share of samples with RR > 1.5.

set.seed(2025)
pred_samples <- inla.posterior.sample(500, model_final, selection = list(Predictor = 0))
mu <- sapply(pred_samples, function(s) exp(s$latent[, 1]))
rm(pred_samples)
gc()

nb_size <- hyper$mean[grepl("size for the nbinomial", hyper$parameter)]

rr_data <- data[, .(council_code, year, cases = Mal_cases, pop = population)]
rr_data[, rate := sum(cases) / sum(pop), by = year]
rr_data[, expected := pop * rate]

group <- paste(rr_data$council_code, rr_data$year)
rr_samples <- rowsum(mu, group) / as.vector(rowsum(rr_data$expected, group))

rr <- rr_data[, .(observed = sum(cases), expected = sum(expected)), by = .(council_code, year)]
rr <- rr[match(rownames(rr_samples), paste(council_code, year))]
rr[, `:=`(
  rr_mean = rowMeans(rr_samples),
  rr_low = apply(rr_samples, 1, quantile, 0.025),
  rr_high = apply(rr_samples, 1, quantile, 0.975),
  p_exceed = rowMeans(rr_samples > 1.5),
  sir = observed / expected
)]
write.csv(rr, file.path(tab_dir, "rr_by_council_year.csv"), row.names = FALSE)
cat("Correlation between model RR and observed SIR:", round(cor(rr$rr_mean, rr$sir), 3), "\n")

rr_map <- council_map |> left_join(rr, by = "council_code")

fig_rr <- ggplot(rr_map) +
  geom_sf(aes(fill = rr_mean), colour = NA) +
  geom_sf(data = ecozone_map, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu")), limits = c(0, 5), oob = scales::squish, name = "RR") +
  facet_wrap(~year, ncol = 4) +
  theme_void()
ggsave(file.path(fig_dir, "fig12_relative_risk_by_year.pdf"), fig_rr, width = 30, height = 25, units = "cm")

fig_exceed <- ggplot(rr_map) +
  geom_sf(aes(fill = p_exceed), colour = NA) +
  geom_sf(data = ecozone_map, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0, 1), name = "P(RR > 1.5)") +
  facet_wrap(~year, ncol = 4) +
  theme_void()
ggsave(file.path(fig_dir, "fig13_exceedance_by_year.pdf"), fig_exceed, width = 30, height = 25, units = "cm")

# Observed vs fitted national MIR with 95% credible interval
national_fit <- rowsum(mu, data$time) / as.vector(rowsum(data$population / 1e5, data$time))
national_obs <- data[, .(obs = sum(Mal_cases) / sum(population) * 1e5), by = time]
national_obs[, `:=`(
  fit = rowMeans(national_fit),
  lo = apply(national_fit, 1, quantile, 0.025),
  hi = apply(national_fit, 1, quantile, 0.975),
  date = as.Date("2015-01-01") %m+% months(time - 13)
)]

fig_fit <- ggplot(national_obs, aes(date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#D37295", alpha = 0.3) +
  geom_line(aes(y = obs, colour = "Observed")) +
  geom_line(aes(y = fit, colour = "Fitted")) +
  scale_colour_manual(values = c(Observed = "#499894", Fitted = "#D37295"), name = NULL) +
  labs(x = NULL, y = "MIR per 100,000") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig14_observed_vs_fitted.pdf"), fig_fit, width = 18, height = 10, units = "cm")


# ==============================================================================
# Part 17. Validation: calibration, holdout and cross-validation
# ==============================================================================

# 17a. Calibration (in sample). Randomised PIT for counts (Czado et al. 2009)
# and coverage of the 95% and 50% predictive intervals.
y <- data$Mal_cases
F_upper <- rowMeans(pnbinom(y, mu = mu, size = nb_size))
F_lower <- rowMeans(pnbinom(y - 1, mu = mu, size = nb_size))
set.seed(1)
pit <- F_lower + runif(length(y)) * (F_upper - F_lower)

set.seed(2)
y_sim <- matrix(rnbinom(length(mu), mu = mu, size = nb_size), nrow = nrow(mu))
q <- apply(y_sim, 1, quantile, c(0.025, 0.25, 0.75, 0.975))

calibration <- data.frame(
  measure = c("Coverage 95% interval", "Coverage 50% interval", "PIT mean (ideal 0.5)", "PIT sd (ideal 0.289)"),
  value = round(c(
    mean(y >= q[1, ] & y <= q[4, ]), mean(y >= q[2, ] & y <= q[3, ]),
    mean(pit), sd(pit)
  ), 3)
)
write.csv(calibration, file.path(tab_dir, "calibration.csv"), row.names = FALSE)
print(calibration)

fig_pit <- ggplot(data.frame(pit), aes(pit)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, 0.05), fill = "#9E9AC8", colour = "white") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Randomised PIT", y = "Density") +
  theme_bw()
ggsave(file.path(fig_dir, "fig15_pit_histogram.pdf"), fig_pit, width = 12, height = 8, units = "cm")

rm(y_sim)
gc()

# 17b. Holdout: hide the last 6 months (all councils), refit, and compare the
# forecast with what was observed.
holdout <- df$T3 > max(df$T3) - 6
df_holdout <- df
df_holdout$Y[holdout] <- NA

holdout_file <- file.path(model_dir, "holdout_model.rds")
if (file.exists(holdout_file)) {
  model_holdout <- readRDS(holdout_file)
} else {
  model_holdout <- inla(final_formula,
    data = df_holdout, family = "nbinomial", offset = log(E),
    control.inla = list(strategy = "simplified.laplace"),
    control.fixed = list(prec.intercept = 1, prec = 1),
    control.predictor = list(link = 1, compute = TRUE)
  )
  saveRDS(model_holdout, holdout_file)
}

size_holdout <- model_holdout$summary.hyperpar[grep("size", rownames(model_holdout$summary.hyperpar)), "mean"]
pred <- data.frame(
  time = data$time[holdout], E = E[holdout],
  observed = y[holdout],
  forecast = model_holdout$summary.fitted.values$mean[holdout]
)

# National 95% prediction interval by simulation
set.seed(99)
sims <- matrix(rnbinom(nrow(pred) * 1000, mu = pred$forecast, size = size_holdout), ncol = 1000)
sims_mir <- rowsum(sims, pred$time) / as.vector(rowsum(pred$E, pred$time))

holdout_national <- pred |>
  group_by(time) |>
  summarise(obs_mir = sum(observed) / sum(E), fit_mir = sum(forecast) / sum(E), E = sum(E)) |>
  mutate(
    lo = apply(sims_mir, 1, quantile, 0.025),
    hi = apply(sims_mir, 1, quantile, 0.975),
    ratio = fit_mir / obs_mir
  )
print(holdout_national)

pred$lo <- qnbinom(0.025, mu = pred$forecast, size = size_holdout)
pred$hi <- qnbinom(0.975, mu = pred$forecast, size = size_holdout)

holdout_metrics <- data.frame(
  measure = c("National MIR MAE", "National MIR RMSE", "National 95% coverage", "Council-month MAE (cases)", "Council-month 95% coverage"),
  value = round(c(
    mean(abs(holdout_national$obs_mir - holdout_national$fit_mir)),
    sqrt(mean((holdout_national$obs_mir - holdout_national$fit_mir)^2)),
    mean(holdout_national$obs_mir >= holdout_national$lo & holdout_national$obs_mir <= holdout_national$hi),
    mean(abs(pred$observed - pred$forecast)),
    mean(pred$observed >= pred$lo & pred$observed <= pred$hi)
  ), 3)
)
write.csv(holdout_metrics, file.path(tab_dir, "holdout_metrics.csv"), row.names = FALSE)
print(holdout_metrics)

fig_holdout <- ggplot(holdout_national, aes(time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#D37295", alpha = 0.3) +
  geom_line(aes(y = obs_mir, colour = "Observed")) +
  geom_line(aes(y = fit_mir, colour = "Forecast")) +
  scale_colour_manual(values = c(Observed = "#499894", Forecast = "#D37295"), name = NULL) +
  labs(x = "Month index", y = "MIR per 100,000") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig16_holdout_forecast.pdf"), fig_holdout, width = 14, height = 9, units = "cm")

# 17c. Leave-one-month-out cross-validation (as in Lowe et al. 2021).
# For each of the 132 months, hide that month in all councils, refit and
# keep the prediction. Each fold is saved, so the loop can be stopped and
# continued later. Hyperparameters start at the final model's mode and are
# fixed there (empirical Bayes) to save time.
if (run_cross_validation) {
  cv_dir <- file.path(model_dir, "cv")
  dir.create(cv_dir, showWarnings = FALSE)
  theta_start <- model_final$mode$theta

  for (k in sort(unique(df$T3))) {
    fold_file <- file.path(cv_dir, sprintf("fold_%03d.rds", k))
    if (file.exists(fold_file)) next
    left_out <- df$T3 == k
    df_k <- df
    df_k$Y[left_out] <- NA
    m_k <- inla(final_formula,
      data = df_k, family = "nbinomial", offset = log(E),
      control.inla = list(strategy = "simplified.laplace", int.strategy = "eb"),
      control.mode = list(theta = theta_start, restart = TRUE),
      control.fixed = list(prec.intercept = 1, prec = 1),
      control.predictor = list(link = 1, compute = TRUE)
    )
    saveRDS(list(rows = which(left_out), pred = m_k$summary.fitted.values$mean[left_out]), fold_file)
    cat("Fold", k, "done\n")
  }

  cv_pred <- rep(NA_real_, nrow(df))
  for (f in list.files(cv_dir, full.names = TRUE)) {
    r <- readRDS(f)
    cv_pred[r$rows] <- r$pred
  }

  in_sample <- rowMeans(mu)
  cv_summary <- data.frame(
    measure = c("MAE in sample (cases)", "MAE cross-validated (cases)", "Correlation cross-validated"),
    value = round(c(
      mean(abs(y - in_sample)),
      mean(abs(y - cv_pred), na.rm = TRUE),
      cor(y, cv_pred, use = "complete.obs")
    ), 3)
  )
  write.csv(cv_summary, file.path(tab_dir, "cross_validation.csv"), row.names = FALSE)
  print(cv_summary)
}
