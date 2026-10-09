# Climate, vegetation and urbanisation as drivers of malaria in Tanzania

R code for a spatio-temporal Bayesian analysis of monthly malaria cases in the
184 councils of mainland Tanzania (2015-2025). The model links malaria to
temperature, rainfall, humidity, drought (PDSI), wind, vegetation greenness
(EVI) and El Nino (ONI), each with a delayed and non-linear effect over 0-6
months, and tests whether urbanisation and vegetation change the effect of
rainfall and greenness.

Everything is in one script, `malaria_climate_urban_tanzania.R`, which runs
from the raw data to the final figures and tables.

## What the script does

| Part | Content |
|------|---------|
| 1-6  | Builds the dataset: TerraClimate climate (population weighted council means), urbanisation from GHSL and WorldPop, MODIS EVI from Google Earth Engine, ONI from NOAA, and the malaria case data |
| 7    | Descriptive figures, correlations and VIF |
| 8    | Lagged exposures and distributed lag non-linear (DLNM) cross-basis matrices |
| 9-12 | Negative binomial models in R-INLA with seasonal (by ecozone), spatial (BYM2, by year) and trend random effects; driver selection, interaction selection, final model and a post-2025 sensitivity model |
| 13-16| Random effects, exposure-lag-response curves, interaction curves, cases attributable to EVI and rainfall, relative risk and exceedance maps |
| 17   | Calibration (PIT, interval coverage), 6-month holdout forecast and leave-one-month-out cross-validation |

## Data

No data are included in this repository.

- **Malaria cases**: monthly confirmed malaria cases at outpatient departments
  by council, from the national HMIS (DHIS2). These data belong to the
  Ministry of Health and cannot be shared here; requests should go to the
  National Malaria Control Programme. The script expects a file
  `Data/tz_council_opd_monthly.csv` with the columns
  `Region, Council, year, month, OPD_cases, OPD_attendance, population`.
- **Council boundaries**: `Data/tz_council_geometry.shp` (columns `Region`, `name`).
- **TerraClimate** (https://www.climatologylab.org/terraclimate.html):
  downloaded by the script.
- **Population**: GPW / WorldPop rasters through the `geodata` package,
  downloaded by the script.
- **GHSL settlement model (SMOD)** 2015 and 2020, 1 km
  (https://human-settlement.emergency.copernicus.eu): put the two GeoTIFFs in `Data/ghsl_smod/`.
- **MODIS EVI (MOD13A3)**: extracted in Google Earth Engine with `rgee`
  (needs an Earth Engine account).
- **ONI**: read directly from NOAA CPC.

## How to run

1. Install R (4.3 or later) and the packages loaded at the top of the script.
   R-INLA is installed from its own repository (see the comment in Part 1).
2. Put the data files in a folder called `Data/` next to the script.
3. Run the script from that folder. Figures go to `outputs/figures/`, tables
   to `outputs/tables/` and fitted models to `outputs/models/`.

Every model is saved after it is fitted, so the script can be stopped and
started again without refitting. Fitting all models takes several hours on a
normal laptop. The cross-validation in Part 17 refits the model 132 times and
is switched off by default (`run_cross_validation <- FALSE`). Once
`Data/analysis_data.csv` exists you can set `build_dataset <- FALSE` to skip
the data preparation.

## Method reference

The modelling framework follows

Lowe R, Lee SA, O'Reilly KM, et al. Combined effects of hydrometeorological
hazards and urbanisation on dengue risk in Brazil: a spatiotemporal modelling
study. *Lancet Planetary Health* 2021; 5: e209-e219.
https://github.com/drrachellowe/hydromet_dengue

## Contact

Lembris Njotto
