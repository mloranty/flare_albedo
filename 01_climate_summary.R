#####################################
# calculate climate vars for study region 
# using ECMWF ERA5-Land
#
# MML 02/18/26
#####################################

rm (list=ls())
# load required packages
library(terra)
library(tidyverse)
library(tidyterra)

# set working directory with options for Mac or PC workstations
ifelse(Sys.info()[1]=="Windows",
       setwd("G:/My Drive/Documents/research/manuscripts/siberia_albedo"),
       setwd("~/Library/CloudStorage/GoogleDrive-mloranty@colgate.edu/My Drive/Documents/research/manuscripts/siberia_albedo"))


clim <- rast("data/ecmwf_era5_land.grib")
crs(clim) <- "EPSG:4326"

fire <- project(vect("data/SiberiaFires2001-2020updated_tree.gpkg"),
                "EPSG:4326")

f.test <- mask(clim,fire)







# 2. Define request
request <- list(
  dataset_short_name = "reanalysis-era5-land",
  product_type = "reanalysis",
  variable = "2m_temperature",
  year = "2023",
  month = "01",
  day = "01",
  time = "12:00",
  format = "netcdf",
  target = "era5_2023_01.nc"
)

# 3. Download data
 file <- wf_request(user = "mloranty@colgate.edu", request = request, transfer = TRUE)
