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
library(dplyr)

# set working directory with options for Mac or PC workstations
ifelse(Sys.info()[1]=="Windows",
       setwd("G:/My Drive/Documents/research/manuscripts/siberia_albedo"),
       setwd("~/Library/CloudStorage/GoogleDrive-mloranty@colgate.edu/My Drive/Documents/research/manuscripts/siberia_albedo"))


# this file contains monthly averages of ERA5-Land climate data
# variables include 2m air temp, snow depth, snow water equivalent, and total precip
# read climate data and define projection
clim <- rast("data/ecmwf_era5_land.grib")
crs(clim) <- "EPSG:4326"

# read fire perimeters to use as mask
# limit to 55-70 latitude
fire <- vect("data/SiberiaFires2001-2020updated_tree.gpkg") %>%
  project("EPSG:4326")%>%
  filter(lat>=55 & 70>=lat) %>%
  buffer(width = 1000000) %>%
  aggregate() %>%
  buffer(width = -1000000)

# get climate variable names
vars <- unique(names(clim))

# create dataframe for summary data
clim.sum <- data.frame(month = 1:12)

#########################
# subset by vars and year
#########################

# 1. Temp
sub <- clim[[which(names(clim)==vars[1])]] 
temp <- sub[[which(year(time(sub))< 2021)]]

m.temp <- tapp(temp, 
               index = month(time(temp)),
               fun = mean, na.rm = T)
# calculate zonal mean temp, convert from k to c, and add to dataframe
clim.sum$temp_C <- t(zonal(m.temp, fire, fun = "mean", na.rm = T)-273.15)

# 2. Snow Depth
#subset by variable and year
sub <- clim[[which(names(clim)==vars[2])]] 
temp <- sub[[which(year(time(sub))< 2021)]]

#ggregate to monthly
m.temp <- tapp(temp, 
               index = month(time(temp)),
               fun = mean, na.rm = T)

# calculate zonal mean temp, convert from m to cm, and add to dataframe
clim.sum$snw_cm <- t(zonal(m.temp, fire, fun = "mean", na.rm = T)*100)

# 3. Snow Water Equivalent
#subset by variable and year
sub <- clim[[which(names(clim)==vars[3])]] 
temp <- sub[[which(year(time(sub))< 2021)]]

#ggregate to monthly
m.temp <- tapp(temp, 
               index = month(time(temp)),
               fun = mean, na.rm = T)

# calculate zonal mean temp, convert from m to cm, and add to dataframe
clim.sum$swe_cm <- t(zonal(m.temp, fire, fun = "mean", na.rm = T)*100)

# 4. Precip
#subset by variable and year
sub <- clim[[which(names(clim)==vars[4])]] 
temp <- sub[[which(year(time(sub))< 2021)]]

#ggregate to monthly
m.temp <- tapp(temp, 
               index = month(time(temp)),
               fun = mean, na.rm = T)

# calculate zonal mean precip, convert from m to cm, and add to dataframe
clim.sum$prcp_cm <- t(zonal(m.temp, fire, fun = "mean", na.rm = T)*100*30)


terra::extract(m.temp,chy)*100
chy <- vect(matrix(c(161.3508, 68.7427), ncol=2), crs="EPSG:4326")
yak <- vect(matrix(c(129.7422, 62.0397), ncol=2), crs="EPSG:4326")

yak <- c(62.0397, 129.7422)
# API code documenting the data request from Climate Data Store

# import cdsapi
# 
# dataset = "reanalysis-era5-land-monthly-means"
# request = {
#   "product_type": ["monthly_averaged_reanalysis"],
#   "variable": [
#     "2m_temperature",
#     "snow_depth",
#     "snow_depth_water_equivalent",
#     "total_precipitation"
#   ],
#   "year": [
#     "2000", "2001", "2002",
#     "2003", "2004", "2005",
#     "2006", "2007", "2008",
#     "2009", "2010", "2011",
#     "2012", "2013", "2014",
#     "2015", "2016", "2017",
#     "2018", "2019", "2020",
#     "2021", "2022", "2023",
#     "2024", "2025", "2026"
#   ],
#   "month": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12"
#   ],
#   "time": ["00:00"],
#   "data_format": "grib",
#   "download_format": "zip"
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()

