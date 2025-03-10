#####################################
# preprocess MODIS albedo for fire perimeters
# from Talucci et al 2022
# albedo provided from Google Earth Engine
# by E. Webb
#
# MML 06/21/22
#####################################

rm (list=ls())
# load required packages
library(terra)
library(lubridate)
library(tidyverse)
#library(ncdf4)

# set working directory with options for Mac or PC workstations
ifelse(Sys.info()[1]=="Windows",
       setwd("G:/My Drive/Documents/research/manuscripts/siberia_albedo"),
       setwd("~/Library/CloudStorage/GoogleDrive-mloranty@colgate.edu/My Drive/Documents/research/manuscripts/siberia_albedo"))

#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#
# READ AND PRE-PROCESS FIRE PERIMETER DATABASE
# add new variables to fire database for analysis (biome, etc.)
# expand to see code, otherwise skip and read updated shapefile
############################################################################

# read fire perimeter data set
fire <- vect("data/SiberiaFires2001-2020/SiberiaFires2001-2020.shp")

##add biome var to fire database
fire$biome <- ifelse(fire$EcoCode == "NST"|fire$EcoCode == "EST","Boreal","Tundra")

## get centroids for fire polygons to add lat/lon
c <- centroids(fire)

#project to lat/lon
d <- project(c,"EPSG:4326")

# get coordinates
l <- crds(d)

# add latitude and longitude to fire database
fire$lat <- l[,2]
fire$lon <- l[,1]
rm(c,d,l)

# create dNBR dummay variable
fire$dnbr <- -9999

## read dNBR data from Anna Talucci and add to shapefile
d <- list.files(path = "data/dnbr", pattern = ".csv", full.names = T)
for(i in 1:length(d))
{
  n <- read.csv(d[i], header = T)
  r <- match(n$UniqueId, fire$UniqueId)
  fire$dnbr[r] <- n$mean
}
rm(d,n,r)

# read CACK radiative kernel data set
# transpose and correct longitude values as well
rf <- rotate(t(rast("data/CACKv1.0/CACKv1.0.nc", "CACK")))
# correct set latitude to correct extent
ext(rf) <- c(-180,180,-90,90)
#crs(rf) <- "EPSG:4326"
# reproject fire to lat/lon in order to extract CACK data
f2 <- project(fire,"EPSG:4326") 
#f3 <- terra::project(rf,fire)
# extract CACK
frf <- terra::extract(rf,f2, "mean", method = "bilinear", bind = T) 

# merge the two vectors, so we have CACK with the original un-transformed vector
fire <- merge(fire, frf)

# read Tree Cover and percent larch 
# tree cover is from Hansen product - Landsat year 2000 (prefire)
# larch is percent larch pixels in fire perimeter based on ESA CCI
# both data sets provided by E. Webb, extracted in GEE
tree <- read.csv("data/fires_treecover_larch.csv", header = T)

# merger tree data with fire database
fire <- merge(fire, tree)

# write updated shapefile
writeVector(fire, "data/SiberiaFires2001-2020/SiberiaFires2001-2020updated_tree.shp", overwrite = T)

# get locations of Retrogressive Thaw Slumps from Runge et al 2022
# not this is not included in analysis

# load rts point locations (https://doi.pangaea.de/10.1594/PANGAEA.941479)
#rts <- vect("data/RTS_NESiberia/RTS_NESiberia.shp")

# reproject rts point locations
#rts <- project(rts,fire)

# use the relate function to find rts within fire perimeters
#j <- relate(rts,fire, "within", pairs = T)
#j <- relate(fire,rts,"contains", pairs = T)

# combine data from each file - not this will be a point vector, which is OK
#rtsf <- cbind(rts[j[,1],1:9],fire[j[,2],1:6])

# calculate when RTS initiated relative to fire
#rtsf$tsf <- rtsf$FirstYear-rtsf$FireYr

# write shapefile of RTS points with fire info
#writeVector(rtsf, "data/RTS_NESiberia/RTS_NESiberia_fire.shp", overwrite = T)
#--------------------------------------------------------------------------#
####################################################################

# read updated fire shapefile with CACK and tree cover data
fire <- vect("data/SiberiaFires2001-2020/SiberiaFires2001-2020updated_tree.shp")
#convert attributes to data frame 
fr <- as.data.frame(fire)
#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#
# READ AND PROCESS MODIS ALBEDO DATA FORE EACH FIRE PERIMETER

# list yearly albedo files
af <- list.files(path = "data/albedo/", pattern = ".csv", full.names = T)

# read first file
alb <- read.csv(af[1], header = T, check.names = F)

# append subsequent files
for(i in 2:length(af))
{
  tmp.alb <- read.csv(af[i], header = T, check.names = F)
  alb <- rbind(alb, tmp.alb)
}

# get date variables
alb$year <- year(strptime(alb$date,format = "%F"))
alb$month <- month(strptime(alb$date,format = "%F"))
alb$day <- day(strptime(alb$date,format = "%F"))

# join variables from fire data set and calculate year since fire
# note using inner join here - there are three fires for which there is no albedo data
alb.fire <- fr %>%
  select(UniqueId, FireYr) %>%
  inner_join(alb,) %>%
  mutate(ysf = year-FireYr)

write.csv(alb.fire, "data/siberia_fire_albedo.csv")

#--------------------------------------------------------------------------#
# I think the code below is better suited for the analyses, 

# calculate prefire albedo
pre <- alb.fire %>%
  select(UniqueId,month,day,ysf,albedo) %>%
  filter(ysf < 0) %>%
  group_by(UniqueId,month,day) %>%
  summarise(pre.alb = mean(albedo, na.rm = T))

# make long CACK data frame for RF calcs
frc <- fr %>%
  select(UniqueId,starts_with("CACK")) %>%
  pivot_longer(cols = starts_with("CACK"),
               names_to = "month",
               names_prefix = "CACK_CM_",
               names_transform = list(month = as.integer),
               values_to = "frc.cack") %>%
  filter(month < 10 & month > 2)

# join all of these into a single data frame?
alb.rf <- alb.fire %>%
  full_join(pre) %>%
  full_join(frc)


