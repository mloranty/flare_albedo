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
library(ggplot2)
library(tidyverse)
library(patchwork)
#library(ncdf4)

# set working directory with options for Mac or PC workstations
ifelse(Sys.info()[1]=="Windows",
       setwd("G:/My Drive/Documents/research/manuscripts/siberia_albedo"),
       setwd("~/Library/CloudStorage/GoogleDrive-mloranty@colgate.edu/My Drive/Documents/research/manuscripts/siberia_albedo"))

# RAW DATA - THIS HAS BEEN PRE-PROCESSED, SEE COLLAPSED CODE BELOW
# read albedo data
#alb <- read.table("data/siberia_fire_albedo.csv", header = T, check.names = F)
# read fire perimeter data set
#fire <- vect("data/SiberiaFires2001-2020/SiberiaFires2001-2020.shp")


# create time since fire date stamp for each fire 
# create new time-since-fire albedo dataframe
# expand to see code, otherwise skip and read alb.tsf file
###########################################################################
# create date/time object
alb$tmstmp <- strptime(alb$date,format = "%F")

# match albedo to fire record
#ref <- match(colnames(alb.s),fire$UniqueId)
ref <- match(colnames(alb),fire$UniqueId)

# which columns in the data frame have matches in fire database
r <- which(is.na(ref)==F)

#fire year
fy <- fire$FireYr[ref[r[1]]]

# year since fire
#ysf <- year(alb.s$tmstmp)-fy
ysf <- year(alb$tmstmp)-fy

# date since fire
#alb.tsf <- data.frame(paste(ysf,month(alb.s$tmstmp), day(alb.s$tmstmp),sep = "/"))
alb.tsf <- data.frame(paste(ysf,month(alb$tmstmp), day(alb$tmstmp),sep = "/"))

# add albedo data
alb.tsf <- cbind(alb.tsf,alb[,r[1]])

#set names
names(alb.tsf) <- c("dsf",names(alb)[r[1]])

#iterate through the remaining columns and merge to a single data frame
for(i in 2:length(r))
{
  fy <- fire$FireYr[ref[r[i]]]
  ysf <- year(alb$tmstmp)-fy
  t <- data.frame(paste(ysf,month(alb$tmstmp), day(alb$tmstmp),sep = "/"))
  t <- cbind(t,alb[,r[i]])
  names(t) <- c("dsf",names(alb)[r[i]])
  alb.tsf <- merge(alb.tsf,t, by = "dsf", all = TRUE)
}

# create separate numeric variables for time since fire date info
t <- strsplit(alb.tsf$dsf, split = "/")

#year since fire
alb.tsf$ysf <- as.numeric(sapply(t,"[[",1))
#month
alb.tsf$month <- as.numeric(sapply(t,"[[",2))
#day 
alb.tsf$day <- as.numeric(sapply(t,"[[",3))

# write the output to a csv file
write.csv(alb.tsf,file = "data/siberia_fire_albedo_tsf.csv", row.names = F, col.names = T)
###########################################################################

#### read data frame for albedo as time since fire ##
alb.tsf <- read.csv("data/siberia_fire_albedo_tsf.csv", header = T)


# add new variables to fire database for analysis (biome, etc.)
# combine retrogressive thaw slump and fire data sets
# expand to see code, otherwise skip and read updated shapefile
###########################################################################

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
writeVector(fire, "data/SiberiaFires2001-2020/SiberiaFires2001-2020updated.shp", overwrite = T)
writeVector(fire, "data/SiberiaFires2001-2020/SiberiaFires2001-2020updated_tree.shp", overwrite = T)

# get locations of Retrogressive Thaw Slumps from Runge et al 2022
# not this is not included in analysis
####################################################################
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
###########################################################################

# read updated fire shapefile this one uses old CACK integration and adds tree data 
# need to do some cleanup above
fire <- vect("data/SiberiaFires2001-2020/SiberiaFires2001-2020updated_tree.shp")

# read RTS shapefiles if analysis will be included...
#rts <- vect("data/RTS_NESiberia/RTS_NESiberia.shp")
#rts.f <- vect("data/RTS_NESiberia/RTS_NESiberia_fire.shp")

## calculate mean pre-fire albedo for the 20 years pre-fire
## and then difference between pre- and post, for RF calculations
## expand to see code, otherwise skip and read relevant csv files
###################################################################
## use dplyr to avoid getting (too) rusty
pre <- alb.tsf %>%
  filter(ysf < 0) %>%
  group_by(month, day) %>%
  summarise(across(2:22089,mean, na.rm = T))
write.csv(pre, file = "data/prefire_mean_albedo.csv", row.names = F)

# filter to create data frame with only post-fire albedo
post <- alb.tsf %>%
  filter(ysf >= 0) %>%
  arrange(ysf,month,day)
write.csv(post, file = "data/postfire_albedo.csv", row.names = F)

# create a reference vector matching rows in post-fire to rows in pre-fire by month and day
r <- match(paste(post$month, post$day, sep = "."),paste(pre$month, pre$day, sep = "."))

# subtract post-fire albedo from pre-fire to get delta
d.alb <- post
d.alb[,2:22089] <- pre[r,3:22090]-post[,2:22089]

# write data fo file
write.csv(d.alb, file = "data/postfire_delta_albedo.csv", row.names = F)
## CALCULATE RADIATIVE FORCING USING CACK DATA
# column names in albedo data frames correspond to unique fire IDs
# but those beginning with a number have an X appended when they are read into R
# need to remove this
z <- ifelse(nchar(colnames(d.alb))==9,substr(colnames(d.alb),2,9),)

# now match fire IDs from shapefile with columns in the delta albedo data frame
r <- match(z[2:22089],fire$UniqueId)

# make a data from from the CACK kernel for March - Sept
y <- data.frame(rbind(fire$CACK_CM_3, fire$CACK_CM_4,fire$CACK_CM_5,
                      fire$CACK_CM_6, fire$CACK_CM_7,fire$CACK_CM_8,fire$CACK_CM_9))

# set colnames as UniqueId
names(y) <- fire$UniqueId

# reorder columns to match albedo data frames
cac <- y[,r]
rm(y)

# add month variable to cack
cac$month <- 3:9

# match month records in delta albedo and cac data frames
r <- match(d.alb$month, cac$month)

# multiply cack x delta albedo to calculate RF values
frc <- d.alb
frc[,2:22089] <- cac[r,1:22088]*d.alb[,2:22089]
write.csv(frc, file = "data/biweekly_radiative_forcing_per_fire.csv", row.names = F)
###################################################################