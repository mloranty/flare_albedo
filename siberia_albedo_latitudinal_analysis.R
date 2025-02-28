#####################################
# analyze MODIS albedo for fire perimeters
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

#read table of post-fire radiative forcing data
frc <- read.csv(file = "data/biweekly_radiative_forcing_per_fire.csv", header = T)
d.alb <- read.csv(file = "data/postfire_delta_albedo.csv", header = T)
pre <- read.csv(file = "data/prefire_mean_albedo.csv", header = T)
post <- read.csv(file = "data/postfire_albedo.csv", header = T)

# fix column header/fire ids (R appends X if they start with a number)
z <- ifelse(nchar(colnames(frc))==9,substr(colnames(frc),2,9),colnames(frc))
colnames(frc) <- z

z <- ifelse(nchar(colnames(d.alb))==9,substr(colnames(d.alb),2,9),colnames(d.alb))
colnames(d.alb) <- z

z <- ifelse(nchar(colnames(pre))==9,substr(colnames(pre),2,9),colnames(pre))
colnames(pre) <- z
rm(z)

z <- ifelse(nchar(colnames(post))==9,substr(colnames(post),2,9),colnames(post))
colnames(post) <- z
rm(z)
###################################################################
# create plots
###################################################################
# summarize area burned by latitude
# latitude of fires ranges from 50.2 to 76.2

# convert fire attribute to data frame
fr <- as.data.frame(fire)

# convert fire year to factor
fr$fy <- as.factor(fr$FireYr)

# create var with latitude rounded to nearest degree
fr$lat.r <- round(fr$lat)

# create var indicating 5 year time periods
# f$tper <- ifelse(f$FireYr < 2006,"2001-2005",
#                  ifelse(f$FireYr > 2005 & f$FireYr<2011, "2006-2010", 
#                         ifelse(f$FireYr>2010 & f$FireYr<2016,"2011-2015","20016-2020")))

# an easier alternative using dplyr
fr <- fr %>% mutate(yr.bin = cut(FireYr, breaks = c(2000,2005,2010,2015,2021)))

# create var with latitudinal bins by 5 degrees
fr <- fr %>% mutate(lat.bin5 = cut(lat.r, breaks = c(49,55,60,65,70,77)))

# group by these bins
f.lat5 <- fr %>%
  group_by(yr.bin, lat.bin5) %>%
  summarise(area = sum(SizeHa), severity = mean(dnbr), tree = mean(treecover2))


# calculate number of fires retained for analysis
nf <- fr %>%
  filter(lat.bin5!="(49,55]" & lat!="(70,77]")

# create and print a barplot of area burned by latitude - 5 degree bins
bp5 <- ggplot(f.lat5, aes(fill=yr.bin, y=area/10^6, x=lat.bin5)) + 
  geom_col(position="dodge") +
  xlab("Latitude") + ylab("Area Burned (MHa)") +
  scale_fill_manual(values = c('#fdbe85','#fd8d3c','#e6550d','#a63603'),
                    labels = c("2001-2005","2006-2010","2011-2015","2016-2020")) +
  scale_x_discrete(labels = c("50-55", "55-60", "60-65", "65-70", "70 < ")) +
  labs(fill = "") +
  theme_bw(base_size = 16)

bp5 + theme(#legend.position = c(0.85, 0.75),
            legend.title=element_blank(),
            legend.position = "top")
ggsave("figures/area_burned_by_latitude.png",
       width = 7, height = 4, units = "in")

# boxplot of mean dNBR in 5 degree latitude bins
sv <- ggplot(fr, aes(x=lat.bin5, y=dnbr/100)) + 
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill = '#e6550d') +
 # coord_cartesian(ylim = quantile(f$dnbr/100, c(0.05, 0.95))) +
  coord_cartesian(ylim = c(0, 4)) +
  xlab("Latitude") + ylab("dNBR") + 
  scale_x_discrete(labels= c("50-55", "55-60", "60-65", "65-70", "70 < ")) +
#  scale_fill_manual(values = c('#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494')) +
  theme_bw(base_size = 18)

ggsave("figures/dnbr_by_latitude.png",
       width = 6, height = 6, units = "in")

# boxplot of mean canopy cover in 5 degree latitude bins
cc <- ggplot(fr, aes(x=lat.bin5, y=treecover2)) + 
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill = '#33a02c') +
  # coord_cartesian(ylim = quantile(f$dnbr/100, c(0.05, 0.95))) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Latitude") + ylab("Canopy Cover (%)") + 
  scale_x_discrete(labels= c("50-55", "55-60", "60-65", "65-70", "70 < ")) +
  theme_bw(base_size = 16) 
  

ggsave("figures/canopy_cover_by_latitude.png",
       width = 6, height = 4, units = "in")

# DELETED - monthly mean albedo by ecozone, biome, and region
# messy and not informative



###################################################################
# plot changes in albedo post-fire
###################################################################

# mean albedo by latitude for each year post-fire
ala <- select(alb.lat, c(-month,-day,-jday)) %>%
  group_by(ysf,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T))

# mean albedo by latitude for each month and year year post-fire
alm <- select(alb.lat, c(-day,-jday)) %>%
  group_by(ysf,month,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T))

# look at years, 1, 2, 3-5, and 5-10 post-fire
b1 <- alb.lat %>%
  filter(2<ysf & ysf<6) %>%
  group_by(month, day, jday, lat) %>%
  summarise(ysf = mean(ysf, na.rm = T),
            albedo = mean(albedo, na.rm = T)) 

b2 <- alb.lat %>%
  filter(5<ysf & ysf<11) %>%
  group_by(month, day, jday, lat) %>%
  summarise(ysf = mean(ysf, na.rm = T),
            albedo = mean(albedo, na.rm = T)) 

b3 <- alb.lat %>%
  filter(10<ysf) %>%
  group_by(month, day, jday, lat) %>%
  summarise(ysf = mean(ysf, na.rm = T),
            albedo = mean(albedo, na.rm = T)) 

alb.lat.bin <- alb.lat %>%
  filter(0>ysf | 0<ysf & ysf<3) %>%
  rbind(b1,b2,b3)


rm(b1,b2,b3)

alb.lat.bin$tsf <- as.factor(alb.lat.bin$ysf)

# mplot post-fire albedo by annual bins, with a separate plot for each latitude bin
a <- ggplot() +
  geom_line(data = alb.lat.bin, aes(x = jday, y = albedo, color = tsf), linewidth = 1) + 
  scale_color_manual(values = c("red", hcl.colors(8,palette = "Blues 3",rev=T)[4:8]),
                     labels = c("Prefire", "1", "2", "3-5", "6-10", "11-20")) +
  theme_bw(base_size = 16) 


ggsave("figures/postfire_seasonal_albedo_latitude_freey.png",
       width = 12, height = 8, units = "in")

# make a plot that excludes southern and northern most latitudinal bins
alb.lat.bin.bor <- alb.lat.bin %>%
  filter(lat!="(49,55]" & lat!="(70,77]")

ab <- ggplot() +
  geom_line(data = alb.lat.bin.bor, aes(x = jday, y = albedo, color = tsf), linewidth = 0.75) + 
  scale_color_manual(values = c("red", hcl.colors(8,palette = "Blues 3",rev=T)[4:8]),
                     labels = c("Prefire", "1", "2", "3-5", "6-10", "11-20")) +
  theme_bw(base_size = 18) 

ab + facet_wrap(vars(lat), scales = "fixed",
               labeller = as_labeller(c("(55,60]" = "55-60",
                                        "(60,65]" = "60-65",
                                        "(65,70]" = "65-70"))) +
  xlab("Julian Day") + 
  ylab("Albedo") + 
  labs(color = "Years\nSince\nFire")

ggsave("figures/postfire_seasonal_albedo_latitude_boreal_fix.png",
       width = 12, height = 4, units = "in")


# plot april and july albedo by year since fire and latitude
alm.ba <- alb.lat.bin.bor%>%
  filter(month == 3) %>%
  group_by(ysf,month,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T))

alm.bj <- alb.lat.bin.bor%>%
  filter(month == 7) %>%
  group_by(ysf,month,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T))

alm.a <- alb.lat.bin.bor%>%
  group_by(ysf,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T))

sz <- 1.5
cl <- c("#fde725", "#21918c", "#440154")
# plot for March
ab1 <- ggplot(data = alm.ba, aes(x = ysf, y = albedo, color = lat)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("March") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

# plot for July
ab2 <- ggplot(data = alm.bj, aes(x = ysf, y = albedo, color = lat)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  labs(color = "Latitude") +
  ggtitle("July") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

# plot for mean 
ab3 <- ggplot(data = alm.a, aes(x = ysf, y = albedo, color = lat)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("March-Sept") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

ab1+ theme(legend.position = "none") +ab2 + theme(legend.position = "none") +ab3

ggsave("figures/postfire_delta_albedo_latitude_boreal.png",
       width = 12, height = 4, units = "in")
###################################################################
# plot changes in radiative forcing post-fire
# mean radiative forcing over March-Sept by year since fire,
# split by latitude as well
###################################################################

# mean RF for each year year post-fire
frc.ann <- select(frc,c(-dsf,-month,-day)) %>%
  filter(ysf > 0) %>%
  group_by(ysf) %>%
  summarise(across(1:22088, ~mean(.x, na.rm = T)))

# look at years, 1, 2, 3-5, and 5-10 post-fire
b1 <- select(frc,-dsf) %>%
  filter(2<ysf & ysf<6) %>%
  group_by(month, day) %>%
  summarise(across(1:22089,~mean(.x, na.rm = T))) 

b2 <- select(frc,-dsf) %>%
  filter(5<ysf & ysf<11) %>%
  group_by(month, day) %>%
  summarise(across(1:22089,~mean(.x, na.rm = T))) 

b3 <- select(frc,-dsf) %>%
  filter(10<ysf) %>%
  group_by(month, day) %>%
  summarise(across(1:22089,~mean(.x, na.rm = T))) 

frc.yr <- select(frc,-dsf) %>%
  filter(0<ysf & ysf<3) %>%
  full_join(b1) %>%
  full_join(b2) %>%
  full_join(b3)

rm(b1,b2,b3)
# test using base R
# b3 <- which(frc$ysf>2 & frc$ysf<6)
# b4 <- whichwhich(frc$ysf>5 & frc$ysf<11)
# b5 <- which(frc$ysf>10)
# tst <- aggregate(frc[b3,2:22089], by = list(frc$month[b3],frc$day[b3]), FUN = "mean", na.rm = T)


# unique ecozones
eco <- unique(fire$EcoCode)

# get fire ids for first ecozone
fid <- fire$UniqueId[which(fire$EcoCode == eco[1])]

# create data frame with first ecozone
frc.eco <- as.data.frame(cbind(frc.ann$ysf,
                          rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T)))

frc.yr.eco <- as.data.frame(cbind(frc.yr$ysf, frc.yr$month, frc.yr$day,
                                  rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.yr)))],na.rm = T)))
# add remaining ecozones
for( i in 2:length(eco))
{
  fid <- fire$UniqueId[which(fire$EcoCode == eco[i])]
  frc.eco <- cbind(frc.eco,rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))
  frc.yr.eco <- cbind(frc.yr.eco, rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.yr)))],na.rm = T))
}

# set column names
names(frc.eco) <- c("ysf",eco)
names(frc.yr.eco) <- c("ysf", "month", "day",eco)

#---------------------------------------------------------------#
# CALCULATE FORCING BY 5 DEGREE LATITUDINAL BANDS
#---------------------------------------------------------------#
l5 <- unique(fr$lat.bin5)

# get fire ids for first band
fid <- fire$UniqueId[which(fr$lat.bin5 == l5[1])]

# create data frame with first ecozone
frc.lat <- as.data.frame(cbind(frc.ann$ysf,
                               rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T)))

frc.yr.lat <- as.data.frame(cbind(frc.yr$ysf, frc.yr$month, frc.yr$day,
                                  rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.yr)))],na.rm = T)))
# add remaining ecozones
for( i in 2:length(l5))
{
  fid <- fire$UniqueId[which(fr$lat.bin5 == l5[i])]
  frc.lat <- cbind(frc.lat,rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))
  frc.yr.lat <- cbind(frc.yr.lat, rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.yr)))],na.rm = T))
}

# set column names
names(frc.lat) <- c("ysf",levels(l5)[l5])
names(frc.yr.lat) <- c("ysf", "month", "day",levels(l5)[l5])

#---------------------------------------------------------------#
# create separate data for biome/region     
#---------------------------------------------------------------#
# add arctic and subarctic
fid <- fire$UniqueId[which(fire$ArcSub == "subarctic")]
frc.reg <- as.data.frame(cbind(frc.ann$ysf,rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T)))
frc.yr.reg <- as.data.frame(cbind(frc.yr$ysf,frc.yr$month, frc.yr$day,
                                  rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T)))
fid <- fire$UniqueId[which(fire$ArcSub == "arctic")]
frc.reg <- cbind(frc.reg,rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))
frc.yr.reg <- cbind(frc.yr.reg, rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))
#add biome
fid <- fire$UniqueId[which(fire$biome == "Boreal")]
frc.reg <- cbind(frc.reg,rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))
frc.yr.reg <- cbind(frc.yr.reg, rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))

fid <- fire$UniqueId[which(fire$biome == "Tundra")]
frc.reg <- cbind(frc.reg,rowMeans(frc.ann[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))
frc.yr.reg <- cbind(frc.yr.reg, rowMeans(frc.yr[,na.omit(match(fid,colnames(frc.ann)))],na.rm = T))

# set column names 
names(frc.reg) <- c("ysf","subarctic", "arctic","boreal","tundra")
names(frc.yr.reg) <- c("ysf","month","day","subarctic", "arctic","boreal","tundra")

# julian day for plotting 
# (note could do this once above, and then use it to aggregate data)
frc.yr.eco$jday <- yday(strptime(paste(frc.yr.eco$month,frc.yr.eco$day,"2010", sep = "-"), 
                                 format = "%m-%e-%Y", tz = ""))
frc.yr.reg$jday <- yday(strptime(paste(frc.yr.reg$month,frc.yr.reg$day,"2010", sep = "-"), 
                                 format = "%m-%e-%Y", tz = ""))
frc.yr.lat$jday <- yday(strptime(paste(frc.yr.lat$month,frc.yr.lat$day,"2010", sep = "-"), 
                                 format = "%m-%e-%Y", tz = ""))

# convert year since fire bins to factor for categorical plotting
frc.yr.eco$ysf <- factor(frc.yr.eco$ysf, levels = c("1", "2", "4", "8", "15.5"))
frc.yr.reg$ysf <- factor(frc.yr.reg$ysf, levels = c("1", "2", "4", "8", "15.5"))
frc.yr.lat$ysf <- factor(frc.yr.lat$ysf, levels = c("1", "2", "4", "8", "15.5"))

levels(frc.yr.reg$ysf) <- c("1", "2", "3-5", "6-10", "11-20")
levels(frc.yr.eco$ysf) <- c("1", "2", "3-5", "6-10", "11-20")
levels(frc.yr.lat$ysf) <- c("1", "2", "3-5", "6-10", "11-20")
# pivot data into long format
frea <- pivot_longer(frc.eco,
                    cols = 2:9,
                    names_to = "eco",
                    values_to = "rf")

frey <- pivot_longer(frc.yr.eco,
                     cols = 4:11,
                     names_to = "eco",
                     values_to = "rf")
# pivot data into long format
frrg <- pivot_longer(frc.reg,
                     cols = 2:5,
                     names_to = "reg",
                     values_to = "rf")

frrgy <- pivot_longer(frc.yr.reg,
                     cols = 4:7,
                     names_to = "reg",
                     values_to = "rf")
# pivot data into long format
frl <- pivot_longer(frc.lat,
                     cols = 2:6,
                     names_to = "lat",
                     values_to = "rf")

frly <- pivot_longer(frc.yr.lat,
                      cols = 4:8,
                      names_to = "lat",
                      values_to = "rf")
#---------------------------------------------------------------#
colMeans(frc.lat)

# plots for mean annual RF 
# create a plot for regions
p1 <- ggplot(data = frrg, aes(x = ysf, y = rf)) + 
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgray", size = 1) +
  geom_line((aes(color = reg)), na.rm = T) + 
  labs(color = "Region") +
  xlab("Year Since Fire") +
  ylab("Radiative Forcing") +
  #scale_color_brewer("Dark2") +
  theme_bw() 


# create a plot for Ecozones
p2 <- ggplot(data = frea, aes(x = ysf, y = rf)) + 
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgray", size = 1) +
  geom_line((aes(color = eco)), na.rm = T) + 
  labs(color = "Ecozone") +
  xlab("Year Since Fire") +
  ylab("Rediative Forcing") +
  #scale_color_brewer("Dark2") +
  theme_bw() 

# create a plot for latitude
p3 <- ggplot(data = frl, aes(x = ysf, y = rf)) + 
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgray", size = 1) +
  geom_line((aes(color = lat)), na.rm = T) + 
  labs(color = "Latitude") +
  xlab("Year Since Fire") +
  ylab("Rediative Forcing") +
  #scale_color_brewer("Dark2") +
  theme_bw() 
#add plots to layout (patchwork package)
p1+p2

# save plot to file
ggsave("figures/annual_rf_eco_reg.png",
       width = 8, height = 4, units = "in")

# plots for mean seasonal RF trajectories, binned by 1, 2, 3-5, 6-10, and 11-20 years post fire
# create a multipanel plot for regions
p <- ggplot(data = frrgy, aes(x = jday, y = rf)) + 
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgray", size = 1) +
  geom_point((aes(color = ysf)), na.rm = T) + 
  scale_color_brewer(palette="Set1") +
  labs(color = "Years Since Fire") +
  #scale_color_brewer("Dark2") +
  theme_bw() 

#p + scale_fill_discrete(limits = c("1", "2", "4", "8", "15.5"))
p + facet_wrap(vars(reg), scales = "free_y") +
  xlab("Day of Year") + ylab("Radiative Forcing")

ggsave("figures/seasonal_rf_reg.png",
       width = 8, height = 4, units = "in")

# create a multipanel plot for ecozones
p <- ggplot(data = frey, aes(x = jday, y = rf)) + 
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgray", size = 1) +
  geom_point((aes(color = ysf)), na.rm = T) + 
#  scale_color_brewer(palette="Set1") +
  labs(color = "Years Since Fire") +
  #scale_color_brewer("Dark2") +
  theme_bw() 

#p + scale_fill_discrete(limits = c("1", "2", "4", "8", "15.5"))
p + facet_wrap(vars(eco), scales = "free_y") +
    xlab("Day of Year") + 
    ylab(expression(paste("Radiative Forcing W",m^-2,sep="")))

ggsave("figures/seasonal_rf_eco.png",
       width = 8, height = 4, units = "in")

# create a multipanel plot for latitudes
p <- ggplot(data = frly, aes(x = jday, y = rf)) + 
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "darkgray", size = 1) +
  geom_point((aes(color = ysf)), na.rm = T) + 
  #  scale_color_brewer(palette="Set1") +
  labs(color = "Years Since Fire") +
  #scale_color_brewer("Dark2") +
  theme_bw() 

#p + scale_fill_discrete(limits = c("1", "2", "4", "8", "15.5"))
p + facet_wrap(vars(lat), scales = "free_y") +
  xlab("Day of Year") + 
  ylab(expression(paste("Radiative Forcing W",m^-2,sep="")))

#--------------------------------------------------#
# boxplot of annual rf by lat bins and time since fire
#--------------------------------------------------#

# maybe need to do this instead with the non-aggregated data set and aggregate to bins of ysf (1,2,3-5,6-10, etc...)
# exclude upper and lower lat bands
frl.bor <- frly %>%
  filter(lat!="(49,55]" & lat!="(70,77]")

bp <- ggplot(frl.bor, aes(x = ysf, y=rf, fill = lat)) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_boxplot(position ="dodge2") +
   #     ylim(-7,7) +
        xlab("Years Since Fire") +
        ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
        ggtitle("March-Sept Radiative Forcing") +
        labs(fill = "Latitude") +
        scale_fill_manual(values = cl,
                          labels = c("55-60", "60-65", "65-70")) +
        theme_bw(base_size = 16)


ggsave("figures/annual_rf_tsf_latitude_boreal.png",
       width = 6, height = 4, units = "in")

#--------------------------------------------------#
# boxplots of spring and summer rf by lat bins and time since fire
#--------------------------------------------------#
# July
bp <- ggplot(frly[which(frly$month==7),], aes(fill=lat, x = ysf, y=rf)) +
  geom_boxplot(position ="dodge2") +
  ylim(-1,3) +
  xlab("Time Since Fire (Years)") +
  ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
  labs(color = "Latitude") +
  scale_color_manual(values = c("red", "orange", "black", "blue", "purple"),
                     labels = c("50-55", "56-60", "61-65", "66-70", "70 <")) 

ggsave("figures/july_rf_tsf_latitude.png",
       width = 6, height = 4, units = "in")

#April
bp <- ggplot(frly[which(frly$month==4),], aes(fill=lat, x = ysf, y=rf)) +
  geom_boxplot(position ="dodge2") +
 # ylim(-1,3) +
  xlab("Time Since Fire (Years)") +
  ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
  labs(color = "Latitude") +
  scale_color_manual(values = c("red", "orange", "black", "blue", "purple"),
                     labels = c("50-55", "56-60", "61-65", "66-70", "70 <"))
  
  ggsave("figures/april_rf_tsf_latitude.png",
         width = 6, height = 4, units = "in")
  
# March-May
bp <- ggplot(frly[which(frly$month<6),], aes(fill=lat, x = ysf, y=rf)) +
  geom_boxplot(position ="dodge2") +
  ylim(-5,5) +
  xlab("Time Since Fire (Years)") +
  ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
  labs(color = "Latitude") +
  scale_color_manual(values = c("red", "orange", "black", "blue", "purple"),
                     labels = c("50-55", "56-60", "61-65", "66-70", "70 <"))

ggsave("figures/spring_rf_tsf_latitude.png",
       width = 6, height = 4, units = "in")