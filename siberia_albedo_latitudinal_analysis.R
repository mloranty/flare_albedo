#####################################
# analyze MODIS albedo for fire perimeters
# from Talucci et al 2022
# albedo provided from Google Earth Engine
# by E. Webb
#
# see albedo_fire_preprocess.R for initial 
# data wrangling
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

#read fire polygons with tree cover and CACK attributes
fire <- vect("data/SiberiaFires2001-2020/SiberiaFires2001-2020updated_tree.shp")

#read tables of albedo and post-fire radiative forcing data
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
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  nrow()

# filter the highest and lowest latitude bands
f.lat5 <- f.lat5 %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") 
  
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
cc <- ggplot(fr, aes(x=lat.bin5, y=treecover2)) + #larch_perc
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill = '#33a02c') +
  # coord_cartesian(ylim = quantile(f$dnbr/100, c(0.05, 0.95))) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Latitude") + ylab("Canopy Cover (%)") + 
  scale_x_discrete(labels= c("50-55", "55-60", "60-65", "65-70", "70 < ")) +
  theme_bw(base_size = 16) 
  
ggsave("figures/canopy_cover_by_latitude.png",
       width = 6, height = 4, units = "in")


###################################################################
# MAKE ALBEDO DATA FRAMES LONGER AND COMBINE WITH FIRE INFO
###################################################################
 pre.alb <- pivot_longer(pre,
                     cols = c(3:22090),
                     names_to = "UniqueId", 
                     values_to = "albedo") 

pre.m <- pre.alb %>%
  group_by(UniqueId, month) %>%
  summarise(albedo = mean(albedo, ne.rm = T)) %>%
  full_join(select(fr,c("UniqueId","dnbr","lat.bin5","treecover2")), by = "UniqueId")


pre.l <- fr %>%
  select(UniqueId,biome,lat, lon, dnbr, yr.bin, lat.bin5) %>%
  full_join(pre.alb, by = "UniqueId")

###################################################################
# delta albedo #
delta <- pivot_longer(d.alb,
                      cols = c(2:22089),
                      names_to = "UniqueId", 
                      values_to = "d.albedo") 

delta.m <- delta %>%
  group_by(UniqueId,ysf,month) %>%
  summarise(d.albedo = mean(d.albedo, na.rm = T)) %>%
  full_join(select(fr,c("UniqueId","dnbr","lat.bin5","treecover2")), by = "UniqueId")

# delta.l <- fr %>%
#   select(UniqueId,biome,lat, lon, dnbr, yr.bin, lat.bin5) %>%
#   full_join(delta, by = "UniqueId")

#plot postfire delta albedo for march, by latitude band
sz <- 1.5
cl <- c("#fde725", "#21918c", "#440154")

delt.mar.plot <- delta.m %>%
  filter(ysf > 0) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(d.albedo = mean(d.albedo, na.rm = T),
            stdev = sd(d.albedo, na.rm = T)) %>%
  filter(month == 3) %>%
  ggplot(aes(x = ysf, y = d.albedo, color = lat.bin5)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab(expression(paste(Delta,"albedo"))) + 
  ggtitle("March") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

#plot postfire delta albedo for july, by latitude band
delt.jul.plot <- delta.m %>%
  filter(ysf > 0) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(d.albedo = mean(d.albedo, na.rm = T),
            stdev = sd(d.albedo, na.rm = T)) %>%
  filter(month == 7) %>%
  ggplot(aes(x = ysf, y = d.albedo, color = lat.bin5)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab(expression(paste(Delta,"albedo"))) + 
  ggtitle("July") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

#plot postfire delta albedo for july, by latitude band
delt.ann.plot <- delta.m %>%
  filter(ysf > 0) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  group_by(lat.bin5, ysf) %>%
  summarise(d.albedo = mean(d.albedo, na.rm = T),
            stdev = sd(d.albedo, na.rm = T)) %>%
  ggplot(aes(x = ysf, y = d.albedo, color = lat.bin5)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab(expression(paste(Delta,"albedo"))) + 
  ggtitle("Mar-Sept") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

delt.mar.plot+ theme(legend.position = "none") +delt.jul.plot + theme(legend.position = "none") +delt.ann.plot

ggsave("figures/postfire_delta_albedo_annual_latitude_boreal.png",
       width = 12, height = 4, units = "in")
###################################################################
# ALBEDO BY 5 DEGREE LATITUDE BINS
###################################################################
# get unique lat bins from fire data set
l5 <- unique(fr$lat.bin5)

# get fire ids for first Lat bin
fid <- fr$UniqueId[which(fr$lat.bin5 == l5[1])]

# create data frame with first Lat bin
pre.lat <- as.data.frame(cbind(pre$month, pre$day,
                               rowMeans(pre[,na.omit(match(fid,colnames(pre)))],na.rm = T)))

post.lat <- as.data.frame(cbind(post$ysf, post$month, post$day,
                                rowMeans(post[,na.omit(match(fid,colnames(post)))],na.rm = T)))
# add remaining bins
for( i in 2:length(l5))
{
  fid <- fr$UniqueId[which(fr$lat.bin5 == l5[i])]
  pre.lat <- cbind(pre.lat,rowMeans(pre[,na.omit(match(fid,colnames(pre)))],na.rm = T))
  post.lat <- cbind(post.lat,rowMeans(post[,na.omit(match(fid,colnames(post)))],na.rm = T))
}

# set column names
names(pre.lat) <- c("month", "day",levels(l5)[l5])
names(post.lat) <- c("ysf","month", "day",levels(l5)[l5])

pre.lat$ysf <- -1

# add julian day for plotting (assume non-leap year)
pre.lat$jday <- yday(strptime(paste(pre.lat$month,pre.lat$day,"2010", sep = "-"), 
                              format = "%m-%e-%Y", tz = ""))

post.lat$jday <- yday(strptime(paste(post.lat$month,post.lat$day,"2010", sep = "-"), 
                               format = "%m-%e-%Y", tz = ""))

# combine to a single dataframe for plotting
alb.l <- rbind(pre.lat[,c(9,1,2,8,6,4,3,5,7)],
               post.lat[,c(1:3,9,7,5,4,6,8)])

alb.lat <- pivot_longer(alb.l,
                        cols = 5:9,
                        names_to = "lat",
                        values_to = "albedo")

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


# plot march and july albedo by year since fire and latitude
alm.ba <- alb.lat.bin.bor%>%
  filter(month == 3) %>%
  group_by(ysf,month,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T),
            sd = sd(albedo, na.rm = T))

alm.bj <- alb.lat.bin.bor%>%
  filter(month == 7) %>%
  group_by(ysf,month,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T),
            sd = sd(albedo, na.rm = T))

alm.a <- alb.lat.bin.bor%>%
  group_by(ysf,lat) %>%
  summarise(albedo = mean(albedo, na.rm = T),
            sd = sd(albedo, na.rm = T))

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