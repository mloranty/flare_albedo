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

#--------------------------------------------------------------------------#
# read fire perimeters and MODIS albedo data
# see albedo_fire_preprocess.R for details
#--------------------------------------------------------------------------#

# read fire polygons with tree cover and CACK attributes
fire <- vect("data/SiberiaFires2001-2020/SiberiaFires2001-2020updated_tree.shp")

# read tables of albedo data
alb <- read.csv(file = "data/siberia_fire_albedo.csv", header = T)


#--------------------------------------------------------------------------#
# summarize latitudinal tree cover and burned area patterns
#--------------------------------------------------------------------------#
############################################################################
# summarize area burned by latitude
# latitude of fires ranges from 50.2 to 76.2
cl <- c("#fde725", "#21918c", "#440154")
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
  summarise(area = sum(SizeHa), severity = mean(dnbr), tree = mean(treecover2)) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") 

# calculate number of fires retained for analysis
nf <- fr %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
#  filter(biome=="Tundra") %>%
  nrow()
  
# create and print a barplot of area burned by latitude - 5 degree bins
bp5 <- f.lat5 %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  ggplot(aes(fill=yr.bin, y=area/10^6, x=lat.bin5)) + 
  geom_col(position="dodge") +
  xlab("Latitude") + ylab("Area Burned (MHa)") +
  scale_fill_manual(values = c('#fdbe85','#fd8d3c','#e6550d','#a63603'),
                    labels = c("2001-2005","2006-2010","2011-2015","2016-2020")) +
  scale_x_discrete(labels= c("55-60", "60-65", "65-70")) +
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
cc <- fr %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  ggplot(aes(x=lat.bin5, y=treecover2)) + #larch_perc
  geom_boxplot(notch = TRUE, outlier.shape = NA, fill = cl) +
  # coord_cartesian(ylim = quantile(f$dnbr/100, c(0.05, 0.95))) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab("Latitude") + ylab("Canopy Cover (%)") + 
  scale_x_discrete(labels= c("55-60", "60-65", "65-70")) +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 16) 
  
ggsave("figures/canopy_cover_by_latitude.png",
       width = 6, height = 4, units = "in")
###################################################################
#-------------------------------------------------------------------------------------------------------------------------#
# prepare data for plotting
#-------------------------------------------------------------------------------------------------------------------------#

# calculate pre-fire albedo
pre <- alb %>%
  select(UniqueId,month,day,ysf,albedo) %>%
  filter(ysf < 0) %>%
  group_by(UniqueId,month,day) %>%
  summarise(pre.alb = mean(albedo, na.rm = T)) 

pre.m <- pre %>%
  group_by(UniqueId, month) %>%
  summarise(albedo = mean(pre.alb, na.rm = T)) %>%
  mutate(ysf = -1) %>%
  full_join(select(fr,c("UniqueId","dnbr","lat.bin5","treecover2")), by = "UniqueId") %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]")

# get cack kernal values for rf calcs
frc <- fr %>%
  select(UniqueId,starts_with("CACK")) %>%
  pivot_longer(cols = starts_with("CACK"),
               names_to = "month",
               names_prefix = "CACK_CM_",
               names_transform = list(month = as.integer),
               values_to = "frc.cack") %>%
  filter(month < 10 & month > 2)

# join these to have a large dataframe
alb.all <- alb %>%
  select(-X) %>% 
  full_join(pre) %>%
  full_join(select(fr,c("UniqueId","dnbr","lat.bin5","treecover2")), by = "UniqueId") %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]")

# add post-fire bins
alb.all$pf.bin <- as.factor(case_match(alb.all$ysf,
                       c(-20:-1)~"pre",
                       0~"FY",
                       c(1:2)~"1-2",
                       c(3:5)~"3-5",
                       c(6:23)~"6-23"))

alb.all$pf.bin <- fct_relevel(alb.all$pf.bin,c("pre", "1-2", "3-5", "6-23", "FY"))
 #levels(alb.all$pf.bin) <- c("pre", "FY", "1-2", "3-5", "6-23")
# add julian day vars
alb.all$jday <- yday(strptime(paste(alb.all$month,alb.all$day,"2010", sep = "-"), 
              format = "%m-%e-%Y", tz = ""))
pre$jday <- yday(strptime(paste(pre$month,pre$day,"2010", sep = "-"), 
                          format = "%m-%e-%Y", tz = ""))

#----------------#
# monthly albedo #
#----------------#

alb.m <- select(alb.all, c("UniqueId", "albedo", "pre.alb", "month", "day", "ysf")) %>%
  filter(ysf > 0) %>%
  mutate(d.alb = pre.alb-albedo) %>%
  group_by(UniqueId, ysf, month) %>%
  summarise(albedo = mean(albedo, na.rm = T),
            pre.alb = mean(pre.alb, na.rm = T),
            d.alb = mean(d.alb, na.rm = T)) %>%
  full_join(select(fr,c("UniqueId","dnbr","lat.bin5","treecover2")), by = "UniqueId")

#-------------------------------------------------------------------------------------------------------------------------#
# calculate summaries for pre-fire tree cover and albedo and perform ANOVAs
#-------------------------------------------------------------------------------------------------------------------------#

# tree cover summary stats by latitude
cc.tbl <- fr %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  group_by(lat.bin5) %>%
  summarise(tree = mean(treecover2, na.rm = T),
            tree.sd = sd(treecover2, na.rm = T),
            tree.sem = sd(treecover2, na.rm = T)/n()) %>%
  write.csv(,file = "results/tree_cover_summary.csv")

# tree cover by latitude ANOVA
cc3 <- fr %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") 

  cc.aov <- aov(cc3$treecover2~cc3$lat.bin5)
  TukeyHSD(cc.aov)
  rm(cc3)

# monthly albedo by latitude    
pre3 <- pre.m %>%
    group_by(lat.bin5, month) %>%
    summarise(alb = mean(albedo/1000, na.rm = T),
              alb.sd = sd(albedo/1000, na.rm = T),
              alb.sem = sd(albedo/1000, na.rm = T)/n()) %>%
    write.csv(,file = "results/prefire_albedo_summary.csv")

pre.3 <- pre.m %>%
  filter(month == 3) 

pre3.aov <- aov(pre.3$albedo/1000~pre.3$lat.bin5)
summary(pre3.aov)
TukeyHSD(pre3.aov)

pre.7 <- pre.m %>%
  filter(month == 7) 

pre7.aov <- aov(pre.7$albedo/1000~pre.7$lat.bin5)
summary(pre7.aov)
TukeyHSD(pre7.aov)
#-------------------------------------------------------------------------------------------------------------------------#
# plot of pre-fire albedo for each latitudinal bin
#-------------------------------------------------------------------------------------------------------------------------#
sz <- 1.5
cl <- c("#fde725", "#21918c", "#440154")

pre.alb.lat <- pre %>%
  full_join(select(fr,c("UniqueId","dnbr","lat.bin5","treecover2")), by = "UniqueId") %>%
  group_by(jday, lat.bin5) %>%
  summarise(alb = mean(pre.alb, na.rm = T),
            stdev = sd(pre.alb, na.rm = T)) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  ggplot(aes(x = jday, y = alb/1000, color = lat.bin5)) + 
  geom_point(size = sz) +
  geom_line(size = sz) +
  #geom_ribbon(aes(ymax = (alb+stdev)/1000, ymin = (alb-stdev)/1000), alpha = 0.1, fill = "gray") 
  xlab("Day of Year") + 
  ylab("Albedo") + 
  ggtitle("Pre-fire Albedo") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

ggsave("figures/prefire_albedo_by_latitude.png",
       width = 6, height = 6, units = "in")

cc + labs(tag = "a") + pre.alb.lat + labs(title = "", tag = "b") 

ggsave("figures/prefire_treecover_albedo.png",
       width = 8, height = 4, units = "in")
#-------------------------------------------------------------------------------------------------------------------------#
# plot postfire delta albedo by latitude band  for march, july, and march-sept
#-------------------------------------------------------------------------------------------------------------------------#

sz <- 1.5
cl <- c("#fde725", "#21918c", "#440154")

delt.mar.plot <- alb.m %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(d.alb = mean(d.alb, na.rm = T),
            stdev = sd(d.alb, na.rm = T)) %>%
  filter(month == 3) %>%
  ggplot(aes(x = ysf, y = d.alb/1000, color = lat.bin5)) + 
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
delt.jul.plot <- alb.m %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(d.alb = mean(d.alb, na.rm = T),
            stdev = sd(d.alb, na.rm = T)) %>%
  filter(month == 7) %>%
  ggplot(aes(x = ysf, y = d.alb/1000, color = lat.bin5)) + 
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
delt.ann.plot <- alb.m %>%
  group_by(lat.bin5, ysf) %>%
  summarise(d.albedo = mean(d.alb, na.rm = T),
            stdev= sd(d.alb, na.rm = T)) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  ggplot(aes(x = ysf, y = d.albedo/1000, color = lat.bin5)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  #geom_ribbon(aes(ymax = stdev, ymin = -stdev), alpha = 0.1, fill = "red") 
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

#-------------------------------------------------------------------------------------------------------------------------#
# plot postfire albedo by latitude band for march, july, and march-sept
#-------------------------------------------------------------------------------------------------------------------------#

sz <- 1.5
cl <- c("#fde725", "#21918c", "#440154")

alb.mar.plot <- alb.m %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(alb = mean(albedo, na.rm = T),
            stdev = sd(albedo, na.rm = T)) %>%
  filter(month == 3) %>%
  ggplot(aes(x = ysf, y = alb/1000, color = lat.bin5)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("March") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

#plot postfire delta albedo for july, by latitude band
alb.jul.plot <- alb.m %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(alb = mean(albedo, na.rm = T),
            stdev = sd(albedo, na.rm = T)) %>%
  filter(month == 7) %>%
  ggplot(aes(x = ysf, y = alb/1000, color = lat.bin5)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("July") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

#plot postfire delta albedo for july, by latitude band
alb.ann.plot <- alb.m %>%
  group_by(lat.bin5, ysf) %>%
  summarise(alb = mean(albedo, na.rm = T),
            stdev = sd(albedo, na.rm = T)) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  ggplot(aes(x = ysf, y = alb/1000, color = lat.bin5)) + 
  geom_point(size = sz*2) +
  geom_line(size = sz) +
  #geom_ribbon(aes(ymax = (alb+stdev)/1000, ymin = (alb-stdev)/1000), alpha = 0.1, fill = "red") 
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("Mar-Sept") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

alb.mar.plot+ theme(legend.position = "none") +alb.jul.plot + theme(legend.position = "none") + alb.ann.plot

ggsave("figures/postfire_albedo_annual_latitude_boreal.png",
       width = 12, height = 4, units = "in")

#-------------------------------------------------------------------------------------------------------------------------#
# plot pre and postfire albedo by latitude band for march, july, and march-sept
#-------------------------------------------------------------------------------------------------------------------------#
sz = 1

p.alb.mar.plot <- alb.m %>%
  bind_rows(pre.m) %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(alb = mean(albedo, na.rm = T),
            stdev = sd(albedo, na.rm = T)) %>%
  filter(month == 3) %>%
  ggplot(aes(x = ysf, y = alb/1000, color = lat.bin5)) + 
  geom_point(size = sz) +
  geom_line(size = sz) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("March") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

#plot postfire delta albedo for july, by latitude band
p.alb.jul.plot <- alb.m %>%
  bind_rows(pre.m) %>%
  group_by(lat.bin5, ysf, month) %>%
  summarise(alb = mean(albedo, na.rm = T),
            stdev = sd(albedo, na.rm = T)) %>%
  filter(month == 7) %>%
  ggplot(aes(x = ysf, y = alb/1000, color = lat.bin5)) + 
  geom_point(size = sz) +
  geom_line(size = sz) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("July") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

#plot postfire delta albedo for july, by latitude band
p.alb.ann.plot <- alb.m %>%
  bind_rows(pre.m) %>%
  group_by(lat.bin5, ysf) %>%
  summarise(alb = mean(albedo, na.rm = T),
            stdev = sd(albedo, na.rm = T)) %>%
  filter(lat.bin5!="(49,55]" & lat.bin5!="(70,77]") %>%
  ggplot(aes(x = ysf, y = alb/1000, color = lat.bin5)) + 
  geom_point(size = sz) +
  geom_line(size = sz) +
  geom_vline(xintercept = 0, linetype = 2) +
  xlab("Years Since Fire") + 
  ylab("Albedo") + 
  ggtitle("Mar-Sept") +
  labs(color = "Latitude") +
  scale_color_manual(values = cl,
                     labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 18)

p.alb.mar.plot+ theme(legend.position = "none") + p.alb.jul.plot + theme(legend.position = "none") + p.alb.ann.plot

ggsave("figures/pre_postfire_albedo_annual_latitude_boreal.png",
       width = 12, height = 4, units = "in")

#-------------------------------------------------------------------------------------------------------------------------#
# seasonal albedo trajectories for binned time since fire, for each lat band
#-------------------------------------------------------------------------------------------------------------------------#

alb.pf.bin <- alb.all %>%
  filter(pf.bin != "FY") %>%
  group_by(lat.bin5, pf.bin, jday) %>%
  summarise(alb = mean(albedo/1000, na.rm = T))

alb.pfp.low <- alb.pf.bin %>%
  filter(lat.bin5 == "(55,60]") %>%
  ggplot(aes(x = jday, y = alb, color = pf.bin)) + 
  #geom_point(size = sz*2) +
  geom_line(size = sz) +
  #geom_vline(xintercept = 0, linetype = 2) +
  xlab("Day of Year") + 
  ylab("Albedo") + 
  ggtitle("55-60") +
  labs(color = "", tag = "a") +
  scale_color_manual(values = c("red", hcl.colors(8,palette = "Blues 3",rev=T)[6:8]),
                     labels = c("Prefire", "1-2", "3-5", "6-23")) +
  theme_bw(base_size = 18)

alb.pfp.mid <- alb.pf.bin %>%
  filter(lat.bin5 == "(60,65]") %>%
  ggplot(aes(x = jday, y = alb, color = pf.bin)) + 
  #geom_point(size = sz*2) +
  geom_line(size = sz) +
  #geom_vline(xintercept = 0, linetype = 2) +
  xlab("Day of Year") + 
  ylab("Albedo") + 
  ggtitle("60-65") +
  labs(color = "", tag = "b") +
  scale_color_manual(values = c("red", hcl.colors(8,palette = "Blues 3",rev=T)[6:8]),
                     labels = c("Prefire", "1-2", "3-5", "6-23")) +
  theme_bw(base_size = 18)

alb.pfp.hi <- alb.pf.bin %>%
  filter(lat.bin5 == "(65,70]") %>%
  ggplot(aes(x = jday, y = alb, color = pf.bin)) + 
  #geom_point(size = sz*2) +
  geom_line(size = sz) +
  #geom_vline(xintercept = 0, linetype = 2) +
  xlab("Day of Year") + 
  ylab("Albedo") + 
  ggtitle("65-70") +
  labs(color = "", tag = "c") +
  scale_color_manual(values = c("red", hcl.colors(8,palette = "Blues 3",rev=T)[6:8]),
                     labels = c("Prefire", "1-2", "3-5", "6-23")) +
  theme_bw(base_size = 18)

alb.pfp.low + theme(legend.position = "none") + 
alb.pfp.mid + theme(legend.position = "none") + 
alb.pfp.hi +
p.alb.mar.plot + theme(legend.position = "none") + labs(tag = "d") + 
p.alb.jul.plot + theme(legend.position = "none") + labs(tag = "e") + 
p.alb.ann.plot + labs(tag = "f") 

ggsave("figures/postfire_albedo_latitude_season.png",
       width = 12, height = 8, units = "in")

#-------------------------------------------------------------------------------------------------------------------------#
# radiative forcing plots - spring, summer, and seasonal, by latitude, with ANOVA
#-------------------------------------------------------------------------------------------------------------------------#
# make a plot that excludes southern and northern most latitudinal bins

alb.frc <- alb.m %>%
  full_join(frc, join_by(UniqueId, month)) %>%
  mutate(rf = (d.alb/1000)*frc.cack,
         pf.bin = as.factor(case_match(ysf,
                                       c(1:2)~"1-2",
                                       c(3:5)~"3-5",
                                       c(6:23)~"6-23")))


#boxplot for march
frc.mar <- alb.frc %>%
  na.omit() %>%
  filter(month == 3) %>%
  ggplot(aes(x = pf.bin, y=rf, fill = lat.bin5)) +
  geom_boxplot(position ="dodge2") +
  geom_hline(yintercept = 0, linetype = 2) +
  ylim(-7,7) +
  xlab("Years Since Fire") +
  ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
  ggtitle("March") +
  labs(fill = "Latitude") +
  scale_fill_manual(values = cl,
                    labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 16)

#boxplot for july
frc.jul <- alb.frc %>%
  na.omit() %>%
  filter(month == 7) %>%
  ggplot(aes(x = pf.bin, y=rf, fill = lat.bin5)) +
  geom_boxplot(position ="dodge2") +
  geom_hline(yintercept = 0, linetype = 2) +
  ylim(-7,7) +
  xlab("Years Since Fire") +
  ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
  ggtitle("July") +
  labs(fill = "Latitude") +
  scale_fill_manual(values = cl,
                    labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 16)

#boxplot for march-sept
frc.ann <- alb.frc %>%
  na.omit() %>%
  ggplot(aes(x = pf.bin, y=rf, fill = lat.bin5)) +
  geom_boxplot(position ="dodge2") +
  geom_hline(yintercept = 0, linetype = 2) +
  ylim(-7,7) +
  xlab("Years Since Fire") +
  ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
  ggtitle("March-Sept") +
  labs(fill = "Latitude") +
  scale_fill_manual(values = cl,
                    labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 16)

frc.mar + theme(legend.position = "none") + labs(tag = "a") +
frc.jul + theme(legend.position = "none") + labs(tag = "b") +
frc.ann +labs(tag = "c") 
  
ggsave("figures/postfire_radiative forcing.png",
       width = 12, height = 4, units = "in")


# average over all years since fire
frc.pf <- alb.frc %>%
  na.omit() %>%
  ggplot(aes(x = lat.bin5, y=rf, fill = lat.bin5)) +
  geom_boxplot(position ="dodge2") +
  geom_hline(yintercept = 0, linetype = 2) +
  ylim(-7,7) +
  xlab("Years Since Fire") +
  ylab(expression(paste("Radiative Forcing W",m^-2,sep=""))) +
  ggtitle("Annual Radiative Forcing") +
  labs(fill = "Latitude") +
  scale_fill_manual(values = cl,
                    labels = c("55-60", "60-65", "65-70")) +
  theme_bw(base_size = 16)


#-------------------------------------------------------------------------------------------------------------------------#
# tree cover vs. albedo/rf figure
#-------------------------------------------------------------------------------------------------------------------------#

# need to decide how for post-fire this should be
# could be 1-2 years, or maybe further out, since albedo keeps increasing
tree.mar <- alb.m %>%
  filter(month == 3 & ysf ==1) %>%
  ggplot(aes(x = treecover2, y = d.alb, color = lat.bin5)) +
  geom_point()



#-------------------------------------------------------------------------------------------------------------------------#
# OLD CODE
#-------------------------------------------------------------------------------------------------------------------------#
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