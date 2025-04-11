######################################################
#
# extract radiative kernels from Bright & O'Halloran
# to use with field obs at various locations
#
# MML 03/19/25
#####################################################

rm (list=ls())
# load required packages
library(terra)
#library(ncdf4)

# for now working from the directory for a siberia fire rf manuscript
# because it contains the CACK files
# set working directory with options for Mac or PC workstations
ifelse(Sys.info()[1]=="Windows",
       setwd("G:/My Drive/Documents/research/manuscripts/siberia_albedo"),
       setwd("~/Library/CloudStorage/GoogleDrive-mloranty@colgate.edu/My Drive/Documents/research/manuscripts/siberia_albedo"))


# read in the cack file
# transpose and correct longitude values as well
rf <- rotate(t(rast("data/CACKv1.0/CACKv1.0.nc", "CACK")))
# correct set latitude to correct extent
ext(rf) <- c(-180,180,-90,90)
#crs(rf) <- "EPSG:4326"

site.crds <- as.data.frame(matrix(c(161.497878, 68.513201, -149.375871, 64.220625), nrow=2, byrow = T))
names(site.crds) <- c("lat", "lon")

# extract the kernel values and average across the 16 years in the CACK data set
site.kern <- terra::extract(rf,site.crds) %>%
  pivot_longer(cols = starts_with("CACK"),
    names_to = c("month", "year"),
    names_sep = "_Year=",
    #names_pattern = "*[0-9]+*[0-9]+",
    values_to = "kernel") %>%
  mutate(site = case_when(
    ID == 1 ~ "pp",
    ID == 2 ~ "afe"), 
    month = as.numeric(substr(month, 12,13))) %>%
  group_by(site, month) %>%
  summarise(kernel = mean(kernel, na.rm = T))

# write output to file
write.csv(site.kern, "results/pleistocene_park_cack.csv")
ext(rf) <- c(-180,180,-90,90)