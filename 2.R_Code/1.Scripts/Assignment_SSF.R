#   Script Details                                                          ####

# Author: Waldemar Ortiz-Calo

# Date:2021-04-10 

# Purpose: The purpose of this code is to 

###############################################################################
#   Library / Functions / Data                                              ####

#      Library                                                              ####
require(sp)
require(raster)
require(lubridate)
require(tidyverse)
require(ggplot2)
require(mapview)
require(maptools)
require(leaflet)
require(xtable)
require(broom)
require(amt)
require(magrittr)
require(cowplot)
library(TMB)
library(glmmTMB)
library(FedData)

#      Functions                                                            ####

#      Data                                                                 ####
wolfGPS <- read.csv("1.Data\\Lab10_data\\wolfGPS.csv")
ggplot(wolfGPS, aes(X_COORD1, Y_COORD1, colour = PACK)) +geom_point()

###############################################################################
#   Data Prep                                                               ####
#      NLCD Data                                                            ####

# Importing Cover Data

land_use <- raster("1.Data\\Lab4_data\\landcover16.tif",values = TRUE)

levels(land_use)
# Creating Categoroes and Reformatting data categories
rat <- levels(land_use)[[1]]

rat$category <- c("",
                  "OpenConifer",
                  "ModerateConifer",
                  "ClosedConifer",
                  "Decidious",
                  "MixedForest",
                  "Regeneration",
                  "Herbaceous",
                  "Shrub",
                  "Water",
                  "Rock-Ice",
                  "Cloud",
                  "BurnForest",
                  "BurnGrassland",
                  "BurnShrub",
                  "AlpineHerb",
                  "AlpineShrub")

levels(land_use) <- rat
levels(land_use)


mapview(land_use)

mapview(land_use) + mapview(SpatOb_85_NAD83, zcol="pack", legend = TRUE)


#      Movement Track for One Indiv                                         ####

# Time Data Formatting

Hour <- paste0(Wolf85$HOUR,":00:00")


DateTime <- paste(Wolf85$DATE,Hour)

Wolf85$DateTime <- mdy_hms(paste(DateTime))

# Creating new dataframe with necessary info
Wolf85_trackdata <- Wolf85 %>% 
  filter(!is.na("LAT")) %>% 
  dplyr::select(x        = "LONG",
                y        = "LAT",
                t        = "DateTime",
                id       = "WOLFUID",
                pack     = "PACK")

Wolf85_trackdata$id <- as.character(Wolf85_trackdata$id)

# Arranging by Date
Wolf85_trackdata <- arrange(Wolf85_trackdata,t)


# Making Movement Track
Wolf85_track <- amt::make_track(Wolf85_trackdata, 
                                .x = x,
                                .y = y,
                                .t = t,
                                id = id,
                                crs = sp::CRS("+init=epsg:4326")) %>% 
  amt::transform_coords(sp::CRS("+init=epsg:2153"))

summarize_sampling_rate(Wolf85_track)

Wolf85_track_BaseSamp <- amt::track_resample(Wolf85_track, rate = hours(1),tolerance = minutes(15)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = FALSE)

str(Wolf85_track_BaseSamp, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

Wolf85_track_4hrResamp <- amt::track_resample(Wolf85_track, rate = hours(4),tolerance = minutes(15)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = FALSE)

str(Wolf85_track_4hrResamp, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

###############################################################################
#   AMT For One Indiv                                                       ####
#      Making Movement Track                                                ####

# Arranging by Date
Indiv85 <- arrange(Indiv85,t)

# Making Movement Track
dat_1 <- amt::make_track(Indiv85, x, y, t, crs = sp::CRS("+init=epsg:4326")) %>% 
  amt::transform_coords(sp::CRS("+init=epsg:3857"))

# Checking Sampling Frequency
summarize_sampling_rate(dat_1)

# Resampling
stps <- amt::track_resample(dat_1, rate = hours(2),tolerance = hours(1)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = FALSE)

# Step Attempt
mapview(land_use)+mapview(SpatOb_85_NAD83)

Grass <- land_use == 12
plot(Grass)

mapview(Grass)+mapview(SpatOb_85_NAD83)
eda1 <- stps %>% extract_covariates(Grass, where ="start") %>% mutate(landuse = factor(Grass, levels = c(0, 1), labels = c("other", "OpenCon")))

extract_covariates(stps,land_use)

head(eda1)


p1 <- eda1 %>% dplyr::select(landuse, tod = tod_end_, sl_, ta_) %>%
  gather(key, val, -landuse, -tod) %>%
  filter(key == "sl_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ landuse, nrow = 2) +
  xlab("Step length [m]") + theme_light() +
  ylab("Density") +
  theme(legend.title = element_blank())

p2 <- eda1 %>% dplyr::select(landuse, tod = tod_end_, sl_, ta_) %>%
  gather(key, val, -landuse, -tod) %>%
  filter(key == "ta_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ landuse, nrow = 2) +
  xlab("Turn angle") + theme_light() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank())

pg1 <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), rel_widths = c(1, 1))


leg <- get_legend(p1)
plot_grid(pg1, leg, rel_widths = c(1, 0.1))

levels(land_use)
