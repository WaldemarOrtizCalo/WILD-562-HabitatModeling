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

# Wolf Data
wolfGPS <- read.csv("1.Data\\Lab10_data\\wolfGPS.csv")
ggplot(wolfGPS, aes(X_COORD1, Y_COORD1, colour = PACK)) +geom_point()

# Fixing Time

Hour <- paste0(wolfGPS$HOUR,":00:00")

DateTime <- paste(wolfGPS$DATE,Hour)

wolfGPS$DateTime <- mdy_hms(paste(DateTime))

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

#        [Conifer Layer]                                                    ####

Conifer <- land_use == 1|land_use == 2|land_use == 3
names(Conifer) <- "Conifer"
plot(Conifer)

#        [Open Layer]                                                       ####

Open <- land_use == 7|land_use == 8|land_use == 10|land_use == 12|land_use == 13 |land_use == 14 |land_use == 15|land_use == 16
names(Open) <- "Open"
plot(Open)

#      Elev Data                                                            ####
#        [Elev]                                                             ####
# Importing Elevation Raster
elev <- raster("1.Data\\Lab1_data\\Elevation2.tif")

#Changing Extent and Name
elev <- projectRaster(elev,land_use)
names(elev) <- "elev"
plot(elev)

# Checking Raster 
mapview(elev)

#        [TRI]                                                              ####
tri <- terrain(elev,opt="tri")
names(tri) <- "tri"

# Checking Raster 
mapview(tri)
###############################################################################
#   AMT For One Individual                                                  ####
#      Movement Track for One Individual                                    ####

# Time Data Formatting

Wolf85 <- subset(wolfGPS, WOLFUID == 85)

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


# Making Basic Movement Track Tibble
Wolf85_track <- amt::make_track(Wolf85_trackdata, 
                                .x = x,
                                .y = y,
                                .t = t,
                                id = id,
                                crs = sp::CRS("+init=epsg:4326")) %>% 
  amt::transform_coords(sp::CRS("+init=epsg:2153"))

summarize_sampling_rate(Wolf85_track)

#      No Temporal Resampling [Doesnt work]                                 ####
#        [Creating Track]                                                   ####

# Base Track with no Temporal Resampling
Wolf85_track_BaseSamp <- amt::track_resample(Wolf85_track, rate = hours(1),tolerance = minutes(15)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = FALSE)

str(Wolf85_track_BaseSamp, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

Wolf85_track_BaseSamp[Wolf85_track_BaseSamp$sl_ == 0]
#        [Extracting Covariates for random points and attaching]            ####

# Doing it for the Base Sampling rate of 1 hour
Wolf85_FinalData <-Wolf85_track_BaseSamp %>% amt::random_steps(n = 9) %>%
  amt::extract_covariates(Conifer) %>%
  amt::time_of_day(include.crepuscule = FALSE) %>%
  mutate(log_sl_ = log(sl_)) -> d1

Wolf85_FinalData$caseF <-as.factor(Wolf85_FinalData$case_)

ggplot(Wolf85_FinalData, aes(x2_ ,y2_, colour = caseF)) + geom_point(aes(size = caseF, colour = caseF))

#        [Models]                                                           ####

m3 <- Wolf85_FinalData %>% amt::fit_issf(case_ ~ Conifer + sl_ + Conifer:tod_end_+ sl_:tod_end_ + strata(step_id_))
m2 <- Wolf85_FinalData %>% amt::fit_issf(case_ ~ Conifer + log_sl_ + Conifer:tod_end_+ log_sl_:tod_end_ + strata(step_id_))
m1 <- Wolf85_FinalData %>% amt::fit_issf(case_ ~ Conifer + log_sl_ + sl_ + Conifer:tod_end_+ log_sl_:tod_end_ + sl_:tod_end_ + strata(step_id_))

#      2hr Temporal Resampling                                              ####

#        [Creating Track]                                                   ####

# Track with Temporal Resampling
Wolf85_track_2hrResamp <- amt::track_resample(Wolf85_track, rate = hours(2),tolerance = minutes(30)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = FALSE)

str(Wolf85_track_4hrResamp, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

#        [Extracting Covariates for random points and attaching]            ####

# Doing it for the Base Sampling rate of 1 hour
Wolf85_2hFinalData <-Wolf85_track_2hrResamp %>% amt::random_steps(n = 9) %>%
  amt::extract_covariates(Conifer) %>%
  amt::time_of_day(include.crepuscule = FALSE) %>%
  mutate(log_sl_ = log(sl_)) -> d1

Wolf85_2hFinalData$caseF <-as.factor(Wolf85_2hFinalData$case_)

ggplot(Wolf85_2hFinalData, aes(x2_ ,y2_, colour = caseF)) + geom_point(aes(size = caseF, colour = caseF))

#        [Models]                                                           ####

# Models
m3 <- Wolf85_2hFinalData %>% amt::fit_issf(case_ ~ Conifer + sl_ + Conifer:tod_end_+ sl_:tod_end_ + strata(step_id_))
m2 <- Wolf85_2hFinalData %>% amt::fit_issf(case_ ~ Conifer + log_sl_ + Conifer:tod_end_+ log_sl_:tod_end_ + strata(step_id_))
m1 <- Wolf85_2hFinalData  %>% amt::fit_issf(case_ ~ Conifer + log_sl_ + sl_ + Conifer:tod_end_+ log_sl_:tod_end_ + sl_:tod_end_ + strata(step_id_))

# AIC
AIC(m1$model, m2$model, m3$model)

# Top Model
summary(m3)

###############################################################################
#   Analysis for Multiple Individuals                                       ####
#      Subsetting Individuals                                               ####

MultiIndiv <- wolfGPS %>%
  filter(!is.na(`X_COORD1`)) %>%
  dplyr::select(x = `LONG`, y = `LAT`,
                t = `DateTime`, id = `WOLFUID`, pack = `PACK`) %>%
  filter(id %in% c(81,87))

MultiIndiv <- MultiIndiv %>% nest(-id)
MultiIndiv

#      Making Tracks                                                        ####

MultiIndiv <- MultiIndiv %>%
  mutate(trk = map(data, function(d) {
    amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:4326")) %>%
      amt::transform_coords(sp::CRS("+init=epsg:2153"))}))

# Sampling Rate
MultiIndiv %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% unnest(cols = c(sr))

# Extracting Covariates
m1 <- MultiIndiv %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = hours(2),tolerance = minutes(30)) %>%
      amt::filter_min_n_burst() %>%
      amt::steps_by_burst() %>% amt::random_steps(n=3) %>%
      amt::extract_covariates(Conifer, where = "both" )%>%
      mutate(Conifer_end = factor(Conifer_end))
  }))



#      Models                                                               ####

# Model
m4 <- m1 %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ Conifer_end + 
                                                       strata(step_id_))))
m4

m4$fit[[1]]$model

# Population Model 
d2 <- m4 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
  dplyr::select(id, coef) %>% unnest(cols = coef) %>%
  mutate(id = factor(id)) %>% group_by(term) %>%
  summarize(
    mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
    ymax = mean + 1.96 * sd(estimate))

d2$x <- 1:nrow(d2)
d2

# Plotting
p1 <- m4 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model, conf.int = T))) %>%
  dplyr::select(id, coef) %>% 
  unnest(cols = coef) %>% 
  mutate(id = factor(id)) %>%
  ggplot(., aes(x = term, y = estimate, group = id, col = id)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Habitat", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops"))

p1

###############################################################################
#   Mixed-Effect cLogit Model                                               ####

# Prepping Data
cLogitData <- MultiIndiv %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = hours(2), tolerance =minutes(30)) %>%
      amt::filter_min_n_burst() %>%
      amt::steps_by_burst() %>% amt::random_steps(n=3) %>%
      amt::extract_covariates(Conifer, where = "both") %>%
      mutate(Conifer_end = factor(Conifer_end))
  })) %>%
  dplyr::select(id, steps) %>%
  unnest()

# ID habitat
cLogitData$landuseName = ifelse(cLogitData$Conifer_end == 1, "Conifer","Other")

# Model
model1 <- glm(case_~ I(Conifer_end), data=cLogitData,family=binomial(link="logit"))
summary(model1)
coef(model1)

# Plot
naive_glm <- broom::tidy(model1) %>% 
  filter(!term=="(Intercept)") 

figNaive <- naive_glm %>%
  ggplot(., aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = estimate - 1.96*std.error, ymax = estimate +1.96*std.error)) +
  labs(x = "Habitat", y = "Relative selection Strength") + 
  theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops")) + geom_hline(yintercept = 0, lty = 2) + ylim(c(-3.75,1.5))

figNaive

fig5 <- plot_grid(p1, figNaive)

fig5

###############################################################################
#   Naive cLogit Model                                                      ####

require(survival)

# Data
NaiveClogit <- cLogitData

# Look at the number of step_id_'s for each id
NaiveClogit %>% group_by(step_id_) %>% summarize(n=n())

# Adding Strata
NaiveClogit$stratum <- paste(cLogitData$id, cLogitData$step_id_)

# Model
clogit1 <- clogit(case_ ~ I(Conifer_end) + strata(stratum), data = NaiveClogit)
summary(clogit1)
coef(clogit1)

# Plots

# tidy up coefficients
clogit_1 <- broom::tidy(clogit1) %>% 
  filter(!term=="(Intercept)") 
## make a figure
figclogit1 <- clogit_1 %>%
  ggplot(., aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = estimate+ 1.96*std.error, ymax = estimate-1.96*std.error)) +
  labs(x = "Habitat", y = "Relative selection Strength") + 
  theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops")) + geom_hline(yintercept = 0, lty = 2) + ylim(c(-3.75,1.5))

figclogit1
###############################################################################
#   Coxme Model                                                             ####
require(coxme)

# Data
CoxmeData <- NaiveClogit

CoxmeData$time_ <- ifelse(CoxmeData$case_ == 0, 2, 1)   #2 for control, 1 for case
table(CoxmeData$time_, CoxmeData$case_)

# Model
clogitM1<- coxme(Surv(time_,case_) ~ I(Conifer_end) + strata(stratum) + (1|id), data=CoxmeData)
AIC(clogitM1)

summary(clogitM1)

###############################################################################
#   glmmTMB Model                                                           ####

# Data
TMBData <- NaiveClogit

# Model
TMBm1 <- glmmTMB(case_~ I(Conifer_end) + (1|step_id_) + (0 + I(Conifer_end)|id),
                 family = poisson, data = TMBData, doFit = FALSE)

TMBm1$parameters$theta[1] <-log(1e3)

TMBm1$mapArg <-list(theta=factor(c(NA, 1:15)))
glmm.TMB.random <- glmmTMB::fitTMB(TMBm1)

###############################################################################
###############################################################################
###############################################################################
#   Final Project                                                           ####
#      Prepping Data                                                        ####

# Wolves to Select
# Cascade 85
# Red Deer 42
# Ranch 65
# Wildhorse NA
# Bow Valley NA

FinalData <- wolfGPS %>%
  filter(!is.na(`X_COORD1`)) %>%
  dplyr::select(x = `LONG`, y = `LAT`,
                t = `DateTime`, id = `WOLFUID`, pack = `PACK`) %>%
  filter(id %in% c(85,42,65))

FinalData <- FinalData %>% nest(-id)
FinalData

#      Making Tracks and Extracting Covariates                              ####

FinalData <- FinalData %>%
  mutate(trk = map(data, function(d) {
    amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:4326")) %>%
      amt::transform_coords(sp::CRS("+init=epsg:2153"))}))

# Sampling Rate
FinalData %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% unnest(cols = c(sr))

# Extracting Covariates
Final_Tracks <- FinalData %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = hours(4),tolerance = minutes(30)) %>%
      amt::filter_min_n_burst() %>%
      amt::steps_by_burst() %>% amt::random_steps(n=9) %>%
      amt::extract_covariates(Open, where = "both" )%>%
      amt::extract_covariates(Conifer, where = "both" )%>%
      amt::extract_covariates(elev, where = "both" )%>%
      amt::extract_covariates(tri, where = "both" )%>%
      mutate(Open_end = factor(Open_end))%>%
      mutate(Conifer_end = factor(Conifer_end))
  }))

#      Models                                                               ####
#        [iSSF]                                                             ####

# Model
iSSF_Model<- Final_Tracks %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ tri_end + Open_end +
                                                       strata(step_id_))))
iSSF_Model

iSSF_Model$fit[[3]]$model

AIC(iSSF_Model$fit[[1]]$model,iSSF_Model$fit[[2]]$model,iSSF_Model$fit[[3]]$model)
# Population Model 
iSSF_Pop <- iSSF_Model %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
  dplyr::select(id, coef) %>% unnest(cols = coef) %>%
  mutate(id = factor(id)) %>% group_by(term) %>%
  summarize(
    mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
    ymax = mean + 1.96 * sd(estimate))

iSSF_Pop$x <- 1:nrow(iSSF_Pop)
iSSF_Pop

# Plotting
iSSF_Plot <- iSSF_Model %>% mutate(coef = map(fit, ~ broom::tidy(.x$model, conf.int = T))) %>%
  dplyr::select(id, coef) %>% 
  unnest(cols = coef) %>% 
  mutate(id = factor(id)) %>%
  ggplot(., aes(x = term, y = estimate, group = id, col = id)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = iSSF_Pop, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = iSSF_Pop, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Habitat Covariates", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels = c("Open", "TRI"))

iSSF_Plot

#        [Mixed Effect Clogit]                                              ####

# Prepping Data
cLogitData <- FinalData %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = hours(4),tolerance = minutes(30)) %>%
      amt::filter_min_n_burst() %>%
      amt::steps_by_burst() %>% amt::random_steps(n=9) %>%
      amt::extract_covariates(Open, where = "both" )%>%
      amt::extract_covariates(Conifer, where = "both" )%>%
      amt::extract_covariates(elev, where = "both" )%>%
      amt::extract_covariates(tri, where = "both" )%>%
      mutate(Conifer_end = factor(Conifer_end)) %>% 
      mutate(Open_end = factor(Open_end))
  })) %>%
  dplyr::select(id, steps) %>%
  unnest()

# ID habitat
cLogitData$landuseName = ifelse(cLogitData$Open_end == 1, "Open","Other")

# Model
model1 <- glm(case_~ I(Open_end) + I(tri_end), data=cLogitData,family=binomial(link="logit"))
summary(model1)
coef(model1)

# Plot
naive_glm <- broom::tidy(model1) %>% 
  filter(!term=="(Intercept)") 

figNaive <- naive_glm %>%
  ggplot(., aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = estimate - 1.96*std.error, ymax = estimate +1.96*std.error)) +
  labs(x = "Habitat", y = "Relative selection Strength") + 
  theme_light() +
  scale_x_discrete(labels = c("Open", "TRI")) + geom_hline(yintercept = 0, lty = 2) + ylim(c(-3.75,1.5))

figNaive

fig5 <- plot_grid(iSSF_Plot, figNaive)

fig5

#        [Naive Clogit]                                                     ####

require(survival)

# Data
NaiveClogit <- cLogitData

# Look at the number of step_id_'s for each id
NaiveClogit %>% group_by(step_id_) %>% summarize(n=n())

# Adding Strata
NaiveClogit$stratum <- paste(cLogitData$id, cLogitData$step_id_)

# Model
clogit1 <- clogit(case_ ~ I(Open_end)+I(tri_end) + strata(stratum), data = NaiveClogit)
summary(clogit1)
coef(clogit1)

# Plots

# tidy up coefficients
clogit_1 <- broom::tidy(clogit1) %>% 
  filter(!term=="(Intercept)") 

## make a figure
figclogit1 <- clogit_1 %>%
  ggplot(., aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = estimate+ 1.96*std.error, ymax = estimate-1.96*std.error)) +
  labs(x = "Habitat", y = "Relative selection Strength") + 
  theme_light() +
  scale_x_discrete(labels = c("Open", "TRI")) + geom_hline(yintercept = 0, lty = 2) + ylim(c(-3.75,1.5))

figclogit1

fig6 <- plot_grid(iSSF_Plot, figclogit1)

fig6

#        [Coxme]                                                            ####

# Data
CoxmeData <- NaiveClogit

CoxmeData$time_ <- ifelse(CoxmeData$case_ == 0, 2, 1)   #2 for control, 1 for case
table(CoxmeData$time_, CoxmeData$case_)

# Model
clogitM1<- coxme(Surv(time_,case_) ~ I(Open_end) + strata(stratum) + (1|id), data=CoxmeData)
clogitM2<- coxme(Surv(time_,case_) ~ I(tri_end) + strata(stratum) + (1|id), data=CoxmeData)
clogitM3<- coxme(Surv(time_,case_) ~ I(Open_end) + I(tri_end) + strata(stratum) + (1|id), data=CoxmeData)


clogitM4<- coxme(Surv(time_,case_) ~ I(Open_end) + I(tri_end) + (Open_end:sl_) + (tri_end:sl_) + strata(stratum) + (1|id), data=CoxmeData)
clogitM5<- coxme(Surv(time_,case_) ~ I(Open_end) + I(tri_end) + strata(stratum) + (1|id), data=CoxmeData)
clogitM6<- coxme(Surv(time_,case_) ~ I(Open_end) + I(tri_end) + strata(stratum) + (1|id), data=CoxmeData)
clogitM7<- coxme(Surv(time_,case_) ~ I(Open_end) + I(tri_end) + strata(stratum) + (1|id), data=CoxmeData)

AIC(clogitM1,clogitM2,clogitM3,clogitM4)

summary(clogitM4)











#        [TMB]                                                              ####

table(cLogitData$Open_end)

TMBm1 <- glmmTMB(case_~ I(Open_end) + I(tri_end) +(1|step_id_) + (0 + I(Open_end)|id), family = poisson, data = cLogitData, doFit = FALSE)


TMBm1$parameters$theta[1] <-log(1e3)

TMBm1$mapArg <-list(theta=factor(c(NA, 1:3)))
glmm.TMB.random <- glmmTMB::fitTMB(TMBm1)
summary(glmm.TMB.random)
