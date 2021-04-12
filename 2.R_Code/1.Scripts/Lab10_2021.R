#title: "Lab 10: Advanced SSF Models - integrated SSF's"
#author: "Mark Hebblewhite and Eric Palm"



## Preliminaries: setting packages
require(sp)
require(raster)
require(lubridate)
require(tidyverse)
require(ggplot2)
require(mapview)
require(maptools)
require(leaflet)
require(broom)
require(amt)
require(magrittr)
require(cowplot)
install.packages('TMB', type = 'source')
library(TMB)
install.packages("glmmTMB", type="source") ## note this really needs to be installed from source
library(glmmTMB)

#A note about tidyverse.  Tidyverse is REALLY useful, but, sometimes has 'hidden' conflicts with other packages.  Check which specific conflicts there are with:
tidyverse_conflicts()


# Animal movement tools (amt)
#Lets get started....
fisher <- read.csv("1.Data\\Lab10_data\\Martes pennanti LaPoint New York.csv")
head(fisher)
str(fisher)
ggplot(fisher, aes(location.long, location.lat, colour = individual.local.identifier)) + geom_point()

## data manipulation
fisher1<-fisher[complete.cases(fisher[4:5]),]
xy <- fisher1[, c(4,5)]
fisherSP <- SpatialPointsDataFrame(coords = xy, data = fisher1, proj4string = sp::CRS("+init=epsg:4326"))
mapview(fisherSP, zcol="individual.local.identifier", legend = TRUE, cex=5, lwd=2, map.type = c("OpenStreetMap.DE", "Esri.WorldShadedRelief"))

#######
#Next we will bringing Fisher data from a single individual, ID 1016, Fisher M1 (known as "RickyT"), into a MOVE object.
dat <- read_csv("1.Data\\Lab10_data\\Martes pennanti LaPoint New York.csv") %>%
   filter(!is.na(`location-lat`)) %>%
   dplyr::select(x = `location-long`, y = `location-lat`,
           t = `timestamp`, id = `tag-local-identifier`) %>%
    filter(id %in% c(1465, 1466, 1072, 1078, 1016, 1469)) # for example 2
   dat_1 <- dat %>% filter(id == 1016)
head(dat_1)


# Make a track just with the data from RickyT. 
dat_1 <- amt::make_track(dat_1, x, y, t, crs = sp::CRS("+init=epsg:4326")) %>%
     amt::transform_coords(sp::CRS("+init=epsg:3857"))

summarize_sampling_rate(dat_1)


## amt::track_resample 
stps <- amt::track_resample(dat_1, rate = minutes(10), tolerance = minutes(1)) %>%
  filter_min_n_burst(min_n = 3) %>% steps_by_burst() %>%
  time_of_day(include.crepuscule = FALSE)

str(stps, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)


## Obtaining the NLCD data from `FedData`
require("devtools")
devtools::install_github("ropensci/FedData")
library(FedData)


# Create mask.raster() extent
## Create a Mask Raster based on the extent of Fisher. Note that I made it arbitrarily larger. 
fisherSP@proj4string
fisherSP2 <-spTransform(fisherSP, CRS("+init=epsg:3857"))
extent(fisherSP2)
mask.raster<-raster()
extent(mask.raster) <- c(xmin=-8240000, xmax=-8160000, ymin=5260000 , ymax=5295000)

res(mask.raster) = 30
projection(mask.raster)<- "+init=epsg:3857"  
#set all values of mask.raster to zero
mask.raster[]<-0
#projection(mask.raster)
plot(mask.raster)
plot(fisherSP2, add = TRUE)

## Import NLCD data
get_nlcd(mask.raster, label = "landuse", year = 2011, dataset = "Land_Cover", extraction.dir = "Data/extract/", force.redo = TRUE)
land_use <- raster("1.Data\\NLCD\\landuse_NLCD_2011_Land_Cover_L48_nlcd.tif",values = TRUE)
str(land_use@data@values)

#Anchoring the Raster Attribute Table (RAT) and tying it to the layer. 
land_use2<-land_use
land_use2[] <-values(land_use) 
land_use3 <-ratify(land_use2)
str(land_use3@data@values)
hist(land_use3@data@values)
unique(land_use3@data@values)

#Create a new Raster Attribute Table (RAT) with these 15 landcover names in it, and then define the levels of our land_use3 as this new RAT table. 
rat <- levels(land_use3)[[1]]
rat$landcover <- c("Open Water", "Developed-Open", "Developed-Low", "Developed- Medium", "Developed-High", "Barren Land", "Deciduous Forest", "Evergreen Forest", "Mixed Forest", "Shrub/Scrub", "Herbaceuous", "Hay/Pasture", "Cultivated Crops", "Woody Wetlands", "Herbaceuous Wetlands")
levels(land_use3) <- rat
levels(land_use3)

## Reproject fisher to NLCD world web mercator projection `"+init=epsg:3857"` 
fisherSP3 <-spTransform(fisherSP2, "+init=epsg:3857") # "+init=epsg:5070"
#Now lets overlay them. 
plot(land_use3)
plot(fisherSP3, add=TRUE, type="p", pch=12, cex = 0.5)

## Mapvie
mapview(fisherSP3, zcol="individual.local.identifier", legend = TRUE, cex=5, lwd=2, map.type = c("OpenStreetMap.DE", "Esri.WorldShadedRelief")) + land_use

## Exploring Landcover for RickyT
wet <- land_use3 == 90
names(wet) <- "wet"
##Lets zoom into Ricky T
rickyT.raster <- raster()
extent(rickyT.raster) <- c(xmin=-8230000, xmax=-8210000, ymin=5270000, ymax=5280000)
plot(wet, ext = rickyT.raster)
plot(fisherSP3, add=TRUE, type="p", pch=12, cex = 0.5, ext = rickyT.raster)

#RickyT selects 'wet' forests. 

## Exploratory Analyses of Step Lengths and Turning Angles
## amt:: extract_covariates
eda1 <- stps %>% extract_covariates(wet, where ="start") %>% mutate(landuse = factor(wet, levels = c(0, 1), labels = c("other", "forested wetland")))
head(eda1)

# Summary plots of step length, turning angle 
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


# To save the figure, use this
#ggsave("Figures/fig_eda_1_animal.pdf", width = 20, height = 18, units = "cm")

################################################################
## Fitting a Step Selection Function


### amt::random_steps
m1 <-stps %>% amt::random_steps(n = 9) %>%
  amt::extract_covariates(wet) %>%
  amt::time_of_day(include.crepuscule = FALSE) %>%
  mutate(log_sl_ = log(sl_)) -> d1

#To see what the random_steps() function did, take a look at the first 18 rows
head(m1, n=10)
#str(m1)


#visualize the SSF points 
m1$caseF <-as.factor(m1$case_)
ggplot(m1, aes(x2_ ,y2_, colour = caseF)) + geom_point(aes(size = caseF, colour = caseF))

### Fitting SSF models using amt::fit_issf()
m3 <- d1 %>% amt::fit_issf(case_ ~ wet + sl_ + wet:tod_end_+ sl_:tod_end_ + strata(step_id_))
m2 <- d1 %>% amt::fit_issf(case_ ~ wet + log_sl_ + wet:tod_end_+ log_sl_:tod_end_ + strata(step_id_))
m1 <- d1 %>% amt::fit_issf(case_ ~ wet + log_sl_ + sl_ + wet:tod_end_+ log_sl_:tod_end_ + sl_:tod_end_ + strata(step_id_))


### Model Selection
AIC(m1$model, m2$model, m3$model)
summary(m1)
s <- summary(m1$model)$coefficients
s


### Examining Movement Statistics
coef(m1)
plot_sl(m1)
sl_m1<-sl_distr(m1)
## note we will extract the scale and shape parameters of the gamma distribution for later simulating the UD
scale <- sl_m1$params$scale
shape <- sl_m1$params$shape

## differences in step lenght between day/night
ggplot(d1, aes(tod_end_, sl_)) +geom_violin()
sl_data <- as_tibble(d1)
sl_data %>% group_by(tod_end_) %>% summarize(median = median(sl_))


## Simulating the Utilization Distributions 

#Crop out wet for a smaller area around just Ricky T using the amt::bbox() function which describes a bounding box. 
wet_c <- crop(wet, amt::bbox(dat_1, spatial = TRUE, buff = 1e3))

## estimate daytime movement kernel
mk <- amt::movement_kernel(scale, shape, wet_c)
plot(mk)

## Habitat kernel 
hk <- amt::habitat_kernel(list(wet = coef(m1)["wet"]), wet_c)
plot(hk)

## Fit simulated steady-state UD 
system.time(ssud_day <- amt::simulate_ud(
  mk, hk,
  as.numeric(stps[1, c("x1_", "y1_")]),
   n = 1e7))
 plot(ssud_day)

## Transient UD 
system.time(tud_day <- amt::simulate_tud(mk, hk, as.numeric(stps[150, c("x1_", "y1_")]), n = 72, n_rep = 5e3))
plot(tud_day)

## Plot all 4 figures
pllog <- list(
  geom_raster(),
   coord_equal(),
  scale_fill_continuous(low = "white", high = "red", tran = "log10", na.value = "white"),
   scale_y_continuous(expand = c(0, 0)),
   scale_x_continuous(expand = c(0, 0)),
   theme_light(),
   theme(legend.position = "none"))

 pl <- list(
   geom_raster(),
   coord_equal(),
   scale_fill_continuous(low = "white", high = "red", na.value = "white"),
 scale_y_continuous(expand = c(0, 0)),
   scale_x_continuous(expand = c(0, 0)),
   theme_light(),
   theme(legend.position = "none"))

r1 <- data.frame(rasterToPoints(mk))
p1 <- ggplot(r1, aes(x, y, fill = d)) + pllog + ggtitle("Movement kernel (day)")

r2 <- data.frame(rasterToPoints(hk))
p2 <- ggplot(r2, aes(x, y, fill = layer)) + pl + ggtitle("Habitat kernel (day)")

r1 <- data.frame(rasterToPoints(tud_day))
 p3 <- ggplot(r1, aes(x, y, fill = layer)) + pllog + ggtitle("Transient UD (day)")

r1 <- data.frame(rasterToPoints(ssud_day))
p5 <- ggplot(r1, aes(x, y, fill = layer)) + pl + ggtitle("Steady state UD (day)")

cowplot::plot_grid(p1, p2, p3, p5, ncol = 2, labels = "AUTO")
#ggsave("fig_one_animal1.pdf", height = 20, width = 24, units = "cm")

######################################################################
######################################################################
# Fitting SSF Models to Multiple Animals

## REad in fisher data again
dat <- read_csv("Data/Martes pennanti LaPoint New York.csv") %>%
  filter(!is.na(`location-lat`)) %>%
  dplyr::select(x = `location-long`, y = `location-lat`,
                t = `timestamp`, id = `tag-local-identifier`) %>%
  filter(id %in% c(1465, 1466, 1072, 1078, 1016, 1469))

dat_all <- dat %>% nest(-id)
dat_all$sex <- c("f", "f", "f", "m", "m", "m")
dat_all
## create a track,  transform the coordinate reference system using the function amt::transform _ coords.
dat_all <- dat_all %>%
  mutate(trk = map(data, function(d) {
    amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:4326")) %>%
      amt::transform_coords(sp::CRS("+init=epsg:3857"))}))

#And summarize sampling rate
dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% unnest(cols = c(sr))


## NLCD Extraction
#Create a reclassification matrix into 5 categories, with water/wet forests as the reference catetgory: water and wetland forests, developed open spaces, other developed areas, forests and shrubs, and crops. 

land_use3@data@attributes

rcl <- cbind(c(11, 12, 21:24, 31, 41:43, 51:52, 71:74, 81:82, 90, 95),c(1, 1, 2, 3, 3, 3, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 1, 1))
# water, dev open, dev, barren, forest, shrub and herb, crops, wetlands
# 1: water, wetlands
# 2: developed (open)
# 3: developed (other)
# 5: forest, herbaceous
# 8: crops
lu <- reclassify(land_use3, rcl, right = NA)
names(lu) <- "landuse"
plot(lu)
plot(fisherSP3, add=TRUE)


## Resample the tracks from 6 individual fishers to 10 minute steps, filter them into consecutive bursts of 3 or more 10 minute locations, create random steps, extract covariate values from both the start and the end, and convert the landuse covariate to a factor all in one step.
```{r, warning = F}
m1 <- dat_all %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = minutes(10), tolerance = seconds(120)) %>%
      amt::filter_min_n_burst() %>%
      amt::steps_by_burst() %>% amt::random_steps() %>%
      amt::extract_covariates(lu, where = "both") %>%
      mutate(landuse_end = factor(landuse_end))
  }))


### fit a simple SSF model to each individual. 
m4 <- m1 %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ landuse_end +strata(step_id_))))
m4


#Lets look at 1 of these models for individual 1
m4$fit[[1]]$model

## Derive marginal, population level averaged coefficient for landuse across our individual animals. 
d2 <- m4 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
   dplyr::select(id, sex, coef) %>% unnest(cols = coef) %>%
  mutate(id = factor(id)) %>% group_by(term) %>%
  summarize(
     mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
     ymax = mean + 1.96 * sd(estimate))

d2$x <- 1:nrow(d2)
d2

# Visualize in some figures
p1 <- m4 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model, conf.int = T))) %>%
  dplyr::select(id, sex, coef) %>% 
  unnest(cols = coef) %>% 
  mutate(id = factor(id)) %>%
  ggplot(., aes(x = term, y = estimate, group = id, col = id, pch = sex)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Habitat", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops"))

p1

## if we want to save
#ggsave("img/fig_all_animals.pdf", width = 24, height = 12, units = "cm")

###############################################################
# Mixed-effect cLogit Models

# compare model interpretations to naive GLM models of the same kind. 
fisher6 <- dat_all %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = minutes(10), tolerance = seconds(120)) %>%
      amt::filter_min_n_burst() %>%
      amt::steps_by_burst() %>% amt::random_steps() %>%
      amt::extract_covariates(lu, where = "both") %>%
      mutate(landuse_end = factor(landuse_end))
  })) %>%
  dplyr::select(id, steps) %>%
  unnest()

fisher6

## Name the NLCD classes names
head(fisher6$landuse_end)
fisher6$landuseName = ifelse(fisher6$landuse_end == 1, "Wet Forests", 
                             ifelse(fisher6$landuse_end == 2, "Developed Open", 
                             ifelse(fisher6$landuse_end == 3, "Developed Other", 
                             ifelse(fisher6$landuse_end == 5, "Natural", "Crops"))))
table(fisher6$landuseName, fisher6$landuse_end)

## Fit a naive GLM
model1 <- glm(case_~ I(landuse_end), data=fisher6,family=binomial(link="logit"))
## I commented out these next few versions of the models fit to landuseName to make comparisons to the previously fit 6 fisher two-step models more comparable, though we have to then keep track of which landovers 2, 3, 5, and 8 are. 
#model1 <- glm(case_~ I(landuseName=="Developed Open") + I(landuseName=="Developed Other") +I(landuseName=="Natural")+I(landuseName=="Crops"), data=fisher6,family=binomial(link="logit"))
summary(model1)

coef(model1)
naive_glm <- broom::tidy(model1) %>% 
  filter(!term=="(Intercept)") 

## Make Figure
figNaive <- naive_glm %>%
  ggplot(., aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = estimate - 1.96*std.error, ymax = estimate +1.96*std.error)) +
  labs(x = "Habitat", y = "Relative selection Strength") + 
    theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops")) + geom_hline(yintercept = 0, lty = 2) + ylim(c(-3.75,1.5))
figNaive

# Comparing the two-stage averaged SSF coefficients from the fisher SSF model above extracted in the object d2.

fig5 <- plot_grid(p1, figNaive)
fig5


## Fitting a 'Naive' cLogit Model with `surv` package

require(survival)
head(fisher6)
## Look at the number of step_id_'s for each id
fisher6 %>% group_by(step_id_) %>% summarize(n=n())
fisher6$stratum <- paste(fisher6$id, fisher6$step_id_)
#We see that there are 66 rows of data for each step ID because there are 6 individuals. 

clogit1 <- clogit(case_ ~ I(landuse_end) + strata(stratum), data = fisher6)
#clogit1 <- clogit(case_ ~ I(landuseName=="Developed Open") + I(landuseName=="Developed Other") +I(landuseName=="Natural")+I(landuseName=="Crops") + strata(stratum), data = fisher6)
summary(clogit1)
coef(clogit1)

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
plot_grid(p1, figclogit1)

## Mixed-effect cLogit Models with `coxme`

require(coxme)

fisher6$time_ <- ifelse(fisher6$case_ == 0, 2, 1)   #2 for control, 1 for case
table(fisher6$time_, fisher6$case_)
head(fisher6)
clogitM1<- coxme(Surv(time_,case_) ~ I(landuse_end) + strata(stratum) + (1|id), data=fisher6)
AIC(clogitM1)
summary(clogitM1)

#clogitM2<- coxme(Surv(time_,case_) ~ I(landuse_end) + strata(stratum) + (1|id) +(I(landuseName=="Developed Other")|id) , data=fisher6)
#summary(clogitM2)

## mclogit
#mclogitTest2 <- mclogit(cbind(used_, stratumN) ~I(landuse_end), random=~ 1|id, data=fisher6)

v1<-model1$coefficients[2:5]
v2<-coef(clogit1)
v3<-coef(clogitM1)
v4<-d2$mean # note that this has female and male two-stage parameter estimates averaged. 
coefSum <- as.data.frame(cbind(v1, v2, v3, v4))
names(coefSum) <- c("Naive", "clogit", "coxme", "two-stage iSSF")
head(coefSum) # 

## Model Selection

AIC(model1, clogit1, clogitM1)

## Fitting iSSF models with glmmTMB

install.packages("glmmTMB", type = 'source')
library(glmmTMB)
head(fisher6)
table(fisher6$landuse_end)
TMBm1 <- glmmTMB(case_~ I(landuse_end) + (1|step_id_) + (0 + I(landuse_end)|id), family = poisson, data = fisher6, doFit = FALSE)

TMBm1$parameters$theta[1] <-log(1e3)

TMBm1$mapArg <-list(theta=factor(c(NA, 1:15)))
glmm.TMB.random <- glmmTMB::fitTMB(TMBm1)
summary(glmm.TMB.random)

## Compare all 5 models now
coefTMB <- fixef(glmm.TMB.random)
v5 <- coefTMB$cond[2:5]
coefSum2 <- as.data.frame(cbind(coefSum, v5))
names(coefSum2) <- c("Naive", "clogit", "coxme", "two-stage iSSF", "glmmTMB")
coefSum2
# recall again that the coefficients are for males and females for two-stage iSSF, which is why there are 8 rows. 
p1

