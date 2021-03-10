#   Script Details                                                          ####

# Author: Waldemar Ortiz-Calo

# Date:2021-03-01 

# Purpose: 

###############################################################################
#   Library / Functions / Data                                              ####

#      Library                                                              ####
library(tidyverse)
library(ggplot2)
library(plotrix)
library(ks)
library(lattice)
library(adehabitatHR)
library(maptools)
library(foreign)
library(rgdal)
library(sp)
library(raster)
library(plot3D)
library(rasterVis)
library(colorRamps)
library(rgeos)
library(mapview)
library(fasterize)
library(effects)
library(magrittr)
library(adehabitatHS)
library(MASS)
library(MuMIn)
library(AICcmodavg)
library(corrgram)
library(GGally)
library(broom)
library(ROCR)
library(VGAM)
library(caret)
library(DescTools)
library(sandwich)
library(lmtest)
library(merTools)
library(ResourceSelection)
library(Hmisc)
library(plotrix)
library(pander)
library(lattice)
library(sjstats)

#      Functions                                                            ####


#      Data                                                                 ####
#        [Wolf]                                                             ####

Wolf <-shapefile("1.Data/Lab1_data/wolfyht.shp")
plot(Wolf)

#        [ELC Habitat]                                                      ####

ELC_Habitat <- shapefile("1.Data/Lab1_data/elc_habitat.shp")
plot(ELC_Habitat)

#           [HSI - Prey]                                                    ####

# Elk 
HSI_elk <- raster("1.Data\\Lab1_data\\elk_w2.tif")
plot(HSI_elk)

# Moose
HSI_moose <- raster("1.Data\\Lab1_data\\moose_w2.tif")
plot(HSI_moose)

# Deer 
HSI_deer<- raster("1.Data\\Lab1_data\\deer_w2.tif")
plot(HSI_deer)

# Goat 
HSI_goat <- raster("1.Data\\Lab1_data\\goat_w2.tif")
plot(HSI_goat)

# Sheep
HSI_sheep <- raster("1.Data\\Lab1_data\\sheep_w2.tif")
plot(HSI_sheep)

#        [Human Access]                                                     ####

HumanAccess <- shapefile("1.Data\\Lab1_data\\humanacess.shp")
plot(HumanAccess)

HumanAccessDistance <- raster("1.Data\\Lab1_data\\DistFromHumanAccess2.tif")
plot(HumanAccessDistance)

HighHumanAccessDistance <- raster("1.Data\\Lab1_data\\DistFromHighHumanAccess2.tif")
plot(HighHumanAccessDistance)

#        [Elevation]                                                        ####

Elevation <- raster("1.Data\\Lab1_data\\Elevation2.tif")
plot(Elevation)

#           [Landscape]                                                     ####
TRI <- terrain(Elevation,opt="TRI")
plot(TRI)

#        [Visual Checking of Data]                                          ####

mapview(Wolf) + mapview(Elevation) + HSI_elk + HighHumanAccessDistance

###############################################################################
#   Home Range Estimates                                                    ####

#      [Data Prep]                                                          ####

#first convert the spatialpointsdataframe to spatial points object
x<-Wolf@data$EASTING
y<-Wolf@data$NORTHING
xy<-cbind(x,y)

all <- data.frame(as.character(Wolf@data$Pack))
coordinates(all) <- xy
proj4string(all) <-  CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#      [Home Range]                                                         ####
#        [Bow Valley]                                                       ####

bv.data<-Wolf[Wolf@data$Pack=="Bow valley",]
x<-bv.data@data$EASTING
y<-bv.data@data$NORTHING
xy<-cbind(x,y)

bv <- data.frame(as.character(bv.data@data$NAME))
coordinates(bv) <- xy
proj4string(bv) <-  CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Fit 95% mpc
cp.bow <- mcp(bv, percent=95)
plot(bv, col="black")
plot(cp.bow[cp.bow@data$id=="63",], col="blue",add=TRUE)
plot(cp.bow[cp.bow@data$id=="87",], col="red",add=TRUE,)
plot(cp.bow[cp.bow@data$id=="44",], col="green",add=TRUE)
plot(bv, col="black", add=TRUE)

#        [Red Deer]                                                         ####

rd.data<-Wolf[Wolf@data$Pack=="Red Deer",]
x<-rd.data@data$EASTING
y<-rd.data@data$NORTHING
xy<-cbind(x,y)

rd <- data.frame(as.character(rd.data@data$NAME))
coordinates(rd) <- xy
proj4string(rd) <-  CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Fit 95% mpc
cp.red <- mcp(rd, percent=95)
plot(rd, col="black")
plot(cp.red[cp.red@data$id=="63",], col="blue",add=TRUE)
plot(cp.red[cp.red@data$id=="87",], col="red",add=TRUE,)
plot(cp.red[cp.red@data$id=="44",], col="green",add=TRUE)
plot(rd, col="black", add=TRUE)

#        [All]                                                              ####

# Fit 95% mpc
cp.all <- mcp(all, percent=95)

# Checking the area
mcp.area(all, percent=seq(50, 100, by=5))

# Creating a Kernel
allUD <- kernelUD(all, grid=30, extent=0.5, same4all=TRUE) # reference grid

# Home Range Contours
homerangeALL <- getverticeshr(allUD)

#      [Plots]                                                              ####

# MCP Plots
plot(Wolf, col="black")
plot(cp.all[cp.all@data$id=="Bow valley",], col="blue",add=TRUE)
plot(cp.all[cp.all@data$id=="Red Deer",], col="green",add=TRUE)
plot(Wolf, col="black", add=TRUE)

# Kernel Plot for both packs
plot(homerangeALL, col=2:3)

###############################################################################
#   RSF Data Prep                                                           ####
#      [Sampling Availability]                                              ####
#        [Extracting Polygons from Pack Home Ranges]                        ####

red.deerPOLY<-homerangeALL[homerangeALL@data$id=="Red Deer",]

bow.valleyPOLY<-homerangeALL[homerangeALL@data$id=="Bow valley",]

#        [Randomly Sampling Available Locations within Home Ranges]         ####

rd.avail<-spsample(red.deerPOLY, 1000, "random")

plot(rd.avail)

bv.avail<-spsample(bow.valleyPOLY, 1000, "random")

plot(bv.avail)

#      [Raster Covariates]                                                  ####

# Raster Stack

RasterStack <- stack(HSI_deer,
                     HSI_elk,
                     HSI_moose,
                     HSI_goat,
                     HSI_sheep,
                     Elevation,
                     HumanAccess,
                     HighHumanAccessDistance,
                     TRI)

plot(RasterStack)
#        [Bow Valley ]                                                      ####

# Used
cov.outBV <- raster::extract(RasterStack, bv.data)

# Available
cov.availBV <- raster::extract(RasterStack, bv.avail)

#        [Red Deer]                                                         ####
 
# Used
cov.outRD <- raster::extract(RasterStack, rd.data)

# Available
cov.availRD <- raster::extract(RasterStack, rd.avail)

#        [Merging Used Data]                                                ####

# Red Deer
rdused <- as.data.frame(cov.outRD)
rdused$pack <- c("Red Deer")

# Bow Valley
bvused <- as.data.frame(cov.outBV)
bvused$pack <- c("Bow Valley")

# Merge used and give a 1 indiciator
Wolf_Used <- merge(rdused, bvused, all.x= TRUE, all.y = TRUE) %>% na.omit
Wolf_Used$used <- 1

#        [Merging Available Data]                                           ####

# Red Deer
rdavail <- as.data.frame(cov.availRD)
rdavail$pack <- c("Red Deer")

# Bow Valley 
bvavail <- as.data.frame(cov.availBV)
bvavail$pack <- c("Bow Valley")

# Merge the two availability samples together and give a 0 indiciator
Wolf_Avail <- rbind(rdavail, bvavail)
Wolf_Avail$used <- 0

#        [Final Dataset]                                                    ####

UsedAvail <- rbind(Wolf_Used, Wolf_Avail)

###############################################################################
#   RSF Models                                                              ####
#      [Graphical Exploration]                                              ####

# Used vs Available Per Pack
table(UsedAvail$used, UsedAvail$pack)

# Used vs Available Per HSI Covariate
table(UsedAvail$used, UsedAvail$deer_w2)
table(UsedAvail$used, UsedAvail$elk_w2)
table(UsedAvail$used, UsedAvail$moose_w2)
table(UsedAvail$used, UsedAvail$goat_w2)
table(UsedAvail$used, UsedAvail$sheep_w2)

# Factor
UsedAvail$usedFactor <- factor(UsedAvail$used, labels=c('0','1'))

#      [Univariate]                                                         ####

# Elevation
RSF_Elev <- glm(used ~ Elevation2, family=binomial(logit), data=UsedAvail)
summary(RSF_Elev)
confint(RSF_Elev)
confint.default(RSF_Elev)

# Predction 
elevBnp = 0:3000 ## creates a new vector elevBnp with ranges from 0 - 3000 in it.
elevPred = predict(RSF_Elev, newdata=data.frame(Elevation2=elevBnp), type = "response") ## uses the predict function to predict Y values given the model object elev
hist(elevPred)
plot(elevBnp, elevPred, type="l", ylim = c(0,1.0), ylab= "Pr(Used)")

#      [MultiVariate]                                                       ####

# Multi-variate model
Multi <- glm(used~Elevation2 + DistFromHighHumanAccess2 , data =UsedAvail, family= binomial(logit))
summary(Multi)

# Checking For Collinearity 
cor.test(UsedAvail$Elevation2, UsedAvail$DistFromHighHumanAccess2)

pairs(~Elevation2+DistFromHighHumanAccess2, data=UsedAvail, main="Scatterplot Matrix")

###############################################################################
#   Final Model                                                             ####
#      [Data Prep]                                                          ####
#        [Colinearity]                                                      ####

pairs(~Elevation2 +
        DistFromHighHumanAccess2 +
        tri +
        deer_w2 +
        elk_w2 +
        moose_w2+
        goat_w2 +
        sheep_w2, data=UsedAvail, main="Scatterplot Matrix")

ggcorrplot <- ggcorr(UsedAvail[1:8], label = TRUE)
ggcorrplot

#      [Univariate Models]                                                  ####
tri <- glm(used~tri , data =UsedAvail, family= binomial(logit))
summary(tri)

ggplot(UsedAvail, aes(x=tri, y=used)) + 
  geom_point() + 
  stat_smooth(method="glm", method.args = list(family="binomial"))

deer_model <- glm(used~ deer_w2 , data =UsedAvail, family= binomial(logit))
summary(deer_model)

ggplot(UsedAvail, aes(x=deer_w2, y=used)) + 
  geom_point() + 
  stat_smooth(method="glm", method.args = list(family="binomial"))

#      [Multi-Variate Models]                                               ####
options(na.action = "na.omit")

prey_model <- glm(used~ deer_w2 + elk_w2 + moose_w2 + sheep_w2 + goat_w2
                    , data =UsedAvail, family= binomial(logit), )
summary(prey_model)

step <- stepAIC(prey_model)
step$anova


#           [Model Selection]                                               ####

# List of Models
m.biotic <- list()
head(m.biotic)

m.biotic[[1]]  <- glm(used ~ 1, family=binomial(logit), data=UsedAvail)
m.biotic[[2]]  <- glm(used ~ elk_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[3]]  <- glm(used ~ deer_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[4]]  <- glm(used ~ moose_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[5]]  <- glm(used ~ sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[6]]  <- glm(used ~ goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[7]]  <- glm(used ~ moose_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[8]]  <- glm(used ~ deer_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[9]]  <- glm(used ~ elk_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[10]] <- glm(used ~ elk_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[11]] <- glm(used ~ deer_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[12]] <- glm(used ~ moose_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[13]] <- glm(used ~ sheep_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[14]] <- glm(used ~ DistFromHighHumanAccess2, family=binomial(logit), data=UsedAvail)
m.biotic[[15]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[16]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[17]] <- glm(used ~ DistFromHighHumanAccess2+sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[18]] <- glm(used ~ DistFromHighHumanAccess2+goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[19]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[20]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[21]] <- glm(used ~ DistFromHighHumanAccess2+elk_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[22]] <- glm(used ~ DistFromHighHumanAccess2+elk_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[23]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[24]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[25]] <- glm(used ~ DistFromHighHumanAccess2+sheep_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[26]] <- glm(used ~ DistFromHighHumanAccess2, family=binomial(logit), data=UsedAvail)
m.biotic[[27]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[28]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[29]] <- glm(used ~ DistFromHighHumanAccess2+sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[30]] <- glm(used ~ DistFromHighHumanAccess2+goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[31]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[32]] <- glm(used ~ tri+DistFromHighHumanAccess2, family=binomial(logit), data=UsedAvail)
m.biotic[[33]] <- glm(used ~ tri+DistFromHighHumanAccess2+deer_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[34]] <- glm(used ~ tri+DistFromHighHumanAccess2+moose_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[35]] <- glm(used ~ tri+DistFromHighHumanAccess2+sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[36]] <- glm(used ~ tri+DistFromHighHumanAccess2+goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[37]] <- glm(used ~ tri+DistFromHighHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[38]] <- glm(used ~ tri+DistFromHighHumanAccess2+deer_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[39]] <- glm(used ~ tri+DistFromHighHumanAccess2+elk_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[40]] <- glm(used ~ tri+DistFromHighHumanAccess2+elk_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[41]] <- glm(used ~ tri+DistFromHighHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[42]] <- glm(used ~ tri+DistFromHighHumanAccess2+moose_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[43]] <- glm(used ~ tri+DistFromHighHumanAccess2+sheep_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[44]] <- glm(used ~ tri+DistFromHighHumanAccess2, family=binomial(logit), data=UsedAvail)
m.biotic[[45]] <- glm(used ~ tri+DistFromHighHumanAccess2+deer_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[46]] <- glm(used ~ tri+DistFromHighHumanAccess2+moose_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[47]] <- glm(used ~ tri+DistFromHighHumanAccess2+sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[48]] <- glm(used ~ tri+DistFromHighHumanAccess2+goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[49]] <- glm(used ~ tri+DistFromHighHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[50]] <- glm(used ~ tri, family=binomial(logit), data=UsedAvail)
m.biotic[[51]] <- glm(used ~ tri+deer_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[52]] <- glm(used ~ tri+moose_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[53]] <- glm(used ~ tri+sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[54]] <- glm(used ~ tri+goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[55]] <- glm(used ~ tri+moose_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[56]] <- glm(used ~ tri+deer_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[57]] <- glm(used ~ tri+elk_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[58]] <- glm(used ~ tri+elk_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[59]] <- glm(used ~ tri+deer_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[60]] <- glm(used ~ tri+moose_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[61]] <- glm(used ~ tri+sheep_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[62]] <- glm(used ~ tri, family=binomial(logit), data=UsedAvail)
m.biotic[[63]] <- glm(used ~ tri+deer_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[64]] <- glm(used ~ tri+moose_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[65]] <- glm(used ~ tri+sheep_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[66]] <- glm(used ~ tri+goat_w2, family=binomial(logit), data=UsedAvail)
m.biotic[[67]] <- glm(used ~ tri+moose_w2 + sheep_w2, family=binomial(logit), data=UsedAvail)

model.names.biotic <-1:67
aictab(cand.set = m.biotic, modnames = model.names.biotic)

# Top Model (Model 41)
top_model <- glm(used ~ tri+DistFromHighHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
summary(top_model)
car::vif(top_model)

# Second Top Model (Model 59)
second_model <- glm(used ~ tri+deer_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
summary(second_model)
car::vif(second_model)

#           [Caterpillar Plots]                                             ####

# GGally Coefficients
ggcoef(top_model, exclude_intercept = TRUE, exponentiate = T, sort = "ascending")


top_model$terms

#      [ROC Plot]                                                           ####


top_model <- glm(used ~ tri+DistFromHighHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=UsedAvail)
summary(top_model)

par(mfrow = c(2,2))
plot(top_model)

UsedAvail <- na.omit(UsedAvail)
UsedAvail$usedFactor <-as.factor(UsedAvail$used) ## make sure usedFactor is a factor

##### Saving predictions manually, an example with the environment model
UsedAvail$fitted.top_model <- fitted(top_model)
#### this is the predicted probability from the model

UsedAvail$residuals.top_model <- residuals(top_model)
## these are the deviations from the predictions for each row (data point)

UsedAvail$rstudent.top_model <- rstudent(top_model)
## This is a standardized residual - the studentized residual

UsedAvail$hatvalues.top_model <- hatvalues(top_model)
#### this is the first of the leverage statistics, the larger hat value is, the bigger the influence on the fitted value

UsedAvail$cooks.distance.top_model <- cooks.distance(top_model)
#### this is the Cooks leverage statistic, the larger hat value is, the bigger the influence on the fitted value

UsedAvail$obsNumber <- 1:nrow(UsedAvail) ## just added a row number for plotting


# Classifcation table
ppused = UsedAvail$fitted.top_model>0.5
table(ppused,UsedAvail$used)

# ROC Plot
#require(ROCR)
pp = predict(top_model,type="response")
pred = prediction(pp, UsedAvail$used)

perf3 <- performance(pred, "sens", x.measure = "cutoff")
plot(perf3)

perf4 <- performance(pred, "spec", x.measure = "cutoff")
plot(perf4)

perfClass <- performance(pred, "tpr","fpr") # change 2nd and/or 3rd arguments for other metrics
fpr <- perfClass@x.values[[1]]
tpr <- perfClass@y.values[[1]]
sum <- tpr + (1-fpr)
index <- which.max(sum)
cutoff <- perfClass@alpha.values[[1]][[index]]
cutoff

plot(perf3, col="blue") # Sensitivity
plot(perf4, add = TRUE) # Specificity
abline(v=cutoff, col="red") ## optimal cutpoint

BMauc <- performance(pred, measure="auc") 
str(BMauc)
auc <- as.numeric(BMauc@y.values)
auc

plot(perfClass, colorize = T, lwd = 5, print.cutoffs.at=seq(0,1,by=0.1),
     text.adj=c(1.2,1.2),
     main = "ROC Curve")
text(0.5, 0.5, "AUC = 0.878")
abline(v=cutoff, col = "red", lwd = 3)

#      [K Folds]                                                            ####
source("2.R_Code\\1.Scripts\\kxv.R", verbose = FALSE)

# Kfolds with a 'fuzz' factor
kxvPrintFlag=FALSE
kxvPlotFlag=TRUE
kxvFuzzFactor = 0.01
kfolds = kxvglm(top_model$formula, data=UsedAvail, k=5, nbin=10)
kfolds

kxvPrintFlag=FALSE
kxvPlotFlag=TRUE
kxvFuzzFactor = 0.01
kfolds2 = kxvglm(top_model$formula, data=UsedAvail, k=5, nbin=10, partition="pack")
kfolds2

#      [Spatial Predictions]                                                ####
par(mfrow = c(1,1))

top_model$coefficients

model.coefs <- top_model$coefficients[c(1:5)]

rast.top.biotic <- exp(model.coefs[1] + 
                       model.coefs[2]*TRI  + 
                       model.coefs[3]*HighHumanAccessDistance +
                       model.coefs[4]*HSI_deer + 
                       model.coefs[5]*HSI_goat) / (1 +exp(model.coefs[1]+
                                                            model.coefs[2]*TRI  + 
                                                            model.coefs[3]*HighHumanAccessDistance+
                                                            model.coefs[4]*HSI_deer + 
                                                            model.coefs[5]*HSI_goat))

plot(rast.top.biotic)


plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=homerangeALL)
plot(homerangeALL, add=TRUE)
plot(Wolf, col='blue', pch = 16, add=TRUE)

# Bow Valley
bv.raster<-raster()
extent(bv.raster) <- c(xmin=570000, xmax=600000, ymin=5665000, ymax=5685000) 
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=bv.raster)
plot(homerangeALL, add=TRUE)
plot(Wolf, col='blue', pch = 16, add=TRUE)

# Red Deer

rd.raster<-raster()
extent(rd.raster) <- c(xmin=540000, xmax=600000, ymin=5700000, ymax=5730000) 
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=rd.raster)
plot(homerangeALL, add=TRUE)
plot(Wolf, col='blue', pch = 16, add=TRUE)
###############################################################################