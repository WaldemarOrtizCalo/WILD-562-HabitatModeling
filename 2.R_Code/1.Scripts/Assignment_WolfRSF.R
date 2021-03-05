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

# Fit 99% mpc
cp.bow <- mcp(bv, percent=99)
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

# Fit 99% mpc
cp.red <- mcp(rd, percent=99)
plot(rd, col="black")
plot(cp.red[cp.red@data$id=="63",], col="blue",add=TRUE)
plot(cp.red[cp.red@data$id=="87",], col="red",add=TRUE,)
plot(cp.red[cp.red@data$id=="44",], col="green",add=TRUE)
plot(rd, col="black", add=TRUE)

#        [All]                                                              ####

# Fit 99% mpc
cp.all <- mcp(all, percent=99)

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
