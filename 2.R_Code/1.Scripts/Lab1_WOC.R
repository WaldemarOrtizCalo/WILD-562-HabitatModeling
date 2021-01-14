#   Script Details                                                          ####

# Author: Waldemar Ortiz-Calo

# Date:2021-01-14 

# Purpose: Script for Lab 1 of Dr. Mark Hebblewhite's Lab 1 Code 

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

#      Functions                                                            ####

#      Data                                                                 ####
#        [Cougar]                                                           ####

# Data Import
cougar <- read.csv("1.Data\\Lab1csvdata\\golfcourse_cougar.csv", header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
head(cougar)
str(cougar)
table(cougar$Use, cougar$UseNonUse)
cougar$factorUSE <- as.factor(cougar$UseNonUse)

#        [Spatial Data]                                                     ####

# reading in shapefiles (raster package)
elc_habitat<-shapefile("1.Data\\Lab1_data\\elc_habitat.shp")
humanaccess<-shapefile("1.Data\\Lab1_data\\humanacess.shp")
mcp2<-shapefile("1.Data\\Lab1_data\\mcp2.shp")
wolfyht<-shapefile("1.Data\\Lab1_data\\wolfyht.shp")

# reading in raster files (raster package)
deer_w<-raster("1.Data\\Lab1_data\\deer_w2.tif")
moose_w<-raster("1.Data\\Lab1_data\\moose_w2.tif") ## missing moose
elk_w<-raster("1.Data\\Lab1_data\\elk_w2.tif")
sheep_w<-raster("1.Data\\Lab1_data\\sheep_w2.tif") ## missing sheep
goat_w<-raster("1.Data\\Lab1_data\\goat_w2.tif")
wolf_w<-raster("1.Data\\Lab1_data\\wolf_w2.tif")#
elevation2<-raster("1.Data\\Lab1_data\\Elevation2.tif") #resampled
disthumanaccess2<-raster("1.Data\\Lab1_data\\DistFromHumanAccess2.tif") #resampled


###############################################################################
#   Cougar Data Breakdown                                                   ####
#      Summary Data                                                         ####

## Cougar summary statistics
summary(cougar)

cougar_df <- as_tibble(cougar)
byUse <- group_by(cougar_df, UseNonUse)
summarise(byUse, slope = mean(Slope))
summarise(byUse, DistTrails = mean(AllTrails))
summarise(byUse, DistCover = mean(CoverDist))
summarise(byUse, DistRoads = mean(Roads))

## Cougar graphing
par(mfrow = c(2,2))
multhist(list(cougar$AllTrails[cougar$factorUSE==1],cougar$AllTrails[cougar$factorUSE==0]), freq = TRUE, main = "Trails")
# I chose to put a legend in the upper right hand graph. 
# That's what the additional arguments in the line below specify.
multhist(list(cougar$CoverDist[cougar$factorUSE==1],cougar$CoverDist[cougar$factorUSE==0]), freq = TRUE, main = "Cover Distance", legend.text = c("Used", "Unused"), args.legend = list(bty = "n"))
multhist(list(cougar$Roads[cougar$factorUSE==1],cougar$Roads[cougar$factorUSE==0]), freq = TRUE, main = "Roads")
multhist(list(cougar$Slope[cougar$factorUSE==1],cougar$Slope[cougar$factorUSE==0]), freq = TRUE, main = "Slope")

### Cougar Boxplots

par(mfrow= c(2,2))
boxplot(AllTrails~factorUSE, ylab="Distance (m)", xlab="Used",main = "Trails", data=cougar)
boxplot(CoverDist~factorUSE, ylab="Distance (m)", xlab="Used", main = "Cover", data=cougar)
boxplot(Roads~factorUSE, ylab="Distance (m)", xlab="Used",main = "Roads", data=cougar)
boxplot(Slope~factorUSE, ylab="Slope", xlab="Used", main = "Slope", data=cougar)

#      Statistical Exploration                                              ####

## Cougar statistical tests

t.test(AllTrails~factorUSE, alternative='two.sided', conf.level=.95, 
       var.equal=FALSE, data=cougar)
t.test(CoverDist~factorUSE, alternative='two.sided', conf.level=.95, 
       var.equal=FALSE, data=cougar)
t.test(Roads~factorUSE, alternative='two.sided', conf.level=.95, 
       var.equal=FALSE, data=cougar)
t.test(Slope~factorUSE, alternative='two.sided', conf.level=.95, 
       var.equal=FALSE, data=cougar)

## Cougar statistical analyses - linear models
USEonly <- subset(cougar, subset=UseNonUse == 1)
# Construct linear models of track counts as a function of 4 covariates
# Distance to trails model
trails <- glm(Use ~ AllTrails, family=gaussian(identity), data=USEonly)
summary(trails)
# Slope model
slope <- glm(Use ~ Slope, family=gaussian(identity), data=USEonly)
summary(slope)
# Distance to cover model
cover <- glm(Use ~ CoverDist, family=gaussian(identity), data=USEonly)
summary(cover)
# Distance to roads model
roads <- glm(Use ~ Roads, family=gaussian(identity), data=USEonly)
summary(roads)

# Visualize with a histogram
par(mfrow= c(1,1))
hist(USEonly$Use, scale="frequency", breaks="Sturges", col="darkgray")
## not really that normal
shapiro.test(USEonly$Use)
# Visualize with a histogram

USEonly$lnUSE <- with(USEonly, log(Use))
hist(USEonly$lnUSE) ## a bit more normal

### Now re-fit the models
# Distance to trails model
ln.trails <- glm(lnUSE ~ AllTrails, family=gaussian(identity), data=USEonly)
summary(ln.trails)
# Slope model
ln.slope <- glm(lnUSE ~ Slope, family=gaussian(identity), data=USEonly)
summary(ln.slope)
# Distance to cover model
ln.cover <- glm(lnUSE ~ CoverDist, family=gaussian(identity), data=USEonly)
summary(ln.cover)
# Distance to roads model
ln.roads <- glm(lnUSE ~ Roads, family=gaussian(identity), data=USEonly)
summary(ln.roads)


shapiro.test(cougar$Use)
hist(cougar$Use, scale="frequency", breaks="Sturges", col="darkgray")

# Trails model
logitTrails <- glm(UseNonUse ~ AllTrails, family=binomial(logit), data=cougar)
summary(logitTrails)
# Slope model
logitSlope <- glm(UseNonUse ~ Slope, family=binomial(logit),   data=cougar)
summary(logitSlope)
# Cover model
logitCover <- glm(UseNonUse~ CoverDist, family=binomial(logit),   data=cougar)
summary(logitCover)
# Roads model
logitRoads <- glm(UseNonUse ~ Roads, family=binomial(logit),  data=cougar)
summary(logitRoads)

## Logistic Regression Visualization - graphing
ggplot(cougar, aes(x=Slope, y=UseNonUse)) + geom_rug() + stat_smooth(method="glm", method.args=list(family="binomial"))
ggplot(cougar, aes(x=CoverDist, y=UseNonUse)) + geom_rug() + stat_smooth(method="glm", method.args=list(family="binomial"))
#   Spatial Data                                                            ####
#        [Shapefile]                                                        ####
# make a very basic plot of shapefile after resetting graphical parameters
par(mfrow= c(1,1))
plot(elc_habitat)
plot(wolfyht)

# look at the class of the shapefile
class(elc_habitat)

# look at structure of shapefile
#str(elc_habitat)

# look at first 20 rows of data for shapefile
head(elc_habitat@data, n=20)

# look at the projection of the shapefile (note the use of "@" instead of "$")
elc_habitat@proj4string@projargs
wolfyht@proj4string@projargs

elc_habitat <- spTransform(elc_habitat, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

# check new projection in geographic coordinate system WGS84
elc_habitat@proj4string@projargs

# reset projection back to what is was previously
elc_habitat <-  spTransform(elc_habitat, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

# look at the projection of the shapefile 
elc_habitat@proj4string@projargs

# another way to change it back to previous 
elc_habitat <-  spTransform(elc_habitat, CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#write shapefile to files
writeOGR(elc_habitat,".","elc_habitat_NEW",driver="ESRI Shapefile", overwrite_layer=TRUE)

## Extents
extent(elc_habitat)
extent(wolfyht)

#        [Raster]                                                           ####

## reset graphical parameters
par(mfrow= c(1,1)) 

# make a very basic plot of raster
plot(deer_w)

# look at the class of the raster
class(deer_w)

# look at basic raster summary
deer_w

# look at raster data
deer_w@data

# look at structure of raster
str(deer_w)

# look at the projection of the raster (note the use of "@" instead of "$")
deer_w@crs@projargs

# look at the spatial extent of the raster
extent(deer_w)

# look at the resolution of the raster
res(deer_w)

#And if we need to change the projection of the raster to another projection (package rgdal) - note for the Tutorial I've masked these out # because they take a LONG time. 

#deer_w2 <- projectRaster(deer_w, crs="+init=epsg:4326")

## or try this the same way we did above

#deer_w2 <- projectRaster(deer_w, crs="+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

#check projection of the raster
deer_w2@crs@projargs

#change it back to what it was using another raster layers projection
#deer_w <- projectRaster(deer_w2, wolf_w)

#check projection of the raster
deer_w@crs@projargs

## Creating a Raster Stack

#One of the GREAT things about doing GIS in R is the ability to create a consistent raster stack of rasters of the same projection, extent, etc.  LEts create a raster stack!
  
all.rasters<-stack(deer_w, moose_w, elk_w, sheep_w, goat_w, elevation2, disthumanaccess2)

plot(all.rasters)

#check class
class(all.rasters)

#learn more about 
#?writeRaster()
writeRaster(deer_w, "1.Data\\Lab1_data\\new_deer.tiff", "GTiff", overwrite=TRUE)

# Mapping with Mapview

mapView(wolfyht2, zcol="NAME", native.crs = FALSE, legend = TRUE, cex=5, lwd=2, map.type = "Esri.DeLorme")

mapView(wolfyht, zcol="NAME", legend = TRUE, cex=5, lwd=2, map.type = "Esri.WorldImagery")

## Open Topo Map
mapView(wolfyht, zcol="NAME", legend = TRUE, cex=5, lwd=2, map.type = "OpenTopoMap")

## Esri World Topo Map
mapView(wolfyht, zcol="NAME", legend = TRUE, cex=5, lwd=2, map.type = "Esri.WorldTopoMap")



