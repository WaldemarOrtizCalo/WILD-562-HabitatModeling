#   Script Details                                                          ####

# Author: Waldemar Ortiz-Calo

# Date:2021-01-21 

# Purpose: 

###############################################################################
#   Library / Functions / Data                                              ####

#      Library                                                              ####
library(ks)
library(plotrix)
library(lattice)
library(adehabitatHR)
library(maptools)
library(mapview)
library(rgdal)
library(sp)
library(raster)
library(ggplot2)
library(colorRamps)
library(rgeos)
library(fasterize)

#      Functions                                                            ####
#      Data                                                                 ####
#        [OBJECTIVE 1]                                                      ####

elc_habitat<-shapefile("1.Data/Lab2_data/elc_habitat.shp")
humanaccess<-shapefile("1.Data/Lab2_data/humanacess.shp")

wolfyht<-shapefile("1.Data/Lab2_data/wolfyht.shp")

# Biological Rasters
deer_w<-raster("1.Data/Lab2_data/deer_w2.tif")
moose_w<-raster("1.Data/Lab2_data/moose_w2.tif")
elk_w<-raster("1.Data/Lab2_data/elk_w2.tif")
sheep_w<-raster("1.Data/Lab2_data/sheep_w2.tif")
goat_w<-raster("1.Data/Lab2_data/goat_w2.tif")
wolf_w<-raster("1.Data/Lab2_data/wolf_w2.tif")
elevation2<-raster("1.Data/Lab2_data/Elevation2.tif")

# Distance to Human Layer
disthumanaccess2<-raster("DistFromHumanAccess2.tif")

# Distance to High Human Access
disthighhumanaccess<-raster("1.Data/Lab2_data/DistFromHighHumanAccess2.tif")

#        [OBJECTIVE 2]                                                      ####
###############################################################################
#   OBJECTIVE 1 - Managing Spatial Raster and Shapefile data                ####
#      [Data Exploration]                                                   ####
 
# Plots

plot(elc_habitat)
plot(humanaccess)



# Note that there are two fields, Easting and Northing which are the X and Y coordinates in UTM zone 11.  We will use these to map it for each PackID
# base plot of wolf packs by color with legend

plot(wolfyht@data$EASTING,wolfyht@data$NORTHING, col=c("red","blue")[wolfyht@data$PackID],ylab="Northing",xlab="Easting")
legend(555000,5742500,unique(wolfyht@data$Pack),col=c("blue","red"),pch=1)

### Habitat Suitability Index Models
elc_habitat@data$id <- rownames(elc_habitat@data)
elc_habitatPoly <- fortify(elc_habitat, region = "id")
elc_habitatDF <- merge(elc_habitatPoly, elc_habitat@data, by = "id")

#construct ggplot2 plot for Moose Winter Habitat
elk_plot<-ggplot(elc_habitatDF, aes(x = long, y = lat,group=group, fill = as.factor(MOOSE_W))) + 
  geom_polygon() + labs(x="Easting",y="Northing") + theme(axis.text.y = element_text(angle = 90, hjust=0.5))

#adjust fill colors of MOOSE_W  (note that I just selected some random colors, but made "7" as blue)
elk_plot2 <- elk_plot + scale_fill_manual(name="MOOSE_W",values=c("gray","gray", "red", "orange", "yellow", "green","darkblue"))
elk_plot2

#      [Bighorn Sheep Winter Habitat Model]                                 ####
#        [Creating Prey Raster Layers]                                      ####

#construct ggplot2 plot for Moose Winter Habitat
sheep_plot<-ggplot(elc_habitatDF, aes(x = long, y = lat,group=group, fill = as.factor(SHEEP_W))) + 
  geom_polygon() + labs(x="Easting",y="Northing") + theme(axis.text.y = element_text(angle = 90, hjust=0.5))

#adjust fill colors of MOOSE_W  (note that I just selected some random colors, but made "7" as blue)
sheep_plot2 <- sheep_plot + scale_fill_manual(name="SHEEP_W",
                                              values=c("gray","gray", "red", "orange", "yellow", "green","darkblue"))
sheep_plot2

## Part 1b - Spatial Raster Operations

# Create a mask raster to use as a template for converting shapefile data to rasters
#create an empty raster
mask.raster <- raster()

#set extent (note that I customized this extent so it covered both elc_habitat and humanacess)
extent(mask.raster) <- c(xmin=443680.6, xmax=650430.4, ymin=5618405, ymax=5789236) 	

#set the resolution to 30 m 
res(mask.raster)<-30

#match projection to elc_habitat shapefile
projection(mask.raster)<- "+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

#set all values of mask.raster to zero
mask.raster[]<-0

# create rasters for DEER_W, MOOSE_W, ELK_W, SHEEP_W, GOAT_W, WOLF_W
## ******  note for today's Lab 2 we have masked out these steps as they take a LONG time. ******
#deer_w<-rasterize(elc_habitat, mask.raster, field="DEER_W")
#moose_w<-rasterize(elc_habitat, mask.raster, field="MOOSE_W")
#elk_w<-rasterize(elc_habitat, mask.raster, field="ELK_W")
#sheep_w<-rasterize(elc_habitat, mask.raster, field="SHEEP_W")
#goat_w<-rasterize(elc_habitat, mask.raster, field="GOAT_W")
#wolf_w<-rasterize(elc_habitat, mask.raster, field="WOLF_W")

#plot result
#plot(wolf_w)

## Fasterize

deer_w <- fasterize(elc_habitat_sf, mask.raster, field = "DEER_W")
# moose_w <- fasterize(elc_habitat_sf, mask.raster, field = "MOOSE_W")
# elk_w <- fasterize(elc_habitat_sf, mask.raster, field = "ELK_W")
# sheep_w <- fasterize(elc_habitat_sf, mask.raster, field = "SHEEP_W")
# goat_w <- fasterize(elc_habitat_sf, mask.raster, field = "GOAT_W")
# wolf_w <- fasterize(elc_habitat_sf, mask.raster, field = "WOLF_W")
# 
# plot(wolf_w)

#resample elevation and humanaccess to match mask.raster
#elevation2<-resample(elevation, mask.raster, method="bilinear")
#disthumaccess2<-resample(disthumaccess, mask.raster, method="bilinear")

#write raster layers to file
#writeRaster(deer_w, "deer_w2.tiff", "GTiff")
#writeRaster(moose_w, "moose_w2.tiff", "GTiff")
#writeRaster(elk_w, "elk_w2.tiff", "GTiff")
#writeRaster(sheep_w, "sheep_w2.tiff", "GTiff")
#writeRaster(goat_w, "goat_w2.tiff", "GTiff")
#writeRaster(wolf_w, "wolf_w2.tiff", "GTiff")
#writeRaster(elevation2, "Elevation2.tiff", "GTiff")
#writeRaster(disthumaccess2, "DistFromHumanAccess2.tiff", "GTiff")


mapview(wolfyht) + deer_w + sheep_w
#        [Distance to Human Layer]                                          ####

## Part 1b - Distance to Human Access Layer _ 

#first create an empty raster
dist.raster <- raster()

#set extent 
extent(dist.raster) <- extent(humanaccess)

#set the resolution to 30 m (Note that this takes a very long time with a 30 m resolution-even on my machine)
res(dist.raster)<-30

#match projection to humanaccess shapefile
projection(dist.raster)<- "+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

#set all values of dist.raster to zero
dist.raster[]<-0

#now rasterize the humanaccess layer and set human features (e.g., roads, trails) to 1
#human.raster<-rasterize(humanaccess, dist.raster, 1)

#calculate distance to human access- NOTE : DO NOT RUN THIS IT TAKES FOREVER; I ENDED UP
#DOING THIS CALCULATION JUST FOR THE POINTS
#accessdist <-system.time(distance(human.raster))

#write raster to file
#writeRaster(accessdist, "DistFromHumanAccess.tiff", "GTiff")
plot(disthumanaccess2)

#        [High Human Access]                                                ####

## Advanced question: how would you calculate distance to high human use?

#first reclassify labels on humanaccess.shp file so they are correct (note: need to bring in humanaccess.shp above)
levels(as.factor(humanaccess$SUM_CLASS))

#[1] "0"         "High"      "HIGH"      "Low"       "LOW"      
#[6] "MEDIUM"    "Moderate"  "Nil"       "NIL"       "VERY HIGH"

#convert humanaccess$SUM_CLASS to a factor
humanaccess$SUM_CLASS<-as.factor(humanaccess$SUM_CLASS)

levels(humanaccess$SUM_CLASS)[1]<-"NIL"
levels(humanaccess$SUM_CLASS)[2]<-"HIGH"
levels(humanaccess$SUM_CLASS)[3]<-"LOW"
levels(humanaccess$SUM_CLASS)[5]<-"MODERATE"
levels(humanaccess$SUM_CLASS)[6]<-"NIL"

#create indicator variable for high or not high human access
highaccess<-humanaccess[humanaccess@data$SUM_CLASS=="HIGH" | humanaccess@data$SUM_CLASS=="VERY HIGH", ]

plot(humanaccess)
plot(highaccess, col="red", add=TRUE)

# Load previously developed layer. 

plot(disthighhumanaccess)

###############################################################################
#   OBJECTIVE 2 - Home Range Analysis                                       ####

rd.data<-wolfyht[wolfyht@data$Pack=="Red Deer",]
x<-rd.data@data$EASTING
y<-rd.data@data$NORTHING
xy<-cbind(x,y)
class(xy)
rd <- data.frame(as.character(rd.data@data$NAME))
coordinates(rd) <- xy
proj4string(rd) <-  CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

class(rd)

# Fit 99% mpc
#cp.rd <- mcp(rd, percent=99)
#note error that one animal does not have at least 5 locations
table(rd.data$NAME)
#42 60 69 70 81 82 84 
#43 25  4 15  3  2  1  
#looks like 4 of the wolves do not have enough locations

#remove these individuals with too few of locations
names(rd)<-"NAME"
rd<-rd[rd@data$NAME!="69" & rd@data$NAME!="81" & rd@data$NAME!="82" & rd@data$NAME!="84",]

#remove unused NAME levels
rd@data$NAME<-factor(rd@data$NAME)


## FIt 99% Minimum Convex Polygon to these cleaned data
# Fit 99% mpc
cp.rd <- mcp(rd, percent=99)
plot(rd, col="black")
plot(cp.rd[cp.rd@data$id=="42",], col="blue", add=TRUE)
plot(cp.rd[cp.rd@data$id=="70",], col="green", add=TRUE)
plot(cp.rd[cp.rd@data$id=="60",], col="red", add=TRUE)
plot(rd, col="black", add=TRUE)

#check area in square meters for each Red Deer wolf pack
as.data.frame(cp.rd)

#calculate area for different percents of MPC
mcp.area(rd, percent=seq(50, 100, by=5))


##### KErnel density estimates

#calculate 99% KDE for Red Deer wolf pack
red.deerUD <- kernelUD(rd, grid=30, extent=0.5, same4all=TRUE) # reference grid
image(red.deerUD)

#get polygons for home ranges
homerangeRD <- getverticeshr(red.deerUD)
as.data.frame(homerangeRD)
class(homerangeRD)
plot(homerangeRD, col=2:4)

#Estimate UD in raster mode
red.deerud <- getvolumeUD(red.deerUD) 
red.deerud

## Set up graphical parameters for the output of getvolumeUD 
par(mar=c(0,0,2,0)) #set margin
image(red.deerud[[1]]) #for first wolf only
title("Red Deer Wolf UD") 
xyzv <- as.image.SpatialGridDataFrame(red.deerud[[1]]) 
contour(xyzv, add=TRUE)

## Example of doing for 1 wolf  
fud <- red.deerud[[1]] #for first wolf only
## store the value of the volume under UD in a vector hr95 
hr95 <- as.data.frame(fud)[,1] 
## if hr95 is <= 95 then the pixel belongs to the home range
## (takes the value 1, 0 otherwise)
hr95 <- as.numeric(hr95 <= 95) 
## Converts into a data frame 
hr95 <- data.frame(hr95) 
## Converts to a SpatialPixelsDataFrame 
coordinates(hr95) <- coordinates(red.deerud[[1]])
gridded(hr95) <- TRUE 
## display the results 
image(hr95)

## Bow Valley wolves 99% MCP

#first convert the spatialpointsdataframe to spatial points object
bv.data<-wolfyht[wolfyht@data$Pack=="Bow valley",]
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

#Red Deer wolf pack
#check area for each Red Deer wolf pack
as.data.frame(cp.bow)

#calculate area for different percents of MPC
mcp.area(bv, percent=seq(50, 100, by=5))

#calculate 99% KDE for Red Deer wolf pack
bow.valleyUD <- kernelUD(bv, grid=30, extent=0.1, same4all=TRUE) # reference grid
image(bow.valleyUD)

#get polygons for home ranges
homerangeBV <- getverticeshr(bow.valleyUD)
as.data.frame(homerangeBV)
class(homerangeBV)
plot(homerangeBV, col=2:4)

#Estimate UD in raster mode
bow.valleyud <- getvolumeUD(bow.valleyUD) 
bow.valleyud

## Set up graphical parameters for the output of getvolumeUD 
par(mar=c(0,0,2,0)) #set margin
image(bow.valleyud[[1]])
title("Bow Valley Pack UD") 
xyzv <- as.image.SpatialGridDataFrame(bow.valleyud[[1]]) 
contour(xyzv, add=TRUE)


## Bow Valley UD

fud <- bow.valleyud[[1]]
## store the value of the volume under 95% UD in a vector hr95 
hr95 <- as.data.frame(fud)[,1] 
## if hr95 is <= 95 then the pixel belongs to the home range
## (takes the value 1, 0 otherwise)
hr95 <- as.numeric(hr95 <= 95) 
## Converts into a data frame 
hr95 <- data.frame(hr95) 
## Converts to a SpatialPixelsDataFrame 
coordinates(hr95) <- coordinates(bow.valleyud[[1]])
gridded(hr95) <- TRUE 
## display the results 
image(hr95)


## Calculate 99% minimum convex polygon for both wolf packs

#first convert the spatialpointsdataframe to spatial points object
x<-wolfyht@data$EASTING
y<-wolfyht@data$NORTHING
xy<-cbind(x,y)

all <- data.frame(as.character(wolfyht@data$Pack))
coordinates(all) <- xy
proj4string(all) <-  CRS("+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Fit 99% mpc
cp.all <- mcp(all, percent=99)

plot(wolfyht, col="black")
plot(cp.all[cp.all@data$id=="Bow valley",], col="blue",add=TRUE)
plot(cp.all[cp.all@data$id=="Red Deer",], col="green",add=TRUE)
plot(wolfyht, col="black", add=TRUE)

#Checking areas, etc 
#check area for each Red Deer wolf pack
as.data.frame(cp.all)

#calculate area for different percents of MPC
mcp.area(all, percent=seq(50, 100, by=5))

#calculate 99% KDE for both wolf packs
allUD <- kernelUD(all, grid=30, extent=0.5, same4all=TRUE) # reference grid
image(allUD)

#Get polygons for home ranges
homerangeALL <- getverticeshr(allUD)
as.data.frame(homerangeALL)
class(homerangeALL)
plot(homerangeALL, col=2:3)

#Estimate UD in raster mode
allud <- getvolumeUD(allUD) 
allud

## Set up graphical parameters for the output of getvolumeUD 
par(mar=c(0,0,2,0)) #set margin
image(allud[[1]]) #for first wolf only
title("Output of getvolumeUD") 
xyzv <- as.image.SpatialGridDataFrame(allud[[1]]) 
contour(xyzv, add=TRUE)

## store the volume under the UD (as computed by getvolumeUD) 
## of the first animal in fud 

fud <- allud[[1]] #for first wolf pack only
## store the value of the volume under UD in a vector hr95 
hr95 <- as.data.frame(fud)[,1] 
## if hr95 is <= 95 then the pixel belongs to the home range
## (takes the value 1, 0 otherwise)
hr95 <- as.numeric(hr95 <= 95) 
## Converts into a data frame 
hr95 <- data.frame(hr95) 
## Converts to a SpatialPixelsDataFrame 
coordinates(hr95) <- coordinates(allud[[1]])
gridded(hr95) <- TRUE 
## display the results 
image(hr95)

###############################################################################
#   OBJECTIVE 3 - Learn How to Sample Availability Within Home Ranges       ####


#subset polygons by wolf pack
red.deerPOLY<-homerangeALL[homerangeALL@data$id=="Red Deer",]
bow.valleyPOLY<-homerangeALL[homerangeALL@data$id=="Bow valley",]

#generate 1000 points from Red Deer wolf pack KDE polygon
rd.avail<-spsample(red.deerPOLY, 1000, "random")
plot(rd.avail)

#generate 1000 points from Bow valley wolf pack KDE polygon
bv.avail<-spsample(bow.valleyPOLY, 1000, "random")
plot(bv.avail)

# lets plot them all together, used and home-range level availability
plot(wolfyht@data$EASTING,wolfyht@data$NORTHING, col=c("red","blue")[wolfyht@data$PackID],ylab="Northing",xlab="Easting")
legend(555000,5742500,unique(wolfyht@data$Pack),col=c("blue","red"),pch=1)
plot(bv.avail, add=TRUE)
plot(rd.avail, add=TRUE)

###############################################################################
#   OBJECTIVE 4 - Extracting GIS covariates For Points                      ####
## 

#IF you have tidyverse loaded, remove it using this command:
rs.unloadPackage("tidyr")
#Note today we are using tidyverse, which depends on tidyr, but I've loaded it separately later in summary statistics below. 
all_rasters<-stack(deer_w, moose_w, elk_w, sheep_w, goat_w, wolf_w,elevation2, disthumanaccess2, disthighhumanaccess)
class(all_rasters)

#Extract covariate values for Red Deer wolf data  
cov.outRD<-extract(all_rasters, rd.data)
head(cov.outRD)

#Extract covariate values for available points
cov.availRD<-extract(all_rasters, rd.avail)

#Extract covariate values for Bow valley wolf data  
cov.outBV<-extract(all_rasters, bv.data)

#Extract covariate values for available points
cov.availBV<-extract(all_rasters, bv.avail)

###############################################################################
#   OBJECTIVE 5 – Exploratory analyses of wolf habitat use with R           ####

rdused <- as.data.frame(cov.outRD)
rdused$pack <- c("Red Deer")

## repeat for Bow Valley pack
bvused <- as.data.frame(cov.outBV)
bvused$pack <- c("Bow Valley")

wolfused <- merge(rdused, bvused, all.x= TRUE, all.y = TRUE)
str(wolfused)
head(wolfused)

## and for next week, lets add a new column for a 1=used 0 = avail
wolfused$used <- 1

### Numerical summaries (review from Lab 1)

summary(wolfused)

#Note that there are 4 missing NA values in the distance to human access layer. Why?

plot(disthumanaccess2)
plot(wolfyht, add = TRUE)

#Note, this replaces just NAs for just the wolf used locations with NA for distance to human access, but note there were 13 missing values, say, for some of the ELC H.S.I models. 
wolfused <- na.omit(wolfused)
summary(wolfused)
dim(wolfused)

### Graphical data visualization
## ungulate HSI models first
par(mfrow = c(2,3))
hist(wolfused$deer_w2)
hist(wolfused$elk_w2)
hist(wolfused$moose_w2)
hist(wolfused$sheep_w2)
hist(wolfused$goat_w2)

#Continuous covariates next, elevation (m), distance (m). 

par(mfrow = c(3,1))
hist(wolfused$Elevation2)
hist(wolfused$DistFromHumanAccess2)
hist(wolfused$DistFromHighHumanAccess2)

### Data exploration by wolf pack

par(mfrow = c(1,1))
# Plot Bow Valley
hist(wolfused$Elevation2[wolfused$pack=="Bow Valley"],breaks=50, xlim = c(1400,2250), probability = TRUE, main="Wolf Habitat Selection", xlab="Elevation") 

#Plot Red Deer
hist(wolfused$Elevation2[wolfused$pack=="Red Deer"],breaks=50, col="darkgray",probability =TRUE, add=TRUE)
# Add legend
legend("topright", c("Bow Valley", "Red Deer"), fill = c("white","darkgray"),border = "black")

# So, the Red Deer wolf pack 'uses' higher elevations than the Bow Valley wolf pack. 

# Now, repeat these steps on your own for all other covariates using the plotrix R package. 

par(mfrow = c(2,1))
multhist(list(wolfused$Elevation2[wolfused$pack=="Bow Valley"],wolfused$Elevation2[wolfused$pack=="Red Deer"]), freq = TRUE, main = "Elevation")
# I chose to put a legend in the lower right hand graph. 
# That's what the additional arguments in the line below specify.
multhist(list(wolfused$DistFromHumanAccess2[wolfused$pack=="Bow Valley"],wolfused$DistFromHumanAccess2[wolfused$pack=="Red Deer"]), freq = TRUE, main = "Distance From Humans", legend.text = c("Bow Valley", "Red Deer"), args.legend = list(bty = "n"))

## Ungulate HSI's 

par(mfrow = c(2,3))
multhist(list(wolfused$elk_w2[wolfused$pack=="Bow Valley"],wolfused$elk_w2[wolfused$pack=="Red Deer"]), freq = TRUE, main = "Elk HSI", legend.text = c("Bow Valley", "Red Deer"), args.legend = list(bty = "n"))
multhist(list(wolfused$deer_w[wolfused$pack=="Bow Valley"],wolfused$deer_w[wolfused$pack=="Red Deer"]), freq = TRUE, main = "Deer HSI", legend.text = c("Bow Valley", "Red Deer"), args.legend = list(bty = "n"))
multhist(list(wolfused$moose_w2[wolfused$pack=="Bow Valley"],wolfused$moose_w2[wolfused$pack=="Red Deer"]), freq = TRUE, main = "Moose HSI", legend.text = c("Bow Valley", "Red Deer"), args.legend = list(bty = "n"))
multhist(list(wolfused$sheep_w2[wolfused$pack=="Bow Valley"],wolfused$sheep_w2[wolfused$pack=="Red Deer"]), freq = TRUE, main = "Sheep HSI", legend.text = c("Bow Valley", "Red Deer"), args.legend = list(bty = "n"))
multhist(list(wolfused$goat_w2[wolfused$pack=="Bow Valley"],wolfused$goat_w2[wolfused$pack=="Red Deer"]), freq = TRUE, main = "Goat HSI", legend.text = c("Bow Valley", "Red Deer"), args.legend = list(bty = "n"))

## Boxplots (from the lattice package)
bwplot(elk_w2 + deer_w2+moose_w2+ sheep_w2+goat_w2~pack, auto.key=TRUE,allow.multiple = TRUE,data=wolfused, outer=TRUE)

bwplot(DistFromHumanAccess2 + DistFromHighHumanAccess2 + Elevation2~pack, auto.key=TRUE,allow.multiple = TRUE,data=wolfused, outer=TRUE)


## Summary Statistics

aggregate(Elevation2 ~ pack, data=wolfused, FUN=mean)
aggregate(DistFromHumanAccess2 ~ pack, data=wolfused, FUN=mean)
aggregate(.~pack, data=wolfused[c("pack", "deer_w2", "elk_w2", "moose_w2", "sheep_w2", "goat_w2")], mean)

sapply(wolfused, median, na.rm=TRUE)


## Summary Statistics with Tibble
library(tidyverse)
wolf_df <- as_tibble(wolfused)

wolf_df %>% group_by(pack) %>% summarise(mean(Elevation2))
wolf_df %>% group_by(pack) %>% summarise(mean(DistFromHumanAccess2))
wolf_df %>% group_by(pack) %>% summarise(mean(DistFromHighHumanAccess2))

wolf_df %>% group_by(pack) %>% summarise(mean(moose_w2))
wolf_df %>% group_by(pack) %>% summarise(mean(elk_w2))
wolf_df %>% group_by(pack) %>% summarise(mean(sheep_w2))
wolf_df %>% group_by(pack) %>% summarise(mean(deer_w2))
wolf_df %>% group_by(pack) %>% summarise(mean(goat_w2))
wolf_df %>% group_by(pack) %>% summarise(mean(wolf_w2))



###############################################################################
