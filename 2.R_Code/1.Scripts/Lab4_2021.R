# title: 'WILD 562 Lab4 : Categorical Covariates'
# author: "Mark Hebblewhite"

setwd("/Users/mark.hebblewhite/Box Sync/Teaching/UofMcourses/WILD562/Spring2021/Labs/lab4")

# Lab 4: Categorical Resource Selection
## Preliminaries - Loading Packages
#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("tidyverse", "adehabitatHS", "adehabitatHR", "mapview", "rgdal", "sp", "raster","ggplot2","colorRamps","rgeos")

#run function to install packages
ipak(packages)
#############################################################
# Importing Landcover Map from an ArcGIS 
## Loading and Manipulating Landcover Data
landcover<-raster("Data/landcover/landcover")
image(landcover, col=rainbow(16))
landcover
#str(landcover)
landcover@data@attributes

landcover@crs@projargs
extent(landcover)

writeRaster(landcover, "Data/landcover16.tif", "GTiff", overwrite = TRUE)
landcover16 <- raster("Data/landcover16.tif") # bringing it back in
res(landcover16)
extent(landcover16)
crs(landcover16)
str(landcover16@data@attributes)
plot(landcover16)

data1 <- data.frame(landcover@data@attributes[[1]][2:3])
names(landcover16@data@attributes[[1]])[2]<-"COUNT"
landcover16@data@attributes <- merge(landcover16@data@attributes, data1, by="COUNT")
landcover16@data@attributes<-landcover16@data@attributes[c(2,1,3)]
str(landcover16@data@attributes)
landcover16@data@attributes

#############################################################
##### Bring in Wolf Data  ###################################

wolfyht<-shapefile("Data/wolfyht.shp")
plot(landcover16, col=rainbow(16))
plot(wolfyht, add=TRUE, type="p", color = "gray25", pch=19, cex = 0.75)
mapview(landcover16, zcol = "HABITATTYPE") + wolfyht

## lets make a second plot zooming into a specific area of the Red Deer pack
yht.raster <- raster()
extent(yht.raster) <- c(xmin=570000, xmax=600000, ymin=5720000, ymax=5740000) 	
plot(landcover16, col=rainbow(16), ext=yht.raster)
legend("topleft", legend = c("Open Conifer", "Mod. Conifer", "Closed Conifer", "Deciduous", "Mixed", "Regen", "Herb", "Shrub", "Water", "Rock-Ice", "Cloud", "Burn-Forest", "Burn-Grassland", "Burn-Shrub", "Alpine Herb", "Alpine Shrub"), fill = rainbow(16), cex=0.75)
plot(wolfyht, add=TRUE, type="p", color = "gray25", pch=19)

#############################################################
# Loading wolfkde dataframe from last lab with landcover extracted
#######

wolfkde <- read.csv("Data/wolfkde.csv")
table(wolfkde$used, wolfkde$pack)
summary(wolfkde)
#wolfkde <- na.omit(wolfkde)
#summary(wolfkde)
#table(wolfkde$used, wolfkde$pack)

ggplot(wolfkde, aes(x=EASTING, y = NORTHING, color=usedFactor)) + geom_point() + stat_density2d() + facet_grid(pack ~ ., scales="free")

# or, Facetting by Used
ggplot(wolfkde, aes(x=EASTING, y = NORTHING)) + geom_point() + stat_density2d() + facet_grid(pack ~ usedFactor, scales="free")

#############################################################
# Univariate Model-fitting 
#############################################################


### First for all packs
elev <- glm(used ~ Elevation2, family=binomial(logit), data=wolfkde)
distacc <- glm(used ~ DistFromHumanAccess2, family=binomial(logit), data=wolfkde)
disthha <- glm(used ~ DistFromHighHumanAccess2, family=binomial(logit), data=wolfkde)
sheep <- glm(used ~ sheep_w2, family=binomial(logit), data=wolfkde)
goat <- glm(used ~ goat_w2, family=binomial(logit), data=wolfkde)
elk <- glm(used ~ elk_w2, family=binomial(logit), data=wolfkde)
moose <- glm(used ~ moose_w2, family=binomial(logit), data=wolfkde)
deer <- glm(used ~ deer_w2, family=binomial(logit), data=wolfkde)

#Creating tables of B, SE. First grab all of the estimates and standard errors

models = rbind(summary(elev)$coefficients[,1:2], summary(disthha)$coefficients[,1:2], summary(distacc)$coefficients[,1:2], summary(sheep)$coefficients[,1:2], summary(goat)$coefficients[,1:2], summary(elk)$coefficients[,1:2], summary(moose)$coefficients[,1:2], summary(deer)$coefficients[,1:2])
# Name your models
modelnames = c("elev","disthha", "distacc", "sheep", "goat", "elk", "moose", "deer")
# Now put all of your estimates in a pretty table with names that you'll remember!
  estimates.all = matrix(models, nrow=2*length(modelnames), ncol=2, dimnames = list(paste(rep(modelnames, each=2),c("intercept", "coefficient")), c("B", "SE")))
estimates.all
plot(estimates.all)

#############################################################
# Categorical Resource Selection
#############################################################

levels(wolfkde$landcover16) ## see, all we have is landcover code

wolfkde$habitatType = ifelse(wolfkde$landcover16 == 0, "NA", 
                             ifelse(wolfkde$landcover16 == 1, "Open Conifer", 
                                    ifelse(wolfkde$landcover16 == 2, "Moderate Conifer", 
                                           ifelse(wolfkde$landcover16 == 3, "Closed Conifer", 
                                                  ifelse(wolfkde$landcover16 == 4, "Deciduous", 
                                                         ifelse(wolfkde$landcover16 == 5, "Mixed", 
                                                                ifelse(wolfkde$landcover16 == 6, "Regen", 
                                                                       ifelse(wolfkde$landcover16 == 7, "Herbaceous",                 
                                                                              ifelse(wolfkde$landcover16 == 8, "Shrub",                       
                                                                                     ifelse(wolfkde$landcover16 == 9, "Water", 
                                                                                            ifelse(wolfkde$landcover16 == 10, "Rock-Ice", 
                                                                                                   ifelse(wolfkde$landcover16 == 11, "Cloud", 
                                                                                                          ifelse(wolfkde$landcover16 == 12, "Burn-Forest",               
                                                                                                                 ifelse(wolfkde$landcover16 == 13, "Burn-Grassland", 
                                                                                                                        ifelse(wolfkde$landcover16 == 14, "Burn-Shrub", 
                                                                                                                               ifelse(wolfkde$landcover16 == 15, "Alpine Herb", "Alpine Shrub"))))))))))))))))

table(wolfkde$landcover16, wolfkde$used)

table(wolfkde$habitatType, wolfkde$usedFactor)
ggplot(wolfkde, aes(x=landcover16, y=..density.., fill = used)) +geom_histogram(binwidth = 1) + facet_grid(used~.)

#What should we do about the NA or Clouds? We will have to discuss what to do with NA's and Cloud? For now, we will decide to remove clouds as missing data
wolfkde2 <- wolfkde[wolfkde$landcover16 != 11, ]
wolfkde3 <-wolfkde2[wolfkde2$landcover16 != 0, ]
table(wolfkde3$habitatType, wolfkde3$usedFactor)

#Next we will create a 'legend' file (names.m) to help us keep track of contrasts
names.m = data.frame(unique(wolfkde3$landcover16),unique(wolfkde3$habitatType))
# Now I put it order
names.m = names.m[order(names.m)[1:15],]
names.m

#Define a factor variable, landcov.f, # the sorted table makes defining the names of your factor level easy!
wolfkde3$landcov.f = factor(wolfkde3$landcover16,labels = names.m$unique.wolfkde3.habitatType)

#Note that there are many alternative ways of defining your landcover/habitattype as a factor. This method seemed most explicit in terms of defining the design matrix for landcover categories. 
table(wolfkde3$landcov.f, wolfkde3$usedFactor)
table(wolfkde3$landcov.f, wolfkde3$landcover16)

#############################################################
# Univariate Selection Ratio's
#############################################################

table(wolfkde3$habitatType, wolfkde3$usedFactor)

## Estimating Proportions

landcovSelection <- table(wolfkde3$habitatType, wolfkde3$usedFactor)
landcovSelection2 <- as.data.frame.matrix(landcovSelection)
colnames(landcovSelection2)[1:2] <- c("avail","used")
## Calculate Proportional Availability
sum(landcovSelection2$used)
landcovSelection2$pUse <- landcovSelection2$used /413
sum(landcovSelection2$avail)
landcovSelection2$pAvail <- landcovSelection2$avail /1996 # note 2000 because of censored cloud and NA's. 
landcovSelection2

## Calculating Selectivity
landcovSelection2$selection <- landcovSelection2$pUse / landcovSelection2$pAvail
plot(landcovSelection2$selection)

#Next, lets compare the calculation of selection ratio from just the # of locations (incorrect). 
landcovSelection2$selectionN <- landcovSelection2$used / landcovSelection2$avail
plot(landcovSelection2$selection, landcovSelection2$selectionN)

landcovSelection2

## Selectivity Coefficient, the Ln-Selection Ratio
Next we take the natural logarithm, ln() which in R is represented by log()
landcovSelection2$lnSelection <- log(landcovSelection2$selection)

## Lets make a new column of habitatType
landcovSelection2$landcoverType <- c("Alpine Herb", "Alpine Shrub", "Burn-Forest", "Burn-Grassland", "Burn-Shrub", "Closed Conifer", "Deciduous", "Herbaceous", "Mixed", "Moderate Conifer" ,"Open Conifer", "Regen", "Rock-Ice", "Shrub", "Water")

## lets make a plot of the Manly (ln) Selectivity Coefficients
ggplot(data=landcovSelection2, aes(x=landcoverType, y = lnSelection)) + geom_point(size=4) + theme(axis.text.x = element_text(angle = 90))

## it might be handy to save this
write.table(landcovSelection2, "Data/wolfselection.csv", sep=",", row.names = TRUE, col.names=TRUE)
#str(landcovSelection2)

## Selection ratio
ggplot(landcovSelection2, aes(x=landcoverType, y = selection)) + geom_bar(stat="Identity") + theme(axis.text.x = element_text(angle = 90))
## Ln-Selection Ratio
ggplot(landcovSelection2, aes(x=landcoverType, y = lnSelection)) + 
geom_bar(stat="Identity") + theme(axis.text.x = element_text(angle = 90))

## Fancier ggplot
ggplot(landcovSelection2, aes(x=selection, y = lnSelection)) + stat_smooth()

#############################################################
# Selection Ratio's in adehabitatHS
#############################################################

## Estimated available proportions on design I data
elk.avail <- c(15, 61, 84, 40)
elk.used <- c(3, 90, 181, 51)
names(elk.used) <- c("0%", "1-25%", "26-75%", ">75%")
names(elk.avail) <- names(elk.used)
## Computation of wi
(wiRatio <- widesI(elk.used, elk.avail, avknown=FALSE))

## plot the values of the selection ratios
plot(wiRatio)

#############################################################
# Categorical Logistic Regression
#############################################################

?contrast
contrasts(wolfkde3$landcov.f) = contr.treatment(15) 
### To see the design matrix assigned
attributes(wolfkde3$landcov.f)
levels(wolfkde3$landcov.f)

levels(wolfkde3$landcov.f)[11:13] = "Burn"
## note this then reduces us from 15 to 13 categories
contrasts(wolfkde3$landcov.f) = contr.treatment(13)
attributes(wolfkde3$landcov.f)

## Incorrectly Treating Landcover as Continuous
naive.nf = glm(used~landcover16,data=wolfkde3, family=binomial(logit))
summary(naive.nf)

oc = glm(used~I(landcov.f=="Open Conifer"),data=wolfkde3, family = binomial(logit))
summary(oc)
#str(summary(oc))
exp(-1.622+0.711*1)/(1+exp(-1.622+0.711*1))
## now compare to the probability of wolf use in non-conifer landcovers ?
exp(-1.622+0.711*0)/(1+exp(-1.622+0.711*0))

landcovSelection2

# Multiple Logistic Regression with Multiple Categories
## with just open conifer and burns
ocb = glm(used~I(landcov.f=="Open Conifer")+I(landcov.f=="Burn"), data = wolfkde3, family = binomial(logit))
summary(ocb)

### and with a few more variables
conif = glm(used~I(landcov.f=="Open Conifer")+I(landcov.f=="Moderate Conifer")
                  +I(landcov.f=="Closed Conifer"), data = wolfkde3, family = binomial(logit))
summary(conif)

##  Full model with all categories considered
# Full model
full = glm(used~I(landcov.f), data=wolfkde3, family = binomial(logit))
summary(full)

#Discussion: What is the intercept? Where did alpine (landcover 15) go? Why did landcover types 4 (decid), 6 (regen) and alpine- herb (12) 'blow' up? Go back and look at this table to undestand

table(wolfkde3$landcov.f, wolfkde3$usedFactor)

#They blew up because there was 0 used observed.  See what its trying to estimate?

exp(-0.974 - 15.592*1)/(1+exp(-0.974 - 15.592*1)) 

## Models Without an Intercept
full.NoInt = glm(used~I(landcov.f) -1, data=wolfkde3, family = binomial(logit))
summary(full.NoInt)

full.model = glm(used~I(landcov.f=="Moderate Conifer")+I(landcov.f=="Closed Conifer") +I(landcov.f=="Deciduous")+I(landcov.f=="Mixed")+I(landcov.f=="Herbaceous")+I(landcov.f=="Regen")+I(landcov.f=="Shrub")+I(landcov.f=="Water")+I(landcov.f=="Rock-Ice") +I(landcov.f=="Burn")+I(landcov.f=="Alpine Herb")+I(landcov.f=="Alpine Shrub"), data = wolfkde3, family = binomial(logit))
summary(full.model)

# Changing the Reference Category in Logistic Regression
## first recheck which # Rock-Ice is
levels(wolfkde3$landcov.f) ## Ok it is # 10

contrasts(wolfkde3$landcov.f) = contr.treatment(13, base = 10)
attributes(wolfkde3$landcov.f)
# and note that rock-ice now is 0. 

rockintercept.model = glm(used~I(landcov.f=="Moderate Conifer") +I(landcov.f=="Closed Conifer") +I(landcov.f=="Deciduous")+I(landcov.f=="Mixed")+I(landcov.f=="Herbaceous") +I(landcov.f=="Regen")+I(landcov.f=="Shrub")+I(landcov.f=="Water")+I(landcov.f=="Open Conifer")+I(landcov.f=="Burn")+I(landcov.f=="Alpine Herb")+I(landcov.f=="Alpine Shrub"), data = wolfkde3, family = binomial(logit))
summary(rockintercept.model)

  ## Manual Dummy (Indicator) Coding
#In practice I find working through the Design matrix coding of R confusing. Instead, I often just create my own 'manual' dummy variables in my data frame, sometimes even beforehand in excel (gasp!). These next commands manually creating 'dummy' variables that replace using the interaction expansion used ~ I.

wolfkde3$closedConif = ifelse(wolfkde3$habitatType == "Closed Conifer", 1, 0)
wolfkde3$modConif = ifelse(wolfkde3$habitatType == "Moderate Conifer", 1, 0)
wolfkde3$openConif = ifelse(wolfkde3$habitatType == "Open Conifer", 1, 0)
wolfkde3$decid = ifelse(wolfkde3$habitatType == "Deciduous", 1, 0)
wolfkde3$regen = ifelse(wolfkde3$habitatType == "Regen", 1, 0)
wolfkde3$mixed = ifelse(wolfkde3$habitatType == "Mixed", 1, 0)
wolfkde3$herb = ifelse(wolfkde3$habitatType == "Herbaceous", 1, 0)
wolfkde3$shrub = ifelse(wolfkde3$habitatType == "Shrub", 1, 0)
wolfkde3$water = ifelse(wolfkde3$habitatType == "Water", 1, 0)
wolfkde3$rockIce = ifelse(wolfkde3$habitatType == "Rock-Ice", 1, 0)
## note here I reclassified all burn = 1 
wolfkde3$burn = ifelse(wolfkde3$habitatType == "Burn-Grassland", 1, ifelse(wolfkde3$habitatType == "Burn-Shrub", 1, ifelse(wolfkde3$habitatType == "Burn-Forest", 1,0 )))
wolfkde3$alpineHerb = ifelse(wolfkde3$habitatType == "Alpine Herb", 1, 0)
wolfkde3$alpineShrub = ifelse(wolfkde3$habitatType == "Alpine Shrub", 1, 0)

head(wolfkde3)

wolfkde3$alpine = wolfkde3$alpineHerb + wolfkde3$alpineShrub

#Refitting model with Open Conifer as the intercept and alpine/burn pooled

oc.intercept.model = glm(used~closedConif + modConif + decid+ regen+mixed+herb+water+rockIce+burn+alpine, data = wolfkde3, family = binomial(logit))
summary(oc.intercept.model)

### refitting model with just Alpine and Rock and Ice as the intercept
rockintercept.alpine.model = glm(used~closedConif + openConif + modConif + decid+ regen+mixed+herb+water+burn+alpine, data = wolfkde3, family = binomial(logit))
summary(rockintercept.alpine.model)

### refitting model manually dropping Decid and Regen - where do they no go?
rock.alpine.regen.decid.intercept.model = glm(used~closedConif + openConif + modConif + mixed+herb+water+burn+alpine, data = wolfkde3, family = binomial(logit))
summary(rock.alpine.regen.decid.intercept.model)

## Comparing Coefficients from Models with Different Intercepts

rockintercept.alpine.model.df <- data.frame(summary(rockintercept.alpine.model)$coefficients[,1:2])
oc.intercept.model.df <- data.frame(summary(oc.intercept.model)$coefficients[,1:2])
coef.table <- rbind(rockintercept.alpine.model.df,oc.intercept.model.df)
coef.table$habitatType <- c(row.names((summary(rockintercept.alpine.model)$coefficients[,1:2])),row.names(summary(oc.intercept.model)$coefficients[,1:2]))
coef.table$habitatType[1] <- "rockIce"
coef.table$habitatType[12] <- "openConif"
coef.table$model <-c(rep("Open Conif Intercept",11),rep( "RockIce Intercept",11))
coef.table

#Now use this table to compare the ABSOLUTE differences say between burn and alpine in both models. In the oc.model the B coefficient for Burn = 2.47 and Alpine = -3.002 an absolute differences of 3.29. In the rock.model the B coefficient for Burn = 2.2 and Alpine is -1.086, an absolute difference of 3.29 - the same! Why is this?

# figure of the Beta coefficients 
ggplot(coef.table, aes(x=habitatType, y=Estimate, colour=model)) + geom_point(size = 5) + theme(axis.text.x = element_text(angle = 90))
