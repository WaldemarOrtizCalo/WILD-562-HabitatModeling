# title: "WILD 562 Lab 6: Evaluating RSF Models"
# author: "Mark Hebblewhite"

  
  # Lab 6: Evaluating RSF Models
  
## 0.1 Preliminaries: setting packages

#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("ROCR", "sp", "raster", "tidyverse", "ks", "mapview", "rgdal", "sp", "raster","colorRamps","rgeos", "VGAM", "AICcmodavg", "MuMIn", "corrgram", "GGally","caret", "DescTools", "car")

#run function to install packages
ipak(packages)

## 0.2 Saving and loading data and shapefile data of Kernel home range from Lab 

wolfkde <- read.csv("1.Data/Lab6_data/wolfkde.csv", header=TRUE, sep = ",", na.strings="NA", dec=".")
wolfkde3 <-na.omit(wolfkde)
wolfkde3$usedFactor <-as.factor(wolfkde3$usedFactor) ## make sure usedFactor is a factor

kernelHR <- shapefile("1.Data/Lab6_data/homerangeALL")
plot(kernelHR)
extent(kernelHR)
kernels <- raster()
extent(kernels) <- c(xmin=546836, xmax=612093, ymin=5662036, ymax=5748911) 


## 0.4 Loading raster's needed for mapping RSF models later
deer_w<-raster("1.Data/Lab6_data/rasters/deer_w2.tif")
moose_w<-raster("1.Data/Lab6_data/rasters/moose_w2.tif")
elk_w<-raster("1.Data/Lab6_data/rasters/elk_w2.tif") # already brought in above
sheep_w<-raster("1.Data/Lab6_data/rasters/sheep_w2.tif")
goat_w<-raster("1.Data/Lab6_data/rasters/goat_w2.tif")
wolf_w<-raster("1.Data/Lab6_data/rasters/wolf_w2.tif")
elevation2<-raster("1.Data/Lab6_data/rasters/Elevation2.tif") #resampled
disthumanaccess2<-raster("1.Data/Lab6_data/rasters/DistFromHumanAccess2.tif") #resampled in lab 4
disthhu2<-raster("1.Data/Lab6_data/rasters/DistFromHighHumanAccess2.tif") #resampled in lab 4
landcover2 <- raster("1.Data/Lab6_data/rasters/landcover2.tif") ## resampled to same extent as lab 4
extent(landcover2)
plot(landcover2)
plot(kernelHR, add=TRUE)


landcover2@data@values <- getValues(landcover2)
table(landcover2@data@values)

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


## 0.5 Creating 'dummy' variable landcover rasters for use in mapping later 

alpine <- mask.raster
burn <- mask.raster
closedConif <- mask.raster
herb <- mask.raster
mixed <- mask.raster
rockIce <- mask.raster
water <- mask.raster
modConif <- mask.raster
decid <- mask.raster

## Create dummy rasters
alpine@data@values <- ifelse(landcover2@data@values== 15 | landcover2@data@values == 16, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
burn@data@values <- ifelse(landcover2@data@values == 12 | landcover2@data@values == 13 | landcover2@data@values == 14, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
closedConif@data@values <- ifelse(landcover2@data@values == 3, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
herb@data@values <- ifelse(landcover2@data@values == 7, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
mixed@data@values <- ifelse(landcover2@data@values == 5, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
rockIce@data@values <- ifelse(landcover2@data@values == 10, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
water@data@values <- ifelse(landcover2@data@values == 9, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
modConif@data@values <- ifelse(landcover2@data@values == 2, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
decid@data@values <- ifelse(landcover2@data@values == 10, 1, ifelse(is.na(landcover2@data@values)==T,NA,0))
plot(rockIce)
plot(closedConif)

# note that open conifer is implicitly being defined as the intercept

## 0.6 Creating a Raster Stack 
#A Raster stack could help with mapping later, but, its not really needed this lab. Recall all rasters must have same extent and resolution, which will help. 

all_rasters <- stack(deer_w, moose_w, elk_w, sheep_w, goat_w, 
                     wolf_w,elevation2, disthumanaccess2, disthhu2, 
                     landcover2, alpine, burn, closedConif, modConif, 
                     herb, mixed, rockIce, water, decid)

plot(all_rasters) ## note limit of plotting 9 layers

#Note, the next two sets of commands give a way to export the processed rasters, all of them, to a new folder in similar format.  We will skip this step today in Lab, but try it 

#names = c("deer_w", "moose_w", "elk_w", "sheep_w", "goat_w", "wolf_w","elevation2", "disthumanaccess2", "disthhu2", "landcover2", "alpine", "burn", "closedConif", "modConif", "herb", "mixed", "rockIce", "water", "decid")

#writeRaster(all_rasters,"lab6Stack.tif", bylayer = TRUE,suffix = 'names', format="GTiff")


# Multiple Logistic Regression 

## top Biotic model was model 41
top.biotic <- glm(used ~ DistFromHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=wolfkde3)
summary(top.biotic)
#double check the VIF for each final model, just to be sure
vif(top.biotic)

#Environmental Model - top model is model 11
top.env <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3)
summary(top.env)
vif(top.env)



# Calculate model selection AIC table

models = list(top.biotic, top.env)
modnames = c("top biotic", "top env")
aictab(cand.set = models, modnames = modnames)
# aictab(models, modnames) ## short form. 

## Evaluating Predictions from Residual Diagnostic Plots
par(mfrow = c(2,2))
plot(top.env)

##### Saving predictions manually, an example with the environment model
wolfkde3$fitted.top.env <- fitted(top.env)
#### this is the predicted probability from the model
wolfkde3$residuals.top.env <- residuals(top.env)
## these are the deviations from the predictions for each row (data point)
wolfkde3$rstudent.top.env <- rstudent(top.env)
## This is a standardized residual - the studentized residual
wolfkde3$hatvalues.top.env <- hatvalues(top.env)
#### this is the first of the leverage statistics, the larger hat value is, the bigger the influence on the fitted value
wolfkde3$cooks.distance.top.env <- cooks.distance(top.env)
#### this is the Cooks leverage statistic, the larger hat value is, the bigger the influence on the fitted value
wolfkde3$obsNumber <- 1:nrow(wolfkde3) ## just added a row number for plotting

#Making manual residual verus predicted plots
ggplot(wolfkde3, aes(fitted.top.env, residuals.top.env)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))
## Manual residual plot

ggplot(wolfkde3, aes(wolfkde3$residuals.top.env, wolfkde3$cooks.distance.top.env)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))
## shows us some points at high cooks values that might be having a big influence

ggplot(wolfkde3, aes(wolfkde3$cooks.distance.top.env, wolfkde3$hatvalues.top.env)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))

#This shows us some points at high cooks values that might be having a big influence. This helps identify some locations that have high leverage that are used (1 - blue) and available (0=black) points. For example, datapoint 16
# e.g., datapoint 16
wolfkde3[16,]

#This is a red deer wolf used point at high elevtion in open conifer far from human access.  The predicted probability is 0.03, yet it was used = 1.  We will return to this point below. 

#Next we will explore evaluating model fit graphically. First lets make a plot of Y against X predictions...
scatterplot(fitted.top.env~Elevation2, reg.line=lm, smooth=TRUE, spread=TRUE, boxplots='xy', span=0.5, xlab="elevation", ylab="residual", cex=1.5, cex.axis=1.4, cex.lab=1.4, data=wolfkde3)

hist(wolfkde3$fitted.top.env, scale="frequency", breaks="Sturges", col="darkgray")

#This last histogram plot is VERY important - it is the predicted probability of a location being a wolf used location given your top model. But how would we divide this into wolf habitat and wolf available? 
ggplot(wolfkde3, aes(x=wolfkde3$fitted.top.env, fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) #+ facet_grid(pack ~ ., scales="free")

#This plot shows that somewhere around 0.25 - 0.40 eyeballing it it looks like we could 'cut' used and available points? 
ggplot(wolfkde3, aes(x=fitted.top.env, y=..density.., fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) + facet_grid(pack ~ ., scales="free")

#But note this 'cut' point looks different for both wolf packs?  

## Evaluating predictions from residual plots for top.biotic model
par(mfrow = c(2,2))
plot(top.biotic)

##### Saving predictions manually, an example with the environment model
wolfkde3$fitted.top.biotic <- fitted(top.biotic)
#### this is the predicted probability from the model
wolfkde3$residuals.top.biotic <- residuals(top.biotic)
## these are the deviations from the predictions for each row (data point)
wolfkde3$rstudent.top.biotic <- rstudent(top.biotic)
## This is a standardized residual - the studentized residual
wolfkde3$hatvalues.top.biotic <- hatvalues(top.biotic)
#### this is the first of the leverage statistics, the larger hat value is, the bigger the influence on the fitted value
wolfkde3$cooks.distance.top.biotic <- cooks.distance(top.biotic)
#### This isthe Cooks leverage statistic

## Making manual residual verus predicted plots
ggplot(wolfkde3, aes(wolfkde3$cooks.distance.top.biotic, wolfkde3$hatvalues.top.biotic)) + geom_point() + geom_text(aes(label = obsNumber, colour = used))

wolfkde3[30,]
#This is another high wolf used point that is being classified as an AVAILABLE location. **Misclassification !!! **

#Next we will evaluate model fit graphically, first making a plot of Y against X predictions... for Elevation. 
scatterplot(fitted.top.biotic~Elevation2, reg.line=lm, smooth=TRUE, spread=TRUE, boxplots='xy', span=0.5, xlab="elevation", ylab="residual", cex=1.5, cex.axis=1.4, cex.lab=1.4, data=wolfkde3)

## next, the histogram of predicted probaiblities
hist(wolfkde3$fitted.top.biotic, scale="frequency", breaks="Sturges", col="darkgray")
#This plot is VERY important - it is the predicted probability of a location being a wolf used location given your top model. But how would we divide this into wolf habitat and wolf available? 

ggplot(wolfkde3, aes(x=wolfkde3$fitted.top.biotic, fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) #+ facet_grid(pack ~ ., scales="free")

ggplot(wolfkde3, aes(x=fitted.top.biotic, y=..density.., fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) + facet_grid(pack ~ ., scales="free")
#### But note this 'cut' point looks different for both wolf packs?

## Comparing the 'Fit' of Predictions From The Biotic and Environmental Models
ggplot(wolfkde3, aes(x=fitted.top.biotic, y=fitted.top.env)) + geom_point() + stat_smooth(method="lm")

#What about the comparison of predictions from these two models when split by wolf packs?
ggplot(wolfkde3, aes(x=fitted.top.biotic, y=fitted.top.env, fill = pack)) + geom_point() + stat_smooth(method="lm") 

#But what if we do not force a linear model through the predictions? 
ggplot(wolfkde3, aes(x=fitted.top.biotic, y=fitted.top.env, fill = pack)) + geom_point() + stat_smooth()


# Classification Tables

## Pseuod R-squared 
#require(DescTools)
PseudoR2(top.biotic, c("McFadden", "CoxSnell", "Nagel"))

## Classification Table for Top Biotic Model

# First we will arbitrarily define the cutpoint between 1 and 0's using p = 0.5
ppused = wolfkde3$fitted.top.biotic > 0.5
table(ppused,wolfkde3$used)

#Next, we will go through Hosmer and Lemeshow Chapter 5 to calculate our classification success for 1's?
167/(167+229)
#So when wolf telemetry locations were known = 1, the model classified 42% as 'used'. *This is pretty terrible, don't you think?* This is also called the **Specificity of a Model, i.e., how well the model specifies TRUTH when = 1**. Another synonym is the True Positive Rate (TPR)

#Now lets do the correct classification rate of the true 0's (i.e., avail)
1655/(1655+67)
#But when the points were really 0s, available, we classified them 96% of the time correctly. **This is called the Sensitivity of a model, the correct classification of 0s.** Otherwise known as the True Negative Rate. 

#Now lets do this manually using cutpoints of 0.25 and 0.1
ppused = wolfkde3$fitted.top.biotic>0.25
table(ppused,wolfkde3$used)
304/(304+92) # specificity
1376 / (1376+346) # sensitivity

#Now lets try a p = 0.10
ppused = wolfkde3$fitted.top.biotic>0.10
table(ppused,wolfkde3$used)
357/(357+39)
1001 / (1001+721)

#Finally, lets try a p = 0.70
ppused = wolfkde3$fitted.top.biotic>0.70
table(ppused,wolfkde3$used)
66/(330+66)
1705 / (1705+17)
                    
# Confusion Matrices
#require(caret)
wolfkde3$pr.top.biotic.used <- ifelse(wolfkde3$fitted.top.biotic>0.5, 1, 0)
xtab1<-table(wolfkde3$pr.top.biotic.used, wolfkde3$used)
xtab1

#?confusionMatrix
confusionMatrix(xtab1) ## Note this is incorrectly specifying what the 1's are. 

confusionMatrix(xtab1, positive = "1") ## Thanks to spring 2019 student Forest HAyes for finding this mistake :)
    
## Receiver Operating Characterstic ROC curves
## Sensitivity and Specificity
#require(ROCR)
pp = predict(top.biotic,type="response")
pred = prediction(pp, wolfkde3$used)

perf3 <- performance(pred, "sens", x.measure = "cutoff")
plot(perf3)

perf4 <- performance(pred, "spec", x.measure = "cutoff")
plot(perf4)

## Estimating the Optimal Cutpoint
perfClass <- performance(pred, "tpr","fpr") # change 2nd and/or 3rd arguments for other metrics
fpr <- perfClass@x.values[[1]]
tpr <- perfClass@y.values[[1]]
sum <- tpr + (1-fpr)
index <- which.max(sum)
cutoff <- perfClass@alpha.values[[1]][[index]]
cutoff
#Where does this value come from?
table(wolfkde3$used)
396/(1722+396)

#Now, lets overlay the sensitivity, specificity, and optimal cutoff curves together. 
plot(perf3, col="blue") # Sensitivity
plot(perf4, add = TRUE) # Specificity
abline(v=cutoff, col="red") ## optimal cutpoint

## ROC Plot 
plot(perfClass)
abline(a=0, b= 1)

BMauc <- performance(pred, measure="auc") 
str(BMauc)
auc <- as.numeric(BMauc@y.values)
auc

#Next, we will Plot ROC Curve with the optimal cut point and AUC
plot(perfClass, colorize = T, lwd = 5, print.cutoffs.at=seq(0,1,by=0.1),
     text.adj=c(1.2,1.2),
     main = "ROC Curve")
text(0.5, 0.5, "AUC = 0.867")
abline(v=cutoff, col = "red", lwd = 3)

acc.perf = performance(pred, measure = "acc")
plot(acc.perf)

## Manually Changing Cutoff Values
### now lets try a p = of our cutoff
ppused = wolfkde3$fitted.top.biotic>cutoff
table(ppused,wolfkde3$used)
#### Now, what is specificity? (i.e., the probability of classifying the 1's correctly?)
320/(320+76)
#### about 80% - Great! But - what happened to our sensitivity (i.e., the probability of classifying the 0's correctly?)
1344 / (1344+378)
#### so our probability of classifying 0's correctly is about 78%

#Lets look at the confusion matrix now for the optimal cutpoint. 
wolfkde3$pr.top.biotic.used2 <- ifelse(wolfkde3$fitted.top.biotic>cutoff, 1, 0)
xtab2<-table(wolfkde3$pr.top.biotic.used2, wolfkde3$used)
xtab2

#?confusionMatrix
confusionMatrix(xtab2, positive = "1")

#Now lets use this cutoff to classify used and avail locations into 1 and 0's, and make a plot of where this cutoff is using geom_vline() in ggplot
ggplot(wolfkde3, aes(x=wolfkde3$fitted.top.biotic, fill=usedFactor)) + geom_histogram(binwidth=0.05, position="identity", alpha=0.7) + xlab("Predicted Probability of Wolf Use") + theme(axis.title.x=element_text(size=16)) + geom_vline(xintercept = cutoff, col="red")

## EXCERCISE 2: Evaluating the top Environmental Model 

# K-folds Cross Validation
## Evaluating the Top Environmental Model with k-folds
source("2.R_Code/1.Scripts/kxv.R", verbose = FALSE)

#Next, we will fit a 5-fold k=5 cross validation of the top.env model
# Kfolds with a 'fuzz' factor
kxvPrintFlag=FALSE
kxvPlotFlag=TRUE
kxvFuzzFactor = 0.01
kfolds = kxvglm(top.env$formula, data=wolfkde3, k=5, nbin=10) ## note we get a lot of ties here and some error messages, this is because of all the categories. Read more about this in the clumpy data.pdf in the Vignette folder.  
kfolds

# Kfolds by each pack with a 'fuzz' factor

kxvPrintFlag=FALSE
kxvPlotFlag=TRUE
kxvFuzzFactor = 0.01
kfolds2 = kxvglm(top.env$formula, data=wolfkde3, k=5, nbin=10, partition="pack")
kfolds2
                                                                               
## Excercise 3: Evaluating the top Biotic Model - On your own for homework - again, a great excercise. 
 
## Manual k-folds Cross-Validation

# 1. Create a vector of random "folds" in this case 5, 1:5
wolfkde3$rand.vec = sample(1:5,nrow(wolfkde3),replace=TRUE)

#2. Run the model for 1 single random subset of the data == 1
top.env.1= glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde3, subset=rand.vec==1) ## note this last subset = rand.vec==1 is just 1 of the 5 subsets. 

# 3. Make predictions for points not used in this random subset (2:5) to fit the model.
pred.prob = predict(top.env.1,newdata=wolfkde3[wolfkde3$rand.vec!=1,],type="response") ## note that != means everything but subset 1. 

# 4. Make quantiles for the predictions - this calculates the 'bin's of the categories of habitat availability
q.pp = quantile(pred.prob,probs=seq(0,1,.1))

# 5. Then for each of 10 bins, put each row of data into a bin
bin = rep(NA,length(pred.prob))
for (i in 1:10){
  bin[pred.prob>=q.pp[i]&pred.prob<q.pp[i+1]] = i
}

## 5. This then makes a count of just the used locations for all other K folds 2:5 
used1 = wolfkde3$used[wolfkde3$rand.vec!=1]

## We then make a table of them
rand.vec.1.table <- table(used1,bin)
rand.vec.1.table

cor.test(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), c(0,2,0,8,6,15,24,50,99,110), method="spearman") 

rand.vec.1.table <- as.data.frame(rand.vec.1.table)
ggplot(rand.vec.1.table, aes(as.numeric(bin), Freq, col = used1)) + geom_point(size=5) + geom_line()

## Excercise: Conduct k-folds cross validation ont the top-biotic models. 

# Mapping Spatial Predictions of 'Top' Model 

## Easy Spatial Prediction just in ggplot
par(mfrow = c(1,1)) # reset graphical parameters

ggplot(wolfkde3, aes(EASTING, NORTHING, col = fitted.top.biotic)) + geom_point(size=5) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')

ggplot(wolfkde3, aes(EASTING, NORTHING, col = fitted.top.env)) + geom_point(size=5) + coord_equal() +  scale_colour_gradient(low = 'yellow', high = 'red')

## Raster Predictions 
par(mfrow = c(1,1))
summary(top.biotic)

biotic.coefs <- top.biotic$coefficients[c(1:4)]
names(all_rasters)

##### Now we will calculate the probability of wolf use using the logistic regression equation y = exp(x) / (1 + exp(x))
rast.top.biotic <- exp(biotic.coefs[1] + biotic.coefs[2]*disthumanaccess2 + biotic.coefs[3]*deer_w + biotic.coefs[4]*goat_w) / (1 +exp(biotic.coefs[1] + biotic.coefs[2]*disthumanaccess2 + biotic.coefs[3]*deer_w + biotic.coefs[4]*goat_w ))
#Note we need to use the names of the raster layers we brought in up above. Note that they are not the same names as stored in the Raster stack.

wolfyht<-shapefile("1.Data/Lab6_data/wolfyht.shp")
# plot predicted raster within extent of kernels raster
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=kernels)
plot(kernelHR, add=TRUE)
plot(wolfyht, col='blue', pch = 16, add=TRUE)

hist(rast.top.biotic@data@values) 

## Bow Valley zoom
bv.raster<-raster()
extent(bv.raster) <- c(xmin=570000, xmax=600000, ymin=5665000, ymax=5685000) 
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=bv.raster)
plot(kernelHR, add=TRUE)
plot(wolfyht, col='blue', pch = 16, add=TRUE)


##Red Deer
rd.raster<-raster()
extent(rd.raster) <- c(xmin=540000, xmax=600000, ymin=5700000, ymax=5730000) 
plot(rast.top.biotic, col=colorRampPalette(c("yellow", "orange", "red"))(255), ext=rd.raster)
plot(kernelHR, add=TRUE)
plot(wolfyht, col='blue', pch = 16, add=TRUE)
