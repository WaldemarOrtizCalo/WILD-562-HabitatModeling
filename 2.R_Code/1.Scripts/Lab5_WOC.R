# title: 'WILD 562 Lab5 : Multiple Logisitic Regression & Model Selection'
# author: "Mark Hebblewhite"

# Lab5: Multiple Logisitic Regression & Model Selectionn


## 0.1 Preliminaries: setting packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load or install these packages:
packages <- c("car", "tidyverse", "MASS", "AICcmodavg", "MuMIn", "corrgram", "GGally", "bootStepAIC", "broom")

#run function to install packages
ipak(packages)


## 0.2 Preliminaries: importing data
wolfkde <- read.csv("Data/wolfkde.csv", header=TRUE, sep = ",", na.strings="NA", dec=".")
head(wolfkde)
table(wolfkde$pack, wolfkde$used)

# Multiple Logistic Regression & Collinearity

## Revisiting Univariate Regressions from Lab 2

## First lets fit Univariate models of these 2 covariates
elev <- glm(used~Elevation2, data =wolfkde, family= binomial(logit))
disthhacc <-  glm(used~DistFromHighHumanAccess2, data =wolfkde, family= binomial(logit))

# Next, fit both in our first multiple logistic regression model
elev.disthhacc <- glm(used~Elevation2 +DistFromHighHumanAccess2 , data =wolfkde, family= binomial(logit))
summary(elev.disthhacc)
# now lets extract coefficients

summary(elev)$coefficients[,1:2]
summary(disthhacc)$coefficients[,1:2]
summary(elev.disthhacc)$coefficients[,1:2]

## lets visually explore differences
disthumanBnp = 0:7000
prDisthhacc <- predict(disthhacc, newdata=data.frame(DistFromHighHumanAccess2=disthumanBnp), type="response")
head(prDisthhacc)

plot(wolfkde$DistFromHighHumanAccess2, wolfkde$used)
lines(disthumanBnp, prDisthhacc, type="l", ylab= "Pr(Used)")

summary(wolfkde$Elevation2)
## ok - lets evaluate the probability of use at 1931 meters from the elev.disthhacc model
medianElev = 1931
prElevMedian.Disthhacc <- predict(elev.disthhacc, newdata=data.frame(DistFromHighHumanAccess2=disthumanBnp, Elevation2=medianElev), type="response")


# Lets take a look at the predicted probabilities from our original univariate regression model - plotted above, and, second, the predicted probabilities at the mean elevation from the multivariate model. 

plot(wolfkde$DistFromHighHumanAccess2, wolfkde$used, xlim=(c(0,10000)))
lines(disthumanBnp, prElevMedian.Disthhacc, type="l", ylab= "Pr(Used)")
lines(disthumanBnp, prDisthhacc, type="l", ylab= "Pr(Used)")

#What is going on?? Why did the coefficient switch sign plotted at the median elevation of 1931m?

##  Partial Regression Coefficients 
newdata <- expand.grid(Elevation2 = pretty(wolfkde$Elevation2, 5), DistFromHighHumanAccess2 = pretty(wolfkde$DistFromHighHumanAccess2, 10))
head(newdata)
newdata$prElev.Disthha <-predict(elev.disthhacc, newdata, type="response")

ggplot(newdata, aes(x = DistFromHighHumanAccess2, y = prElev.Disthha)) + geom_line() + facet_wrap(~Elevation2)

## Collinearity and Correlations
#Why are Elevation and DistHighHumanUse changing? They are CORRELATED! Lets formally test the new assumption of multiple linear regression using Pearsons correlation test.  
cor.test(wolfkde$Elevation2, wolfkde$DistFromHighHumanAccess2)


#Second, we will fit a linear model of distance as a function of elevation to see the regression coefficient between the two, and finally, we will plot elevation and distance using two approaches. 
elev_disthha <- lm(DistFromHighHumanAccess2~Elevation2, data=wolfkde)
summary(elev_disthha)

plot(wolfkde$Elevation2,wolfkde$DistFromHighHumanAccess2, type="p")
abline(lm(DistFromHighHumanAccess2~Elevation2, data=wolfkde), col="red")

#Can we hold effects of X1 constant while varying X1? Lets use the pairs() to visualize the correlation with elevation as X and Y. 
pairs(~Elevation2+DistFromHighHumanAccess2, data=wolfkde, main="Scatterplot Matrix")

#Lets examine our 'study design'
wolfkde.Used <- subset(wolfkde, wolfkde$used==1)
wolfkde.Used$elev2F <- cut(wolfkde.Used$Elevation2, 2)
wolfkde.Used$DistFromHighHumanAccess2.2F <- cut(wolfkde.Used$DistFromHighHumanAccess2, 2)
table(wolfkde.Used$elev2F, wolfkde.Used$DistFromHighHumanAccess2.2F)

#Note the 2 elevaition classes, low (1.4e+03,1.81e+03]) and high, are the ROWS, and the 2 distance classes are the COLUMNS (close to humans, far from humans).
#What this table tells us is that most of our wolf 'used' locations occurs at low elevations in areas close to human access (n=336), and that we have very few locations at high elevations, far from humans, and the combination. This visualizes the experimental confound that is present in our data. This question is unresolvable with the current data. 

# However, the guidelines range from a $\rho$ of 0.3 - 0.7 according to various sources. 

#The key reason why there is no hard and fast guidelines for a threshold that is set and stone for collinearity is because collinearity is really a sign of the bigger experimental and statistical design flaw, confounding. 

## Confounding
summary(elev)$coefficients[,1:2]
summary(disthhacc)$coefficients[,1:2]
summary(elev.disthhacc)$coefficients[,1:2]

#We note that the sign for Elevation is quite stable when alone and when in the presence of distance to high human access. Alone, it is -0.0055, together, it is  -0.0071.  This is a magnitude of change of about 27% from alone to when in a multiple logistic model.

#Conversely, the coefficient for distance varies dramatically from  -0.0002  when alone in a univariate model, to 0.00023 when in the presence of elevation.  This is a change of - 115%, a complete sign flip in the opposite direction.  The coefficient for elevation is what I call a _stable_ coefficient. In contrast, the coefficient for distance is _unstable_ when in the presence of other covariates. This is a sign of both collinearity and confounding, and points to the experimental design challenge of being able to say anything about distance to high human access. 

# Lets next test 2 other variables, elk and deer...
## lets test 2 other variables, elk and deer...
deer <- glm(used~deer_w2, data =wolfkde, family= binomial(logit))
elk <-  glm(used~elk_w2, data =wolfkde, family= binomial(logit))

# Next, fit both in our first multiple logistic regression model
deer.elk <- glm(used~deer_w2 + elk_w2, data =wolfkde, family= binomial(logit))

# now lets extract coefficients
summary(deer)$coefficients[,1:2]
summary(elk)$coefficients[,1:2]
summary(deer.elk)$coefficients[,1:2]

#Note this time the sign's didn't flip, but, they significantly changed, weakening in the presence of each other. The coefficient for deer was 1.17 alone, 0.74 together, a 36% change!  For elk, it was 1.12 alone, 0.633 together, representing a 44% change. Both are dramatic changes, suggestive of a correlation.    Which we confirm:
cor.test(wolfkde$deer_w2, wolfkde$elk_w2)
plot(wolfkde$deer_w2,wolfkde$elk_w2, type="p")
abline(lm(elk_w2~deer_w2, data=wolfkde), col="red")

#We confirm that the deer H.S.I. and elk H.S.I. scores are correlated with a $\rho$ = 0.86, representing essentially the exact same variable.  But doing these collinearity screens one at a time is very tedious.  In the next section we will advance to checking all covariates in a dataset at the same time. 

#  Screening for Collinearity

## Continuous variables

plot(wolfkde$Elevation2 ,wolfkde$goat_w2, type="p")
abline(lm(goat_w2~Elevation2, data=wolfkde), col="red")
## graphically examining collinearity

## Scatterplot Matrices

pairs(~Elevation2+DistFromHumanAccess2+DistFromHighHumanAccess2, data=wolfkde, main="Scatterplot Matrix")
pairs(~deer_w2+moose_w2+elk_w2+sheep_w2+goat_w2+Elevation2, data=wolfkde, main="Scatterplot Matrix")
#Here again we see the now familiar wedge shaped, triangular distribution that is indicative of our poor experimental design again with respect to the sampling of wolf territories across an elevation and distance from human gradient.  

### ScatterplotMatrix from Library car
## using car library
scatterplotMatrix(~deer_w2+moose_w2+elk_w2+sheep_w2+goat_w2+Elevation2, data=wolfkde, main="Scatterplot Matrix")

scatterplotMatrix(~Elevation2+DistFromHumanAccess2+DistFromHighHumanAccess2, data=wolfkde, main="Scatterplot Matrix")

scatterplotMatrix(~deer_w2+moose_w2+elk_w2+sheep_w2+goat_w2+Elevation2+DistFromHumanAccess2+DistFromHighHumanAccess2, data=wolfkde, main="Scatterplot Matrix")

### Using Library corrgram
#This makes a scatterplot matrix with the upper panel defined as the Pearsons correlation coefficient expressed in Pie graph form (where r = 1.0), red = negative correlation and blue equals a positive correlation. Similarly, the bottom panel displays the strength of the correlation in shaded color. 
corrgram(wolfkde[1:9], order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Correlations in the Wolf Data")

corrgram(wolfkde[1:9], order=TRUE, lower.panel=panel.pts,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Correlations in the Wolf Data")

corrgram(wolfkde[1:9], order=TRUE, lower.panel=panel.ellipse,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Correlations in the Wolf Data")

### ggcorr
## using the ggcorr package
ggcorrplot <- ggcorr(wolfkde[1:9], label = TRUE)
ggcorrplot
## GGally package with ggpairs()
ggpairplot<-ggpairs(wolfkde[1:9])
ggpairplot

### Multicollinearity Function
cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X, use="complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R
}

cor.prob(as.matrix(wolfkde[,c("deer_w2","elk_w2", "moose_w2", "sheep_w2", "goat_w2", "Elevation2", "DistFromHumanAccess2", "DistFromHighHumanAccess2")]))

cor.prob2 <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X, use="complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  Rstar = ifelse(R[above]<0.05, "***", "NS")
  R[above]=paste(R[above],Rstar)
  R
}

cor.prob2(as.matrix(wolfkde[,c("deer_w2","elk_w2", "moose_w2", "sheep_w2", "goat_w2", "Elevation2", "DistFromHumanAccess2", "DistFromHighHumanAccess2")]))

#So which covariates have the highest correlations??  Deer, Elk, and Moose all have correlation coefficients > 0.65; Sheep and Goats are correlated > 0.4; elevation is inversely correlated with an R of -0.75 with deer, elk , moose

### Variance Inflation Factors

?vif()

full.model = glm(used~deer_w2 + elk_w2 +moose_w2 +sheep_w2+goat_w2+Elevation2+DistFromHumanAccess2+DistFromHighHumanAccess2, data =wolfkde, family= binomial(logit))
summary(full.model)

vif(full.model)

## Collinearity amongst Categorical Variables
cor.test(wolfkde$alpine, wolfkde$Elevation2)
cor.test(wolfkde$burn, wolfkde$Elevation2)
cor.test(wolfkde$closedConif, wolfkde$Elevation2)
cor.test(wolfkde$herb, wolfkde$Elevation2)
cor.test(wolfkde$goat_w2, wolfkde$Elevation2)

#Lets do all variables, continuous and categorical together. 
cor.prob(as.matrix(wolfkde[,c("Elevation2", "DistFromHumanAccess2", "openConif", "closedConif", "modConif", "burn", "herb", "decid", "burn", "alpine")]))

#Lets plot the correlations between Elevation [7] and landcover types [18:29]
corrgram(wolfkde[c(7, 18:29)], order=TRUE, lower.panel=panel.fill,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Landcover Correlations with Elevation")
#So nothing too egregious except, unsurprisingly, Rock and Ice and Elevaiton. 

#Next, lets test for correlations between human access[8], and the landcover dummy variables. 

corrgram(wolfkde[c(8, 18:29)], order=TRUE, lower.panel=panel.cor,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Landcover Correlations with Distance from Human Access")
#Again, only issue is Rock and Ice but even then its not a huge effect. We can essentially see this here: 
boxplot(Elevation2~landcov.f, ylab="Elevation (m)", data=wolfkde, las=3)
boxplot(DistFromHumanAccess2~landcov.f, ylab="Elevation (m)", data=wolfkde, las=3)
#So collinearity is not as important for categorical variables but it becomes important if we start to assess categorical interactions with continuous factors. 

## Interaction Between Categorical Factors and Continuous
wolfkde$closed = 0
wolfkde$closed <- wolfkde$closedConif + wolfkde$modConif + wolfkde$openConif + wolfkde$decid + wolfkde$mixed + wolfkde$burn
## note I considered burn here as 'closed' - could change. 

wolfkde$closedFactor <-as.factor(wolfkde$closed)

ggplot(wolfkde, aes(x=DistFromHighHumanAccess2, y=used, fill=closedFactor)) + stat_smooth(method="glm", method.args = list(family="binomial"), level=0.95) #+ facet_grid(closed~.)

boxplot(DistFromHighHumanAccess2~closed, ylab="Distance from High Human (m)", data=wolfkde)
cor.test(wolfkde$closed, wolfkde$DistFromHighHumanAccess2)
cor.test(wolfkde$closed, wolfkde$Elevation2)

#So yes, you can only get far away from humans evidently in open landcover types (rock / ice) but this isnt that big a problem, based on the correlation coefficient of $\rho$= -0.35. 

#
ggplot(wolfkde, aes(x=Elevation2, y=used, fill=closedFactor)) + stat_smooth(method="glm", method.args = list(family="binomial"), level=0.95) 

#In our final step, lets fit the model now
disthha.cover <-  glm(used~closed + DistFromHighHumanAccess2 + closed*DistFromHighHumanAccess2, data =wolfkde, family= binomial(logit))
summary(disthha.cover)
boxplot(DistFromHighHumanAccess2~closedFactor, ylab="Distance From High Human Access (m)", data=wolfkde)
#So yes, you can only get far away from humans evidently in open landcover types (rock / ice) but this is not that big a problem. 
#Lets try again with distance to human access

ggplot(wolfkde, aes(x=DistFromHumanAccess2, y=used, fill=closedFactor)) + stat_smooth(method="glm", method.args = list(family="binomial"), level=0.95) #+ facet_grid(closed~.)
#Note, there is a bit stronger evidence of an interaction here (the lines cross), which brings us back to our original observation above. 
boxplot(DistFromHumanAccess2~closedFactor, ylab="Distance From High Human Access (m)", data=wolfkde)
distha.cover <-  glm(used~closed + DistFromHumanAccess2 + closed*DistFromHumanAccess2, data =wolfkde, family= binomial(logit))
summary(distha.cover)
#While it is tempting to think that we have succesfully used the interaction between elevation and distance to high human access to 'separate' the confounding between elevation and distance to hha, the problem of collinearity remains problematic for these 2 variables. They are, fundamentally, too correlated at r = 0.53, confounded, and we have too few observaitons of wolves at high elevations far from human access to be able to reliable interpret the interaction here.  Interactions offer a way to 'break' collinearity, but not solve confounding.  I would not be tempted here for risk of introducing a spurious interaction to the model. 

# Model Building and Model Selection
## Cleaning Up Missing Data
#But first we have to clean up msising data for all model selection steps, as missing data that are unbalanced between covariates will result in different degrees of freedom, or, certainty, between models. And models with different degrees of freedom or rows of data cannot be compared directly.  For example, in our dataset we have different NAs in elevation and other variables. 

length(wolfkde$Elevation2)
wolfkde2 <- na.omit(wolfkde)
length(wolfkde2$Elevation2)
#Note there were 252 NA's for some covariates.  Check that there are no other NA's in any of the other fields (on your own for homework).

# Akaike Information Critera
# Manual Model Selection Using AIC
cover <-  glm(used~closed, data =wolfkde2, family= binomial(logit))
## Estimating AIC manually
logLik(cover) ## this is the log likelihood
2*(length(cover$coefficients)) ## this is the number of parameters
#Note that we can calcualte AIC manually using -2* LL + 2*K where LL is logLik and K is # of parameters (without considering the small sample size correction)
-2*as.numeric(logLik(cover))+2*(length(cover$coefficients))
## Note we don't have to do this all manually, in the model str(cover) we see
cover$aic
#Lets use using AIC to select interactions...
distha <-  glm(used~DistFromHumanAccess2, data =wolfkde2, family= binomial(logit))
distha.cover <-  glm(used~closed + DistFromHumanAccess2, data =wolfkde2, family= binomial(logit)) ## Main effects only
disthaXcover <-  glm(used~closed + DistFromHumanAccess2 + closed*DistFromHumanAccess2, data =wolfkde2, family= binomial(logit))
##
AIC(cover, distha, distha.cover, disthaXcover)

#Lets redo with distance to high human access
disthha <-  glm(used~DistFromHighHumanAccess2, data =wolfkde2, family= binomial(logit))
disthha.cover <-  glm(used~closed + DistFromHighHumanAccess2, data =wolfkde2, family= binomial(logit)) ## Main effects only
disthhaXcover <-  glm(used~closed + DistFromHighHumanAccess2 + closed*DistFromHighHumanAccess2, data =wolfkde2, family= binomial(logit))
AIC(cover, disthha, disthha.cover, disthhaXcover)
#Again, , here there is STRONG evidence that model disthhaXcover is much better than the model with the additive main effects. 

## Stepwise Model Selection
# Lets review the full.model again
full.model = glm(used~deer_w2 + elk_w2 +moose_w2 +sheep_w2+goat_w2+Elevation2+DistFromHumanAccess2+DistFromHighHumanAccess2 +closed + closed*DistFromHighHumanAccess2, data =wolfkde2, family= binomial(logit))

#There are two ways to consider stepwise model selection.  Backwards, starting from the most complex 'global' model (full.model) above:
 ## Backwards selection
stepAIC(full.model, direction = "backward")
top.backwards = glm(used ~ deer_w2 + elk_w2 + moose_w2 + sheep_w2 + goat_w2 + Elevation2 + DistFromHumanAccess2 + DistFromHighHumanAccess2, data=wolfkde2,family=binomial(logit))
summary(top.backwards)

#and second, Forward Stepwise model selection, which starts from the null model and proceeds to the maximum # of parameters specified in the full.model
# Forwards selection - First define a NULL model as the starting place
null.model = glm(used~1,data=wolfkde2,family=binomial(logit))
### This time with output from stepAIC supressed 
stepAIC(null.model, scope=list(upper=full.model, lower= null.model),direction="forward")
## lots of output supressed in Rmarkdown
top.forward <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + 
                     moose_w2 + elk_w2 + goat_w2 + DistFromHumanAccess2 + sheep_w2 + 
                     deer_w2 + closed, family = binomial(logit), data = wolfkde2)
summary(top.forward)
#Very similar - identical - to the backwards top model.  Right ! we are done??? are we? Lets screen for collinearity. 
vif(top.forward)
vif(top.backwards)

#But there are a bunch of collinear variables in the model, moose/elk/deer, human/human/human. Basically everything is being retained, not much kicked out. 
#*This is our first experience that you can throw garbage into a model selection algorithm and get garbage out. * Model selection using anything, AIC, should never replace careful consideration of all the variables in the top model, their collinearity, confounding, and interactions. 

#Now what about using AIC to whittle down landcover (with rock Ice as the intercept)??
full.model.landcov = glm(used~ closedConif +modConif+openConif+decid+regen+mixed+herb+shrub+water+burn+alpine, data =wolfkde2, family= binomial(logit))
stepAIC(full.model.landcov, direction = "backward")
#note that AIC helped us trim our landcover categories down from 12 categories to 9, kicking out what - alpine (which is now lumped with Rock/Ice - makes sense), decid and regen - which, if you recall, were both categories with very very few observaitons. This is encouraging, as we hope this is what AIC is supposed to do.

top.model.landcov = glm(used~openConif+modConif+closedConif+mixed+herb+shrub+water+burn, data =wolfkde2, family= binomial(logit))
summary(top.model.landcov)
vif(top.model.landcov)
#Despite the evidence for potential collinearity amongst our categorical variables, even though its harder to discern, lets use this combination of Landcover covariates next as the BEST top model selected using stepAIC later.

#  Model selection using the AICcmodavg package
## Biotic Model List
#Model set 1: Biotic interactions, deer/elk/moose all too correlated to put in the same model, sheep and goat are OK. 
m.biotic <- list()
head(m.biotic)

#lets fit our a-priori list of models 
## Model set 1: Biotic
m.biotic[[1]] <- glm(used ~ 1, family=binomial(logit), data=wolfkde2)
m.biotic[[2]] <- glm(used ~ elk_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[3]] <- glm(used ~ deer_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[4]] <- glm(used ~ moose_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[5]] <- glm(used ~ sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[6]] <- glm(used ~ goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[7]] <- glm(used ~ moose_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[8]] <- glm(used ~ deer_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[9]] <- glm(used ~ elk_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[10]] <- glm(used ~ elk_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[11]] <- glm(used ~ deer_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[12]] <- glm(used ~ moose_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[13]] <- glm(used ~ sheep_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[14]] <- glm(used ~ DistFromHighHumanAccess2, family=binomial(logit), data=wolfkde2)
m.biotic[[15]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[16]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[17]] <- glm(used ~ DistFromHighHumanAccess2+sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[18]] <- glm(used ~ DistFromHighHumanAccess2+goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[19]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[20]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[21]] <- glm(used ~ DistFromHighHumanAccess2+elk_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[22]] <- glm(used ~ DistFromHighHumanAccess2+elk_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[23]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[24]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[25]] <- glm(used ~ DistFromHighHumanAccess2+sheep_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[26]] <- glm(used ~ DistFromHighHumanAccess2, family=binomial(logit), data=wolfkde2)
m.biotic[[27]] <- glm(used ~ DistFromHighHumanAccess2+deer_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[28]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[29]] <- glm(used ~ DistFromHighHumanAccess2+sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[30]] <- glm(used ~ DistFromHighHumanAccess2+goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[31]] <- glm(used ~ DistFromHighHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[32]] <- glm(used ~ DistFromHumanAccess2, family=binomial(logit), data=wolfkde2)
m.biotic[[33]] <- glm(used ~ DistFromHumanAccess2+deer_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[34]] <- glm(used ~ DistFromHumanAccess2+moose_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[35]] <- glm(used ~ DistFromHumanAccess2+sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[36]] <- glm(used ~ DistFromHumanAccess2+goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[37]] <- glm(used ~ DistFromHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[38]] <- glm(used ~ DistFromHumanAccess2+deer_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[39]] <- glm(used ~ DistFromHumanAccess2+elk_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[40]] <- glm(used ~ DistFromHumanAccess2+elk_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[41]] <- glm(used ~ DistFromHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[42]] <- glm(used ~ DistFromHumanAccess2+moose_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[43]] <- glm(used ~ DistFromHumanAccess2+sheep_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[44]] <- glm(used ~ DistFromHumanAccess2, family=binomial(logit), data=wolfkde2)
m.biotic[[45]] <- glm(used ~ DistFromHumanAccess2+deer_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[46]] <- glm(used ~ DistFromHumanAccess2+moose_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[47]] <- glm(used ~ DistFromHumanAccess2+sheep_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[48]] <- glm(used ~ DistFromHumanAccess2+goat_w2, family=binomial(logit), data=wolfkde2)
m.biotic[[49]] <- glm(used ~ DistFromHumanAccess2+moose_w2 + sheep_w2, family=binomial(logit), data=wolfkde2)
                   
## then name our models .
## note you can name your models with a command like this
# model.names <-  ("null", "disthha", "distacc", "sheepwi", "goatwin", "elkwint", "moosewin", "deerwin") but in this case there were 49 models
model.names.biotic <-c("m0","m1","m2","m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","m13","m14","m15","m16","m17","m18","m19","m20","m21","m22","m23","m24","m25","m26","m27","m28","m29","m30","m31","m32","m33","m34","m35","m36","m37","m38","m39","m40","m41","m42","m43","m44", "m45","m46","m47","m48")
model.names.biotic <-1:49
aictab(cand.set = m.biotic, modnames = model.names.biotic)
## OK so the top model was model 41
top.biotic <- glm(used ~ DistFromHumanAccess2+deer_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
summary(top.biotic)
vif(top.biotic)
#So - not too badly collinear. 
#And the 2nd ranked top biotic model was  model 40
second.biotic <- glm(used ~ DistFromHumanAccess2+elk_w2 + goat_w2, family=binomial(logit), data=wolfkde2)
summary(second.biotic)
vif(second.biotic)

## Model set 2: Environmental Covariates Only
m.env <- list()
head(m.env)

## Model set 1: Biotic
m.env[[1]] <- glm(used ~ 1, family=binomial(logit), data=wolfkde2)
m.env[[2]] <- glm(used ~ Elevation2, family=binomial(logit), data=wolfkde2)
m.env[[3]] <- glm(used ~ DistFromHighHumanAccess2, family=binomial(logit), data=wolfkde2)
m.env[[4]] <- glm(used ~ DistFromHumanAccess2, family=binomial(logit), data=wolfkde2)
m.env[[5]] <- glm(used ~ openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde2)
m.env[[6]] <- glm(used ~ Elevation2 + DistFromHumanAccess2, family=binomial(logit), data=wolfkde2)
m.env[[7]] <- glm(used ~ DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde2)
m.env[[8]] <- glm(used ~ DistFromHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde2)
m.env[[9]] <- glm(used ~ Elevation2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde2)
m.env[[10]] <- glm(used ~ Elevation2 + DistFromHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde2)
m.env[[11]] <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde2)
m.env[[12]] <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + closed + closed*DistFromHighHumanAccess2, family=binomial(logit), data=wolfkde2)
m.env[[13]] <- glm(used ~ Elevation2 + DistFromHumanAccess2 + closed + closed*DistFromHumanAccess2, family=binomial(logit), data=wolfkde2)
m.env[[14]] <- glm(used ~ DistFromHighHumanAccess2 + closed + closed*DistFromHighHumanAccess2, family=binomial(logit), data=wolfkde2)
m.env[[15]] <- glm(used ~ DistFromHumanAccess2 + closed + closed*DistFromHumanAccess2, family=binomial(logit), data=wolfkde2)

model.names.env <-1:15
aictab(cand.set = m.env, modnames = model.names.env)

#OK - top model is model 11
top.env <- glm(used ~ Elevation2 + DistFromHighHumanAccess2 + openConif+modConif+closedConif+mixed+herb+shrub+water+burn, family=binomial(logit), data=wolfkde2)
summary(top.env)
vif(top.env)

#Now - which 'set' of covariates is best? Env? or Biotic?
AIC(top.env, top.biotic)

## Environmental model HANDS DOWN. 
## now go back and compare 'top' model to top model selected by AIC
AIC(top.forward, top.biotic, second.biotic, top.env)

#This reveals that all model selection methods, especially stepwise, will overfit models and does not penalize for collinearity.  This is a crucial lesson from today. *Model selection methods will not screen out collinear variables for you!*


# Model Selection using the MuMIn Package

# re-run FULL logistic regression model
top.forward = glm(used ~ deer_w2 + elk_w2 + moose_w2 + sheep_w2 + goat_w2 + Elevation2 + DistFromHumanAccess2 + DistFromHighHumanAccess2 + closed + DistFromHighHumanAccess2*closed, data=wolfkde2,family=binomial(logit), na.action ="na.fail")
summary(top.forward)
#install and load MuMIn package
require(MuMIn)

#use dredge function to get all possible models
x1<-dredge(top.forward)

#x1 looks at ALL possible model combinations here, there are over 1000 models fit! 10! models = ? models.  Dredge has fit XX models in total out of this candidate set of 19 candidate variables.
head(x1, n = 10) ## only shows top 10 models fit
plot(x1)
#get top models with AICc <2
top.models<-get.models(x1, subset=delta<2)
#model average covariate effects
x6<-model.avg(top.models)
summary(x6)
#This is a model averaged set of coefficients across the top model set. Given we did not have much model selection uncertainty here, this model does not differ that much from our top model(s). 

#Next lets do landcover using dredge to whittle down our landcover model list. This is a COMMON approach I take with trying to reduce landcover. 
top.dredge.lc = glm(used ~ openConif+modConif+closedConif+mixed+herb+shrub+water+burn+decid+regen+alpine, data=wolfkde2,family=binomial(logit), na.action ="na.fail")
x2<-dredge(top.dredge.lc)
head(x2, n=10)
top.lc <- glm(used ~ openConif+modConif+closedConif+mixed+herb+shrub+water+burn, data=wolfkde2,family=binomial(logit))
summary(top.lc)
#compare to full landcover model
AIC(top.lc, top.dredge.lc)
#So the top landcover model has 2 fewer landcover types, and a dAIC of 3.6 'better' than the full model.  

## Manipulating Coefficients from Dredge
coefficients(top.models[[1]])
top.model.coef <- lapply(top.models, coefficients)
#str(top.model.coef)
str(top.model.coef) makes a horrendous list. You can collapse those into a single data frame by doing something like this with `ldply`:
require(plyr)
ldply(top.models, function(l) as.data.frame(t(coefficients(l))))
summary(top.models[[1]])$coefficients
tidyList1<- ldply(top.models, function(l) as.data.frame(t(summary(l)$coefficients[,4])))
head(tidyList1)

#Next we will use the `broom` package, which has a `tidy` command which might automate everything for you.  Check this out:
require(broom)
tidy(top.models[[1]])
tidyList2 <- ldply(top.models, tidy)
head(tidyList2)
#What you actually want is confidence intervals.  You can get those out of the `tidy` output with a little piping %>% (with magrittr, which is included in tidyverse)
(CI.table <- ldply(top.models, tidy) %>%  mutate(CI.low = estimate - 2*std.error, CI.high = estimate + 2*std.error))
#This thing you can quickly ggplot, which is nice:
ggplot(mutate(CI.table, model = .id), aes(model, estimate)) + 
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high)) + 
  facet_wrap(.~term, scales = "free_y") + theme(axis.text.x = element_text(angle = 90))

#  Model Selection using BIC
## First manually
AIC(top.forward, top.biotic, second.biotic, top.env)
BIC(top.forward, top.biotic, second.biotic, top.env)
#OK - so not much difference in top models using BIC and AIC in this dataset.

## BIC with Dredge 
x1.bic<-dredge(top.forward, rank=BIC) ## note this now ranks using BIC
plot(x1.bic)
## x1.bic - look at all 

head(x1.bic, n = 10) ## only shows top 10 models fit
# lets compare the top model from AIC and BIC
head(x1.bic, n = 1) ## only shows top 1 models fit with BIC
head(x1, n = 1) ## only shows top 1 models fit with AIC
#So AIC is overfitting here potentially, selecting a model with 11 parameters versus 7 parameters with BIC. Take note.  This is a theme. 

#get top models with BIC <2
top.models.bic<-get.models(x1.bic, subset=delta<2)
top.models.bic 
#Note there is only 1 top model here using BIC, thus, there is no need to model average covariate effects using this command. 
x.top.bic<-model.avg(top.models.bic) ## only 1 top model, so this doesnt work

## Lets run the 'top' model selected using BIC for next week
top.model.bic = glm(used ~ DistFromHighHumanAccess2 + DistFromHumanAccess2+Elevation2+elk_w2+goat_w2+moose_w2, data=wolfkde2,family=binomial(logit), na.action ="na.fail")
summary(top.model.bic)
## compare to top AIC model
summary(top.forward)

# Caterpillar plots of coefficients
# run logistic regression model
summary(full.model)
B<-summary(full.model)$coefficient[1:length(summary(full.model)$coefficient[,1]),1]
#create margin of error (ME) for 95% CI
ME <- summary(full.model)$coefficient[1:length(summary(full.model)$coefficient[,1]),2]*1.96
lower<-B - ME
upper<-B + ME

# bundle into data frame
logisData<-data.frame(B, lower, upper, names(summary(full.model)$coefficient[,2]))
names(logisData) <- c("Coefficient", "lower.ci", "upper.ci", "Variable")
levels(logisData$Variable)[1] <- "Intercept"
#logisData$Variable <- relevel(logisData$Variable, ref="Intercept")
## Lets make nicer labels for graphing of the covariate oders that I pulled out of logisData
figLabels = c("B0", "Closed", "Deer", "DHHA", "D:C", "DHA", "Elev", "Elk", "Goat", "Moose", "Sheep")
pd <- position_dodge(0.6) # move them .05 to the left and right
x1<-ggplot(data=logisData, aes(x=Variable,y=Coefficient)) +
  geom_errorbar(data=logisData,aes(ymin=lower.ci, ymax=upper.ci), width=.4,position=pd,size=1) + geom_point(size=3, col="blue") 

p6<-x1+theme(axis.text.y = element_text(size=14, family="Times"),axis.text.x = element_text(size=14, family="Times", angle = 90, vjust = 0.5),text = element_text(size=16, family="Times"),axis.title.x=element_text(size=16, family="Times"),axis.title.y=element_text(size=16, family="Times",vjust=1))

p7<-p6+theme(axis.line.x = element_line(color="black", size = 0.25),
             axis.line.y = element_line(color="black", size = 0.25),legend.title=element_blank(),legend.text=element_text(size=16, family="Times"))+ylab("Estimate") + xlab("Coefficient") + scale_x_discrete(labels = figLabels)

p7

tiff("coefPlot.tiff", res=600, compression = "lzw", height=5, width=7, units="in")
p7
dev.off()

##  Caterpillar Plots with the GGally Package   
ggcoef(full.model)

#Note how the intercept is off the charts, so lets remove that, and play around with some other options such as sorting, exponentiating the coefficients which then displays them as Odds ratio's from a logistic regression model, etc. 

ggcoef(full.model, exclude_intercept = TRUE, exponentiate = FALSE, sort = "ascending")

ggcoef(full.model, exclude_intercept = TRUE, exponentiate = TRUE, sort = "ascending")

# Variable reduction using PCA - Principle Components Analysis
head(wolfkde2)
pcawolf <-princomp(na.omit(wolfkde2[1:9]), cor=TRUE)
summary(pcawolf)
loadings(pcawolf)
plot(pcawolf, type="lines")
biplot(pcawolf, xlim =c(-0.06, 0.04))

wolfkde2$Comp.1 <- -0.406*wolfkde2$deer_w2 - 0.370*wolfkde2$moose_w2 - 0.402*wolfkde2$elk_w2 +0.182*wolfkde2$goat_w2 - 0.415*wolfkde2$wolf_w2 + 0.408*wolfkde2$Elevation2 + 0.318*wolfkde2$DistFromHumanAccess2 + 0.233*wolfkde2$DistFromHighHumanAccess2

wolf_comp1 <- glm(used ~ Comp.1, family=binomial (logit), data=wolfkde2)
wolfkde2$fitted1 <- fitted(wolf_comp1)
hist(wolfkde2$fitted1)
plot(wolfkde2$fitted1, wolfkde2$Comp.1)

figPCA <- ggplot(wolfkde2, aes(x=Comp.1, y=used)) + stat_smooth(method="glm", method.args = list(family="binomial"))
x.axis = "-0.41*deer - 0.37*moose - 0.4*elk +0.18*goat - 0.42*wolf + 0.41*Elev + 0.32*DistHum + 0.23*DistHighHum"
figPCA2 <- figPCA + xlab(x.axis)
figPCA2