#######################################
##
##       Leave-one-out cross
##       validation for GLMM
##
##       Journal: http://onlinelibrary.wiley.com/doi/10.1890/10-0751.1/full
##       Code: http://www.esapubs.org/archive/ecol/E092/051/appendix-C.htm
##
#######################################

# Packages
require(lme4)

#### START HERE ###
#======================================
#=         Plotting Function          =
#======================================
### Change xlim and ylim as needed to 
### capture all values as described  
### across n plots

gof<-function(fits, obs, title)        # plotting goodness-of-fit plot for a binomial model
{
  Precision <- sum((fits-obs)^2)
  Bias <- sum(fits-obs)
  o <- cbind(Precision, Bias)
  smoothScatter(fits,obs, main=title,  
                xlab="Standardized estimates (n-1)", ylab="Standardized estimates (with-held animal)", nrpoints=500
                ,xlim=c(1.0e-05, 9e-05), ylim=c(0,0.0003) ### can use this to put all graphs on same scale
  )
  abline(0,1)
  return(o)
}

gof2 <- function(fits, obs)        # plotting goodness-of-fit plot for a binomial model
{
  o <- list()
  o$values <- cbind(expected = fits, observed = obs)
  Precision <- sum((fits-obs)^2)
  Bias <- sum(fits-obs)
  o$pb <- cbind(Precision, Bias)
  return(o)
}

tableLOO<- function(looo) {
 o <- lapply(looo, function(x) x$pb)
 ou <- c()
 for (l in 1:length(o)) {
   ou <- rbind(ou, o[[l]])
 }
 row.names(ou) <- names(o)
 return(ou)
}

plotLOO <- function(looo) {
  xmin <- min(unlist(lapply(looo, function (x) min(x$values[,1]))))
  xmax <- max(unlist(lapply(looo, function (x) max(x$values[,1]))))
  xlim <- c(xmin, xmax)
  
  ymin <- min(unlist(lapply(looo, function (x) min(x$values[,2]))))
  ymax <- max(unlist(lapply(looo, function (x) max(x$values[,2]))))
  ylim <- c(ymin, ymax)
  
  for (l in 1:length(looo)) {
    t <- looo[[l]]$values
    smoothScatter(t[,1],t[,2], main=paste(names(looo)[l], 'with-held'),  
                  xlab="Standardized estimates (n-1)", ylab="Standardized estimates (with-held animal)", nrpoints=500,
                  xlim=xlim, ylim=ylim)
    abline(0,1)                  
    par(ask=T)
  }
  par(ask=F)
}

# 2.LEAVE-ONE-OUT VALIDATION ===================================================
loo <- function(dat, formMix, formNonMix, thin=10, rE = 'ID') 
{
  # Reduced sample size
  # Used to increase speed of processing time
  aces<-length(dat$Used[dat$Used>0])
  ids<-c(seq(1,aces, thin), (aces+1):length(dat[,1]))
  subda<-dat[ids,]   # creates a sparser data set for fitting
  
  ### Dealing with factors 
  ### Creates dummy codes for each category
  ### Not using this currently; see cited paper above 
  ### for details
  #lkup<-list(1, 2, 3, 4, 5, 6) ### alpha-numeric representation of categories
  #labs<-c("glv", "dlv", "gbv", "sbv", "lbv", "blv")
  
  #for(k in 1:length(lkup))
  #{
  #  ids<-which(subda[,21] %in% lkup[[k]])
  #  colNew<-rep(0, length(subda[,21]))
  #  colNew[ids]<-1
  #  ifelse (k==1, cols<-data.frame(colNew), cols<-cbind(cols, colNew))
  #}
  #names(cols)<-labs
  #subda<-cbind(subda,cols)
  
  # fitting the model with one missing animal
  uE <- unique(subda[, rE])
  out <- list()
  for (id in uE) {
    datless <- subset(subda, ID!=id)
    datonly <- subset(subda, ID==id)
    
    modTrain <- glmer(formMix, datless, family=binomial, verbose=F)
    bestLev <- glm(formNonMix, datonly, family=binomial) ### best GLM model for with-held animal
    bestPredict <- predict.glm(bestLev, type="link")
    
    preds <- model.matrix(terms(modTrain), datonly) %*% fixef(modTrain)
    #rec <- gof(plogis(preds)/sum(plogis(preds)), plogis(bestPredict)/sum(plogis(bestPredict)), 
    #           paste(id, "with-held"))
    rec <- gof2(plogis(preds)/sum(plogis(preds)), plogis(bestPredict)/sum(plogis(bestPredict)))
    
    rec$pb <- data.frame(rec$pb, Correlation = cor(preds, bestPredict))
    out[[as.character(id)]] <- rec
  }
  return(out)
}



