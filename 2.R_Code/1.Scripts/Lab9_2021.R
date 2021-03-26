#Title: "WILD 562: Intro to Step Selection Functions"
#author: "Mark Hebblewhite and Eric Palm"
# Introduction to Step Selection Functions

# Loading necessary packages
require(raster)
require(rgdal)
require(amt)
require(mapview)
require(tidyverse)
require(survival)
require(sjPlot)
require(lme4)
#

# Biased-Correlated Random Walks

#_From Fagan and Calabrese - Bulletin of the Ecological Society of America 2014_

# The origins of Step Selection Functions stem from the earliest days of animal ecology with Skellam's classic 1950's paper.  But the field of movement ecology didnt really get started until field ecologists started working with theoretical physicists and ecologists to understand how to conceptualize movement based on ideal gas law theories.  These seminal studies were summarized by Peter Turchin in his landmark book, Quantitative Analysis of Animal Movement (Turchin 1999).   
# 
# More than 30 years ago, an early, sturdy bridge between field data and spatial ecological theory was built when the article “Analyzing insect movement as a correlated random walk” was published in Oecologia. This paper, which represented a collaboration between ecologist Peter Kareiva and mathematician Nanako Shigesada, is a milestone along the Paper Trail because it marks a critical link between the abstract world of ecological theory and the hands-on way in which ecologists actually collect data on individual animals.
# 
# This correlated random walk model was comprised of steps and turns, and Kareiva and Shigesada showed how one could estimate these distributions, and, make them functions of spatial or temporal covariates through field data.  The biased correlated random walk emerged, and represents the cornerstone of the step dplyr::selection function concept. And it links the mechanistic movement models of Moorcroft, Lewis and Barnett to field data approaches commonly collected with GPS data. 
# 
# In this first excercise, we will explore how different 'parameters' of movement (step, turns) and bias towards a centroid influence the spatial pattern of movement.  In essence, the BCRW is the driver of the movement kernel distribution in the iSSF models we will use from the package amt. 
# 
# First, we make a function that draws random movements based on 3 parameters, a, b, rho (the degree of correaltion), and an attraction to a home range activity center.  a and b are parameters of the step lenght distribution, fit as a Weibull distribution. Larger values of a or b represent more or less longer step lenghts.  Rho is the degree of directional persistence or 'bias' in the correlation in direction between steps, and the attraction is equivalent to the mathematical advection term in Moorcroft and Barnett.  
# 
# Here, we will compare just 3 types of fits to explore, but I encourage you to play around with the paramters on your own to get a feel for unbiased and biased correlated random walks. 
# crw}

#### Correlated Random Walks
BCRW <- function(a = 2, b = 1, rho = 0.7, Z.center = 0, attraction = 0.5, n = 50, Z0 = 0){
  require(CircStats)
  
  Z <- c(Z0, rep(NA,n-1))
  phi <- runif(1, -pi, pi)
  for(i in 2:n)
  {
    # relative orientation to center of attraction
    chi <- Arg(Z.center - Z[i-1])
    dphi <- chi-phi
    if(abs(chi - phi) > pi) dphi <- (chi - phi) - pi
    
    # adjust the location 
    location <- phi + attraction * dphi
    
    # pick a new absolute direction ... but MUST BE BETWEEN -pi and pi
    phi <- rwrpcauchy(1, location, rho) - 2*pi
    if(phi > pi) phi <- phi-2*pi
    if(phi < -pi) phi <- phi+2*pi
    
    Z[i] <- Z[i-1] + complex(arg = phi, mod = rweibull(1, a, b))
  }
  return(Z)
}

BCRW(a = 2, b = 1, rho = 0, Z.center = 10, attraction = 0.25, n = 2000) %>% 
  plot(type="o", asp=1, pch = 21, bg= grey(seq(0,1,length = 2000)),
       main = "a = 2, b = 1, rho = 0.2, attraction = 0")

BCRW(a = 2, b = 1, rho = 0.5, Z.center = 10, attraction = 0.25, n = 2000) %>% 
  plot(type="o", asp=1, pch = 21, bg= grey(seq(0,1,length = 200)),
       main = "a = 2, b = 1, rho = 0.7, attraction = 0.5")

BCRW(a = 2, b = 1, rho = 0.7, Z.center = 10, attraction = 0.25, n = 2000) %>% 
  plot(type="o", asp=1, pch = 21, bg= grey(seq(0,1,length = 200 )),
       main = "a = 2, b = 1, rho = 0.7, attraction = 0.8")

#
# Loading and importing data

#, warning = F}
elev<-raster("1.Data/Lab9_data/elev.tif")
slope<-raster("1.Data/Lab9_data/slope.tif")
d_human<-raster("1.Data/Lab9_data/d_human.tif")
d_high_human <-raster("1.Data/Lab9_data/d_high_human.tif")
habitat_stack<-stack(elev, slope, d_human, d_high_human)
plot(habitat_stack)
habitat_stack@layers


#And then our elk telemetry data:
  #}
elk_df <- read_csv("1.Data/Lab9_data/elk_df.csv")
#
#For some reason the `read_csv` function didn't parse the timestamp column as "datetime" format, so we'll manually convert it to POSIXct format, which is the date-time format that `amt` likes:
elk_df$timestamp <- as.POSIXct(elk_df$timestamp, format = "%m/%d/%y %H:%M")
elk_df
#
# Data visualization and exploration 
elk_sp <- SpatialPointsDataFrame(coords = elk_df[, c("lon","lat")], data = elk_df, proj4string = sp::CRS("+init=epsg:4326"))
elk_sp_UTM <- spTransform(elk_sp, habitat_stack@crs)

# map
mapview(elk_sp, zcol="id", legend = TRUE, cex=5, lwd=2, map.type = "Esri.DeLorme")
# We can overlay our elk telemetry locations on the elevation raster:
plot(habitat_stack$elev)
points(elk_sp_UTM, pch=20, col=c("blue", "red", "green", "purple", "navy", "darkgreen")[as.factor(elk_sp_UTM$id)])

#To get an idea of how many locations we have per individual:
table(elk_df$id)

# Creating and nesting an `amt` track

#}
elk_trk <- amt::make_track(elk_df, .x=lon, .y=lat, .t=timestamp, id=id, crs = sp::CRS("+init=epsg:4326")) %>%
  amt::transform_coords(habitat_stack@crs)
elk_trk
#
elk_trk_nested <-
  elk_trk %>% 
  nest(-id)
#
head(elk_trk_nested$data[[1]])

# Check for duplicated time stamps and complete cases.
all(complete.cases(elk_trk))
any(duplicated(elk_trk$ts))

#time of day
elk_trk <- time_of_day(elk_trk)
head(elk_trk)
table(elk_trk$tod_, elk_trk$id)

## `amt` Data Summaries and Visualization
  elk_trk %>% 
  nest(-id) %>% 
  mutate(lags = map(data, summarize_sampling_rate))
## lets take a look at some of the summary statistics
elk_trk2 <-elk_trk_nested %>% 
  mutate(lags = map(data, summarize_sampling_rate))
elk_trk2$lags

# we can now "unnest" it to see our sampling rate summary.  Let's keep the "id" column too so we can see our animal ids.
elk_trk %>% 
  nest(data = -id) %>% 
  mutate(lags = map(data, summarize_sampling_rate)) %>%
  dplyr::select(-data) %>%
  unnest(lags)
#

# amt summary statistics
elk_trk_stats <- 
  elk_trk %>% 
  nest(data = -id) %>% 
  mutate(speed = map(data, speed),
         step_length = map(data, step_lengths),
         tod_ = map(data, time_of_day), 
         turn_angle = map(data, ~direction_rel(.) %>% as_degree(.))) %>%   dplyr::select(-data) %>% 
  unchop(speed:turn_angle)

# You can see there are NAs in "speed", "step_length" and "turn_angle":
summary(elk_trk_stats)

#Something strange happed with the summary statistics of time of day. 
head(elk_trk_stats)
str(elk_trk_stats$tod_)

## Summary plots
#, warning=F}
elk_trk_stats %>% 
  ggplot(., aes(x = turn_angle, fill=id)) +  
  geom_histogram(breaks = seq(-180,180, by=10))+
  theme_classic() + 
  ylab("Count") + 
  ggtitle("Relative turn angles") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=60),
                     labels = seq(-180, 180, by=60)) +
  facet_wrap(~id, scales="free") +
  theme(legend.position = "none")
#

#, warning=F}
elk_trk_stats %>% 
  ggplot(., aes(x = turn_angle)) +  
  geom_histogram(breaks = seq(-180,180, by=10))+
  theme_classic() + 
  ylab("Count") + 
  ggtitle("Relative turn angles") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=60),
                     labels = seq(-180, 180, by=60))

# Polar coordinate plots
elk_trk_stats %>% 
  ggplot(., aes(x = turn_angle, y = ..density..)) +  
  geom_histogram(breaks = seq(-180,180, by=20))+
  coord_polar(start = 0)+
  theme_classic() + 
  ylab("Count") + 
  ggtitle("Relative turn angles") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=60), labels = seq(-180, 180, by=60))


#} Polar by individuals
elk_trk_stats %>% 
  ggplot(., aes(x = turn_angle, y = ..density.., fill = id)) +  
  geom_histogram(breaks = seq(-180,180, by=20))+
  coord_polar(start = 0)+
  theme_classic() + 
  ylab("Count") + 
  ggtitle("Relative turn angles") + 
  scale_x_continuous("", limits = c(-180, 180), breaks = seq(-180, 180, by=60), labels = seq(-180, 180, by=60)) + facet_wrap( ~ id)
#


#Next we can plot histograms of step lengths faceted by individual.
#, warning=F}
elk_trk_stats %>% 
  ggplot(., aes(x = step_length, fill=id)) +  
  geom_histogram(breaks = seq(0,4000, by=250))+
  theme_classic() + 
  ylab("Count") + 
  ggtitle("Step lengths (m)") + 
  scale_x_continuous("", limits = c(0, 4000), breaks = seq(0, 4000, by=1000),
                     labels = seq(0, 4000, by=1000)) +
  facet_wrap(~id, scales="free") +
  theme(legend.position = "none")
#

#Next we can plot histograms of step lengths faceted by individual.
#, warning=F}
elk_trk_stats %>% 
  ggplot(., aes(x = log(step_length), fill=id)) +  
  geom_histogram(breaks = seq(0,10, by=0.5))+
  theme_classic() + 
  ylab("Count") + 
  ggtitle("Step lengths (m)") + 
  scale_x_continuous("", limits = c(0, 10), breaks = seq(0, 10, by=1),
                     labels = seq(0, 10, by=1)) +
  facet_wrap(~id, scales="free") +
  theme(legend.position = "none")
#

#So, movement rate is approximately log-normal. Remember this - there is never ever ever anything normal about movement rate parameters. 

#How about we see if animals move faster during the day versus at night.
#, warning=F}
ggplot(elk_trk_stats, aes(x = tod_[[4]], y = speed, fill=tod_[[4]])) + 
  geom_violin() +
  theme_bw() +
  facet_wrap(~id, scales = "free") +
  theme(legend.position = "none") +
  ylab("speed (m/s)") +
  xlab("time of day")
#

# Log of speed
ggplot(elk_trk_stats, aes(x = tod_[[4]], y = log(speed), fill=tod_[[4]])) + 
  geom_violin() +
  theme_bw() +
  facet_wrap(~id, scales = "free") +
  theme(legend.position = "none") +
  ylab("log(speed)") +
  xlab("time of day")
# Seems reasonable that they move a bit faster during the day.

# Prepare SSF data frame by individual

elk_trk_id <- 
  elk_df %>% 
  nest(-id) %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, lon, lat, timestamp, crs = sp::CRS("+init=epsg:4326")) %>%
      transform_coords(habitat_stack@crs)
  }))
#

#Now we've made six tracks, one for each individual.
elk_trk_id
elk_trk_id$trk[[1]]
#

 
## Create available steps and extract covariates
ssf_2_hr <- elk_trk_id %>%
  mutate(steps_2_hr = map(trk, function(x) {
    x %>%
      track_resample(rate = minutes(120), tolerance = minutes(20)) %>%
      filter_min_n_burst(min_n = 3) %>%
      steps_by_burst(diff_time_units = "hours") %>%
      random_steps(n = 3) %>% 
      extract_covariates(habitat_stack, where = "end")
  })) %>%
  dplyr::select(id, steps_2_hr) %>%
  unnest() 
head(ssf_2_hr)
#
print(ssf_2_hr, width=Inf)
#There seems to be a lot of zeros in the "d_high_human" column, which makes me think something might not be right with that layer, so I'm just not going to include it in our simple model below.  
                                                                     #Visualizing the SSF sampling for a zoomed in portion of our study area:
ggplot(ssf_2_hr, aes(x=x2_, y= y2_, colour = case_)) + geom_point() + geom_path()+ xlim(550000, 560000) + ylim(5700000, 5710000)
                                                                      
## Univariate Plotting of Used versus Available Points
ggplot(ssf_2_hr, aes(x = case_, y = slope, fill=case_)) + geom_violin() +   theme_bw() +
  theme(legend.position = "none") +
  ylab("slope (m)") +
  xlab("")
#
# Elevation
ggplot(ssf_2_hr, aes(x = case_, y = elev, fill=case_)) + 
  geom_violin() +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("slope (m)") +
  xlab("")
#

# human
ggplot(ssf_2_hr, aes(x = case_, y = d_human, fill=case_)) + 
  geom_violin() +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("slope (m)") +
  xlab("")
#

# Running an SSF model in `amt`
ssf_2_hr$stratum <- paste(ssf_2_hr$id, ssf_2_hr$burst_, ssf_2_hr$step_id_)
head(ssf_2_hr$stratum)

#}
ssf_2_hr_raw <-
  ssf_2_hr %>%
  dplyr::select(id, case_, t2_, elev, slope, d_human, d_high_human, step_id_, stratum)


#Then we can scale and center our variables 

ssf_2_hr_scaled <-
  ssf_2_hr_raw %>%
  mutate(
    elev = as.numeric(scale(elev)),
    slope = as.numeric(scale(slope)),
    d_human = as.numeric(scale(d_human)),
    d_high_human = as.numeric(scale(d_high_human))
  )

## Fitting SSFs with clogit in R
`?clogit`

ssf_model <- clogit(case_ ~ elev + slope + d_human + strata(stratum) + cluster(id), method = "approximate", data = ssf_2_hr_scaled)

summary(ssf_model)
#
plot_model(ssf_model, title="SSF Coefficients", transform = NULL)
#

## Comparing to a Naive Logistic Regression
glmm_model = glmer(case_ ~ elev + slope + d_human + (1|id), data=ssf_2_hr_scaled, family=binomial(link="logit"))
summary(glmm_model)
plot_model(ssf_model, transform = NULL)
plot_model(glmm_model, transform = NULL)
coef(ssf_model)
fixef(glmm_model)


## Interpreting Clogit Models

  exp(coef(ssf_model))
#

## Predicting
ssf_2_hr_scaled$predSSF <- predict(ssf_model, type = "expected")
hist(ssf_2_hr_scaled$predSSF)
plot(ssf_2_hr_scaled$elev, ssf_2_hr_scaled$predSSF)
ggplot(ssf_2_hr_scaled, aes(x=elev, y = predSSF, colour = id)) + stat_smooth(method="glm", method.args = list(family="binomial"))

#Note the Y axis here is the relative predicted probabilitiy of selection 
ggplot(ssf_2_hr_scaled, aes(x=slope, y = predSSF, colour = id)) + stat_smooth(method="glm", method.args = list(family="binomial"))

ggplot(ssf_2_hr_scaled, aes(x=d_human, y = predSSF, colour = id)) + stat_smooth(method="glm", method.args = list(family="binomial"))


## Model selection in SSF models
ssf_model1 <- clogit(case_ ~ elev + strata(stratum) + cluster(id), method = "approximate", data = ssf_2_hr_scaled)

ssf_model2 <- clogit(case_ ~ elev + slope + strata(stratum) + cluster(id), method = "approximate", data = ssf_2_hr_scaled)

ssf_model3 <- clogit(case_ ~ elev + d_human + strata(stratum) + cluster(id), method = "approximate", data = ssf_2_hr_scaled)

ssf_model4 <- clogit(case_ ~ d_human + strata(stratum) + cluster(id), method = "approximate", data = ssf_2_hr_scaled)

AIC(ssf_model, ssf_model1, ssf_model2, ssf_model3, ssf_model4)
#
AIC(glmm_model, ssf_model)
