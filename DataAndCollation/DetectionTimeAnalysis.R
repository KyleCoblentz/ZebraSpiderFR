###############################################################################
### Detection time -- Stage and Sex
###############################################################################

### load packages 

library(ggplot2); library(dplyr); library(cowplot); library(chron); library(tidyr); library(brms); library(RColorBrewer)

### load data

detection <- read.csv('DetectionTimeData.csv')

### use chron package for times to calculate differences in time

detection$TimeStart <- chron(times = detection$TimeStart)

detection$TimeStop <- chron(times = detection$TimeStop)

detection <- detection %>% mutate(TimeDiff = as.numeric(TimeStop - TimeStart))

detection <- detection %>%  mutate(TempStart = replace(TempStart, TempUnit == 'F', (5/9)*(TempStart[TempUnit == 'F']-32)))  

detection <- detection %>% filter(PreySize != is.na(PreySize))

################################################################################
### fit regression in Bayesian framework
################################################################################


fit <- brm(formula = log(TimeDiff) ~ log(SpiderLength) + PreySize + log(TempStart),
               data = detection, seed = 666)

summary(fit)

### model coefficients
# B0 = 0.10
# B1 = -0.71
# B2 = -0.76
# B3 = -1.37
# B4 = -0.67

























