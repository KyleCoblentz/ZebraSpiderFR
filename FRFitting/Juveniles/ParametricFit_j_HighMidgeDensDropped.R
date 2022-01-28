################################################################################
### Juveniles -- Only midge densities less than 4
################################################################################

### load packages

library(dplyr); library(ggplot2); library(cowplot); library(rstan); library(loo)

### set some options for Stan

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

### load data

juveniles <- read.csv('FeedingData_Manipulated_j.csv', stringsAsFactors = TRUE)

### filter data

juveniles <- juveniles %>% filter(NumberMidges > 0 & NumberMidges < 4)

juveniles <- juveniles %>% mutate(FeedingRate = ProportionFeed/(DetectionEst*24))


################################################################################
### Type IV FR
################################################################################

quadtemp_sepint_typeiv_data <- list(N = nrow(juveniles),
                                    y = juveniles$NumberFeed,
                                    n_t = juveniles$NumberPred,
                                    d = juveniles$DetectionEst*24,
                                    R = juveniles$NumberMidges,
                                    C = juveniles$JuvenileDens,
                                    F = juveniles$FemaleDens,
                                    M = juveniles$MaleDens,
                                    Temp = juveniles$TempStart,
                                    Temp2 = juveniles$TempStart^2)

quadtemp_sepint_typeiv_fit <- stan(file = 'Juvenile_tempquad_BDall_TypeIV.stan', data = quadtemp_sepint_typeiv_data)

samples_quadtemp_sepint_typeiv <- rstan::extract(quadtemp_sepint_typeiv_fit)

### get summary of model fit

summary_quadtemp_sepint_typeiv <- summary(quadtemp_sepint_typeiv_fit)

summary_quadtemp_sepint_typeiv$summary[1:7,]

quantile(samples_quadtemp_sepint_typeiv$ca, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint_typeiv$ba, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint_typeiv$qa, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint_typeiv$g, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint_typeiv$gF, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint_typeiv$gM, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint_typeiv$w, probs = c(0.05, 0.5, 0.95))

