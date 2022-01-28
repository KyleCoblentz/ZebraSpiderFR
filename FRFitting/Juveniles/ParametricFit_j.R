################################################################################
### Parametric fit of juvenile functional response to data
################################################################################

### load packages

library(dplyr); library(ggplot2); library(cowplot); library(rstan); library(loo)

### set some options for Stan

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

### load data

juveniles <- read.csv('FeedingData_Manipulated_j.csv', stringsAsFactors = TRUE)

### manipulate data

juveniles <- juveniles %>% mutate(FeedingRate = ProportionFeed/(DetectionEst*24),
                                  TotalSpiderDens = FemaleDens + JuvenileDens + MaleDens)


################################################################################
### fit full model -- quadratic temperature, separate interference
################################################################################

quadtemp_sepint_data <- list(N = nrow(juveniles),
                             y = juveniles$NumberFeed,
                             n_t = juveniles$NumberPred,
                             d = juveniles$DetectionEst*24,
                             R = juveniles$NumberMidges,
                             C = juveniles$JuvenileDens,
                             F = juveniles$FemaleDens,
                             M = juveniles$MaleDens,
                             Temp = juveniles$TempStart,
                             Temp2 = juveniles$TempStart^2)

quadtemp_sepint_fit <- stan(file = 'Juvenile_tempquad_BDall.stan', data = quadtemp_sepint_data)

samples_quadtemp_sepint <- rstan::extract(quadtemp_sepint_fit)

### get summary of model fit

summary_quadtemp_sepint <- summary(quadtemp_sepint_fit)

summary_quadtemp_sepint$summary[1:6,]

quantile(samples_quadtemp_sepint$ca, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint$ba, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint$qa, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint$g, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint$gF, probs = c(0.05, 0.5, 0.95))

quantile(samples_quadtemp_sepint$gM, probs = c(0.05, 0.5, 0.95))


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

