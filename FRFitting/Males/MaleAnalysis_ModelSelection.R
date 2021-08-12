################################################################################
### Male functional response analysis -- Model selection
################################################################################

### load required packages

library(dplyr); library(ggplot2); library(cowplot); library(rstan); library(loo)

### set some options for Stan

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

### load data

spider <- read.csv('FeedingData_Manipulated_m.csv')

### note that sometimes a chain gets stuck in a weird place fitting the temperature dependent models
### this will lead to divergent transitions and Rhat values that are high or NA. Running the 
### model again (or a few times) usually works.

################################################################################
### fit full model -- quadratic temperature, separate interference
################################################################################

quadtemp_sepint_data <- list(N = nrow(spider),
                             y = spider$NumberFeed,
                             n_t = spider$NumberPred,
                             d = spider$DetectionEst,
                             R = spider$NumberMidges,
                             C = spider$MaleDens,
                             F = spider$FemaleDens,
                             J = spider$JuvDens,
                             Temp = spider$TempStart,
                             Temp2 = spider$TempStart^2)

quadtemp_sepint_fit <- stan(file = 'Male_tempquad_BDall.stan', data = quadtemp_sepint_data)

samples_quadtemp_sepint <- rstan::extract(quadtemp_sepint_fit)

### get summary of model fit

summary_quadtemp_sepint <- summary(quadtemp_sepint_fit)

summary_quadtemp_sepint$summary[1:6,]

log_lik_quadtemp_sepint <- extract_log_lik(quadtemp_sepint_fit, merge_chains = FALSE)

waic_quadtemp_sepint <-waic(log_lik_quadtemp_sepint)

################################################################################
### same model but with interference combined to a single interference parameter
################################################################################

quadtemp_combint_data <- list(N = nrow(spider),
                             y = spider$NumberFeed,
                             n_t = spider$NumberPred,
                             d = spider$DetectionEst,
                             R = spider$NumberMidges,
                             C = spider$FemaleDens + spider$MaleDens + spider$JuvDens,
                             Temp = spider$TempStart,
                             Temp2 = spider$TempStart^2)

quadtemp_combint_fit <- stan(file = 'Male_tempquad_BDcomb.stan', data = quadtemp_combint_data)

samples_quadtemp_combint <- rstan::extract(quadtemp_combint_fit)

### get summary of model fit

summary_quadtemp_combint <- summary(quadtemp_combint_fit)

summary_quadtemp_combint$summary[1:4,]

log_lik_quadtemp_combint <- extract_log_lik(quadtemp_combint_fit, merge_chains = FALSE)

waic_quadtemp_combint <- waic(log_lik_quadtemp_combint)

loo_compare(waic_quadtemp_sepint, waic_quadtemp_combint)

################################################################################
### quadratic temperature -- no interference
################################################################################

quadtemp_noint_data <- list(N = nrow(spider),
                              y = spider$NumberFeed,
                              n_t = spider$NumberPred,
                              d = spider$DetectionEst,
                              R = spider$NumberMidges,
                              Temp = spider$TempStart,
                              Temp2 = spider$TempStart^2)

quadtemp_noint_fit <- stan(file = 'Male_tempquad_BDnone.stan', data = quadtemp_noint_data)

samples_quadtemp_noint <- rstan::extract(quadtemp_noint_fit)

### get summary of model fit

summary_quadtemp_noint <- summary(quadtemp_noint_fit)

summary_quadtemp_noint$summary[1:4,]

log_lik_quadtemp_noint <- extract_log_lik(quadtemp_noint_fit, merge_chains = FALSE)

waic_quadtemp_noint <-waic(log_lik_quadtemp_noint)

loo_compare(waic_quadtemp_sepint, waic_quadtemp_combint, waic_quadtemp_noint)


################################################################################
### no temperature -- interference separate
################################################################################

notemp_sepint_data <- list(N = nrow(spider),
                            y = spider$NumberFeed,
                            n_t = spider$NumberPred,
                            d = spider$DetectionEst,
                            R = spider$NumberMidges,
                            C = spider$MaleDens,
                            F = spider$FemaleDens,
                            J = spider$JuvDens)

notemp_sepint_fit <- stan(file = 'Male_notemp_BDall.stan', data = notemp_sepint_data)

samples_notemp_sepint <- rstan::extract(notemp_sepint_fit)

### get summary of model fit

summary_notemp_sepint <- summary(notemp_sepint_fit)

summary_notemp_sepint$summary[1:5,]

log_lik_notemp_sepint <- extract_log_lik(notemp_sepint_fit, merge_chains = FALSE)

waic_notemp_sepint <-waic(log_lik_notemp_sepint)

loo_compare(waic_quadtemp_sepint, waic_quadtemp_combint, waic_quadtemp_noint, waic_notemp_sepint)

################################################################################
### no temperature -- interference combined
################################################################################

notemp_combint_data <- list(N = nrow(spider),
                           y = spider$NumberFeed,
                           n_t = spider$NumberPred,
                           d = spider$DetectionEst,
                           R = spider$NumberMidges,
                           C = spider$FemaleDens + spider$MaleDens + spider$JuvDens)

notemp_combint_fit <- stan(file = 'Male_notemp_BDcomb.stan', data = notemp_combint_data)

samples_notemp_combint <- rstan::extract(notemp_combint_fit)

### get summary of model fit

summary_notemp_combint <- summary(notemp_combint_fit)

summary_notemp_combint$summary[1:5,]

log_lik_notemp_combint <- extract_log_lik(notemp_combint_fit, merge_chains = FALSE)

waic_notemp_combint <-waic(log_lik_notemp_combint)

loo_compare(waic_quadtemp_sepint, waic_quadtemp_combint, waic_quadtemp_noint, waic_notemp_sepint,
            waic_notemp_combint)

################################################################################
### Type II functional response -- no temp, no interference
################################################################################

notemp_noint_data <- list(N = nrow(spider),
                            y = spider$NumberFeed,
                            n_t = spider$NumberPred,
                            d = spider$DetectionEst,
                            R = spider$NumberMidges)

notemp_noint_fit <- stan(file = 'Male_notemp_BDnone.stan', data = notemp_noint_data)

samples_notemp_noint <- rstan::extract(notemp_noint_fit)

### get summary of model fit

summary_notemp_noint <- summary(notemp_noint_fit)

summary_notemp_noint$summary[1:5,]

log_lik_notemp_noint <- extract_log_lik(notemp_noint_fit, merge_chains = FALSE)

waic_notemp_noint <- waic(log_lik_notemp_noint)

loo_compare(waic_quadtemp_sepint, waic_quadtemp_combint, waic_quadtemp_noint, waic_notemp_sepint,
            waic_notemp_combint, waic_notemp_noint)

################################################################################
### null model
################################################################################

null_data <- list(N = nrow(spider),
                  y = spider$NumberFeed,
                  n_t = spider$NumberPred,
                  d = spider$DetectionEst)

null_fit <- stan(file = 'Male_null.stan', data = null_data)

samples_null <- rstan::extract(null_fit)


### get summary of model fit

summary_null <- summary(null_fit)

summary_null$summary[1:5,]

log_lik_null <- extract_log_lik(null_fit, merge_chains = FALSE)

waic_null <-waic(log_lik_null)

print(loo_compare(waic_quadtemp_sepint, waic_quadtemp_combint, waic_quadtemp_noint, waic_notemp_sepint,
            waic_notemp_combint, waic_notemp_noint, waic_null), simplify = FALSE)

### get weights for the models

waics <- c(
  waic_quadtemp_sepint$estimates["elpd_waic", 1],
  waic_quadtemp_combint$estimates["elpd_waic", 1],
  waic_quadtemp_noint$estimates["elpd_waic", 1],
  waic_notemp_sepint$estimates["elpd_waic", 1],
  waic_notemp_combint$estimates["elpd_waic", 1],
  waic_notemp_noint$estimates["elpd_waic", 1],
  waic_null$estimates["elpd_waic", 1]
)

waic_wts <- exp(waics) / sum(exp(waics))
round(waic_wts, 5)

### save objects containing MCMC samples for model averaging

### make MCMC sample objects male specific

samples_quadtemp_sepint_m <- samples_quadtemp_sepint

samples_quadtemp_combint_m <- samples_quadtemp_combint

samples_quadtemp_noint_m <- samples_quadtemp_noint

samples_notemp_sepint_m <- samples_notemp_sepint

samples_notemp_combint_m <- samples_notemp_combint

samples_notemp_noint_m <- samples_notemp_noint

samples_null_m <- samples_null

### save waic weights as male specific

waic_wts_m <- waic_wts

### save all of the info for model averaging 

save(samples_quadtemp_sepint_m, samples_quadtemp_combint_m, samples_quadtemp_noint_m,
     samples_notemp_sepint_m, samples_notemp_combint_m, samples_notemp_noint_m,
     samples_null_m, waic_wts_m, file = 'MaleModelAvgData.RData')


















