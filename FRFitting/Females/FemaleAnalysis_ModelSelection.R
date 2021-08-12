################################################################################
### Female functional response analysis -- Model selection
################################################################################

### load required packages

library(dplyr); library(ggplot2); library(cowplot); library(rstan); library(loo)

### set some options for Stan

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

### load data

spider <- read.csv('FeedingData_Manipulated_f.csv')

### drop surveys with no midges

spider <- spider %>% filter(NumberMidges > 0)

### note that sometimes one or more of the chains in the temperature models 
### will get stuck in a weird place and cause divergent transitions and high or NA Rhat values. 
### Running the model again (or a couple times) usually fixes the problem

################################################################################
### fit full model -- quadratic temperature, separate interference
################################################################################

quadtemp_sepint_data <- list(N = nrow(spider),
                             y = spider$NumberFeed,
                             n_t = spider$NumberPred,
                             d = spider$DetectionEst,
                             R = spider$NumberMidges,
                             C = spider$FemaleDens,
                             M = spider$MaleDens,
                             J = spider$JuvDens,
                             Temp = spider$TempStart,
                             Temp2 = spider$TempStart^2)

quadtemp_sepint_fit <- stan(file = 'Female_tempquad_BDall.stan', data = quadtemp_sepint_data)

samples_quadtemp_sepint <- rstan::extract(quadtemp_sepint_fit)

### get summary of model fit

summary_quadtemp_sepint <- summary(quadtemp_sepint_fit)

summary_quadtemp_sepint$summary[1:6,]

log_lik_quadtemp_sepint <- extract_log_lik(quadtemp_sepint_fit, merge_chains = FALSE)

rel_eff_quadtemp_sepint <- relative_eff(log_lik_quadtemp_sepint)

loo_quadtemp_sepint <- loo(log_lik_quadtemp_sepint, r_eff = rel_eff_quadtemp_sepint)

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

quadtemp_combint_fit <- stan(file = 'Female_tempquad_BDcomb.stan', data = quadtemp_combint_data)

samples_quadtemp_combint <- rstan::extract(quadtemp_combint_fit)

### get summary of model fit

summary_quadtemp_combint <- summary(quadtemp_combint_fit)

summary_quadtemp_combint$summary[1:4,]

log_lik_quadtemp_combint <- extract_log_lik(quadtemp_combint_fit, merge_chains = FALSE)

rel_eff_quadtemp_combint <- relative_eff(log_lik_quadtemp_combint)

loo_quadtemp_combint <- loo(log_lik_quadtemp_combint, r_eff = rel_eff_quadtemp_combint)

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

quadtemp_noint_fit <- stan(file = 'Female_tempquad_BDnone.stan', data = quadtemp_noint_data)

samples_quadtemp_noint <- rstan::extract(quadtemp_noint_fit)

### get summary of model fit

summary_quadtemp_noint <- summary(quadtemp_noint_fit)

summary_quadtemp_noint$summary[1:4,]

log_lik_quadtemp_noint <- extract_log_lik(quadtemp_noint_fit, merge_chains = FALSE)

rel_eff_quadtemp_noint <- relative_eff(log_lik_quadtemp_noint)

loo_quadtemp_noint <- loo(log_lik_quadtemp_noint, r_eff = rel_eff_quadtemp_noint)

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
                            C = spider$FemaleDens,
                            M = spider$MaleDens,
                            J = spider$JuvDens)

notemp_sepint_fit <- stan(file = 'Female_notemp_BDall.stan', data = notemp_sepint_data)

samples_notemp_sepint <- rstan::extract(notemp_sepint_fit)

### get summary of model fit

summary_notemp_sepint <- summary(notemp_sepint_fit)

summary_notemp_sepint$summary[1:5,]

log_lik_notemp_sepint <- extract_log_lik(notemp_sepint_fit, merge_chains = FALSE)

rel_eff_notemp_sepint <- relative_eff(log_lik_notemp_sepint)

loo_notemp_sepint <- loo(log_lik_notemp_sepint, r_eff = rel_eff_notemp_sepint)

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

notemp_combint_fit <- stan(file = 'Female_notemp_BDcomb.stan', data = notemp_combint_data)

samples_notemp_combint <- rstan::extract(notemp_combint_fit)

### get summary of model fit

summary_notemp_combint <- summary(notemp_combint_fit)

summary_notemp_combint$summary[1:5,]

log_lik_notemp_combint <- extract_log_lik(notemp_combint_fit, merge_chains = FALSE)

rel_eff_notemp_combint <- relative_eff(log_lik_notemp_combint)

loo_notemp_combint <- loo(log_lik_notemp_combint, r_eff = rel_eff_notemp_combint)

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

notemp_noint_fit <- stan(file = 'Female_notemp_BDnone.stan', data = notemp_noint_data)

samples_notemp_noint <- rstan::extract(notemp_noint_fit)

### get summary of model fit

summary_notemp_noint <- summary(notemp_noint_fit)

summary_notemp_noint$summary[1:5,]

log_lik_notemp_noint <- extract_log_lik(notemp_noint_fit, merge_chains = FALSE)

rel_eff_notemp_noint <- relative_eff(log_lik_notemp_noint)

loo_notemp_noint <- loo(log_lik_notemp_noint, r_eff = rel_eff_notemp_noint)


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

null_fit <- stan(file = 'Female_null.stan', data = null_data)

samples_null <- rstan::extract(null_fit)

### get summary of model fit

summary_null <- summary(null_fit)

summary_null$summary[1:5,]

log_lik_null <- extract_log_lik(null_fit, merge_chains = FALSE)

rel_eff_null <- relative_eff(log_lik_null)

loo_null <- loo(log_lik_null, r_eff = rel_eff_null)


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

### save each of the objects containing the MCMC samples

# name each of the objects containing MCMC samples to be female specific

samples_quadtemp_sepint_f <- samples_quadtemp_sepint

samples_quadtemp_combint_f <- samples_quadtemp_combint

samples_quadtemp_noint_f <- samples_quadtemp_noint

samples_notemp_sepint_f <- samples_notemp_sepint

samples_notemp_combint_f <- samples_notemp_combint

samples_notemp_noint_f <- samples_notemp_noint

# save WAIC weights for model averaging

waic_wts_f <- round(waic_wts, 5)

# save all of these to a .RData file

save(samples_quadtemp_sepint_f, samples_quadtemp_combint_f, samples_quadtemp_noint_f,
     samples_notemp_sepint_f, samples_notemp_combint_f, samples_notemp_noint_f, 
     waic_wts_f, file = 'FemaleModelAvgData.RData')





