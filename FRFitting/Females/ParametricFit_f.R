################################################################################
### Parametric functional response fit -- females
################################################################################

### load packages

library(dplyr); library(ggplot2); library(cowplot); library(rstan);

### set some options for Stan

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

### load data

females <- read.csv('FeedingData_Manipulated_f.csv', stringsAsFactors = TRUE)

### filter data 

females <- females %>% filter(NumberMidges > 0)

females <- females %>% mutate(TotalPredDens = FemaleDens + JuvDens + MaleDens)

females <- females %>% mutate(FeedingRate = ProportionFeed/(DetectionEst*24))

### fit parametric model to data

################################################################################
### no temperature -- interference combined
################################################################################

notemp_combint_data <- list(N = nrow(females),
                            y = females$NumberFeed,
                            n_t = females$NumberPred,
                            d = females$DetectionEst*24,
                            R = females$NumberMidges,
                            C = females$FemaleDens + females$MaleDens + females$JuvDens)

notemp_combint_fit <- stan(file = 'Female_notemp_BDcomb.stan', data = notemp_combint_data)

samples_notemp_combint <- rstan::extract(notemp_combint_fit)

### get summary of model fit

summary_notemp_combint <- summary(notemp_combint_fit)

summary_notemp_combint$summary[1:5,]

quantile(samples_notemp_combint$a, probs = c(0.05, 0.5, 0.95))

quantile(samples_notemp_combint$g, probs = c(0.05, 0.5, 0.95))



