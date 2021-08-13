################################################################################
### Code to do model averaging using the posterior distributions
### Code for interference effects -- perhaps just females here?
################################################################################

### load libraries

library(dplyr); library(ggplot2); library(cowplot); library(viridis); library(RColorBrewer)

### load the RData files with the MCMC samples for all of the sexes/stages

load('FemaleModelAvgData.RData')

### also load spider survey data 

females <- read.csv('FeedingData_Manipulated_f.csv')

### model averaged functional responses -- females

### data required to make predictions

AvgPrey <- mean(females$NumberMidges)

AvgTemp <- mean(females$TempStart)

AvgDetection <- mean(females$DetectionEst)

AvgFemale <- 0

AvgMale <- 0

AvgJuvenile <- 0

FemaleRange <- seq(from = 0, to = 0.83, length.out = 100)

MaleRange <- seq(from = 0, to = 0.43, length.out = 100)

JuvenileRange <- seq(from = 0, to = 0.36, length.out = 100)

################################################################################
### Predictions from model with temperature, separate interference
################################################################################

samples_quadtemp_sepint_f <-  as.data.frame(samples_quadtemp_sepint_f[1:6])

### use apply statement to get feeding rate predictions

### feeding rates with variation in female densities

FR_pred_temp_sepint_f <- apply(samples_quadtemp_sepint_f, 1, function(x) 
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgPrey/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*AvgPrey + x['g']*FemaleRange + x['gM']*AvgMale + x['gJ']*AvgJuvenile))

FR_pred_temp_sepint_m <- apply(samples_quadtemp_sepint_f, 1, function(x) 
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgPrey/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*AvgPrey + x['g']*AvgFemale + x['gM']*MaleRange + x['gJ']*AvgJuvenile))

FR_pred_temp_sepint_j <- apply(samples_quadtemp_sepint_f, 1, function(x) 
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgPrey/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*AvgPrey + x['g']*AvgFemale + x['gM']*AvgMale + x['gJ']*JuvenileRange))

################################################################################
### Predictions from model with temperature and combined interference terms
################################################################################

samples_quadtemp_combint_f <- as.data.frame(samples_quadtemp_combint_f[1:4])

### use apply statement to get feeding rate predictions

FR_pred_temp_combint_f <- apply(samples_quadtemp_combint_f, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgPrey/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*AvgPrey + x['g']*FemaleRange + x['g']*AvgMale + x['g']*AvgJuvenile))

FR_pred_temp_combint_m <- apply(samples_quadtemp_combint_f, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgPrey/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*AvgPrey + x['g']*AvgFemale + x['g']*MaleRange + x['g']*AvgJuvenile))

FR_pred_temp_combint_j <- apply(samples_quadtemp_combint_f, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgPrey/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*AvgPrey + x['g']*AvgFemale + x['g']*AvgMale + x['g']*JuvenileRange))

################################################################################
### Predictions from model with temperature and no interference
################################################################################

samples_quadtemp_noint_f <-  as.data.frame(samples_quadtemp_noint_f[1:3])

### use apply statement to get feeding rate predictions

FR_pred_temp_noint_f <- apply(samples_quadtemp_noint_f, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgPrey/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*AvgPrey))

FR_pred_temp_noint_f <- matrix(rep(FR_pred_temp_noint_f, 100), ncol = 4000, byrow = TRUE)

################################################################################
### Predictions from model with no temperature and separate interference
################################################################################

samples_notemp_sepint_f <- as.data.frame(samples_notemp_sepint_f[1:4])

### use apply statement to get feeding rate predictions

FR_pred_notemp_sepint_f <- apply(samples_notemp_sepint_f, 1, function(x)
  x['a']*AvgPrey/(1 + x['a']*AvgPrey*AvgDetection + x['g']*FemaleRange + x['gM']*AvgMale + x['gJ']*AvgJuvenile))

FR_pred_notemp_sepint_m <- apply(samples_notemp_sepint_f, 1, function(x)
  x['a']*AvgPrey/(1 + x['a']*AvgPrey*AvgDetection + x['g']*AvgFemale + x['gM']*MaleRange + x['gJ']*AvgJuvenile))

FR_pred_notemp_sepint_j <- apply(samples_notemp_sepint_f, 1, function(x)
  x['a']*AvgPrey/(1 + x['a']*AvgPrey*AvgDetection + x['g']*AvgFemale + x['gM']*AvgMale + x['gJ']*JuvenileRange))

################################################################################
### Predictions from model with no temperature and combined interference
################################################################################

samples_notemp_combint_f <- as.data.frame(samples_notemp_combint_f[1:2])

### use apply statement to get feeding rate predictions

FR_pred_notemp_combint_f <- apply(samples_notemp_combint_f, 1, function(x)
  x['a']*AvgPrey/(1 + x['a']*AvgPrey*AvgDetection + x['g']*FemaleRange + x['g']*AvgMale + x['g']*AvgJuvenile))

FR_pred_notemp_combint_m <- apply(samples_notemp_combint_f, 1, function(x)
  x['a']*AvgPrey/(1 + x['a']*AvgPrey*AvgDetection + x['g']*AvgFemale + x['g']*MaleRange + x['g']*AvgJuvenile))

FR_pred_notemp_combint_j <- apply(samples_notemp_combint_f, 1, function(x)
  x['a']*AvgPrey/(1 + x['a']*AvgPrey*AvgDetection + x['g']*AvgFemale + x['g']*AvgMale + x['g']*JuvenileRange))


################################################################################
### Predictions from model with no temperature and no interference
################################################################################

samples_notemp_noint_f <- as.data.frame(samples_notemp_noint_f[1])

### use apply statement to get feeding rate predictions

FR_pred_notemp_noint_f <- apply(samples_notemp_noint_f, 1, function(x)
  x['a']*AvgPrey/(1 + x['a']*AvgPrey*AvgDetection))

FR_pred_notemp_noint_f <- matrix(rep(FR_pred_notemp_noint_f, 100), ncol = 4000, byrow = TRUE)

################################################################################
### need to put the predictions together using Bayesian model averaging
################################################################################

round(waic_wts_f * 4000)

################################################################################
### Interference from females
################################################################################

### samples from + temp sep interference

FR_pred_temp_sepint_f <- apply(FR_pred_temp_sepint_f, 1, function(x) sample(x, 82))

### samples from + temp comb interference

FR_pred_temp_combint_f <- apply(FR_pred_temp_combint_f, 1, function(x) sample(x, 161))

### samples from + temp no interference

FR_pred_temp_noint_f <- apply(FR_pred_temp_noint_f, 1, function(x) sample(x, 60))

### samples from no temp sep interference

FR_pred_notemp_sepint_f <- apply(FR_pred_notemp_sepint_f, 1, function(x) sample(x, 1707))

### samples from no temp comb interference

FR_pred_notemp_combint_f <- apply(FR_pred_notemp_combint_f, 1, function(x) sample(x, 1490))

### samples from no temp no interference

FR_pred_notemp_noint_f <- apply(FR_pred_notemp_noint_f, 1, function(x) sample(x, 501))

### now we can put all of these together into a single matrix

FR_int_f <- rbind(FR_pred_temp_sepint_f, FR_pred_temp_combint_f, FR_pred_temp_noint_f,
                   FR_pred_notemp_sepint_f, FR_pred_notemp_combint_f, FR_pred_notemp_noint_f)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

FR_FemaleInterference_DataFrame <- data.frame(FemaleDensity = FemaleRange,
                                Q5 = apply(FR_int_f, 2, function(x) quantile(x, probs = c(0.05))),
                                Q25 = apply(FR_int_f, 2, function(x) quantile(x, probs = c(0.25))),
                                Q50 = apply(FR_int_f, 2, function(x) quantile(x, probs = c(0.50))),
                                Q75 = apply(FR_int_f, 2, function(x) quantile(x, probs = c(0.75))),
                                Q95 = apply(FR_int_f, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = FR_FemaleInterference_DataFrame, aes(x = FemaleRange, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

################################################################################
### Now have to put together the interference predictions for males
################################################################################

### samples from + temp sep interference

FR_pred_temp_sepint_m <- apply(FR_pred_temp_sepint_m, 1, function(x) sample(x, 82))

### samples from + temp comb interference

FR_pred_temp_combint_m <- apply(FR_pred_temp_combint_m, 1, function(x) sample(x, 161))

### samples from + temp no interference

# FR_pred_temp_noint_f <- apply(FR_pred_temp_noint_f, 1, function(x) sample(x, 60))

### samples from no temp sep interference

FR_pred_notemp_sepint_m <- apply(FR_pred_notemp_sepint_m, 1, function(x) sample(x, 1707))

### samples from no temp comb interference

FR_pred_notemp_combint_m <- apply(FR_pred_notemp_combint_m, 1, function(x) sample(x, 1490))

### samples from no temp no interference

# FR_pred_notemp_noint_f <- apply(FR_pred_notemp_noint_f, 1, function(x) sample(x, 501))

### now we can put all of these together into a single matrix

FR_int_m <- rbind(FR_pred_temp_sepint_m, FR_pred_temp_combint_m, FR_pred_temp_noint_f,
                  FR_pred_notemp_sepint_m, FR_pred_notemp_combint_m, FR_pred_notemp_noint_f)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

FR_MaleInterference_DataFrame <- data.frame(MaleDensity = MaleRange,
                                              Q5 = apply(FR_int_m, 2, function(x) quantile(x, probs = c(0.05))),
                                              Q25 = apply(FR_int_m, 2, function(x) quantile(x, probs = c(0.25))),
                                              Q50 = apply(FR_int_m, 2, function(x) quantile(x, probs = c(0.50))),
                                              Q75 = apply(FR_int_m, 2, function(x) quantile(x, probs = c(0.75))),
                                              Q95 = apply(FR_int_m, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = FR_MaleInterference_DataFrame, aes(x = MaleDensity, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

################################################################################
### Now, lastly need to do the same for juveniles
################################################################################

### samples from + temp sep interference

FR_pred_temp_sepint_j <- apply(FR_pred_temp_sepint_j, 1, function(x) sample(x, 82))

### samples from + temp comb interference

FR_pred_temp_combint_j <- apply(FR_pred_temp_combint_j, 1, function(x) sample(x, 161))

### samples from + temp no interference

# FR_pred_temp_noint_f <- apply(FR_pred_temp_noint_f, 1, function(x) sample(x, 60))

### samples from no temp sep interference

FR_pred_notemp_sepint_j <- apply(FR_pred_notemp_sepint_j, 1, function(x) sample(x, 1707))

### samples from no temp comb interference

FR_pred_notemp_combint_j <- apply(FR_pred_notemp_combint_j, 1, function(x) sample(x, 1490))

### samples from no temp no interference

# FR_pred_notemp_noint_f <- apply(FR_pred_notemp_noint_f, 1, function(x) sample(x, 501))

### now we can put all of these together into a single matrix

FR_int_j <- rbind(FR_pred_temp_sepint_j, FR_pred_temp_combint_j, FR_pred_temp_noint_f,
                  FR_pred_notemp_sepint_j, FR_pred_notemp_combint_j, FR_pred_notemp_noint_f)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

FR_JuvenileInterference_DataFrame <- data.frame(JuvenileDensity = JuvenileRange,
                                            Q5 = apply(FR_int_j, 2, function(x) quantile(x, probs = c(0.05))),
                                            Q25 = apply(FR_int_j, 2, function(x) quantile(x, probs = c(0.25))),
                                            Q50 = apply(FR_int_j, 2, function(x) quantile(x, probs = c(0.50))),
                                            Q75 = apply(FR_int_j, 2, function(x) quantile(x, probs = c(0.75))),
                                            Q95 = apply(FR_int_j, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = FR_JuvenileInterference_DataFrame, aes(x = JuvenileDensity, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

### put all of the data frames together into a single data frame

Female_Int_DataFrame <- bind_rows(FR_FemaleInterference_DataFrame, FR_MaleInterference_DataFrame, 
                                  FR_JuvenileInterference_DataFrame)

Female_Int_DataFrame <- Female_Int_DataFrame %>% mutate(Density = coalesce(FemaleDensity, MaleDensity, JuvenileDensity))

Female_Int_DataFrame$SexStage <- rep(c('Female', 'Male', 'Juvenile'), each = 100) 

# ggplot(data = Female_Int_DataFrame, aes(x = Density, y = Q50/24, color = SexStage)) + geom_line(size = 1) + 
#   geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24, fill = SexStage), alpha = 0.25, color = NA) + 
#   geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24, fill = SexStage), alpha = 0.5, color = NA) + 
#   theme_cowplot() + xlab('Density of Sex/Stage') + ylab('Feeding Rate \n(Midges per hour)') + labs(fill = 'Sex/Stage') +
#   labs(color = 'Sex/Stage') + scale_color_brewer(palette = 'Dark2') + scale_fill_brewer(palette = 'Dark2') + ggtitle('Females')
# 
# save(Female_Int_DataFrame, file = 'Female_Int_Plot.RData')
# 













