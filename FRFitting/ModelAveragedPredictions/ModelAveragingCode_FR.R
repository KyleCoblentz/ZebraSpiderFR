################################################################################
### Code to do model averaging using the posterior distributions
### Code for the functional response comparsion
################################################################################

### load libraries

library(dplyr); library(ggplot2); library(cowplot); library(RColorBrewer);

### load the RData files with the MCMC samples for all of the sexes/stages

load('FemaleModelAvgData.RData')

load('JuvenileModelAvgData.RData')

load('MaleModelAvgData.RData')

### also load spider survey data 

females <- read.csv('FeedingData_Manipulated_f.csv')

males <- read.csv('FeedingData_Manipulated_m.csv')

juveniles <- read.csv('FeedingData_Manipulated_j.csv')

### model averaged functional responses -- females

### data required to make predictions

PreyRange <- seq(from = 0, to = 7, length.out = 100)

AvgTemp <- mean(females$TempStart)

AvgFemale <- mean(females$FemaleDens)

AvgMale <- mean(females$MaleDens)

AvgJuvenile <- mean(females$JuvDens)

AvgDetection <- mean(females$DetectionEst)

################################################################################
### Predictions from model with temperature and sep interference terms
################################################################################

samples_quadtemp_sepint_f <-  as.data.frame(samples_quadtemp_sepint_f[1:6])

### use apply statement to get feeding rate predictions

FR_pred_temp_sepint_f <- apply(samples_quadtemp_sepint_f, 1, function(x) 
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange + x['g']*AvgFemale + x['gM']*AvgMale + x['gJ']*AvgJuvenile))

################################################################################
### Predictions from model with temperature and combined interference terms
################################################################################

samples_quadtemp_combint_f <- as.data.frame(samples_quadtemp_combint_f[1:4])

### use apply statement to get feeding rate predictions

FR_pred_temp_combint_f <- apply(samples_quadtemp_combint_f, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange + x['g']*AvgFemale + x['g']*AvgMale + x['g']*AvgJuvenile))

################################################################################
### Predictions from model with temperature and no interference
################################################################################

samples_quadtemp_noint_f <-  as.data.frame(samples_quadtemp_noint_f[1:3])

### use apply statement to get feeding rate predictions

FR_pred_temp_noint_f <- apply(samples_quadtemp_noint_f, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange))

################################################################################
### Predictions from model with no temperature and separate interference
################################################################################

samples_notemp_sepint_f <- as.data.frame(samples_notemp_sepint_f[1:4])

### use apply statement to get feeding rate predictions

FR_pred_notemp_sepint_f <- apply(samples_notemp_sepint_f, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection + x['g']*AvgFemale + x['gM']*AvgMale + x['gJ']*AvgJuvenile))

################################################################################
### Predictions from model with no temperature and combined interference
################################################################################

samples_notemp_combint_f <- as.data.frame(samples_notemp_combint_f[1:2])

### use apply statement to get feeding rate predictions

FR_pred_notemp_combint_f <- apply(samples_notemp_combint_f, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection + x['g']*AvgFemale + x['g']*AvgMale + x['g']*AvgJuvenile))

################################################################################
### Predictions from model with no temperature and no interference
################################################################################

samples_notemp_noint_f <- as.data.frame(samples_notemp_noint_f[1])

### use apply statement to get feeding rate predictions

FR_pred_notemp_noint_f <- apply(samples_notemp_noint_f, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection))

################################################################################
### now need to put all of the predictions together using Bayesian model
### averaging
################################################################################

### need to make a new matrix that has samples at each of the prey densities
### in proportion to the WAIC weights

### number of samples that should be included from each of the models

round(waic_wts_f * 4000)

### sample from temp + sep interference

FR_pred_temp_sepint_f <- apply(FR_pred_temp_sepint_f, 1, function(x) sample(x, 82))

### sample from temp + combined interference

FR_pred_temp_combint_f <- apply(FR_pred_temp_combint_f, 1, function(x) sample(x, 161))

### sample from temp + no interference

FR_pred_temp_noint_f <- apply(FR_pred_temp_noint_f, 1, function(x) sample(x, 60))

### sample from no temp + sep interference

FR_pred_notemp_sepint_f <- apply(FR_pred_notemp_sepint_f, 1, function(x) sample(x, 1707))

### sample from no temp + combined interference 

FR_pred_notemp_combint_f <- apply(FR_pred_notemp_combint_f, 1, function(x) sample(x, 1490))

### sample from no temp + no interference

FR_pred_notemp_noint_f <-  apply(FR_pred_notemp_noint_f, 1, function(x) sample(x, 501))

### now we can put all of these together into a single matrix

FR_pred_f <- rbind(FR_pred_temp_sepint_f, FR_pred_temp_combint_f, FR_pred_temp_noint_f,
                   FR_pred_notemp_sepint_f, FR_pred_notemp_combint_f, FR_pred_notemp_noint_f)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

FR_Female_DataFrame <- data.frame(MidgeDensity = seq(from = 0, to = 7, length.out = 100),
                                  Q5 = apply(FR_pred_f, 2, function(x) quantile(x, probs = c(0.05))),
                                  Q25 = apply(FR_pred_f, 2, function(x) quantile(x, probs = c(0.25))),
                                  Q50 = apply(FR_pred_f, 2, function(x) quantile(x, probs = c(0.50))),
                                  Q75 = apply(FR_pred_f, 2, function(x) quantile(x, probs = c(0.75))),
                                  Q95 = apply(FR_pred_f, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = FR_Female_DataFrame, aes(x = MidgeDensity, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

################################################################################
################################################################################
### Now need to do the same thing for males
################################################################################
################################################################################

### data required to make predictions

PreyRange <- seq(from = 0, to = 7, length.out = 100)

AvgTemp <- mean(males$TempStart)

AvgFemale <- mean(males$FemaleDens)

AvgMale <- mean(males$MaleDens)

AvgJuvenile <- mean(males$JuvDens)

AvgDetection <- mean(males$DetectionEst)

################################################################################
### Predictions from model with temperature and sep interference terms
################################################################################

samples_quadtemp_sepint_m <-  as.data.frame(samples_quadtemp_sepint_m[1:6])

### use apply statement to get feeding rate predictions

FR_pred_temp_sepint_m <- apply(samples_quadtemp_sepint_m, 1, function(x) 
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange + x['g']*AvgMale + x['gF']*AvgFemale + x['gJ']*AvgJuvenile))

################################################################################
### Predictions from model with temperature and combined interference terms
################################################################################

samples_quadtemp_combint_m <- as.data.frame(samples_quadtemp_combint_m[1:4])

### use apply statement to get feeding rate predictions

FR_pred_temp_combint_m <- apply(samples_quadtemp_combint_m, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange + x['g']*AvgFemale + x['g']*AvgMale + x['g']*AvgJuvenile))

################################################################################
### Predictions from model with temperature and no interference
################################################################################

samples_quadtemp_noint_m <-  as.data.frame(samples_quadtemp_noint_m[1:3])

### use apply statement to get feeding rate predictions

FR_pred_temp_noint_m <- apply(samples_quadtemp_noint_m, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange))

################################################################################
### Predictions from model with no temperature and separate interference
################################################################################

samples_notemp_sepint_m <- as.data.frame(samples_notemp_sepint_m[1:4])

### use apply statement to get feeding rate predictions

FR_pred_notemp_sepint_m <- apply(samples_notemp_sepint_m, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection + x['g']*AvgMale + x['gF']*AvgFemale + x['gJ']*AvgJuvenile))

################################################################################
### Predictions from model with no temperature and combined interference
################################################################################

samples_notemp_combint_m <- as.data.frame(samples_notemp_combint_m[1:2])

### use apply statement to get feeding rate predictions

FR_pred_notemp_combint_m <- apply(samples_notemp_combint_m, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection + x['g']*AvgFemale + x['g']*AvgMale + x['g']*AvgJuvenile))

################################################################################
### Predictions from model with no temperature and no interference
################################################################################

samples_notemp_noint_m <- as.data.frame(samples_notemp_noint_m[1])

### use apply statement to get feeding rate predictions

FR_pred_notemp_noint_m <- apply(samples_notemp_noint_m, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection))

################################################################################
### feeding rate predictions from the null model
################################################################################

samples_null_m <- as.data.frame(samples_null_m[1])

FR_pred_null_m <- matrix(rep(samples_null_m$f, 100), ncol = 4000, byrow = TRUE)

################################################################################
### now need to put all of the predictions together using Bayesian model
### averaging
################################################################################

### need to make a new matrix that has samples at each of the prey densities
### in proportion to the WAIC weights

### number of samples that should be included from each of the models

round(waic_wts_m * 4000)

### sample from temp + sep interference

FR_pred_temp_sepint_m <- apply(FR_pred_temp_sepint_m, 1, function(x) sample(x, 56))

### sample from temp + combined interference

FR_pred_temp_combint_m <- apply(FR_pred_temp_combint_m, 1, function(x) sample(x, 201))

### sample from temp + no interference

FR_pred_temp_noint_m <- apply(FR_pred_temp_noint_m, 1, function(x) sample(x, 734))

### sample from no temp + sep interference

FR_pred_notemp_sepint_m <- apply(FR_pred_notemp_sepint_m, 1, function(x) sample(x, 175))

### sample from no temp + combined interference 

FR_pred_notemp_combint_m <- apply(FR_pred_notemp_combint_m, 1, function(x) sample(x, 431))

### sample from no temp + no interference

FR_pred_notemp_noint_m <-  apply(FR_pred_notemp_noint_m, 1, function(x) sample(x, 2295))

### sample from null model

FR_pred_null_m <- apply(FR_pred_null_m, 1, function(x) sample(x, 108))

### now we can put all of these together into a single matrix

FR_pred_m <- rbind(FR_pred_temp_sepint_m, FR_pred_temp_combint_m, FR_pred_temp_noint_m,
                   FR_pred_notemp_sepint_m, FR_pred_notemp_combint_m, FR_pred_notemp_noint_m,
                   FR_pred_null_m)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

FR_Male_DataFrame <- data.frame(MidgeDensity = seq(from = 0, to = 7, length.out = 100),
                                  Q5 = apply(FR_pred_m, 2, function(x) quantile(x, probs = c(0.05))),
                                  Q25 = apply(FR_pred_m, 2, function(x) quantile(x, probs = c(0.25))),
                                  Q50 = apply(FR_pred_m, 2, function(x) quantile(x, probs = c(0.50))),
                                  Q75 = apply(FR_pred_m, 2, function(x) quantile(x, probs = c(0.75))),
                                  Q95 = apply(FR_pred_m, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = FR_Male_DataFrame, aes(x = MidgeDensity, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

################################################################################
################################################################################
### Now need to do the same for Juveniles
################################################################################
################################################################################

### data required to make predictions

PreyRange <- seq(from = 0, to = 7, length.out = 100)

AvgTemp <- mean(juveniles$TempStart)

AvgFemale <- mean(juveniles$FemaleDens)

AvgMale <- mean(juveniles$MaleDens)

AvgJuvenile <- mean(juveniles$JuvenileDens)

AvgDetection <- mean(juveniles$DetectionEst)

################################################################################
### Predictions from model with temperature and sep interference terms
################################################################################

samples_quadtemp_sepint_j <-  as.data.frame(samples_quadtemp_sepint_j[1:6])

### use apply statement to get feeding rate predictions

FR_pred_temp_sepint_j <- apply(samples_quadtemp_sepint_j, 1, function(x) 
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange + x['g']*AvgJuvenile + x['gF']*AvgFemale + x['gM']*AvgMale))

################################################################################
### Predictions from model with temperature and combined interference terms
################################################################################

samples_quadtemp_combint_j <- as.data.frame(samples_quadtemp_combint_j[1:4])

### use apply statement to get feeding rate predictions

FR_pred_temp_combint_j <- apply(samples_quadtemp_combint_j, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange + x['g']*AvgFemale + x['g']*AvgMale + x['g']*AvgJuvenile))

################################################################################
### Predictions from model with temperature and no interference
################################################################################

samples_quadtemp_noint_j <-  as.data.frame(samples_quadtemp_noint_j[1:3])

### use apply statement to get feeding rate predictions

FR_pred_temp_noint_j <- apply(samples_quadtemp_noint_j, 1, function(x)
  x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*PreyRange/(1 + x['ca']*exp(x['ba']*AvgTemp + x['qa']*AvgTemp^2)*AvgDetection*PreyRange))

################################################################################
### Predictions from model with no temperature and separate interference
################################################################################

samples_notemp_sepint_j <- as.data.frame(samples_notemp_sepint_j[1:4])

### use apply statement to get feeding rate predictions

FR_pred_notemp_sepint_j <- apply(samples_notemp_sepint_j, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection + x['g']*AvgJuvenile + x['gF']*AvgFemale + x['gM']*AvgMale))

################################################################################
### Predictions from model with no temperature and combined interference
################################################################################

samples_notemp_combint_j <- as.data.frame(samples_notemp_combint_j[1:2])

### use apply statement to get feeding rate predictions

FR_pred_notemp_combint_j <- apply(samples_notemp_combint_j, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection + x['g']*AvgFemale + x['g']*AvgMale + x['g']*AvgJuvenile))

################################################################################
### Predictions from model with no temperature and no interference
################################################################################

samples_notemp_noint_j <- as.data.frame(samples_notemp_noint_j[1])

### use apply statement to get feeding rate predictions

FR_pred_notemp_noint_j <- apply(samples_notemp_noint_j, 1, function(x)
  x['a']*PreyRange/(1 + x['a']*PreyRange*AvgDetection))

################################################################################
### now need to put all of the predictions together using Bayesian model
### averaging
################################################################################

### need to make a new matrix that has samples at each of the prey densities
### in proportion to the WAIC weights

### number of samples that should be included from each of the models

round(waic_wts_j * 4000)

### sample from temp + sep interference

FR_pred_temp_sepint_j <- apply(FR_pred_temp_sepint_j, 1, function(x) sample(x, 847))

### sample from temp + combined interference

FR_pred_temp_combint_j <- apply(FR_pred_temp_combint_j, 1, function(x) sample(x, 589))

### sample from temp + no interference

FR_pred_temp_noint_j <- apply(FR_pred_temp_noint_j, 1, function(x) sample(x, 323))

### sample from no temp + sep interference

FR_pred_notemp_sepint_j <- apply(FR_pred_notemp_sepint_j, 1, function(x) sample(x, 880))

### sample from no temp + combined interference 

FR_pred_notemp_combint_j <- apply(FR_pred_notemp_combint_j, 1, function(x) sample(x, 901))

### sample from no temp + no interference

FR_pred_notemp_noint_j <-  apply(FR_pred_notemp_noint_j, 1, function(x) sample(x, 460))

### now we can put all of these together into a single matrix

FR_pred_j <- rbind(FR_pred_temp_sepint_j, FR_pred_temp_combint_j, FR_pred_temp_noint_j,
                   FR_pred_notemp_sepint_j, FR_pred_notemp_combint_j, FR_pred_notemp_noint_j)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

FR_Juvenile_DataFrame <- data.frame(MidgeDensity = seq(from = 0, to = 7, length.out = 100),
                                Q5 = apply(FR_pred_j, 2, function(x) quantile(x, probs = c(0.05))),
                                Q25 = apply(FR_pred_j, 2, function(x) quantile(x, probs = c(0.25))),
                                Q50 = apply(FR_pred_j, 2, function(x) quantile(x, probs = c(0.50))),
                                Q75 = apply(FR_pred_j, 2, function(x) quantile(x, probs = c(0.75))),
                                Q95 = apply(FR_pred_j, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = FR_Juvenile_DataFrame, aes(x = MidgeDensity, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

################################################################################
### put all of the data frames together
################################################################################

FR_DataFrame <- bind_rows(FR_Female_DataFrame, FR_Male_DataFrame, FR_Juvenile_DataFrame)

FR_DataFrame$SexStage <- rep(c('Female', 'Male', 'Juvenile'), each = 100) 

# FR_plot <- ggplot(data = FR_DataFrame, aes(x = MidgeDensity, y = Q50/24, color = SexStage)) + geom_line() + 
#   geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24, fill = SexStage), alpha = 0.25, color = NA) + 
#   geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24, fill = SexStage), alpha = 0.5, color = NA) + 
#   theme_cowplot() + xlab('Density of Midges') + ylab('Feeding Rate \n(Midges per hour)') + labs(fill = 'Sex/Stage') +
#   labs(color = 'Sex/Stage') + scale_color_brewer(palette = 'Dark2') + scale_fill_brewer(palette = 'Dark2')
# 
# 
# save_plot(filename = 'FunctionalResponseFigure.png', plot = FR_plot)



































