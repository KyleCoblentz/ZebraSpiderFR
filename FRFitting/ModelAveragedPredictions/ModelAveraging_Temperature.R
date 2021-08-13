################################################################################
### model averaging -- temperature effects on space clearance rate
################################################################################

library(dplyr); library(ggplot2); library(cowplot); 

### load the RData files with the MCMC samples for all of the sexes/stages

load('FemaleModelAvgData.RData')

load('JuvenileModelAvgData.RData')

load('MaleModelAvgData.RData')

### also load spider survey data 

females <- read.csv('FeedingData_Manipulated_f.csv')

males <- read.csv('FeedingData_Manipulated_m.csv')

juveniles <- read.csv('FeedingData_Manipulated_j.csv')

### model averaged temperature responses -- females

### should just need the temperature range to get the effects of temperature 
### on the space clearance rates

TempRange <- seq(from = 15.5, to = 45, length.out = 100)

### a predictions, temp and sep interference

samples_quadtemp_sepint_f <- as.data.frame(samples_quadtemp_sepint_f[1:6])

a_pred_temp_sepint_f <- apply(samples_quadtemp_sepint_f, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, temp and combined interference

samples_quadtemp_combint_f <- as.data.frame(samples_quadtemp_combint_f[1:4])

a_pred_temp_combint_f <- apply(samples_quadtemp_combint_f, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, temp and no interference

samples_quadtemp_noint_f <-  as.data.frame(samples_quadtemp_noint_f[1:3])

a_pred_temp_noint_f <- apply(samples_quadtemp_noint_f, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, no temp sep interference

samples_notemp_sepint_f <- as.data.frame(samples_notemp_sepint_f[1:4])

a_pred_notemp_sepint_f <- matrix(rep(samples_notemp_sepint_f$a, 100), ncol = 4000, byrow = TRUE)

### a predictions, no temp combined interference

samples_notemp_combint_f <- as.data.frame(samples_notemp_combint_f[1:2])

a_pred_notemp_combint_f <- matrix(rep(samples_notemp_combint_f$a, 100), ncol = 4000, byrow = TRUE)

### a predictions, no temp no interference

samples_notemp_noint_f <- as.data.frame(samples_notemp_noint_f[1])

a_pred_notemp_noint_f <- matrix(rep(samples_notemp_noint_f$a, 100), ncol = 4000, byrow = TRUE)

### sample predictions in accordance to waic weights

### number of samples that should be included from each of the models

round(waic_wts_f * 4000)

### sample from temp + sep interference

a_pred_temp_sepint_f <- apply(a_pred_temp_sepint_f, 1, function(x) sample(x, 82))

### sample from temp + combined interference

a_pred_temp_combint_f <- apply(a_pred_temp_combint_f, 1, function(x) sample(x, 161))

### sample from temp + no interference

a_pred_temp_noint_f <- apply(a_pred_temp_noint_f, 1, function(x) sample(x, 60))

### sample from no temp + sep interference

a_pred_notemp_sepint_f <- apply(a_pred_notemp_sepint_f, 1, function(x) sample(x, 1707))

### sample from no temp + combined interference 

a_pred_notemp_combint_f <- apply(a_pred_notemp_combint_f, 1, function(x) sample(x, 1490))

### sample from no temp + no interference

a_pred_notemp_noint_f <-  apply(a_pred_notemp_noint_f, 1, function(x) sample(x, 501))

### now we can put all of these together into a single matrix

a_pred_f <- rbind(a_pred_temp_sepint_f, a_pred_temp_combint_f, a_pred_temp_noint_f,
                   a_pred_notemp_sepint_f, a_pred_notemp_combint_f, a_pred_notemp_noint_f)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

a_Female_DataFrame <- data.frame(Temp = TempRange,
                                  Q5 = apply(a_pred_f, 2, function(x) quantile(x, probs = c(0.05))),
                                  Q25 = apply(a_pred_f, 2, function(x) quantile(x, probs = c(0.25))),
                                  Q50 = apply(a_pred_f, 2, function(x) quantile(x, probs = c(0.50))),
                                  Q75 = apply(a_pred_f, 2, function(x) quantile(x, probs = c(0.75))),
                                  Q95 = apply(a_pred_f, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = a_Female_DataFrame, aes(x = Temp, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

################################################################################
### now males 
################################################################################

### a predictions, temp and sep interference

samples_quadtemp_sepint_m <- as.data.frame(samples_quadtemp_sepint_m[1:6])

a_pred_temp_sepint_m <- apply(samples_quadtemp_sepint_m, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, temp and combined interference

samples_quadtemp_combint_m <- as.data.frame(samples_quadtemp_combint_m[1:4])

a_pred_temp_combint_m <- apply(samples_quadtemp_combint_m, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, temp and no interference

samples_quadtemp_noint_m <-  as.data.frame(samples_quadtemp_noint_m[1:3])

a_pred_temp_noint_m <- apply(samples_quadtemp_noint_m, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, no temp sep interference

samples_notemp_sepint_m <- as.data.frame(samples_notemp_sepint_m[1:4])

a_pred_notemp_sepint_m <- matrix(rep(samples_notemp_sepint_m$a, 100), ncol = 4000, byrow = TRUE)

### a predictions, no temp combined interference

samples_notemp_combint_m <- as.data.frame(samples_notemp_combint_m[1:2])

a_pred_notemp_combint_m <- matrix(rep(samples_notemp_combint_m$a, 100), ncol = 4000, byrow = TRUE)

### a predictions, no temp no interference

samples_notemp_noint_m <- as.data.frame(samples_notemp_noint_m[1])

a_pred_notemp_noint_m <- matrix(rep(samples_notemp_noint_m$a, 100), ncol = 4000, byrow = TRUE)

### a predictions, null model -- there is no a, so it is zero?

a_pred_null_m <- matrix(data = 0, nrow = 100, ncol = 4000)


### sample predictions in accordance to waic weights

### number of samples that should be included from each of the models

round(waic_wts_m * 4000)

### sample from temp + sep interference

a_pred_temp_sepint_m <- apply(a_pred_temp_sepint_m, 1, function(x) sample(x, 56))

### sample from temp + combined interference

a_pred_temp_combint_m <- apply(a_pred_temp_combint_m, 1, function(x) sample(x, 201))

### sample from temp + no interference

a_pred_temp_noint_m <- apply(a_pred_temp_noint_m, 1, function(x) sample(x, 734))

### sample from no temp + sep interference

a_pred_notemp_sepint_m <- apply(a_pred_notemp_sepint_m, 1, function(x) sample(x, 175))

### sample from no temp + combined interference 

a_pred_notemp_combint_m <- apply(a_pred_notemp_combint_m, 1, function(x) sample(x, 431))

### sample from no temp + no interference

a_pred_notemp_noint_m <-  apply(a_pred_notemp_noint_m, 1, function(x) sample(x, 2295))

### sample from null

a_pred_null_m <- apply(a_pred_null_m, 1, function(x) sample(x, 108))


### now we can put all of these together into a single matrix

a_pred_m <- rbind(a_pred_temp_sepint_m, a_pred_temp_combint_m, a_pred_temp_noint_m,
                  a_pred_notemp_sepint_m, a_pred_notemp_combint_m, a_pred_notemp_noint_m,
                  a_pred_null_m)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

a_Male_DataFrame <- data.frame(Temp = TempRange,
                                 Q5 = apply(a_pred_m, 2, function(x) quantile(x, probs = c(0.05))),
                                 Q25 = apply(a_pred_m, 2, function(x) quantile(x, probs = c(0.25))),
                                 Q50 = apply(a_pred_m, 2, function(x) quantile(x, probs = c(0.50))),
                                 Q75 = apply(a_pred_m, 2, function(x) quantile(x, probs = c(0.75))),
                                 Q95 = apply(a_pred_m, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = a_Male_DataFrame, aes(x = Temp, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

################################################################################
### Last, the juveniles
################################################################################

### a predictions, temp and sep interference

samples_quadtemp_sepint_j <- as.data.frame(samples_quadtemp_sepint_j[1:6])

a_pred_temp_sepint_j <- apply(samples_quadtemp_sepint_j, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, temp and combined interference

samples_quadtemp_combint_j <- as.data.frame(samples_quadtemp_combint_j[1:4])

a_pred_temp_combint_j <- apply(samples_quadtemp_combint_j, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, temp and no interference

samples_quadtemp_noint_j <-  as.data.frame(samples_quadtemp_noint_j[1:3])

a_pred_temp_noint_j <- apply(samples_quadtemp_noint_j, 1, function(x) 
  x['ca']*exp(x['ba']*TempRange + x['qa']*TempRange^2))

### a predictions, no temp sep interference

samples_notemp_sepint_j <- as.data.frame(samples_notemp_sepint_j[1:4])

a_pred_notemp_sepint_j <- matrix(rep(samples_notemp_sepint_j$a, 100), ncol = 4000, byrow = TRUE)

### a predictions, no temp combined interference

samples_notemp_combint_j <- as.data.frame(samples_notemp_combint_j[1:2])

a_pred_notemp_combint_j <- matrix(rep(samples_notemp_combint_j$a, 100), ncol = 4000, byrow = TRUE)

### a predictions, no temp no interference

samples_notemp_noint_j <- as.data.frame(samples_notemp_noint_j[1])

a_pred_notemp_noint_j <- matrix(rep(samples_notemp_noint_j$a, 100), ncol = 4000, byrow = TRUE)

### number of samples that should be included from each of the models

round(waic_wts_j * 4000)

### sample from temp + sep interference

a_pred_temp_sepint_j <- apply(a_pred_temp_sepint_j, 1, function(x) sample(x, 847))

### sample from temp + combined interference

a_pred_temp_combint_j <- apply(a_pred_temp_combint_j, 1, function(x) sample(x, 589))

### sample from temp + no interference

a_pred_temp_noint_j <- apply(a_pred_temp_noint_j, 1, function(x) sample(x, 323))

### sample from no temp + sep interference

a_pred_notemp_sepint_j <- apply(a_pred_notemp_sepint_j, 1, function(x) sample(x, 880))

### sample from no temp + combined interference 

a_pred_notemp_combint_j <- apply(a_pred_notemp_combint_j, 1, function(x) sample(x, 901))

### sample from no temp + no interference

a_pred_notemp_noint_j <-  apply(a_pred_notemp_noint_j, 1, function(x) sample(x, 460))

### now we can put all of these together into a single matrix

a_pred_j <- rbind(a_pred_temp_sepint_j, a_pred_temp_combint_j, a_pred_temp_noint_j,
                  a_pred_notemp_sepint_j, a_pred_notemp_combint_j, a_pred_notemp_noint_j)

### create a data frame that contains midge densities, median, upper and lower quantiles for feeding rates

a_Juvenile_DataFrame <- data.frame(Temp = TempRange,
                               Q5 = apply(a_pred_j, 2, function(x) quantile(x, probs = c(0.05))),
                               Q25 = apply(a_pred_j, 2, function(x) quantile(x, probs = c(0.25))),
                               Q50 = apply(a_pred_j, 2, function(x) quantile(x, probs = c(0.50))),
                               Q75 = apply(a_pred_j, 2, function(x) quantile(x, probs = c(0.75))),
                               Q95 = apply(a_pred_j, 2, function(x) quantile(x, probs = c(0.95))))


ggplot(data = a_Juvenile_DataFrame, aes(x = Temp, y = Q50/24)) + geom_line() + 
  geom_line(aes(y = Q5/24)) + geom_line(aes(y = Q25/24)) + geom_line(aes(y = Q75/24)) + 
  geom_line(aes(y = Q95/24))

### now need to put all of the data frames together 

# Female_Temp_plot <- ggplot(data = a_Female_DataFrame, aes(x = Temp, y = Q50/24)) + geom_line() + 
#   geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24), alpha = 0.25, color = NA, fill = 'black') + 
#   geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24), alpha = 0.5, color = NA, fill = 'black') + 
#   theme_cowplot() + xlab('Temperature') + ylab('Space Clearance Rate \n(m^2 per hour)') + 
#   ggtitle('Females')
# 
# Male_Temp_plot <- ggplot(data = a_Male_DataFrame, aes(x = Temp, y = Q50/24)) + geom_line() + 
#   geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24), alpha = 0.25, color = NA, fill = 'black') + 
#   geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24), alpha = 0.5, color = NA, fill = 'black') + 
#   theme_cowplot() + xlab('Temperature') + ylab('Space Clearance Rate \n(m^2 per hour)') +
#   ggtitle('Males')
# 
# Juvenile_Temp_plot <- ggplot(data = a_Juvenile_DataFrame, aes(x = Temp, y = Q50/24)) + geom_line() + 
#   geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24), alpha = 0.25, color = NA, fill = 'black') + 
#   geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24), alpha = 0.5, color = NA, fill = 'black') + 
#   theme_cowplot() + xlab('Temperature') + ylab('Space Clearance Rate \n(m^2 per hour)') + 
#   ggtitle('Juveniles')
# 
# Combined_plot <- plot_grid(Female_Temp_plot, Male_Temp_plot, Juvenile_Temp_plot, nrow = 1, ncol = 3)
# 
# save_plot(filename = 'Temp_a_plot.png', plot = Combined_plot, nrow = 1, ncol = 3, base_asp = 1.3)
















