################################################################################
### Prey Abundance Corrections and calculating detection times
################################################################################

### load libraries

library(dplyr); library(tidyr);

### load correction data

corr_data <- read.csv('PreyCorrections_Simp.csv')

colnames(corr_data)[1] <- 'Observer'

### need to add up all of the observations then divide to get proportions
### of each of the transition types

corr_data <- corr_data %>% group_by(Observer) %>% select(-Survey) %>% summarise_all(.funs = sum)

corr_data <- corr_data %>% group_by(Observer) %>% mutate(Prop_S2M = SmalltoMedium/TotalSmall, 
                                                         Prop_M2L = MediumtoLarge/TotalMedium)

corr_data$Prop_M2L[2] <- 0

### Need to load the prey data to apply the corrections to the data before 6.6

prey_abundance <- read.csv('PreyAbundance.csv')

### manipulate to get just sizes and numbers for each size

prey_abundance <- prey_abundance %>% group_by(Date, SurveyCode, PreySize) %>% 
  summarise(Number = sum(Number))

### need to add observer to the surveys

### first load survey info.

Surv_info <- read.csv('SurveyInfo.csv', stringsAsFactors = TRUE)

Surv_info <- Surv_info %>% select(SurveyCode, OBS)

### join prey_abundance and Surv_info.

prey_abundance <- left_join(prey_abundance, Surv_info, by = 'SurveyCode')

colnames(prey_abundance)[5] <- 'Observer'

### split into surveys before and after the correction

before_corr <- prey_abundance %>% filter(Date < 6.6)

after_corr <- prey_abundance %>% filter(Date >= 6.6)

### make these into wide form

before_corr <- before_corr %>% pivot_wider(names_from = PreySize, values_from = Number, values_fill = 0)

after_corr <- after_corr %>% pivot_wider(names_from = PreySize, values_from = Number, values_fill = 0)

### now we need to apply the correction to the before_corr data


### join corr_data and before_corr by Observer after dropping some unnecessary stuff

corr_data <- corr_data %>% select(Observer, Prop_S2M, Prop_M2L)

before_corr <- left_join(before_corr, corr_data, by = 'Observer')

### now we can correct the prey abundances by subtracting proportions from the small 
### and medium and adding them to the medium and large 

before_corr <- before_corr %>% mutate(Small_change = Small*Prop_S2M, Med_change = Medium*Prop_M2L)

before_corr <- before_corr %>% mutate(Small = Small- Small_change, Medium = Medium + Small_change - Med_change, 
                                      Large = Large + Med_change)

### join the before_corr and after_corr data back together

### drop columns in before_corr that aren't in after_corr

before_corr <- before_corr %>% select(-c(Prop_S2M, Prop_M2L, Small_change, Med_change))

prey_abundance <- bind_rows(before_corr, after_corr)

###  add in the prey that were being eaten

spider_surv <- read.csv('IndSpiderSurveys.csv')

spider_surv <- spider_surv %>% group_by(SurveyCode, Date, Building, Wall, PreySize) %>% filter(PreySize != is.na(PreySize))  %>% summarise(Number = n()) %>% 
  pivot_wider(names_from = PreySize, values_from = Number, values_fill = 0)

### now join this with the prey abundance data frame

prey_abundance <- full_join(prey_abundance, spider_surv, by = 'SurveyCode')

### now want to add the numbers of prey together

### need to replace some NA's with zeros

prey_abundance[,c(4:6, 10:12)][is.na(prey_abundance[,c(4:6, 10:12)])] <- 0

### now can add the .x's and .y's together and drop a bunch of unnecessary columns

prey_abundance <- prey_abundance %>% mutate(Medium = Medium.x + Medium.y, 
                                            Small = Small.x + Small.y, 
                                            Large = Large.x + Large.y) %>% 
  mutate(Date.x = coalesce(Date.x, Date.y)) %>% 
  select(Date.x, SurveyCode, Medium, Small, Large)

### average proportions of each prey type

prey_average <- prey_abundance %>% group_by(SurveyCode) %>% mutate(Total = Medium + Small + Large) %>% summarise(Medium = mean(Medium)/mean(Total),
                                                                      Small = mean(Small)/mean(Total),
                                                                      Large = mean(Large)/mean(Total))
med_avg <- mean(prey_average$Medium)

small_avg <- mean(prey_average$Small)

large_avg <- mean(prey_average$Large)


### This gives us the total corrected prey densities! Whew.

###############################################################################################
### Use corrected abundances to calculate average detection times for each survey
###############################################################################################

### Now we want to use these corrected prey densities to get to mean detection time 
### estimates for each individual in the surveys if it were able to feed on all of the 
### individuals

### going to need to pair the individual spider survey data with the estimated 
### midge abundances

ind_spiders <- read.csv('IndSpiderSurveys.csv')

### only zebra jumpers

ind_spiders <- ind_spiders %>% filter(SpiderSpecies == 'Zebra')

### replace lengths of spiders with NA lengths with the mean for each sex/stage

mean_length_female <- ind_spiders %>% filter(SpiderSex == 'F') %>% select(SpiderLength)

mean_length_female <- mean(mean_length_female$SpiderLength, na.rm = TRUE)

mean_length_male <- ind_spiders %>% filter(SpiderSex == 'M') %>% select(SpiderLength)

mean_length_male <- mean(mean_length_male$SpiderLength, na.rm = TRUE)

mean_length_juvenile <- ind_spiders %>% filter(SpiderSex == 'J') %>% select(SpiderLength)

mean_length_juvenile <- mean(mean_length_juvenile$SpiderLength, na.rm = TRUE)

ind_spiders <- ind_spiders %>% mutate(SpiderLength_mean = ifelse(is.na(SpiderLength) & SpiderSex == 'F', mean_length_female,
                                                                 ifelse(is.na(SpiderLength) & SpiderSex == 'M', mean_length_male,
                                                                        ifelse(is.na(SpiderLength) & SpiderSex == 'J', mean_length_juvenile,
                                                                               SpiderLength))))
  
### now we join to add the prey densities

ind_spiders <- left_join(ind_spiders, prey_abundance, by = 'SurveyCode')


### need to add temperatures from the survey data

Surv_info <- read.csv('SurveyInfo.csv')

### drop all columns except for survey code and temperature

Surv_info <- Surv_info %>% select(SurveyCode, TempStart)

### now join the two data frames

ind_spiders <- left_join(x = ind_spiders, y = Surv_info, by = 'SurveyCode')

### now we need to calculate the detection times 
### the estimated detection time for a spider feeding on all of the different prey available
### is just a weighted average. (Nsmall*detection on small + Nmed*detection on med + NLarge * detection on large)/(Ntotal)

### need the coefficients of the model to figure out what detections for each spider should be

### model is E(log(detection)) = B0 + B1*log(SpiderLength) + B2*PreyMedium + B3*PreySmall + B4*log(Temperature)

### coefficients are:

B0 = 0.10; B1 = -0.71; B2 = -0.76; B3 = -1.37; B4 = -0.67;

ind_spiders <- ind_spiders %>% mutate(detect_small = exp(B0 + B1*log(SpiderLength_mean) + B3 + B4 * log(TempStart)) ,
                                      detect_med = exp(B0 + B1*log(SpiderLength_mean) + B2 + B4 * log(TempStart)),
                                      detect_large = exp(B0 + B1*log(SpiderLength_mean) + B4 * log(TempStart)))

### for each spider, we can take the weighted mean using the detection times and the prey abundances

ind_spiders <- ind_spiders %>% mutate(detection_est = (detect_small*Small + detect_med*Medium + detect_large*Large)/(Small + Medium + Large))

### take the mean within each survey

detection_predict <- ind_spiders %>% group_by(SurveyCode) %>% summarise(detection_est = mean(detection_est, na.rm = TRUE))

hist(detection_predict$detection_est)

### write this to a .csv file

write.csv(detection_predict, file = 'DetectionAvgAcrossInd.csv')

### seperate by sexes and then take averages

detection_predict_f <- ind_spiders %>% filter(SpiderSex == 'F') %>% group_by(SurveyCode) %>% summarise(detection_est = mean(detection_est))

detection_predict_m <- ind_spiders %>% filter(SpiderSex == 'M') %>% group_by(SurveyCode) %>% summarise(detection_est = mean(detection_est))                                                                                                       
                                                                                                       
detection_predict_j  <- ind_spiders %>% filter(SpiderSex == 'J') %>% group_by(SurveyCode) %>% summarise(detection_est = mean(detection_est))                                                                                                      

### write this to a .csv file

write.csv(detection_predict_f, file = 'DetectionAvgAcrossInd_f.csv')

write.csv(detection_predict_m, file = 'DetectionAvgAcrossInd_m.csv')

write.csv(detection_predict_j, file = 'DetectionAvgAcrossInd_j.csv')

### some histograms

hist(detection_predict_f$detection_est)

hist(detection_predict_m$detection_est)

hist(detection_predict_j$detection_est)
