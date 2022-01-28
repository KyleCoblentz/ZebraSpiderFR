################################################################################
### GAM Analysis of feeding rates -- Juveniles
################################################################################

### load libraries

library(dplyr); library(ggplot2); library(cowplot); library(mgcv); library(gratia);

### load data

juveniles <- read.csv('FeedingData_Manipulated_j.csv', stringsAsFactors = TRUE)

################################################################################
### fit GAMs to get general relationship between feeding rates and 
### and potential functional response variables. Will fit a parametric 
### functional response afterwards to get FR parameter estimates
################################################################################

### manipulate data to get feeding rates and a Total Predator Density Variable

juveniles <- juveniles %>% mutate(FeedingRate = ProportionFeed/(DetectionEst*24), 
                                  TotalPredDens = JuvenileDens + FemaleDens + MaleDens)

### full model -- temperature, separate interference from each stage/sex

full_fit <- gam(FeedingRate ~  s(NumberMidges) + s(FemaleDens) + s(JuvenileDens) + 
                  s(MaleDens) + s(TempStart) + s(BuildingWall, bs = 're'), 
                weights = juveniles$NumberPred/mean(juveniles$NumberPred), data = juveniles, method = 'REML')

summary(full_fit)

draw(full_fit, residuals = TRUE)

AIC(full_fit)

### interference combined

intcomb_fit <- gam(FeedingRate ~  s(NumberMidges) + s(TotalPredDens) + s(TempStart) + s(BuildingWall, bs = 're'), 
                   weights = juveniles$NumberPred/mean(juveniles$NumberPred), data = juveniles, method = 'REML')

summary(intcomb_fit)

draw(intcomb_fit, residuals = TRUE)

AIC(full_fit, intcomb_fit)

### no interference, temperature

noint_fit <- gam(FeedingRate ~ s(NumberMidges) + s(TempStart) + s(BuildingWall, bs = 're'), 
                 weights = juveniles$NumberPred/mean(juveniles$NumberPred), data = juveniles, method = 'REML')

summary(noint_fit)

draw(noint_fit, residuals = TRUE)

AIC(full_fit, intcomb_fit, noint_fit)

### no temperature, interference separate

notemp_fit <- gam(FeedingRate ~  s(NumberMidges) + s(FemaleDens) + s(JuvenileDens) + 
                    s(MaleDens) +  s(BuildingWall, bs = 're'), 
                  weights = juveniles$NumberPred/mean(juveniles$NumberPred), data = juveniles, method = 'REML')

summary(notemp_fit)

draw(notemp_fit, residuals = TRUE)

AIC(full_fit, intcomb_fit, noint_fit, notemp_fit)

### no temperature + interference combined

notemp_intcomb_fit <- gam(FeedingRate ~  s(NumberMidges) + s(TotalPredDens) +  s(BuildingWall, bs = 're'),
                          weights = juveniles$NumberPred/mean(juveniles$NumberPred), data = juveniles, method = 'REML')

summary(notemp_intcomb_fit)

appraise(notemp_intcomb_fit)

AIC(full_fit, intcomb_fit, noint_fit, notemp_fit, notemp_intcomb_fit)

### no temperature, no interference

notemp_noint_fit <- gam(FeedingRate ~  s(NumberMidges)  +  s(BuildingWall, bs = 're'),
                        weights = juveniles$NumberPred/mean(juveniles$NumberPred), data = juveniles, method = 'REML')

summary(notemp_noint_fit)

draw(notemp_noint_fit, residuals = TRUE)

AIC(full_fit, intcomb_fit, noint_fit, notemp_fit, notemp_intcomb_fit, notemp_noint_fit)

### null model

null <- gam(FeedingRate ~ 1 + s(BuildingWall, bs = 're'), 
            weights = juveniles$NumberPred/mean(juveniles$NumberPred), data = juveniles, method = 'REML')

summary(null)

AIC(null)

### make a table of the AIC results

f_AicTable <- as.data.frame(AIC(full_fit, intcomb_fit, noint_fit, notemp_fit, notemp_intcomb_fit, notemp_noint_fit, null))

f_AicTable$Model <- c("Full", "IntComb", "NoInt", "NoTemp", "NoTemp-IntComb", "NoTemp-NoInt", "Null")

f_AicTable <- f_AicTable %>% mutate(DeltaAIC = round(AIC, 2) - 220.94)

f_AicTable <- f_AicTable %>% mutate(AIC_weight = round(exp(-0.5 * DeltaAIC), 4))

f_AicTable <- f_AicTable %>% mutate(AIC_weight = AIC_weight/sum(AIC_weight))

f_AicTable <- f_AicTable %>% arrange(desc(AIC_weight))

f_AicTable


