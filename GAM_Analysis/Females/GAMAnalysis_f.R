################################################################################
### GAM analysis of feeding rates -- Females
################################################################################

### load libraries

library(mgcv); library(dplyr); library(ggplot2); library(cowplot); library(gratia);

### load data

females <- read.csv('FeedingData_Manipulated_f.csv', stringsAsFactors = TRUE)

################################################################################
### Goal is to fit GAM models to the data and use model selection to pick 
### which model best describes the feeding rates of females. Will then fit 
### a parametric model that matches the GAM results
################################################################################

### first need to calculate the feeding rates for each of the surveys

females <- females %>% mutate(FeedingRate = ProportionFeed/(DetectionEst*24))

### create a Total Predator Density variable 

females <- females %>% mutate(TotalPredDens = FemaleDens + MaleDens + JuvDens)

### drop surveys with no midges

females <- females %>% filter(NumberMidges > 0)


### full model -- temperature, separate interference from each stage/sex

full_fit <- gam(FeedingRate ~  s(NumberMidges) + s(FemaleDens) + s(JuvDens) +
                s(MaleDens) + s(TempStart) + s(BuildingWall, bs = 're'), 
                weights = females$NumberPred/mean(females$NumberPred), data = females, method = 'REML')

summary(full_fit)

AIC(full_fit)

### interference combined 

intcomb_fit <- gam(FeedingRate ~ s(NumberMidges) + s(TotalPredDens) + s(TempStart) + s(BuildingWall, bs = 're'),
                   weights = females$NumberPred/mean(females$NumberPred), data = females, method = 'REML')

summary(intcomb_fit)

AIC(full_fit, intcomb_fit)

### temperature no interference

noint_fit <- gam(FeedingRate ~ s(NumberMidges) + s(TempStart) + s(BuildingWall, bs = 're'), 
                 weights = females$NumberPred/mean(females$NumberPred), data = females, method = 'REML')

summary(noint_fit)

AIC(full_fit, intcomb_fit, noint_fit)

### no temperature + interference separate

notemp_fit <- gam(FeedingRate ~ s(NumberMidges) + s(FemaleDens) + s(JuvDens) + 
                  s(MaleDens) + s(BuildingWall, bs = 're'), 
                  weights = females$NumberPred/mean(females$NumberPred), data = females, method = 'REML')

AIC(full_fit, intcomb_fit, noint_fit, notemp_fit)

### no temperature + interference combined

notemp_intcomb_fit <- gam(FeedingRate ~ s(NumberMidges) + s(TotalPredDens) + s(BuildingWall, bs = 're'),
                          weights = females$NumberPred/mean(females$NumberPred), data = females, method = 'REML')

summary(notemp_intcomb_fit)

AIC(full_fit, intcomb_fit, notemp_fit, notemp_intcomb_fit)

### no temperature + no interference

notemp_noint_fit <- gam(FeedingRate ~ s(NumberMidges) + s(BuildingWall, bs = 're'), 
                        weights = females$NumberPred/mean(females$NumberPred), data = females, method = 'REML')

summary(notemp_noint_fit)

AIC(full_fit, intcomb_fit, notemp_fit, notemp_intcomb_fit, notemp_noint_fit)

### can we do a null model?

null <- gam(FeedingRate ~ 1 + s(BuildingWall, bs = 're'), 
            weights = females$NumberPred/mean(females$NumberPred), data = females, method = 'REML')

summary(null)

AICnull <- AIC(null)

### as one might expect for the female data

f_AicTable <- as.data.frame(AIC(full_fit, intcomb_fit, noint_fit, notemp_fit, notemp_intcomb_fit, notemp_noint_fit, null))

f_AicTable$Model <- c("Full", "IntComb", "NoInt", "NoTemp", "NoTemp-IntComb", "NoTemp-NoInt", "Null")

f_AicTable <- f_AicTable %>% mutate(DeltaAIC = round(AIC, 2) - 301.42)
f_AicTable <- f_AicTable %>% mutate(AIC_weight = round(exp(-0.5 * DeltaAIC), 4))

f_AicTable <- f_AicTable %>% mutate(AIC_weight = AIC_weight/sum(AIC_weight))

f_AicTable <- f_AicTable %>% arrange(desc(AIC_weight))

f_AicTable




