#################################################################################
### Data manipulation
#################################################################################

### need to pull together data into a dataframe we can use for the analysis

### load packages

library(dplyr); library(tidyr); library(ggplot2); library(cowplot);

### upload feeding data

feed <- read.csv('IndSpiderSurveys.csv', stringsAsFactors = TRUE)

### drop some unnecessary columns

feed <- feed[,-c(9:12,14:17)]

### split data by stage and sex. Zebra Jumpers only.

feed_f <- feed %>% filter(SpiderSex == 'F', SpiderSpecies == 'Zebra')

feed_m <- feed %>% filter(SpiderSex == 'M', SpiderSpecies == 'Zebra')

feed_j <- feed %>% filter(SpiderSex == 'J', SpiderSpecies == 'Zebra')

### need to summarise data to one line per survey with all of the info
### that we will need for analysis

### lets first just create a column that has a 1 if individual is feeding on a midge 
### and a zero if it isn't

feed_f <- feed_f %>% mutate(FeedMidge = ifelse(is.na(PreySize), 0, 1))

feed_m <- feed_m %>% mutate(FeedMidge = ifelse(is.na(PreySize), 0, 1))

feed_j <- feed_j %>% mutate(FeedMidge = ifelse(is.na(PreySize), 0, 1))


### summarise data to number feeding on midges, total number of female spiders, 
### the proportion feeding on midges

feed_f <- feed_f %>% group_by(SurveyCode, Date, Building, Wall) %>% summarise(NumberFeed = sum(FeedMidge),
                                                                                NumberPred = n(),
                                                                                ProportionFeed = NumberFeed/NumberPred,
                                                                                SpiderLength = mean(SpiderLength, na.rm = TRUE))

feed_m <- feed_m %>% group_by(SurveyCode, Date, Building, Wall) %>% summarise(NumberFeed = sum(FeedMidge),
                                                                              NumberPred = n(),
                                                                              ProportionFeed = NumberFeed/NumberPred,
                                                                              SpiderLength = mean(SpiderLength, na.rm = TRUE))

feed_j <- feed_j %>% group_by(SurveyCode, Date, Building, Wall) %>% summarise(NumberFeed = sum(FeedMidge),
                                                                              NumberPred = n(),
                                                                              ProportionFeed = NumberFeed/NumberPred,
                                                                              SpiderLength = mean(SpiderLength, na.rm = TRUE))

### we have 150 surveys for females, 113 surveys for juveniles, and 96 surveys for males

### need to replace some NaN's in the spider length data

feed_f$SpiderLength <-  ifelse(is.na(feed_f$SpiderLength), mean(feed_f$SpiderLength, na.rm = TRUE), feed_f$SpiderLength)

feed_m$SpiderLength <-  ifelse(is.na(feed_m$SpiderLength), mean(feed_m$SpiderLength, na.rm = TRUE), feed_m$SpiderLength)

feed_j$SpiderLength <-  ifelse(is.na(feed_j$SpiderLength), mean(feed_j$SpiderLength, na.rm = TRUE), feed_j$SpiderLength)

### Now want to add columns for number of males, number of each other spider species,
### and  the temperature for each survey

SpiderNumbers <- feed %>% group_by(SurveyCode) %>% summarise(ZebraMale = sum(SpiderSpecies == 'Zebra' & SpiderSex == 'M'),
                                                                 ZebraFemale = sum(SpiderSpecies == 'Zebra' & SpiderSex == 'F'),
                                                                 ZebraJuvenile = sum(SpiderSpecies == 'Zebra' & SpiderSex == 'J'),
                                                                 Bold = sum(SpiderSpecies == 'Bold'),
                                                                 Unknown = sum(SpiderSpecies == 'Unknown'), 
                                                                 Tan = sum(SpiderSpecies == 'Tan'),
                                                                 BHLB = sum(SpiderSpecies == 'BHLB'),
                                                                 TotalOtherSpider = sum(Bold, Unknown, Tan, BHLB))

feed_f <- left_join(feed_f, SpiderNumbers, by = 'SurveyCode')

feed_m <- left_join(feed_m, SpiderNumbers, by = 'SurveyCode')

feed_j <- left_join(feed_j, SpiderNumbers, by = 'SurveyCode')

### Need to add temperature to surveys

SurveyInfo <- read.csv('SurveyInfo.csv', stringsAsFactors = TRUE)

### drop everything except temperature and survey code then join to feed_nm

SurveyInfo <- SurveyInfo %>% mutate(TempStart = ifelse(is.na(TempStartC), (5/9)*(TempStartF-32), TempStartC))  

SurveyInfo <- SurveyInfo %>% select(SurveyCode, TempStart)

feed_f <- left_join(feed_f, SurveyInfo, by = "SurveyCode")

feed_m <- left_join(feed_m, SurveyInfo, by = "SurveyCode")

feed_j <- left_join(feed_j, SurveyInfo, by = "SurveyCode")

### Now need to get wall sizes to get densities

wall_sizes <- read.csv("WallInfo.csv")

### drop bathhouse from data becuase of unknown area, create buildingwall variable

feed_f <- feed_f %>% filter(Building != 'BH') %>% mutate(BuildingWall = paste0(Building, Wall))

feed_m <- feed_m %>% filter(Building != 'BH') %>% mutate(BuildingWall = paste0(Building, Wall))

feed_j <- feed_j %>% filter(Building != 'BH') %>% mutate(BuildingWall = paste0(Building, Wall))

### create building wall in wall_size

wall_sizes <- wall_sizes %>% mutate(BuildingWall = paste0(Abbreviation, Wall.Direction))

### Change name of wall area

wall_sizes <- wall_sizes %>% rename(TotalArea_m2 = Total.Area..m.2.)

wall_sizes$TotalArea_m2 <- as.numeric(paste(wall_sizes$TotalArea_m2))

### add wall area to data

wall_sizes <- wall_sizes %>% select(BuildingWall, TotalArea_m2)

### change a couple of the names to match the spider survey data

wall_sizes$BuildingWall[which(wall_sizes$BuildingWall == 'DBNW')] <- 'DHNW'

wall_sizes$BuildingWall[which(wall_sizes$BuildingWall == 'TBNW')] <- 'TBN'

feed_f <- left_join(feed_f, wall_sizes, by = 'BuildingWall')

feed_m <- left_join(feed_m, wall_sizes, by = "BuildingWall")

feed_j <- left_join(feed_j, wall_sizes, by = "BuildingWall")

### load detection time data, join to spider survey data

detectiondata_f <- read.csv('DetectionAvgAcrossInd_f.csv', stringsAsFactors = TRUE)

feed_f <- left_join(feed_f, detectiondata_f, by = "SurveyCode")

detectiondata_m <- read.csv('DetectionAvgAcrossInd_m.csv', stringsAsFactors = TRUE)

feed_m <- left_join(feed_m, detectiondata_m, by = 'SurveyCode')

detectiondata_j <- read.csv('DetectionAvgAcrossInd_j.csv', stringsAsFactors = TRUE)

feed_j <- left_join(feed_j, detectiondata_j, by = 'SurveyCode')

### set NA detection times to the mean across all surveys 

meanDetection_f <- mean(detectiondata_f$detection_est, na.rm = TRUE)

feed_f <- feed_f %>% mutate(DetectionEst = ifelse(is.na(detection_est), meanDetection_f, detection_est))

meanDetection_m <- mean(detectiondata_m$detection_est, na.rm = TRUE)

feed_m <- feed_m %>% mutate(DetectionEst = ifelse(is.na(detection_est), meanDetection_m, detection_est))

meanDetection_j <- mean(detectiondata_j$detection_est, na.rm = TRUE)

feed_j <- feed_j %>% mutate(DetectionEst = ifelse(is.na(detection_est), meanDetection_j, detection_est))

### Now add the prey abundances

### load the prey abundance data

preydata <- read.csv('PreyAbundance.csv', stringsAsFactors = TRUE)

### manipulate prey data to sum up all of the midges

preydata <- preydata %>% group_by(SurveyCode) %>% summarise(NumberMidges = sum(Number))

### join with the data

feed_f <- left_join(feed_f, preydata, by = 'SurveyCode')

feed_m <- left_join(feed_m, preydata, by = "SurveyCode")

feed_j <- left_join(feed_j, preydata, by = "SurveyCode")

### need to change the NA's to the number that are being eaten by the predators 
### add the number being eaten by the predators to the number not being eaten

feed <- feed %>% mutate(FeedMidge = ifelse(is.na(PreySize), 0, 1)) %>% filter(Building != 'BH')

feed_nfeed <- feed %>% group_by(SurveyCode) %>% summarise(TotalSpiderFeed = sum(FeedMidge))

### join to add the total number of spiders feeding 

feed_f <- left_join(feed_f, y = feed_nfeed, by = 'SurveyCode')

feed_m <- left_join(feed_m, y = feed_nfeed, by = 'SurveyCode')

feed_j <- left_join(feed_j, y = feed_nfeed, by = 'SurveyCode')

### add the number of midges being eaten 

feed_f <- feed_f %>% mutate(NumberMidges = ifelse(is.na(NumberMidges), 0, NumberMidges))

feed_f <- feed_f %>% mutate(NumberMidges = NumberMidges + TotalSpiderFeed)

feed_m <- feed_m %>% mutate(NumberMidges = ifelse(is.na(NumberMidges), 0, NumberMidges))

feed_m <- feed_m %>% mutate(NumberMidges = NumberMidges + TotalSpiderFeed)

feed_j <- feed_j %>% mutate(NumberMidges = ifelse(is.na(NumberMidges), 0, NumberMidges))

feed_j <- feed_j %>% mutate(NumberMidges = NumberMidges + TotalSpiderFeed)

### convert all of the numbers to densities

feed_f <- feed_f %>% mutate(PredDens = (NumberPred)/TotalArea_m2, 
                            JuvDens = ZebraJuvenile/TotalArea_m2, FemaleDens = (ZebraFemale-1)/TotalArea_m2,
                              MaleDens = ZebraMale/TotalArea_m2, Bold = Bold/TotalArea_m2,
                              Unknown = Unknown/TotalArea_m2, Tan = Tan/TotalArea_m2, 
                              BHLB = BHLB/TotalArea_m2, TotalOtherSpider = TotalOtherSpider/TotalArea_m2,
                              NumberMidges = NumberMidges/TotalArea_m2)

feed_m <- feed_m %>% mutate(MaleDens = (ZebraMale)/TotalArea_m2, 
                            FemaleDens = ZebraFemale/TotalArea_m2,
                            JuvDens = ZebraJuvenile/TotalArea_m2, Bold = Bold/TotalArea_m2,
                            Unknown = Unknown/TotalArea_m2, Tan = Tan/TotalArea_m2, 
                            BHLB = BHLB/TotalArea_m2, TotalOtherSpider = TotalOtherSpider/TotalArea_m2,
                            NumberMidges = NumberMidges/TotalArea_m2)

feed_j <- feed_j %>% mutate(JuvenileDens = (ZebraJuvenile)/TotalArea_m2, 
                            FemaleDens = ZebraFemale/TotalArea_m2,
                            MaleDens = ZebraMale/TotalArea_m2, Bold = Bold/TotalArea_m2,
                            Unknown = Unknown/TotalArea_m2, Tan = Tan/TotalArea_m2, 
                            BHLB = BHLB/TotalArea_m2, TotalOtherSpider = TotalOtherSpider/TotalArea_m2,
                            NumberMidges = NumberMidges/TotalArea_m2)

### save to .csv files

 write.csv(feed_f, 'FeedingData_Manipulated_f.csv')

 write.csv(feed_m, 'FeedingData_Manipulated_m.csv')

 write.csv(feed_j, 'FeedingData_Manipulated_j.csv')

