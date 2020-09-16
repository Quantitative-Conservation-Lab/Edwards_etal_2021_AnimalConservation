#-----------------Set working directory--------------------
rm(list = ls())

#setwd("/Users/hannahedwards/Documents/Calgary/CZ/Incubation Study")
setwd("C:/Users/HannahE/OneDrive - The Calgary Zoological Society/My Documents/Incubation Study/")

### Set the library path (just in case, as it seems a bit unpredictable at times)
#.libPaths("C:/R/Library")


#----------------Load libraries.-------------
require(tidyverse)
library(stringr)
library(xts)

library(plyr)
library(dplyr)
library(gdata)

# Import the data for each egg (no data collected in 2018).
###PWRCC 17
egg1 <- read.csv("Patuxent/2017/P1_R34no1.csv")
egg2 <- read.csv("Patuxent/2017/P1_S17no2.csv")
egg3 <- read.csv("Patuxent/2017/P2_B5no6.csv")
egg4 <- read.csv("Patuxent/2017/P3_R34no2.csv")
egg5 <- read.csv("Patuxent/2017/P3_S17no1.csv")
egg6 <- read.csv("Patuxent/2017/P3_Y16no2.csv")
egg7 <- read.csv("Patuxent/2017/P4_Y16no1.csv")
egg8 <- read.csv("Patuxent/2017/P5_S32no2.csv")
egg9 <- read.csv("Patuxent/2017/P6_Y16no1g.csv")
egg10 <- read.csv("Patuxent/2017/P7_S32no1.csv")
###PWRCC 16
egg11 <- read.csv("Patuxent/2016/P1_S18no.4.2.csv")
egg12 <- read.csv("Patuxent/2016/P1_Y27no.2.csv")
egg13 <- read.csv("Patuxent/2016/P2_S18no.3.csv")
egg14 <- read.csv("Patuxent/2016/P3_Y17no.1b.csv")
egg15 <- read.csv("Patuxent/2016/P4_B6no.3.csv")
egg16 <- read.csv("Patuxent/2016/P6_Y17no.1g.csv")
egg17 <- read.csv("Patuxent/2016/P7_Y17no.2b.csv")
egg18 <- read.csv("Patuxent/2016/P9_Y27no.1.csv")
egg19 <- read.csv("Patuxent/2016/T3_B6no.4.csv")
###PWRCC 15
egg20 <- read.csv("Patuxent/2015/P1_T7no.2.csv")
egg21 <- read.csv("Patuxent/2015/P2_R34no.4.csv")
egg22 <- read.csv("Patuxent/2015/P3_R34no.3.csv")
egg23 <- read.csv("Patuxent/2015/P3_S2no.5.csv")
egg24 <- read.csv("Patuxent/2015/P4_Y34no.3.csv")
egg25 <- read.csv("Patuxent/2015/P5_S17no.3.csv")
egg26 <- read.csv("Patuxent/2015/P5_Y42no.4.csv")
egg27 <- read.csv("Patuxent/2015/P6_R30no.2.csv")
egg28 <- read.csv("Patuxent/2015/P8_R6no.6.csv")
egg29 <- read.csv("Patuxent/2015/P8_S17no.1.csv")
egg30 <- read.csv("Patuxent/2015/P9_R6no.5.csv")
egg31 <- read.csv("Patuxent/2015/P9_S2no.2.csv")
egg32 <- read.csv("Patuxent/2015/P9_S17no.4.csv")
egg33 <- read.csv("Patuxent/2015/P9_Y34no.4.csv")
egg34 <- read.csv("Patuxent/2015/P9_Y42no.5.csv")
egg35 <- read.csv("Patuxent/2015/P10_R30no.1.csv")
egg36 <- read.csv("Patuxent/2015/P10_S2no.6.csv")
egg37 <- read.csv("Patuxent/2015/P10_S17no.2.csv")
###PWRCC 14
egg38 <- read.csv("Patuxent/2014/P3_R23.1_10 days long.csv")
egg39 <- read.csv("Patuxent/2014/P7_R23.4.3.csv")
egg40 <- read.csv("Patuxent/2014/P9_R23.2_11 days long.csv")
egg41 <- read.csv("Patuxent/2014/P9_R33.2.3_8 days long.csv")
egg42 <- read.csv("Patuxent/2014/P9_S18.2.csv")
egg43 <- read.csv("Patuxent/2014/P10_S1.2_4 days long.csv")
#egg44 <- read.csv("Patuxent/2014/P10_S18.4.csv") removed because it received SHC and Petersime treatment
####CZ 19
egg45<-read.csv("Calgary Zoo/2019/CZ-3_PS2-2-19.csv")
#egg46<-read.csv("Calgary Zoo/2019/CZ-2_BL1-1-19.csv")#Not included because it failed in the wild, so outcome in captivity is unclear
####CZ 17
egg47 <- read.csv("Calgary Zoo/2017/CZ-1_IN3-1-17.csv")
egg48 <- read.csv("Calgary Zoo/2017/CZ-4_BL3-1-17.csv")
####CZ 16
egg49 <- read.csv("Calgary Zoo/2016/CZ1_MD2-2-16.csv")
egg50 <- read.csv("Calgary Zoo/2016/CZ2_A1-2-16.csv")
egg51 <- read.csv("Calgary Zoo/2016/CZ2_MD2-1-16.csv")
egg52 <- read.csv("Calgary Zoo/2016/CZ3_N2-1-16.csv")
egg53 <- read.csv("Calgary Zoo/2016/CZ3_PS2-2-16.csv")
egg54 <- read.csv("Calgary Zoo/2016/CZ4_DO1-1-16_16 days long.csv")
egg55 <- read.csv("Calgary Zoo/2016/CZ4_DO1-2-16.csv")
egg56 <- read.csv("Calgary Zoo/2016/CZ4_MD1-2-16_14 days long.csv")
egg57 <- read.csv("Calgary Zoo/2016/CZ5_PS2-1-16_15 days long.csv")
####CZ 15
egg58 <- read.csv("Calgary Zoo/2015/CZ1_N1-1-15_14 days long.csv")
egg59 <- read.csv("Calgary Zoo/2015/CZ3_N4-1-15.csv")
egg60 <- read.csv("Calgary Zoo/2015/CZ3_PS3-2-15_14 days long.csv")
egg61 <- read.csv("Calgary Zoo/2015/CZ4_N4-2-15.csv")
####CZ 14
egg62 <- read.csv("Calgary Zoo/2014/CZ2_MD2-1-14.csv")
egg63 <- read.csv("Calgary Zoo/2014/CZ2_MD2-2-14.csv")
egg64 <- read.csv("Calgary Zoo/2014/CZ2_N2-3-14_3 days long.csv")
####ICF 17
egg65 <- read.csv("ICF/2017/I-8_13-144-04-17.csv")
egg66 <- read.csv("ICF/2017/I-10_13-116-06-17.csv")
egg67 <- read.csv("ICF/2017/P-11_07-117-01-17.csv")
####ICF 16
egg68 <- read.csv("ICF/2016/I8_13-46-04-16.csv")
egg69 <- read.csv("ICF/2016/I14_13-35-05-16.csv")
egg70 <- read.csv("ICF/2016/I19_07-117-05-16.csv")
####ICF 15
egg71 <- read.csv("ICF/2015/I8+P11_07-114-02-15.csv")
egg72 <- read.csv("ICF/2015/I8_07-117-03-15.csv")
egg73 <- read.csv("ICF/2015/I8_07-117-04-15.csv")
egg74 <- read.csv("ICF/2015/I14_07-117-05-15.csv")
egg75 <- read.csv("ICF/2015/I19_13-46-03-15.csv")

egg_IDs<-read.csv("Summary data/egg_IDs.csv")

# Combine all the data into one list for analysis.
eggs_all <- list(egg1, egg2,egg3,egg4,egg5,egg6,egg7,egg8,egg9,egg10,
                 egg11, egg12,egg13,egg14,egg15,egg16,egg17,egg18,egg19,egg20,
                 egg21, egg22,egg23,egg24,egg25,egg26,egg27,egg28,egg29,egg30,
                 egg31, egg32,egg33,egg34,egg35,egg36,egg37,egg38,egg39,egg40,
                 egg41, egg42,egg43,egg45,egg47,egg48,egg49,egg50,
                 egg51, egg52,egg53,egg54,egg55,egg56,egg57,egg58,egg59,egg60,
                 egg61, egg62,egg63,egg64,egg65,egg66,egg67,egg68,egg69,egg70,
                 egg71, egg72,egg73, egg74, egg75)

#For each egg create a new dataset with the 40th datapoint every 
#24 hours starting from the first data point recorded each day and put it in egg_data1 
egg_data1 = list()
for(i in 1:length(eggs_all)){
  tmp <- eggs_all[[i]]
  tmp$sequential_day <- tmp$Day - min(tmp$Day) 
  tmp$time_in_minutes <- (tmp$sequential_day * 24 * 60) +(tmp$Hour*60) + tmp$Min
  
  first_time <- tmp[1, "time_in_minutes"]
  lags_40min <-  cumsum(c(first_time, rep(40, (24*60 *(max(tmp$Day)-min(tmp$Day))))))
  lags_40min <- lags_40min[which(lags_40min <= max(tmp$time_in_minutes))]
  lags <- cumsum(c(first_time, rep(24*60, max(tmp$Day)-min(tmp$Day))))
  lags <- lags[which(lags <= max(tmp$time_in_minutes))]
  
  #S17#1, S17#2,S2#2, Y42#4, Y42#5-have a small time gap (couple hours) so do not have data every 40 mins during that time
  if(length(which(tmp$time_in_minutes %in% lags))== length(lags)){
    tmp[which(tmp$time_in_minutes %in% lags), "lags_24_hours"] <- lags
  } else {
    tmp[which(tmp$time_in_minutes %in% lags), "lags_24_hours"] <- lags[which(lags %in% tmp$time_in_minutes)]
  }
  tmp <- tmp %>% fill(lags_24_hours, .direction="down")
  tmp <- tmp[tmp$time_in_minutes %in% lags_40min, ]
  Day <- tmp[1, "Day"]
  egg_data1[[i]] <- tmp
} 

#Calculate the mean and variance 
egg_data2 = list()
for(j in 1:length(egg_data1)){
  egg_summs <- egg_data1[[j]]
  #Create temp mean and var across the three RH sensors
  egg_summs$Temp1.<-as.numeric(as.character(egg_summs$Temp1.))
  egg_summs$Temp2.<-as.numeric(as.character(egg_summs$Temp2.))
  egg_summs$Temp3.<-as.numeric(as.character(egg_summs$Temp3.))
    
  #Then from the thinned data (egg data1) calc for each egg the mean and var 
  #of each environ param for each day
  egg_summs$temp_mean <- apply(cbind(egg_summs$Temp1., egg_summs$Temp2., egg_summs$Temp3.), 1, mean)
  
  #Create RH mean and var across the two RH sensors
  egg_summs$RH1. <- as.numeric(as.character(gsub("%", "", egg_summs$RH1.)))
  egg_summs$RH2. <- as.numeric(as.character(gsub("%", "", egg_summs$RH2.)))
  # Take the average of the two humidity readings for each minute
  egg_summs$rh_mean <- apply(cbind(egg_summs$RH1., egg_summs$RH2.), 1, mean)

  #Convert position numeric
  egg_summs$Position<-as.numeric(as.character(egg_summs$Position))
  #Calc the min difference between the current position and next position
  nrows <- nrow(egg_summs)[[1]]
  egg_summs$position_change[2:nrows] <- pmin( abs(egg_summs$Position[2:nrows] - egg_summs$Position[1:nrows-1]),
  360-abs(egg_summs$Position[2:nrows] - egg_summs$Position[1:nrows-1]))
  egg_summs2 <- egg_summs
  #Calc the mean and var per day for each environ variable.
  if(any(names(egg_summs)=="EggID")){
    egg_summs$Egg.ID <- egg_summs$EggID
  }
  Daily_temp_mean<-aggregate(temp_mean~Day+Egg.ID, data=egg_summs, FUN='mean')
  Daily_temp_var<-aggregate(temp_mean~Day, data=egg_summs, FUN='var')
  Daily_temp_var<-rename(Daily_temp_var, replace= c("temp_mean"="temp_var"))
  Daily_rh_mean<-aggregate(rh_mean~Day, data=egg_summs, FUN='mean')
  Daily_rh_var<-aggregate(rh_mean~Day, data=egg_summs, FUN='var')
  Daily_rh_var<-rename(Daily_rh_var, replace= c("rh_mean"="rh_var"))
  Daily_rot_mean<-aggregate(position_change~Day, data=egg_summs2, FUN='mean')
  Daily_rot_mean<-rename(Daily_rot_mean, replace= c("position_change"="rot_mean"))
  Daily_rot_var<-aggregate(position_change~Day, data=egg_summs2, FUN='var')
  Daily_rot_var<-rename(Daily_rot_var, replace= c("position_change"="rot_var"))
  #Merge
  Daily_rot<-merge(Daily_rot_mean, Daily_rot_var, by="Day")
  egg_data_summs<-merge(Daily_rot, Daily_temp_mean, by="Day", all=TRUE)
  egg_data_summs<-cbind(egg_data_summs, Daily_temp_var$temp_var, Daily_rh_mean$rh_mean, 
                        Daily_rh_var$rh_var)
  #Rename
  egg_data_summs<-rename(egg_data_summs, c("Daily_temp_var$temp_var"="temp_var", "Daily_rh_mean$rh_mean"="rh_mean",
                                           "Daily_rh_var$rh_var"="rh_var"))
  
  egg_data2[[j]] <- egg_data_summs
  
  print(j)
}

  # Turn the lists into a single dataframe after the loop and write as csv
  egg_data2 = do.call(rbind,egg_data2)
  
  #Add paired date from another summary datasheet-not included in individual csvs
  Pair<-read.csv("Summary data/pair day.csv")
  Paired.date<-subset(Pair, select=c(Egg.ID, Paired.date))
  egg_data2<-merge(egg_data2, Paired.date, by="Egg.ID")
  
  #Add unpaired date column to egg_data2
  egg_data2$Unpaired.date<- NA
  #Move unpaired data from eggs_all to egg_data2
  for(k in 1:length(eggs_all)){
    nrows <- dim(eggs_all[[k]])[[1]]
    egg <- eggs_all[[k]]$Egg.ID[1] 
    index <- which(match(egg_data2$Egg.ID, egg) == 1) 
    
    # Extract the unpaired date and add it to the dataframe
    Unpaired.date <- eggs_all[[k]]$Unpaired.date[1] 
    egg_data2$Unpaired.date[index] <- as.character(Unpaired.date)
  }
  
  #write the csv
  write.csv(egg_data2, file = "Egg_summary_survival_dailydata_40min.csv", row.names = F)

#------ Initialize a column in our predictor dataframe for the things we will add for each egg
# As we calculate summary data below for each egg we will add it to this data frame.
egg_summary <- egg_IDs # initialize as the egg outcomes
colnames(egg_summary)[[1]] <- 'Egg.ID'
  
egg_summary$Facility<- NA
egg_summary$Fate <- NA
egg_summary$Year <- NA
egg_summary$Treatment <- NA
egg_summary$Logger.ID <- NA
egg_summary$Lay.date <- NA
egg_summary$Pair.ID <- NA
egg_summary$Unpaired.date<-NA
egg_summary$Last.alive.date<-NA
egg_summary$Hatch.fail.date<-NA

#Move all the relevant data from eggs_all to summary dataframe
for(l in 1:length(eggs_all)){
  # Identify the number of rows, egg ID for this egg (i.e. egg #i), and index for egg summary.
  nrows <- dim(eggs_all[[l]])[[1]]
  egg <- eggs_all[[l]]$Egg.ID[1] # egg ID will be the same for all rows so can just take the first one
  index <- which(match(egg_summary$Egg.ID, egg) == 1) # which row in our summary dataframe is for the name of the eggID that we are dealing with here in iteration 'i'
  
  Facility <- eggs_all[[l]]$Facility[1] # Facility (an all variables below) will be the same for all rows so can just take the first one
  egg_summary$Facility[index] <- as.character(Facility) 
  
  # Extract the fate and add it to the summary dataframe
  Fate <- eggs_all[[l]]$Fate[1] 
  egg_summary$Fate[index] <- as.character(Fate)
  
  # Extract the year and add it to the summary dataframe
  Year <- eggs_all[[l]]$Year[1] 
  egg_summary$Year[index] <- as.character(Year)
  
  # Extract the treatment and add it to the summary dataframe
  Treatment <- eggs_all[[l]]$Treatment[1] 
  egg_summary$Treatment[index] <- as.character(Treatment)
  
  # Extract the logger ID and add it to the summary dataframe
  Logger.ID <- eggs_all[[l]]$Logger.ID[1] 
  egg_summary$Logger.ID[index] <- as.character(Logger.ID)
  
  # Extract the pair ID and add it to the summary dataframe
  Pair.ID <- eggs_all[[l]]$Pair.ID[1] 
  egg_summary$Pair.ID[index] <- as.character(Pair.ID)
  
  # Extract the lay date and add it to the summary dataframe
  Lay.date <- eggs_all[[l]]$Lay.date[1] 
  egg_summary$Lay.date[index] <- as.character(Lay.date)
  
  # Extract the unpaired date and add it to the summary dataframe
  Unpaired.date <- eggs_all[[l]]$Unpaired.date[1] 
  egg_summary$Unpaired.date[index] <- as.character(Unpaired.date)
  
  # Extract the last alive date (only for unhatched eggs) and add it to the summary dataframe
  Last.alive.date <- eggs_all[[l]]$Last.alive.date[1] 
  egg_summary$Last.alive.date[index] <- as.character(Last.alive.date)
  
  # Extract the hatch/fail date and add it to the summary dataframe
  Hatch.fail.date <- eggs_all[[l]]$Hatch.date[1] 
  egg_summary$Hatch.fail.date[index] <- as.character(Hatch.fail.date)
  
  print(l)
}

#Add paired date from another summary datasheet-not included in individual csvs
egg_summary_predictors<-merge(egg_summary, Paired.date, by="Egg.ID")

#Remove white space from SHC-SHC group
egg_summary_predictors$Treatment<-gsub('\\s+', '', egg_summary_predictors$Treatment)

#write csv
write.csv(egg_summary_predictors, file = "Egg_summary_survival_predictors.csv", row.names = F)


