#-----------------Set working directory--------------------
rm(list = ls())

#setwd("/Users/hannahedwards/Documents/Calgary/CZ/Incubation Study")
#setwd("C:/Users/HannahE/OneDrive - The Calgary Zoological Society/My Documents/Incubation Study/")

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
loc_path <- "data/daily_data"
all_f <- list.files(path=loc_path, full.names=T)
for(i in 1:length(all_f)){
  assign(paste0("egg", i), read.csv(all_f[i], header = T))
}
egg_IDs<-read.csv("data/egg_IDs.csv")

# Combine all the data into one list for analysis.
#######
#issue! egg2 has no date recoreded
#egg2$Date <- rep("19/04/2017", dim(egg2)[1])
eggs_all <- list(egg1, egg2, egg3)

#add a column for each egg that creates a time stamp
eggs_all<- lapply(1:length(eggs_all), function(x){
  nrows <- dim(eggs_all[[x]])[[1]]
  Min2<-str_pad(eggs_all[[x]]$Min, 2, pad = "0")#Add zero to single digits
  eggs_all[[x]]$Time <- paste(eggs_all[[x]]$Hour,":",Min2,":00",sep="")
  eggs_all[[x]]$Date <- as.character(eggs_all[[x]]$Date)
  eggs_all[[x]]$Dateformat <- as.Date(eggs_all[[x]]$Date, format="%d/%m/%Y")
  timestamp<- paste(eggs_all[[x]]$Dateformat, eggs_all[[x]]$Time, sep=" ")
  eggs_all[[x]]$timestamp <- as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S") 
  return(eggs_all[[x]])
  
})

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
  if(length(which(tmp$time_in_minutes %in% lags))== length(lags)){
    tmp[which(tmp$time_in_minutes %in% lags), "lags_24_hours"] <- lags
  } else {
    tmp[which(tmp$time_in_minutes %in% lags), "lags_24_hours"] <- lags[which(lags %in% tmp$time_in_minutes)]
  }
  # If the last day does not cover the 24 hours delete it
  last_day <- dim(tmp)[1] - which(tmp$lags == max(lags))
  if(last_day < 1440){
    tmp <- tmp[-c(which(tmp$lags == max(lags)):dim(tmp)[1]), ]
  }
  
  tmp <- tmp %>% fill(lags_24_hours, .direction="down")
  tmp <- tmp[tmp$time_in_minutes %in% lags_40min, ]
  Day <- tmp[1, "Day"]
  if(is.null(Day)){
    tmp[, "new_day"] <- rep(NA, dim(tmp)[1])
    warning(paste0("ID ", unique(tmp$Egg.ID), " has no recorded Day. Is this true or do we have a typo in column name?"))
  } else{
    tmp[, "new_day"] <- Day:(Day+dim(tmp)[1]-1)
  }
  egg_data1[[i]] <- tmp
} 

egg_data2 = list()
for(j in 1:length(egg_data1)){
  egg_summs <- egg_data1[[j]]
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
  #Drop the first observation of the day in position change col
  #egg_summs2<- droplevels( egg_summs[-which(egg_summs$Time == "6:40:00"), ] )
  # I am not sure what this is about! therefore I leave as is
  egg_summs2 <- egg_summs
  #Calc the mean and var per day for each environ variable.
  if(any(names(egg_summs)=="EggID")){
    egg_summs$Egg.ID <- egg_summs$EggID
  }
  Daily_temp_mean<-aggregate(temp_mean~new_day+Egg.ID, data=egg_summs, FUN='mean')
  Daily_temp_var<-aggregate(temp_mean~new_day, data=egg_summs, FUN='var')
  Daily_temp_var<-rename(Daily_temp_var, replace= c("temp_mean"="temp_var"))
  Daily_rh_mean<-aggregate(rh_mean~new_day, data=egg_summs, FUN='mean')
  Daily_rh_var<-aggregate(rh_mean~new_day, data=egg_summs, FUN='var')
  Daily_rh_var<-rename(Daily_rh_var, replace= c("rh_mean"="rh_var"))
  Daily_rot_mean<-aggregate(position_change~new_day, data=egg_summs2, FUN='mean')
  Daily_rot_mean<-rename(Daily_rot_mean, replace= c("position_change"="rot_mean"))
  Daily_rot_var<-aggregate(position_change~new_day, data=egg_summs2, FUN='var')
  Daily_rot_var<-rename(Daily_rot_var, replace= c("position_change"="rot_var"))
  #Merge
  Daily_rot<-merge(Daily_rot_mean, Daily_rot_var, by="new_day")
  egg_data_summs<-merge(Daily_rot, Daily_temp_mean, by="new_day", all=TRUE)
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
  Pair<-read.csv("data/pair_day.csv")
  Paired.date<-subset(Pair, select=c(Egg.ID, Paired.date))
  egg_data2<-merge(egg_data2, Paired.date, by="Egg.ID")
  
  #Add unpaired date column to egg_data2
  egg_data2$Unpaired.date<- NA
  #Move unpaired data from eggs_all to egg_data2
  for(k in 1:length(eggs_all)){
    nrows <- dim(eggs_all[[k]])[[1]]
    egg <- eggs_all[[k]]$Egg.ID[1] # egg ID will be the same for all rows so can just take the first one
    index <- which(match(egg_data2$Egg.ID, egg) == 1) # which row in egg_data2 is for the name of the eggID that we are dealing with here in iteration 'i'
    
    # Extract the unpaired date and add it to the summary dataframe
    if(!is.null(eggs_all[[k]]$Unpaired.date[1])){
      Unpaired.date <- eggs_all[[k]]$Unpaired.date[1] 
      egg_data2$Unpaired.date[index] <- as.character(Unpaired.date)
    } else {
      egg_data2$Unpaired.date[index] <- NA
      warning(paste(" egg ", k, " has no recorded Unpaired.date"))
    }
  }
  
  #Remove the fist and last day the egg was paired-not a full day of data collection-can exclude code
  #if we calc the mean and variance within a 24 hour period from the time it was paired.
  #egg_data2$Unpaired.date<-as.integer(as.character(egg_data2$Unpaired.date))
  #egg_data2$Paired.date<-as.integer(as.character(egg_data2$Paired.date))
  #egg_data2$Day<-as.integer(as.character(egg_data2$Day))
  
  #egg_data2<- egg_data2[!(egg_data2$Day>=egg_data2$Unpaired.date),]
  #egg_data2<- egg_data2[!(egg_data2$Day<=egg_data2$Paired.date),]
  
  write.csv(egg_data2, file = "data/Egg_summary_survival_dailydata_40min.csv", row.names = F)

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
  if(!is.null(Last.alive.date)){
    egg_summary$Last.alive.date[index] <- as.character(Last.alive.date)
  } else {
    warning(paste("egg ", l, " has no Last.alive.date"))
  }
  # Extract the hatch/fail date and add it to the summary dataframe
  Hatch.fail.date <- eggs_all[[l]]$Hatch.date[1] 
  egg_summary$Hatch.fail.date[index] <- as.character(Hatch.fail.date)
  
  print(l)
}

#Add paired date from another summary datasheet-not included in individual csvs
egg_summary_predictors<-merge(egg_summary, Paired.date, by="Egg.ID")

#Remove white space from SHC-SHC group
egg_summary_predictors$Treatment<-gsub('\\s+', '', egg_summary_predictors$Treatment)

write.csv(egg_summary_predictors, file = "data/Egg_summary_survival_predictors2.csv", row.names = F)

#S17#1, S17#2,S2#2, Y42#4, Y42#5-added data from a different time point because
#of data gaps when data was offloaded during the study period
