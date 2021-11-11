# Script to read and combine data of egg survival and
# data from the data loggers placed at the nest.
# Use the .Rproj default working directory to run the script

rm(list = ls())

#----------------Load libraries.-------------
require(reshape2)
library(dplyr)

if(!dir.exists("data/derived_data/")) dir.create("data/derived_data/")

#read in the pairing, hatch and non-environmental predictor data
eggs.new <- read.csv("data/eggs.new.csv")
unique(cbind(eggs.new$Facility, as.numeric(as.factor(eggs.new$Facility))))
eggs.new$Facility <- as.numeric(as.factor(eggs.new$Facility))
eggs.new$Year <- as.numeric(as.factor(eggs.new$Year))
#Treatment is coded so level one has the largest sample size
eggs.new$Treatment2<-rep(0,length(eggs.new$Treatment))

eggs.new$Treatment2[which(eggs.new$Treatment=="GQF1")] <- 1
eggs.new$Treatment2[which(eggs.new$Treatment=="GQF2")] <- 2
eggs.new$Treatment2[which(eggs.new$Treatment=="P")] <- 3
eggs.new$Treatment2[which(eggs.new$Treatment=="SHC")] <- 4
eggs.new$Treatment2[which(eggs.new$Treatment=="B")] <- 5
eggs.new$Treatment2[which(eggs.new$Treatment=="WC")] <- 6
#random effect hyp 1: egg type nested in incubation treatment
eggs.new$Egg.treatment <- as.integer(as.factor(eggs.new$Egg.treatment))
#random effect hyp 2
eggs.new$Egg.type <- as.integer(as.factor(eggs.new$Egg.type))

#start all nests at time 1  
eggs.new$start <- eggs.new$Paired.date - (eggs.new$Paired.date) + 1
#calc time paired
eggs.new$unpaired <- eggs.new$Unpaired.date - (eggs.new$Paired.date) + 1
#calc time between last alive and paired
eggs.new$last.alive <- eggs.new$Last.alive.date - (eggs.new$Paired.date) + 1
#calc time between paired date and hatch/fail date
eggs.new$term <- eggs.new$Hatch.fail.date - (eggs.new$Paired.date) + 1

#build encounter histories 
#make a empty matrix of the number of days between paired and hatch/fail (term) for each egg
enc.hist <- matrix(nrow=nrow(eggs.new),ncol=max(as.numeric(eggs.new$term)))
enc.hist[] <- NA
#between the paired date and the hatch fail date for each egg, give a one for each day it was 
#alive and a NA when we're unsure if it was alive
#for some failed Patuxent and one failed ICF egg this is the whole time it was paired
#as the last time it was noted alive was the time it was paired with the logger
first <- last <- rep(NA,nrow(eggs.new))
for(i in 1:nrow(enc.hist)){
  enc.hist[i,((eggs.new$start[i]):(eggs.new$last.alive[i]))] <- 1    
  if(eggs.new$outcome[i] == 0){enc.hist[i,(eggs.new$term[i])] <- 0}
  first[i] <- as.integer(eggs.new$start[i])
  last[i] <- as.integer(eggs.new$term[i])
}
nind <- nrow(enc.hist)#No. of individuals

#### -> matrices of the predictors
#### -> the same dimensions as enc.hist (rows = eggs, cols = days)
#### -> NAs for missing predictor data  
#### -> matrices start on the day of placement with logger egg and end on hatch/fail date
#read in the daily environmental data
eggs.predict <- read.csv("data/Egg_summary_survival_dailydata_40min.csv")

#Create a variable that accounts for the first day of data being dropped as full day wasnt collected
eggs.new$Paired.date_1<-eggs.new$Paired.date+1
eggs.new$Paired.date_relative<-eggs.new$Paired.date_1-eggs.new$Paired.date

nam <- unique(eggs.predict$Egg.ID)
tres <- vector("list", length(nam))
names(tres) <- nam

vart <- eggs.predict$rot_mean
eggs.predict$rot_mean <- round((vart - mean(vart, na.rm=T))/ sd(vart, na.rm = T), 2)

vart <- eggs.predict$rot_var
eggs.predict$rot_var <- round((vart - mean(vart, na.rm=T))/ sd(vart, na.rm = T), 2)

vart <- eggs.predict$rh_mean
  eggs.predict$rh_mean <- round((vart - mean(vart, na.rm=T))/ sd(vart, na.rm = T), 2)

vart <- eggs.predict$rh_var
  eggs.predict$rh_var <- round((vart - mean(vart, na.rm=T))/ sd(vart, na.rm = T), 2)

vart <- eggs.predict$temp_mean
  eggs.predict$temp_mean <- round((vart - mean(vart, na.rm=T))/ sd(vart, na.rm = T), 2)

vart <- eggs.predict$temp_var
  eggs.predict$temp_var <-round( (vart - mean(vart, na.rm=T))/ sd(vart, na.rm = T), 2)


df_rot_mean <- df_rot_var <-df_temp_mean<-df_temp_var<-df_rh_mean<-df_rh_var<- data.frame(matrix(NA, nrow=73, ncol= 33))
df_rot_mean$Egg.ID <- df_rot_var$Egg.ID<-df_temp_mean$Egg.ID<-df_temp_var$Egg.ID<-df_rh_mean$Egg.ID<-df_rh_var$Egg.ID <- nam

for(i in 1:length(nam)){
  pred <- subset(eggs.predict, Egg.ID == nam[i], 
                 select = c(Egg.ID, Day, rot_mean, rot_var, temp_mean, temp_var, rh_mean, rh_var))
  info <- data.frame(Egg.ID = rep(nam[i], 33), Paired.date = rep(0,33), Paired.date_relative =1:33)
  
  real <- subset(eggs.new, Egg.ID ==nam[i],
                 select =c(Egg.ID, Paired.date, Paired.date_relative))
  info[real$Paired.date_relative:33, "Paired.date"] <- real$Paired.date:(real$Paired.date + 33 - real$Paired.date_relative)
  
  temp <- merge(pred, info, by.x = "Day", by.y = "Paired.date", all.y = T)
  
  if(length(df_rot_mean[df_rot_mean$Egg.ID == nam[i], 1:33]) != length(t(temp$rot_mean))) 
    print(paste("id= ", nam[i], " and i= ", i, " lengths are different"))
  
  df_rot_mean[df_rot_mean$Egg.ID == nam[i], 1:33] <- t(temp$rot_mean)
  df_rot_var[df_rot_var$Egg.ID == nam[i], 1:33]   <- t(temp$rot_var)
  df_temp_mean[df_temp_mean$Egg.ID == nam[i], 1:33]   <- t(temp$temp_mean)
  df_temp_var[df_temp_var$Egg.ID == nam[i], 1:33]   <- t(temp$temp_var)
  df_rh_mean[df_rh_mean$Egg.ID == nam[i], 1:33]   <- t(temp$rh_mean)
  df_rh_var[df_rh_var$Egg.ID == nam[i], 1:33]   <- t(temp$rh_var)
} 

#Drop ID col and convert to matrix
dfs <- list(df_rh_mean, df_rh_var, df_rot_mean, df_rot_var, df_temp_mean, df_temp_var)
dfs<-lapply(dfs, function(x) x[-34])
df_rh_mean <- matrix(unlist(dfs[1]), ncol = 33)
df_rh_var <- matrix(unlist(dfs[2]), ncol = 33) 
df_rot_mean <- matrix(unlist(dfs[3]), ncol = 33) 
df_rot_var <- matrix(unlist(dfs[4]), ncol = 33) 
df_temp_mean <- matrix(unlist(dfs[5]), ncol = 33) 
df_temp_var <- matrix(unlist(dfs[6]), ncol = 33) 

#Code a simpler egg ID for both files 
eggs.predict$Egg <- as.numeric(as.factor(eggs.predict$Egg.ID)) 
eggs.new$Egg <- as.numeric(as.factor(eggs.new$Egg.ID)) 

#this object binds together encounter histories, the predictor data (temp mean in this case) and first/last 
#for the purposes of examining the data in a useful way  
#dimensions are row (enc.hist), col(enc.hist+1) - to allow for inclusion of first/last, and 2 
#print all.data.array[x,,] for any x to look at a summary of data for the full nest 
#the first column will be enc.hist, the second column will be temp_mean, and the last row will be first/last
all.data.array <- array(NA,dim=c((dim(enc.hist)[1]), (dim(enc.hist)[2]+1), 2))
all.data.array[,1:dim(enc.hist)[2],1] <- enc.hist
all.data.array[,1:dim(enc.hist)[2],2] <- df_temp_mean
all.data.array[,dim(enc.hist)[2]+1,1] <- first
all.data.array[,dim(enc.hist)[2]+1,2] <- last
            
write.csv(enc.hist, file="data/enc.hist.csv", row.names=F)
write.csv(df_temp_mean, file="data/df_temp_mean.csv", row.names=F)
write.csv(df_temp_var, file="data/df_temp_var.csv", row.names=F)
write.csv(df_rh_mean, file="data/df_rh_mean.csv", row.names=F)
write.csv(df_rh_var, file="data/df_rh_var.csv", row.names=F)
write.csv(df_rot_mean, file="data/df_rot_mean.csv", row.names=F)
write.csv(df_rot_var, file="data/df_rot_var.csv", row.names=F)
write.csv(first, file="data/first.csv", row.names=F)
write.csv(last, file="data/last.csv", row.names=F)


###### create dummy variable for contrasts

eggs.new$B <- ifelse(eggs.new$Treatment=="B", 1, 0)
eggs.new$GQF1 <- ifelse(eggs.new$Treatment=="GQF1", 1, 0)
eggs.new$GQF2 <- ifelse(eggs.new$Treatment=="GQF2", 1, 0)
eggs.new$P <- ifelse(eggs.new$Treatment=="P", 1, 0)
eggs.new$SHC <- ifelse(eggs.new$Treatment=="SHC", 1, 0)
eggs.new$WC <- ifelse(eggs.new$Treatment=="WC", 1, 0)
eggs.new$CZ <- ifelse(eggs.new$Facility==1, 1, 0)
eggs.new$ICF <- ifelse(eggs.new$Facility==2, 1, 0)     
eggs.new$Patuxent <- ifelse(eggs.new$Facility==3, 1, 0)
eggs.new$notCZ <- ifelse(eggs.new$Facility!=1, 1, 0)

last.enc <- NULL
for(i in 1:73){
  last.enc[i] <- enc.hist[i, last[i]]
}

################################DATASET FOR HYP 1 (temp mean)##########################
dataset <- list(nind=nind,first=first,last=last,enc.hist=enc.hist,pair.ID=eggs.new$Pair.ID,
                npair=length(unique(eggs.new$Pair.ID)), ntreat=length(unique(eggs.new$Egg.treatment)),
                facility=eggs.new$Facility, egg.treatment=eggs.new$Egg.treatment,temp.mn=df_temp_mean,
                temp.var = df_temp_var, rot.mn = df_rot_mean, rot.var = df_rot_var,
                rh.mn = df_rh_mean, rh.var = df_rh_var,
                cumsumm = c(0, cumsum(last-1))
                )
save(dataset, file = "data/derived_data/dataset_hyp_1_40min.RData")
##########################DATASET FOR HYP 2#########################################
dataset2 <- list(nind=nind,first=first,last=last,enc.hist=enc.hist,pair.ID=eggs.new$Pair.ID,
                 npair=length(unique(eggs.new$Pair.ID)), ntype=length(unique(eggs.new$Egg.type)),
                 facility=eggs.new$Facility, species=eggs.new$Egg.type, treatment=eggs.new$Treatment2,
                 treatment.range=length(unique(eggs.new$Egg.treatment)),
                 cumsumm = c(0, cumsum(last-1)),
                 last.enc = last.enc,
                 ntreat_by_species=length(unique(eggs.new$Egg.treatment)),
                 treat_by_species=eggs.new$Egg.treatment, 
                 treat_by_species1 = ifelse(eggs.new$Egg.treatment ==1, 1, 0),
                 treat_by_species2 = ifelse(eggs.new$Egg.treatment ==2, 1, 0),
                 treat_by_species3 = ifelse(eggs.new$Egg.treatment ==3, 1, 0),
                 treat_by_species4 = ifelse(eggs.new$Egg.treatment ==4, 1, 0),
                 treat_by_species5 = ifelse(eggs.new$Egg.treatment ==5, 1, 0),
                 treat_by_species6 = ifelse(eggs.new$Egg.treatment ==6, 1, 0),
                 treat_by_species7 = ifelse(eggs.new$Egg.treatment ==7, 1, 0),
                 treat_by_species8 = ifelse(eggs.new$Egg.treatment ==8, 1, 0),
                 
                 treatment_B = eggs.new$B,
                 treatment_GQF1 = eggs.new$GQF1, treatment_GQF2 = eggs.new$GQF2, 
                 treatment_P = eggs.new$P, treatment_SHC = eggs.new$SHC, treatment_WC = eggs.new$WC,         
                 facility_CZ = eggs.new$CZ, facility_notCZ =  eggs.new$notCZ, facility_ICF = eggs.new$ICF, facility_Patuxent = eggs.new$Patuxent
)
dataset3 <- list(nind=nind,first=first,last=last,enc.hist=enc.hist,pair.ID=eggs.new$Pair.ID,
                 npair=length(unique(eggs.new$Pair.ID)), ntype=length(unique(eggs.new$Egg.type)),
                 facility=eggs.new$Facility, species=eggs.new$Egg.type, treatment=eggs.new$Treatment2,
                 treatment_factor=eggs.new$Treatment,
                 treatment.range=length(unique(eggs.new$Egg.treatment)),
                 cumsumm = c(0, cumsum(last-1)),
                             treatment_B = eggs.new$B,
                             treatment_GQF1 = eggs.new$GQF1, treatment_GQF2 = eggs.new$GQF2, 
                             treatment_P = eggs.new$P, treatment_SHC = eggs.new$SHC, treatment_WC = eggs.new$WC,         
facility_CZ = eggs.new$CZ, facility_notCZ =  eggs.new$notCZ, facility_ICF = eggs.new$ICF, facility_Patuxent = eggs.new$Patuxent
)
save(dataset2, file = "data/derived_data/dataset_hyp_2_40min.RData")
save(dataset3, file = "data/derived_data/dataset3_hyp_2_40min.RData")
