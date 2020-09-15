#-----------------Set working directory--------------------
rm(list = ls())

#setwd("C:/Users/HannahE/OneDrive - The Calgary Zoological Society/My Documents/Incubation Study/")

#----------------Load libraries.-------------
require(reshape2)
library(dplyr)
library(here)

#read in the pairing, hatch and non-environmental predictor data
eggs.new <- read.csv(here("data","eggs.new.csv"))

eggs.new$Facility <- as.numeric(as.factor(eggs.new$Facility))
eggs.new$Year <- as.numeric(as.factor(eggs.new$Year))
eggs.new$Treatment2<-rep(0,length(eggs.new$Treatment))#Coded so level one has the largest sample size
eggs.new$Treatment2[which(eggs.new$Treatment=="B")] <- 5
eggs.new$Treatment2[which(eggs.new$Treatment=="P")] <- 2
eggs.new$Treatment2[which(eggs.new$Treatment=="GQF1")] <- 3
eggs.new$Treatment2[which(eggs.new$Treatment=="GQF2")] <- 4
eggs.new$Treatment2[which(eggs.new$Treatment=="SHC")] <- 1
eggs.new$Treatment2[which(eggs.new$Treatment=="WC")] <- 6
eggs.new$Egg.treatment <- as.integer(as.factor(eggs.new$Egg.treatment))#random effect hyp 1-egg type nested in incubation treatment
eggs.new$Egg.type <- as.integer(as.factor(eggs.new$Egg.type))#random effect hyp 2

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

#### -> here you need to develop matrices of the predictors
#### -> these will need to be the same dimensions as enc.hist (rows = eggs, cols = days)
#### -> if the data are missing, put in NAs 
#### -> these matrices will start on the day of placement with logger egg and end hatch/fail date
#read in the daily environmental data
eggs.predict <- read.csv(here("data","Egg_summary_survival_dailydata.csv"))

#Create a variable that accounts for the first day of data being dropped as full day wasnt collected
eggs.new$Paired.date_1<-eggs.new$Paired.date+1
eggs.new$Paired.date_relative<-eggs.new$Paired.date_1-eggs.new$Paired.date

nam <- unique(eggs.predict$Egg.ID)
tres <- vector("list", length(nam))
names(tres) <- nam

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

#day.1 = rep(0, 73), day.2 = rep(0, 73))

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

St.temp <- df_temp_mean
for(i in 1:nrow(df_temp_mean)){
  St.temp[i,1] <- mean(df_temp_mean[i,],na.rm=T)
}
for(i in 1:nrow(df_temp_mean)){
  for(j in 2:ncol(df_temp_mean)){
    if(is.na(St.temp[i,j])){ 
      St.temp[i,j] <- St.temp[i,j-1]
    }
  }
}
St.temp[!is.na(St.temp)] <- NA
for(i in 1:nrow(df_temp_mean)){
  St.temp[i,1:(first[i]-1)] <- NA
  #St.temp[i,(last[i]+1):ncol(df_temp_mean)] <- NA
}

################################DATASET FOR HYP 1 (temp mean)##########################
dataset <- list(nind=nind,first=first,last=last,enc.hist=enc.hist,pair.ID=eggs.new$Pair.ID,
                npair=length(unique(eggs.new$Pair.ID)), ntreat=length(unique(eggs.new$Egg.treatment)),
                facility=eggs.new$Facility, egg.treatment=eggs.new$Egg.treatment,temp.mn=df_temp_mean)

#initial values
inits <- function(){
  list (beta.fac = c(runif(2),NA))
}

#parameters
parameters <- c("beta.int","beta.fac")

######################### CREATE MODEL FILE HYP 1 #########################
#Daily hatch/fail~environ param + environ param^2 (maybe?) + facility + random = pair ID + egg type/treatment

cat("
model {
    
#####LIKELIHOOD
for(i in 1:nind){  
  for(t in (first[i]+1):last[i]){
    
    enc.hist[i,t] ~ dbern(S[i,t-1]*enc.hist[i,t-1])   
    
    logit(S[i,t-1]) <- beta.int + alpha.Pr[pair.ID[i]] + alpha.Treat[egg.treatment[i]]
    
    + beta.fac[facility[i]] + beta.temp.mn*temp.mn[i,t-1] 

    temp.mn[i,t] ~ dnorm(mn.temp.mn[i,t],tau.temp.mn[egg.treatment[i]])

    mn.temp.mn[i,t] <- mu.mn.temp + rho.mn.temp[egg.treatment[i]]*(temp.mn[i,t-1]-mu.mn.temp)

  }
}
    
#PARAMETER PRIORS 
beta.int ~ dnorm(0,0.001)
beta.temp.mn ~ dnorm(0,0.001)

for(f in 1:8){
  tau.temp.mn[f] <- pow(sigma.temp.mn[f],-2)
  sigma.temp.mn[f] ~ dunif(0,5)
  rho.mn.temp[f] ~ dnorm(0,0.001)
}

mu.mn.temp ~ dnorm(0,0.001)


for(f in 1:2){
  beta.fac[f] ~ dnorm(0,0.001)
}
beta.fac[3] <- 0

#RANDOM EFFECT
for(i in 1:npair){
  alpha.Pr[i] ~ dnorm(0,tau.pair)
}
tau.pair <- pow(sigma.pair,-2)
sigma.pair ~ dunif(0,25) 

for(i in 1:ntreat){
  alpha.Treat[i] ~ dnorm(0,tau.treat)
}
tau.treat <- pow(sigma.treat,-2)
sigma.treat ~ dunif(0,25) 

}
",file="nestmodel.txt")

######################### RUN BUGS #########################

start <- Sys.time()

nc <- 3
nb <- 1000
nt <- 1
ni <- 10000

library("jagsUI")
jagsfit.1 <- jags(data=dataset, inits=inits, parameters.to.save=parameters, n.chains=nc, n.burnin = nb, n.iter=ni, n.thin=nt, model.file="nestmodel.txt", parallel=TRUE)

end <- Sys.time() - start
end

jagsfit.1$summary

##########################DATASET FOR HYP 2#########################################
dataset2 <- list(nind=nind,first=first,last=last,enc.hist=enc.hist,pair.ID=eggs.new$Pair.ID,
                npair=length(unique(eggs.new$Pair.ID)), ntype=length(unique(eggs.new$Egg.type)),
                facility=eggs.new$Facility, egg.type=eggs.new$Egg.type, treatment=eggs.new$Treatment2)

#initial values
inits2 <- function(){
  list (beta.fac = c(runif(2),NA), beta.treat=c(runif(5), NA))
}

#parameters
parameters2 <- c("beta.int","beta.fac", "beta.treat")

#for treatment predictions
treatment.range <- seq(1,6, by = 1)

######################### CREATE MODEL FILE HYP 2 #########################
#Daily hatch/fail~treatment + facility + random = pair ID + egg type

cat("
model {
    
#####LIKELIHOOD
for(i in 1:nind){  
  for(t in first[i]:(last[i]-1)){
    
    enc.hist[i,(t+1)] ~ dbern(S[i,t]*enc.hist[i,t])   
    
    logit(S[i,t]) <- beta.int + alpha.Pr[pair.ID[i]] + alpha.Type[egg.type[i]]
    
    +beta.fac[facility[i]]  + beta.treat[treatment[i]]

  }
}
    
#PARAMETER PRIORS 
beta.int ~ dnorm(0,0.001)

for(f in 1:2){
  beta.fac[f] ~ dnorm(0,0.001)
}
beta.fac[3] <- 0

for(j in 2:6) {
    beta.treat[j] ~ dnorm(0, 0.001)
    }
beta.treat[1] <- 0
    
#RANDOM EFFECT
for(i in 1:npair){
  alpha.Pr[i] ~ dnorm(0,tau.pair)
}
tau.pair <- pow(sigma.pair,-2)
sigma.pair ~ dunif(0,25) 

for(i in 1:ntype){
  alpha.Type[i] ~ dnorm(0,tau.type)
}
tau.type <- pow(sigma.type,-2)
sigma.type ~ dunif(0,25) 

#TREATMENT PREDICTIONS
 for(i in 1:length(treatment.range)) {
  logit(pred.S[i]) <- 
  b1[i] + #treatment-specific
  b2[1] + #Predict for Calgary (for cost-benefit comparison)
  b0 +           # common intercept
  0                   # random effect = 0
    } # loop units

}
",file="nestmodel2.txt")

######################### RUN BUGS #########################

start <- Sys.time()

nc <- 3
nb <- 1000
nt <- 1
ni <- 10000

jagsfit.2 <- jags(data=dataset2, inits=inits2, parameters.to.save=parameters2, n.chains=nc, n.burnin = nb, n.iter=ni, n.thin=nt, model.file="nestmodel.txt", parallel=TRUE)

end <- Sys.time() - start
end

jagsfit.2$summary
