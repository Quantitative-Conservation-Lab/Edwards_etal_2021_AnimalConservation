#This code runs the model for hypothesis three described in Edwards et al.   
#Written by Hannah Edwards

# Scripts to run the models in JAGS with the jagsUI package.
# These models address the hypothesis that incubation treatment affects  
# the temp/humidity/rotation of the egg during incubation.

#----------------Load libraries.-------------
library(dplyr)
library(ggplot2) # for graphing
library(jagsUI) # for the JAGS model
library(loo) # for the JAGS model
library(rjags) # for the JAGS model
library(reshape2) # for prior to graphing with ggplot

#---------------Load the data------------------
data_environ <- read.csv("Egg_summary_survival_dailydata_40min_hyp2.csv")#daily averages and variances
data_predict<- read.csv("eggs.new.csv")

#---------------Code variables-----------------
nind <- length(data_predict$Egg.ID)

tempmean<-data_environ$temp_mean
tempvar<-data_environ$temp_var
rhmean<-data_environ$rh_mean
rhvar<-data_environ$rh_var
rotmean<-data_environ$rot_mean
rotvar<-data_environ$rot_var

treatment<-rep(0,length(data_predict$Treatment))#Coded so that the treatment with the best hatch success is contrasted with all other treatments
treatment[which(data_predict$Treatment=="B")] <- 5
treatment[which(data_predict$Treatment=="P")] <- 2
treatment[which(data_predict$Treatment=="GQF1")] <- 3
treatment[which(data_predict$Treatment=="GQF2")] <- 4
treatment[which(data_predict$Treatment=="SHC")] <- 1
treatment[which(data_predict$Treatment=="WC")] <- 6

treatment2<-rep(0,length(data_predict$Treatment2))
treatment2[which(data_predict$Treatment2=="NI")] <- 2
treatment2[which(data_predict$Treatment2=="AI")] <- 1

pair <- data_predict$Pair.ID
Npair <- length(levels(as.factor(data_predict$Pair.ID)))

egg.type<-rep(0,length(data_predict$Egg.type))
egg.type[which(data_predict$Egg.type=="SHC")] <- 1
egg.type[which(data_predict$Egg.type=="WC")] <- 2 
Neggtype <- length(levels(as.factor(data_predict$Egg.type)))

facility<-rep(0,length(data_predict$Facility))
facility[which(data_predict$Facility=="CZ")] <- 3 
facility[which(data_predict$Facility=="ICF")] <- 2
facility[which(data_predict$Facility=="Patuxent")] <- 1 

#-------------------------Define and run the model------------------------------------

##Random intercepts model (equivalent to tempmean~treatment+facility+egg type)
#Specify normal distribution: http://biometry.github.io/APES//LectureNotes/StatsCafe/Linear_models_jags.html
#Temp mean
treatment.range <- seq(1,6, by = 1)

model.string.1<-"model {

     for(i in 1:nind){
     #Likelihood function
     tempmean[i] ~ dnorm(mu[i], tau)
     mu[i] <- b0 + b1[treatment[i]] + b2[facility[i]] + b3[egg.type[i]]
     }

     #First level priors
     b0 ~ dnorm(0, 0.01)
     b1[1] <- 0
     for(i in 2:6){
     b1[i] ~ dnorm(0,0.01)
     }
     b2[1] <- 0
     for(i in 2:3){
     b2[i] ~ dnorm(0,0.01)
     }

     #Second level priors
     tau <- 1/(sigma*sigma)
     sigma ~ dunif(0,100)

    b3[1] <- 0 
    b3[2] ~ dnorm(0,0.01)

   }
" cat(file = "tempmean_treatment.jags ", model.string.1)

##Temp var
model.string.2<- "model { 
    
    for(i in 1:nind){
    
    #Likelihood function
    
    tempvar[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1[treatment[i]] + b2[facility[i]] + b3[egg.type[i]]
    
    } 
    
    #First level priors
    b0 ~ dnorm(0, 0.01)
    b1[1] <- 0
    for(i in 2:6){
    b1[i] ~ dnorm(0,0.01)
    }
    b2[1] <- 0
    for(i in 2:3){
    b2[i] ~ dnorm(0,0.01)
    }

    #Second level priors
    tau <- 1/(sigma*sigma)
    sigma ~ dunif(0,100)

    b3[1] <- 0 
    b3[2] ~ dnorm(0,0.01)
    
    }
    "

cat(file = "tempvar_treatment.jags ", model.string.2)

tempvar_treatment.jags <-"tempvar_treatment.jags "

#RH mean
model.string.3<- "model { 
    
    for(i in 1:nind){
    
    #Likelihood function
    
    rhmean[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1[treatment[i]] + b2[facility[i]] + b3[egg.type[i]]
    
    } 
    
    #First level priors
    b0 ~ dnorm(0, 0.01)
    b1[1] <- 0
    for(i in 2:6){
    b1[i] ~ dnorm(0,0.01)
    }
    b2[1] <- 0
    for(i in 2:3){
    b2[i] ~ dnorm(0,0.01)
    }

    #Second level priors
    tau <- 1/(sigma*sigma)
    sigma ~ dunif(0,100)

    b3[1] <- 0 
    b3[2] ~ dnorm(0,0.01)
    
    }
    "
cat(file = "rhmean_treatment.jags ", model.string.3)

rhmean_treatment.jags <-"rhmean_treatment.jags "

##Rh var
model.string.4<-"model { 
    
    for(i in 1:nind){
    
    #Likelihood function
    
    rhvar[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1[treatment[i]] + b2[facility[i]] + b3[egg.type[i]]
    
    } 
    
    #First level priors
    b0 ~ dnorm(0, 0.01)
    b1[1] <- 0
    for(i in 2:6){
    b1[i] ~ dnorm(0,0.01)
    }
    b2[1] <- 0
    for(i in 2:3){
    b2[i] ~ dnorm(0,0.01)
    }

    #Second level priors
    tau <- 1/(sigma*sigma)
    sigma ~ dunif(0,100)

    b3[1] <- 0 
    b3[2] ~ dnorm(0,0.01)
    
    }
    "
cat(file = "rhvar_treatment.jags ", model.string.4)

rhvar_treatment.jags <-"rhvar_treatment.jags "

#Rot mean
model.string.5<-"model { 
    
    for(i in 1:nind){
    
    #Likelihood function
    
    rotmean[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1[treatment[i]] + b2[facility[i]] + b3[egg.type[i]]
    
    } 
    
    #First level priors
    b0 ~ dnorm(0, 0.01)
    b1[1] <- 0
    for(i in 2:6){
    b1[i] ~ dnorm(0,0.01)
    }
    b2[1] <- 0
    for(i in 2:3){
    b2[i] ~ dnorm(0,0.01)
    }

    #Second level priors
    tau <- 1/(sigma*sigma)
    sigma ~ dunif(0,100)

    b3[1] <- 0 
    b3[2] ~ dnorm(0,0.01)
    
    }
    "
cat(file = "rotmean_treatment.jags ", model.string.5)

rotmean_treatment.jags <-"rotmean_treatment.jags "

#Rot var
model.string.6<-"model { 
    
    for(i in 1:nind){
    
    #Likelihood function
    
    rotvar[i] ~ dnorm(mu[i], tau)
    
    mu[i] <- b0 + b1[treatment[i]] + b2[facility[i]] + b3[egg.type[i]]
    
    } 
    
    #First level priors
    b0 ~ dnorm(0, 0.01)
    b1[1] <- 0
    for(i in 2:6){
    b1[i] ~ dnorm(0,0.01)
    }
    b2[1] <- 0
    for(i in 2:3){
    b2[i] ~ dnorm(0,0.01)
    }

    #Second level priors
    tau <- 1/(sigma*sigma)
    sigma ~ dunif(0,100)

    b3[1] <- 0 
    b3[2] ~ dnorm(0,0.01)
     
    }
    "
cat(file = "rotvar_treatment.jags ", model.string.6)

rotvar_treatment.jags <-"rotvar_treatment.jags "

#Complile JAGS data
jags.data <- list(tempmean=tempmean,nind=nind,treatment=treatment,facility=facility, egg.type=egg.type,Neggtype = Neggtype, treatment.range=treatment.range)
jags.data_eggfixed <- list(tempmean=tempmean,nind=nind,treatment=treatment,facility=facility, egg.type=egg.type,Neggtype = Neggtype, treatment.range=treatment.range, egg.range=egg.range)
jags.data2 <- list(tempvar=tempvar,nind=nind,treatment=treatment,facility=facility, egg.type=egg.type,Neggtype = Neggtype, treatment.range=treatment.range)
jags.data3 <- list(rhmean=rhmean,nind=nind,treatment=treatment,facility=facility, egg.type=egg.type,Neggtype = Neggtype, treatment.range=treatment.range)
jags.data4 <- list(rhvar=rhvar,nind=nind,treatment=treatment,facility=facility, egg.type=egg.type,Neggtype = Neggtype, treatment.range=treatment.range)
jags.data5 <- list(rotmean=rotmean,nind=nind,treatment=treatment,facility=facility, egg.type=egg.type,Neggtype = Neggtype, treatment.range=treatment.range)
jags.data6 <- list(rotvar=rotvar,nind=nind,treatment=treatment,facility=facility, egg.type=egg.type,Neggtype = Neggtype, treatment.range=treatment.range)

# Initial values 
inits <- function(){list(b0=runif(1), b1=c(NA, runif(5)), b2=c(NA, runif(2)), b3=c(NA, runif(1)))} 
inits_eggfixed <- function(){list(b0=runif(1), b1=c(NA, runif(5)), b2=c(NA, runif(2)), b3=c(NA, runif(1)))} 

# Parameters monitored
parameters <- c("b0","b1","b2","b3")

# MCMC settings
ni <- 300000; nt <- 5; nb <- ni-80000; nc <- 4

# Call JAGS from R (jagsUI), check convergence and summarize posteriors
egg.model <- jags(jags.data, inits, parameters, tempmean_treatment.jags, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
egg.model_eggfixed <- jags(jags.data_eggfixed, inits_eggfixed, parameters, tempmean_treatment_eggfixed.jags, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
egg.model2 <- jags(jags.data2, inits, parameters, tempvar_treatment.jags, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
egg.model3 <- jags(jags.data3, inits, parameters, rhmean_treatment.jags, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
egg.model4 <- jags(jags.data4, inits, parameters, rhvar_treatment.jags, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
egg.model5 <- jags(jags.data5, inits, parameters, rotmean_treatment.jags, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)
egg.model6 <- jags(jags.data6, inits, parameters, rotvar_treatment.jags, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE)

# See https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/, for WAIC comparison
save(egg.model, file ="Model output/Hyp3_Tempmean.Rdata")
load(file ="Model output/Hyp3_Tempmean.Rdata")
print(egg.model)

save(egg.model2, file ="Model output/Hyp3_Tempvar.Rdata")
load(file ="Model output/Hyp3_Tempvar.Rdata")
print(egg.model2)

save(egg.model3, file ="Model output/Hyp3_RHmean.Rdata")
load(file ="Model output/Hyp3_RHmean.Rdata")
print(egg.model3)

save(egg.model4, file ="Model output/Hyp3_RHvar.Rdata")
load(file ="Model output/Hyp3_RHvar.Rdata")
print(egg.model4)

save(egg.model5, file ="Model output/Hyp3_Rotmean.Rdata")
load(file ="Model output/Hyp3_Rotmean.Rdata")
print(egg.model5)

save(egg.model6, file ="Model output/Hyp3_Rotvar.Rdata")
load(file ="Model output/Hyp3_Rotvar.Rdata")
print(egg.model6)

#Traceplots model 1
op <- par(mfrow = c(2, 1))
jagsUI::traceplot(egg.model, parameters = c("b0"))
jagsUI::traceplot(egg.model, parameters = c("b1"))
plot(egg.model)#if you want to look at den plots

#model 2
op <- par(mfrow = c(2, 1))
jagsUI::traceplot(egg.model2, parameters = c("b0"))
jagsUI::traceplot(egg.model2, parameters = c("b1"))

#model 3
op <- par(mfrow = c(2, 1))
jagsUI::traceplot(egg.model3, parameters = c("b0"))
jagsUI::traceplot(egg.model3, parameters = c("b1"))

#model 4
op <- par(mfrow = c(2, 1))
jagsUI::traceplot(egg.model4, parameters = c("b0"))
jagsUI::traceplot(egg.model4, parameters = c("b1"))

#model 5
op <- par(mfrow = c(2, 1))
jagsUI::traceplot(egg.model5, parameters = c("b0"))
jagsUI::traceplot(egg.model5, parameters = c("b1"))

#model 6
op <- par(mfrow = c(2, 1))
jagsUI::traceplot(egg.model6, parameters = c("b0"))
jagsUI::traceplot(egg.model6, parameters = c("b1"))

###Visualize relationships
#Mean temp
ggplot(data, aes(x=Treatment, y=temp_mean)) + 
  geom_boxplot()

#Variation in temp
ggplot(data, aes(x=Treatment, y=temp_var)) + 
  geom_boxplot()

#Mean RH
ggplot(data, aes(x=Treatment, y=rh_mean)) + 
  geom_boxplot()

#Variation in RH
ggplot(data, aes(x=Treatment, y=rh_var)) + 
  geom_boxplot()

#Rotation
ggplot(data, aes(x=Treatment, y=rot)) + 
  geom_boxplot()
