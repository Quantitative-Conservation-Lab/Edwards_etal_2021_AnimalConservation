#-----------------------------------------------------


# Scripts to run the models in JAGS with the jagsUI package.
# These models address the hypothesis that the environmental predictor rh.var measured in the nest 
# affect the survival of the eggs.
# Models are autoregressive 
# JAGS scripts are saved as txt files within the folder scripts/jags/jags_models/rh.var

# At the end we check convergency of the chains, do basic model diagnostics and
# compare model weights.

# All models are autoregressive and loop through each individual egg.
# The m_only_intercepts model is the NULL hypothesis and includes egg.treatment (species * incubation treatment of the eggs)
# and facility for multilevel intercepts.
# The environmental predictor is included as an exploratory variable and also as a quadratic term.


# Output is "sink" to file "model_outputs/jags/compare_rh.var.txt" 
# Plots are saved in model_outputs/jags/rh.var
#-----------------------------------------------------

#----------------Load libraries.-------------# 
library("jagsUI")
library(coda)
library(gdata)
library(loo)
library(stringr)
library(HDInterval)
library(forestplot)
library(viridis)

if(file.exists("data/derived_data/dataset_hyp_1_40min.RData")){
  load(file = "data/derived_data/dataset_hyp_1_40min.RData")
} else{
  # if data/derived_data/dataset_list.RData is not available first run
  source("scripts/getData_for_jags_40min.R")
  rm(list=ls())
  load(file = "data/derived_data/dataset_hyp_1_40min.RData")
}
if(!dir.exists("model_outputs/")){
  dir.create("model_outputs/")
}
if(!dir.exists("model_outputs/jags/")) {
     dir.create("model_outputs/jags/")
}
if(!dir.exists("model_outputs/jags/jags_models/")) {
  dir.create("model_outputs/jags/jags_models/")
}
if(!dir.exists("model_outputs/jags/mean_hdi/")){
  dir.create("model_outputs/jags/mean_hdi/")
}
if(!dir.exists("scripts/jags/jags_models/rh.var/")){
    dir.create("scripts/jags/jags_models/rh.var/")
}
if(!dir.exists("model_outputs/jags/plots/")){
  dir.create("model_outputs/jags/plots/")
}
test<- FALSE
#test <- TRUE

if(test==TRUE){
  # run short chains for testing
  nc <- 2
  nt <- 1#5 #thinning
  ni <- 500#50000 #iteration
  nb <- 50#ni - 40000 #burn-in
  na <- 500
} else{
  nc <- 4
  nt <- 5 #thinning
  ni <- 300000 #iteration
  nb <- ni - 80000 #burn-in
  na <- 100000 # adaptation
  (ni-nb)/5
  if(nb<=0) stop("burn in is negative")
}


#initial values
inits <- function(){
  list (beta.fac = c(runif(2),NA))
}

#parameters

#-----------------We want to monitor the intercept and all the predictors, and do we want to calculate log-likelihood (and thus need to monitor that)
#-----------------to compare this model to a null model and a model with a quadratic term for the environmental parameter?
parameters_quadratic <- c("log.like.id", "beta.int","beta.fac","beta.rh.var","beta.rh.var2", "mu.mn.rh.var","rho.mn.rh.var","sigma.rh.var",
                         "alpha.Treat", "sigma.treat")
parameters_linear <- c("log.like.id","beta.int","beta.fac","beta.rh.var","mu.mn.rh.var","rho.mn.rh.var","sigma.rh.var",
                        "alpha.Treat", "sigma.treat")
parameters_null <- c("log.like.id", "beta.int","beta.fac", "mu.mn.rh.var","rho.mn.rh.var","sigma.rh.var",
                      "alpha.Treat", "sigma.treat")
######################### CREATE MODEL FILE HYP 1 #########################
# quadratic
cat("
model {
    
#####LIKELIHOOD
for(i in 1:nind){  
  log.like2[i, 1] <- 0
  for(t in (first[i]+1):last[i]){
    
    enc.hist[i,t] ~ dbern(S[i,t-1]*enc.hist[i,t-1])   
    
    #model for daily survival  
    logit(S[i,t-1]) <- alpha.Treat[egg.treatment[i]]
    
    + beta.fac[facility[i]] + beta.rh.var*rh.var[i,t-1] + beta.rh.var2*rh.var[i,t-1]*rh.var[i,t-1]

    rh.var[i,t] ~ dnorm(mn.rh.var[i,t],tau.rh.var[egg.treatment[i]])

    mn.rh.var[i,t] <- mu.mn.rh.var + rho.mn.rh.var*(rh.var[i,t-1]-mu.mn.rh.var)

   #log-likelihood calc
    #log.like[((cumsumm[i] +t)-1)]<-logdensity.bin(enc.hist[i,t], S[i,t-1]*enc.hist[i,t-1] , 1) 
       log.like2[i, t] <- logdensity.bin(enc.hist[i,t], S[i,t-1]*enc.hist[i,t-1] , 1) 
  }
}
log.like2[1, 32] <- 0;log.like2[1, 33] <- 0;log.like2[2, 33] <- 0;log.like2[3, 33] <- 0;log.like2[4, 32] <- 0;log.like2[4, 33] <- 0;log.like2[5, 32] <- 0;log.like2[5, 33] <- 0;log.like2[6, 32] <- 0;log.like2[6, 33] <- 0;log.like2[7, 30] <- 0;log.like2[7, 31] <- 0;log.like2[7, 32] <- 0;log.like2[7, 33] <- 0;log.like2[8, 28] <- 0;log.like2[8, 29] <- 0;log.like2[8, 30] <- 0;log.like2[8, 31] <- 0;log.like2[8, 32] <- 0;log.like2[8, 33] <- 0;log.like2[9, 20] <- 0;log.like2[9, 21] <- 0;log.like2[9, 22] <- 0;log.like2[9, 23] <- 0;log.like2[9, 24] <- 0;log.like2[9, 25] <- 0;log.like2[9, 26] <- 0;log.like2[9, 27] <- 0;log.like2[9, 28] <- 0;log.like2[9, 29] <- 0;log.like2[9, 30] <- 0;log.like2[9, 31] <- 0;log.like2[9, 32] <- 0;log.like2[9, 33] <- 0;log.like2[10, 30] <- 0;log.like2[10, 31] <- 0;log.like2[10, 32] <- 0;log.like2[10, 33] <- 0;log.like2[11, 30] <- 0;log.like2[11, 31] <- 0;log.like2[11, 32] <- 0;log.like2[11, 33] <- 0;log.like2[12, 31] <- 0;log.like2[12, 32] <- 0;log.like2[12, 33] <- 0;log.like2[13, 27] <- 0;log.like2[13, 28] <- 0;log.like2[13, 29] <- 0;log.like2[13, 30] <- 0;log.like2[13, 31] <- 0;log.like2[13, 32] <- 0;log.like2[13, 33] <- 0;log.like2[14, 29] <- 0;log.like2[14, 30] <- 0;log.like2[14, 31] <- 0;log.like2[14, 32] <- 0;log.like2[14, 33] <- 0;log.like2[15, 23] <- 0;log.like2[15, 24] <- 0;log.like2[15, 25] <- 0;log.like2[15, 26] <- 0;log.like2[15, 27] <- 0;log.like2[15, 28] <- 0;log.like2[15, 29] <- 0;log.like2[15, 30] <- 0;log.like2[15, 31] <- 0;log.like2[15, 32] <- 0;log.like2[15, 33] <- 0;log.like2[16, 32] <- 0;log.like2[16, 33] <- 0;log.like2[17, 18] <- 0;log.like2[17, 19] <- 0;log.like2[17, 20] <- 0;log.like2[17, 21] <- 0;log.like2[17, 22] <- 0;log.like2[17, 23] <- 0;log.like2[17, 24] <- 0;log.like2[17, 25] <- 0;log.like2[17, 26] <- 0;log.like2[17, 27] <- 0;log.like2[17, 28] <- 0;log.like2[17, 29] <- 0;log.like2[17, 30] <- 0;log.like2[17, 31] <- 0;log.like2[17, 32] <- 0;log.like2[17, 33] <- 0;log.like2[18, 25] <- 0;log.like2[18, 26] <- 0;log.like2[18, 27] <- 0;log.like2[18, 28] <- 0;log.like2[18, 29] <- 0;log.like2[18, 30] <- 0;log.like2[18, 31] <- 0;log.like2[18, 32] <- 0;log.like2[18, 33] <- 0;log.like2[19, 30] <- 0;log.like2[19, 31] <- 0;log.like2[19, 32] <- 0;log.like2[19, 33] <- 0;log.like2[20, 16] <- 0;log.like2[20, 17] <- 0;log.like2[20, 18] <- 0;log.like2[20, 19] <- 0;log.like2[20, 20] <- 0;log.like2[20, 21] <- 0;log.like2[20, 22] <- 0;log.like2[20, 23] <- 0;log.like2[20, 24] <- 0;log.like2[20, 25] <- 0;log.like2[20, 26] <- 0;log.like2[20, 27] <- 0;log.like2[20, 28] <- 0;log.like2[20, 29] <- 0;log.like2[20, 30] <- 0;log.like2[20, 31] <- 0;log.like2[20, 32] <- 0;log.like2[20, 33] <- 0;log.like2[21, 31] <- 0;log.like2[21, 32] <- 0;log.like2[21, 33] <- 0;log.like2[22, 30] <- 0;log.like2[22, 31] <- 0;log.like2[22, 32] <- 0;log.like2[22, 33] <- 0;log.like2[24, 31] <- 0;log.like2[24, 32] <- 0;log.like2[24, 33] <- 0;log.like2[25, 31] <- 0;log.like2[25, 32] <- 0;log.like2[25, 33] <- 0;log.like2[26, 31] <- 0;log.like2[26, 32] <- 0;log.like2[26, 33] <- 0;log.like2[27, 30] <- 0;log.like2[27, 31] <- 0;log.like2[27, 32] <- 0;log.like2[27, 33] <- 0;log.like2[28, 22] <- 0;log.like2[28, 23] <- 0;log.like2[28, 24] <- 0;log.like2[28, 25] <- 0;log.like2[28, 26] <- 0;log.like2[28, 27] <- 0;log.like2[28, 28] <- 0;log.like2[28, 29] <- 0;log.like2[28, 30] <- 0;log.like2[28, 31] <- 0;log.like2[28, 32] <- 0;log.like2[28, 33] <- 0;log.like2[29, 31] <- 0;log.like2[29, 32] <- 0;log.like2[29, 33] <- 0;log.like2[30, 15] <- 0;log.like2[30, 16] <- 0;log.like2[30, 17] <- 0;log.like2[30, 18] <- 0;log.like2[30, 19] <- 0;log.like2[30, 20] <- 0;log.like2[30, 21] <- 0;log.like2[30, 22] <- 0;log.like2[30, 23] <- 0;log.like2[30, 24] <- 0;log.like2[30, 25] <- 0;log.like2[30, 26] <- 0;log.like2[30, 27] <- 0;log.like2[30, 28] <- 0;log.like2[30, 29] <- 0;log.like2[30, 30] <- 0;log.like2[30, 31] <- 0;log.like2[30, 32] <- 0;log.like2[30, 33] <- 0;log.like2[31, 21] <- 0;log.like2[31, 22] <- 0;log.like2[31, 23] <- 0;log.like2[31, 24] <- 0;log.like2[31, 25] <- 0;log.like2[31, 26] <- 0;log.like2[31, 27] <- 0;log.like2[31, 28] <- 0;log.like2[31, 29] <- 0;log.like2[31, 30] <- 0;log.like2[31, 31] <- 0;log.like2[31, 32] <- 0;log.like2[31, 33] <- 0;log.like2[32, 20] <- 0;log.like2[32, 21] <- 0;log.like2[32, 22] <- 0;log.like2[32, 23] <- 0;log.like2[32, 24] <- 0;log.like2[32, 25] <- 0;log.like2[32, 26] <- 0;log.like2[32, 27] <- 0;log.like2[32, 28] <- 0;log.like2[32, 29] <- 0;log.like2[32, 30] <- 0;log.like2[32, 31] <- 0;log.like2[32, 32] <- 0;log.like2[32, 33] <- 0;log.like2[33, 16] <- 0;log.like2[33, 17] <- 0;log.like2[33, 18] <- 0;log.like2[33, 19] <- 0;log.like2[33, 20] <- 0;log.like2[33, 21] <- 0;log.like2[33, 22] <- 0;log.like2[33, 23] <- 0;log.like2[33, 24] <- 0;log.like2[33, 25] <- 0;log.like2[33, 26] <- 0;log.like2[33, 27] <- 0;log.like2[33, 28] <- 0;log.like2[33, 29] <- 0;log.like2[33, 30] <- 0;log.like2[33, 31] <- 0;log.like2[33, 32] <- 0;log.like2[33, 33] <- 0;log.like2[34, 14] <- 0;log.like2[34, 15] <- 0;log.like2[34, 16] <- 0;log.like2[34, 17] <- 0;log.like2[34, 18] <- 0;log.like2[34, 19] <- 0;log.like2[34, 20] <- 0;log.like2[34, 21] <- 0;log.like2[34, 22] <- 0;log.like2[34, 23] <- 0;log.like2[34, 24] <- 0;log.like2[34, 25] <- 0;log.like2[34, 26] <- 0;log.like2[34, 27] <- 0;log.like2[34, 28] <- 0;log.like2[34, 29] <- 0;log.like2[34, 30] <- 0;log.like2[34, 31] <- 0;log.like2[34, 32] <- 0;log.like2[34, 33] <- 0;log.like2[35, 14] <- 0;log.like2[35, 15] <- 0;log.like2[35, 16] <- 0;log.like2[35, 17] <- 0;log.like2[35, 18] <- 0;log.like2[35, 19] <- 0;log.like2[35, 20] <- 0;log.like2[35, 21] <- 0;log.like2[35, 22] <- 0;log.like2[35, 23] <- 0;log.like2[35, 24] <- 0;log.like2[35, 25] <- 0;log.like2[35, 26] <- 0;log.like2[35, 27] <- 0;log.like2[35, 28] <- 0;log.like2[35, 29] <- 0;log.like2[35, 30] <- 0;log.like2[35, 31] <- 0;log.like2[35, 32] <- 0;log.like2[35, 33] <- 0;log.like2[36, 24] <- 0;log.like2[36, 25] <- 0;log.like2[36, 26] <- 0;log.like2[36, 27] <- 0;log.like2[36, 28] <- 0;log.like2[36, 29] <- 0;log.like2[36, 30] <- 0;log.like2[36, 31] <- 0;log.like2[36, 32] <- 0;log.like2[36, 33] <- 0;log.like2[37, 20] <- 0;log.like2[37, 21] <- 0;log.like2[37, 22] <- 0;log.like2[37, 23] <- 0;log.like2[37, 24] <- 0;log.like2[37, 25] <- 0;log.like2[37, 26] <- 0;log.like2[37, 27] <- 0;log.like2[37, 28] <- 0;log.like2[37, 29] <- 0;log.like2[37, 30] <- 0;log.like2[37, 31] <- 0;log.like2[37, 32] <- 0;log.like2[37, 33] <- 0;log.like2[38, 30] <- 0;log.like2[38, 31] <- 0;log.like2[38, 32] <- 0;log.like2[38, 33] <- 0;log.like2[39, 11] <- 0;log.like2[39, 12] <- 0;log.like2[39, 13] <- 0;log.like2[39, 14] <- 0;log.like2[39, 15] <- 0;log.like2[39, 16] <- 0;log.like2[39, 17] <- 0;log.like2[39, 18] <- 0;log.like2[39, 19] <- 0;log.like2[39, 20] <- 0;log.like2[39, 21] <- 0;log.like2[39, 22] <- 0;log.like2[39, 23] <- 0;log.like2[39, 24] <- 0;log.like2[39, 25] <- 0;log.like2[39, 26] <- 0;log.like2[39, 27] <- 0;log.like2[39, 28] <- 0;log.like2[39, 29] <- 0;log.like2[39, 30] <- 0;log.like2[39, 31] <- 0;log.like2[39, 32] <- 0;log.like2[39, 33] <- 0;log.like2[40, 28] <- 0;log.like2[40, 29] <- 0;log.like2[40, 30] <- 0;log.like2[40, 31] <- 0;log.like2[40, 32] <- 0;log.like2[40, 33] <- 0;log.like2[41, 29] <- 0;log.like2[41, 30] <- 0;log.like2[41, 31] <- 0;log.like2[41, 32] <- 0;log.like2[41, 33] <- 0;log.like2[42, 25] <- 0;log.like2[42, 26] <- 0;log.like2[42, 27] <- 0;log.like2[42, 28] <- 0;log.like2[42, 29] <- 0;log.like2[42, 30] <- 0;log.like2[42, 31] <- 0;log.like2[42, 32] <- 0;log.like2[42, 33] <- 0;log.like2[43, 28] <- 0;log.like2[43, 29] <- 0;log.like2[43, 30] <- 0;log.like2[43, 31] <- 0;log.like2[43, 32] <- 0;log.like2[43, 33] <- 0;log.like2[44, 26] <- 0;log.like2[44, 27] <- 0;log.like2[44, 28] <- 0;log.like2[44, 29] <- 0;log.like2[44, 30] <- 0;log.like2[44, 31] <- 0;log.like2[44, 32] <- 0;log.like2[44, 33] <- 0;log.like2[45, 28] <- 0;log.like2[45, 29] <- 0;log.like2[45, 30] <- 0;log.like2[45, 31] <- 0;log.like2[45, 32] <- 0;log.like2[45, 33] <- 0;log.like2[46, 6] <- 0;log.like2[46, 7] <- 0;log.like2[46, 8] <- 0;log.like2[46, 9] <- 0;log.like2[46, 10] <- 0;log.like2[46, 11] <- 0;log.like2[46, 12] <- 0;log.like2[46, 13] <- 0;log.like2[46, 14] <- 0;log.like2[46, 15] <- 0;log.like2[46, 16] <- 0;log.like2[46, 17] <- 0;log.like2[46, 18] <- 0;log.like2[46, 19] <- 0;log.like2[46, 20] <- 0;log.like2[46, 21] <- 0;log.like2[46, 22] <- 0;log.like2[46, 23] <- 0;log.like2[46, 24] <- 0;log.like2[46, 25] <- 0;log.like2[46, 26] <- 0;log.like2[46, 27] <- 0;log.like2[46, 28] <- 0;log.like2[46, 29] <- 0;log.like2[46, 30] <- 0;log.like2[46, 31] <- 0;log.like2[46, 32] <- 0;log.like2[46, 33] <- 0;log.like2[47, 28] <- 0;log.like2[47, 29] <- 0;log.like2[47, 30] <- 0;log.like2[47, 31] <- 0;log.like2[47, 32] <- 0;log.like2[47, 33] <- 0;log.like2[48, 24] <- 0;log.like2[48, 25] <- 0;log.like2[48, 26] <- 0;log.like2[48, 27] <- 0;log.like2[48, 28] <- 0;log.like2[48, 29] <- 0;log.like2[48, 30] <- 0;log.like2[48, 31] <- 0;log.like2[48, 32] <- 0;log.like2[48, 33] <- 0;log.like2[49, 32] <- 0;log.like2[49, 33] <- 0;log.like2[50, 32] <- 0;log.like2[50, 33] <- 0;log.like2[51, 27] <- 0;log.like2[51, 28] <- 0;log.like2[51, 29] <- 0;log.like2[51, 30] <- 0;log.like2[51, 31] <- 0;log.like2[51, 32] <- 0;log.like2[51, 33] <- 0;log.like2[52, 24] <- 0;log.like2[52, 25] <- 0;log.like2[52, 26] <- 0;log.like2[52, 27] <- 0;log.like2[52, 28] <- 0;log.like2[52, 29] <- 0;log.like2[52, 30] <- 0;log.like2[52, 31] <- 0;log.like2[52, 32] <- 0;log.like2[52, 33] <- 0;log.like2[53, 24] <- 0;log.like2[53, 25] <- 0;log.like2[53, 26] <- 0;log.like2[53, 27] <- 0;log.like2[53, 28] <- 0;log.like2[53, 29] <- 0;log.like2[53, 30] <- 0;log.like2[53, 31] <- 0;log.like2[53, 32] <- 0;log.like2[53, 33] <- 0;log.like2[54, 27] <- 0;log.like2[54, 28] <- 0;log.like2[54, 29] <- 0;log.like2[54, 30] <- 0;log.like2[54, 31] <- 0;log.like2[54, 32] <- 0;log.like2[54, 33] <- 0;log.like2[55, 24] <- 0;log.like2[55, 25] <- 0;log.like2[55, 26] <- 0;log.like2[55, 27] <- 0;log.like2[55, 28] <- 0;log.like2[55, 29] <- 0;log.like2[55, 30] <- 0;log.like2[55, 31] <- 0;log.like2[55, 32] <- 0;log.like2[55, 33] <- 0;log.like2[56, 30] <- 0;log.like2[56, 31] <- 0;log.like2[56, 32] <- 0;log.like2[56, 33] <- 0;log.like2[57, 27] <- 0;log.like2[57, 28] <- 0;log.like2[57, 29] <- 0;log.like2[57, 30] <- 0;log.like2[57, 31] <- 0;log.like2[57, 32] <- 0;log.like2[57, 33] <- 0;log.like2[58, 25] <- 0;log.like2[58, 26] <- 0;log.like2[58, 27] <- 0;log.like2[58, 28] <- 0;log.like2[58, 29] <- 0;log.like2[58, 30] <- 0;log.like2[58, 31] <- 0;log.like2[58, 32] <- 0;log.like2[58, 33] <- 0;log.like2[59, 31] <- 0;log.like2[59, 32] <- 0;log.like2[59, 33] <- 0;log.like2[60, 31] <- 0;log.like2[60, 32] <- 0;log.like2[60, 33] <- 0;log.like2[61, 30] <- 0;log.like2[61, 31] <- 0;log.like2[61, 32] <- 0;log.like2[61, 33] <- 0;log.like2[62, 31] <- 0;log.like2[62, 32] <- 0;log.like2[62, 33] <- 0;log.like2[63, 31] <- 0;log.like2[63, 32] <- 0;log.like2[63, 33] <- 0;log.like2[64, 32] <- 0;log.like2[64, 33] <- 0;log.like2[65, 31] <- 0;log.like2[65, 32] <- 0;log.like2[65, 33] <- 0;log.like2[66, 32] <- 0;log.like2[66, 33] <- 0;log.like2[67, 31] <- 0;log.like2[67, 32] <- 0;log.like2[67, 33] <- 0;log.like2[68, 27] <- 0;log.like2[68, 28] <- 0;log.like2[68, 29] <- 0;log.like2[68, 30] <- 0;log.like2[68, 31] <- 0;log.like2[68, 32] <- 0;log.like2[68, 33] <- 0;log.like2[69, 32] <- 0;log.like2[69, 33] <- 0;log.like2[70, 29] <- 0;log.like2[70, 30] <- 0;log.like2[70, 31] <- 0;log.like2[70, 32] <- 0;log.like2[70, 33] <- 0;log.like2[71, 30] <- 0;log.like2[71, 31] <- 0;log.like2[71, 32] <- 0;log.like2[71, 33] <- 0;log.like2[72, 29] <- 0;log.like2[72, 30] <- 0;log.like2[72, 31] <- 0;log.like2[72, 32] <- 0;log.like2[72, 33] <- 0;log.like2[73, 31] <- 0;log.like2[73, 32] <- 0;log.like2[73, 33] <- 0;
for(i in 1:nind){  
      log.like.id[i]<- sum(log.like2[i, ]) 
  }
    
#PARAMETER PRIORS 
beta.int ~ dnorm(0,0.01)
beta.rh.var ~ dnorm(0,0.01)
beta.rh.var2 ~ dnorm(0,0.01)
for(f in 1:8){
  tau.rh.var[f] <- pow(sigma.rh.var[f],-2)
  sigma.rh.var[f] ~ dunif(0,5)
}

mu.mn.rh.var ~ dnorm(0,0.01)
rho.mn.rh.var ~ dnorm(0,0.01)

for(f in 1:2){
  beta.fac[f] ~ dnorm(0,0.01)
}
beta.fac[3] <- 0



for(i in 1:ntreat){
  alpha.Treat[i] ~ dnorm(beta.int,tau.treat)
}
tau.treat <- pow(sigma.treat,-2)
sigma.treat ~ dunif(0,25) 

}
",file="scripts/jags/jags_models/rh.var/nestmodel_quadratic.txt")


#########################################################################
# linear
cat("
model {
    
#####LIKELIHOOD
for(i in 1:nind){  
  log.like2[i, 1] <- 0
  for(t in (first[i]+1):last[i]){
    
    enc.hist[i,t] ~ dbern(S[i,t-1]*enc.hist[i,t-1])   
    
    #model for daily survival  
    logit(S[i,t-1]) <- alpha.Treat[egg.treatment[i]]
    
    + beta.fac[facility[i]] + beta.rh.var*rh.var[i,t-1] 

    rh.var[i,t] ~ dnorm(mn.rh.var[i,t],tau.rh.var[egg.treatment[i]])

    mn.rh.var[i,t] <- mu.mn.rh.var + rho.mn.rh.var*(rh.var[i,t-1]-mu.mn.rh.var)
   #log-likelihood calc
  #log.like[((cumsumm[i] +t)-1)]<-logdensity.bin(enc.hist[i,t], S[i,t-1]*enc.hist[i,t-1] , 1) 
     log.like2[i, t] <- logdensity.bin(enc.hist[i,t], S[i,t-1]*enc.hist[i,t-1] , 1) 
  }
}
log.like2[1, 32] <- 0;log.like2[1, 33] <- 0;log.like2[2, 33] <- 0;log.like2[3, 33] <- 0;log.like2[4, 32] <- 0;log.like2[4, 33] <- 0;log.like2[5, 32] <- 0;log.like2[5, 33] <- 0;log.like2[6, 32] <- 0;log.like2[6, 33] <- 0;log.like2[7, 30] <- 0;log.like2[7, 31] <- 0;log.like2[7, 32] <- 0;log.like2[7, 33] <- 0;log.like2[8, 28] <- 0;log.like2[8, 29] <- 0;log.like2[8, 30] <- 0;log.like2[8, 31] <- 0;log.like2[8, 32] <- 0;log.like2[8, 33] <- 0;log.like2[9, 20] <- 0;log.like2[9, 21] <- 0;log.like2[9, 22] <- 0;log.like2[9, 23] <- 0;log.like2[9, 24] <- 0;log.like2[9, 25] <- 0;log.like2[9, 26] <- 0;log.like2[9, 27] <- 0;log.like2[9, 28] <- 0;log.like2[9, 29] <- 0;log.like2[9, 30] <- 0;log.like2[9, 31] <- 0;log.like2[9, 32] <- 0;log.like2[9, 33] <- 0;log.like2[10, 30] <- 0;log.like2[10, 31] <- 0;log.like2[10, 32] <- 0;log.like2[10, 33] <- 0;log.like2[11, 30] <- 0;log.like2[11, 31] <- 0;log.like2[11, 32] <- 0;log.like2[11, 33] <- 0;log.like2[12, 31] <- 0;log.like2[12, 32] <- 0;log.like2[12, 33] <- 0;log.like2[13, 27] <- 0;log.like2[13, 28] <- 0;log.like2[13, 29] <- 0;log.like2[13, 30] <- 0;log.like2[13, 31] <- 0;log.like2[13, 32] <- 0;log.like2[13, 33] <- 0;log.like2[14, 29] <- 0;log.like2[14, 30] <- 0;log.like2[14, 31] <- 0;log.like2[14, 32] <- 0;log.like2[14, 33] <- 0;log.like2[15, 23] <- 0;log.like2[15, 24] <- 0;log.like2[15, 25] <- 0;log.like2[15, 26] <- 0;log.like2[15, 27] <- 0;log.like2[15, 28] <- 0;log.like2[15, 29] <- 0;log.like2[15, 30] <- 0;log.like2[15, 31] <- 0;log.like2[15, 32] <- 0;log.like2[15, 33] <- 0;log.like2[16, 32] <- 0;log.like2[16, 33] <- 0;log.like2[17, 18] <- 0;log.like2[17, 19] <- 0;log.like2[17, 20] <- 0;log.like2[17, 21] <- 0;log.like2[17, 22] <- 0;log.like2[17, 23] <- 0;log.like2[17, 24] <- 0;log.like2[17, 25] <- 0;log.like2[17, 26] <- 0;log.like2[17, 27] <- 0;log.like2[17, 28] <- 0;log.like2[17, 29] <- 0;log.like2[17, 30] <- 0;log.like2[17, 31] <- 0;log.like2[17, 32] <- 0;log.like2[17, 33] <- 0;log.like2[18, 25] <- 0;log.like2[18, 26] <- 0;log.like2[18, 27] <- 0;log.like2[18, 28] <- 0;log.like2[18, 29] <- 0;log.like2[18, 30] <- 0;log.like2[18, 31] <- 0;log.like2[18, 32] <- 0;log.like2[18, 33] <- 0;log.like2[19, 30] <- 0;log.like2[19, 31] <- 0;log.like2[19, 32] <- 0;log.like2[19, 33] <- 0;log.like2[20, 16] <- 0;log.like2[20, 17] <- 0;log.like2[20, 18] <- 0;log.like2[20, 19] <- 0;log.like2[20, 20] <- 0;log.like2[20, 21] <- 0;log.like2[20, 22] <- 0;log.like2[20, 23] <- 0;log.like2[20, 24] <- 0;log.like2[20, 25] <- 0;log.like2[20, 26] <- 0;log.like2[20, 27] <- 0;log.like2[20, 28] <- 0;log.like2[20, 29] <- 0;log.like2[20, 30] <- 0;log.like2[20, 31] <- 0;log.like2[20, 32] <- 0;log.like2[20, 33] <- 0;log.like2[21, 31] <- 0;log.like2[21, 32] <- 0;log.like2[21, 33] <- 0;log.like2[22, 30] <- 0;log.like2[22, 31] <- 0;log.like2[22, 32] <- 0;log.like2[22, 33] <- 0;log.like2[24, 31] <- 0;log.like2[24, 32] <- 0;log.like2[24, 33] <- 0;log.like2[25, 31] <- 0;log.like2[25, 32] <- 0;log.like2[25, 33] <- 0;log.like2[26, 31] <- 0;log.like2[26, 32] <- 0;log.like2[26, 33] <- 0;log.like2[27, 30] <- 0;log.like2[27, 31] <- 0;log.like2[27, 32] <- 0;log.like2[27, 33] <- 0;log.like2[28, 22] <- 0;log.like2[28, 23] <- 0;log.like2[28, 24] <- 0;log.like2[28, 25] <- 0;log.like2[28, 26] <- 0;log.like2[28, 27] <- 0;log.like2[28, 28] <- 0;log.like2[28, 29] <- 0;log.like2[28, 30] <- 0;log.like2[28, 31] <- 0;log.like2[28, 32] <- 0;log.like2[28, 33] <- 0;log.like2[29, 31] <- 0;log.like2[29, 32] <- 0;log.like2[29, 33] <- 0;log.like2[30, 15] <- 0;log.like2[30, 16] <- 0;log.like2[30, 17] <- 0;log.like2[30, 18] <- 0;log.like2[30, 19] <- 0;log.like2[30, 20] <- 0;log.like2[30, 21] <- 0;log.like2[30, 22] <- 0;log.like2[30, 23] <- 0;log.like2[30, 24] <- 0;log.like2[30, 25] <- 0;log.like2[30, 26] <- 0;log.like2[30, 27] <- 0;log.like2[30, 28] <- 0;log.like2[30, 29] <- 0;log.like2[30, 30] <- 0;log.like2[30, 31] <- 0;log.like2[30, 32] <- 0;log.like2[30, 33] <- 0;log.like2[31, 21] <- 0;log.like2[31, 22] <- 0;log.like2[31, 23] <- 0;log.like2[31, 24] <- 0;log.like2[31, 25] <- 0;log.like2[31, 26] <- 0;log.like2[31, 27] <- 0;log.like2[31, 28] <- 0;log.like2[31, 29] <- 0;log.like2[31, 30] <- 0;log.like2[31, 31] <- 0;log.like2[31, 32] <- 0;log.like2[31, 33] <- 0;log.like2[32, 20] <- 0;log.like2[32, 21] <- 0;log.like2[32, 22] <- 0;log.like2[32, 23] <- 0;log.like2[32, 24] <- 0;log.like2[32, 25] <- 0;log.like2[32, 26] <- 0;log.like2[32, 27] <- 0;log.like2[32, 28] <- 0;log.like2[32, 29] <- 0;log.like2[32, 30] <- 0;log.like2[32, 31] <- 0;log.like2[32, 32] <- 0;log.like2[32, 33] <- 0;log.like2[33, 16] <- 0;log.like2[33, 17] <- 0;log.like2[33, 18] <- 0;log.like2[33, 19] <- 0;log.like2[33, 20] <- 0;log.like2[33, 21] <- 0;log.like2[33, 22] <- 0;log.like2[33, 23] <- 0;log.like2[33, 24] <- 0;log.like2[33, 25] <- 0;log.like2[33, 26] <- 0;log.like2[33, 27] <- 0;log.like2[33, 28] <- 0;log.like2[33, 29] <- 0;log.like2[33, 30] <- 0;log.like2[33, 31] <- 0;log.like2[33, 32] <- 0;log.like2[33, 33] <- 0;log.like2[34, 14] <- 0;log.like2[34, 15] <- 0;log.like2[34, 16] <- 0;log.like2[34, 17] <- 0;log.like2[34, 18] <- 0;log.like2[34, 19] <- 0;log.like2[34, 20] <- 0;log.like2[34, 21] <- 0;log.like2[34, 22] <- 0;log.like2[34, 23] <- 0;log.like2[34, 24] <- 0;log.like2[34, 25] <- 0;log.like2[34, 26] <- 0;log.like2[34, 27] <- 0;log.like2[34, 28] <- 0;log.like2[34, 29] <- 0;log.like2[34, 30] <- 0;log.like2[34, 31] <- 0;log.like2[34, 32] <- 0;log.like2[34, 33] <- 0;log.like2[35, 14] <- 0;log.like2[35, 15] <- 0;log.like2[35, 16] <- 0;log.like2[35, 17] <- 0;log.like2[35, 18] <- 0;log.like2[35, 19] <- 0;log.like2[35, 20] <- 0;log.like2[35, 21] <- 0;log.like2[35, 22] <- 0;log.like2[35, 23] <- 0;log.like2[35, 24] <- 0;log.like2[35, 25] <- 0;log.like2[35, 26] <- 0;log.like2[35, 27] <- 0;log.like2[35, 28] <- 0;log.like2[35, 29] <- 0;log.like2[35, 30] <- 0;log.like2[35, 31] <- 0;log.like2[35, 32] <- 0;log.like2[35, 33] <- 0;log.like2[36, 24] <- 0;log.like2[36, 25] <- 0;log.like2[36, 26] <- 0;log.like2[36, 27] <- 0;log.like2[36, 28] <- 0;log.like2[36, 29] <- 0;log.like2[36, 30] <- 0;log.like2[36, 31] <- 0;log.like2[36, 32] <- 0;log.like2[36, 33] <- 0;log.like2[37, 20] <- 0;log.like2[37, 21] <- 0;log.like2[37, 22] <- 0;log.like2[37, 23] <- 0;log.like2[37, 24] <- 0;log.like2[37, 25] <- 0;log.like2[37, 26] <- 0;log.like2[37, 27] <- 0;log.like2[37, 28] <- 0;log.like2[37, 29] <- 0;log.like2[37, 30] <- 0;log.like2[37, 31] <- 0;log.like2[37, 32] <- 0;log.like2[37, 33] <- 0;log.like2[38, 30] <- 0;log.like2[38, 31] <- 0;log.like2[38, 32] <- 0;log.like2[38, 33] <- 0;log.like2[39, 11] <- 0;log.like2[39, 12] <- 0;log.like2[39, 13] <- 0;log.like2[39, 14] <- 0;log.like2[39, 15] <- 0;log.like2[39, 16] <- 0;log.like2[39, 17] <- 0;log.like2[39, 18] <- 0;log.like2[39, 19] <- 0;log.like2[39, 20] <- 0;log.like2[39, 21] <- 0;log.like2[39, 22] <- 0;log.like2[39, 23] <- 0;log.like2[39, 24] <- 0;log.like2[39, 25] <- 0;log.like2[39, 26] <- 0;log.like2[39, 27] <- 0;log.like2[39, 28] <- 0;log.like2[39, 29] <- 0;log.like2[39, 30] <- 0;log.like2[39, 31] <- 0;log.like2[39, 32] <- 0;log.like2[39, 33] <- 0;log.like2[40, 28] <- 0;log.like2[40, 29] <- 0;log.like2[40, 30] <- 0;log.like2[40, 31] <- 0;log.like2[40, 32] <- 0;log.like2[40, 33] <- 0;log.like2[41, 29] <- 0;log.like2[41, 30] <- 0;log.like2[41, 31] <- 0;log.like2[41, 32] <- 0;log.like2[41, 33] <- 0;log.like2[42, 25] <- 0;log.like2[42, 26] <- 0;log.like2[42, 27] <- 0;log.like2[42, 28] <- 0;log.like2[42, 29] <- 0;log.like2[42, 30] <- 0;log.like2[42, 31] <- 0;log.like2[42, 32] <- 0;log.like2[42, 33] <- 0;log.like2[43, 28] <- 0;log.like2[43, 29] <- 0;log.like2[43, 30] <- 0;log.like2[43, 31] <- 0;log.like2[43, 32] <- 0;log.like2[43, 33] <- 0;log.like2[44, 26] <- 0;log.like2[44, 27] <- 0;log.like2[44, 28] <- 0;log.like2[44, 29] <- 0;log.like2[44, 30] <- 0;log.like2[44, 31] <- 0;log.like2[44, 32] <- 0;log.like2[44, 33] <- 0;log.like2[45, 28] <- 0;log.like2[45, 29] <- 0;log.like2[45, 30] <- 0;log.like2[45, 31] <- 0;log.like2[45, 32] <- 0;log.like2[45, 33] <- 0;log.like2[46, 6] <- 0;log.like2[46, 7] <- 0;log.like2[46, 8] <- 0;log.like2[46, 9] <- 0;log.like2[46, 10] <- 0;log.like2[46, 11] <- 0;log.like2[46, 12] <- 0;log.like2[46, 13] <- 0;log.like2[46, 14] <- 0;log.like2[46, 15] <- 0;log.like2[46, 16] <- 0;log.like2[46, 17] <- 0;log.like2[46, 18] <- 0;log.like2[46, 19] <- 0;log.like2[46, 20] <- 0;log.like2[46, 21] <- 0;log.like2[46, 22] <- 0;log.like2[46, 23] <- 0;log.like2[46, 24] <- 0;log.like2[46, 25] <- 0;log.like2[46, 26] <- 0;log.like2[46, 27] <- 0;log.like2[46, 28] <- 0;log.like2[46, 29] <- 0;log.like2[46, 30] <- 0;log.like2[46, 31] <- 0;log.like2[46, 32] <- 0;log.like2[46, 33] <- 0;log.like2[47, 28] <- 0;log.like2[47, 29] <- 0;log.like2[47, 30] <- 0;log.like2[47, 31] <- 0;log.like2[47, 32] <- 0;log.like2[47, 33] <- 0;log.like2[48, 24] <- 0;log.like2[48, 25] <- 0;log.like2[48, 26] <- 0;log.like2[48, 27] <- 0;log.like2[48, 28] <- 0;log.like2[48, 29] <- 0;log.like2[48, 30] <- 0;log.like2[48, 31] <- 0;log.like2[48, 32] <- 0;log.like2[48, 33] <- 0;log.like2[49, 32] <- 0;log.like2[49, 33] <- 0;log.like2[50, 32] <- 0;log.like2[50, 33] <- 0;log.like2[51, 27] <- 0;log.like2[51, 28] <- 0;log.like2[51, 29] <- 0;log.like2[51, 30] <- 0;log.like2[51, 31] <- 0;log.like2[51, 32] <- 0;log.like2[51, 33] <- 0;log.like2[52, 24] <- 0;log.like2[52, 25] <- 0;log.like2[52, 26] <- 0;log.like2[52, 27] <- 0;log.like2[52, 28] <- 0;log.like2[52, 29] <- 0;log.like2[52, 30] <- 0;log.like2[52, 31] <- 0;log.like2[52, 32] <- 0;log.like2[52, 33] <- 0;log.like2[53, 24] <- 0;log.like2[53, 25] <- 0;log.like2[53, 26] <- 0;log.like2[53, 27] <- 0;log.like2[53, 28] <- 0;log.like2[53, 29] <- 0;log.like2[53, 30] <- 0;log.like2[53, 31] <- 0;log.like2[53, 32] <- 0;log.like2[53, 33] <- 0;log.like2[54, 27] <- 0;log.like2[54, 28] <- 0;log.like2[54, 29] <- 0;log.like2[54, 30] <- 0;log.like2[54, 31] <- 0;log.like2[54, 32] <- 0;log.like2[54, 33] <- 0;log.like2[55, 24] <- 0;log.like2[55, 25] <- 0;log.like2[55, 26] <- 0;log.like2[55, 27] <- 0;log.like2[55, 28] <- 0;log.like2[55, 29] <- 0;log.like2[55, 30] <- 0;log.like2[55, 31] <- 0;log.like2[55, 32] <- 0;log.like2[55, 33] <- 0;log.like2[56, 30] <- 0;log.like2[56, 31] <- 0;log.like2[56, 32] <- 0;log.like2[56, 33] <- 0;log.like2[57, 27] <- 0;log.like2[57, 28] <- 0;log.like2[57, 29] <- 0;log.like2[57, 30] <- 0;log.like2[57, 31] <- 0;log.like2[57, 32] <- 0;log.like2[57, 33] <- 0;log.like2[58, 25] <- 0;log.like2[58, 26] <- 0;log.like2[58, 27] <- 0;log.like2[58, 28] <- 0;log.like2[58, 29] <- 0;log.like2[58, 30] <- 0;log.like2[58, 31] <- 0;log.like2[58, 32] <- 0;log.like2[58, 33] <- 0;log.like2[59, 31] <- 0;log.like2[59, 32] <- 0;log.like2[59, 33] <- 0;log.like2[60, 31] <- 0;log.like2[60, 32] <- 0;log.like2[60, 33] <- 0;log.like2[61, 30] <- 0;log.like2[61, 31] <- 0;log.like2[61, 32] <- 0;log.like2[61, 33] <- 0;log.like2[62, 31] <- 0;log.like2[62, 32] <- 0;log.like2[62, 33] <- 0;log.like2[63, 31] <- 0;log.like2[63, 32] <- 0;log.like2[63, 33] <- 0;log.like2[64, 32] <- 0;log.like2[64, 33] <- 0;log.like2[65, 31] <- 0;log.like2[65, 32] <- 0;log.like2[65, 33] <- 0;log.like2[66, 32] <- 0;log.like2[66, 33] <- 0;log.like2[67, 31] <- 0;log.like2[67, 32] <- 0;log.like2[67, 33] <- 0;log.like2[68, 27] <- 0;log.like2[68, 28] <- 0;log.like2[68, 29] <- 0;log.like2[68, 30] <- 0;log.like2[68, 31] <- 0;log.like2[68, 32] <- 0;log.like2[68, 33] <- 0;log.like2[69, 32] <- 0;log.like2[69, 33] <- 0;log.like2[70, 29] <- 0;log.like2[70, 30] <- 0;log.like2[70, 31] <- 0;log.like2[70, 32] <- 0;log.like2[70, 33] <- 0;log.like2[71, 30] <- 0;log.like2[71, 31] <- 0;log.like2[71, 32] <- 0;log.like2[71, 33] <- 0;log.like2[72, 29] <- 0;log.like2[72, 30] <- 0;log.like2[72, 31] <- 0;log.like2[72, 32] <- 0;log.like2[72, 33] <- 0;log.like2[73, 31] <- 0;log.like2[73, 32] <- 0;log.like2[73, 33] <- 0;
  for(i in 1:nind){  
      log.like.id[i]<- sum(log.like2[i, ]) 
  }
    
#PARAMETER PRIORS 
beta.int ~ dnorm(0,0.01)
beta.rh.var ~ dnorm(0,0.01)

for(f in 1:8){
  tau.rh.var[f] <- pow(sigma.rh.var[f],-2)
  sigma.rh.var[f] ~ dunif(0,5)
}

mu.mn.rh.var ~ dnorm(0,0.01)
rho.mn.rh.var ~ dnorm(0,0.01)

for(f in 1:2){
  beta.fac[f] ~ dnorm(0,0.01)
}
beta.fac[3] <- 0

for(i in 1:ntreat){
  alpha.Treat[i] ~ dnorm(beta.int,tau.treat)
}
tau.treat <- pow(sigma.treat,-2)
sigma.treat ~ dunif(0,25) 

}
",file="scripts/jags/jags_models/rh.var/nestmodel_linear.txt")

#########################################################################
# null
cat("
model {
    
#####LIKELIHOOD
for(i in 1:nind){  
  log.like2[i, 1] <- 0
  for(t in (first[i]+1):last[i]){
    
    enc.hist[i,t] ~ dbern(S[i,t-1]*enc.hist[i,t-1])   
    
    #model for daily survival  
    logit(S[i,t-1]) <- alpha.Treat[egg.treatment[i]]
    
    + beta.fac[facility[i]] 

 
   #log-likelihood calc
  # log.like[((cumsumm[i] +t)-1)]<-logdensity.bin(enc.hist[i,t], S[i,t-1]*enc.hist[i,t-1] , 1) 
   log.like2[i, t] <- logdensity.bin(enc.hist[i,t], S[i,t-1]*enc.hist[i,t-1] , 1) 
  }
 
}
log.like2[1, 32] <- 0;log.like2[1, 33] <- 0;log.like2[2, 33] <- 0;log.like2[3, 33] <- 0;log.like2[4, 32] <- 0;log.like2[4, 33] <- 0;log.like2[5, 32] <- 0;log.like2[5, 33] <- 0;log.like2[6, 32] <- 0;log.like2[6, 33] <- 0;log.like2[7, 30] <- 0;log.like2[7, 31] <- 0;log.like2[7, 32] <- 0;log.like2[7, 33] <- 0;log.like2[8, 28] <- 0;log.like2[8, 29] <- 0;log.like2[8, 30] <- 0;log.like2[8, 31] <- 0;log.like2[8, 32] <- 0;log.like2[8, 33] <- 0;log.like2[9, 20] <- 0;log.like2[9, 21] <- 0;log.like2[9, 22] <- 0;log.like2[9, 23] <- 0;log.like2[9, 24] <- 0;log.like2[9, 25] <- 0;log.like2[9, 26] <- 0;log.like2[9, 27] <- 0;log.like2[9, 28] <- 0;log.like2[9, 29] <- 0;log.like2[9, 30] <- 0;log.like2[9, 31] <- 0;log.like2[9, 32] <- 0;log.like2[9, 33] <- 0;log.like2[10, 30] <- 0;log.like2[10, 31] <- 0;log.like2[10, 32] <- 0;log.like2[10, 33] <- 0;log.like2[11, 30] <- 0;log.like2[11, 31] <- 0;log.like2[11, 32] <- 0;log.like2[11, 33] <- 0;log.like2[12, 31] <- 0;log.like2[12, 32] <- 0;log.like2[12, 33] <- 0;log.like2[13, 27] <- 0;log.like2[13, 28] <- 0;log.like2[13, 29] <- 0;log.like2[13, 30] <- 0;log.like2[13, 31] <- 0;log.like2[13, 32] <- 0;log.like2[13, 33] <- 0;log.like2[14, 29] <- 0;log.like2[14, 30] <- 0;log.like2[14, 31] <- 0;log.like2[14, 32] <- 0;log.like2[14, 33] <- 0;log.like2[15, 23] <- 0;log.like2[15, 24] <- 0;log.like2[15, 25] <- 0;log.like2[15, 26] <- 0;log.like2[15, 27] <- 0;log.like2[15, 28] <- 0;log.like2[15, 29] <- 0;log.like2[15, 30] <- 0;log.like2[15, 31] <- 0;log.like2[15, 32] <- 0;log.like2[15, 33] <- 0;log.like2[16, 32] <- 0;log.like2[16, 33] <- 0;log.like2[17, 18] <- 0;log.like2[17, 19] <- 0;log.like2[17, 20] <- 0;log.like2[17, 21] <- 0;log.like2[17, 22] <- 0;log.like2[17, 23] <- 0;log.like2[17, 24] <- 0;log.like2[17, 25] <- 0;log.like2[17, 26] <- 0;log.like2[17, 27] <- 0;log.like2[17, 28] <- 0;log.like2[17, 29] <- 0;log.like2[17, 30] <- 0;log.like2[17, 31] <- 0;log.like2[17, 32] <- 0;log.like2[17, 33] <- 0;log.like2[18, 25] <- 0;log.like2[18, 26] <- 0;log.like2[18, 27] <- 0;log.like2[18, 28] <- 0;log.like2[18, 29] <- 0;log.like2[18, 30] <- 0;log.like2[18, 31] <- 0;log.like2[18, 32] <- 0;log.like2[18, 33] <- 0;log.like2[19, 30] <- 0;log.like2[19, 31] <- 0;log.like2[19, 32] <- 0;log.like2[19, 33] <- 0;log.like2[20, 16] <- 0;log.like2[20, 17] <- 0;log.like2[20, 18] <- 0;log.like2[20, 19] <- 0;log.like2[20, 20] <- 0;log.like2[20, 21] <- 0;log.like2[20, 22] <- 0;log.like2[20, 23] <- 0;log.like2[20, 24] <- 0;log.like2[20, 25] <- 0;log.like2[20, 26] <- 0;log.like2[20, 27] <- 0;log.like2[20, 28] <- 0;log.like2[20, 29] <- 0;log.like2[20, 30] <- 0;log.like2[20, 31] <- 0;log.like2[20, 32] <- 0;log.like2[20, 33] <- 0;log.like2[21, 31] <- 0;log.like2[21, 32] <- 0;log.like2[21, 33] <- 0;log.like2[22, 30] <- 0;log.like2[22, 31] <- 0;log.like2[22, 32] <- 0;log.like2[22, 33] <- 0;log.like2[24, 31] <- 0;log.like2[24, 32] <- 0;log.like2[24, 33] <- 0;log.like2[25, 31] <- 0;log.like2[25, 32] <- 0;log.like2[25, 33] <- 0;log.like2[26, 31] <- 0;log.like2[26, 32] <- 0;log.like2[26, 33] <- 0;log.like2[27, 30] <- 0;log.like2[27, 31] <- 0;log.like2[27, 32] <- 0;log.like2[27, 33] <- 0;log.like2[28, 22] <- 0;log.like2[28, 23] <- 0;log.like2[28, 24] <- 0;log.like2[28, 25] <- 0;log.like2[28, 26] <- 0;log.like2[28, 27] <- 0;log.like2[28, 28] <- 0;log.like2[28, 29] <- 0;log.like2[28, 30] <- 0;log.like2[28, 31] <- 0;log.like2[28, 32] <- 0;log.like2[28, 33] <- 0;log.like2[29, 31] <- 0;log.like2[29, 32] <- 0;log.like2[29, 33] <- 0;log.like2[30, 15] <- 0;log.like2[30, 16] <- 0;log.like2[30, 17] <- 0;log.like2[30, 18] <- 0;log.like2[30, 19] <- 0;log.like2[30, 20] <- 0;log.like2[30, 21] <- 0;log.like2[30, 22] <- 0;log.like2[30, 23] <- 0;log.like2[30, 24] <- 0;log.like2[30, 25] <- 0;log.like2[30, 26] <- 0;log.like2[30, 27] <- 0;log.like2[30, 28] <- 0;log.like2[30, 29] <- 0;log.like2[30, 30] <- 0;log.like2[30, 31] <- 0;log.like2[30, 32] <- 0;log.like2[30, 33] <- 0;log.like2[31, 21] <- 0;log.like2[31, 22] <- 0;log.like2[31, 23] <- 0;log.like2[31, 24] <- 0;log.like2[31, 25] <- 0;log.like2[31, 26] <- 0;log.like2[31, 27] <- 0;log.like2[31, 28] <- 0;log.like2[31, 29] <- 0;log.like2[31, 30] <- 0;log.like2[31, 31] <- 0;log.like2[31, 32] <- 0;log.like2[31, 33] <- 0;log.like2[32, 20] <- 0;log.like2[32, 21] <- 0;log.like2[32, 22] <- 0;log.like2[32, 23] <- 0;log.like2[32, 24] <- 0;log.like2[32, 25] <- 0;log.like2[32, 26] <- 0;log.like2[32, 27] <- 0;log.like2[32, 28] <- 0;log.like2[32, 29] <- 0;log.like2[32, 30] <- 0;log.like2[32, 31] <- 0;log.like2[32, 32] <- 0;log.like2[32, 33] <- 0;log.like2[33, 16] <- 0;log.like2[33, 17] <- 0;log.like2[33, 18] <- 0;log.like2[33, 19] <- 0;log.like2[33, 20] <- 0;log.like2[33, 21] <- 0;log.like2[33, 22] <- 0;log.like2[33, 23] <- 0;log.like2[33, 24] <- 0;log.like2[33, 25] <- 0;log.like2[33, 26] <- 0;log.like2[33, 27] <- 0;log.like2[33, 28] <- 0;log.like2[33, 29] <- 0;log.like2[33, 30] <- 0;log.like2[33, 31] <- 0;log.like2[33, 32] <- 0;log.like2[33, 33] <- 0;log.like2[34, 14] <- 0;log.like2[34, 15] <- 0;log.like2[34, 16] <- 0;log.like2[34, 17] <- 0;log.like2[34, 18] <- 0;log.like2[34, 19] <- 0;log.like2[34, 20] <- 0;log.like2[34, 21] <- 0;log.like2[34, 22] <- 0;log.like2[34, 23] <- 0;log.like2[34, 24] <- 0;log.like2[34, 25] <- 0;log.like2[34, 26] <- 0;log.like2[34, 27] <- 0;log.like2[34, 28] <- 0;log.like2[34, 29] <- 0;log.like2[34, 30] <- 0;log.like2[34, 31] <- 0;log.like2[34, 32] <- 0;log.like2[34, 33] <- 0;log.like2[35, 14] <- 0;log.like2[35, 15] <- 0;log.like2[35, 16] <- 0;log.like2[35, 17] <- 0;log.like2[35, 18] <- 0;log.like2[35, 19] <- 0;log.like2[35, 20] <- 0;log.like2[35, 21] <- 0;log.like2[35, 22] <- 0;log.like2[35, 23] <- 0;log.like2[35, 24] <- 0;log.like2[35, 25] <- 0;log.like2[35, 26] <- 0;log.like2[35, 27] <- 0;log.like2[35, 28] <- 0;log.like2[35, 29] <- 0;log.like2[35, 30] <- 0;log.like2[35, 31] <- 0;log.like2[35, 32] <- 0;log.like2[35, 33] <- 0;log.like2[36, 24] <- 0;log.like2[36, 25] <- 0;log.like2[36, 26] <- 0;log.like2[36, 27] <- 0;log.like2[36, 28] <- 0;log.like2[36, 29] <- 0;log.like2[36, 30] <- 0;log.like2[36, 31] <- 0;log.like2[36, 32] <- 0;log.like2[36, 33] <- 0;log.like2[37, 20] <- 0;log.like2[37, 21] <- 0;log.like2[37, 22] <- 0;log.like2[37, 23] <- 0;log.like2[37, 24] <- 0;log.like2[37, 25] <- 0;log.like2[37, 26] <- 0;log.like2[37, 27] <- 0;log.like2[37, 28] <- 0;log.like2[37, 29] <- 0;log.like2[37, 30] <- 0;log.like2[37, 31] <- 0;log.like2[37, 32] <- 0;log.like2[37, 33] <- 0;log.like2[38, 30] <- 0;log.like2[38, 31] <- 0;log.like2[38, 32] <- 0;log.like2[38, 33] <- 0;log.like2[39, 11] <- 0;log.like2[39, 12] <- 0;log.like2[39, 13] <- 0;log.like2[39, 14] <- 0;log.like2[39, 15] <- 0;log.like2[39, 16] <- 0;log.like2[39, 17] <- 0;log.like2[39, 18] <- 0;log.like2[39, 19] <- 0;log.like2[39, 20] <- 0;log.like2[39, 21] <- 0;log.like2[39, 22] <- 0;log.like2[39, 23] <- 0;log.like2[39, 24] <- 0;log.like2[39, 25] <- 0;log.like2[39, 26] <- 0;log.like2[39, 27] <- 0;log.like2[39, 28] <- 0;log.like2[39, 29] <- 0;log.like2[39, 30] <- 0;log.like2[39, 31] <- 0;log.like2[39, 32] <- 0;log.like2[39, 33] <- 0;log.like2[40, 28] <- 0;log.like2[40, 29] <- 0;log.like2[40, 30] <- 0;log.like2[40, 31] <- 0;log.like2[40, 32] <- 0;log.like2[40, 33] <- 0;log.like2[41, 29] <- 0;log.like2[41, 30] <- 0;log.like2[41, 31] <- 0;log.like2[41, 32] <- 0;log.like2[41, 33] <- 0;log.like2[42, 25] <- 0;log.like2[42, 26] <- 0;log.like2[42, 27] <- 0;log.like2[42, 28] <- 0;log.like2[42, 29] <- 0;log.like2[42, 30] <- 0;log.like2[42, 31] <- 0;log.like2[42, 32] <- 0;log.like2[42, 33] <- 0;log.like2[43, 28] <- 0;log.like2[43, 29] <- 0;log.like2[43, 30] <- 0;log.like2[43, 31] <- 0;log.like2[43, 32] <- 0;log.like2[43, 33] <- 0;log.like2[44, 26] <- 0;log.like2[44, 27] <- 0;log.like2[44, 28] <- 0;log.like2[44, 29] <- 0;log.like2[44, 30] <- 0;log.like2[44, 31] <- 0;log.like2[44, 32] <- 0;log.like2[44, 33] <- 0;log.like2[45, 28] <- 0;log.like2[45, 29] <- 0;log.like2[45, 30] <- 0;log.like2[45, 31] <- 0;log.like2[45, 32] <- 0;log.like2[45, 33] <- 0;log.like2[46, 6] <- 0;log.like2[46, 7] <- 0;log.like2[46, 8] <- 0;log.like2[46, 9] <- 0;log.like2[46, 10] <- 0;log.like2[46, 11] <- 0;log.like2[46, 12] <- 0;log.like2[46, 13] <- 0;log.like2[46, 14] <- 0;log.like2[46, 15] <- 0;log.like2[46, 16] <- 0;log.like2[46, 17] <- 0;log.like2[46, 18] <- 0;log.like2[46, 19] <- 0;log.like2[46, 20] <- 0;log.like2[46, 21] <- 0;log.like2[46, 22] <- 0;log.like2[46, 23] <- 0;log.like2[46, 24] <- 0;log.like2[46, 25] <- 0;log.like2[46, 26] <- 0;log.like2[46, 27] <- 0;log.like2[46, 28] <- 0;log.like2[46, 29] <- 0;log.like2[46, 30] <- 0;log.like2[46, 31] <- 0;log.like2[46, 32] <- 0;log.like2[46, 33] <- 0;log.like2[47, 28] <- 0;log.like2[47, 29] <- 0;log.like2[47, 30] <- 0;log.like2[47, 31] <- 0;log.like2[47, 32] <- 0;log.like2[47, 33] <- 0;log.like2[48, 24] <- 0;log.like2[48, 25] <- 0;log.like2[48, 26] <- 0;log.like2[48, 27] <- 0;log.like2[48, 28] <- 0;log.like2[48, 29] <- 0;log.like2[48, 30] <- 0;log.like2[48, 31] <- 0;log.like2[48, 32] <- 0;log.like2[48, 33] <- 0;log.like2[49, 32] <- 0;log.like2[49, 33] <- 0;log.like2[50, 32] <- 0;log.like2[50, 33] <- 0;log.like2[51, 27] <- 0;log.like2[51, 28] <- 0;log.like2[51, 29] <- 0;log.like2[51, 30] <- 0;log.like2[51, 31] <- 0;log.like2[51, 32] <- 0;log.like2[51, 33] <- 0;log.like2[52, 24] <- 0;log.like2[52, 25] <- 0;log.like2[52, 26] <- 0;log.like2[52, 27] <- 0;log.like2[52, 28] <- 0;log.like2[52, 29] <- 0;log.like2[52, 30] <- 0;log.like2[52, 31] <- 0;log.like2[52, 32] <- 0;log.like2[52, 33] <- 0;log.like2[53, 24] <- 0;log.like2[53, 25] <- 0;log.like2[53, 26] <- 0;log.like2[53, 27] <- 0;log.like2[53, 28] <- 0;log.like2[53, 29] <- 0;log.like2[53, 30] <- 0;log.like2[53, 31] <- 0;log.like2[53, 32] <- 0;log.like2[53, 33] <- 0;log.like2[54, 27] <- 0;log.like2[54, 28] <- 0;log.like2[54, 29] <- 0;log.like2[54, 30] <- 0;log.like2[54, 31] <- 0;log.like2[54, 32] <- 0;log.like2[54, 33] <- 0;log.like2[55, 24] <- 0;log.like2[55, 25] <- 0;log.like2[55, 26] <- 0;log.like2[55, 27] <- 0;log.like2[55, 28] <- 0;log.like2[55, 29] <- 0;log.like2[55, 30] <- 0;log.like2[55, 31] <- 0;log.like2[55, 32] <- 0;log.like2[55, 33] <- 0;log.like2[56, 30] <- 0;log.like2[56, 31] <- 0;log.like2[56, 32] <- 0;log.like2[56, 33] <- 0;log.like2[57, 27] <- 0;log.like2[57, 28] <- 0;log.like2[57, 29] <- 0;log.like2[57, 30] <- 0;log.like2[57, 31] <- 0;log.like2[57, 32] <- 0;log.like2[57, 33] <- 0;log.like2[58, 25] <- 0;log.like2[58, 26] <- 0;log.like2[58, 27] <- 0;log.like2[58, 28] <- 0;log.like2[58, 29] <- 0;log.like2[58, 30] <- 0;log.like2[58, 31] <- 0;log.like2[58, 32] <- 0;log.like2[58, 33] <- 0;log.like2[59, 31] <- 0;log.like2[59, 32] <- 0;log.like2[59, 33] <- 0;log.like2[60, 31] <- 0;log.like2[60, 32] <- 0;log.like2[60, 33] <- 0;log.like2[61, 30] <- 0;log.like2[61, 31] <- 0;log.like2[61, 32] <- 0;log.like2[61, 33] <- 0;log.like2[62, 31] <- 0;log.like2[62, 32] <- 0;log.like2[62, 33] <- 0;log.like2[63, 31] <- 0;log.like2[63, 32] <- 0;log.like2[63, 33] <- 0;log.like2[64, 32] <- 0;log.like2[64, 33] <- 0;log.like2[65, 31] <- 0;log.like2[65, 32] <- 0;log.like2[65, 33] <- 0;log.like2[66, 32] <- 0;log.like2[66, 33] <- 0;log.like2[67, 31] <- 0;log.like2[67, 32] <- 0;log.like2[67, 33] <- 0;log.like2[68, 27] <- 0;log.like2[68, 28] <- 0;log.like2[68, 29] <- 0;log.like2[68, 30] <- 0;log.like2[68, 31] <- 0;log.like2[68, 32] <- 0;log.like2[68, 33] <- 0;log.like2[69, 32] <- 0;log.like2[69, 33] <- 0;log.like2[70, 29] <- 0;log.like2[70, 30] <- 0;log.like2[70, 31] <- 0;log.like2[70, 32] <- 0;log.like2[70, 33] <- 0;log.like2[71, 30] <- 0;log.like2[71, 31] <- 0;log.like2[71, 32] <- 0;log.like2[71, 33] <- 0;log.like2[72, 29] <- 0;log.like2[72, 30] <- 0;log.like2[72, 31] <- 0;log.like2[72, 32] <- 0;log.like2[72, 33] <- 0;log.like2[73, 31] <- 0;log.like2[73, 32] <- 0;log.like2[73, 33] <- 0;
for(i in 1:nind){  
      log.like.id[i]<- sum(log.like2[i, ]) 
  }
#PARAMETER PRIORS 
beta.int ~ dnorm(0,0.01)

for(f in 1:2){
  beta.fac[f] ~ dnorm(0,0.01)
}
beta.fac[3] <- 0

#RANDOM EFFECT
for(i in 1:ntreat){
  alpha.Treat[i] ~ dnorm(beta.int,tau.treat)
}
tau.treat <- pow(sigma.treat,-2)
sigma.treat ~ dunif(0,25) 

}
",file="scripts/jags/jags_models/nestmodel_null.txt")
######################### RUN BUGS #########################

start <- Sys.time()


jagsfit.1.quadratic <- jags(data=dataset, inits=inits, parameters.to.save=parameters_quadratic, 
                  n.chains=nc, n.burnin = nb, n.iter=ni, n.thin=nt, n.adapt = na, 
                  model.file="scripts/jags/jags_models/rh.var/nestmodel_quadratic.txt", parallel=TRUE)

jagsfit.1.linear <- jags(data=dataset, inits=inits, parameters.to.save=parameters_linear, 
                          n.chains=nc, n.burnin = nb, n.iter=ni, n.thin=nt, n.adapt = na,
                          model.file="scripts/jags/jags_models/rh.var/nestmodel_linear.txt", parallel=TRUE)

if(file.exists("model_outputs/jags/nestmodel_null.RData")){
  load("model_outputs/jags/nestmodel_null.RData")
} else {
  jagsfit.1.null <- jags(data=dataset, inits=inits, parameters.to.save=parameters_null, 
                          n.chains=nc, n.burnin = nb, n.iter=ni, n.thin=nt, n.adapt = na,
                          model.file="scripts/jags/jags_models/nestmodel_null.txt", parallel=TRUE)
}
if(!file.exists("model_outputs/jags/nestmodel_null.RData")){
  save(jagsfit.1.null, file = "model_outputs/jags/nestmodel_null.RData")
}
end <- Sys.time() - start
end
#######################################################################################
#################### MODEL COMPARISON
if(file.exists("model_outputs/jags/compare_rh.var.txt")) unlink("model_outputs/jags/compare_rh.var.txt")
sink(file = "model_outputs/jags/compare_rh.var.txt")
print("--------------------------------------------------------------------")
print("Log likelihood by ID")
print("--------------------------------------------------------------------")

log.lik.jagsfit.1.null <- jagsfit.1.null$sims.list$log.like.id
waic_only_intercepts <- waic(log.lik.jagsfit.1.null)
loo_only_intercepts <- loo(log.lik.jagsfit.1.null)

log.lik.jagsfit.1.linear <- jagsfit.1.linear$sims.list$log.like.id
waic_linear <- waic(log.lik.jagsfit.1.linear)
loo_linear <- loo(log.lik.jagsfit.1.linear)

log.lik.jagsfit.1.quadratic <- jagsfit.1.quadratic$sims.list$log.like.id
waic_quadratic <- waic(log.lik.jagsfit.1.quadratic)
loo_quadratic <- loo(log.lik.jagsfit.1.quadratic)

waics <- c(
  waic_only_intercepts$estimates["elpd_waic", 1],
  waic_linear$estimates["elpd_waic", 1],
  waic_quadratic$estimates["elpd_waic", 1]
)

lpd_point <- cbind(
  loo_only_intercepts$pointwise[,"elpd_loo"], 
  loo_linear$pointwise[,"elpd_loo"],
  loo_quadratic$pointwise[,"elpd_loo"]
)
####################################################################
print("compare null model with linear and quadratic for variable rh.var ")
print("model 1=null; model2=linear; model3=quadratic")
print("waic")
print(loo_compare(waic_only_intercepts,waic_linear, waic_quadratic), simplify = F)
print("loo")
print(loo_compare(loo_only_intercepts,loo_linear, loo_quadratic), simplify = F)

print("Compare null with linear for variable rh.var")
print("model 1=null; model2=linear")
print("waic")
print(loo_compare(waic_only_intercepts,waic_linear), simplify = F)
print("loo")
print(loo_compare(loo_only_intercepts,loo_linear), simplify = F)

print("Model weights")
print(" 1) WAIC weights, 2) Pseudo-BMA weights without Bayesian bootstrap, \n
3) Pseudo-BMA+ weights with Bayesian bootstrap, and 4) Bayesian stacking weights.")

waic_wts <- exp(waics) / sum(exp(waics))
pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)
print(round(cbind(waic_wts, pbma_wts, pbma_BB_wts, stacking_wts), 2))

print("###################################################################")
print("###################################################################")

print("null model")
print("loo")
print(loo_only_intercepts)
print("waic")
print(waic_only_intercepts)


print("###################################################################")
print("###################################################################")

print("Linear model")
print("loo")
print(loo_linear)

print("waic")
print(waic_linear)



print("###################################################################")
print("###################################################################")

print("Quadratic model")
print("loo")
print(loo_quadratic)

print("waic")
print(waic_quadratic)


print("--------------------------------------------------------------------")
print("Model diagnostics")
print("--------------------------------------------------------------------")

# Model diagnostics
rn <- row.names(jagsfit.1.quadratic$summary)
para_quadratic <- rn[str_detect(rn, "log.like", negate = T)]
rn <- row.names(jagsfit.1.linear$summary)
para_linear <- rn[str_detect(rn, "log.like", negate = T)]
rn <- row.names(jagsfit.1.null$summary)
para_null <- rn[str_detect(rn, "log.like", negate = T)]


mean_hdi_quadratic <- t(round(rbind(jagsfit.1.quadratic$summary[para_quadratic,"mean"],
                                    hdi(jagsfit.1.quadratic, credMass = 0.95)[,para_quadratic]),2))
print("mean hdi quadratic")
print(mean_hdi_quadratic)

print("Summary quadratic model rh.var")
print(round(jagsfit.1.quadratic$summary[para_quadratic,], 2))


print("###################################################################")
print("###################################################################")
mean_hdi_linear <- t(round(rbind(jagsfit.1.linear$summary[para_linear,"mean"],
                                 hdi(jagsfit.1.linear, credMass = 0.95)[,para_linear]),2))
print("mean hdi linear")
print(mean_hdi_linear)

print("Summary linear model rh.var")
print(round(jagsfit.1.linear$summary[para_linear,], 2))
print("###################################################################")
print("###################################################################")

print("mean hdi null model")
mean_hdi_null <- t(round(rbind(jagsfit.1.null$summary[para_null,"mean"],
                               hdi(jagsfit.1.null, credMass = 0.95)[,para_null]),2))

print(mean_hdi_null)
print("Summary null model")
print(round(jagsfit.1.null$summary[para_null,], 2))



print("###################################################################")
print("###################################################################")
################# Diagnostic plots

pdf(file="model_outputs/jags/plots/jagsfit.rh.var_quadratic.pdf")
plot(jagsfit.1.quadratic$samples[,para_quadratic])
dev.off()
pdf(file="model_outputs/jags/plots/jagsfit.rh.var_linear.pdf")
plot(jagsfit.1.linear$samples[,para_linear])
dev.off()
if(!file.exists("model_outputs/jags/plots/jagsfit.null.pdf")){
  pdf(file="model_outputs/jags/plots/jagsfit.null.pdf")
  plot(jagsfit.1.null$samples[,para_null])
  dev.off()
}
print("#############################################################mo")
save(mean_hdi_quadratic, mean_hdi_linear, mean_hdi_null, file = "model_outputs/jags/mean_hdi/means_hdi_rh.var.RData")
print("###################################################################")
print("###################################################################")
# Create linear plots
load(file = "model_outputs/jags/mean_hdi/means_hdi_rh.var.RData")
load(file = "data/derived_data/dataset_hyp_1_40min.RData")
nvalues = 2000
expit <- function (x) exp(x)/(1 + exp(x))
b_length_new <- seq(min(dataset$rh.var,na.rm=T), max(dataset$rh.var,na.rm=T), 
                    length.out=nvalues)
lower = expit(mean_hdi_linear["beta.rh.var", 2] * b_length_new)
upper =expit(mean_hdi_linear["beta.rh.var", 3] * b_length_new)
means = expit(mean_hdi_linear["beta.rh.var", 1] * b_length_new)

means_quadratic = expit(mean_hdi_quadratic["beta.rh.var", 1] * b_length_new +
                          mean_hdi_quadratic["beta.rh.var2", 1] * b_length_new * b_length_new)

# Plot linear and quadratic functions
pdf(file="model_outputs/jags/plots/regression_rh.var_linear_versus-quadratic.pdf")
plot(x=dataset$rh.var, y=jitter(dataset$enc.hist, 0.2), ylab ="Egg survival", 
     xlab = "Humidity Variance", pch=20)

lines(y = means, x=b_length_new, col=viridis(2)[1])
lines(y=means_quadratic, x = b_length_new, col=viridis(2)[2])
legend(x=min(dataset$rh.var, na.rm=T)-0.15*min(dataset$rh.var, na.rm=T), y=0.3,  
       col=viridis(2), legend=c("linear", "quadratic"), lty=c(1,1))
dev.off()

# Plot linear functions, including hdi
pdf(file="model_outputs/jags/plots/regression_rh.var_linear_2.pdf")
plot(x=dataset$rh.var, y=jitter(dataset$enc.hist, 0.2), ylab ="Egg survival", 
     xlab = "Humidity Variance", pch=20)
polygon(c(b_length_new,rev(b_length_new)),c(lower,rev(upper)),
        col=viridis(1, alpha=0.4), border=NA)
lines(y = means, x=b_length_new, col=viridis(1)[1])

legend(x=min(dataset$rh.var, na.rm=T)-0.15*min(dataset$rh.var, na.rm=T), y=0.3,  
       col=c(viridis(1), viridis(1, alpha=0.4)), legend=c("mean", "confidence"), 
       lty=c(1,1), lwd=c(1,10))
dev.off()


fac_par <- row.names(mean_hdi_linear)[str_detect(row.names(mean_hdi_linear), "beta.fac", negate = F)]
fac_int <- mean_hdi_linear[fac_par,]
treat_par <- row.names(mean_hdi_linear)[str_detect(row.names(mean_hdi_linear), "alpha.Treat", negate = F)]
treat_int <- mean_hdi_linear[treat_par,]


ymax = expit(max(fac_int[,1], treat_int[,1]) -mean_hdi_linear["beta.int",1] + mean_hdi_linear["beta.rh.var", 1] * b_length_new)
ymin = expit(min(fac_int[,1], treat_int[,1]) -mean_hdi_linear["beta.int",1] + mean_hdi_linear["beta.rh.var", 1] * b_length_new)
# plot linear function including variation between facilities and treatments
pdf(file="model_outputs/jags/plots/regression_rh.var_linear_1.pdf")
plot(x=dataset$rh.var, y=jitter(dataset$enc.hist, 0.2), ylab ="Egg survival", 
     xlab = "Humidity Variance", pch=20)
polygon(c(b_length_new,rev(b_length_new)),c(ymin,rev(ymax)),
        col=viridis(1, alpha=0.4), border=NA)
lines(y = means, x=b_length_new, col=viridis(1)[1])
legend(x=min(dataset$rh.var, na.rm=T)-0.15*min(dataset$rh.var, na.rm=T),
       y=0.5,  
       col=c(viridis(1), viridis(1, alpha=0.4)), legend=c("mean", 
      "Variation between facilities, incubation treatments and egg type"), 
       lty=c(1,1), lwd=c(1,10), bty="n", cex=0.8)
dev.off()
print("###################################################################")
print("###################################################################")
print(paste("test equal ", test))
print(paste("cores = ", nc))
print(paste("thinning = ", nt))
print(paste("iterations = ", ni))
print(paste("burn-in = ", nb))

info <- sessionInfo()


print(info)


sink(NULL)

print("the end")
