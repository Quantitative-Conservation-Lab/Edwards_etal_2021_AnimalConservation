#-----------------------------------------------------

# Scripts to run the models in JAGS with the jagsUI package.
# These models address the hypothesis that incubation treatment affects  
# the survival of the eggs.
# Models are autoregressive 
# JAGS scripts are saved as txt files within the folder scripts/jags/jags_models
# MCMC chains are saved as RData files within the folder model_outputs/jags/ for further analyses
# The m_specis model is the NULL hypothesis and includes species as multilevel intercepts.
# Incubation treatment is a categorical variable with 6 levels: B, WC, GQF1, GQF2, P and SHC
# Script also does check convergency of the chains, do basic model diagnostics and
# compare model weights.
# The m_null model is the NULL hypothesis and includes species as multilevel intercepts.
# Pair ID is a random effect nested into species.
# Facility is a random effect with three levels: CZ, ICF and Patuxent
# Basic output is "sink" to file "model_outputs/jags/compare_jags_hyp2.txt"


#----------------Load libraries.-------------# 

if(file.exists("data/derived_data/dataset3_hyp_2_40min.RData")){
  load(file = "data/derived_data/dataset_hyp_2_40min.RData")
  load(file = "data/derived_data/dataset3_hyp_2_40min.RData")
} else{
  # if data/derived_data/dataset_list.RData is not available first run
  source("scripts/getData_for_jags_40min.R")
  rm(list=ls())
  load(file = "data/derived_data/dataset_list.RData")
  load(file = "data/derived_data/dataset3_hyp_2_40min.RData")
}
if(!file.exists("model_outputs/")) dir.create("model_outputs/")
if(!file.exists("model_outputs/jags")) dir.create("model_outputs/jags")
if(!file.exists("model_outputs/jags/plots")) dir.create("model_outputs/jags/plots")
if(!file.exists("model_outputs/jags/mean_hdi/")) dir.create("model_outputs/jags/mean_hdi/")
# 
library("jagsUI")
library(coda)
library(gdata)
library(loo)
library(stringr)
library(HDInterval)
library(forestplot)
test <- FALSE
#test <- TRUE
if(test==TRUE){
  nc <- 4
  nt <- 1#5 #thinning
  ni <- 500#50000 #iteration
  nb <- 50#ni - 40000 #burn-in
  na <- 500 # adaptation
} else{
  set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion")
  nc <- 4
  nt <- 5 #thinning
  ni <- 300000 #iteration
  nb <- ni - 80000 #burn-in
  na <- 100000 # adaptation
  (ni-nb)/5
  if(nb<=0) stop("burn in is negative")
}
#initial values
inits <- NULL
# For nested random effects see example in https://stackoverflow.com/questions/50471437/three-level-nested-mixed-effects-model

#parameters
par_m_treatment <- c("log.like","treats", "int" #alpha.B",
                             
)  


######################### RUN BUGS #########################
####### m_treatment
start <- Sys.time()
jagsfit.m_treatment <- jags(data=dataset2, inits=inits, parameters.to.save=par_m_treatment, n.chains=nc, 
                                    n.burnin = nb, n.iter=ni, n.thin=nt, n.adapt = na, 
                            model.file="scripts/jags/jags_models/m_treatment_no_random.txt", parallel=TRUE)

end <- Sys.time() - start
end
# save(jagsfit.m_treatment, file="model_outputs/jags/jagsfit.m_treatment.RData")
# rm(jagsfit.m_treatment)
######################################################################

# Calculate observed probabilities from the data
prob <- NULL
for(i in 1:6){
  zeros <- length(which(dataset2$enc.hist[dataset2$treatment==i,]==0))
  ones <- sum(dataset2$enc.hist[dataset2$treatment==i,], na.rm = T)
  prob[i] <- round(ones/(zeros+ones),2)
}
prob_by_egg <- NULL
for(i in 1:6){
  tmp <- NULL
  dtmp <- dataset2$enc.hist[dataset2$treatment==i,]
  for(n in 1:nrow(dtmp)){
    if(any(dtmp[n,]==0, na.rm=T)){
      tmp[n] <- 0
    } else {
      tmp[n] <- 1
    }
  }
  prob_by_egg[i] <- round(sum(tmp)/length(tmp), 2)
  
}



######################### RUN BUGS #########################
if(file.exists("model_outputs/jags/compare_jags_hyp2_no_random.txt")) unlink("model_outputs/jags/compare_jags_hyp2_no_random.txt")
sink(file = "model_outputs/jags/compare_jags_hyp2_no_random.txt")

print("--------------------------------------------------------------------")
print(" Log likelihood by ID")
print("--------------------------------------------------------------------")


log.lik.m_treatment <- jagsfit.m_treatment$sims.list$log.like#list of simulated params model 1
waic_m_treatment <- waic(log.lik.m_treatment)
loo_m_treatment <- loo(log.lik.m_treatment)
print("\n waic")
print(waic_m_treatment)

print("\n loo")
print(loo_m_treatment)


print("#-------------------------------------------------------------------------#")
print("Model diagnostics")
print("#")
print("#")
print("#-------------------------------------------------------------------------#")

rn <- row.names(jagsfit.m_treatment$summary)
para <- rn[str_detect(rn, "log.like", negate = T)]
print("Summary jagsfit.m_treatment")
print(round(jagsfit.m_treatment$summary[para,], 2))

mean_hdi_m_treatment <- t(round(rbind(jagsfit.m_treatment$summary[para,"mean"],
                                              hdi(jagsfit.m_treatment, credMass = 0.95)[,para]),2))
print("mean hdi jagsfit.m_treatment")
print(mean_hdi_m_treatment)

pdf(file="model_outputs/jags/plots/jagsfit.m_treatment_no_random.pdf")
plot(jagsfit.m_treatment$samples[,para])
dev.off()


print("###########################################################################")
print("###########################################################################")
print("###########################################################################")

save(jagsfit.m_treatment, 
     file = "model_outputs/jags/mean_hdi/jagsfit.m_treatment_no_random.RData")
save(mean_hdi_m_treatment, 
     file = "model_outputs/jags/mean_hdi/means_hdi_hyp2_no_random.RData")
print("###########################################################################")
# Compare to NULL model

par_m_null <- c("log.like", "int")

######################### RUN BUGS #########################
######################################################################
jagsfit.m_null <- jags(data=dataset2, inits=inits, parameters.to.save=par_m_null, n.chains=nc,
                          n.burnin = nb, n.iter=ni, n.thin=nt, n.adapt = na,
                       model.file="scripts/jags/jags_models/m_int_no_random.txt", parallel=TRUE)



log.lik.m_null <- jagsfit.m_null$sims.list$log.like#list of simulated params model 1
waic_m_null <- waic(log.lik.m_null)
loo_m_null <- loo(log.lik.m_null)


print("compare to null model")
print("loo")
print(loo_compare(loo_m_null, loo_m_treatment), simplify = F)
######################################################################
print("###########################################################################")
print("#######################################################")
print("The treatments coefficients are generated in two steps: first we get a set of temp values drawn from broad normal distributions, 
then we subtract the mean from each.")
# see http://mikemeredith.net/blog/2017/Categories_in_JAGS.htm


print("check that the treatments coefficients sum to 0 or very near to zero:")

tst <- rowSums(jagsfit.m_treatment$sims.list$treats)
print(range(tst))


print("The individual treatment coefficients tell us the difference in success for each treatment from the mean of all treatments 
(on the logit scale).")
print("The mean is given by the intercept. So to recover estimates of each treatment we need to add the intercept:")

logitT6 <- jagsfit.m_treatment$sims.list$treats + jagsfit.m_treatment$sims.list$int
T6 <- plogis(logitT6)
#hdi_T6 <- round(plogis(hdi(logitT6, credMass = 0.95)), 2)
hdi_T6 <- round(hdi(plogis(logitT6), credMass = 0.95), 2)


backtransformed <- round(data.frame(mean=round(colMeans(T6), 2), lower=t(hdi_T6)[, 1], upper=t(hdi_T6)[, 2]), 2)




untrasformed <- data.frame(par= c("int", "treats[1]", "treats[2]", 
                  "treats[3]", 
                  "treats[4]", 
                  "treats[5]", 
                  "treats[6]"), 
           mean =round(c(mean(jagsfit.m_treatment$sims.list$int), colMeans(jagsfit.m_treatment$sims.list$treats)),2),
           lower=round(c(hdi(jagsfit.m_treatment$sims.list$int)[1],hdi(jagsfit.m_treatment$sims.list$treats)[1,]),2),
           upper=round(c(hdi(jagsfit.m_treatment$sims.list$int)[2],hdi(jagsfit.m_treatment$sims.list$treats)[2,]),2))


print("###########################################################################")
print("###########################################################################")
print("#######################################################")

expit <- function (x) exp(x)/(1 + exp(x))



order_fac <- c("GQF1", "GQF2", "P", "SHC", "B", "WC")
backtransformed$names <- order_fac
backtransformed$sample_size <- table(dataset3$treatment_factor)[order_fac]
backtransformed$observed <- prob_by_egg
backtransformed <- backtransformed[order(backtransformed$observed, decreasing = T),]
print("The actual estimates are:")
print(backtransformed)
load(file="data/derived_data/dataset3_hyp_2_40min.RData")


#colnames(values_fore) <- c("mean", "lower", "upper")
x<- backtransformed$sample_size
ss <- c("Sample size", x[1], "", x[2], "", x[3], "", x[4], "", x[5], "", x[6], "")
x <- backtransformed$observed
pe <- c("Obs prob",x[1], "", x[2], "", x[3], "", x[4], "", x[5], "", x[6], "")
x <- backtransformed$names
nn <- c("Treatment", x[1], "", x[2], "", x[3], "", x[4], "", x[5], "", x[6], "")
tabletext<-cbind(
  nn,
  ss, c("", rep(c("observed", "predicted"), 6)))


png(file="model_outputs/jags/plots/forestplot_treatments_no_random.png")
forestplot(tabletext, 
           mean  = c(NA, backtransformed$observed[1], backtransformed[1,"mean"],
                     backtransformed$observed[2], backtransformed[2,"mean"],
                     backtransformed$observed[3], backtransformed[3,"mean"],
                     backtransformed$observed[4], backtransformed[4,"mean"],
                     backtransformed$observed[5], backtransformed[5,"mean"],
                     backtransformed$observed[6], backtransformed[6,"mean"]), 
           lower = c(NA, backtransformed$observed[1], backtransformed[1,"lower"],
                     backtransformed$observed[2], backtransformed[2,"lower"],
                     backtransformed$observed[3], backtransformed[3,"lower"],
                     backtransformed$observed[4], backtransformed[4,"lower"],
                     backtransformed$observed[5], backtransformed[5,"lower"],
                     backtransformed$observed[6], backtransformed[6,"lower"]),
           upper = c(NA, backtransformed$observed[1], backtransformed[1,"upper"],
                     backtransformed$observed[2], backtransformed[2,"upper"],
                     backtransformed$observed[3], backtransformed[3,"upper"],
                     backtransformed$observed[4], backtransformed[4,"upper"],
                     backtransformed$observed[5], backtransformed[5,"upper"],
                     backtransformed$observed[6], backtransformed[6,"upper"]),
           hrzl_lines=list("5" = gpar(lwd=85, lineend="butt", columns=c(1:4), col="#99999922"), 
                           "9" = gpar(lwd=85, lineend="butt", columns=c(1:4), col="#99999922"),
                           "13" = gpar(lwd=85, lineend="butt", columns=c(1:4), col="#99999922")),
           col=fpColors(box="black", lines="black", zero = "gray50"),
           xlab = "hatching success", xticks = c(0, 0.25, 0.50, 0.75, 1),
           boxsize=c(0.1,0.1))

dev.off()


##########################################################################
tabletext<-cbind(
  c("Mean", "Difference from the mean",  "GQF1", "GQF2", "P", "SHC", "B", "WC"),
  #c("", "", "GQF1", "GQF2", "P", "SHC", "B", "WC"),
  c("", "Sample size", 
    table(dataset3$treatment_factor)[order_fac]))

print("the untrasformed expectations are")
print(untrasformed)

png(file="model_outputs/jags/plots/forestplot_treatments_difference_from_the_mean_no_random.png")
forestplot(tabletext, 
           mean  = c(untrasformed[1,"mean"], NA, untrasformed[2:7,"mean"]), 
           lower = c(untrasformed[1,"lower"], NA, untrasformed[2:7,"lower"]),
           upper = c(untrasformed[1,"upper"], NA, untrasformed[2:7,"upper"]),
                    col=fpColors(box="black", lines="black", zero = "gray50"),
           xticks = c(0, 0.25, 0.50, 0.75, 1),xlab = "hatching success",
           boxsize=c(0.1,0.1))       

dev.off()


print("###########################################################################")
print("###########################################################################")
print("#######################################################")
print(paste("test equal ", test))
print(paste("cores = ", nc))
print(paste("thinning = ", nt))
print(paste("iterations = ", ni))
print(paste("burn-in = ", nb))
info <- sessionInfo()
print(info)
#######################################
print("the end")
sink(file = NULL)

