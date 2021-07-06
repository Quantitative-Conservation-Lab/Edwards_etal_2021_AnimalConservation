##########################################
#This code runs the simulation described in Edwards et al.   
#First drafted 12/2/2020
#Written by Sarah Converse with edits by Hannah Edwards

#this library makes the data objects available to the simulation code 
library(here)

###SET UP CONTROL PARAMETERS 
#capacity of the artificial incubator
inc.spaces <- 25
#WC pairs
pairs <- 10
#minimum number of days before a pair can recycle after it has laid or has been incubating an egg
recycle.days <- 10
#the maximum number of eggs a pair can lay
max.eggs.pair <- 6
#length of the annual laying period, e.g., April 1 through June 1 
lay.days <- 60 
#if an egg is laid within this number of days before end of lay days, the WC incubate it themselves
season.end.threshold <- 10 
#number of reps to run - run something like 10000 or more  
reps <- 10000

#inc.method is either "WC" = whooping cranes, "ART" = artificial incubator, or "SHC" = sandhill cranes
#if you set inc.method to "SHC" you must also set the number of SHC pairs (no.SHC)
#all other control parameters in the function are defined above 

##########################################
####EXPECTED OUTCOMES - CHICKS HATCHED#### 
##########################################
#This function simulates the incubation process when the incubation method is: 
#WC = the whooping cranes incubate their own eggs
#SHC = there are sandhill cranes to incubate WC eggs
#ART = there is an artificial incubator to incubate the eggs 

# load the posteriors and format them for the analyses
# these are the MCMC samples from the posterior for the effect of treatment on hatching success
load(here("data","jagsfit_treatment_no_random.RData"))

#intercept 
int <- jagsfit.m_treatment$sims.list$int
#treatment effects 
treats <- jagsfit.m_treatment$sims.list$treats 
trt <- treats 
#back transform posterior samples to the probability scale for each of the 6 treatments 
for(i in 1:6){
  trt[,i] <- 1/(1+exp(-(int + treats[,i])))
}
#get a daily rate for each treatment 
trt.daily.all <- trt^(1/30)
#remove treatments that are not of interest (i.e., suboptimal incubators)
trt.daily <- trt.daily.all[,-c(2,5,6)]
colnames(trt.daily) <- c("WC","GQF1","SHC")

#put the predictions for daily incubation success in a data frame
#this is of dimension `MCMC samples` by 3 (number of treatments)
trt.daily <- as.data.frame(trt.daily)

#function to run the simulation 
inc.eval <- function(inc.method = NA, no.SHC = NA, inc.spaces = NA, pairs = NA, recycle.days = NA, max.eggs.pair = NA, lay.days = NA, season.end.threshold = NA, reps = NA){
  
  #this transposed object is WC, ART, and SHC values in rows, MCMC samples in columns 
  daily.surv.all <- t(trt.daily)
  
  #management rule - WC, SHC, or artificial
  #set incubator spaces available under SHC or artificial incubator
  #eliminate the other value for daily egg survival for either SHC or artificial incubator
  if(inc.method == "WC"){spaces = 0; daily.surv.all <- daily.surv.all} #here you will only use the first row of daily.surv.all but we need it to be a matrix 
  if(inc.method == "SHC"){spaces <- no.SHC; daily.surv.all <- daily.surv.all[-c(2),] }#get rid of artificial incubator values
  if(inc.method == "ART"){spaces = inc.spaces; daily.surv.all <- daily.surv.all[-c(3),] }#get rid of SHC values
  
  #calculate the max length of incubation season
  inc.days <- lay.days + 30
  
  #set up objects to hold all results 
  incubate <- array(0,dim=c(reps,pairs,inc.days))
  artificial <- array(0,dim=c(reps,max(1,spaces),inc.days))
  eggs.laid <- matrix(0,nrow=reps,ncol=pairs)
  successes.crane <- matrix(0,nrow=reps,ncol=pairs)  
  successes.artif <- matrix(0,nrow=reps,ncol=pairs)  
  
  #loop through reps
  for(r in 1:reps){
    #choose a sample from the MCMC chains for the daily egg survival values 
    choose.sample <- sample(c(1:ncol(daily.surv.all)),1)
    daily.surv <- daily.surv.all[,choose.sample]
    
    #loop through pairs
    for(p in 1:pairs){ 
      #initialize the pair at the beginning of the season
      egg.start <- 0
      #a pair can only lay if the total number of eggs it has laid is < max eggs per pair
      while (eggs.laid[r,p] < max.eggs.pair){
        #get date of egg laying - once a pair is able to lay, it will lay in the following 10 days
        egg.date <- egg.start + which(rmultinom(1,1,prob=rep(0.1,10))==1)
        #make sure a pair doesn't lay outside of the laying days
        if(egg.date > lay.days){break}
        else{
          #note egg laid by pair
          incubate[r,p,egg.date] <- 1 
          #increment eggs laid by 1 
          eggs.laid[r,p] <- eggs.laid[r,p] + 1
          #determine whether to divert to non-WC incubator 
          if(inc.method == "WC"){divert.1 <- 0}else(divert.1 <- 1)
          #let WCs incubate any eggs that are laid within 'season.end.threshold' days of the end of the laying season
          if(egg.date > lay.days - season.end.threshold){divert.2 <- 0}else(divert.2 <- 1)
          #determine whether non-WC incubator space is available
          if(length(which(artificial[r,,egg.date]==0))==0){divert.3 <- 0}else(divert.3 <- 1)
          #if all conditions are met, divert to non-WC incubator  
          divert <- divert.1*divert.2*divert.3
          #for eggs going to non-WC incubator 
          if(divert == 1){
            #find out what space they will be in (note that 'artificial' can be either SHC or artificial incubator)
            avail.space <- min(which(artificial[r,,egg.date]==0))
            #put egg in the space
            artificial[r,avail.space,egg.date] <- 1
            #daily incubation success 
            for(t in (egg.date+1):(egg.date+30)){
              artificial[r,avail.space,t] <- rbinom(1,1,daily.surv[2])*artificial[r,avail.space,t-1]
            }
            #how long did it incubate and was it successful - assign success to the pair that laid it 
            days.artificial <- sum(artificial[r,avail.space,egg.date:(egg.date+30)])
            successes.artif[r,p] <- successes.artif[r,p] + ifelse(days.artificial==31,1,0)
            #for eggs staying with the WC pair for incubation   
          }else if(divert == 0){  
            #daily incubation success 
            for(t in (egg.date+1):(egg.date+30)){
              incubate[r,p,t] <- rbinom(1,1,daily.surv[1])*incubate[r,p,t-1]
            }
            #how long did it incubate and was it successful?
            days.incubate <- sum(incubate[r,p,egg.date:(egg.date+30)])
            successes.crane[r,p] <- successes.crane[r,p] + ifelse(days.incubate==31,1,0) 
          }
        }
        #assume a break before a pair is able to lay again 
        egg.start <- max(which(incubate[r,p,]==1))+recycle.days
      }
    }
  }
  
  #summarize outputs 
  eggs.laid.total <- apply(eggs.laid,1,sum)
  successes.artif.total <- apply(successes.artif,1,sum)
  successes.crane.total <- apply(successes.crane,1,sum)
  successes.total <- successes.artif.total + successes.crane.total
  success.rate.total <- successes.total/eggs.laid.total
  
  #output from function 
  list(eggs.laid = eggs.laid.total, successes.crane = successes.crane.total, successes.artif = successes.artif.total, successes = successes.total, success.rate = success.rate.total)
  
}


#Run the function for multiple scenarios - these are the scenarios in the paper
#WC incubation 
out.WC <- inc.eval(inc.method = "WC", no.SHC = NA, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
#artificial incubation
out.ART <- inc.eval(inc.method = "ART", no.SHC = NA, inc.spaces = inc.spaces, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
#5 SHC pairs
out.SHC5 <- inc.eval(inc.method = "SHC", no.SHC = 5, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
#10 SHC pairs
out.SHC10 <- inc.eval(inc.method = "SHC", no.SHC = 10, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
#15 SHC pairs
out.SHC15 <- inc.eval(inc.method = "SHC", no.SHC = 15, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)

#summarize the results for these scenarios 
mean(out.WC$eggs.laid)
quantile(out.WC$eggs.laid,probs=c(0.025,0.975))
mean(out.WC$successes)
quantile(out.WC$successes,probs=c(0.025,0.975))
mean(out.WC$success.rate)
quantile(out.WC$success.rate,probs=c(0.025,0.975))

mean(out.ART$eggs.laid)
quantile(out.ART$eggs.laid,probs=c(0.025,0.975))
mean(out.ART$successes)
quantile(out.ART$successes,probs=c(0.025,0.975))
mean(out.ART$success.rate)
quantile(out.ART$success.rate,probs=c(0.025,0.975))

mean(out.SHC5$eggs.laid)
quantile(out.SHC5$eggs.laid,probs=c(0.025,0.975))
mean(out.SHC5$successes)
quantile(out.SHC5$successes,probs=c(0.025,0.975))
mean(out.SHC5$success.rate)
quantile(out.SHC5$success.rate,probs=c(0.025,0.975))

mean(out.SHC10$eggs.laid)
quantile(out.SHC10$eggs.laid,probs=c(0.025,0.975))
mean(out.SHC10$successes)
quantile(out.SHC10$successes,probs=c(0.025,0.975))
mean(out.SHC10$success.rate)
quantile(out.SHC10$success.rate,probs=c(0.025,0.975))

mean(out.SHC15$eggs.laid)
quantile(out.SHC15$eggs.laid,probs=c(0.025,0.975))
mean(out.SHC15$successes)
quantile(out.SHC15$successes,probs=c(0.025,0.975))
mean(out.SHC15$success.rate)
quantile(out.SHC15$success.rate,probs=c(0.025,0.975))
