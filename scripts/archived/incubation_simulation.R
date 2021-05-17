##########################################
####EXPECTED OUTCOMES - CHICKS HATCHED#### 
##########################################
#This function simulates the incubation process when the incubation method is: 
#WC = the whooping cranes incubate their own eggs
#SHC = there are sandhill cranes to incubate WC eggs
#ART = there is an artificial incubator to incubate the eggs 

#Current draft 12/2/2020
#Written by sconver@uw.edu 
#Updated by HannahE@calgaryzoo.com on 12/4/2020

#library(stringr)
library(here)

# load the posteriors and format them for the analyses
load(here("data","jagsfit_treatment_no_random.RData"))
#rn <- row.names(jagsfit.m_treatment$summary)
#para <- rn[str_detect(rn, "log.like", negate = T)]

#expit <- function (x) round(exp(x)/(1 + exp(x)), 2)

int <- jagsfit.m_treatment$sims.list$int
treats <- jagsfit.m_treatment$sims.list$treats 
trt <- treats 
for(i in 1:6){
  trt[,i] <- 1/(1+exp(-(int + treats[,i])))
}
trt.daily.all <- trt^(1/30)
trt.daily <- trt.daily.all[,-c(2,5,6)]
colnames(trt.daily) <- c("WC","GQF1","SHC")

trt.daily <- as.data.frame(trt.daily)

#test <- FALSE
##test <- TRUE
#if(test==TRUE){
#  nc <- 4
#  nt <- 1#5 #thinning
#  ni <- 500#50000 #iteration
#  nb <- 50#ni - 40000 #burn-in
#  na <- 500 # adaptation
#} else{
#  set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion")
#  nc <- 4
#  nt <- 5 #thinning
#  ni <- 300000 #iteration
#  nb <- ni - 80000 #burn-in
#  na <- 100000 # adaptation
#  (ni-nb)/5
#  if(nb<=0) stop("burn in is negative")
#}

#samples_per_chain <- ((ni-nb)/nt)

#posteriors <- data.frame(GQF1 = numeric(nc*samples_per_chain), WC = numeric(nc*samples_per_chain), SHC = numeric(nc*samples_per_chain))

#for(i in 1:nc){
#  posteriors[(1+(samples_per_chain*(i-1))):(samples_per_chain*i), "GQF1"] <- expit(jagsfit.m_treatment$samples[, "int"][[i]] +
#                                                                                     jagsfit.m_treatment$samples[, "treats[1]"][[i]])
#  posteriors[(1+(samples_per_chain*(i-1))):(samples_per_chain*i), "WC"] <- expit(jagsfit.m_treatment$samples[, "int"][[i]] + 
#                                                                                   jagsfit.m_treatment$samples[, "treats[6]"][[i]])
#  posteriors[(1+(samples_per_chain*(i-1))):(samples_per_chain*i), "SHC"] <- expit(jagsfit.m_treatment$samples[, "int"][[i]] + 
#                                                                                    jagsfit.m_treatment$samples[, "treats[4]"][[i]])
#}

inc.eval <- function(inc.method = NA, no.SHC = NA, inc.spaces = NA, pairs = NA, recycle.days = NA, max.eggs.pair = NA, lay.days = NA, season.end.threshold = NA, reps = NA){

#this object is WC, ART, and SHC values in rows, MCMC samples in columns 
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

eggs.laid.total <- apply(eggs.laid,1,sum)
successes.artif.total <- apply(successes.artif,1,sum)
successes.crane.total <- apply(successes.crane,1,sum)
successes.total <- successes.artif.total + successes.crane.total
success.rate.total <- successes.total/eggs.laid.total

list(eggs.laid = eggs.laid.total, successes.crane = successes.crane.total, successes.artif = successes.artif.total, successes = successes.total, success.rate = success.rate.total)

}
