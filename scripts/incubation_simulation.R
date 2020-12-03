##########################################
####EXPECTED OUTCOMES - CHICKS HATCHED#### 
##########################################
#This function simulates the incubation process when the incubation method is either 
#WC = the whooping cranes incubate their own eggs
#SHC = there are sandhill cranes to incubate WC eggs
#ART = there is an artificial incubator to incubate the eggs 
#Current draft 12/2/2020
#Written by sconver@uw.edu 

inc.eval <- function(inc.method = "WC", no.SHC = NA, inc.spaces = NA, pairs = NA, recycle.days = NA, max.eggs.pair = NA, lay.days = NA, season.end.threshold = NA, reps = NA){
  
#######################################################################################################################
#####HANNAH: need to bring in parametric uncertainty here - replace the short vector here with the full vector of posterior samples#####  
#######################################################################################################################
daily.surv.all <- matrix(0,nrow=3,ncol=2) 
daily.surv.all[1,] <- c(0.98,0.981) #pretend samples from WC daily incubation success
daily.surv.all[2,] <- c(0.97,0.971) #pretend samples from artificial incubator daily incubation success
daily.surv.all[3,] <- c(0.99,0.991) #pretend samples from SHC daily incubation success

#sample from the MCMC chains 
choose.sample <- sample(c(1:ncol(daily.surv.all)),1)
daily.surv <- daily.surv.all[1:3,choose.sample]

#management rule - WC, SHC, or artificial
#set incubator spaces available under SHC or artificial incubator
#select the necessary values for daily nest survival
if(inc.method == "WC"){spaces = 0}
if(inc.method == "SHC"){spaces <- no.SHC; daily.surv <- daily.surv[-c(2)] }
if(inc.method == "ART"){spaces = inc.spaces; daily.surv <- daily.surv[-c(3)] }

#calculate the total days of incubation season
inc.days <- lay.days + 30

#set up objects to hold all results 
incubate <- array(0,dim=c(reps,pairs,inc.days))
artificial <- array(0,dim=c(reps,max(1,spaces),inc.days))
eggs.laid <- matrix(0,nrow=reps,ncol=pairs)
successes.crane <- matrix(0,nrow=reps,ncol=pairs)  
successes.artif <- matrix(0,nrow=reps,ncol=pairs)  

#loop through reps
for(r in 1:reps){
  #loop through pairs
  for(p in 1:pairs){ 
    #initialize the pair at the beginning of the season
    egg.start <- 0
    #a pair can only lay if the total number of eggs laid is < max eggs per pair
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
        #determine whether to attempt to divert to non-WC incubator 
        if(inc.method == "WC"){divert.1 <- 0}else(divert.1 <- 1)
        #let the WCs incubate any eggs that are laid within X days of the end of the season
        if(egg.date > lay.days - season.end.threshold){divert.2 <- 0}else(divert.2 <- 1)
        #determine whether non-WC incubator space is available
        if(length((which(artificial[r,,egg.date]==0)))==0){divert.3 <- 0}else(divert.3 <- 1)
        divert <- divert.1*divert.2*divert.3
        #for eggs going to non-WC incubator 
        if(divert == 1){
          #find out what space they will be in
          avail.space <- min(which(artificial[r,,egg.date]==0))
          artificial[r,avail.space,egg.date] <- 1
          #daily incubation success 
          for(t in (egg.date+1):(egg.date+30)){
            artificial[r,avail.space,t] <- rbinom(1,1,daily.surv[2])*artificial[r,avail.space,t-1]
          }
          #how long did it incubate and was it successful - assign success to the pair that laid it 
          days.artificial <- sum(artificial[r,avail.space,egg.date:(egg.date+30)])
          successes.artif[r,p] <- successes.artif[r,p] + ifelse(days.artificial==31,1,0)
        }else if(divert == 0){  
          #incubation success 
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

eggs.laid.tot <- apply(eggs.laid,1,sum)
successes.artif.tot <- apply(successes.artif,1,sum)
successes.crane.tot <- apply(successes.crane,1,sum)
successes.tot <- successes.artif.tot + successes.crane.tot
success.rate.tot <- successes.tot/eggs.laid.tot

list(eggs.laid = eggs.laid.tot, successes.crane = successes.crane.tot, successes.artif = successes.artif.tot, successes = successes.tot, success.rate = success.rate.tot)

}