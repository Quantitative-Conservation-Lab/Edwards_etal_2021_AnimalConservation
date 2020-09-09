#read in data to clean up some stuff and re-output 
eggs <- read.csv("Egg_summary_survival_predictors.csv")

#get a simpler Egg ID to use  
eggs$Egg <- as.numeric(as.factor(egg.data.all$Egg.ID))

#get a simpler Logger ID to use 
eggs$Logger <- as.numeric(as.factor(egg.data.all$Logger.ID))

#get outcome of nest as 0/1 = "fail" or "failed"/"hatch" or "hatched" 
eggs$outcome <- rep(0,nrow(egg.data.all))
eggs$outcome[which(egg.data.all$Fate == "Failed")] <- 0 
eggs$outcome[which(egg.data.all$Fate == "Hatched")] <- 1 

keep <- c("Egg","Logger","Treatment","Facility","Year","Pair.ID","Paired.date","Unpaired.date","Last.alive.date","Hatch.fail.date","outcome")
eggs.new <- eggs[,keep]

write.csv(eggs.new,"eggs.new.csv",row.names=FALSE)

#########################################################################################
#I read the data out and read it back in as this might be the cleaner dataset to publish
eggs.new <- read.csv("eggs.new.csv")

eggs.new$Facility <- as.numeric(as.factor(eggs.new$Facility))
eggs.new$Year <- as.numeric(as.factor(eggs.new$Year))
eggs.new$Treatment <- as.numeric(as.factor(eggs.new$Treatment))

#start all nests at time 1  
eggs.new$start <- eggs.new$Paired.date - (eggs.new$Paired.date) + 1
eggs.new$unpaired <- eggs.new$Unpaired.date - (eggs.new$Paired.date) + 1
eggs.new$last.alive <- eggs.new$Last.alive.date - (eggs.new$Paired.date) + 1
eggs.new$term <- eggs.new$Hatch.fail.date - (eggs.new$Paired.date) + 1

#build encounter histories 
enc.hist <- matrix(nrow=nrow(eggs.new),ncol=max(as.numeric(eggs.new$term)))
enc.hist[] <- NA
first <- last <- rep(NA,nrow(eggs.new))
for(i in 1:nrow(enc.hist)){
  enc.hist[i,((eggs.new$start[i]):(eggs.new$last.alive[i]))] <- 1    
  if(eggs.new$outcome[i] == 0){enc.hist[i,(eggs.new$term[i])] <- 0}
  first[i] <- as.integer(eggs.new$start[i])
  last[i] <- as.integer(eggs.new$term[i])
}
nind <- nrow(enc.hist)










#data
dataset <- list(nind=nind,first=first,last=last,enc.hist=enc.hist,pair.ID=eggs.new$Pair.ID,npair=length(unique(eggs.new$Pair.ID)),facility=eggs.new$Facility)

#initial values
inits <- function(){
  list (beta.fac = c(runif(2),NA))
}

#parameters
parameters <- c("beta.int","beta.fac")

######################### CREATE MODEL FILE #########################

cat("
model {
    
#####LIKELIHOOD
for(i in 1:nind){  
  for(t in first[i]:(last[i]-1)){
    enc.hist[i,(t+1)] ~ dbern(S[i,t]*enc.hist[i,t])   
    logit(S[i,t]) <- beta.int + beta.fac[facility[i]] + alpha.Pr[pair.ID[i]]
  }
}
    

#####PARAMETER PRIORS
    
#PARAMETER PRIORS 
beta.int ~ dnorm(0,0.001)

for(f in 1:2){
  beta.fac[f] ~ dnorm(0,0.001)
}
beta.fac[3] <- 0

#RANDOM EFFECTS AND ATTEMPT
for(i in 1:npair){
  alpha.Pr[i] ~ dnorm(0,tau.pair)
}
tau.pair <- pow(sigma.pair,-2)
sigma.pair ~ dunif(0,25) 

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

