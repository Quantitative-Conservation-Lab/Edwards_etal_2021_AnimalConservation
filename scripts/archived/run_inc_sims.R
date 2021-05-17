##########################################
#This code runs the incubation_simulation script 
#Current draft 12/2/2020
#Written by sconver@uw.edu 

library(here)

source(here("scripts","incubation_simulation.r"))

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
#number of reps to run - maybe something like 10000 ultimately 
reps <- 100

#inc.method is either "WC" = whooping cranes, "ART" = artificial incubator, or "SHC" = sandhill cranes
#if you set inc.method to "SHC" you must also set the number of SHC pairs (no.SHC)
#all other control parameters in the function are defined above 

#Run multiple scenarios
out.1 <- inc.eval(inc.method = "WC", no.SHC = NA, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.2 <- inc.eval(inc.method = "ART", no.SHC = NA, inc.spaces = inc.spaces, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.3 <- inc.eval(inc.method = "SHC", no.SHC = 3, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.4 <- inc.eval(inc.method = "SHC", no.SHC = 4, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.5 <- inc.eval(inc.method = "SHC", no.SHC = 5, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.6 <- inc.eval(inc.method = "SHC", no.SHC = 6, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.7 <- inc.eval(inc.method = "SHC", no.SHC = 7, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.8 <- inc.eval(inc.method = "SHC", no.SHC = 8, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)
out.15 <- inc.eval(inc.method = "SHC", no.SHC = 15, inc.spaces = NA, pairs = pairs, recycle.days = recycle.days, max.eggs.pair = max.eggs.pair, lay.days = lay.days, season.end.threshold = season.end.threshold, reps = reps)

mean(out.1$eggs.laid)
quantile(out.1$eggs.laid,probs=c(0.025,0.975))
mean(out.1$successes)
quantile(out.1$successes,probs=c(0.025,0.975))
mean(out.1$success.rate)
quantile(out.1$success.rate,probs=c(0.025,0.975))

mean(out.2$eggs.laid)
quantile(out.2$eggs.laid,probs=c(0.025,0.975))
mean(out.2$successes)
quantile(out.2$successes,probs=c(0.025,0.975))
mean(out.2$success.rate)
quantile(out.2$success.rate,probs=c(0.025,0.975))

mean(out.3$eggs.laid)
quantile(out.3$eggs.laid,probs=c(0.025,0.975))
mean(out.3$successes)
quantile(out.3$successes,probs=c(0.025,0.975))
mean(out.3$success.rate)
quantile(out.3$success.rate,probs=c(0.025,0.975))

mean(out.4$eggs.laid)
quantile(out.4$eggs.laid,probs=c(0.025,0.975))
mean(out.4$successes)
quantile(out.4$successes,probs=c(0.025,0.975))
mean(out.4$success.rate)
quantile(out.4$success.rate,probs=c(0.025,0.975))

mean(out.5$eggs.laid)
quantile(out.5$eggs.laid,probs=c(0.025,0.975))
mean(out.5$successes)
quantile(out.5$successes,probs=c(0.025,0.975))
mean(out.5$success.rate)
quantile(out.5$success.rate,probs=c(0.025,0.975))

mean(out.6$eggs.laid)
quantile(out.6$eggs.laid,probs=c(0.025,0.975))
mean(out.6$successes)
quantile(out.6$successes,probs=c(0.025,0.975))
mean(out.6$success.rate)
quantile(out.6$success.rate,probs=c(0.025,0.975))

mean(out.7$eggs.laid)
quantile(out.7$eggs.laid,probs=c(0.025,0.975))
mean(out.7$successes)
quantile(out.7$successes,probs=c(0.025,0.975))
mean(out.7$success.rate)
quantile(out.7$success.rate,probs=c(0.025,0.975))

mean(out.8$eggs.laid)
quantile(out.8$eggs.laid,probs=c(0.025,0.975))
mean(out.8$successes)
quantile(out.8$successes,probs=c(0.025,0.975))
mean(out.8$success.rate)
quantile(out.8$success.rate,probs=c(0.025,0.975))

mean(out.15$eggs.laid)
quantile(out.15$eggs.laid,probs=c(0.025,0.975))
mean(out.15$successes)
quantile(out.15$successes,probs=c(0.025,0.975))
mean(out.15$success.rate)
quantile(out.15$success.rate,probs=c(0.025,0.975))
