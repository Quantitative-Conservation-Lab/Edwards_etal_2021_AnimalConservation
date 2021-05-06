
data <- sample(c(0:359),1000,replace=TRUE)
data <- matrix(data,nrow=100,ncol=10)

diff <- matrix(NA,nrow=nrow(data),ncol=ncol(data)-1)
for(i in 1:nrow(data)){
  for(j in 1:(ncol(data)-1)){
    diff[i,j] <-  min( abs(data[i,(j+1)] - data[i,j]), 360-abs(data[i,(j+1)] - data[i,j]))
  }
}
mean.diff <- apply(diff,1,mean)
var.diff <- apply(diff,1,var)