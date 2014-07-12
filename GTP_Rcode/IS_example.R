####
#### Example for Importance Sampling
####

#### Input functions for calculating target density
source("IS_fn.R")

#Parameters
alpha=c(3,0.5)
sigma=c(5,1)

#Plot target distribution
x=seq(0.001,0.999,by=0.001)
par(mfrow=c(1,2))
logf=dWF(x,alpha,sigma,log=T)
plot(x,logf,type="l",lwd=2,xlab="x",ylab="log(f)")

plot(x,exp(logf-max(logf)),type="l",lwd=2,xlab="x",ylab="f")

##Importance Sampling: from uniform proposal

N=10000 ##number of samples
x=runif(N) ##generate from proposal
##calculate weights
w=dWF(x,alpha,sigma)/dunif(x)

##Example usess
logChat=log(sum(w)/N)##estimate of log of normalising constant
w=w/sum(w) ##normalised weights

#######################Alternative IS approach -- using log-weights to avoid numerical instability
logw=dWF(x,alpha,sigma,log=T)-dunif(x,log=T) 
maxw=max(logw)
w=exp(logw-maxw)
sumw=sum(w)
w=w/sumw ##normalised weights
logChat=maxw+log(sumw/N) ##estimate of log of normalising constant

####ESS calculation -- measure of accuracy
ESS=sum(w)^2/sum(w^2) ##Effective Sample Size
### Plot weights
hist(w)

####Use of weighted sample
sum(w*x) ##estimate of mean frequency
sum(w*(x>0.99)) ##estimate of probability frequency > 0.99

###Resampling to get approximate (unweighted) sample
M=1000 ##sample size
index=sample(1:N,prob=w,size=M,rep=T) ##index of particles to keep
xsam=x[index]
hist(xsam,breaks=seq(0,1,by=0.01))



