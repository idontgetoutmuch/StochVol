###############################
#Example ABC:
#
# Inference for alpha-stable
###############################

library(stabledist) ##load library

###simulate data to analyse
alpha=1.2 ##shape
beta=0.1 ##skew
gamma=1 ##scale
delta=2 ##location

n=100

y=rstable(n,alpha,beta,gamma,delta) ##data to analyse

######
# Assume delta unknown; with N(0,p.sd^2) prior

p.sd=5


###define summary statistics
S=function(y){
 return(mean(y))
}

Sobs=S(y) ##observed summary statistics

d=length(Sobs) ##dimension of summary statistics
A=matrix(0,nrow=d,ncol=d); diag(A)=1 ##matrix for distance -- equal weight for all components

############ABC REJECTION SAMPLER
#
# Rather than pre-specifying threshold for distance: store 
#    (delta,dist) pairs; and post-process these for different thresholds
#

N=10000 ##number of samples
ABC.out=matrix(0,nrow=N,ncol=2) ##matrix to store (delta,dist) pairs

for(i in 1:N){ ##LOOP 
  delta.prop=rnorm(1,0,p.sd)
  ysim=rstable(n,alpha,beta,gamma,delta.prop) ##simulated data
  Ssim=S(ysim) #summaries
  dist=matrix(Sobs-Ssim,nrow=1) %*% A %*% matrix(Sobs-Ssim,ncol=1)
  ABC.out[i,]=c(delta.prop,dist)
}

h=quantile(ABC.out[,2],0.02) ##threshold for ABC
ABCsam=ABC.out[ABC.out[,2]<h,1] ##sample of delta values
plot(density(ABCsam))
cat("ABC posterior mean",mean(ABCsam),"\tABC posterior variance",var(ABCsam),"\n")


########################ABC MCMC ALGORITHM
#
#
# Need to pre-specify threshold; and random walk variance


N=10000 ##number of samples
ABC.out=rep(0,N)
h=0.1
prop.sd=0.1

##initialise 
delta.cur=median(y)

for(i in 1:N){
  delta.prop=rnorm(1,delta.cur,prop.sd)
  ysim=rstable(n,alpha,beta,gamma,delta.prop) ##simulated data
  Ssim=S(ysim) #summaries
  dist=matrix(Sobs-Ssim,nrow=1) %*% A %*% matrix(Sobs-Ssim,ncol=1)
  if(dist<h){
    acc.p=dnorm(delta.prop,0,p.sd)/dnorm(delta.cur,0,p.sd)
    if(runif(1)<acc.p){
      delta.cur=delta.prop
    }
  }
  ABC.out[i]=delta.cur
}
plot.ts(ABC.out) ##trace plot

ABCsam=ABC.out[-(1:(N/10)) ##remove burn-in
plot(density(ABC.sam))
cat("ABC posterior mean",mean(ABCsam),"\tABC posterior variance",var(ABCsam),"\n")


  
