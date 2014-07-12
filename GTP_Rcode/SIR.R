#######
#SIR Filter for SV model
#
# x_t=phi x_{t-1}+N(0,sigma^2)
# y_t= N(0,exp\{2(gamma+ x_t)\})
#
# input y: observation
#      par: vector of parameters (phi,sigma,gamma)
#      prior: vector of prior parameters (prior mn, prior sd) for x_1
#      N: Number of particles
#      truth: true state -- used for plots
#      PLOT: produce plots of filtered density?
#
#  Output: matrix of approximate filtered mean and variance of X_t.
#
#######

SIRfilterSV=function(N,y,par,prior,truth=NULL,PLOT=T)
{
   n=length(y) ##number of observations
   d=1 ##dimension
   
   ##(1)simulate particles
  particles=matrix(0,nrow=N,ncol=d) ## particles: X_t^{(i)}
  for(j in 1:d) particles[,j]=rnorm(N,prior[j],prior[j+d])

   ##(2) calculate weights
   w=dnorm(y[1],0,exp(par[3]+particles[,1])) 
   
  ##(3) log of estimate of likelihood
  logl=log(mean(w))

   ##OPTIONAL -- STORAGE OF FEATURES OF FILTER -- MEAN AND VAR
   M.st=matrix(0,nrow=n,ncol=d)
   V.st=matrix(0,nrow=n,ncol=d)
   w=w/sum(w)
   for(j in 1:d){
     M.st[1,]=sum(w*particles[,j])
     V.st[1,]=sum(w*particles[,j]^2)-M.st[1,]^2
   }
   ##PLOT
   if(PLOT){
   par(mfrow=c(1,2))
   plot(range(particles[,1]),c(0,max(w)),xlab="X_t",ylab="Weight",type="n")
   for(j in 1:N) lines(c(particles[j,1],particles[j,1]),c(0,w[j]))
   if(length(truth>=1)) abline(v=truth[1],col=2)
    plot(density(particles[,1],weights=w),xlab="X_t",ylab="Density",main="")
 
    if(length(truth>=1)) abline(v=truth[1],col=2)
 }
     for(i in 2:n){##LOOP
  ##(1) Resample 
      index=sample(1:N,prob=w,size=N,rep=T)
      particles=particles[index,]
      if(d==1) particles=matrix(particles,ncol=1)
          
  ##(2) Propagate particles
    #x_t=phi x_{t-1}+N(0,sigma^2); phi=par[1]; sigma=par[2]
    particles[,1]=particles[,1]*par[1]+rnorm(N,0,par[2])
      
 
  ##(3) Weight
      w=dnorm(y[i],0,exp(par[3]+particles[,1])) 
  
  ##(4) update log of estimate of likelihood
   logl=logl+log(mean(w))

  ##OPTIONAL -- STORE FEATURES OF THE PARTICLE APPROXIMATION
    #normalise weights
    w=w/sum(w)
    for(j in 1:d){
     M.st[i,]=sum(w*particles[,j])
     V.st[i,]=sum(w*particles[,j]^2)-M.st[i,]^2
   }
   ##PLOT
      if(PLOT){
   par(mfrow=c(1,1))
     plot(density(particles[,1],weights=w),xlab="X_t",ylab="Density",main="") 
   if(length(truth>=i)) abline(v=truth[i],col=2)
    }
    }
  return(list(mean=M.st,var=V.st,l=logl)) ##output summaries 
 }
############################################################################
#######
#SIR Filter for SV model -- unknown parameters
#
# x_t=phi x_{t-1}+N(0,sigma^2)
# y_t= N(0,exp\{2(gamma+ x_t)\})
#
# STATE = (X_t,logit(phi))
#
# where if z=logit(phi) then phi= exp(z)/(1+exp(z))
#
# Parameters: tau=sigma/sqrt(1-rho^2); gamma
#
# input y: observation
#      prior: vector of prior parameters (prior mn, prior sd) for
#             (x_1, logit(phi))
#      par: vector of parameters (tau,gamma)
#      N: Number of particles
#      h: for KDE methods, we jitter by addining N(0,hSig^2) where
#         Sig^2 is the estimate of the posterior variance
#     LW: Use LiuWest shrinkage?
#
#  Output: matrix of approximate filtered mean and variance of STATE.
#
#######

SIRfilterSVpar=function(N,y,prior,par,h=0,LW=F,PLOT=T)
{
   n=length(y) ##number of observations
   d=2 ##dimension -- STATE=(X_t,logit(phi))
   
   ##(1)simulate particles
  particles=matrix(0,nrow=N,ncol=d) ## particles
  for(j in 1:d) particles[,j]=rnorm(N,prior[j],prior[j+d])

   ##(2) calculate weights y_1~N(0,exp(gamma+x_1)^2)
   w=dnorm(y[1],0,exp(par[2]+particles[,1])) 
   
  ##(3) log of estimate of likelihood
  logl=log(mean(w))

   ## STORAGE OF FEATURES OF FILTER -- MEAN AND VAR
   M.st=matrix(0,nrow=n,ncol=d)
   V.st=matrix(0,nrow=n,ncol=d)
   w=w/sum(w)
   for(j in 1:d){
     M.st[1,j]=sum(w*particles[,j])
     V.st[1,j]=sum(w*particles[,j]^2)-M.st[1,j]^2
   }
   ##PLOT
   if(PLOT){
   par(mfrow=c(1,1))
   plot(density(particles[,2],kernel="gaussian",weights=w,bw=0.001+sqrt(h*V.st[1,2])),xlab="logit(Phi)",ylab="Density",main="Time: 1")

   }
     for(i in 2:n){##LOOP
       ##(0) Estimate Mean and Variance --need if using KDE: i.e. h>0
       ## Mean  Not needed as we have stored these as M.st[i-1,]
       ## Var    Not needed as we have stored these as V.st[i-1,]
 
       
  ##(1) Resample 
      index=sample(1:N,prob=w,size=N,rep=T)
      particles=particles[index,]
      if(d==1) particles=matrix(particles,ncol=1)
          
  ##(2) Propagate particles
      #First update parameters
      if(h>0 && V.st[i-1,2]>0){##JITTER
        if(LW){#Liu-West -- shrinkage
          lambda=sqrt(1-h) ##
           particles[,2]=lambda*particles[,2]+(1-lambda)*M.st[i-1,2]
        }
        #Jitter --for both LW and non LW 
         particles[,2]=particles[,2]+rnorm(N,0,sqrt(h*V.st[i-1,2]))
          
      }
  #x_t=phi x_{t-1}+N(0,sigma^2);
   ## phi=exp(particles[,2])/(1+exp(particles[,2]); sigma=par[1]*sqrt(1-phi^2)
   phi=exp(particles[,2])/(1+exp(particles[,2]))
      sigma=par[1]*sqrt(1-phi^2)
     particles[,1]=particles[,1]*phi+rnorm(N,0,sigma)

  ##(3) Weight y_t~N(0,exp(gamma+x_t)^2) gamma=par[2]
      w=dnorm(y[i],0,exp(par[2]+particles[,1])) 
  
  ##(4) update log of estimate of likelihood
   logl=logl+log(mean(w))

  ##OPTIONAL -- STORE FEATURES OF THE PARTICLE APPROXIMATION
    #normalise weights
    w=w/sum(w)
    for(j in 1:d){
     M.st[i,j]=sum(w*particles[,j])
     V.st[i,j]=sum(w*particles[,j]^2)-M.st[i,j]^2
   }
   ##PLOT
      if(PLOT){
          plot(density(particles[,2],kernel="gaussian",weights=w,bw=0.001+sqrt(h*V.st[i,2])),xlab="logit(Phi)",ylab="Density",main=paste("Time",i))
        }
    }
  return(list(mean=M.st,var=V.st,l=logl)) ##output summaries 
 }








###########################################################################
#####
#SV simulation
#   n is number of time-points
#   par = (phi,sigma,gamma)
#
#   uses stationary distribution of state-process as the prior
#
###

SVsim=function(n,par){

  x=rep(0,n)
  x[1]=rnorm(1,0,par[2]/sqrt(1-par[1]^2)) ##stationary distribution
  for(i in 2:n) x[i]=x[i-1]*par[1]+rnorm(1,0,par[2])
  y=rnorm(n,0,exp(par[3]+x))
  return(list(x=x,y=y))
}


#######
#SIR Filter for SV model
#
# x_t=phi x_{t-1}+N(0,sigma^2)
# y_t= N(0,exp\{2(gamma+ x_t)\})
#
# input y: observation
#      par: vector of parameters (phi,sigma,gamma)
#      prior: vector of prior parameters (prior mn, prior sd) for x_1
#      N: Number of particles
#      truth: true state -- used for plots
#      PLOT: produce plots of filtered density?
#
#  Output: matrix of approximate filtered mean and variance of X_t.
#
#######

SIRfilterNoisyRV=function(N,y,par,prior,truth=NULL,PLOT=T)
{
   n=length(y) ##number of observations
   d=1 ##dimension

   ##(1)simulate particles
  particles=matrix(0,nrow=N,ncol=d) ## particles: X_t^{(i)}
  for(j in 1:d) particles[,j]=rnorm(N,prior[j],prior[j+d])

   ##(2) calculate weights
   w=dnorm(y[1],0,exp(par[3]+particles[,1]))

  ##(3) log of estimate of likelihood
  logl=log(mean(w))

   ##OPTIONAL -- STORAGE OF FEATURES OF FILTER -- MEAN AND VAR
   M.st=matrix(0,nrow=n,ncol=d)
   V.st=matrix(0,nrow=n,ncol=d)
   w=w/sum(w)
   for(j in 1:d){
     M.st[1,]=sum(w*particles[,j])
     V.st[1,]=sum(w*particles[,j]^2)-M.st[1,]^2
   }
   ##PLOT
   if(PLOT){
   par(mfrow=c(1,2))
   plot(range(particles[,1]),c(0,max(w)),xlab="X_t",ylab="Weight",type="n")
   for(j in 1:N) lines(c(particles[j,1],particles[j,1]),c(0,w[j]))
   if(length(truth>=1)) abline(v=truth[1],col=2)
    plot(density(particles[,1],weights=w),xlab="X_t",ylab="Density",main="")

    if(length(truth>=1)) abline(v=truth[1],col=2)
 }
     for(i in 2:n){##LOOP
  ##(1) Resample
      index=sample(1:N,prob=w,size=N,rep=T)
      particles=particles[index,]
      if(d==1) particles=matrix(particles,ncol=1)

  ##(2) Propagate particles
    #x_t=phi x_{t-1}+N(0,sigma^2); phi=par[1]; sigma=par[2]
    particles[,1]=particles[,1]*par[1]+rnorm(N,0,par[2])


  ##(3) Weight
      w=dnorm(y[i],0,exp(par[3]+particles[,1]))

  ##(4) update log of estimate of likelihood
   logl=logl+log(mean(w))

  ##OPTIONAL -- STORE FEATURES OF THE PARTICLE APPROXIMATION
    #normalise weights
    w=w/sum(w)
    for(j in 1:d){
     M.st[i,]=sum(w*particles[,j])
     V.st[i,]=sum(w*particles[,j]^2)-M.st[i,]^2
   }
   ##PLOT
      if(PLOT){
   par(mfrow=c(1,1))
     plot(density(particles[,1],weights=w),xlab="X_t",ylab="Density",main="")
   if(length(truth>=i)) abline(v=truth[i],col=2)
    }
    }
  return(list(mean=M.st,var=V.st,l=logl)) ##output summaries
 }

NoisyRVsim=function(n,sigmaRV,sigmaObs){

  x=rep(0,n)
  x[1]=rnorm(1,0,sigmaRV)
  for(i in 2:n) x[i]=x[i-1]
  y=rnorm(n,x,sigmaObs)
  return(list(x=x,y=y))
}
