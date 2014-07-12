#################
# Example of SIR filter for Noisy Observation of a Single Random Variable Model
#
# x_t = x_{t-1}
# y_t = x_t + N(0,sigma^2)
#
#
################
source("SIR.R")

##First simulate data:
n=100 ##number of data points

##parameters

sigma=1.0

## par=c(phi,sigma,gamma) ##parameter vector

data=NoisyRVsim(n,sigma,sigma)
y=data$y ##observations
truth=data$x ##true state

prior=c(0,sigma/sqrt(1-phi^2)) ##prior for SIR filter
N=100 ##number of particles
PLOT=T
par(ask=T) ## this will produce a prompt after each plot
## par(ask=F)## use to remove the prompt after each plot
out=SIRfilterSV(N,y,par,prior,truth,PLOT=T)

##Can plot the output
plot(c(0,n),c(-3,3),type="n",xlab="t",ylab="X_t")
lines(1:n,truth,lwd=2) ##true state
lines(1:n,out$m,lwd=2,col=2) ##filtered estimate of state
lines(1:n,out$m+2*sqrt(out$v),lwd=1,col=2,lty=2) ##upper error band
lines(1:n,out$m-2*sqrt(out$v),lwd=1,col=2,lty=2) ##lower error band

##One way to test accuracy is to do multiple runs -- and look at the variability
##across runs (for the SAME data set)

##Use large N to get "accurate" estimate of posterior
out.accurate=SIRfilterSV(1e4,y,par,prior,truth,PLOT=F)
out=SIRfilterSV(N,y,par,prior,truth,PLOT=F)
plot(c(0,n),c(-2,2),type="n",xlab="t",ylab="X_t")
lines(1:n,truth,lwd=2) ##true state
##Accurate estimate in Blue
lines(1:n,out.accurate$m,lwd=2,col=4) ##filtered estimate of state
lines(1:n,out.accurate$m+2*sqrt(out$v),lwd=1,col=4,lty=4) ##upper error band
lines(1:n,out.accurate$m-2*sqrt(out$v),lwd=1,col=4,lty=4) ##lower error band
##Standard estimate in Red
lines(1:n,out$m,lwd=2,col=2) ##filtered estimate of state
lines(1:n,out$m+2*sqrt(out$v),lwd=1,col=2,lty=2) ##upper error band
lines(1:n,out$m-2*sqrt(out$v),lwd=1,col=2,lty=2) ##lower error band

##############################################################################
#
#SV model with unknown parameters
#
##############################################################################
n=2000
phi=0.9
sigma=sqrt(1-phi^2) ##this choice fixes the marginal variance of X_t at stationarity to 1
gamma=1

data=SVsim(n,c(phi,sigma,gamma))
y=data$y ##observations
truth=data$x ##true state

prior=c(0,1.5, 1,1) ##prior means and sds for X_1,logit(phi)
par=c(1,gamma) ## parameter for tau=sigma/sqrt(1-phi^2), and gamma.
N=1000 ##number of particles
PLOT=T 
#par(ask=T) ## this will produce a prompt after each plot
par(ask=F)## use to remove the prompt after each plot
out=SIRfilterSVpar(N,y,prior,par,h=0,LW=F,PLOT=PLOT)
outJ=SIRfilterSVpar(N,y,prior,par,h=0.01,LW=F,PLOT=PLOT)
outLW=SIRfilterSVpar(N,y,prior,par,h=0.1,LW=T,PLOT=PLOT)

###PLOT
par(mfrow=c(1,1))
j=2
plot(c(0,n),c(prior[j]-3*prior[j+2],prior[j]+3*prior[j+2]),type="n",xlab="Time",ylab="logit(phi)")
lines(1:n,rep(log(phi/(1-phi)),n),lwd=2)
lines(1:n,out$m[,j],col=2) ##Red
lines(1:n,outJ$m[,j],col=3) ##Jitter - Green
lines(1:n,outLW$m[,j],col=4) ##Liu-West Shrinkage - Blue
lines(1:n,out$m[,j]+2*sqrt(out$v[,j]),lty=2,col=2)
lines(1:n,out$m[,j]-2*sqrt(out$v[,j]),lty=2,col=2)

lines(1:n,outJ$m[,j]+2*sqrt(outJ$v[,j]),lty=2,col=3)
lines(1:n,outJ$m[,j]-2*sqrt(outJ$v[,j]),lty=2,col=3)

lines(1:n,outLW$m[,j]+2*sqrt(outLW$v[,j]),lty=2,col=4)
lines(1:n,outLW$m[,j]-2*sqrt(outLW$v[,j]),lty=2,col=4)
###LOOK AT POSTERIOR SDS
plot.ts(sqrt(outJ$v[,2]),col=3,ylim=c(0,sqrt(max(outJ$v[,2]))))
lines.ts(sqrt(out$v[,2]),col=2)
lines.ts(sqrt(outLW$v[,2]),col=4)
