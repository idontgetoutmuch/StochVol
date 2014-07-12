#######################################################################
#
# PseudoMarginal/Particle MCMC example
#
source("SIR.R") ##Load in function SIRfilterSV; and SVsim

##simulate data -- remember observations are data$y

data=SVsim(100,c(0.9, sqrt(1-0.9^2), 1)) ##parameters are (phi,sigma,gamma)

###PRIOR on logit(phi); log(tau); and gamma are normal; means and sd given in vector
##  p.m and p.sd respectively

p.m=c(3,-1,1);p.sd=c(1,1,1)

### random walk sds given by prop.sd
prop.sd=c(0.5,0.1,0.1)

M=5000 #number of MCMC iterations
N=100 #number of particles

###Storage
state.st=matrix(0,nrow=M,ncol=4)
acc.n=rep(0,3) ##accept prob for each update

##Initialise
state.cur=rep(0,4)  ##State is logit(phi),log(tau),gamma, log(Z)
state.cur[1:3]=rnorm(3,p.m,p.sd) 
###Initially just assume phi is unknown 
state.cur[3]=1 ##gamma=1
state.cur[2]=0 ##log-tau=0
phi=exp(state.cur[1])/(1+exp(state.cur[1]))
sigma=exp(state.cur[2])*sqrt(1-phi^2)
par=c(phi,sigma,state.cur[3] ) ##parameters needed for input to PF
prior=c(0,exp(state.cur[2])) ##Prior for X_1 is N(0,tau^2)
state.cur[4]=SIRfilterSV(N,data$y,par,prior,PLOT=F)$l 
for(i in 1:M){
	##Update logit(phi)
	lp.prop=rnorm(1,state.cur[1],prop.sd[1])
	##Generate estimate of log(p(data|parameters))
	phi=exp(lp.prop)/(1+exp(lp.prop))
	sigma=exp(state.cur[2])*sqrt(1-phi^2)
	par=c(phi,sigma,state.cur[3] ) ##parameters needed for input to PF
	prior=c(0,exp(state.cur[2])) ##Prior for X_1 is N(0,tau^2)
	Z.prop=SIRfilterSV(N,data$y,par,prior,PLOT=F)$l 
	##Acceptance probability: again proposal ratio cancels
	acc.p=exp(Z.prop-state.cur[4]+ sum(dnorm(c(lp.prop,state.cur[2:3]),p.m,p.sd,log=T))-sum(dnorm(state.cur[1:3],p.m,p.sd,log=T)) )
	if(runif(1)<acc.p){
		state.cur[1]=lp.prop
		state.cur[4]=Z.prop
		acc.n[1]=acc.n[1]+1
	}
	
	##update log(tau)
	
	##update gamma
		
    state.st[i,]=state.cur
}
