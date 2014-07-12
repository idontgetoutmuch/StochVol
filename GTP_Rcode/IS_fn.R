###
###Functions for Importance Sampling Example/Exercise
###
###

###These functions evaluate the (un-normalised density)
###for the Wright-Fisher with selection model of allele frequencies.

###1-d model
###
### f(x) \propto x^{alpha_1-1}*(1-x)^{alpha_2-1} * exp\{-sigma_1 x -sigma_2 x^2\}

###input is:
###      x -- vector of points at which to evaluate the density function
###      alpha -- vector of (alpha_1,alpha_2)
###      sigma -- vector of (sigma_1,sigma_2)
###      log -- return the log density?
dWF=function(x,alpha,sigma,log=F)
  {
    if(min(alpha)<=0 || length(alpha)<2){
      cat("\n Error in alpha vector\n")
      return()
    }
    if(length(sigma)<2){
      cat("\n Error in sigma vector\n")
      return()
    }
    ##calculate log density
    lf=(alpha[1]-1)*log(x)+(alpha[2]-1)*log(1-x)-sigma[1]*x-sigma[2]*x^2
    if(log){
     return(lf) 
    }else{
      return(exp(lf))
    }
  }
