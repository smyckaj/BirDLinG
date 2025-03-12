###############################
#This script is for handling the probabilistic aspect of predicting number of offspring species from a tip estimate of diversification rates. The core function is generating generalized Rsquare vaĺue of the tip rates when we know real number of offspring species, e.g. from a fossil phylogeny. Lower-level function perform a respective regression models and code a probability density of distribution of number of offspring species after time t given speciation and extinction.


###############################
#probability density (actually mass) function from Luke Harmons book https://lukejharmon.github.io/pcm/pd/phylogeneticComparativeMethods.pdf, but with typo corrected according to Raup 1985 Paleobiology https://www.jstor.org/stable/2400422
dbdnnv=function(n,lamb,mu,t,extantonly=F, log=F,jittervar=0.000001){
  
  if(n!=as.integer(n)){
    print("some of n values are not integers, function is untested in this context")
  }
  
  #calculate transformed parameters alpha and beta
  al=(mu*((exp((lamb-mu)*t))-1))/((lamb*exp((lamb-mu)*t))-mu)
  be=(lamb*((exp((lamb-mu)*t))-1))/((lamb*exp((lamb-mu)*t))-mu)
  
  #the sppecification for extantonly==T is experimental 
  if(extantonly){
    #n=0 cannot be observed in extant phylogeny, higher values are normalized to satisfy summing across n to 1
    if(n==0){
      d=0
    }else{
      d=((1-al)*(1-be)*(be^(n-1)))/(1-al)
    }
    
  } else{
    #branching equation exactly as in the Harmons book  
    if(n==0){
      d=al
    }else{
      d=(1-al)*(1-be)*(be^(n-1))
    }
  }


  #some values throw NaN, for instance lamb==0 or lamb==mu, ad truncated normal jitter to them and recursively rerun
  if(is.nan(d)){
    lamb=abs(rnorm(1,lamb,jittervar))
    mu=abs(rnorm(1,mu,jittervar))
    print(paste("lambda and mu had to be jittered to", lamb, mu, sep=" "))
    d=dbdnnv(n=n, lamb=lamb, mu=mu, t=t,extantonly=extantonly, log=F, jittervar=jittervar)
  }

  #negative lamb, mu and t have probability mass of 0
  if(lamb<0 | mu<0 | t<0 | n<0){
    d=0
  }
  
  if(d<0 | d>1){
    print(paste("something is wrong, detected probability is", d, sep=" "))
    d=0
  }
  
  #return log values in needed
  if(log){
    d=log(d)
  }
  
  return(d)
}


#' Probability density function of number of species from a birth-death model
#'
#' The function is taken from Luke Harmons book https://lukejharmon.github.io/pcm/pd/phylogeneticComparativeMethods.pdf, but with typo corrected according to Raup 1985 Paleobiology https://www.jstor.org/stable/2400422. All parameters are fully vectorized.
#' @param n Number of observed species
#' @param lamb Speciation rate of a birth-death process
#' @param mu Extinction rate of a birth-death process
#' @param t Time period after which species are observed
#' @param extantonly Condition the distribution for non extinction, default is F
#' @param log Return log-probability, default is F
#' @param jittervar Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation
#' @return Probability of number of species for given alpha and mu
#' @examples 
#' 
#' #test for writing it well, with mu=1, t=1, the distribution should equal to geometric for n-1
#' dbdn(1,2,0,1)
#' dbdn(1,3,3,1)
#' dbdn(0,0,3,1)
#' dbdn(2,0,3,1)
#' dbdn(2,-1,3,1)
#' dbdn(2,0,3,5)
#' dbdn(2,5,5,5)
#' dbdn(1,2,0,1)
#' dgeom(1-1,exp(-2))
#' 
#' dbdn(10,2,0,1)
#' dgeom(10-1,exp(-2))
#' 
#' dbdn(1,12,0,1)
#' dgeom(1-1,exp(-12))
#' 
#' #tests for vectorisation
#' dbdn(1:5,1:5,0,1,log=c(F,T,T,T,T))
#' dbdn(2,2,0,1,F)
#' 
#' #negative values
#' dbdn(-2,12,0,1)
#' dbdn(1,-12,0,1)
#' dbdn(1,12,-10,1)
#' dbdn(0.1,12,-10,1)
#'
#' #the distribution sums to 1
#' sum(dbdn(0:100000,1,2,1))
#' sum(dbdn(0:100000,1,0,1))
#' sum(dbdn(0:100000,1,0.5,10))
#'
#' #extantonly testing
#' dbdn(1,1,2,1,extantonly = F)
#' dbdn(1,1,2,1, extantonly = T)
#' dbdn(0,1,2,1, extantonly = F)
#' dbdn(0,1,2,1, extantonly = T)
#' sum(dbdn(1:100000,1,2,1,extantonly = T))
#' sum(dbdn(1:100000,1,2,1,extantonly = F))
#' sum(dbdn(0,1,2,1,extantonly = F))
#' @export
dbdn=function(n,lamb,mu,t,extantonly=F,log=F, jittervar=0.000001){
  #it is technically possibly to vectorize also along extantonly, log and jittervar, why not:-)
  d=mapply(dbdnnv,n = n, lamb = lamb, mu = mu, t=t,extantonly=extantonly, log = log, jittervar=jittervar, SIMPLIFY = T)
  return(d)
}



###############################
#regression log likelihood, params are in order speciation intercept, speciation slope, extinction intercept, extinction slope
#data must be a data frame with columns named n, tipspec, tipext, and with rows representing tips species richness after time t, and tip estimates of speciation and extinction 

#' Log-likelihood of linear transformation of tip-estimated lambda and mu using dbdn distribution
#'
#' Calculates a log-likelihood as a sum of individual log probabilities for individual data points with real lambdas and mus being the linenar functions of observed tip lambdas and mus 
#' @param params A vector of length 4, with regression parameters ordered c(intercept_lambda, slope_lambda, intercept_mu, slope_mu)
#' @param data A data frame with columns n, tipspec and tipext, referring to number of species resulting from every tip after time t, estimated tip speciation and tip extinction rates 
#' @param t Time period after which species are observed
#' @param extantonly Condition the distribution for non extinction, default is F
#' @param jittervar Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation
#' @return Log-likelihood value
#' @examples 
#' #dataset without extinct tips
#' testdata=data.frame(n=1:10,tipspec=1:10,tipext=rep(0,10))
#' bdnregloglik(c(0,1,0,0),testdata,1)
#' bdnregloglik(c(0,1,1,0),testdata,1)
#' bdnregloglik(c(0,1,0,1),testdata,1)
#' bdnregloglik(c(0,1,1,1),testdata,1)
#' 
#' #dataset with extinct tips
#' testdata0=data.frame(n=0:10,tipspec=0:10,tipext=rep(0,11))
#' bdnregloglik(c(0,1,1,0),testdata0,1, extantonly = T)
#' bdnregloglik(c(0,1,1,0),testdata0,1, extantonly = F)
#' bdnregloglik(c(0,1,1,0),testdata,1, extantonly = T)
#' bdnregloglik(c(0,1,1,0),testdata,1, extantonly = F)
#' @export
bdnregloglik=function(params,data,t,extantonly=F, jittervar=0.000001){
  #equations for expected speciation and extinction rates based on tip rate metrics
  lambhat=params[1]+params[2]*data$tipspec
  muhat=params[3]+params[4]*data$tipext
  
  #likelihood parsing 
  llh=sum(dbdn(n = data$n, lamb = lambhat, mu = muhat, t=t, extantonly=extantonly, log = T, jittervar = jittervar))
  return(llh)
}



###############################
#saturated likelihood
#first explore behaviour of dbdn for boundary conditions where it shows NaN
#put likelihood 1 for n=1
# dbdn(1,log(1),0,1)
# dbdn(1,log(1)+0.1,0.1,1)
# dbdn(1,0.000000000000001+0.1,0.1,1)
# dbdn(1,0.000000000000001+0.2,0.2,1)
# dbdn(1,0.0000000000000000000000000000000001,0,1)
#
#put likelihood 1 for n=0
# dbdn(0,0.0000000000000000000000000000000001,2,1)
# dbdn(0,0.0000000000000000000000000000000001,200,1)
# dbdn(0,0.1,2,1)
#
#look at which parameter combinations are most probable for n>1, it is the ones where mu=0 and lamb=log(n), see deterministic equation in Harmons book
# dbdn(2,log(2),0,1)
# dbdn(2,log(2)+0.1,0.1,1)
# dbdn(2,log(2),0.001,1)
# dbdn(50,log(50),0,1)
# dbdn(50,log(50)-0.001,0,1)

#this function calculates a saturated probability mass for every count of species, i.e. maximum of probability mass function
#it also works around jittering for NaN values, by inputting asymptotic values in branching
dbdnsatnv=function(n,mu=0,t=1,extantonly=F,log=F, jittervar=0.000001,
                   extinctmarginal=F, startval=0.1,
                   optimpar_plex = list(reltolpar=1e-04, 
                                        reltolval=1e-05, 
                                        abstolpar=1e-07,
                                        maxiter=10000, 
                                        num_cycles=1)){
  
  if(extinctmarginal){
    # if we pick marginalisation over extinction, the function will find optimal lambda for given extinction numerically
    fn=function(lambda){
      dbdn(n=n, lamb=lambda, mu=mu, t=t, extantonly=extantonly,jittervar=jittervar)
    }
    
    d=optimglwrap(optimmethod = "subplex",
                  optimpar_plex = optimpar_plex, 
                  fc = fn,
                  trparsopt = startval)$fvalues
    
    
  } else {
    # if we let extinction free, the optimum is this branching, with mu=0 for all except n=0
    if(n==0){
      if(extantonly){
        d=0 #cannot happen
      } else{
        d=1 #for mu->Inf
      }
    }else if(n==1){
      d=1 #for mu=0 and lamb=0
    }else{
      d=dbdn(n,lamb=log(n)/t, mu=0, t=t,extantonly=extantonly,log=F,jittervar=jittervar) #deterministic equation
    }
  }

  #log if needed
  if(log){
    d=log(d)
  }
  return(d)
}


#' Saturated likelihood of lambda and mu values for individual tip counts
#'
#' Calculates a saturated likelihoods for individual species counts, i.e. the maximum values of dbdn function across lambda and mu. All parameters are fully vectorized.
#' @param n Number of observed species
#' @param mu Assumed extinction rate of a birth-death process in case we use extinction marginalisation
#' @param t Time period after which species are observed
#' @param extantonly Condition the distribution for non extinction, default is F
#' @param log Return log-probability, default is F
#' @param jittervar Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation
#' @param extinctmarginal Do we marginalize the distribution over the assumed mu values, default is F
#' @param startval Starting values of the subplex optimisation of lambda in case we marginalize over assumed mu
#' @param optimpar_plex The control list of subplex optimisation of lambda in case we marginalize over assumed mu
#' @return Probability of number of species for maximum lambda and mu
#' @examples 
#' #test if it matches nonsaturated lhs
#' dbdnsat(1:10)
#' dbdn(1:10, log(1:10),0, t=1)
#'
#' #test marginalized version
#' dbdnsatnv(5,extinctmarginal = F)
#' dbdnsatnv(5,extinctmarginal = T)
#' dbdnsatnv(5,mu=5,extinctmarginal = T)
#' curve(dbdn(n=5,lamb=x,mu=0,t=1),0,10)
#' curve(dbdn(n=5,lamb=x,mu=5,t=1),add=T)
#' dbdnsat(5,extinctmarginal = F)
#' dbdnsat(0:5,extinctmarginal = F)
#' dbdnsat(0:5,extinctmarginal = T)
#' dbdnsat(0:5,mu=c(10,0,0,0,0,0),extinctmarginal = T)
#' dbdnsat(0:5,mu=c(1,1,1,1,1,1),extinctmarginal = c(F,T,T,F,F,F))
#' dbdnsat(0:5,mu=0,extinctmarginal = T,startval=c(0,0,0,0,-10,-10))
#' @export
dbdnsat=function(n,mu=0,t=1,extantonly=F, log=F, jittervar=0.000001,
                 extinctmarginal=F, startval=0.1,
                 optimpar_plex = list(reltolpar=1e-04, 
                                      reltolval=1e-05, 
                                      abstolpar=1e-07,
                                      maxiter=10000, 
                                      num_cycles=1)){
  d=mapply(dbdnsatnv,n = n, mu=mu, t=t, extantonly=extantonly, log = log, jittervar=jittervar,
           extinctmarginal=extinctmarginal, startval=startval, 
           MoreArgs = list(optimpar_plex=optimpar_plex),
           SIMPLIFY = T)
  return(d)
}


  
#' Saturated log-likelihood for the whole dataset
#'
#' Calculates a saturated log-likelihood as a sum of individual saturated log probabilities for individual data points, i.e. the maximum likelihood if lambda and mu are left free for every data point. Allows for marginalization of distributions using supplied mu values.
#' @param data A data frame with columns n, tipspec and tipext, referring to number of species resulting from every tip after time t, estimated tip speciation and tip extinction rates 
#' @param mu Assumed extinction rates of a birth-death process in case we use extinction marginalisation
#' @param t Time period after which species are observed
#' @param extantonly Condition the distribution for non extinction, default is F
#' @param jittervar Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation
#' @param extinctmarginal Do we marginalize the distribution over the assumed mu values, default is F
#' @param startval Starting values of the subplex optimisation of lambda in case we marginalize over assumed mu
#' @param optimpar_plex The control list of subplex optimisation of lambda in case we marginalize over assumed mu
#' @return Saturated log-likelihood value
#' @examples 
#' testdata=data.frame(n=1:10,tipspec=1:10,tipext=rep(0,10))
#' bdnsatloglik(testdata,t=1,extinctmarginal = F)
#' bdnsatloglik(testdata,t=1,extinctmarginal = T)
#' bdnsatloglik(testdata,mus=1,t=1,extinctmarginal = T)
#' bdnsatloglik(testdata,mus=rep(1,10),t=1,extinctmarginal = T)
#' bdnsatloglik(testdata,t=1,extinctmarginal = T, startval = c(0,0,0,0,0,0,0,0,0,0))
#' bdnsatloglik(testdata,t=1,extinctmarginal = T, startval = c(-10,0,0,0,0,0,0,0,0,0))
#' bdnsatloglik(testdata,t=1,extinctmarginal = T, startval = c(0,0,0,0,0,0,0,0,0,-10))
#' @export
bdnsatloglik=function(data, mus=0, t=1, extantonly=F, jittervar=0.000001,
                      extinctmarginal=F, startval=0.1,
                      optimpar_plex = list(reltolpar=1e-04, 
                                           reltolval=1e-05, 
                                           abstolpar=1e-07,
                                           maxiter=10000, 
                                           num_cycles=1)){
  llhsat=sum(dbdnsat(n = data$n, mu=mus,t=t, extantonly=extantonly, jittervar=jittervar, log = TRUE, 
                     extinctmarginal=extinctmarginal, startval=startval, optimpar_plex=optimpar_plex))
  return(llhsat)
}




#' Saturated log-likelihood for the whole dataset (with linear conditioning of extinction)
#'
#' Calculates a saturated log-likelihood as a sum of individual saturated log probabilities for individual data points, with marginalization of extinction as a maximum likelihood fitted linear function fo tip observed extinction values.
#' @param data A data frame with columns n, tipspec and tipext, referring to number of species resulting from every tip after time t, estimated tip speciation and tip extinction rates 
#' @param t Time period after which species are observed
#' @param startvals Starting values of the subplex optimisation of lambdas, mu intercept and mu slope in this order, vector of length n+2
#' @param extantonly Condition the distribution for non extinction, default is F
#' @param jittervar Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation
#' @param optimpar_plex The control list of subplex optimisation of lambdas and mu parameters
#' @param onlyfout Return just functional value or full output of optimisation, default is T
#' @return Saturated log-likelihood value
#' @examples 
#' testdata=data.frame(n=0:9,tipspec=1:10,tipext=rep(1,10))
#' bdnsatloglik(testdata,t=1,extinctmarginal = F)
#' bdnsatloglik(testdata,mus=testdata$tipext,t=1,extinctmarginal = T)
#' dbdnsatlinear(testdata, t=1, startvals=c(rep(0.1,length(testdata$n)),0.1,0.1), onlyfout = F)
#' @export
dbdnsatlinear=function(data, t, 
                       startvals, 
                       extantonly=F, jittervar=0.0001,
                       optimpar_plex = list(reltolpar=1e-04, 
                                            reltolval=1e-05, 
                                            abstolpar=1e-07,
                                            maxiter=10000, 
                                            num_cycles=1),
                       onlyfout=T) 
              
{
  fn=function(params){
    #break down the
    lambdas=params[1:length(data$n)]
    mualpha=params[length(data$n)+1]
    mubeta=params[length(data$n)+2]
    
    #linear restriction on extinction
    muhat=mualpha+mubeta*data$tipext
    
    #likelihood parsing 
    llh=sum(dbdn(n = data$n, lamb = lambdas, mu = muhat, t=t, extantonly=extantonly, log = T, jittervar = jittervar))
    return(llh)
  }
  
  out=optimglwrap(optimmethod = "subplex",
                optimpar_plex = optimpar_plex, 
                fc = fn,
                trparsopt = startvals)
  
  if(onlyfout){
    out2=out$fvalues
  } else{
    out2=list(lambdas=out$par[1:length(data$n)],mualpha=out$par[length(data$n)+1], mubeta=out$par[length(data$n)+2], llh=out$fvalues)
  }
  return(out2)
  
}



###############################
#maximum likelihood under different free parameters

#' Maximum likelihood estimation of linear regression using dbdn distribution
#'
#' Estimates maximum likelihood values of linear regression parameters of a model where number of species after time t is drawn from a dbdn distribution with lambda and mu estimated as linear functions of tip observed lambda and mus.  
#' @param data A data frame with columns n, tipspec and tipext, referring to number of species resulting from every tip after time t, estimated tip speciation and tip extinction rates 
#' @param t Time period after which species are observed
#' @param startval Starting or fixed values for the optimisation of model parameters, lambda intercept, lambda slope, mu intercept and mu slope in this order, vector of length 4
#' @param freeparid Indices of model parameters startval that are free to be be optimised
#' @param extantonly Condition the distribution for non extinction, default is F
#' @param jittervar Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation
#' @param optimmethod Optimisation algorithm to use, default local "subplex" was sufficient in our testing, but "simplex"  or global "GenSA" or "DEoptim" can also be used
#' @param optimpar_plex List of optimisation parameters passed to local "subplex" or "simplex" routines if these are used
#' @param optimpar_DE List of optimisation parameters passed to "DE_optim" if this algorithm used
#' @param optimpar_GenSA List of optimisation parameters passed to "GenSA" if this algorithm used
#' @param convdetails Return convergence diagnostics, default is F
#' @return A list of all parameter values par, free parameter values parexp, maximum log-likelihood llh and optionally convergence details 
#' @examples 
#' testdata=data.frame(n=1:10,tipspec=1:10,tipext=rep(0,10))
#' 
#' bdnmaxlikval=bdnmaxlik(data=testdata, t=1, startval=c(0,1,0,0),freeparid = 1:4)
#' bdnmaxlikval$par
#' bdnmaxlikval$llh
#' 
#' bdnsatloglikval=bdnsatloglik(testdata,1)
#' bdnsatloglikval
#' 
#' testdata2=data.frame(n=1:10,tipspec=1:10,tipext=20:11)
#' 
#' bdnmaxlikval=bdnmaxlik(data=testdata2, t=1, startval=c(0,1,0,0),freeparid = 1:4)
#' bdnmaxlikval$par
#' bdnmaxlikval$llh
#' 
#' bdnsatloglikval=bdnsatloglik(testdata,1)
#' bdnsatloglikval
#' 
#' restbdnmaxlikval=bdnmaxlik(data=testdata2, t=1, startval=c(0,0,0,0),freeparid = c(1,3,4))
#' restbdnmaxlikval
#' @export
bdnmaxlik=function(data, t, 
                   startval, freeparid,
                   extantonly=F, jittervar=0.0001,
                   optimmethod = "subplex",
                   optimpar_plex = list(reltolpar=1e-04, 
                                        reltolval=1e-05, 
                                        abstolpar=1e-07,
                                        maxiter=10000, 
                                        num_cycles=1), 
                   optimpar_DE = list(lower=-100,
                                      upper=100,
                                      VTR=-Inf, 
                                      strategy = 2, 
                                      bs = FALSE, 
                                      NP = 100,
                                      itermax = 200, 
                                      CR = 0.5, 
                                      F = 0.8, 
                                      trace = T, 
                                      initialpop = NULL,
                                      storepopfrom = 1001, 
                                      storepopfreq = 1, 
                                      p = 0.2, 
                                      c = 0, 
                                      reltol=sqrt(.Machine$double.eps), 
                                      steptol=200,
                                      parallelType="none",
                                      cluster=NULL, 
                                      packages = c(), 
                                      parVar = c(),
                                      foreachArgs = list(), 
                                      parallelArgs = NULL), 
                   optimpar_GenSA = list(lower = -100,
                                         upper = 100,
                                         maxit = 5000,
                                         nb.stop.improvement = 1e+06,
                                         smooth = TRUE,
                                         max.call = 1e+07,
                                         max.time = 3.154e+07,
                                         temperature = 5230,
                                         visiting.param = 2.62,
                                         acceptance.param = -5,
                                         simple.function = FALSE,
                                         trace.mat = TRUE,
                                         seed = -100377),
                   convdetails=F)
{
  if(length(freeparid)>0) {
    #regression with all four parameters relaxed
    fn=function(parsexp){
      parsin=startval
      parsin[freeparid]=parsexp
      return(bdnregloglik(params=parsin,data=data,t=t, extantonly=extantonly,jittervar=jittervar))
    }
    
    #optimize it
    out=optimglwrap(optimmethod = optimmethod,
                    optimpar_plex = optimpar_plex,
                    optimpar_DE=optimpar_DE,
                    optimpar_GenSA=optimpar_GenSA,
                    fc = fn,
                    trparsopt = startval[freeparid])
    
    parsout2=startval
    parsout2[freeparid]=out$par
    
    out2=list(par=parsout2,parsexp=out$par, llh=out$fvalues)
    if(convdetails){
      out2$convdetails=out$details
    }
  } else{
    out2=list(par=startval,
              parsexp=NULL, 
              llh=bdnregloglik(params=startval,data=data,t=t, extantonly=extantonly, jittervar=jittervar),
              convdetails=NULL)
    
  }
  
  return(out2)
}



###############################
#generalized r square via deviance ratio

#' Generalized Rsquare of linear regression using dbdn distribution
#'
#' Generalized Rsquare of linear regression parameters of a model where number of species after time t is drawn from a dbdn distribution with lambda and mu estimated as linear functions of tip observed lambda and mus.  
#' @param data A data frame with columns n, tipspec and tipext, referring to number of species resulting from every tip after time t, estimated tip speciation and tip extinction rates 
#' @param t Time period after which species are observed
#' @param startvalalt Starting or fixed values for the alternative model, lambda intercept, lambda slope, mu intercept and mu slope in this order, vector of length 4
#' @param freeparidalt Indices of alternative model parameters startvalalt that are free to be be optimised
#' @param startvalnul Starting or fixed values for the null model, lambda intercept, lambda slope, mu intercept and mu slope in this order, vector of length 4
#' @param freeparidnul Indices of alternative model parameters startvalalt that are free to be be optimised, default is c(1,3), that is lambda and mu intercepts
#' @param extantonly Condition the distribution for non extinction, default is F
#' @param jittervar Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation
#' @param extmargtype Type of extinction marginalisation, default is "no" i.e. mu is freely optimised, but can be "sats_as_alt" where extinction estimates froma lternative model are used, "sat1_as_alt_sat2_as_null" where extinction estimates from both alternative and null model are used in the numerator and denominator of R square formula, or "linear" where mus are optimised but constrained to be linear function of tip observed mus
#' @param startvalsat Starting value of extinction marginalisation procedure, default is "as_alt" i.e. taking them from alternative model estimates, but can be user supplied too
#' @param optimmethod Optimisation algorithm to use, default local "subplex" was sufficient in our testing, but "simplex"  or global "GenSA" or "DEoptim" can also be used
#' @param optimpar_plex List of optimisation parameters passed to local "subplex" or "simplex" routines if these are used
#' @param optimpar_DE List of optimisation parameters passed to "DE_optim" if this algorithm used
#' @param optimpar_GenSA List of optimisation parameters passed to "GenSA" if this algorithm used
#' @param convdetails Return convergence diagnostics, default is F
#' @return A list of generalized R square genRsquare,deviance of alternative model devalt, deviance of null model devnull, log-likelihood of alternative model loglikalt, log-likelihood of null model logliknul, saturated log-likelihood in the numerator of Rsquare formula logliksat1, saturated log-likelihood in the numerator of Rsquare formula logliksat2, regression parameters of alternative model paralt, regression parameters of null model parnul, and optionally convergence details 
#' @examples 
#' testdata2=data.frame(n=1:10,tipspec=1:10,tipext=20:11)
#' r=bdngenRsquare(data=testdata2, t=1,
#'                 startvalalt=c(0.5,1,0.5,1), freeparidalt=c(1,2,3,4),
#'                 startvalnul=c(0.5,0,0.5,0), freeparidnul=c(1,3))
#' r
#' 
#' re=bdngenRsquare(data=testdata2, t=1, jittervar = 10,
#'                  startvalalt=c(0.5,1,0.5,1), freeparidalt=c(1,2,3,4),
#'                  startvalnul=c(0.5,0,0.5,0), freeparidnul=c(1,3),
#'                  extantonly = T)
#' re
#' 
#' rm=bdngenRsquare(data=testdata2, t=1,
#'                 startvalalt=c(0.5,1,0.5,1), freeparidalt=c(1,2,3,4),
#'                 startvalnul=c(0.5,0,0.5,0), freeparidnul=c(1,3),
#'                 extmargtype = "linear")
#' rm
#' @export
bdngenRsquare=function(data, t, 
                       startvalalt, freeparidalt,
                       startvalnul, freeparidnul=c(1,3),
                       extantonly=F, jittervar=0.0001,
                       extmargtype="no", startvalsat="as_alt", 
                       optimmethod = "subplex",
                       optimpar_plex = list(reltolpar=1e-08, 
                                                reltolval=1e-10, 
                                                abstolpar=1e-14,
                                                maxiter=10000, 
                                                num_cycles=1), 
                       optimpar_DE = list(lower=-100,
                                              upper=100,
                                              VTR=-Inf, 
                                              strategy = 2, 
                                              bs = FALSE, 
                                              NP = 100,
                                              itermax = 200, 
                                              CR = 0.5, 
                                              F = 0.8, 
                                              trace = T, 
                                              initialpop = NULL,
                                              storepopfrom = 1001, 
                                              storepopfreq = 1, 
                                              p = 0.2, 
                                              c = 0, 
                                              reltol=sqrt(.Machine$double.eps), 
                                              steptol=200,
                                              parallelType="none",
                                              cluster=NULL, 
                                              packages = c(), 
                                              parVar = c(),
                                              foreachArgs = list(), 
                                              parallelArgs = NULL), 
                       optimpar_GenSA = list(lower = -100,
                                                 upper = 100,
                                                 maxit = 1000,
                                                 nb.stop.improvement = 1e+06,
                                                 smooth = TRUE,
                                                 max.call = 1e+07,
                                                 max.time = 3.154e+07,
                                                 temperature = 5230,
                                                 visiting.param = 2.62,
                                                 acceptance.param = -5,
                                                 simple.function = FALSE,
                                                 trace.mat = TRUE,
                                                 seed = -100377),
                       convdetails=F){
  
  #calculation of likelihood of alternative model
  loglikalt=bdnmaxlik(data=data, t=t,
                      startval=startvalalt,freeparid = freeparidalt,
                      extantonly=extantonly, jittervar=jittervar,
                      optimmethod = optimmethod,
                      optimpar_plex = optimpar_plex, 
                      optimpar_DE = optimpar_DE, 
                      optimpar_GenSA = optimpar_GenSA, 
                      convdetails=convdetails)
  if(loglikalt$llh==-Inf) print("The loglikelihood of alternative model is -Inf. Either your model structure or your starting values cannot result in the observed data.")
  
  
  
  #calculation of likelihood of null model
  logliknul=bdnmaxlik(data=data, t=t,
                      startval=startvalnul,freeparid = freeparidnul,
                      extantonly=extantonly, jittervar=jittervar,
                      optimmethod = optimmethod,
                      optimpar_plex = optimpar_plex, 
                      optimpar_DE = optimpar_DE, 
                      optimpar_GenSA = optimpar_GenSA,
                      convdetails=convdetails)
  
  if(logliknul$llh==-Inf) print("The loglikelihood of null model is -Inf. Either your model structure or your starting values cannot result in the observed data.")
  
  
  
  
  # branching across various options of extinction marginalization when calculating saturated likelihood
  if(extmargtype=="sats_as_alt"){ #using alternative model extinction fit for saturated model 
      muhat=loglikalt$par[3]+loglikalt$par[4]*data$tipext
      if(startvalsat[1]=="as_alt"){ #use starting values from alternative model
        startvalsat=loglikalt$par[1]+loglikalt$par[2]*data$tipspec
      }
      logliksat1=bdnsatloglik(data=data, mu=muhat,t=t, extantonly=extantonly,
                              extinctmarginal=T, startval=startvalsat,
                              optimpar_plex = optimpar_plex)
      logliksat2=logliksat1

      
    } else if(extmargtype=="sat1_as_alt_sat2_as_null"){ #using alternative and null model extinctions fit for saturated models
      muhatalt=loglikalt$par[3]+loglikalt$par[4]*data$tipext
      muhatnul=logliknul$par[3]+logliknul$par[4]*data$tipext
      if(startvalsat[1]=="as_alt"){ #use starting values from alternative model
        startvalsat1=loglikalt$par[1]+loglikalt$par[2]*data$tipspec
        startvalsat2=logliknul$par[1]+logliknul$par[2]*data$tipspec
      } else {
        startvalsat1=startvalsat
        startvalsat2=startvalsat
      }
      logliksat1=bdnsatloglik(data=data, mu=muhatalt,t=t, extantonly=extantonly,
                              extinctmarginal=T, startval=startvalsat1,
                              optimpar_plex = optimpar_plex)
      logliksat2=bdnsatloglik(data=data, mu=muhatnul,t=t, extantonly=extantonly,
                              extinctmarginal=T, startval=startvalsat2,
                              optimpar_plex = optimpar_plex)
      
      
    } else if(extmargtype=="linear"){ #make a separate extinction fit for saturated model using linear regression
      if(startvalsat[1]=="as_alt"){  #use starting values from alternative model
        lambdastart=loglikalt$par[1]+loglikalt$par[2]*data$tipspec
        mustart=c(loglikalt$par[3],loglikalt$par[4])
        startvalsat=c(lambdastart, mustart)      
      }

      logliksat1=dbdnsatlinear(data=data, t=t, 
                               extantonly=extantonly, jittervar=jittervar, 
                               startvals=startvalsat, 
                               optimpar_plex=optimpar_plex)
      logliksat2=logliksat1
      
      
    } else if(extmargtype=="no"){ #use no extinction marginalisation
      logliksat1=bdnsatloglik(data=data, t=t, extantonly=extantonly,
                             extinctmarginal=F, startval=startvalsat[1],
                             optimpar_plex = optimpar_plex)
      logliksat2=logliksat1
      
      
    } else{
      print("This marginalizazion type is to be implemented, crashing")
    }
  
  
  if(logliksat1==-Inf) print("The loglikelihood of saturated model 1 is -Inf. Either your model structure or your starting values cannot result in the observed data.")
  
  if(logliksat2==-Inf) print("The loglikelihood of saturated model 2 is -Inf. Either your model structure or your starting values cannot result in the observed data.")
  
  
  
  #deviance calculation
  devalt=(2*logliksat1)-(2*loglikalt$llh)
  devnul=(2*logliksat2)-(2*logliknul$llh)
  
  
  
  #generalized R square calculation
  genRsquare=1-(devalt/devnul)
  
  
  
  #parsing the output list
  res=list(genRsquare=genRsquare,
           devalt=devalt, devnul=devnul, 
           loglikalt=loglikalt$llh, logliknul=logliknul$llh, logliksat1=logliksat1, logliksat2=logliksat2,
           paralt=loglikalt$par,parnul=logliknul$par)
  
  #optionally printing the convergence details
  if(convdetails){
    res$convdetailsalt=loglikalt$details
    res$convdetailsnul=logliknul$details
  }
  
  return(res)
  
}


###############################
#optimizer with local simplex, subplex or global DEoptim and GenSA algorithms, the global ones currently overwhelm the memory stack 

#' Wrapper for various local and global optimiser routines
#'
#' Wrapper function for subplex, simplex, differential evolution (DEoptim) and generalized simulated annealing (GenSA) optimisation routines.
#' @param optimmethod Optimisation algorithm to use, default is local "subplex", but "simplex"  or global "GenSA" or "DEoptim" can also be used
#' @param optimpar_plex List of optimisation parameters passed to local "subplex" or "simplex" routines if these are used, relative tolerance of parameters, relative tolerance of functional values, absolute tolerance of parameters, number of iterations and number of cycles for which the routine is run from the start can be set
#' @param optimpar_DE List of optimisation parameters passed to "DE_optim" if this algorithm used, see DEoptim package for details
#' @param optimpar_GenSA List of optimisation parameters passed to "GenSA" if this algorithm used, see GenSA package for details
#' @param fc Function to be optimised
#' @param trparsopt Starting values used in subplex, simplex or GenSA
#' @return A list of generalized R square genRsquare,deviance of alternative model devalt, deviance of null model devnull, log-likelihood of alternative model loglikalt, log-likelihood of null model logliknul, saturated log-likelihood in the numerator of Rsquare formula logliksat1, saturated log-likelihood in the numerator of Rsquare formula logliksat2, regression parameters of alternative model paralt, regression parameters of null model parnul, and optionally convergence details 
#' @export
optimglwrap=function(optimmethod = 'simplex',
                     optimpar_plex = list(reltolpar=1e-04, 
                                          reltolval=1e-05, 
                                          abstolpar=1e-07,
                                          maxiter=10000, 
                                          num_cycles=1), 
                     optimpar_DE = list(lower=0,
                                        upper=1,
                                        VTR=-Inf, 
                                        strategy = 2, 
                                        bs = FALSE, 
                                        NP = 100,
                                        itermax = 200, 
                                        CR = 0.5, 
                                        F = 0.8, 
                                        trace = T, 
                                        initialpop = NULL,
                                        storepopfrom = 1001, 
                                        storepopfreq = 1, 
                                        p = 0.2, 
                                        c = 0, 
                                        reltol=sqrt(.Machine$double.eps), 
                                        steptol=200,
                                        parallelType="none",
                                        cluster=NULL, 
                                        packages = c(), 
                                        parVar = c(),
                                        foreachArgs = list(), 
                                        parallelArgs = NULL), 
                     optimpar_GenSA = list(lower = 0,
                                           upper = 1,
                                           maxit = 1000,
                                           nb.stop.improvement = 1e+06,
                                           smooth = TRUE,
                                           max.call = 1e+07,
                                           max.time = 3.154e+07,
                                           temperature = 5230,
                                           visiting.param = 2.62,
                                           acceptance.param = -5,
                                           simple.function = FALSE,
                                           trace.mat = TRUE,
                                           seed = -100377),
                     fc,
                     trparsopt,
                     ...)
{
  
  
  #### Simplex
  if(optimmethod == 'simplex'){
    max_cycles=optimpar_plex[5]
    cy=1
    fvalue <- rep(-Inf,max_cycles)
    out <- NULL
    while(cy <= max_cycles){
      # Print cycle number
      if(max_cycles > 1) cat(paste('Cycle ',cy,'\n',sep =''))
      # Simplex
      outnew <- suppressWarnings(DDD::simplex(fun = fc,
                                         trparsopt = trparsopt,
                                         optimpars = c(as.numeric(optimpar_plex[1]),
                                                       as.numeric(optimpar_plex[2]),
                                                       as.numeric(optimpar_plex[3]),
                                                       as.numeric(optimpar_plex[4])),
                                         ...))
      outnew$details=outnew
      
      # Break if any NAs
      if(cy > 1 & (any(is.na(outnew$par)) | any(is.nan(outnew$par)) | is.na(outnew$fvalues) | is.nan(outnew$fvalues) | outnew$conv != 0))
      {cat('The last cycle failed; second last cycle result is returned.\n')
        return(out)} else{
          out <- outnew
          trparsopt <- out$par
          fvalue[cy] <- out$fvalues}
      
      # Break if we get under absolute tolerance
      if(cy > 1){
        if(abs(fvalue[cy] - fvalue[cy - 1]) < optimpar_plex[3]){
          if(cy < max_cycles) cat('No more cycles needed.\n')
          cy <- max_cycles} else if(cy == max_cycles){
            cat('More cycles in optimization recommended.\n')}}
      
      cy <- cy + 1
    }
  }
  
  ### Subplex
  if(optimmethod == 'subplex'){
    max_cycles=optimpar_plex[5]
    cy=1
    fvalue <- rep(-Inf,max_cycles)
    out <- NULL
    while(cy <= max_cycles){
      # Print cycle number
      if(max_cycles > 1) cat(paste('Cycle ',cy,'\n',sep =''))
      # Invert function
      minfc <- function(par){return(-fc(par))}
      # Subplex
      outnew <- suppressWarnings(subplex::subplex(par = trparsopt,
                                                  fn = minfc,
                                                  control = list(abstol = optimpar_plex[3],
                                                                 reltol = optimpar_plex[1],
                                                                 maxit = optimpar_plex[4]),
                                                  ...))
      outnew <- list(par = outnew$par, fvalues = -outnew$value, conv = outnew$convergence, details=outnew)
      
      # Break if any NAs
      if(cy > 1 & (any(is.na(outnew$par)) | any(is.nan(outnew$par)) | is.na(outnew$fvalues) | is.nan(outnew$fvalues) | outnew$conv != 0))
      {cat('The last cycle failed; second last cycle result is returned.\n')
        return(out)} else{
          out <- outnew
          trparsopt <- out$par
          fvalue[cy] <- out$fvalues}
      
      # Break if we get under absolute tolerance
      if(cy > 1){
        if(abs(fvalue[cy] - fvalue[cy - 1]) < optimpar_plex[3]){
          if(cy < max_cycles) cat('No more cycles needed.\n')
          cy <- max_cycles} else if(cy == max_cycles){
            cat('More cycles in optimization recommended.\n')}}
      
      cy <- cy + 1
    }
  }
  
  
  ### DEoptim
  if(optimmethod == 'DEoptim'){
    # Invert function
    minfc <- function(par){return(-fc(par))}
    outnew <- suppressWarnings(DEoptim::DEoptim(fn = minfc,
                                                lower = rep(as.numeric(optimpar_DE[1]), length(trparsopt)),
                                                upper = rep(as.numeric(optimpar_DE[2]), length(trparsopt)),
                                                control = optimpar_DE[3:23],
                                                ...))
    out <- list(par = outnew$optim$bestmem, fvalues = -outnew$optim$bestval, conv = 0, details=outnew)
    
  }
  
  ### GenSA
  if(optimmethod == "GenSA"){
    # Invert function
    minfc <- function(par){return(-fc(par))}
    # GenSA
    outnew <- suppressWarnings(GenSA::GenSA(par = trparsopt,
                                            fn = minfc,
                                            lower = rep(as.numeric(optimpar_GenSA[1]), length(trparsopt)),
                                            upper = rep(as.numeric(optimpar_GenSA[2]), length(trparsopt)),
                                            control = optimpar_GenSA[3:13],
                                            ...))
    out <- list(par = outnew$par, fvalues = -outnew$value, conv = 0, details=outnew)
    
    
  }
  
  return(out)
}



###############################
#Simulation function using bd simulated subtrees
#we need a simulation procedure that allows to simulate trees from root, not crown
# plot(TreeSim::sim.bd.age(age=1,numbsim=1,lambda=20,mu=10)[[1]])
# axisPhylo()


#' Simulation procedure for input data of bdn linear models
#'
#' A function generating numbers data frame with numbers of species n resulting from a birth-death processes with given lambdas, mus and ts. All these three parameters are vectorized over.
#' @param lambdas Vector od speciation rates lambda
#' @param mus Vector of extinction rates mu
#' @param ts Vector of waiting times t
#' @param extantonly Iterate until a surviving tree is produced
#' @return A data frame with columns n, tipspec, tipext and time, referring to number of species resulting from every tip after given  tip speciation, tip extinction rates, and waiting times 
#' @examples
#' bdntreesim(rep(1,1000),rep(0,1000),rep(1,1000))
#' mean(bdntreesim(rep(1,1000),rep(0,1000),rep(1,1000))$n)
#' exp(1)
#' 
#' mean(bdntreesim(rep(2,1000),rep(0,1000),rep(1,1000))$n)
#' exp(2)
#' 
#' mean(bdntreesim(rep(2,1000),rep(1,1000),rep(1,1000))$n)
#' exp(1)
#' @export
bdntreesim=function(lambdas,mus,ts,extantonly=F){
  fn=function(lambda, mu, t, extantonly){
    repeat{
      tr=TreeSim::sim.bd.age(age=t,numbsim=1,lambda=lambda,mu=mu)[[1]]
      
      if(geiger::is.phylo(tr)){
        num=length(tr$tip.label)-length(geiger::is.extinct(tr))
      } else {
        num=tr
      }
      if(!extantonly | num>=1){
        break
      }
    }
    return(num)

  } 
  
  sizevec=mapply(fn,lambda=lambdas , mu=mus, t=ts, extantonly=extantonly, SIMPLIFY = T)
  d=data.frame(n=sizevec,
               tipspec=rep_len(lambdas,length(sizevec)),
               tipext=rep_len(mus,length(sizevec)),
               time=rep_len(ts,length(sizevec))
  )
  
  if(extantonly){
    d=d[which(d$n>0),]
  }
  return(d)
}



