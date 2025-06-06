% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdngenRsquare.R
\name{bdngenRsquare}
\alias{bdngenRsquare}
\title{Generalized Rsquare of linear regression using dbdn distribution}
\usage{
bdngenRsquare(
  data,
  t,
  startvalalt,
  freeparidalt,
  startvalnul,
  freeparidnul = c(1, 3),
  extantonly = F,
  jittervar = 1e-04,
  extmargtype = "no",
  startvalsat = "as_alt",
  optimmethod = "subplex",
  optimpar_plex = list(reltolpar = 1e-08, reltolval = 1e-10, abstolpar = 1e-14, maxiter =
    10000, num_cycles = 1),
  optimpar_DE = list(lower = -100, upper = 100, VTR = -Inf, strategy = 2, bs = FALSE, NP
    = 100, itermax = 200, CR = 0.5, F = 0.8, trace = T, initialpop = NULL, storepopfrom =
    1001, storepopfreq = 1, p = 0.2, c = 0, reltol = sqrt(.Machine$double.eps), steptol =
    200, parallelType = "none", cluster = NULL, packages = c(), parVar = c(), foreachArgs
    = list(), parallelArgs = NULL),
  optimpar_GenSA = list(lower = -100, upper = 100, maxit = 1000, nb.stop.improvement =
    1e+06, smooth = TRUE, max.call = 1e+07, max.time = 31540000, temperature = 5230,
    visiting.param = 2.62, acceptance.param = -5, simple.function = FALSE, trace.mat =
    TRUE, seed = -100377),
  convdetails = F
)
}
\arguments{
\item{data}{A data frame with columns n, tipspec and tipext, referring to number of species resulting from every tip after time t, estimated tip speciation and tip extinction rates}

\item{t}{Time period after which species are observed}

\item{startvalalt}{Starting or fixed values for the alternative model, lambda intercept, lambda slope, mu intercept and mu slope in this order, vector of length 4}

\item{freeparidalt}{Indices of alternative model parameters startvalalt that are free to be be optimised}

\item{startvalnul}{Starting or fixed values for the null model, lambda intercept, lambda slope, mu intercept and mu slope in this order, vector of length 4}

\item{freeparidnul}{Indices of alternative model parameters startvalalt that are free to be be optimised, default is c(1,3), that is lambda and mu intercepts}

\item{extantonly}{Condition the distribution for non extinction, default is F}

\item{jittervar}{Standard deviation of a normal distribution by which lamb and mu are jittered in boundary situation}

\item{extmargtype}{Type of extinction marginalisation, default is "no" i.e. mu is freely optimised, but can be "sats_as_alt" where extinction estimates froma lternative model are used, "sat1_as_alt_sat2_as_null" where extinction estimates from both alternative and null model are used in the numerator and denominator of R square formula, or "linear" where mus are optimised but constrained to be linear function of tip observed mus}

\item{startvalsat}{Starting value of extinction marginalisation procedure, default is "as_alt" i.e. taking them from alternative model estimates, but can be user supplied too}

\item{optimmethod}{Optimisation algorithm to use, default local "subplex" was sufficient in our testing, but "simplex"  or global "GenSA" or "DEoptim" can also be used}

\item{optimpar_plex}{List of optimisation parameters passed to local "subplex" or "simplex" routines if these are used}

\item{optimpar_DE}{List of optimisation parameters passed to "DE_optim" if this algorithm used}

\item{optimpar_GenSA}{List of optimisation parameters passed to "GenSA" if this algorithm used}

\item{convdetails}{Return convergence diagnostics, default is F}
}
\value{
A list of generalized R square genRsquare,deviance of alternative model devalt, deviance of null model devnull, log-likelihood of alternative model loglikalt, log-likelihood of null model logliknul, saturated log-likelihood in the numerator of Rsquare formula logliksat1, saturated log-likelihood in the numerator of Rsquare formula logliksat2, regression parameters of alternative model paralt, regression parameters of null model parnul, and optionally convergence details
}
\description{
Generalized Rsquare of linear regression parameters of a model where number of species after time t is drawn from a dbdn distribution with lambda and mu estimated as linear functions of tip observed lambda and mus.
}
\examples{
testdata2=data.frame(n=1:10,tipspec=1:10,tipext=20:11)
r=bdngenRsquare(data=testdata2, t=1,
                startvalalt=c(0.5,1,0.5,1), freeparidalt=c(1,2,3,4),
                startvalnul=c(0.5,0,0.5,0), freeparidnul=c(1,3))
r

re=bdngenRsquare(data=testdata2, t=1, jittervar = 10,
                 startvalalt=c(0.5,1,0.5,1), freeparidalt=c(1,2,3,4),
                 startvalnul=c(0.5,0,0.5,0), freeparidnul=c(1,3),
                 extantonly = T)
re

rm=bdngenRsquare(data=testdata2, t=1,
                startvalalt=c(0.5,1,0.5,1), freeparidalt=c(1,2,3,4),
                startvalnul=c(0.5,0,0.5,0), freeparidnul=c(1,3),
                extmargtype = "linear")
rm
}
