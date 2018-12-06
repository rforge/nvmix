### d/p/rnvmixcop() ############################################################


##' Density function of a Multivariate Normal Variance Mixture Copula
##' @param u (n,d) matrix of evaluation points. Have to be in (0,1)
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale matrix).
##' @param control see ?pnvmixcop()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @param log logical indicating whether the logarithmic density is to be computed
##' @param ... see ?pnvmix
##' @author Erik Hintz and Marius Hofert
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)

dnvmixcop <- function(u, qmix, scale = diag(d), control = list(), 
                      verbose = FALSE, log = FALSE, ...){
  
  ## Most arguments are checked by qnvmix() and pnvmix()
  if(!is.matrix(u)) u <- rbind(u)
  d <- ncol(u) # for 'scale'
  
  ## Obtain quantiles. Note that qnvmix() takes in and returns a vector
  qu <- qnvmix(as.vector(u), qmix = qmix, control = control, 
               verbose = verbose, q.only =  FALSE, ...) # length n*d 

  ## log f_{X, scale} (F_{X1}^{-1}(u_{j1}),...,F_X1 ^{-1}(u_{jd})), j = 1,...,n
  num <- dnvmix(matrix(qu$q, ncol = d), qmix = qmix, scale = scale, 
                control = control, verbose = verbose, log = TRUE, ...)# length n
  
  ## sum_{i=1}^d log f_{X1}( F_{X1}^{-1}(u_{ji})), j = 1,..,n
  ## Note that the log-density values are already calculated by qnvmix() 
  denom <- rowSums(matrix(qu$log.density, ncol = d)) # length n
  
  if(log) num - denom else exp(num - denom)
}



##' Distribution function of a Multivariate Normal Variance Mixture Copula
##' @param u (n,d) matrix of evaluation points. Have to be in (0,1)
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale matrix).
##' @param control see ?pnvmixcop()
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @author Erik Hintz and Marius Hofert
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)

pnvmixcop <- function(u, qmix, scale = diag(d), control = list(), 
                      verbose = FALSE, ...){
  
  ## Most arguments are checked by qnvmix() and pnvmix()
  if(!is.matrix(u)) u <- rbind(u)
  d <- ncol(u) # for 'scale'
  ## Obtain quantiles. Note that qnvmix() returns a vector
  qu <- matrix(qnvmix(as.vector(u), qmix = qmix, control = control, 
                      verbose = verbose, q.only = TRUE, ...), ncol = d)
  pnvmix(qu, qmix = qmix, scale = scale, control = control, 
         verbose = verbose, ...)
}


##' Distribution function of a Multivariate Normal Variance Mixture Copula
##' @param u (n,d) matrix of evaluation points. Have to be in (0,1)
##' @param qmix see ?pnvmix
##' @param scale (d, d)-covariance matrix (scale matrix).
##' @param control see ?pnvmixcop()
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @author Erik Hintz and Marius Hofert
##' @return (n, d)-matrix with NVM(0, scale)-copula samples 

rnvmixcop <- function(n, qmix, scale = diag(2), factor = NULL,
                      method = c("PRNG", "sobol", "ghalton"), skip = 0, 
                      control = list(), verbose = TRUE, ...)
{               
  d <- dim(scale)[1]
  stopifnot(dim(scale) == c(d,d))
  scale <- cov2cor(scale) # only need correlation matrix
  
  ## Sample from the nvmix dist'n
  sample.nvmix <- rnvmix(n = n, qmix = qmix, scale = scale, factor = factor,
                method = method, skip = skip, ...)
  ## Apply univariate margins
  ## Need (n,1) matrix as input so that pnvmix() gets the dimension right:
  sample.nvmixcop <- pnvmix(upper = matrix(sample.nvmix, ncol = 1), qmix = qmix, 
                            scale = matrix(1), standardized = TRUE, 
                            control = control, verbose = verbose, ...)
  ## Get dimensions correct and return
  matrix(sample.nvmixcop, ncol = d)
}




# 
# library(mvtnorm)
# d <- 2 # dimension
# rho <- 0.7 # off-diagonal entry of the correlation matrix P
# P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
# diag(P) <- 1
# n <- 10000
# 
# set.seed(64)
# k <- 100
# u <- matrix(runif(k*d), ncol = d) # generate two random evaluation points
# 
# nu = 2
# df = nu
# qmix. <- function(u) 1/qgamma(u, shape = df/2, rate = df/2)
# tc <- tCopula(rho, dim = d, df = nu)
# 
# 
# 
# 
# pnvmixcop(u, qmix = "inverse.gamma", scale = P, df = nu)
# dnvmixcop(u, qmix = "inverse.gamma", scale = P, df = nu)
# 
# pnvmixcop(u, qmix = qmix., scale = P)
# dnvmixcop(u, qmix = qmix., scale = P)
# 
# 
# pCopula(u, copula = tc) # value of the copula at u
# dCopula(u, copula = tc) # value of the copula at u
# 
# 
# sample <- rnvmixcop(n, qmix = "inverse.gamma", df = nu, scale = P)
# plot(sample)
# # 












