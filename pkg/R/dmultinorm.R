### dmultinorm() ###################################################################

##' @title Density of the Multivariate t distribution
##' @param x (n, d)-matrix of evaluation points
##' @param loc location vector of dimension d
##' @param scale covariance matrix of dimension (d, d)
##' @param factor factorization matrix of the covariance matrix scale;
##'        caution: this has to be an *upper triangular* matrix R
##'        such that R^T R = scale here (otherwise det(scale) not computed correctly)
##' @param abserr numeric and non-negative. Absolute precision required. If abserr = 0, algorithm will run until total number of function evaluations exceeds Nmax (see also Nmax)
##' @param gam Monte Carlo confidence multiplier. Algorithm runs until gam * (estimated standard error) < abserr. gam = 3.3 means that one can expect
##'        that in 99.9% of the cases the actual absolute error is less than abserr.
##' @param Nmax Total number of function evaluations allowed.
##' @param B Number of randomizations to get error estimate.
##' @param n_init First loop uses n_init function evaluations. Any positive integer allowed, powers or at least multiples of 2 are recommended.
##' @param method Character string indicating method to be used. Allowed are
##'         - "sobol" for a Sobol sequence.
##'         - "ghalton" for a generalized Halton sequence.
##'         - "prng" for a pure Monte Carlo approach.
##' @return n-vector with N(loc, scale) density values
##' @author Marius Hofert and Erik Hintz

dmultinorm <- function(x, loc = rep(0, d), scale,
                     factor = tryCatch(factorize(scale), error = function(e) e), # needs to be triangular!
                     log = FALSE){
  
  if(!is.matrix(x)) x <- rbind(x)
  d <- ncol(x)
  
  dnvmix(x = x, loc = loc, scale = scale, mix = "constant", factor = factor, log = log)
  
}


