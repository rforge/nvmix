### pmultinorm() ###################################################################

##' @title Distribution Function of the Multivariate Normal Distribution
##' @param upper vector of upper limits
##' @param lower vector of lower limits
##' @param scale covariance matrix of dimension (q, q)
##' @param standardized logical. If TRUE, scale is assumed to be a correlation matrix; if FALSE (default), lower, upper and scale will be normalized.
##' @param abserr numeric and non-negative. Absolute precision required. If abserr = 0, algorithm will run until total number of function evaluations exceeds Nmax (see also Nmax)
##' @param gam Monte Carlo confidence multiplier. Algorithm runs until gam * (estimated standard error) < abserr. gam = 3.3 means that one can expect 
##'        that in 99.9% of the cases the actual absolute error is less than abserr.
##' @param Nmax Total number of function evaluations allowed. 
##' @param N Number of randomizations to get error estimate. 
##' @param n_init First loop uses n_init function evaluations. Any positive integer allowed, powers or at least multiples of 2 are recommended. 
##' @param precond Logical. If TRUE (recommended), variable reordering as described in [genzbretz2002] pp. 955-956 is performed. Variable reordering can lead to a significant variance
##'        reduction and decrease in computational time. 
##' @param method Character string indicating method to be used. Allowed are
##'         - "sobol" for a Sobol sequence.
##'         - "ghalton" for a generalized Halton sequence.
##'         - "prng" for a pure Monte Carlo approach. 
##' @author Erik Hintz
pmultinorm <- function(upper, lower = rep(-Inf, length(upper)), mean = rep(0, length(upper)), scale, standardized = FALSE, gam = 3.3, abserr = 0.001, Nmax = 1e8, N = 12, n_init = 2^6, precond = TRUE, method = "sobol")
{
  pnvmix(upper = upper, lower = lower, shift = mean, scale = scale, mix = "constant", standardized = standardized, gam = gam, abserr = abserr, Nmax = Nmax, N = N, n_init = n_init, precond = precond, method = method)
}