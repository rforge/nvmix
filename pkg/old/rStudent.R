### rStudent() ###################################################################

##' @title Random Number Generator for the Multivariate Student t Distribution
##' @param n sample size
##' @param loc location vector of dimension d
##' @param scale covariance matrix of dimension (d, d)
##' @param df degrees of freedom (positive real or Inf in which case the density
##'        of a N(loc, scale) is evaluated)
##' @param factor factorization matrix of the covariance matrix scale; a matrix
##'        R such that R^T R = scale. R is multiplied to the (n, d)-matrix of
##'        independent N(0,1) random variates in the construction from the *right*
##'        (hence the notation).
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Marius Hofert, Erik Hintz

rStudent <- function(n, loc = rep(0, d), scale, df, factor = factorize(scale))
{
  d <- nrow(as.matrix(factor))
  
  rnvmix(n = n, mix = "inverse.gamma", loc = loc, scale = scale, factor = factor, df = df)
}
