### rnvmix() ###################################################################

##' @title Random Number Generator for Multivariate Normal Variance Mixtures
##' @param n sample size
##' @param df degrees of freedom (positive real or Inf in which case samples
##'        from N(loc, scale) are drawn).
##' @param loc location vector of dimension d
##' @param scale covariance matrix of dimension (d, d)
##' @param factor factorization matrix of the covariance matrix scale; a matrix
##'        R such that R^T R = scale. R is multiplied to the (n, d)-matrix of
##'        independent N(0,1) random variates in the construction from the *right*
##'        (hence the notation).
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Marius Hofert
rnvmix <- function(n, df, loc = rep(0, d), scale, factor = factorize(scale))
{
    ## Checks
    d <- nrow(as.matrix(factor))
    stopifnot(n >= 1, df > 0)
    ## Generate Z ~ N(0, I)
    Z <- matrix(rnorm(n * d), ncol = d) # (n, d)-matrix of N(0, 1)
    ## Generate Y ~ N(0, scale)
    Y <- Z %*% factor # (n, d) %*% (d, k) = (n, k)-matrix of N(0, scale); allows for different k
    ## Generate Y ~ t_nu(0, scale)
    ## Note: sqrt(W) for W ~ df/rchisq(n, df = df) but rchisq() calls rgamma(); see ./src/nmath/rchisq.c
    ##       => W ~ 1/rgamma(n, shape = df/2, rate = df/2)
    if(is.finite(df)) {
        df2 <- df/2
        Y <- Y / rgamma(n, shape = df2, rate = df2) # also fine for different k
    }
    ## Generate X ~ t_nu(loc, scale)
    sweep(Y, 2, loc, "+") # X
}
