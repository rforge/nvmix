### rnvmix() ###################################################################

##' @title Random Number Generator for Multivariate Normal Variance Mixtures
##' @param n sample size
##' @param mix specification of the (mixture) distribution of W. This can be:
##'        1) character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "r", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) function being interpreted as the quantile function F_W^-.
##'        4) n-vector containing a random sample from W.
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param factor factor R of the covariance matrix 'scale' such that R^T R
##'        = 'scale'
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Marius Hofert
##' @note - For the Student t distribution, W ~ df/rchisq(n, df = df) but
##'         rchisq() simply calls rgamma(); see ./src/nmath/rchisq.c
##'         => W ~ 1/rgamma(n, shape = df/2, rate = df/2)
##'       - For a generalized inverse Gaussian distribution one could use:
##'         + "Runuran": faster if n large and parameters fixed; based on density
##'         + "GIGrvg":  faster if n small and often called with several parameters
##'         see examples of 'GIGrvg' for both methods
rnvmix <- function(n, mix, loc = rep(0, d), scale = diag(2),
                   factor = factorize(scale), ...)
{
    ## Checks
    d <- nrow(as.matrix(factor))
    stopifnot(n >= 1, df > 0)
    ## Generate W (do that first, so that results are reproducible when user
    ## provides realizations of W; see examples on ?rnvmix)
    W <- if(is.character(mix)) { # 'mix' is a character vector specifying supported mixture distributions (utilizing '...')
             mix <- match.arg(mix, choices = c("constant", "inverse.gamma"))
             switch(mix,
                    "constant" = {
                        rep(1, n)
                    },
                    "inverse.gamma" = {
                        if(hasArg(df)) df <- list(...)$df else
                            stop("'mix = \"inverse.gamma\"' requires 'df' to be provided.")
                        ## Still allow df = Inf (normal distribution)
                        stopifnot(is.numeric(df), length(df) == 1, df > 0)
                        if(is.finite(df)) {
                            df2 <- df/2
                            1 / rgamma(n, shape = df2, rate = df2)
                        } else {
                            rep(1, n)
                        }
                    },
                    stop("Currently unsupported 'mix'"))
         } else if(is.list(mix)) { # 'mix' is a list of the form (<character string>, <parameters>)
             stopifnot(length(mix) >= 1, is.character(distr <- mix[[1]]))
             rmix <- paste0("r", distr)
             if(!existsFunction(rmix))
                 stop("No function named '", rmix, "'.")
             do.call(rmix, c(n, mix[-1]))
         } else if(is.function(mix)) { # 'mix' is interpreted as the quantile function F_W^- of the mixture distribution F_W of W
             mix(runif(n))
         } else if(is.numeric(mix) && length(mix) == n && all(mix >= 0)) { # 'mix' is the vector of realizations of W
             mix
         } else stop("'mix' must be a character string, list, quantile function or n-vector of non-negative random variates.")
    ## Generate Z ~ N(0, I)
    Z <- matrix(rnorm(n * d), ncol = d) # (n, d)-matrix of N(0, 1)
    ## Generate Y ~ N(0, scale)
    Y <- Z %*% factor # (n, d) %*% (d, k) = (n, k)-matrix of N(0, scale); allows for different k
    ## Generate X ~ M_k(0, Sigma, LS[F_W])
    X <- sqrt(W) * Y # also fine for different k
    ## Generate X ~ M_k(mu, Sigma, LS[F_W])
    sweep(X, 2, loc, "+")
}
