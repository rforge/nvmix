### rnvmix() ###################################################################

##' @title Random Number Generator for Multivariate Normal Variance Mixtures
##'        and Gamma mixtures
##' @param n see ?rnvmix()
##' @param rmix ?rnvmix()
##' @param qmix ?rnvmix()
##' @param loc ?rnvmix()
##' @param scale ?rnvmix()
##' @param factor ?rnvmix()
##' @param method ?rnvmix()
##' @param skip ?rnvmix()
##' @param which string, either "nvmix" for generating NVM(loc, scale, F_W) samples
##'        or "maha2" for generating samples from the squared maha distances of
##'        NVM(loc, scale, F_W)
##' @param ... ?rnvmix()
##' @return ?rnvmix()
##' @author Marius Hofert and Erik Hintz
##' @note - For the Student t distribution, W ~ df/rchisq(n, df = df) but
##'         rchisq() simply calls rgamma(); see ./src/nmath/rchisq.c
##'         => W ~ 1/rgamma(n, shape = df/2, rate = df/2)
##'       - For a generalized inverse Gaussian distribution one could use:
##'         + "Runuran": faster if n large and parameters fixed; based on density
##'         + "GIGrvg":  faster if n small and often called with several parameters
##'         see examples of 'GIGrvg' for both methods
##'       - user friendly wrappers are provided in 'rnvmix()' and 'rgammamix()'
rnvmix_ <- function(n, rmix, qmix, loc = rep(0, d), scale = diag(2), factor = NULL,
                    method = c("PRNG", "sobol", "ghalton"), skip = 0, 
                    which = c("nvmix", "maha2"), ...)
{
    ## Basic checks
    stopifnot(n >= 1)
    method <- match.arg(method)
    ## Set [q/r]mix to NULL if not provided (needed for 'get_mix_()' below)
    if(!hasArg(qmix)) qmix <- NULL
    if(!hasArg(rmix)) rmix <- NULL
    ## Deal with 'factor' (more general here than in dnvmix() and pnvmix1())
    if(is.null(factor)) { # => let 'factor' (internally here) be an *upper* triangular matrix
        factor <- chol(scale) # *upper* triangular; by this we avoid t() internally here and below around Z
        d <- nrow(factor)
        k <- d # => factor a square matrix here
    } else { # => 'factor' is a provided (d, k)-matrix (factor %*% factor^T = (d, d)-scale)
        d <- nrow(factor <- as.matrix(factor)) # (d, k)-matrix here...
        k <- ncol(factor)
        factor <- t(factor) # ... converted to a (k, d)-matrix to avoid t() below around Z
    }
    ## Determine if inversion is to be used
    inversion <- FALSE
    ## This is the case if the method used is "sobol" or "ghalton"
    if(method != "PRNG") inversion <- TRUE
    ## Or if the method is  "PRNG" but 'rmix' was not provivded
    if(method == "PRNG" && is.null(rmix)) inversion <- TRUE
    ## Check if 'qmix' provided for inversion based methods
    if(inversion & is.null(qmix)) 
       stop("'qmix' needs to be provided for methods 'sobol' and 'ghalton'")
    ## Logical if supplied 'rmix' is a vector of realizations
    rmix.sample <- (is.numeric(rmix) & length(rmix) == n & all(rmix > 0))
    if(!rmix.sample){
       mix_list      <- get_mix_(qmix = qmix, rmix = rmix, 
                                 callingfun = "rnvmix", ... ) 
       mix_          <- mix_list[[1]] # function(u) or function(n)
       special.mix   <- mix_list[[2]] # string or NA
       use.q         <- mix_list$use.q # logical if 'mix_' is a quantile function
    }
    ## Obtain n realizations of the mixing rv 
    W <- if(inversion){
       ## Get low discrepancy pointset
       ## For 'which = "nvmix"' need k-dim normal distribution later => k+1 uniforms;
       ## otherwise 2 (one for 'W', one for the gamma)
       dim. <- if(which == "nvmix") k + 1 else 2
       U <- switch(method,
                   "sobol" = {
                      qrng::sobol(n, d = dim., randomize = "digital.shift", 
                                  skip = skip)
                   },
                   "ghalton" = {
                      qrng::ghalton(n, d = dim., method = "generalized")
                   },
                   "PRNG" = {
                      matrix(runif(n* dim.), ncol = dim.)
                   })  # (n, dim.) matrix
       ## Get quasi realizations of W via 'mix_()' (quantile function here)
       mix_(U[, 1])
    } else if(rmix.sample){
       ## A sample was provided
       rmix
    } else {
       ## Call provided RNG 
       mix_(n)
    } 
    
    ## Generate normals or gamma variates 
    if(which == "nvmix") {
        ## Generate Z ~ N(0, I_k)
        Z <- if(!inversion) {
                 matrix(rnorm(n * k), ncol = k) # (n, k)-matrix of N(0, 1)
             } else {
                 qnorm(U[, 2:(k+1)]) # (n, k)-matrix of N(0, 1)
             }
        ## Generate Y ~ N_d(0, scale)
        ## Recall that factor had been transposed, i.e. factor is (k,d)
        Y <- Z %*% factor # (n, k) %*% (k, d) = (n, d)-matrix of N(0, scale)
        ## Generate X ~ NVM_d(0, Sigma, F_W)
        X <- sqrt(W) * Y # also fine for different k
        ## Generate X ~ NVM_d(mu, Sigma, F_W)
        sweep(X, 2, loc, "+")
    } else {
        ## Generate Z^2 ~ chi^2_d
        Zsq <- if(!inversion) {
                   rgamma(n, shape = d/2, scale = 2)
               } else {
                   qgamma(U[, 2], shape = d/2, scale = 2)
               }
        W * Zsq
    }
}


### rnvmix() ###################################################################

##' @title Random Number Generator for Multivariate Normal Variance Mixtures
##' @param n sample size
##' @param rmix specification of random number generator of the  (mixture) distribution
##'        of W. This can be:
##'        1) character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "r", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) function being interpreted as a random number generator of W.
##'           Additional arguments can be passed via '...'
##'        4) n-vector containing a random sample from W.
##' @param qmix specification of the quantile function of the  (mixture) distribution
##'        of W. This needs to be supplied for the methods "sobol" and "ghalton".This can be:
##'        1) character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "qmix" random number generator.
##'        3) function being interpreted as the quantile function F_W^-.
##'           Additional arguments can be passed via '...'
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param factor (d, k)-matrix such that factor %*% t(factor) = scale;
##'        internally determined via chol() (and then an upper triangular
##'        matrix) if not provided
##' @param method character string indicating the method to be used:
##'         - "PRNG":    pure Monte Carlo
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         Note: For the methods "sobol" and "ghalotn", qmix() must be provided
##'         and rmix() is ignored. For the method "PRNG", either qmix() or rmix()
##'         needs to be provided. If both are provided, qmix() is ignored and
##'         rmix() is used.
##' @param skip numeric integer. How many points should be skipped when method='sobol'?
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return (n, d)-matrix with NVM(loc,scale, F_W) samples if 'which == "nvmix"'
##'         or n-vector with samples of the squared mahalanobis distance of
##'         NVM(loc,scale, F_W)
##' @author Marius Hofert and Erik Hintz

rnvmix <- function(n, rmix, qmix, loc = rep(0, d), scale = diag(2),
                   factor = NULL, method = c("PRNG", "sobol", "ghalton"),
                   skip = 0, ...)
{
    ## Get 'd':
    d <- if(is.null(factor)) dim(scale)[1] else nrow(factor <- as.matrix(factor))
    ## Call internal 'rnvmix_()' 
    rnvmix_(n, rmix = rmix, qmix = qmix, loc = loc, scale = scale,
            factor = factor, method = method, skip = skip, which = "nvmix", ...)
}