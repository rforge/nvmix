### dnvmix() ###################################################################
### I changed the code to account for general nvmix distributions. The old code is commented out below.


##' @title Density of Multivariate Normal Variance Mixtures
##' @param x (n, d)-matrix of evaluation points
##' @param loc location vector of dimension d
##' @param scale covariance matrix of dimension (d, d)
##' @param mix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param factor factorization matrix of the covariance matrix scale;
##'        caution: this has to be an *upper triangular* matrix R
##'        such that R^T R = scale here (otherwise det(scale) not computed correctly)
##' @param abserr numeric and non-negative. Absolute precision required. If abserr = 0, algorithm will run until total number of function evaluations exceeds Nmax (see also Nmax)
##' @param gam Monte Carlo confidence multiplier. Algorithm runs until gam * (estimated standard error) < abserr. gam = 3.3 means that one can expect
##'        that in 99.9% of the cases the actual absolute error is less than abserr.
##' @param Nmax Total number of function evaluations allowed.
##' @param B Number of randomizations to get error estimate.
##' @param n_init First loop uses n_init function evaluations. Any positive integer allowed, powers or at least multiples of 2 are recommended.
##' @param precond Logical. If TRUE (recommended), variable reordering as described in [genzbretz2002] pp. 955-956 is performed. Variable reordering can lead to a significant variance
##'        reduction and decrease in computational time.
##' @param method Character string indicating method to be used. Allowed are
##'         - "sobol" for a Sobol sequence.
##'         - "ghalton" for a generalized Halton sequence.
##'         - "prng" for a pure Monte Carlo approach.
##' @return n-vector with NVM(loc,scale,mix) density values
##' @author Marius Hofert and Erik Hintz
dnvmix <- function(x, loc = rep(0, d), scale, mix,
                   factor = tryCatch(factorize(scale), error = function(e) e), # needs to be triangular!
                   abserr = 0.001, gam = 3.3, Nmax = 1e8, B = 12, n_init = 2^6, method = "sobol",
                   log = FALSE, ...)
{
  ## Logicals if we are dealing with a multivariate normal or multivariate t
  const <- FALSE
  inv.gam <- FALSE
  
  ## Define the quantile function of the mixing variable.
  ## If mix is either "constant" or "inverse.gamma", we don't need the quantile function as there is a closed
  ## formula for the density in these cases.
  if(is.character(mix)){
    mix <- match.arg(mix, choices = c("constant", "inverse.gamma"))
    switch(mix,
           "constant" = {
             const <- TRUE
           },
           "inverse.gamma" = {
             if(hasArg(df)) df <- list(...)$df else stop("'mix = \"inverse.gamma\"' requires 'df' to be provided.")
             ## Still allow df = Inf (normal distribution)
             stopifnot(is.numeric(df), length(df) == 1, df > 0)
             if(is.finite(df)) {
               inv.gam <- TRUE
             } else {
               const <- TRUE
             }
           },
           stop("Currently unsupported 'mix'"))
  } else {
    ## In all the other cases we do need a quantile function, which we define now:
    W <- if(is.list(mix)) { # 'mix' is a list of the form (<character string>, <parameters>)
      stopifnot(length(mix) >= 1, is.character(distr <- mix[[1]]))
      qmix <- paste0("q", distr)
      if(!existsFunction(qmix))
        stop("No function named '", qmix, "'.")
      function(u){
        return(  do.call(qmix, c(u, mix[-1])))
      }
    } else if(is.function(mix)) { # 'mix' is interpreted as the quantile function F_W^- of the mixture distribution F_W of W
      function(u){
        return(mix(u, ...))
      }
    } else stop("'mix' must be a character string, list or quantile function.")
  }
  
  if(!is.matrix(x)) x <- rbind(x)
  n <- nrow(x)
  d <- ncol(x)
  stopifnot(length(loc) == d)
  notNA <- apply(!is.na(x), 1, all)
  lres <- rep(-Inf, n)
  lres[!notNA] <- NA
  x <- x[notNA,] # available points
  tx <- t(x) # (d, n)-matrix
  if(inherits(factor, "error") || is.null(factor)) {
    lres[notNA & (colSums(tx == loc) == d)] <- Inf
  } else {
    ## Solve R^T * z = x - mu for z, so z = (R^T)^{-1} * (x - mu) (a (d, d)-matrix)
    ## => z^2 (=> componentwise) = z^T z = (x - mu)^T * ((R^T)^{-1})^T (R^T)^{-1} (x - mu)
    ##                           = z^T z = (x - mu)^T * R^{-1} (R^T)^{-1} (x - mu)
    ##                           = (x - mu)^T * (R^T R)^{-1} * (x - mu)
    ##                           = (x - mu)^T * scale^{-1} * (x - mu) = quadratic form
    z <- backsolve(factor, tx - loc, transpose = TRUE)
    maha2 <- colSums(z^2) # = sum(z^T z); squared Mahalanobis distance from x to mu w.r.t. scale
    ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
    ## = log(prod(diag(R))) = sum(log(diag(R)))
    lrdet <- sum(log(diag(factor)))
    ## First we catch the case of a mulitvariate t/multivariate normal
    if(inv.gam){
      df.d.2 <- (df + d) / 2
      lres[notNA] <- lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) - lrdet - df.d.2 * log1p(maha2 / df)
    } else if(const){
      lres[notNA] <- -(d/2) * log(2 * pi) - lrdet - maha2/2
    } else{
      
      ## If we are not dealing with a multivariate t/normal, we use a RQMC procedure as in pnvmix
      
      gam <- gam / sqrt(B) # instead of dividing sigma by sqrt(B) each time
      n. <- n_init # initial n
      T. <- matrix(0, ncol = n, nrow = B) # matrix to store RQMC estimates
      err <- abserr + 42 # initialize err to something bigger than abserr so that we can enter the while loop
      N. <- 0 # N. will count the total number of function evaluations
      i. <- 0 # initialize counter; this will count the number of iterations in the while loop
      useskip <- 0 # will need that because the first iteration is a little different from all the others
      denom <- 1
      
      maha2 <- as.matrix(maha2, nrow = 1, ncol = n) # 1 x n matrix 
      
      if(method == "sobol") {
        if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
        seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used.
      }
      
      while(err > abserr && N. < Nmax)
      {
        
        if(method == "sobol") .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
        
        ## Get B RQCM estimates
        for(l in 1:B){
          
          # Get the pointset
          U <- switch(method,
                      "sobol"   = {
                        qrng::sobol(n = n., d = 1, randomize = TRUE, skip = (useskip * n.))
                      },
                      "gHalton" = {
                        qrng::ghalton(n = n., d = 1, method = "generalized")
                      },
                      "prng"    = {
                        matrix(runif( n. ), ncol = 1)
                      })


          U <- as.matrix( W(U), ncol = 1, nrow = n.)

          b <- - (d/2) * matrix(log(2 * pi * U), nrow = n., ncol = n) - lrdet - (1/U) %*% t(maha2) / 2 #n. x n matrix, each column corresponds to "one x"
          
          bmax <- apply(b, 2, max) # n vector
          
          T.[l,] <- ( T.[l,] - log(n.) + bmax + log( colSums( exp (b - matrix(bmax, ncol = n, nrow = n., byrow = TRUE) ) ) ) )/denom
        }
          
        ## Update the total number of function evaluations
        N. <- N. + B * n.
        
        ## Change denom and useksip. This is done exactly once, namely in the first iteration.
        if(i. == 0){
          
          denom <- 2
          useskip <- 1
          
        } else {
          
          ## Increase sample size n. This is done in all iterations except for the first two.
          n. <- 2 * n.
          
        }  
        
        sig <- max( apply(T., 2, sd)) # get standard deviation of the column with the largest standard deviation
        err <- gam * sig # update error. Note that this gam is actually gamma/sqrt(N)
        i. <- i. + 1 # update counter 
      }
      lres[notNA] <- apply(T., 2, mean)
    }
  }
  if(log) lres else exp(lres) # also works with NA, -Inf, Inf
}




## OLD CODE:
# dnvmix <- function(x, df, loc = rep(0, d), scale,
#                    factor = tryCatch(factorize(scale), error = function(e) e), # needs to be triangular!
#                    log = FALSE)
# {
#   if(!is.matrix(x)) x <- rbind(x)
#   n <- nrow(x)
#   d <- ncol(x)
#   stopifnot(df > 0, length(loc) == d)
#   notNA <- apply(!is.na(x), 1, all)
#   lres <- rep(-Inf, n)
#   lres[!notNA] <- NA
#   x <- x[notNA,] # available points
#   tx <- t(x) # (d, n)-matrix
#   if(inherits(factor, "error") || is.null(factor)) {
#     lres[notNA & (colSums(tx == loc) == d)] <- Inf
#   } else {
#     ## Solve R^T * z = x - mu for z, so z = (R^T)^{-1} * (x - mu) (a (d, d)-matrix)
#     ## => z^2 (=> componentwise) = z^T z = (x - mu)^T * ((R^T)^{-1})^T (R^T)^{-1} (x - mu)
#     ##                           = z^T z = (x - mu)^T * R^{-1} (R^T)^{-1} (x - mu)
#     ##                           = (x - mu)^T * (R^T R)^{-1} * (x - mu)
#     ##                           = (x - mu)^T * scale^{-1} * (x - mu) = quadratic form
#     z <- backsolve(factor, tx - loc, transpose = TRUE)
#     maha2 <- colSums(z^2) # = sum(z^T z); squared Mahalanobis distance from x to mu w.r.t. scale
#     ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
#     ## = log(prod(diag(R))) = sum(log(diag(R)))
#     lrdet <- sum(log(diag(factor)))
#     lres[notNA] <- if(is.finite(df)) {
#       df.d.2 <- (df + d) / 2
#       lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) - lrdet - df.d.2 * log1p(maha2 / df)
#     } else {
#       - (d/2) * log(2 * pi) - lrdet - maha2/2
#     }
#   }
#   if(log) lres else exp(lres) # also works with NA, -Inf, Inf
# }
