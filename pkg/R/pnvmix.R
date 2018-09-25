### pnvmix() ###################################################################

##' @title Re-order Variables According to their Expected Integration Limits
##'        (Precondition) for the case of a multivariate t distribution. See [genzbretz2002, p. 957].
##' @param a vector, lower integration limit
##' @param b vector, upper integration limit
##' @param C Cholesky (lower triangular) factor of the sccale matrix R
##' @param R scale matrix
##' @param q dimension of the problem
##' @param meansqrtmix mean of sqrt(W)
##' @return list a, b, R, C of integration limits/scale matrix and Cholesky factor after reordering has been performed
##' @author Erik Hintz
precond <- function(a, b, R, C, q, meansqrtmix)
{
  y <- rep(0, q - 1)

  ## Find i = argmin_j { <expected length of interval> }
  i <- which.min( apply( pnorm( cbind(b, a) / (meansqrtmix * sqrt(diag(R)) ) ), 1, diff) )

  if(i != 1) {
    ## Swap 1 and i
    tmp <- swap(a = a, b = b, R = R, i = 1, j = i)
    a <- tmp$a
    b <- tmp$b
    R <- tmp$R
  }

  ## Store y1
  y[1] <- - ( dnorm(b[1]/meansqrtmix) - dnorm(a[1]/meansqrtmix)) / (pnorm(b[1]/meansqrtmix) - pnorm(a[1]/meansqrtmix))

  ## Update Cholesky
  C[1, 1] <- sqrt( R[1, 1] )
  C[2:q, 1] <- as.matrix( R[2:q, 1] / C[1, 1] )

  for(j in 2:(q-1)){

    denom <- sqrt( diag(R)[j:q] - rowSums( as.matrix( C[j:q, 1:(j-1)] )^2 ) )
    c <- as.matrix( C[j:q, 1:j-1] ) %*% y[1:(j-1)]

    # Find i = argmin { <expected length of interval j> }
    i <- which.min( pnorm( (b[j:q] / meansqrtmix - c) / denom ) - pnorm( (a[j:q] / meansqrtmix - c) / denom ) ) + j - 1

    if(i != j){
      ## Swap i and j
      tmp <- swap(a = a, b = b, R = R, i = i, j = j)
      a <- tmp$a
      b <- tmp$b
      R <- tmp$R

      C[c(i,j), ]    <- as.matrix(C[c(j,i), ])
      C[j, (j+1):i] <- as.matrix(0, ncol = i - j, nrow = 1)
    }

    ## Update Cholesky
    C[j,j] <- sqrt(R[j,j] - sum( C[j, 1:(j-1)]^2 ) )
    if(j< (q-1)) C[(j+1):q, j] <- ( R[(j+1):q, j] - as.matrix( C[(j+1):q, 1:(j-1)] ) %*% C[j, 1:(j-1)] ) / C[j, j]
    else C[(j+1):q, j] <- ( R[(j+1):q, j] - C[(j+1):q, 1:(j-1)] %*% C[j, 1:(j-1)] ) / C[j, j]

    ## Get yj
    ajbj <-  c( (a[j] / meansqrtmix - C[j, 1:(j-1)] %*% y[1:(j-1)] ) , (b[j] / meansqrtmix - C[j, 1:(j-1)] %*% y[1:(j-1)]) ) / C[j,j]
    y[j] <- ( dnorm(ajbj[1]) - dnorm(ajbj[2]) ) / ( pnorm(ajbj[2]) - pnorm(ajbj[1]) )
  }

  C[q,q] <- sqrt(R[q, q] - sum(C[q, 1:(q-1)]^2))
  list(a = a, b = b, R = R, C = C)
}

##' @title Distribution Function of the Multivariate t Distribution
##' @param upper vector of upper limits
##' @param lower vector of lower limits
##' @param loc location vector (corresponds to mu in the paper)
##' @param scale covariance matrix of dimension (q, q)
##' @param mix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param meansqrtmix Mean of sqrt(W)
##' @param standardized logical. If TRUE, scale is assumed to be a correlation matrix; if FALSE, lower,upper and scale will be normalized
##' @param ... additional arguments containing parameters of
##'        mixing distributions
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
##' @author Erik Hintz
pnvmix <- function(upper, lower = rep(-Inf, length(upper)), loc = rep(0, length(upper)), scale, mix, meansqrtmix = NA, standardized = FALSE,
                     gam = 3.3, abserr = 0.001, Nmax = 1e8, B = 12, n_init = 2^6, precond = TRUE, method = "sobol", ... )
{


  if(!is.matrix(scale)) scale <- as.matrix(scale)
  q <- dim(scale)[1] # dimension of the problem

  ## Checks
  if( method != "sobol" && method != "ghalton" && method != "prng") stop("Only sobol, ghalton and prng are allowed as methods")
  if( length(lower) != length(upper) ) stop("Lenghts of lower and upper differ")
  if( any(lower >= upper) ) stop("lower needs to be smaller than upper (componentwise)")
  if( q != length(lower) ) stop("Dimension of scale does not match dimension of lower")
  if( q != length(loc) ) stop("Dimension of loc does not match dimension of scale")

  ## Find infinite limits
  infina  <-  (lower == -Inf)
  whichfina <- which(!infina)
  infinb  <-  (upper == Inf)
  whichfinb <- which(!infinb)
  infinab <-  infina * infinb

  ## Remove double infinities
  if( sum(infinab) > 0 )
  {
    whichdoubleinf <- which(infinab == 1) # indicides with limits [-Inf, Inf]
    lower <- lower[ -whichdoubleinf ]
    upper <- upper[ -whichdoubleinf ]
    scale <- scale[ -whichdoubleinf, -whichdoubleinf ]

    ## If all but one element is deleted, still want scale to be a matrix
    if(!is.matrix(scale)) scale <- as.matrix(scale)

    ## Update dimension
    q <- dim(scale)[1]
  }

  ## Subtract shift if necessary:
  if(any(loc != 0)){
    lower <- lower - loc
    upper <- upper - loc
  }

  ## Standardize if necessary:
  if(!standardized){
    Dinv <- diag(1/sqrt(diag(scale)))
    scale <- Dinv %*% scale %*% Dinv
    lower[whichfina] <- as.vector(Dinv[whichfina, whichfina] %*% lower[whichfina]) # only works for those values which are not +/- Inf
    upper[whichfinb] <- as.vector(Dinv[whichfinb, whichfinb] %*% upper[whichfinb])
  }

  ## Logicals if we are dealing with a multivariate normal or multivariate t
  const <- FALSE
  inv.gam <- FALSE

  ## Define the quantile function of the mixing variable:
  W <- if(is.character(mix)) { # 'mix' is a character vector specifying supported mixture distributions (utilizing '...')
    mix <- match.arg(mix, choices = c("constant", "inverse.gamma"))
    switch(mix,
           "constant" = {
             const <- TRUE
             function(u){
               return(1)
             }
           },
           "inverse.gamma" = {
             if(hasArg(df)) df <- list(...)$df else stop("'mix = \"inverse.gamma\"' requires 'df' to be provided.")
             ## Still allow df = Inf (normal distribution)
             stopifnot(is.numeric(df), length(df) == 1, df > 0)
             if(is.finite(df)) {
               inv.gam <- TRUE
               df2 <- df / 2
               meansqrtmix <- sqrt(df) * gamma(df2) / ( sqrt(2) * gamma( (df+1) / 2 ) ) # mean of sqrt(W) in this case, will be used for preconditioning
               function(u){
                 return(1 / qgamma(u, shape = df2, rate = df2))
               }
             } else {
               const <- TRUE
               meansqrtmix <- 1 # mean of sqrt(W) in this case, will be used for preconditioning
               function(u){
                 return(1)
               }
             }
           },
           stop("Currently unsupported 'mix'"))
  } else if(is.list(mix)) { # 'mix' is a list of the form (<character string>, <parameters>)
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


  ## If q = 1 and const or mvt is true, we are dealing with a univariate normal or t distribution
  ## These cases can be dealt with pnorm() and pt()

  if(q == 1){

    if(const){

      Prob <- pnorm(upper) - pnorm(lower)
      return(list(Prob = Prob, N = 0, i = 0, ErrEst = 0, Var = 0))
    }

    if(inv.gam){
      Prob <- pt(upper, df = df) - pt(lower, df = df)
      return(list(Prob = Prob, N = 0, i = 0, ErrEst = 0, Var = 0))
    }

  }

  ## Get Cholesky factor (lower triangular)
  C <- t(chol(scale))

  ## Precondtioning (resorting the limits (cf precond); only for q > 2)
  if(precond && q>2) {

    if(is.na(meansqrtmix)){
      ## If meansqrtmix was not supplied, we approximate it
      meansqrtmix <- mean(sqrt(W(qrng::sobol(n = 5000, d = 1, randomize = TRUE))))

    } else if(meansqrtmix <= 0) {

      stop("Meansqrtmix has to be positive.")
    }

    temp <- precond(a = lower, b = upper, R = scale, C = C, q = q, meansqrtmix = meansqrtmix)
    lower <- temp$a
    upper <- temp$b
    C <- temp$C
  }

  gam <- gam / sqrt(B) # instead of dividing sigma by sqrt(B) each time
  n. <- n_init # initial n
  T. <- rep(0, B) # vector to store RQMC estimates

  ONE <- 1-.Machine$double.neg.eps
  ZERO <- .Machine$double.eps

  if(method == "sobol") {
      if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
      seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used.
  }

  err <- abserr + 42 # initialize err to something bigger than abserr so that we can enter the while loop
  N. <- 0 # N. will count the total number of function evaluations
  i. <- 0 # initialize counter; this will count the number of iterations in the while loop

  # useskip <- 0 # will need that because the first iteration is a little different from all the others
  # denom <- 1

  while(err > abserr && N. < Nmax)
  {

    if(method == "sobol") .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )

    ## Get N RQCM estimates
    for(l in 1:B){

      ## Get the pointset
      ## If const = TRUE, we only need q - 1 (quasi) random numbers
      ## The case const = TRUE and q = 1 has already been dealt with

      if(const){
        U <- switch(method,
                    "sobol"   = {
                      # qrng::sobol(n = n., d = q - 1, randomize = TRUE, skip = (useskip * n.))
                      qrng::sobol(n = n., d = q - 1, randomize = TRUE, skip = (i. * n.) )
                    },
                    "ghalton" = {
                      qrng::ghalton(n = n., d = q - 1, method = "generalized")
                    },
                    "prng"    = {
                      matrix(runif( n. * (q - 1)), ncol = q - 1)
                    })

        ## First and last column contain 1s corresponding to simulated values from sqrt(mix)
        U <- cbind( rep(1, n.), U, rep(1, n.))

      } else {

        U <- switch(method,
                    "sobol"   = {
                      # qrng::sobol(n = n., d = q, randomize = TRUE, skip = (useskip * n.))
                      qrng::sobol(n = n., d = q, randomize = TRUE, skip = (i. * n.))
                    },
                    "ghalton" = {
                      qrng::ghalton(n = n., d = q, method = "generalized")
                    },
                    "prng"    = {
                      matrix(runif( n. * q), ncol = q)
                    })

        ## Case q = 1 somewhat special again:

        if( q == 1){

          U <- cbind( sqrt( W(U) ), sqrt( W( 1 - U) ))

        } else {

          ## Column 1: sqrt(mix), Columns 2 - q: unchanged (still uniforms), Column q + 1: anitithetic realization of sqrt(mix)

          U <- cbind(sqrt( W(U[, 1]) ), U[,2:q], sqrt( W( 1 - U[, 1]) ))
        }
      }

      ## Evaluate the integrand at the point-set and save it in T[] (calls C Code).
      ## Both, T.[l] and the new estimate are based on n. evaluations, so we can just average them unless we are in the first iteration.
      ## In this case, denom = 1 and T.[l] = 0.

      if(q == 1){

        ## Case of dimension 1 seperate: Here, we do not need to approximate the multivariate normal cdf and can just use pnorm
        ## The case of dimension 1 for a normal / t distribution has already been dealt with

        #T.[l] <- (T.[l] + mean( ( pnorm( upper / U[,1]) - pnorm( lower / U[,1]) + pnorm( upper / U[, q + 1]) - pnorm( lower / U[, q + 1])) / 2) ) / denom
        T.[l] <- (i. * T.[l] + mean( ( pnorm( upper / U[,1]) - pnorm( lower / U[,1]) + pnorm( upper / U[, q + 1]) - pnorm( lower / U[, q + 1])) / 2) ) / (i. + 1)

      } else {

        ## TODO: why keep this code? give a reason
        ## T.[l] <- (T.[l] + .Call("eval_nvmix_integral",
        ##                a    = as.double(lower),
        ##                b    = as.double(upper),
        ##                U    = as.double(U),
        ##                n    = as.integer(n.),
        ##                q    = as.integer(q),
        ##                C    = as.double(C),
        ##                ZERO = as.double(ZERO),
        ##                ONE  = as.double(ONE)) )/denom

          T.[l] <- (i. * T.[l] + .Call("eval_nvmix_integral",
                                a    = as.double(lower),
                                b    = as.double(upper),
                                U    = as.double(U),
                                n    = as.integer(n.),
                                q    = as.integer(q),
                                C    = as.double(C),
                                ZERO = as.double(ZERO),
                                ONE  = as.double(ONE)) ) / (i. + 1)

      }
    } # end for(l in 1:N)

    ## Update the total number of function evaluations; mutliplied by 2 since antithetic variates are being used in eval_nvmix_integral
    N. <- N. + 2 * B * n.

    ## TODO why keep this code? give a reason
    ## ## Change denom and useksip. This is done exactly once, namely in the first iteration.
    ## if(i. == 0){
    ##
    ##   denom <- 2
    ##   useskip <- 1
    ##
    ## } else {
    ##
    ##   ## Increase sample size n. This is done in all iterations except for the first two.
    ##   n. <- 2 * n.
    ##
    ## }

    sig <- sd(T.) # get standard deviation of the estimator
    err <- gam * sig # update error. Note that this gam is actually gamma/sqrt(N)
    i. <- i. + 1 # update counter
  }

  ## Calculate the RQMC estimator:
  T <- mean(T.)

  ## Get the variance
  var <- (sig/sqrt(B))^2

  ## Print a warning if precision level not reached:
  if(err > abserr) warning("Precision level abserr not reached; consider increasing Nmax.")

  ## And return
  list(Prob = T, N = N., i = i., ErrEst = err, Var = var)
}

