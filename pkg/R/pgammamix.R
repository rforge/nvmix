# ### pgammamix() ###################################################################

##' @title Distribution function of a chi-sq.-mixture distribution
##' @param m n-vector of evaluation points
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param lower.tail logical        
##' @param d degree-of-freedom parameter 
##' @param control list of control arguments; see ?pnvmix 
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return n-vector with computed cdf values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert

pgammamix <- function(m, qmix, d, lower.tail = TRUE, 
                      control = list(), verbose = TRUE, ...)
{
  ## Checks
  if(!is.vector(m)) m <- as.vector(m) 
  n <- length(m) # length of input
  ## Deal with algorithm parameters, see also get.set.parameters():
  ## get.set.parameters() also does argument checking, so not needed here.
  control <- nvmix:::get.set.parameters(control)
  ## 1 Define the quantile function of the mixing variable ###################
  qW <- if(is.character(qmix)) { # 'qmix' is a character vector
    qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
    switch(qmix,
           "constant" = {
             function(u) 1
           },
           "inverse.gamma" = {
             if(hasArg(df)) {
               df <- list(...)$df
             } else {
               stop("'qmix = \"inverse.gamma\"' requires 'df' to be provided.")
             }
             ## Still allow df = Inf (normal distribution)
             stopifnot(is.numeric(df), length(df) == 1, df > 0)
             if(is.finite(df)) {
               df2 <- df / 2
               function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
             } else {
               function(u) 1
             }
           },
           stop("Currently unsupported 'qmix'"))
  } else if(is.list(qmix)) { # 'mix' is a list of the form (<character string>, <parameters>)
    stopifnot(length(qmix) >= 1, is.character(distr <- qmix[[1]]))
    qmix. <- paste0("q", distr)
    if(!existsFunction(qmix.))
      stop("No function named '", qmix., "'.")
    function(u)
      do.call(qmix., append(list(u), qmix[-1]))
  } else if(is.function(qmix)) { # 'mix' is the quantile function F_W^- of F_W
    function(u)
      qmix(u, ...)
  } else stop("'qmix' must be a character string, list or quantile function.")
  ## Build result object 
  pres <- rep(0, n) # n-vector of results
  notNA <- which(!is.na(m)) 
  pres[!notNA] <- NA
  m <- m[notNA, drop = FALSE] # non-missing data (rows)
    ## Counter
  numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
  ## 1 Basics ##################################################################
  ## Define various quantites:
  dblng           <- (control$increment == "doubling")
  B               <- control$B # number of randomizations
  current.n       <- control$fun.eval[1] #initial sample size
  numiter         <- 0 # counter for the number of iterations
  total.fun.evals <- 0
  ZERO            <- .Machine$double.neg.eps
  ## Absolte/relative precision?
  if(is.na(control$pgammamix.reltol)){
    ## Use absolute error
    tol <- control$pgammamix.abstol
    do.reltol <- FALSE
  } else {
    ## Use relative error
    tol <- control$pgammamix.reltol
    do.reltol <- TRUE
  }
  ## Store seed if 'sobol' is used to get the same shifts later:
  if(control$method == "sobol") {
    if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
    seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
  }
  ## Additional variables needed if the increment chosen is "dblng"
  if(dblng) {
    if(control$method == "sobol") useskip <- 0
    denom <- 1
  }
  ## Matrix to store RQMC estimates
  rqmc.estimates <- matrix(0, ncol = n, nrow = B)
  ## Will be needed a lot:
  CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
  ## Initialize 'max.error' to > tol so that we can enter the while loop:
  max.error <- tol + 42
  ## 2 Main loop ###############################################################
  ## while() runs until precision abstol is reached or the number of function
  ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
  ## the desired log-densities are calculated.
  while(max.error > tol && numiter < control$max.iter.rqmc &&
        total.fun.evals < control$fun.eval[2])
  {
    ## Reset seed to have the same shifts in sobol( ... )
    if(control$method == "sobol" && numiter > 0)
      .Random.seed <<- seed # reset seed to have the same shifts in sobol( ... )
    for(b in 1:B){
      ## 2.1 Get the point set ###########################################
      U <- switch(control$method,
                       "sobol" = {
                         if(dblng) {
                           qrng::sobol(n = current.n, d = 1,
                                       randomize = TRUE,
                                       skip = (useskip * current.n))
                         } else {
                           qrng::sobol(n = current.n, d = 1,
                                       randomize = TRUE,
                                       skip = (numiter * current.n))
                         }
                       },
                       "ghalton" = {
                         qrng::ghalton(n = current.n, d = 1,
                                       method = "generalized")
                       },
                       "PRNG" = {
                         runif(current.n)
                       }) # sorted for later!
      ## 2.2 Evaluate the integrand at the (next) point set #############
      W <- pmax(qW(U), ZERO) # realizations of the mixing variable
      next.estimate <- .colMeans(pchisq( 
        sapply(m, function(i) i/W), df = d, lower.tail = lower.tail), 
        current.n, n, 0)
      ## 2.3 Update RQMC estimates #######################################
      rqmc.estimates[b,] <-
        if(dblng) {
          ## In this case both, rqmc.estimates[b,] and
          ## next.estimate depend on n.current points
          (rqmc.estimates[b,] + next.estimate)/denom
        } else {
          ## In this case, rqmc.estimates[b,] depends on
          ## numiter * n.current points whereas next.estimate
          ## depends on n.current points
          (numiter * rqmc.estimates[b,] + next.estimate)/(numiter + 1)
        }
    } # end for(b in 1:B)
    ## Update of various variables
    ## Double sample size and adjust denominator in averaging as well as useskip
    if(dblng) {
      ## Change denom and useksip (exactly once, in the first iteration)
      if(numiter == 0){
        denom <- 2
        useskip <- 1
      } else {
        ## Increase sample size n. This is done in all iterations
        ## except for the first two
        current.n <- 2 * current.n
      }
    }
    ## Total number of function evaluations:
    total.fun.evals <- total.fun.evals + B * current.n
    numiter <- numiter + 1
    ## Update error. The following is slightly faster than 'apply(..., 2, var)' 
    pres <- .colMeans(rqmc.estimates, B, n , 0)
    vars <- .colMeans((rqmc.estimates - rep(pres, each = B))^2, B, n, 0)
    errors <- if(!do.reltol){
      sqrt(vars)*CI.factor.sqrt.B
    } else {
      sqrt(vars)/pres*CI.factor.sqrt.B
    }
    max.error <- max(errors)
  } # while()
  if(verbose && max.error > tol)
    warning("Tolerance not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
  ## 4 Return ##################################################################
  attr(pres, "error")   <- errors 
  attr(pres, "numiter") <- numiter
  pres
}
  
  
 