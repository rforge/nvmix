##' @title Extrapolate log-densitiy of a NVM dist'n 
##' @param maha2.2 vector of squared mahalanobis distances (sorted)
##' @param ldensities log-densities of NVM density evaluated at maha2.2
##' @param errors estimated errors of ldensities (from dnvmix.int)
##' @param control list of control arguments, see ?dnvmix
##' @author Erik Hintz

dnvmix.extrapolate <- function(maha2.2, ldensities, errors, control){
  n <- length(maha2.2)
  ## logical vector if precision reached
  notRchd <- (errors >= control$dnvmix.abstol)
  ## First not reached in the *tail* (ignore first ones, if any)
  fnotrchd <- if(!notRchd[1]) which.max(notRchd)[1] else {
    if(prod(notRchd)) 1 else {
      firstFalse <- which.min(notRchd)[1]
      firstFalse - 1 + which.max(notRchd[firstFalse:n]) 
    }
  }
  ## Indices used to fit:
  ind.fit <- if(fnotrchd - control$dnvmix.extrap.num.fit > 1){
    ## Do we have enough to fit?
    (fnotrchd - control$dnvmix.extrap.num.fit - 1):fnotrchd
  } else {
    ## If not, just take the first 'dnvmix.extrap.num.fit' ones
    1:control$dnvmix.extrap.num.fit 
  }
  ## Indices used to test:
  ind.test <- if(fnotrchd + control$dnvmix.extrap.num.test - 1 <= n){
    ## Do we have enough to test?
    fnotrchd:(fnotrchd + control$dnvmix.extrap.num.test - 1)
  } else {
    ## If not, just take the 'dnvmix.extrap.num.test' ones after the fitting ones
    (n - control$dnvmix.extrap.num.test + 1):n
  }
  
  ## Maha distances + logdensities used to fit:
  maha.fit  <- maha2.2[ind.fit]
  ldens.fit <- ldensities[ind.fit]
  ## ldensities used to test:
  ldens.test <- maha2.2[ind.test]
  
  ## Set up optimizer function:
  log.approx <- function(maha, dens, linear = FALSE){
    err <- if(linear){
      start <- c(dens[1], -1)
      low <- c(-Inf, -Inf)
      upp <- c(Inf, 0)
      function(coeff){
        sum((coeff[1] + coeff[2]*maha - dens)^2) 
      }
    } else {
      start <- c(dens[1], 0, 0, 1, -1, 1)
      low <- c(-Inf, -Inf, -1, 0 , -Inf, -Inf)
      upp <- c(Inf, 0, Inf, Inf, Inf, Inf)
      function(coeff){
        sum((coeff[1] + coeff[2]*log(1 + coeff[3] + coeff[4]*maha)
             + coeff[5]*maha^coeff[6] - dens)^2) 
      } 
    }
    o <- optim( start,
                err,
                lower = low,
                upper = upp,
                method = c("L-BFGS-B"),
                control = list(maxit = 500))
    o$par
  }
  
  ## Call optimizer, once linear, once general:
  coeff_linear <- log.approx(maha.fit, ldens.fit, linear = TRUE)
  coeff_notLin <- log.approx(maha.fit, ldens.fit, linear = FALSE)
  
  ## Approximate densities:
  dens.approx_linear <- coeff_linear[1] + coeff_linear[2]*maha2.2
  dens.approx_notLin <- coeff_notLin[1] + coeff_notLin[2]*log(coeff_notLin[3] + coeff_notLin[4]*maha2.2) +
    coeff_notLin[5]*maha2.2^coeff_notLin[6] 
  
  ## Which one is better? 
  err <- c( sum(dens.approx_linear[ind.test] - ldens.test)^2, 
            sum(dens.approx_notLin[ind.test] - ldens.test)^2)
  
  ## Return correspondingly
  ## Set errors for extrapolated values to 'NA' 
  errors[fnotrchd:n] <- NA 
  
  ldensities[fnotrchd:n] <-  if(err[1] < err[2]){
    dens.approx_linear[fnotrchd:n] 
  } else {
    dens.approx_notLin[fnotrchd:n]
  }
  ## Could probably return more later, e.g. the fitted parameters
  return(list(ldensities = ldensities, errors = errors))
}



### dnvmix() ###################################################################

##' @title Density of a Multivariate Normal Variance Mixture - Internal version 
##'        (not exported)
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lrdet log(sqrt(det(scale))) where 'scale' is the scale matrix of 
##'        the normal variance mixture distribution. 
##' @param U0 vector of first input uniforms. Length has to be multiple of B. 
##' @param d dimension of the Normal Variance Mixture        
##' @param method see details in ?pnvmix
##' @param increment see details in ?pnvmix
##' @param abstol see details in ?pnvmix
##' @param CI.factor see details in ?pnvmix
##' @param max.iter.rqmc see details in ?pnvmix
##' @param B see details in ?pnvmix
##' @param seed value of .Random.seed when U0 was generated; needed to get 
##'        the same shifts in the randomized Sobol approach. 
##' @return List of three:
##'         $ldensities n-vector with computed log-density values 
##'         $error error estimate (largest error in ldensities)
##'         $numiter number of iterations needed 
##' @author Erik Hintz and Marius Hofert
dnvmix.int <- function(qW, maha2.2, lrdet, U0, d, 
                       method = c("sobol", "ghalton", "PRNG"),
                       increment = "doubling",
                       abstol = 0.001, reltol = NA, CI.factor = 3.3, 
                       fun.eval = c(2^6, 1e8), max.iter.rqmc, B, seed)
{
  
  ## 1 Basics ##################################################################
  ## 
  ## Note: Most checking was done in dnvmix()
  n <- length(maha2.2) # sample size 
  numiter <- 1 # counter for the number of iterations (one iteration is done here)
  CI.factor.sqrt.B <- CI.factor / sqrt(B) # instead of dividing by sqrt(B) each time
  rqmc.estimates <- matrix(0, ncol = n, nrow = B) # matrix to store RQMC estimates
  
  ## Absolute or relative precision required?
  if(is.na(reltol)){
    ## Use absolute error
    tol <- abstol/CI.factor.sqrt.B
    do.reltol <- FALSE
  } else {
    ## Use relative error
    tol <- reltol/CI.factor.sqrt.B
    do.reltol <- TRUE
  }
  
  ## 2 First point-set #########################################################
  
  ## First pointset that was passed as vector: U has length B * current.n
  ## Realizations of l'th shift are elements (l-1)*current.n + (1:current.n)
  current.n <- length(U0)/B 
  W <- qW(U0)
  
  for(l in 1:B){
    ## Grab realizations corresponding to l'th shift and use exp-log trick 
    ## The underlying C function "eval_dnvmix_integrand" needs both, maha2_2 
    ## and W to be sorted in increasing order.     
    rqmc.estimates[l,] <- .Call("eval_dnvmix_integrand", 
                                W          = as.double(sort(W[(l-1)*current.n + (1:current.n)])),
                                maha2_2    = as.double(maha2.2),
                                current_n  = as.integer(current.n),
                                n          = as.integer(n),
                                d          = as.integer(d),
                                k          = as.integer(d), 
                                lrdet      = as.double(lrdet))
  }
  
  ## Error meausre. Note that 'tol' was divided by 'CI.factor' 
  errors <- if(!do.reltol){
    apply(rqmc.estimates, 2, sd)
  } else {
    apply(rqmc.estimates, 2, sd)/abs(.colMeans(rqmc.estimates, B, n, 0))
  }
  error <- max(errors)
  total.fun.evals <- B*current.n
  
  ## 3 Main loop ###############################################################
  while(error > tol && numiter < max.iter.rqmc && total.fun.evals < fun.eval[2]) {
    
    if(method == "sobol") .Random.seed <<- seed # reset seed to have the same shifts in sobol( ... )
    
    ## Get B RQCM estimates
    for(b in 1:B) {
      ## Get the point set
      U <- switch(method,
                  "sobol"   = {
                    qrng::sobol(current.n, d = 1, randomize = TRUE, skip = numiter*current.n)
                  },
                  "gHalton" = {
                    qrng::ghalton(current.n, d = 1, method = "generalized")
                  },
                  "prng"    = {
                    cbind(runif(current.n)) # 1-column matrix
                  })
      
      ## Exp-log trick 
      ## The underlying C function "eval_dnvmix_integrand" needs both, maha2_2 
      ## and W to be sorted in increasing order.     
      next.estimate <- .Call("eval_dnvmix_integrand", 
                             W          = as.double(sort(qW(U))),
                             maha2_2    = as.double(maha2.2),
                             current_n  = as.integer(current.n),
                             n          = as.integer(n),
                             d          = as.integer(d),
                             k          = as.integer(d), 
                             lrdet      = as.double(lrdet)) 
      
      rqmc.estimates[b,] <- 
        if(increment == "doubling"){
          ## In this case both, rqmc.estimates[b] and
          ## next.estimate depend on n.current points
          (rqmc.estimates[b, ] + next.estimate) / 2
        } else {
          ## In this case, rqmc.estimates[b] depends on
          ## numiter * n.current points whereas next.estimate
          ## depends on n.current points
          (numiter * rqmc.estimates[b, ] + next.estimate) / (numiter + 1)
        }
      
    } # end for(b in 1:B)
    
    ## Update various variables
    ## Number of function evaluations:
    total.fun.evals <- total.fun.evals + B*current.n
    ## Number of iterations:
    numiter <- numiter + 1 # update counter
    ## Next sample size:
    current.n <- if(increment == "doubling") 2*current.n else current.n 
    ## Error measure
    errors <- if(!do.reltol){
      apply(rqmc.estimates, 2, sd)
    } else {
      apply(rqmc.estimates, 2, sd)/abs(.colMeans(rqmc.estimates, B, n, 0))
    }
    error <- max(errors)
  } # while()
  
  ## 4 Return ##################################################################
  ## Get estimates for log-density:
  lres <- .colMeans(rqmc.estimates, B, n, 0)
  ## Get correct error estimate:
  errors <- errors*CI.factor.sqrt.B
  return(list(ldensities = lres, numiter = numiter, error = errors))
}

##' @title Density of a Multivariate Normal Variance Mixture
##' @param x (n, d)-matrix of evaluation points
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param loc d-vector (location vector)
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param factor Cholesky factor (lower triangular matrix) of 'scale';
##'        important here so that det(scale) is computed correctly!
##' @param method character string indicating the method to be used:
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         - "prng":    pure Monte Carlo
##' @param abstol numeric >= 0 providing the absolute precision required.
##'        If abstol = 0, algorithm will run until total number of function
##'        evaluations exceeds fun.eval[2].
##' @param CI.factor Monte Carlo confidence interval multiplier. Algorithm runs
##'        CI.factor * (estimated standard error) < abstol. If CI.factor = 3.3
##'        (default), one can expect the actual absolute error to be less than
##'        abstol in 99.9% of the cases
##' @param fun.eval 2-vector giving the initial function evaluations (in the
##'        first loop; typically powers of 2) and the maximal number of
##'        function evaluations
##' @param max.iter.rqmc maximum number of iterations in the RQMC approach        
##' @param B number of randomizations to get error estimates.
##' @param log logical indicating whether the logarithmic density is to be computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return n-vector with computed density values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert
dnvmix <- function(x, qmix, loc = rep(0, d), scale = diag(d),
                   factor = NULL, # needs to be lower triangular!
                   control = list(), log = FALSE, verbose = TRUE,...)
{
  ## Checks 
  if(!is.matrix(x)) x <- rbind(x)
  d <- ncol(x) # dimension
  if(!is.matrix(scale)) scale <- as.matrix(scale)
  stopifnot(length(loc) == d, dim(scale) == c(d, d))
  
  
  ## Deal with algorithm parameters, see also get.set.parameters():
  ## get.set.parameters() also does argument checking, so not needed here. 
  control <- get.set.parameters(control)
  
  ## Grab method, increment and B 
  method    <- control$method
  increment <- control$increment
  B         <- control$B
  
  ## If factor is not provided, determine it here as a *lower* triangular matrix
  if(is.null(factor)) factor <- t(chol(scale)) # lower triangular
  
  ## 1 Define the quantile function of the mixing variable ###################
  
  ## If 'mix' is "constant" or "inverse.gamma", we use the analytical formulas
  is.const.mix <- FALSE # logical indicating whether we have a multivariate normal
  inv.gam <- FALSE # logical indicating whether we have a multivariate t
  qW <- if(is.character(qmix)) { # 'qmix' is a character vector
    qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
    switch(qmix,
           "constant" = {
             is.const.mix <- TRUE
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
               inv.gam <- TRUE
               df2 <- df / 2
               mean.sqrt.mix <- sqrt(df) * gamma(df2) / (sqrt(2) * gamma((df+1)/2)) # used for preconditioning
               function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
             } else {
               is.const.mix <- TRUE
               mean.sqrt.mix <- 1 # used for preconditioning
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
  
  ## Build result object (log-density)
  n <- nrow(x)
  lres <- rep(-Inf, n) # n-vector of results
  notNA <- rowSums(is.na(x)) == 0
  lres[!notNA] <- NA
  x <- x[notNA,, drop = FALSE] # non-missing data (rows)
  
  ## 2 Actual computation ######################################################
  
  ## Recall that 'scale' is *lower triangular*. For short, let 'scale' = L
  ## Solve L * z = x_i - mu for z, so z = L^{-1} * (x_i - mu)   (d vector)
  ## => z^2 (=> componentwise) = z^T z = (x_i - mu)^T * (L^{-1})^T L^{-1} (x_i - mu)
  ##                           = z^T z = (x_i - mu)^T * (L L^T )^{-1} (x_i - mu)
  ##                           = (x_i - mu)^T * scale^{-1} * (x_i - mu)
  ##                           = quadratic form
  ## Now do this for *all* x_i simultaneously using that L is lower triangular:
  ## Forwardsolve: "right-hand sides" of equation must be in the *columns*, thus t(x)
  z <- forwardsolve(factor, t(x) - loc, transpose = FALSE)
  maha2 <- colSums(z^2) # = sum(z^T z); n-vector of squared Mahalanobis distances from x to mu w.r.t. scale
  ## Note: could probably be done with mahalanobis() but unclear how we would
  ##       get det(scale) then.
  ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
  ##                       = log(prod(diag(R))) = sum(log(diag(R)))
  lrdet <- sum(log(diag(factor)))
  if(!is.finite(lrdet)) stop(paste("Density not defined for singular 'scale' "))
  
  ## Counter
  numiter <- 0 # initialize counter; this will count the number of iterations in the while loop
  
  ## Deal with the different distributions
  if(inv.gam) { # multivariate t
    df.d.2 <- (df + d) / 2
    lres[notNA] <- lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) -
      lrdet - df.d.2 * log1p(maha2 / df)
    error <- 0
  } else if(is.const.mix) { # multivariate normal
    lres[notNA] <- -(d/2) * log(2 * pi) - lrdet - maha2/2
    error <- 0
  } else { 
    ## General case of a multivariate normal variance mixture (RQMC)
    ## Prepare inputs for dnvmix.internal 
    if(method == "sobol") {
      if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
      seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
    }
    ## Initial point-set as vector 
    U0 <- switch(method,
                 "sobol"   = {
                   as.vector(sapply(1:B, function(i) 
                     sobol(control$fun.eval[1], d = 1, randomize = TRUE)))
                 },
                 "gHalton" = {
                   as.vector(sapply(1:B, function(i) 
                     ghalton(control$fun.eval[1], d = 1, method = "generalized")))
                 },
                 "prng"    = {
                   runif(control$fun.eval[1]*B)
                 })
    
    ## Sort maha-distance and divide by 2; store ordering to recover "original
    ## ordering" later:
    ordering.maha <- order(maha2)
    maha2.2 <- maha2[ordering.maha]/2
    
    ## Call internal dnvix (which itself calls C-Code)
    ests <- dnvmix.int(qW, maha2.2 = maha2.2, lrdet = lrdet, U0 = U0, d = d, 
                       method = method, increment = control$increment, 
                       abstol = control$dnvmix.abstol, 
                       reltol = control$dnvmix.reltol,
                       CI.factor = control$CI.factor, fun.eval = control$fun.eval, 
                       max.iter.rqmc = control$max.iter.rqmc, B = B, seed = seed)
    error <- ests$error
    numiter <- ests$numiter 
    ## Finalize
    max.error <- max(error)
    
    ## Extrapolation?
    if(max.error > control$dnvmix.abstol && control$dnvmix.extrap){
      extrap.obj <- dnvmix.extrapolate(maha2.2, ldensities = ests$ldensities,
                                       errors = error, control = control)
      ## Sort ldensities according to original input ordering
      lres[notNA] <- extrap.obj$ldensities[order(ordering.maha)]
      error <- if(log){
        extrap.obj$errors[order(ordering.maha)]
      } else {
        extrap.obj$errors[order(ordering.maha)]*max(exp(max(lres[notNA])), 1)
      }
      if(verbose) warning("'abstol' not reached; extrapolation used")
    } else {
      ## No extrapolation, grab results, sort them, print warning if needed
      lres[notNA] <- ests$ldensities[order(ordering.maha)]
      error <- if(log){
        error[order(ordering.maha)]
      } else {
        error[order(ordering.maha)]*max(exp(max(lres[notNA])), 1)
      }
      if(verbose && max.error > control$dnvmix.abstol) warning("'abstol' not reached; consider increasing 'maxiter.rqmc'")
    }
  }
  ## Return
  attr(lres, "error")   <- error
  attr(lres, "numiter") <- numiter
  if(log) lres else exp(lres)
}
