### dnvmix() ###################################################################


##' @title Extrapolate log-densitiy of a NVM dist'n 
##' @param maha2.2 vector of squared mahalanobis distances (sorted)
##' @param ldensities log-densities of NVM density evaluated at maha2.2
##' @param errors estimated errors of ldensities (from dnvmix.int)
##' @param control list of control arguments, see ?dnvmix
##' @author Erik Hintz

dnvmix.extrapolate <- function(maha2.2, ldensities, errors, tol, control){
  n <- length(maha2.2)
  ## logical vector if precision reached
  notRchd <- (errors >= tol)
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



##' @title Density of a Multivariate Normal Variance Mixture - Internal version 
##'        (not exported)
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lrdet log(sqrt(det(scale))) where 'scale' is the scale matrix of 
##'        the normal variance mixture distribution. 
##' @param d dimension of the Normal Variance Mixture        
##' @param control see ?dnvmix
##' @param verbose see ?dnvmix
##' @return List of three:
##'         $ldensities n-vector with computed log-density values 
##'         $error n-vector of error estimates for log-densities; either relative
##'         error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $numiter number of iterations needed 
##' @author Erik Hintz and Marius Hofert
dnvmix.int <- function(qW, maha2.2, lrdet, d, control, verbose)
{
  ## 1 Basics ##################################################################
  ## Define various quantites:
  ## Note: Most checking was done in dnvmix()
  dblng           <- (control$increment == "doubling")
  B               <- control$B # number of randomizations
  n               <- length(maha2.2) # sample size 
  current.n       <- control$fun.eval[1] #initial sample size 
  numiter         <- 0 # counter for the number of iterations 
  total.fun.evals <- 0
  ## Store seed if 'sobol' is used to get the same shifts later:
  if(control$method == "sobol") {
    if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
    seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
  }
  ## Additional variables needed if the increment chosen is "doubling"
  if(dblng) {
    if(control$method == "sobol") useskip <- 0
    denom <- 1
  }
  ## Matrix to store RQMC estimates
  rqmc.estimates <- matrix(0, ncol = n, nrow = B) 
  ## Absolute or relative precision required?
  CI.factor.sqrt.B <- control$CI.factor / sqrt(B) 
  if(is.na(control$dnvmix.reltol)){
    ## Use absolute error
    tol <- control$dnvmix.abstol/CI.factor.sqrt.B
    do.reltol <- FALSE
  } else {
    ## Use relative error
    tol <- control$dnvmix.reltol/CI.factor.sqrt.B
    do.reltol <- TRUE
  }
  ## Initialize 'max.error' to > tol so that we can enter the while loop:
  max.error <- tol + 42 
  
  ## 2 Main loop ###############################################################
  
  ## while() runs until precision abstol is reached or the number of function
  ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
  ## the desired log-densities are calculated.
  while(max.error > tol && numiter < control$max.iter.rqmc && 
        total.fun.evals < control$fun.eval[2]) {
    
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
                  })
      
      ## 2.2 Evaluate the integrand at the (next) point set #############
      
      W <- qW(U) # realizations of the mixing variable
      c <- - (d/2) * log(2 * pi) - d/2 * log(W) - lrdet - outer(1/W, maha2.2)
      cmax <- apply(c, 2, max)
      next.estimate <- -log(current.n) + cmax + log(colSums(exp(c - rep(cmax, each = current.n))))
      
      # Note: Problem in C function...
      # next.estimate <- .Call("eval_dnvmix_integrand", 
      #                        W          = as.double(sort(W)),
      #                        maha2_2    = as.double(maha2.2),
      #                        current_n  = as.integer(current.n),
      #                        n          = as.integer(n),
      #                        d          = as.integer(d),
      #                        k          = as.integer(d), 
      #                        lrdet      = as.double(lrdet))
      
      ## 2.3 Update RQMC estimates #######################################
      
      rqmc.estimates[b,] <-
        if(dblng) {
          ## In this case both, rqmc.estimates[b,] and
          ## next.estimate depend on n.current points
          (rqmc.estimates[b,] + next.estimate) / denom
        } else {
          ## In this case, rqmc.estimates[b,] depends on
          ## numiter * n.current points whereas next.estimate
          ## depends on n.current points
          (numiter * rqmc.estimates[b,] + next.estimate) / (numiter + 1)
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
    ## Update error. Note that 'tol' was divided by 'CI.factor.sqrt.B'
    errors <- if(!do.reltol){
      apply(rqmc.estimates, 2, sd)
    } else {
      apply(rqmc.estimates, 2, sd)/abs(.colMeans(rqmc.estimates, B, n, 0))
    }
    max.error <- max(errors)
  } # while()
  
  ## Finalize 
  ldensities <- .colMeans(rqmc.estimates, B, n, 0)
  ## Get correct error estimate:
  errors <- errors*CI.factor.sqrt.B
  ## Tolerance not reached?
  if(max.error > tol){
    if(control$dnvmix.extrap){
      ## Extrapolation:
      extrap.obj  <- dnvmix.extrapolate(maha2.2, ldensities = ldensities,
                                        errors = errors, tol = 1.05*tol, 
                                        control = control)
      ldensities  <- extrap.obj$ldensities
      errors      <- extrap.obj$errors
      ## Print waring that extrapolation was used:
      if(verbose){
        if(do.reltol){
          warning("'reltol' not reached for all inputs; extrapolation used")
        } else {
          warning("'abstol' not reached for all inputs; extrapolation used")
        }
      }
    } else if(verbose){
      ## Tolerance not reached, no extrapolation but warning:
      if(do.reltol){
        warning("'reltol' not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
      } else {
        warning("'abstol' not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
      }
    }
  }
  ## 4 Return ##################################################################
  return(list(ldensities = ldensities, numiter = numiter, error = errors))
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
##'         - "PRNG":    pure Monte Carlo
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
  control <- nvmix:::get.set.parameters(control)
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
  numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
  ## Deal with the different distributions
  if(inv.gam) { # multivariate t
    df.d.2 <- (df + d) / 2
    lres[notNA] <- lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) -
      lrdet - df.d.2 * log1p(maha2 / df)
    if(!log) lres <- exp(lres) # already exponentiate
    error <- rep(0, length(z))
  } else if(is.const.mix) { # multivariate normal
    lres[notNA] <- -(d/2) * log(2 * pi) - lrdet - maha2/2
    if(!log) lres <- exp(lres) # already exponentiate
    error <- rep(0, length(z))
  } else { 
    ## General case of a multivariate normal variance mixture (RQMC)
    ## Prepare inputs for dnvmix.int
    ## Sort maha-distance and divide by 2; store ordering to recover original
    ## ordering later:
    ordering.maha <- order(maha2)
    maha2.2 <- maha2[ordering.maha]/2
    ## Call internal dnvix (which itself calls C-Code)
    ## Note that 'dnvmix.int' calls 'dnvmix.extrapolate' (if necessary)
    ests <- dnvmix.int(qW, maha2.2 = maha2.2, lrdet = lrdet, d = d,
                       control = control, verbose = verbose)
    ## Grab results, correct 'error' and 'lres' if 'log = FALSE'
    lres[notNA] <- ests$ldensities[order(ordering.maha)]
    error <- if(log){
      ests$error[order(ordering.maha)]
    } else {
      lres <- exp(lres)
      ests$error[order(ordering.maha)]*pmax(lres[notNA], 1)
    }
    numiter <- ests$numiter 
  }
  ## Return
  ## Note that 'lres' was exponentiated already if necessary. 
  attr(lres, "error")   <- error
  attr(lres, "numiter") <- numiter
  lres 
}

