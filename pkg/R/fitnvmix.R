##### fitnvmix() #################################################################

## 1. Functions to estimate weights ##############################################

#' Weights for fitnvmix() (corresponds to delta_ki in the paper)
#' for the special case of a t dist'n
#'
#' @param tx t(x) with x as in fitnvmix()
#' @param nu degree of freedom parameter
#' @param loc location vector
#' @param scale scale matrix
#' @return vector of length nrow(tx) with corresponding weights
#' @author Erik Hintz
get.weights.t <- function(tx, nu, loc, scale, ...){
   d <- nrow(tx)
   factor <- t(chol(scale))
   z <- forwardsolve(factor, tx - loc, transpose = FALSE) # use the full sample!
   maha2 <- colSums(z^2)
   (nu + d) / (nu + maha2)
}
#' Weights for fitnvmix() (corresponds to delta_ki in the paper)
#' for the special case of a t dist'n
#'
#' @param maha2.2 
#' @param nu degree of freedom parameter
#' @param d dimension
#' @param ... nothing (will be ignored)
#' @return vector with the same lenght as 'maha2.2' with corresponding weights
#' @author Erik Hintz
get.weights.maha2.2 <- function(maha2.2, nu, d,...){
   (nu + d) / (nu + maha2.2*2)
}






#' Estimate weights for fitnvmix() (corresponds to delta_ki in the paper)
#' (internal function)
#'
#' @param maha2.2 squared maha distances divided by 2
#' @param qW see ?fitnvmix() ('qmix' there)
#' @param nu parameter (vector) nu of W
#' @param lrdet log(sqrt(det(scale)))
#' @param d dimension 
#' @param U initial point-set of uniforms to estimate weights
#' @param control see ?fitnvmix()
#' @param seed seed to get the same shifts if 'method = "sobol"'
#' @param verbose see ?fitnvmix
#' @return vector of length nrow(tx) with corresponding weights
#' @author Erik Hintz
get.weights <- function(maha2.2, qW, nu, lrdet, d, U, control, seed, verbose)
{
   ## Define various quantities needed later:
   B                 <- control$B
   CI.factor.sqrt.B  <- control$CI.factor/sqrt(B)
   reltol            <- control$weights.reltol
   abstol            <- control$weights.abstol
   n                 <- length(maha2.2)
   method            <- control$method
   ## Initialize RQMC procedure to estimate the weights
   ## Matrix to store RQMC estimates for weights
   rqmc.estimates.logweights  <- matrix(0, ncol = n, nrow = B)
   total.fun.evals            <- 0
   current.n                  <- length(U)/B
   ## First pointset that was passed as vector: U has length B * current.n
   ## Realizations of l'th shift are elements (l-1)*current.n + (1:current.n)
   W <- qW(U, nu = nu)
   for(b in 1:B){
      ## Grab realizations corresponding to l'th shift and use exp-log trick
      ## Note: The C function eval_dnvmix_integrand performs the following
      ##  b <- - (d/2) * log(2 * pi) - k/2 * log(W) - lrdet - outer(1/W, maha2.2)
      ##  bmax <- apply(b, 2, max)
      ##  log(current_n) + bmax + log(colSums(exp(b - rep(bmax, each = current_n))))
      W.current.sorted <- sort(W[(b-1)*current.n + (1:current.n)])
      condexp <- .Call("eval_dnvmix_integrand",
                       W          = as.double(W.current.sorted),
                       maha2_2    = as.double(maha2.2),
                       current_n  = as.integer(current.n),
                       n          = as.integer(n),
                       d          = as.integer(d),
                       k          = as.integer(d + 2), # note k = d+2 here!
                       lrdet      = as.double(lrdet))
      ldens <- .Call("eval_dnvmix_integrand",
                     W          = as.double(W.current.sorted),
                     maha2_2    = as.double(maha2.2),
                     current_n  = as.integer(current.n),
                     n          = as.integer(n),
                     d          = as.integer(d),
                     k          = as.integer(d),
                     lrdet      = as.double(lrdet))
      
      rqmc.estimates.logweights[b, ] <- condexp - ldens
   }
   ## Get error estimates
   lweights       <- nvmix:::logsumexp(rqmc.estimates.logweights) - log(B) # performs better than .colMeans
   vars           <- .colMeans((rqmc.estimates.logweights - rep(lweights, each = B))^2, B, n, 0)
   errors         <- sqrt(vars)/abs(lweights)*CI.factor.sqrt.B
   max.rel.error  <- max(errors)
   
   ## Update counters 
   total.fun.evals <- B * current.n
   numiter <- 1
   ## Main loop
   while(max.rel.error > reltol && total.fun.evals < control$fun.eval[2] &&
         numiter <= control$max.iter.rqmc)
   {
      if(method == "sobol") .Random.seed <<- seed # reset seed to have the same shifts in sobol( ... )
      ## Get next point-set
      U.next <- switch(method,
                       "sobol"   = {
                          as.vector(sapply(1:B, function(i)
                             sobol(current.n, d = 1, randomize = TRUE, skip = current.n)))
                       },
                       "gHalton" = {
                          as.vector(sapply(1:B, function(i)
                             ghalton(current.n, d = 1, method = "generalized")))
                       },
                       "prng"    = {
                          runif(current.n*B)
                       })
      ## Realizations of W
      W <- qW(U.next, nu = nu)
      for(b in 1:B){
         ## Get realizations of W
         W.current.sorted <- sort(W[(b-1)*current.n + (1:current.n)])
         ## See above for details on eval_dnvmix_integrand
         
         condexp <- .Call("eval_dnvmix_integrand",
                          W          = as.double(W.current.sorted),
                          maha2_2    = as.double(maha2.2),
                          current_n  = as.integer(current.n),
                          n          = as.integer(n),
                          d          = as.integer(d),
                          k          = as.integer(d+2),
                          ##  Note k = d+2 here!
                          lrdet      = as.double(lrdet))
         ldens <- .Call("eval_dnvmix_integrand",
                        W          = as.double(W.current.sorted),
                        maha2_2    = as.double(maha2.2),
                        current_n  = as.integer(current.n),
                        n          = as.integer(n),
                        d          = as.integer(d),
                        k          = as.integer(d),
                        lrdet      = as.double(lrdet))
         
         rqmc.estimates.logweights[b, ] <-
            .Call("logsumexp2",
                  a = as.double(rqmc.estimates.logweights[b, ]),
                  b = as.double(condexp - ldens),
                  n = as.integer(n)) - log(2)
         ## was  '(rqmc.estimates.logweights[b, ] + condexp - ldens)/2' 
      }
      ## Updates:
      total.fun.evals  <- total.fun.evals + B * current.n
      numiter          <- numiter + 1
      current.n        <- 2 * current.n
      ## As in dnvmix():
      ## Update error. The following is slightly faster than 'apply(..., 2, var)' 
      lweights       <- nvmix:::logsumexp(rqmc.estimates.logweights) - log(B) # performs better than .colMeans
      vars           <- .colMeans((rqmc.estimates.logweights - rep(lweights, each = B))^2, B, n, 0)
      errors         <- sqrt(vars)/abs(lweights)*CI.factor.sqrt.B
      max.rel.error  <- max(errors)
   }
   ## Return
   weights <- exp(lweights)
   if(max.rel.error > reltol && verbose) warning('weights.abstol in get.weights not reached')
   weights
}

#' Estimate weights for fitnvmix() (corresponds to delta_ki in the paper)
#' (internal function)
#'
#' @param maha2.2 squared maha distances divided by 2 (length n)
#' @param qW see ?fitnvmix() ('qmix' there)
#' @param nu parameter (vector) nu of W
#' @param lrdet log(sqrt(det(scale)))
#' @param d dimension 
#' @param control see ?fitnvmix()
#' @param verbose see ?fitnvmix
#' @param return.all logical; if true, matrix (U, qW(U)) also returned.
#' @return List of three:
#'         $weights n-vector with computed log-density values
#'         $numiter numeric, number of iterations needed
#'         $error n-vector of error estimates for log-densities; either relative
#'         error or absolte error depending on is.na(control$dnvmix.reltol)
#'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms (only if return.all = TRUE)
#' @author Erik Hintz
weights.internal.RQMC <- function(maha2.2, qW, nu, lrdet, d, max.iter.rqmc,
                                  control, verbose, return.all)
{
   ## 1 Basics ##################################################################
   ## Define various quantites:
   dblng           <- TRUE
   B               <- control$B # number of randomizations
   n               <- length(maha2.2) # sample size
   current.n       <- control$fun.eval[1] #initial sample size
   numiter         <- 0 # counter for the number of iterations
   ZERO            <- .Machine$double.neg.eps
   total.fun.evals <- 0
   ## Absolte/relative precision?
   if(is.na(control$weights.reltol)){
      ## Use absolute error
      tol <- control$weights.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error
      tol <- control$weights.reltol
      do.reltol <- TRUE
   }
   ## Store seed if 'sobol' is used to get the same shifts later:
   if(control$method == "sobol") {
      useskip <- 0 # to get correct shifts if method = "sobol"
      if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
      seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
   }
   denom <- 1 
   ## Matrix to store RQMC estimates
   rqmc.estimates.lweights <- matrix(-Inf, ncol = n, nrow = B)
   ## Will be needed a lot:
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ## Initialize 'max.error' to > tol so that we can enter the while loop:
   max.error <- tol + 42
   ## Matrix to store U, W values => nrows = maximal number of funevals
   if(return.all){
      max.nrow <- current.n*B*2^(max.iter.rqmc-1)
      UsWs <- matrix(NA, ncol = 2, nrow = max.nrow)
      curr.lastrow <- 0 # will count row-index additional points are being inserted after
   }
   ## 2 Main loop ###############################################################
   ## while() runs until precision abstol is reached or the number of function
   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
   ## the desired log-densities are calculated.
   while(max.error > tol && numiter < max.iter.rqmc &&
         total.fun.evals < control$fun.eval[2])
   {
      ## Reset seed to have the same shifts in sobol( ... )
      if(control$method == "sobol" && numiter > 0)
         .Random.seed <<- seed # reset seed to have the same shifts in sobol( ... )
      for(b in 1:B){
         ## 2.1 Get the point set ###########################################
         U <- sort(switch(control$method,
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
                          })) # sorted for later!
         ## 2.2 Evaluate the integrand at the (next) point set #############
         W <- qW(U, nu = nu) # realizations of the mixing variable; sorted!
         ## Need to replace values < ZERO by ZERO. W is *sorted*, so check using
         ## loop instd of 'pmax' (more efficient)
         for(ind in 1:current.n) if(W[ind] < ZERO) W[ind] <- ZERO else break
         ## Update 'UsWs'
         if(return.all){
            UsWs[(curr.lastrow + 1) : (curr.lastrow + current.n), ] <- cbind(U, W)
            curr.lastrow <- curr.lastrow + current.n
         }
         next.est.condexp <- .Call("eval_dnvmix_integrand",
                                   W          = as.double(W),
                                   maha2_2    = as.double(maha2.2),
                                   current_n  = as.integer(current.n),
                                   n          = as.integer(n),
                                   d          = as.integer(d),
                                   k          = as.integer(d + 2), # k=d+2 here!
                                   lrdet      = as.double(lrdet))
         next.est.ldens <- .Call("eval_dnvmix_integrand",
                                 W          = as.double(W),
                                 maha2_2    = as.double(maha2.2),
                                 current_n  = as.integer(current.n),
                                 n          = as.integer(n),
                                 d          = as.integer(d),
                                 k          = as.integer(d), # k=d+2 here!
                                 lrdet      = as.double(lrdet))
         
         ## 2.3 Update RQMC estimates #######################################
         rqmc.estimates.lweights[b, ] <-
            .Call("logsumexp2",
                  a = as.double(rqmc.estimates.lweights[b, ]),
                  b = as.double(next.est.condexp - next.est.ldens),
                  n = as.integer(n)) - log(denom)
         
      } # end for(b in 1:B)
      ## Update of various variables
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(numiter == 0){
         ## Change denom and useksip (exactly once, in the first iteration)
         denom <- 2
         useskip <- 1
      } else {
         ## Increase sample size n. This is done in all iterations
         ## except for the first two
         current.n <- 2 * current.n
      }
      ## Total number of function evaluations:
      total.fun.evals <- total.fun.evals + B * current.n
      numiter <- numiter + 1
      ## Update error for 'weights' (NOT log values!)
      rqmc.estimates.weights <- exp(rqmc.estimates.lweights)
      weights  <- colMeans(rqmc.estimates.weights) 
      vars     <- .colMeans((rqmc.estimates.weights - rep(weights, each = B))^2, B, n, 0)
      errors   <- if(!do.reltol){
         sqrt(vars)*CI.factor.sqrt.B
      } else {
         sqrt(vars)/abs(weights)*CI.factor.sqrt.B
      }
      max.error <- max(errors)
   } # while()
   ## 4 Return ##################################################################
   ret.obj <- if(return.all){
      list(weights = weights, numiter = numiter, error = errors,
           UsWs = UsWs)
   } else {
      list(weights = weights, numiter = numiter, error = errors)
   }
   ret.obj
}

#' Estimate weights for fitnvmix() (corresponds to delta_ki in the paper)
#' (internal function)
#'
#' @param maha2.2 squared maha distances divided by 2 (length n)
#' @param qW see ?fitnvmix() ('qmix' there)
#' @param nu parameter (vector) nu of W
#' @param lrdet log(sqrt(det(scale)))
#' @param d dimension 
#' @param control see ?fitnvmix()
#' @param verbose see ?fitnvmix
#' @return List of three:
#'         $weights n-vector with computed log-density values
#'         $numiter numeric, number of iterations needed
#'         $error n-vector of error estimates for log-densities; either relative
#'         error or absolte error depending on is.na(control$dnvmix.reltol)
#'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms (only if return.all = TRUE)
#' @author Erik Hintz
weights.internal <- function(maha2.2, qW, nu, lrdet, d, control,
                             verbose)
{
   ## Absolte/relative precision?
   if(is.na(control$weights.reltol)){
      tol <- control$weights.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error
      tol <- control$weights.reltol
      do.reltol <- TRUE
   }
   
   ## Call RQMC procedure without any stratification
   rqmc.obj <- nvmix:::weights.internal.RQMC(maha2.2, qW = qW, nu = nu, lrdet = lrdet, 
                                             d = d, max.iter.rqmc = control$dnvmix.max.iter.rqmc.pilot,
                                             control = control, verbose = TRUE, 
                                             return.all = TRUE)
   ## Extract results
   weights <- rqmc.obj$weights
   numiter <- rep(rqmc.obj$numiter, length(maha2.2))
   error   <- rqmc.obj$error
   if(any(error > tol)){
      ## Accuracy not reached for at least one 'maha2.2' value
      ## => Use adaptive approach for those
      notRchd <- which(error > tol)
      qW. <- function(u) qW(u, nu = nu)
      ldens.obj <- nvmix:::dnvmix.internal.adaptRQMC(qW., maha2.2 = maha2.2[notRchd], lrdet = lrdet, 
                                                     d = d, UsWs = rqmc.obj$UsWs, 
                                                     control = control)
      lcond.obj <- nvmix:::dnvmix.internal.adaptRQMC(qW., maha2.2 = maha2.2[notRchd], lrdet = lrdet, 
                                                     d = d, k = d + 2, UsWs = rqmc.obj$UsWs, 
                                                     control = control)
      weights[notRchd]   <- exp(lcond.obj$ldensities - ldens.obj$ldensities)
      ## Which weights cannot be reliably estimated?
      which.errorNA     <- which(is.na(lcond.obj$error) | is.na(ldens.obj$error))
      if(any(which.errorNA)){
         if(verbose) warning('Some weights cannot be reliably estimated')
         ## Check if 'weights' decreasing in 'maha2.2'
         n <- length(weights)
         for(i in 1:n){
            if(i <= 2) next else if(weights[i] <= weights[i-1]) next else {
               ## In this case, weights[i] > weights[i-1]
               ## Case 1: Is there any weight beyond i which is smaller than weights[i-1]?
               smallerthani <- which(weights[i:n] <= weights[i-1]) + i - 1
               if(length(smallerthani) > 0){
                  ## If that's the case, interpolate all weights between i 
                  ## and the first one smaller than weights[i-1]
                  firstsmaller <- smallerthani[1]
                  slope <- (weights[firstsmaller] - weights[i-1])/(maha2.2[firstsmaller] - maha2.2[i-1])
                  weights[i:(firstsmaller-1)] <- weights[i-1] + slope*
                     (maha2.2[i:(firstsmaller-1)] - maha2.2[i-1])
               } else {
                  ## In this case we extrapolate
                  ## Following extrapolates linearly based on weights, maha2.2
                  ## slope <- (weights[i-1] - weights[i-2])/(maha2.2[i-1] - maha2.2[i-2])
                  ## weights[i:n] <- pmax(weights[i-1] + slope*
                  ##   (maha2.2[i:n] - maha2.2[i-1]), 1e-16)
                  ##   
                  ## Behavior of the function suggests log-log extrapolation:
                  slope <- log(weights[i-1]/weights[i-2])/log(maha2.2[i-1]/maha2.2[i-2])
                  weights[i:n] <- weights[i-1]*exp(slope*log(maha2.2[i:n]/maha2.2[i-1]))
               }
            }
         }
      }
      #weights.est[errorNA] <- get.weights.maha2.2(maha2.2[notRchd][errorNA], nu = nu, d = d)
      #weights[notRchd]  <- weights.est 
      # numiter[notRchd]  <- numiter[notRchd] + rqmc.obj$numiter
      # error[notRchd]    <- rqmc.obj$error
      # ## Handle warnings:
      # if(verbose){
      #    if(any(is.na(error))){
      #       ## At least one error is NA
      #       warning("Estimation unreliable, corresponding error estimate NA")
      #    }
      #    whichNA <- which(is.na(error))
      #    if(any(error[setdiff(1:length(error), whichNA)] > tol)) # 'setdiff' needed if 'whichNA' is empty
      #       warning("Tolerance not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
      # }
      # ## Transform error back to *absolute* errors:
      # if(do.reltol) error <- error * abs(ldens)
   }
   list(weights = weights, numiter = numiter, error = error)
}

## 2. Functions to estimate wnu given 'loc', 'scale'############################

# Estimate nu using copula density (internal function)#'
#' @param U matrix of pseudo obs
#' @param qW quantile function of W; must be function(u, nu)
#' @param init.nu initial estimate of nu
#' @param factor cholesky factor of the scale matrix
#' @param control see ?fitnvmix()
#' @param control.optim passed to optim; see ?optim
#' @param mix.param.bounds see ?fitnvmix()
#' @param verbose see ?fitnvmix()
#' @param inv.gam logical indicating if W is inv.gamma (special case)
#' @param seed current seed.
#' @return list of two: $nu.est (scalar of vector of length init.nu; MLE estimate of nu)
#'                      $max.ll (negative log-likelihood at nu.est)
#' @author Erik Hintz
estim.nu.cop <- function(U, qW, init.nu, factor, control, control.optim,
                         mix.param.bounds, verbose, inv.gam, seed){
   ## Set up -log likelihood as a function of nu (given P), based on copula
   neg.log.likelihood.nu <- if(inv.gam){
      function(nu){
         -sum(dnvmixcop(U, qmix = "inverse.gamma", factor = factor,
                        control = control, verbose = verbose,
                        log = TRUE, df = nu))}
   } else {
      function(nu){
         .Random.seed <<- seed
         qmix. <- function(u) qW(u, nu = nu) # function of u
         - sum(dnvmixcop(U, qmix = qmix., factor = factor, control = control,
                         verbose = verbose, log = TRUE))}
   }
   ## Optimize neg.log.likelihood over nu
   opt.obj <- optim(init.nu, fn = neg.log.likelihood.nu,
                    lower = mix.param.bounds[, 1],
                    upper = mix.param.bounds[, 2],
                    method = "L-BFGS-B", control = control.optim)
   list(nu.est = opt.obj$par,
        max.ll = opt.obj$value)
}




#' Estimate nu given loc, scale by maximizing log-likelihood (internal function)
#'
#' @param tx t(x) where x is as in ?fitnmvix()
#' @param qW quantile function of W; must be function(u, nu)
#' @param init.nu initial estimate of nu
#' @param factor cholesky factor of the scale matrix
#' @param control see ?fitnvmix()
#' @param control.optim passed to optim; see ?optim
#' @param mix.param.bounds see ?fitnvmix()
#' @param inv.gam logical indicating if W is inv.gamma (special case)
#' @param U0 vector of uniforms, see also nvmix:::dnvmix.int
#' @param seed seed used to produce U0
#' @param verbose see ?fitnvmix
#' @return list of two: $nu.est (scalar of vector of length init.nu; MLE estimate of nu)
#'                      $max.ll (negative log-likelihood at nu.est)
#' @author Erik Hintz
estim.nu <- function(tx, qW, init.nu, loc, scale, control, control.optim,
                     mix.param.bounds, inv.gam = FALSE, seed, verbose)
{
   factor <- t(chol(scale))
   if(inv.gam){ ## in this case, dnvmix() uses analytical formula for density
      neg.log.likelihood.nu <- function(nu){
         -sum(dnvmix(t(tx), qmix = "inverse.gamma", loc = loc, scale = scale,
                     df = nu, log = TRUE, verbose = verbose))
      }
   } else {
      ## Get various quantitites passed to 'dnvmix.internal'
      z <- forwardsolve(factor, tx - loc, transpose = FALSE)
      maha2.2 <- sort(colSums(z^2)/2)
      lrdet <- sum(log(diag(factor)))
      d <- ncol(factor)
      ## Set up -loglikelihood as a function of 'nu'
      neg.log.likelihood.nu <- function(nu){
         .Random.seed <<- seed # reset seed => monotonicity (not bc of sobol shifts!)
         qmix. <- function(u) qW(u, nu = nu) # function of u only
         ## Call nvmix:::dnvmix.int which by default returns the log-density
         ldens.obj <- nvmix:::dnvmix.internal(qW = qmix., maha2.2 = maha2.2, lrdet = lrdet,
                                              d = d, control = control, verbose = verbose)
         ## Return -log-density
         -sum(ldens.obj$ldensities)
      }
   }
   ## Optimize neg.log.likelihood over nu
   opt.obj <- optim(init.nu, fn = neg.log.likelihood.nu,
                    lower = mix.param.bounds[, 1],
                    upper = mix.param.bounds[, 2],
                    method = "L-BFGS-B", control = control.optim)
   list(nu.est = opt.obj$par,
        max.ll = opt.obj$value,
        num.llevals = opt.obj$counts[1])
}


## 3. Main function ('fitnvmix') that is exported ###############################

##' @title Fitting Multivariate Normal Variance Mixtures
##'        (only Student-t right now)
##' @param x (n,d) data matrix
##' @param qmix character string ("constant", "inverse.gamma") or function. If
##'        function, it *has* to be qmix(u, nu) and the length of nu has
##'        to be provided in mix.param.length.
##' @param mix.param.length
##' @param mix.param.bounds either NA or (mix.param.length, 2) matrix, where
##'         1st/2nd column corresponds to lower/upper limits for i'th element
##'         of nu
##'         (eg if W~exp(nu), then mix.param.length = 1 and mix.param.bounds =
##'         c(0, Inf))
##'         Can have NAs for unrestricted parameters
##'         TODO: Re-think if to use NA or Inf etc
##'         TODO: Replace '0' by 'zero' internally?
##' @param size.subsample
##' @param resample
##' @param ECMEstep
##' @param control.optim
##' @param control
##' @param verbose
##' @return list of three (if qmix = "constant"), otherwise five:
##'         $nu: estimate for nu (omitted if qmix = "constant")
##'         $loc: estimate for the location vector
##'         $scale: estimate for scale matrix
##'         $iter.ECME: of EM iterations
##'         $max.ll: log-likelihood at (nu, loc, scale)
##' TODO    include option to give names to parameters etc
##' TODO    maybe hide some arguments (CI.factor etc) or have them passed via ...
##' @author Erik Hintz
fitnvmix <- function(x, qmix,
                     nu.init = NA, init.size.subsample = min(n, 50),
                     mix.param.length = 1, mix.param.bounds = NA,
                     size.subsample = n, resample = FALSE,
                     ECMEstep = TRUE, control.optim = list(),  control = list(),
                     verbose = TRUE,
                     useCop = FALSE)
{
   ## 0: Initialize various quantities: #######################################
   control <- nvmix:::get.set.parameters(control)
   ## Get quantile function:
   ## If 'mix' is "constant" or "inverse.gamma", we use the analytical formulas
   is.const.mix  <- FALSE # logical indicating whether we have a multivariate normal
   inv.gam       <- FALSE # logical indicating whether we have a multivariate t
   inv.gam.weights <- FALSE
   ## Set up qW as function(u, nu)
   qW <- if(is.character(qmix)) {# 'qmix' is a character vector
      qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
      switch(qmix,
             "constant" = {
                is.const.mix <- TRUE
                function(u) 1
             },
             "inverse.gamma" = {
                inv.gam <- TRUE
                mix.param.length <- 1
                mix.param.bounds <- matrix(c(0.5, NA), ncol = 2)
                function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
             },
             stop("Currently unsupported 'qmix'"))
   } else if(is.list(qmix)) { # 'mix' is a list of the form (<character string>, <parameters>)
      #TODO: Do this!
      #stopifnot(length(qmix) >= 1, is.character(distr <- qmix[[1]]))
      #qmix. <- paste0("q", distr)
      #if(!existsFunction(qmix.))
      #  stop("No function named '", qmix., "'.")
      #function(u)
      # do.call(qmix., append(list(u), qmix[-1]))
   } else if(is.function(qmix)) { # 'mix' is the quantile function F_W^- of F_W
      function(u, nu)
         qmix(u, nu)
   } else stop("'qmix' must be a character string, list or quantile function.")
   ## Case of MVN: MLEs are sample mean and sample cov matrix
   if(is.const.mix){
      loc.est <- colMeans(x)
      ## TODO: Do this better (as.matrix to remove attributes, can be done better)
      scale.est <- as.matrix(nearPD(cov(x))$mat) # sample covariance matrix
      return(list(loc = loc.est, scale = scale.est, iter = 0))
   }
   ## Check inputs, get dimensions
   ## TODO: More checking (INFs, ...)
   if(!is.matrix(x)) x <- rbind(x)
   notNA <- rowSums(is.na(x)) == 0
   x     <- x[notNA,] # non-missing data (rows)
   tx    <- t(x)
   n     <- nrow(x)
   d     <- ncol(x)
   ## Use only sub-sample to estimate nu?
   if(size.subsample < n){
      sampled.ind <- sample(n, size.subsample)
      x.sub       <- x[sampled.ind,]
      tx.sub      <- tx[, sampled.ind]
   } else {
      x.sub  <- x
      tx.sub <- tx
   }
   ## Check/define parameter bounds on nu
   if(!any(!is.na(mix.param.bounds))){
      mix.param.bounds <- cbind(rep(-Inf, mix.param.length), rep(Inf, mix.param.length))
   } else{
      stopifnot(all.equal( dim(mix.param.bounds), c(mix.param.length, 2)))
   }
   ## Get initial pointset that is being reused again and again:
   if(!exists(".Random.seed")) runif(1)
   seed <- .Random.seed
   U0 <- switch(control$method,
                "sobol"   = {
                   as.vector(sapply(1:control$B, function(i)
                      sobol(control$fun.eval[1], d = 1, randomize = TRUE)))
                },
                "gHalton" = {
                   as.vector(sapply(1:control$B, function(i)
                      ghalton(control$fun.eval[1], d = 1, method = "generalized")))
                },
                "prng"    = {
                   runif(control$fun.eval[1]*control$B)
                })
   if(control$addReturns){
      ## Matrix storing all 'nu' estimates with corresponding likelihoods 
      nu.Ests <- matrix(NA, ncol = mix.param.length + 1, nrow = control$ECME.maxiter + 2)
      rownames(nu.Ests) <- c("Initial", 
                             sapply(1:control$ECME.maxiter, function(i) paste0("ECME-iteration ", i)),
                             "Laststep")
      colnames(nu.Ests) <- c(sapply(1:mix.param.length, function(i) paste0("nu[", i, "]")),
                             "Log-likelihood")
      current.iter.total <- 1
   }
   
   ## 1: Initial estimates for nu, loc, scale: ##################################
   ## Unbiased estimator for 'loc' based on full sample:
   loc.est <- colMeans(x)
   ## Sample covariance matrix based on full sample:
   SCov    <- as.matrix(nearPD(cov(x))$mat) # TODO do this smarter
   # ## Determine maha distance and determinant of 'SCov'
   # chol.SCov <- t(chol(SCov))
   # z         <- forwardsolve(chol.SCov, tx.sub - loc.est, transpose = FALSE)
   # maha2.2   <- sort(colSums(z^2))/2 # sorted for nvmix:::dnvmix.int
   # lrdet     <- sum(log(diag(chol.SCov)))
   ## Check if init.nu was provided. If so, calculate 'scale' as (1/E(W))*SCov
   if(!is.na(nu.init)){
      nu.est <- nu.init
      scale.est <- 1/mean(qW(runif(100000), nu.est))* SCov
      scale.est <- SCov
   } else if(!useCop) {
      ## Optimize log-likelihood:
      ## -loglikelihood as function of param=(nu, c) of length mix.param.length + 1
      neg.log.likelihood.init <- function(param){
         if(inv.gam){
            ## In case of inv.gam, a closed formula for the density exists:
            return(-sum(dnvmix(x.sub, qmix = "inverse.gamma", 
                               loc = loc.est, scale = param[2] * SCov, 
                               df = param[1], log = TRUE)))
         } else {
            ## Define a 'qmix' function of u only that can be passed to dnvmix():
            qmix. <- function(u) qW(u, nu = param[1:mix.param.length]) 
            .Random.seed <<- seed # for monotonicity
            ## Return - loglikelihood
            -sum(dnvmix(x.sub, qmix = qmix., loc = loc.est, 
                        scale = param[mix.param.length + 1] * SCov,
                        control = control, verbose = verbose, log = TRUE))
         }
      }
      ## Optimize -log.likelihood over (nu, c)
      init.param <- c(1, 1)
      opt.obj <- optim(init.param, fn = neg.log.likelihood.init,
                       lower = c(mix.param.bounds[, 1], 0.1),
                       upper = c(mix.param.bounds[, 2], NA),
                       method = "L-BFGS-B", control = control.optim)
      ## Grab estimate for 'nu' as well as for the 'scale' matrix
      nu.est    <- opt.obj$par[1:mix.param.length]
      scale.est <- opt.obj$par[mix.param.length + 1] * SCov
      max.ll    <- opt.obj$value
   } else {
      ## Here ('useCop = TRUE'), take a subsample of size 'init.size.subsample',
      ## transfrom to pseudo-Obs, optimize copula log-likelihood over 'nu'
      pObs <-  copula::pobs(x) # pseudo obs: full sample
      ## Correlation matrix based on full pseudo-sample
      corr.Pobs <- as.matrix(nearPD(cor(pObs))$mat) 
      ## Subset of 'pObs' used to determine likelihood
      pObs.sub <- pObs[sample(n, init.size.subsample), ]
      ## Set up -log likelihood as a function of nu (given P), based on copula
      neg.log.cop.likelihood.nu <- if(inv.gam){
         function(nu){
            -sum(dnvmixcop(pObs.sub, qmix = "inverse.gamma", scale = corr.Pobs,
                           control = control, verbose = verbose,
                           log = TRUE, df = nu))}
      } else {
         function(nu){
            .Random.seed <<- seed
            qmix. <- function(u) qW(u, nu = nu) # function of u
            - sum(dnvmixcop(pObs.sub, qmix = qmix., scale = corr.Pobs, control = control,
                            verbose = verbose, log = TRUE))}
      }
      ## Optimize neg.log.likelihood over nu
      opt.obj <- optim(rowMeans(mix.param.bounds), fn = neg.log.cop.likelihood.nu,
                       lower = mix.param.bounds[, 1],
                       upper = mix.param.bounds[, 2],
                       method = "L-BFGS-B", control = control.optim)
      ## Obtain optimal nu
      nu.est <- opt.obj$par
      scale.est <- 1/mean(qW(runif(500), nu.est))* SCov
      list(nu.est = opt.obj$par,
           max.ll = opt.obj$value)
   }
   if(control$addReturns){
      ## Matrix storing all 'nu' estimates with corresponding likelihoods 
      ll <- sum(dnvmix(x, qmix = qW, loc  = loc.est, scale = scale.est, 
                       nu = nu.est, control = control, verbose = verbose, log = TRUE))
      nu.Ests[current.iter.total, ] <- c(nu.est, ll)
      current.iter.total <- current.iter.total + 1
      iter.converged <- control$ECME.maxiter + 2
   }
   
   ## 2: ECME step: #############################################################
   
   if(ECMEstep){
      
      ## Initialize various quantities
      iter.ECME            <- 0
      converged            <- FALSE
      while(iter.ECME < control$ECME.maxiter && !converged){
         
         ## 2.1: 'loc.est' and 'scale.est' updates #############################
         
         converged.locscale   <- FALSE
         iter.locscaleupdate  <- 1
         
         ## Update 'scale.est' and 'loc.est' given current estimate of 'nu.est'
         ## until convergence. 
         
         while(!converged.locscale && 
               iter.locscaleupdate < control$max.iter.locscaleupdate){
            ## Get new 'weights'
            ## First, get new maha distances (with current 'loc.est' and 'scale.est')
            factor               <- t(chol(scale.est))
            lrdet                <- sum(log(diag(factor)))
            z                    <- forwardsolve(factor, tx - loc.est, transpose = FALSE) # use the full sample!
            maha2.2.new          <- colSums(z^2)/2
            order.maha2.2.new    <- order(maha2.2.new)
            maha2.2.new          <-  maha2.2.new[order.maha2.2.new] # sorted increasingly
            if(iter.locscaleupdate == 1){
               ## Only in the first iteration do we approximate all weights by RQMC.
               ## Get weights 
               if(inv.gam.weights){
                  weights <- nvmix:::get.weights.maha2.2(maha2.2.new, nu = nu.est, d = d)
               } else {
                  #weights <- nvmix:::get.weights(maha2.2.new, qW = qW, nu = nu.est,
                  #                               lrdet = lrdet, d = d, U = U0, 
                  #                               control = control, seed = seed, 
                  #                               verbose = verbose)
                  weights <- nvmix:::weights.internal(maha2.2.new, qW = qW, nu = nu.est,
                                                      lrdet = lrdet, d = d, control = control,
                                                      verbose = verbose)$weights
               }
               weights.new <- weights[order(order.maha2.2.new)]
               maha2.2     <- maha2.2.new # need to store maha-distances for interpolation
               length.maha  <- n # store length of 'maha2.2' and 'weights'
            } else {
               ## Linearly interpolate 'weights' to get new weights
               weights.new          <- rep(NA, n)
               curr.index           <- 1 # index to look for close values in 'maha2.2'
               notInterpol          <- rep(NA, n)
               notInterpolcounter   <- 1
               for(ind in 1:n){
                  curr.maha2.2 <- maha2.2.new[ind]
                  if(curr.maha2.2 < maha2.2[1] || curr.maha2.2 > maha2.2[length.maha]){
                     notInterpol[notInterpolcounter] <- ind
                     notInterpolcounter <- notInterpolcounter + 1
                  } else {
                     ## Start looking for close maha values in 'maha2.2'
                     found <- FALSE
                     while(!found && curr.index < length.maha){
                        if(maha2.2[curr.index] <= curr.maha2.2 && 
                           curr.maha2.2 <= maha2.2[curr.index+1]){
                           found <- TRUE
                           ## Now check if we can interpolate (ie rel.error small)
                           if(abs(weights[curr.index+1] - weights[curr.index])/
                              weights[curr.index+1] < control$weights.interpol.reltol){
                              weights.new[ind] <- weights[curr.index] + 
                                 (curr.maha2.2 - maha2.2[curr.index])*
                                 (weights[curr.index+1] - weights[curr.index])/
                                 (maha2.2[curr.index+1]-maha2.2[curr.index])
                           } else {
                              ## if not, will use 'get.weights' for this maha.
                              notInterpol[notInterpolcounter] <- ind
                              notInterpolcounter <- notInterpolcounter + 1
                           }
                        } else {
                           curr.index <- curr.index + 1
                        }
                     }
                  }
               }
               ## Now need to approximate weights for those maha in 'notInterpol'
               if(notInterpolcounter > 1){
                  notInterpol <- notInterpol[1:(notInterpolcounter-1)]
                  weights.new[notInterpol] <- if(inv.gam.weights){
                     nvmix:::get.weights.maha2.2(maha2.2.new[notInterpol], 
                                                 nu = nu.est, d = d)
                  } else {
                     #nvmix:::get.weights(maha2.2.new[notInterpol], qW = qW, 
                     #                    nu = nu.est, lrdet = lrdet, d = d, U = U0, 
                     #                    control = control, seed = seed, 
                     #                    verbose = verbose)
                     nvmix:::weights.internal(maha2.2.new[notInterpol], qW = qW, nu = nu.est,
                                              lrdet = lrdet, d = d, control = control,
                                              verbose = verbose)$weights
                  }
                  ## Add estimated weights to 'weights' and corresponding 
                  ## 'maha2.2.new' to 'maha2.2' so that they can be reused
                  maha2.2 <- c(maha2.2, maha2.2.new[notInterpol])
                  temp.ordering <- order(maha2.2) # only needed here
                  maha2.2 <- maha2.2[temp.ordering] # sort 'maha2.2' again
                  weights <- c(weights, weights.new[notInterpol])[temp.ordering] # and weights accordingly
                  length.maha <- length.maha + notInterpolcounter - 1
               }
               ## Recover original ordering and set negative weights to zero.
               weights.new <- weights.new[order(order.maha2.2.new)] ## TODO: Omit?
               
               #                OLD:
               #               
               #                ## First: We *avoid* EXTRApolation and estimate weights for those maha's 
               #                ## that are out of range of the previous maha's
               #                anyTooSmall <- (maha2.2.new[1] < maha2.2[1])
               #                anyTooBig   <- (maha2.2.new[n] > maha2.2[n])
               #                if(anyTooBig || anyTooSmall){
               #                   ## Find indices of too small and too big 'maha2.2new' values
               #                   if(anyTooSmall){
               #                      lastTooSmall <- which.min((maha2.2.new < maha2.2[1]))[1]-1
               #                      length.maha  <- length.maha + lastTooSmall
               #                      whichTooSmall <- 1:lastTooSmall
               #                   } else {
               #                      whichTooSmall <- c()
               #                   }
               #                   if(anyTooBig){
               #                      firstTooBig <- which.max((maha2.2.new > maha2.2[length.maha]))[1]
               #                      length.maha <- length.maha + n - firstTooBig + 1
               #                      whichTooBig <- firstTooBig:n
               #                   } else {
               #                      whichTooBig <- c()
               #                   }
               #                   ## Update 'maha2.2' so that those values can be used again
               #                   maha2.2 <- c(maha2.2.new[whichTooSmall], maha2.2, maha2.2.new[whichTooBig])
               #                   weights.new[c(whichTooSmall, whichTooBig)] <- 
               #                      if(inv.gam){
               #                         nvmix:::get.weights.maha2.2(maha2.2.new[c(whichTooSmall, whichTooBig)], 
               #                                                     nu = nu.est, d = d)
               #                      } else {
               #                         nvmix:::get.weights(maha2.2.new[c(whichTooSmall, whichTooBig)], qW = qW, 
               #                                             nu = nu.est, lrdet = lrdet, d = d, U = U0, 
               #                                             control = control, seed = seed, 
               #                                             verbose = verbose)
               #                      }
               #                   ## Store those additionally estimated weights in 'weights'
               #                   weights <- c(weights.new[whichTooSmall], weights, weights.new[whichTooBig])
               #                   whichInRange <- setdiff(1:n, c(whichTooSmall, whichTooBig))
               #                } else {
               #                   whichInRange <- 1:n
               #                }
               #                ## Second: Interpolation for those weights whose maha's are in the range
               #                curr.index <- 1 # current index to look up values in 'maha2.2'
               #                for(ind in whichInRange){
               #                   curr.maha2.2 <- maha2.2.new[ind]
               #                   found <- FALSE
               #                   while(!found && curr.index < length.maha){
               #                      if(maha2.2[curr.index] <= curr.maha2.2 && curr.maha2.2 <= maha2.2[curr.index+1]){
               #                         found <- TRUE
               #                         weights.new[ind] <- weights[curr.index] + (curr.maha2.2 - maha2.2[curr.index])*
               #                            (weights[curr.index+1] - weights[curr.index])/
               #                            (maha2.2[curr.index+1]-maha2.2[curr.index])
               #                      } else {
               #                         curr.index <- curr.index + 1
               #                      }
               #                   }
               #                }
            }
            ## Get new 'scale.est': 1/n * sum_{i=1}^n weights_i (x_i-mu)(x_i-mu)^T
            ## where 'mu' corresponds to current 'loc.est'
            scale.est.new <- crossprod(sqrt(weights.new)*sweep(x, 2, loc.est, check.margin = FALSE))/n
            

            
            ## Get new 'loc.est': sum_{i=1}^n weights_i x_i / (sum weights)
            ## as.vector because we need 'loc.est' as a vector, not (d, 1) matrix
            loc.est.new <- as.vector(crossprod(x, weights.new)/sum(weights.new))
            ## Did we converge?
            scale.est.rel.diff <- abs((scale.est - scale.est.new)/scale.est)
            loc.est.rel.diff   <- abs((loc.est - loc.est.new)/loc.est)
            converged.locscale <- (max(loc.est.rel.diff) < control$ECME.rel.conv.tol[2]) &&
               (max(scale.est.rel.diff) < control$ECME.rel.conv.tol[2])
            
            ## Update counter
            iter.locscaleupdate <- iter.locscaleupdate + 1
            
            ## Update 'loc.est' and 'scale.est'
            loc.est     <- loc.est.new
            scale.est   <- scale.est.new
         }
         
         ## 2.2: Update 'nu.est', if desired/necessary: ########################
         if(control$ECMEstep.do.nu){
            ## New subsample used for this 'nu' update?
            if(resample && size.subsample < n){
               ## TODO: Is this sketchy? 
               rm(".Random.seed") # destroy the reseted seed 
               runif(1) # get a new seed
               sampled.ind <- sample(n, size.subsample)
               tx.sub      <- tx[,sampled.ind]
            }
            ## Optimize neg.log.likelihood over nu
            est.obj <- nvmix:::estim.nu(tx.sub, qW = qW, init.nu = nu.est,
                                        loc = loc.est, scale = scale.est,
                                        control = control, control.optim = control.optim,
                                        mix.param.bounds = mix.param.bounds, inv.gam = inv.gam,
                                        seed = seed, verbose = verbose)
            nu.est.new        <- est.obj$nu.est
            nu.est.rel.diff   <- abs((nu.est.new - nu.est)/nu.est)
            nu.est            <- nu.est.new
            max.ll            <- est.obj$max.ll
         } else {
            nu.est.rel.diff <- 0
         }
         ## Did we converge?
         converged <- (abs(nu.est.rel.diff) < control$ECME.rel.conv.tol[3])
         ## Update counter and 'nu.Ests'
         iter.ECME <- iter.ECME + 1
         if(control$addReturns){
            ## Store new 'nu.est' along with log-likelihood
            nu.Ests[current.iter.total, ] <- c(nu.est, -max.ll)
            ## If 'converged', set all future iteration values to current one
            if(converged){
               nu.Ests[current.iter.total:(dim(nu.Ests)[1]-1), ] <- 
                  matrix(c(nu.est, -max.ll), ncol = mix.param.length+1, 
                         nrow = dim(nu.Ests)[1]-current.iter.total,
                         byrow = TRUE)
               iter.converged <- current.iter.total
            }
            current.iter.total <- current.iter.total + 1
         }
      } # end while()
   } #end if(ECMEstep)
   
   ## 3: Another last 'nu.est' with *full* sample? #############################   
   if(control$laststep.do.nu){
      ## One last nu update with the *full* sample.
      est.obj <- nvmix:::estim.nu(tx, qW = qW, init.nu = nu.est, loc = loc.est,
                                  scale = scale.est, control = control, 
                                  control.optim = control.optim,
                                  mix.param.bounds = mix.param.bounds, 
                                  inv.gam = inv.gam,  seed = seed, verbose = verbose)
      nu.est <- est.obj$nu.est
      max.ll <- est.obj$max.ll
      if(control$addReturns){
         ## Store new 'nu.est' along with log-likelihood
         nu.Ests[dim(nu.Ests)[1], ] <- c(nu.est, -max.ll)
         current.iter.total <- current.iter.total + 1
      }
   }
   
   ## 4: Return ################################################################
   if(control$addReturns){
      return(list(nu = nu.est, loc = loc.est, scale = scale.est, iter = iter.ECME,
                  max.ll = -max.ll, nu.Ests = nu.Ests, iter.converged = iter.converged))
   } else {
      return(list(nu = nu.est, loc = loc.est, scale = scale.est, iter = iter.ECME,
                  max.ll = -max.ll))
   }
}