### pnvmix() ###################################################################

##' @title Swap Variables i and j in a, b and R
##' @param i variable to be switched with j
##' @param j variable to be switched with i
##' @param lower d-vector of lower evaluation limits
##' @param upper d-vector of upper evaluation limits
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @return list with lower, upper and scale after components/rows/columns
##'         i and j have been switched
##' @author Erik Hintz and Marius Hofert
swap <- function(i, j, lower, upper, scale)
{
  ## Build vector
  ij <- c(i,j)
  ji <- c(j,i)
  ## Reorder lower, upper
  lower[ij] <- lower[ji]
  upper[ij] <- upper[ji]
  ## Reorder scale
  wo.ij <- setdiff(seq_len(nrow(scale)), ij)
  temp_i <- scale[wo.ij,i, drop = FALSE]
  temp_j <- scale[wo.ij,j, drop = FALSE]
  temp_ii <- scale[i,i]
  scale[wo.ij,i] <- temp_j
  scale[wo.ij,j] <- temp_i
  scale[i,wo.ij] <- temp_j
  scale[j,wo.ij] <- temp_i
  scale[i,i] <- scale[j,j]
  scale[j,j] <- temp_ii
  ## Return
  list(lower = lower, upper = upper, scale = scale)
}

##' @title Preconditioning (Reordering Variables According to their Expected Integration Limits)
##' @param lower d-vector of lower evaluation limits
##' @param upper d-vector of upper evaluation limits
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param cholScale Cholesky factor (lower triangular matrix) of 'scale'
##' @param mean.sqrt.mix E(sqrt(W)) or NULL; the latter if not available
##'        in which case it is estimated by QMC
##' @return list with reordered integration limits, scale matrix and Cholesky factor
##' @author Erik Hintz and Marius Hofert
##' @note See Genz and Bretz (2002, p. 957)
precond <- function(lower, upper, scale, cholScale, mean.sqrt.mix)
{
  d <- length(lower)
  y <- rep(0, d - 1)
  
  for(j in 1:(d-1)) {
    
    ## Case j = 1 somewhat special:
    if(j == 1){
      denom <- sqrt(diag(scale))
      c <- 0
    } else {
      denom <- sqrt(diag(scale)[j:d] - rowSums(cholScale[j:d, 1:(j-1), drop = FALSE]^2))
      c <- cholScale[j:d, 1:(j-1), drop = FALSE] %*% y[1:(j-1)]
    }
    
    ## Find i = argmin { <expected length of interval j> }
    i <- which.min(pnorm((upper[j:d] / mean.sqrt.mix - c) / denom) -
                     pnorm((lower[j:d] / mean.sqrt.mix - c) / denom)) + j - 1
    
    
    ## Swap i and j if they are different:
    if(i != j){
      tmp <- swap(i = i, j = j, lower = lower, upper = upper, scale = scale)
      lower <- tmp$lower
      upper <- tmp$upper
      scale <- tmp$scale
      
      ## If j>1 and an actual swap has occured, need to reorder cholesky factor:
      if(j > 1){
        cholScale[c(i,j),]   <- cholScale[c(j,i),, drop = FALSE]
        cholScale[j,(j+1):i] <- matrix(0, ncol = i - j, nrow = 1)
      }
    }
    
    ## Update Cholesky factor
    if(j == 1){
      cholScale[1, 1] <- sqrt(scale[1, 1])
      cholScale[2:d, 1] <- scale[2:d, 1, drop = FALSE] / cholScale[1, 1]
      ## Store y1
      y[1] <- -(dnorm(upper[1]/mean.sqrt.mix) - dnorm(lower[1]/mean.sqrt.mix)) /
        (pnorm(upper[1]/mean.sqrt.mix) - pnorm(lower[1]/mean.sqrt.mix))
    } else {
      cholScale[j,j] <- sqrt(scale[j,j] - sum(cholScale[j,1:(j-1)]^2))
      if(j < d-1) {
        cholScale[(j+1):d, j] <- (scale[(j+1):d, j] - cholScale[(j+1):d, 1:(j-1), drop = FALSE] %*%
                                    cholScale[j, 1:(j-1)]) / cholScale[j, j]
      } else {
        cholScale[(j+1):d, j] <- (scale[(j+1):d, j] - cholScale[(j+1):d, 1:(j-1)] %*%
                                    cholScale[j, 1:(j-1)]) / cholScale[j, j]
      }
      ## Get yj
      low.j.up.j <- c(lower[j] / mean.sqrt.mix - cholScale[j, 1:(j-1)] %*% y[1:(j-1)],
                      upper[j] / mean.sqrt.mix - cholScale[j, 1:(j-1)] %*% y[1:(j-1)]) / cholScale[j, j]
      y[j] <- (dnorm(low.j.up.j[1]) - dnorm(low.j.up.j[2])) / (pnorm(low.j.up.j[2]) - pnorm(low.j.up.j[1]))
    }
  }
  cholScale[d, d] <- sqrt(scale[d, d] - sum(cholScale[d, 1:(d-1)]^2))
  
  ## Return
  list(lower = lower, upper = upper, scale = scale, cholScale = cholScale)
}

##' @title Distribution Function of a Multivariate Normal Variance Mixture for a Single Observation
##' @note Internal function being called by pnvmix. Unless a parameter has a description, 
##'        it is the same as in pnvmix() 
##' @param upper d vector
##' @param lower d vector 
##' @param qW quantile function of W. This has already been built and 'checked' in pnvmix()
##' @param mean.sqrt.mix
##' @param loc
##' @param scale
##' @param cholScale lower triangular cholesky factor of scale. If NULL, it will be calculated.
##' @param method
##' @param precond
##' @param abstol
##' @param CI.factor
##' @param fun.eval
##' @param increment 
##' @param const logical, TRUE if qW is constant (=> dealing with normal dist'n)
##' @param B
##' @param ...
##' @return list of length 3:
##'         - value: computed probability
##'         - error: error estimate
##'         - numiter: number of iterations needed
##' @author Erik Hintz and Marius Hofert
pnvmix1 <- function(upper, lower = rep(-Inf, d), qW = NULL, mean.sqrt.mix,
                    loc = rep(0, d), scale = diag(d), cholScale = NULL, 
                    method = c("sobol", "ghalton", "PRNG"), precond = TRUE, 
                    abstol = 1e-3, CI.factor = 3.3, fun.eval = c(2^6, 1e8), 
                    increment = c("doubling", "num.init"), B = 12, 
                    const = FALSE, ...)
{
  
  ## (Only) basic check; most checking and building was done in pnvmix()
  d <- length(upper)
  stopifnot(length(lower) == d, is.function(qW))
  if(any(lower == upper))
    return(list(value = 0, error = 0, numiter = 0))
  
  ## Get Cholesky factor if not provided (lower triangular)
  ## This is only needed if the internal function pnvmix1() is called directly,
  ## if pnvmix() is used, it was determined there. 
  if(is.null(cholScale)) cholScale <- t(chol(scale)) # get Cholesky factor if not provided (lower triangular)
  
  ## Preconditioning (resorting the limits; only for d > 2)    
  if(precond && d > 2) {
    ## Note that mean.sqrt.mix has already been calculated in pnvmix() 
    temp <- precond(lower = lower, upper = upper, scale = scale,
                    cholScale = cholScale, mean.sqrt.mix = mean.sqrt.mix)
    lower <- temp$lower
    upper <- temp$upper
    cholScale <- temp$cholScale
  }
  
  
  ## Define some basics needed for the while loop below:
  ##
  ## Error is calculated as CI.factor * sd( estimates) / sqrt(B);
  ## replace CI.factor by CI.factor / sqrt(B) to avoid calculating sqrt(B) everytime
  CI.factor <- CI.factor / sqrt(B)
  
  ## Grab the number of sample points for the first iteration
  current.n <- fun.eval[1] #
  
  ## Vector to store the B RQMC estimates
  rqmc.estimates <- rep(0, B)
  
  ## Initialize error to something bigger than abstol so that we can enter the while loop below
  error <- abstol + 42
  ## Initialize the total number of function evaluations
  total.fun.evals <- 0
  ## Initialize a variable that counts the number of iterations in the while loop below
  numiter <- 0
  
  ## It may happen that qnorm(u) for u too close to 1 (or 0) is evaluated; in those
  ## cases, u will be replaced by ONE and ZERO which is the largest (smallest) number
  ## different from 1 (0) such that qnorm(u) is not +/- Inf
  ZERO <- .Machine$double.eps # for symmetry reasons (-8/+8), use this instead of .Machine$double.xmin
  ONE <- 1-.Machine$double.neg.eps
  
  ## If method == sobol, we want the same random shifts in each iteration below,
  ## this is accomplished by reseting to the "original" seed
  if(method == "sobol") {
    if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
    seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used.
  }
  
  ## Additional variables needed if the increment chosen is "doubling"
  if(increment == "doubling") {
    if(method == "sobol") useskip <- 0
    denom <- 1
  }
  
  ## Now enter the while loop, which runs until precision abstol is reached or
  ## the number of function evaluations exceed fun.eval[2]. In each iteration of
  ## the while loop, B RQMC estimates of the desired probability are calculated.
  
  while(error > abstol && total.fun.evals < fun.eval[2])
  {
    if(method == "sobol" && numiter > 0)
      .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
    
    ## Get B RQCM estimates
    for(b in 1:B)
    {
      ## Get the point set:
      
      ## If const = TRUE, we only need (d - 1) (quasi) random numbers
      ## (the case const = TRUE and d = 1 has already been dealt with)
      U <- if(const){
        U <- switch(method,
                    "sobol"   = {
                      if(increment == "doubling") {
                        qrng::sobol(n = current.n, d = d - 1, randomize = TRUE, skip = (useskip * current.n))
                      } else {
                        qrng::sobol(n = current.n, d = d - 1, randomize = TRUE, skip = (numiter * current.n) )
                      }
                    },
                    "ghalton" = {
                      qrng::ghalton(n = current.n, d = d - 1, method = "generalized")
                    },
                    "PRNG"    = {
                      matrix(runif( current.n * (d - 1)), ncol = d - 1)
                    })
        ## First and last column contain 1s corresponding to "simulated" values from sqrt(mix)
        cbind( rep(1, current.n), U, rep(1, current.n))
      } else {
        U <- switch(method,
                    "sobol"   = {
                      if(increment == "doubling") {
                        qrng::sobol(n = current.n, d = d, randomize = TRUE, skip = (useskip * current.n))
                      } else {
                        qrng::sobol(n = current.n, d = d, randomize = TRUE, skip = (numiter * current.n))
                      }
                    },
                    "ghalton" = {
                      qrng::ghalton(n = current.n, d = d, method = "generalized")
                    },
                    "PRNG"    = {
                      matrix(runif( current.n * d), ncol = d)
                    })
        
        ## Case d = 1 somewhat special again:
        if(d == 1){
          cbind( sqrt(qW(U)), sqrt(qW(1 - U)) )
        } else {
          ## Column 1:sqrt(mix), Columns 2--d: unchanged (still uniforms),
          ## Column d + 1: antithetic realization of sqrt(mix)
          cbind(sqrt(qW(U[, 1])), U[, 2:d], sqrt(qW(1 - U[, 1])))
        }
      }
      
      ## Evaluate the integrand at the (next) point set
      
      if(d == 1) {
        ## Case of dimension 1: Don't need to approximate the multivariate
        ##                      normal df and can just use pnorm()
        ## Note that d = 1 for a pure normal or t df has already been addressed
        
        next.estimate <- mean( (pnorm(upper/U[,1])   - pnorm(lower/U[,1]) +
                                  pnorm(upper/U[,d+1]) - pnorm(lower/U[,d+1]) ) / 2)
        
      } else {
        
        next.estimate <- .Call("eval_nvmix_integral",
                               lower     = as.double(lower),
                               upper     = as.double(upper),
                               U         = as.double(U),
                               n         = as.integer(current.n),
                               d         = as.integer(d),
                               cholScale = as.double(cholScale),
                               ZERO      = as.double(ZERO),
                               ONE       = as.double(ONE))
      }
      
      ## Update RQMC estimates
      if(increment == "doubling") {
        ## In this case both, rqmc.estimates[b] and next.estimate depend on n.current points
        rqmc.estimates[b] <- ( rqmc.estimates[b] + next.estimate ) / denom
      } else {
        ## In this case, rqmc.estimates[b] depends on numiter*n.current points whereas
        ## next.estimate depends on n.current points
        rqmc.estimates[b] <- ( numiter * rqmc.estimates[b] + next.estimate ) / (numiter + 1)
      }
      
    } # end for(b in 1:B)
    
    ## Update various variables
    
    ## Number of function evaluations: (* 2 since antithetic variates are used in eval_nvmix_integral())
    total.fun.evals <- total.fun.evals + 2 * B * current.n
    
    ## Double sample size and adjust denominator in averaging as well as useskip
    if(increment == "doubling") {
      ## Change denom and useksip. This is done exactly once, namely in the first iteration.
      if(numiter == 0){
        denom <- 2
        useskip <- 1
      } else {
        ## Increase sample size n. This is done in all iterations except for the first two.
        current.n <- 2 * current.n
      }
    }
    
    ## Update error; note that this CI.factor is actually CI.factor/sqrt(B) (see above)
    error <- CI.factor * sd(rqmc.estimates)
    numiter <- numiter + 1 # update counter
  }
  
  ## Finalize
  value <- mean(rqmc.estimates) # calculate the RQMC estimator
  
  ## Return
  list(value = value, error = error, numiter = numiter)
}

##' @title Distribution Function of a Multivariate Normal Variance Mixture
##' @param upper (n,d) matrix of upper evaluation limits
##' @param lower (n,d) matrix of lower evaluation limits (<= upper)
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param mean.sqrt.mix E(sqrt(W)) or NULL; the latter if not available
##'        in which case it is estimated by QMC
##' @param loc d-vector (location vector)
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param method character string indicating the method to be used:
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         - "prng":    pure Monte Carlo
##' @param precond logical; if TRUE (recommended), variable reordering
##'        similar to Genz and Bretz (2002, pp. 955--956) is performed.
##'        Variable reordering can lead to a significant variance reduction
##'        and decrease in computational time.
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
##' @param increment character string indicating how the sample size should
##'        be increased in each iteration:
##'        - "doubling": next iteration has as many sample points as all the previous
##'          iterations combined
##'        - "num.init": all iterations use an additional fun.eval[1] many points
##' @param B numeric >= 2 providing number of randomizations to get error estimates
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'abstol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return numeric vector with the computed probabilities and attributes "error" 
##'         (error estimate of the RQMC estimator) and "numiter" (number of iterations)
##' @author Erik Hintz and Marius Hofert
pnvmix <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d), qmix, mean.sqrt.mix = NULL,
                   loc = rep(0, d), scale = diag(d), standardized = FALSE,
                   method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
                   abstol = 1e-3, CI.factor = 3.3, fun.eval = c(2^6, 1e8),
                   increment = c("doubling", "num.init"), B = 12, verbose = TRUE, ...)
{
  ## Checks
  if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
  n <- nrow(upper) # number of evaluation points  
  d <- ncol(upper) # dimension
  if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
  
  if(!is.matrix(scale)) scale <- as.matrix(scale)
  stopifnot(dim(lower) == c(n, d), length(loc) == d, # note: 'mean.sqrt.mix' is tested in pnvmix1()
            dim(scale) == c(d, d), is.logical(standardized), is.logical(precond),
            abstol >= 0, CI.factor >= 0, length(fun.eval) == 2, fun.eval >= 0, B >= 1)
  method <- match.arg(method)
  increment <- match.arg(increment)
  
  ## Setting abstol to a negative value will ensure that algorithm runs until fun.eval[2] function
  ## evaluations are done.
  if(is.null(abstol)) abstol <- -42
  
  ## Define the quantile function of the mixing variable.
  const <- FALSE # logical indicating whether we have a multivariate normal
  inv.gam <- FALSE # logical indicating whether we have a multivariate t
  qW <- if(is.character(qmix)) { # 'qmix' is a character vector specifying supported mixture distributions (utilizing '...')
    qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
    switch(qmix,
           "constant" = {
             const <- TRUE
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
               mean.sqrt.mix <- sqrt(df) * gamma(df2) / ( sqrt(2) * gamma( (df+1) / 2 ) ) # used for preconditioning
               function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
             } else {
               const <- TRUE
               mean.sqrt.mix <- 1 # used for preconditioning
               function(u) 1
             }
           },
           stop("Currently unsupported 'qmix'"))
  } else if(is.list(qmix)) { # 'qmix' is a list of the form (<character string>, <parameters>)
    stopifnot(length(qmix) >= 1, is.character(distr <- qmix[[1]]))
    qmix. <- paste0("q", distr)
    if(!existsFunction(qmix.))
      stop("No function named '", qmix., "'.")
    function(u)
      do.call(qmix., append(list(u), qmix[-1])) 
  } else if(is.function(qmix)) { # 'mix' is interpreted as the quantile function F_W^- of the mixture distribution F_W of W
    function(u)
      qmix(u, ...)
  } else stop("'qmix' must be a character string, list or quantile function.")
  
  ## Grab / approximate mean.sqrt.mix, which will be needed for preconditioning in pnvmix1
  ## This only depends on 'qmix', hence it is done (once) here in pnvmix.
  if(precond && d > 2){
    if(is.null(mean.sqrt.mix))
      mean.sqrt.mix <- mean(sqrt(qW(qrng::sobol(n = 2^12, d = 1, randomize = TRUE))))
    ## Check if provided/approximated mean.sqrt.mix is strictly positive 
    if(mean.sqrt.mix <= 0)
      stop("'mean.sqrt.mix' has to be positive (possibly after being generated in pnvmix())")
  }
  
  ## Build temporary result list
  res1 <- vector("list", length = n) # results from calls of pnvmix1()
  
  ## Deal with NA
  NAs <- apply(is.na(lower) | is.na(upper), 1, any) # at least one NA => use NA => nothing left to do
  
  ## Standardize 'scale' if necessary; limits will be standardized later
  if(!standardized){
    Dinv <- diag(1/sqrt(diag(scale))) 
    scale <- Dinv %*% scale %*% Dinv
  }
  ## Get (lower triangular) Cholesky factor of 'scale';
  cholScale <- tryCatch(   t(chol(scale)), error = function(error) error)
  if(inherits(cholScale, "error")) stop(paste("Error occured when calling 'chol(scale)':", cholScale$message))
  
  
  
  ## Loop over observations
  reached <- rep(TRUE, n) # indicating whether 'abstol' has been reached in the ith integration bounds (has to have default TRUE!)
  for(i in seq_len(n)) {
    
    if(NAs[i]) {
      res1[[i]] <- list(value = NA, error = 0, numiter = 0)
      next # => 'reached' already has the right value
    }
    
    ## Pick out ith row of lower and upper
    low <- lower[i,]
    up  <- upper[i,]
    
    ## Deal with equal bounds (result is 0 as integral over null set)
    if(any(low == up)) {
      res1[[i]] <- list(value = 0, error = 0, numiter = 0)
      next # => 'reached' already has the right value
    }
    
    ## Deal with case where components of both low and up are Inf
    lowFin <- is.finite(low)
    upFin  <- is.finite(up)
    lowupFin <- lowFin | upFin # at least one finite limit
    if(any(!lowupFin)) {
      ## Update low, up
      low <- low[lowupFin]
      up  <- up [lowupFin]
      ## Grab (new) dimension. If 0, then all upper are +Inf, all lower are -Inf
      ## => Return 0
      d <- length(low) # Update dimension 
      if(d==0){
        res1[[i]] <- list(value = 1, error = 0, numiter = 0)
        next 
      }
      ## Update scale, Dinv etc
      scale <- scale[lowupFin, lowupFin, drop = FALSE] # update scale
      Dinv <- Dinv[lowupFin, lowupFin, drop = FALSE]  # update Dinv
      cholScaleFin <- t(chol(scale)) # Cholesky factor changes!
    } else {
      ## If no such component exists, set Choleksy factor correctly
      cholScaleFin <- cholScale 
    }
    
    ## Standardize the ith row of lower and upper if necessary
    ## Shift
    if(any(loc != 0)) {
      low <- low - loc
      up  <- up  - loc
    }
    ## Scale
    ## Note that 'scale' has already been standardized 
    if(!standardized) {
      low[lowFin] <- as.vector(Dinv[lowFin, lowFin] %*% low[lowFin]) # only works for those values which are not +/- Inf
      up [upFin]  <- as.vector(Dinv[upFin,  upFin]  %*% up [upFin])
    }
    
    ## If d = 1, deal with multivariate normal or t via pnorm() and pt()
    ## Note that everything has been standardized.  s 
    if(d == 1){
      if(const){
        value <- pnorm(up) - pnorm(low)
        res1[[i]] <- list(value = value, error = 0, numiter = 0)
        next
      }
      if(inv.gam){
        value <- pt(up, df = df) - pt(low, df = df)
        res1[[i]] <- list(value = value, error = 0, numiter = 0)
        next
      }
    }
    
    ## Compute result for ith row of lower and upper (in essence,
    ## because of the preconditioning, one has to treat each such
    ## row separately)
    res1[[i]] <- pnvmix1(up, lower = low, qW = qW, mean.sqrt.mix = mean.sqrt.mix,
                         loc = loc, scale = scale, cholScale = cholScaleFin,
                         method = method, precond = precond, abstol = abstol,
                         CI.factor = CI.factor, fun.eval = fun.eval, 
                         increment = increment, B = B, 
                         const = const, inv.gam = inv.gam, ...)
    
    ## Check if desired precision was reached
    reached[i] <- res1[[i]]$error <= abstol
    if(verbose >= 2 && !reached[i])
      warning(sprintf("'abstol' not reached for pair %d of integration bounds; consider increasing 'fun.eval[2]'", i))
  }
  
  if(verbose == 1 && any(!reached)) { # <=> verbose == 1
    ii <- which(!reached) # (beginning of) indices
    strng <- if(length(ii) > 6) paste0(paste(head(ii), collapse = ", "), ",...") else paste(ii, collapse = ", ")
    warning("'abstol' not reached for pair(s) ",strng," of integration bounds; consider increasing 'fun.eval[2]'")
  }
  
  ## Return
  res <- vapply(res1, function(r) r$value, NA_real_)
  attr(res, "error")   <- vapply(res1, function(r) r$error, NA_real_)
  attr(res, "numiter") <- vapply(res1, function(r) r$numiter, NA_real_)
  res
}
