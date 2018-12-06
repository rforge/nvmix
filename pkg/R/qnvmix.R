### qnvmix() ##################################################################
##'
##' @param u vector of probabilities 
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param stored.values matrix with 3 columns of the form [x, F(x), logf(x)] where
##'           F and logf are the df and log-density of the distribution specified
##'           in 'qmix'. 
##'           If provided it will be used to determine starting values for 
##'           the internal newton proceudure. Only very basic checking is done.        
##' @param CI.factor see details in ?qnvmix()
##' @param n0 size of initial pointset
##' @param B see ?pnvmix()
##' @param abstol.cdf abstol to estimate F(x)
##' @param abstol.logdensity abstol to estimate logf(x) for the internal Newton procedure
##' @param abstol.newton abstol for newton
##' @param max.iter.newton max number of iterations for newton *per element*
##' @param max.iter.rqmc max number of iterations to estimate F(x) and logf(x) *per element*
##' @param verbose 
##' @param q.only if TRUE, only quantiles will be returned, ow additional quantites (see return )
##' @param ... see ?pnvmix()
##' @return if q.only is TRUE, vector of same length as u with entries F^{-1}(u)
##'         if q.only is FALSE, a list of 
##'         - $q: Vector of quantiles
##'         - $log.density: Vector log-density values at q
##'         - $computed.values: matrix with 3 columns [x, F(x), logf(x)] 
##'         - $newton.iterations: Vector, element i gives nb of 
##'            Newton iterations needed for u[i]
##' @author Erik Hintz, Marius Hofert and Christiane Lemieux  
##' @note - If only the quantiles are needed, abstol.logdensity does not need to be as small.
##'       - *TODO* Maybe do something smart in the tails
##'       - *TODO* print better warnings if accuracies not reached and give better
##'       'verbose' argument (as is done in pnvmix, telling the user which u[i] is
##'       the one that causes problems). 
##'       - *TODO* Maybe one can make better use of stored values if one changes
##'       the order in which the quantiles are evaluated. 
##'       - *TODO* provide rmix option. (=> later, not important)
##'       - *TODO* Is it worth having another argument, specifying if the actual
##'       density shall be returned (instead of logdensity) if  q.only = FALSE?


qnvmix <- function(u, qmix, control = list(),
                   verbose = TRUE, q.only = FALSE, stored.values = NULL, ...) {
  
  ## 1 Checking and definitions ################################################
  
  ## Basic input checking (stored.values is checked below)
  stopifnot(!any(u>=1), !any(u<=0), is.logical(verbose), is.logical(q.only))
  
  ## Deal with algorithm parameters, see also get.set.parameters():
  ## get.set.parameters() also does argument checking, so not needed here. 
  control <- get.set.parameters(control)
  
  ## Grab method, B and n0
  method    <- control$method
  B         <- control$B
  n0        <- control$fun.eval[1]
  
  ## Define the quantile function of the mixing variable: 
  is.const.mix <- FALSE # logical indicating whether we have a multivariate normal
  inv.gam <- FALSE # logical indicating whether we have a multivariate t
  W <- if(is.character(qmix)) { # 'qmix' is a character vector
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
  
  ## Special case of normal and t distribution:
  if(is.const.mix || inv.gam){
    q <- if(inv.gam) qt(u, df = df) else qnorm(u)
    if(q.only) return(q) else{
      l.dens <- if(inv.gam) dt(q, df = df, log = TRUE) else dnorm(q, log = TRUE)
      return(list(q = q, 
                  log.density = l.dens, 
                  computed.values = cbind(q, u, l.dens, deparse.level = 0),
                  newton.iterations = rep(0, length(u))))
    }
  }
  
  ## 2 Set up est.cdf.dens() ## ################################################
  ## Initialize first pointset needed for rqmc approach:
  if(method == "sobol"){
    if(!exists(".Random.seed")) runif(1)
    seed <- .Random.seed
  } 
  
  ## Get realizations of W and sqrt(W)
  ## Initial point-set with B columns (each column = one shift)
  U0 <- switch(method,
               "sobol"   = {
                 sapply(1:B, function(i) 
                   sobol(n0, d = 1, randomize = TRUE))
               },
               "gHalton" = {
                 sapply(1:B, function(i) 
                   ghalton(n0, d = 1, method = "generalized"))
               },
               "prng"    = {
                 matrix(runif(B*n0), ncol = B)
               }) # (n0, B) matrix
  
  mixings <- W(U0)  # (n0, B) matrix
  sqrt.mixings <- sqrt(mixings) # (n0, B) matrix
  
  ## Set up various quantities for est.cdf.dens():
  CI.factor <- control$CI.factor/sqrt(B) # instead of dividing by sqrt(B) all the time
  current.n <- n0 # will give ncol(mixings) 
  
  ## Function that estimates F(x) and logf(x) for *scalar* input x. Similar to
  ## pnvmix() and dnvmix() with increment = "num.init", tailored to the 
  ## one - dimensional case. Previous realizations of the mixing variable are
  ## reused (=> arguments mixings and sqrt.mixings)
  est.cdf.dens <- function(x, mixings, sqrt.mixings){
    ## Define various quantities:
    xx2 <- x*x/2 
    rqmc.estimates.log.density <- rep(NA, B) #  result vector
    rqmc.estimates.cdf <- rep(NA, B)   # result vector
    current.n <- dim(mixings)[1]
    
    ## First we use 'mixings' and 'sqrt.mixings' that are already available.
    for(l in 1:B){
      ## Grab realizations corresponding to l'th shift and use exp-log trick 
      log.dens <- -1/2 * log(2* pi * mixings[,l]) - 
        xx2 / mixings[,l] # length current.n
      log.dens.max <- max(log.dens)
      rqmc.estimates.log.density[l] <- -log(current.n) + log.dens.max + 
        log(sum(exp(log.dens - log.dens.max)))
      rqmc.estimates.cdf[l] <- mean( pnorm(x/sqrt.mixings[,l]) )
    }
    ## Check if precisions are reached
    error <- c(sd(rqmc.estimates.cdf), sd(rqmc.estimates.log.density))*CI.factor
    precision.reached <- (error[1] <= control$newton.df.abstol && 
                            error[2] <= control$newton.logdens.abstol)
    
    if(!precision.reached){
      ## Set up while loop
      iter.rqmc <- 1
      while(!precision.reached && iter.rqmc < control$max.iter.rqmc){
        ## Reset seed and get n0 realizations:
        if(method == "sobol") .Random.seed <- seed
        
        U.next <- switch(method,
                         "sobol"   = {
                           sapply(1:B, function(i) 
                             sobol(n0, d = 1, randomize = TRUE))
                         },
                         "gHalton" = {
                           sapply(1:B, function(i) 
                             ghalton(n0, d = 1, method = "generalized"))
                         },
                         "prng"    = {
                           matrix(runif(B*n0), ncol = B)
                         })
        mixings.next <- W(U.next) # (n0, B) matrix
        sqrt.mixings.next <- sqrt(mixings.next ) # (n0, B) matrix
        
        ## Update RQMC estimators
        for (l in 1:B) {
          ## Grab realizations corresponding to l'th shift and use exp-log trick
          log.dens <- -1/2 * log(2*pi*mixings.next[, l]) - 
            xx2/mixings.next[, l] # length n0
          log.dens.max <- max(log.dens)
          ## Previos estimate based on current.n samples, new one based on n0 samples
          
          rqmc.estimates.log.density[l] <- (current.n*rqmc.estimates.log.density[l] +
                                              n0*(-log(n0) + log.dens.max + log(sum(exp(log.dens - log.dens.max)))))/(current.n + n0) 
          rqmc.estimates.cdf[l] <- (current.n * rqmc.estimates.cdf[l] + 
                                      sum(pnorm(x/sqrt.mixings.next[,l])))/(current.n + n0)
        }
        ## Update mixings and sqrt.mixings so that they can be reused:
        mixings <- rbind(mixings, mixings.next)
        sqrt.mixings <- rbind(sqrt.mixings, sqrt.mixings.next)
        current.n <- current.n + n0
        ## Update iteration number and error(s)
        iter.rqmc <- iter.rqmc + 1
        error <- c(sd(rqmc.estimates.cdf), sd(rqmc.estimates.log.density))*CI.factor
        precision.reached <- (error[1] <= control$newton.df.abstol
                              && error[2] <= control$newton.logdens.abstol)
      }
    }
    if(verbose){
      if(error[1] > control$newton.df.abstol) warning("'abstol.cdf' not reached; consider increasing 'max.iter.rqmc'")
      if(error[2] > control$newton.logdens.abstol) warning("'abstol.logdensity' not reached; consider increasing 'max.iter.rqmc'")
    }
    ## Return
    return(list(estimates = c(mean(rqmc.estimates.cdf), mean(rqmc.estimates.log.density)),
                mixings = mixings, # return 'new' mxinings (= old mixings and new mixings)
                sqrt.mixings = sqrt.mixings))
  }
  
  ## 3 Actual computation of quantiles (Newton's method) #######################
  
  ## Only compute quantiles for u>=0.5 (use symmetry in the other case)
  lower <- (u < 0.5)
  u[lower] <- 1 - u[lower]
  ## Sort u. Ordering needed to return the quantiles in the correct order later
  ordering <- order(u)
  u.sorted <- u[ordering]
  n <- length(u.sorted)
  ## Build result vectors
  quantiles <- rep(NA, n)
  log.density <- rep(NA, n)
  num.iter.newton <- rep(0, n)
  
  if(is.null(stored.values)){
    ## Matrix of the form [x, F(x), logf(x)] to store df and density evaluations
    ## First element corresponds to F(x) = 0.5, so that x = 0 
    ## Note: est.cdf.dens(0)[1] = 0.5 (i.e. no error here), but need log-density
    cdf.dens.mixings <- est.cdf.dens(0, mixings, sqrt.mixings)
    ## Store these values in stored.values
    stored.values <- matrix(c(0, cdf.dens.mixings$estimates), nrow = 1)
    ## Update mixings and sqrt.mixings
    mixings       <- cdf.dens.mixings$mixings
    sqrt.mixings  <- cdf.dens.mixings$sqrt.mixings
  } else {
    ## Some very basic checking if stored.values was provided
    stopifnot(is.matrix(stored.values), dim(stored.values)[2] == 3,
              !any(stored.values[,2] > 1 || stored.values[,2] < 0))
  }
  
  ## Main loop for Newton's method:
  for (i in 1:n) {
    ## Initialize error and counter for Newton
    error <- control$newton.conv.abstol+ 42
    iter.newton <- 0
    ## Grab current u
    current.u <- u.sorted[i]
    ## Did we already have that u? (could happen if, e.g., u = c(0.4,0.6))
    if(i > 1 && current.u == u.sorted[i-1]){
      quantiles[i]   <- quantiles[i-1]
      log.density[i] <- log.density[i-1]
      next
    }
    ## Get starting value: x in stored.values sth F(x) is close to u
    closest.ind <- which.min(abs(stored.values[, 2] - current.u))
    current.qu <- stored.values[closest.ind, 1]
    current.funvals <- stored.values[closest.ind, 2:3]
    
    ## Main loop for Newton procedure
    while (error > control$newton.conv.abstol && iter.newton < control$max.iter.newton) {
      ## Update quantile:
      diff.qu <-  current.qu - 
        (current.qu <- current.qu - sign(current.funvals[1]-current.u)*
           exp(log( abs(current.funvals[1]-current.u)) - current.funvals[2]))
      ## Call est.cdf.dens:
      cdf.dens.mixings <- est.cdf.dens(current.qu, mixings, sqrt.mixings)
      current.funvals  <-  cdf.dens.mixings$estimates
      ## Store these values in stored.values
      stored.values    <- rbind( stored.values, c(current.qu, current.funvals),
                                 deparse.level = 0)
      ## Update mixings and sqrt.mixings
      mixings <- cdf.dens.mixings$mixings
      sqrt.mixings <- cdf.dens.mixings$sqrt.mixings
      ## Update error and increase counter
      error <- abs(diff.qu)
      iter.newton <- iter.newton + 1
    }
    ## Safe result
    quantiles[i]       <- current.qu
    log.density[i]     <- current.funvals[2]
    num.iter.newton[i] <- iter.newton
    if(verbose && error > control$newton.conv.abstol) 
      warning("'abstol.newton' not reached; consider increasing 'max.iter.newton'")
  }
  
  ## 4 Clean-up and return   ###################################################  
  
  ## Order results according to original ordering of input u
  quantiles <- quantiles[order(ordering)]
  ## Use symmetry for those u which were < 0.5
  quantiles[lower] <- -quantiles[lower]
  ## Return
  if(q.only) return(quantiles) else{
    num.iter.newton <- num.iter.newton[order(ordering)]
    log.density     <- log.density[order(ordering)]
    return(list(q = quantiles, 
                log.density = log.density,
                computed.values = stored.values,
                newton.iterations = num.iter.newton))
  }
}