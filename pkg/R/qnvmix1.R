### qnvmix1() ##################################################################
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
##' @param CI.factor see ?pnvmix()
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


qnvmix1 <- function(u, qmix, 
                    stored.values = NULL,
                    control = list(),
                    #CI.factor = 3.3, n0 = 2^7, B = 8,
                    #abstol.cdf = 1e-4, abstol.logdensity = 1e-2, abstol.newton = 1e-4,
                    #max.iter.newton = 30, max.iter.rqmc = 15, 
                    verbose = TRUE, q.only = FALSE, ...) {
  ## 1 Checking and definitions ################################################
  
  ## Basic input checking (stored.values is checked below)
  stopifnot(!any(u>=1), !any(u<=0), is.logical(verbose), is.logical(q.only))
  
  ## Deal with algorithm parameters: The idea of the following is taken from
  ## how optim() handles parameters.
  ## Default algorithm parameters:   
  control.int <- list(method = "sobol",
                      mean.sqrt.mix = NULL,
                      precond = TRUE,
                      abstol = 1e-3,
                      CI.factor = 3.3,
                      fun.eval = c(2^7, 1e8),
                      max.iter.rqmc = 15,
                      max.iter.newton = 40,
                      increment = "doubling",
                      B = 8,
                      n0 = 2^7,
                      abstol.newton = 1e-4,
                      abstol.logdensity = 1e-2,
                      abstol.cdf = 1e-4)
  names.control <- names(control.int)
  ## Overwrite those defaults by the ones that the user provided (if applicable)
  control.int[(names.provided <- names(control))] <- control
  ## Did they provide something that is not used?
  if (length(unmatched <- names.provided[!names.provided %in% names.control])) 
    warning("unknown names in control: ", paste(unmatched, collapse = ", "))
  
  ## Now some checkings if the arguments provided make sense
  stopifnot(is.logical(control.int$precond), control.int$abstol >= 0,
            control.int$CI.factor >= 0, length(control.int$fun.eval) == 2,
            control.int$fun.eval >= 0, control.int$max.iter.rqmc > 0,
            control.int$max.iter.newton > 0, control.int$B > 1,
            control.int$abstol.newton >= 0, control.int$abstol.logdensity.newton >= 0,
            control.int$abstol.cdf.newton >= 0)
  
  ## Grab method, increment and B
  method    <- match.arg(control.int$method, choices = c("sobol", "ghalton", "PRNG"))
  #increment <- match.arg(control.int$increment, choices = c("doubling", "num.init"))
  B         <- control.int$B
  n0        <- control.int$n0
  
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
  
  ## Initialize first pointset needed for rqmc approach:
  if(method == "sobol"){
    if(!exists(".Random.seed")) runif(1)
    seed <- .Random.seed
  } 
  
  ## Get realizations of W and sqrt(W)
  ## Initial point-set as vector 
  U0 <- switch(method,
               "sobol"   = {
                 as.vector(sapply(1:B, function(i) 
                   sobol(control.int$fun.eval[1], d = 1, randomize = TRUE)))
               },
               "gHalton" = {
                 as.vector(sapply(1:B, function(i) 
                   ghalton(control.int$fun.eval[1], d = 1, method = "generalized")))
               },
               "prng"    = {
                 runif(control.int$fun.eval[1]*B)
               })
  
  mixings <- W(U0)  #length B*fun.eval[1]
  sqrt.mixings <- sqrt(mixings)
  CI.factor <- control.int$CI.factor/sqrt(B) # instead of dividing by sqrt(B) all the time
  
  ## Function that estimates F(x) and logf(x) for *scalar* input x
  est.cdf.dens <- function(x){
    ## Set up result vectors
    rqmc.estimates.log.density <- rep(NA, B)
    rqmc.estimates.cdf <- rep(NA, B)
    ## "realizations" of log density
    log.dens <- - 1/2 * log(2 * pi * mixings) - x*x / (2*mixings) # length B*no
    ## "realizations" of cdf
    cdf <- pnorm(x/sqrt.mixings)
    ## Grab realizations corresponding to l'th shift and use exp-log trick 
    for(l in 1:B){
      log.dens.max <- max( log.dens[ (l-1)*n0 + (1:n0)])
      rqmc.estimates.log.density[l] <- -log(n0) + log.dens.max + 
        log(sum(exp(log.dens[ (l-1)*n0 + (1:n0)] - log.dens.max)))
      rqmc.estimates.cdf[l] <- mean( cdf[(l-1)*n0 + (1:n0)] )
    }
    ## Check if precisions are reached
    error <- c(sd(rqmc.estimates.cdf), sd(rqmc.estimates.log.density))*CI.factor
    precision.reached <- (error[1] <= control.int$abstol.cdf && 
                            error[2] <= control.int$abstol.logdensity)
      
    if(!precision.reached){
      ## Set up while loop
      iter.rqmc <- 0
      current.n <- n0
      
      while(!precision.reached && iter.rqmc < control.int$max.iter.rqmc){
        ## Reset seed and get realizations
        if(method == "sobol") .Random.seed <- seed
        
        U0 <- switch(method,
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
        
        mixings.new <- W(U0) # length B*current.n
        ## "realizations" of log density
        log.dens <- - 1/2 * log(2 * pi * mixings.new) - x*x / (2*mixings.new) 
        ## "realizations" of cdf
        cdf <- pnorm(x/sqrt(mixings.new))
        ## Update RQMC estimators
        for(l in 1:B){
          ## Grab realizations corresponding to l'th shift and use exp-log trick 
          log.dens.max <- max( log.dens[ (l-1)*current.n + (1:current.n)])
          rqmc.estimates.log.density[l] <- (rqmc.estimates.log.density[l] - 
            log(current.n) + log.dens.max + 
            log(sum(exp(log.dens[ (l-1)*current.n + (1:current.n)] - log.dens.max))))/2
          rqmc.estimates.cdf[l] <- (rqmc.estimates.cdf[l] + 
            mean( cdf[(l-1)*current.n + (1:current.n)] ) )/2
        }
        ## Update sample size, iteration number and error(s)
        current.n <- 2*current.n
        iter.rqmc <- iter.rqmc + 1
        error <- c(sd(rqmc.estimates.cdf), sd(rqmc.estimates.log.density))*CI.factor
        precision.reached <- (error[1] <= control.int$abstol.cdf 
                              && error[2] <= control.int$abstol.logdensity)
      }
    }
    if(verbose){
      if(error[1] > control.int$abstol.cdf) warning("'abstol.cdf' not reached; consider increasing 'max.iter.rqmc'")
      if(error[2] > control.int$abstol.logdensity) warning("'abstol.logdensity' not reached; consider increasing 'max.iter.rqmc'")
    }
    ## Return
    c(mean(rqmc.estimates.cdf), mean(rqmc.estimates.log.density))
  }
  
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
    ## Note: est.cdf.dens(0)[1] = 0.5 (i.e. no error here)
    stored.values <- matrix(c(0, est.cdf.dens(0)), nrow = 1)
  } else {
    ## Some very basic checking if stored.values was provided
    stopifnot(is.matrix(stored.values), dim(stored.values)[2] == 3,
              !any(stored.values[,2] > 1 || stored.values[,2] < 0))
  }
  
  ## 2 Main Loop ###############################################################
  for (i in 1:n) {
    ## Initialize error and counter for Newton
    error <- control.int$abstol.newton + 42
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
    while (error > control.int$abstol.newton && iter.newton < control.int$max.iter.newton) {
      ## Update quantile:
      diff.qu <-  current.qu - 
        (current.qu <- current.qu - sign(current.funvals[1]-current.u)*
           exp(log( abs(current.funvals[1]-current.u)) - current.funvals[2]))
      ## Calculate F and logf and store those values
      current.funvals <- est.cdf.dens(current.qu)
      stored.values   <- rbind( stored.values, c(current.qu, current.funvals),
                               deparse.level = 0)
      ## Update error and increase counter
      error <- abs(diff.qu)
      iter.newton <- iter.newton + 1
    }
    ## Safe result
    quantiles[i]       <- current.qu
    log.density[i]     <- current.funvals[2]
    num.iter.newton[i] <- iter.newton
    if(verbose && error > control.int$abstol.newton) warning("'abstol.newton' not reached; consider increasing 'max.iter.newton'")
  }
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