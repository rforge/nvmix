### dgnvmix() ###################################################################

##' @title Density of a Grouped  Normal Variance Mixture for restricted W
##' @param qW function of one variable specifying the quantile functions of W.
##' @param x (n, d) matrix of evaluation points 
##' @param scale.inv (d, d) matrix scale^{-1} 
##' @param lrdet log(sqrt(det(scale)))
##' @param u.left numeric in (0,1)
##' @param u.right numeric in (0,1), > u.left. Density will be estimated
##'         conditional on W being between its 'u.left' and 'u.right' quantile.
##' @param groupings see ?dgnvmix()         
##' @param max.iter.rqmc maximum number of iterations
##' @param return.all logical; if true, matrix (U, qW(U)) also returned.
##' @param control see ?get_set_param() 
##' @return List of three:
##'         $ldensities n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either
##'         relative error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms
##'         (only if return.all = TRUE)
##' @author Erik Hintz and Marius Hofert
densgmix_rqmc <- function(qW, x, scale.inv, lrdet = 0, u.left = 0, u.right = 1, 
                          groupings = 1:d, max.iter.rqmc = 10, return.all = FALSE,
                          control = list())
{
   ## 1. Checks and set-up  ####################################################
   stopifnot(is.function(qW)) # sanity check
   if(!is.matrix(x)) x <- rbind(x) # 1-row matrix if 'x' is a vector
   d <- ncol(x)
   n <- nrow(x) # number of evaluation points
   numgroups <- length(unique(groupings)) # number of groups 
   stopifnot(dim(scale.inv) == c(d, d)) 
   ## Is there an integration region? 
   if(u.left == u.right){
      return(list(ldens = -Inf, error = 0, numiter = 0))
   }
   ## Define various quantities 
   control <- get_set_param(control)
   dblng   <- (control$increment == "doubling")
   B       <- control$B # number of randomizations
   current.n <- control$fun.eval[1] # initial sample size
   numiter   <- 0 # counter for the number of iterations
   total.fun.evals <- 0 # counter for the number of fun evaluations 
   ## Absolte/relative precision?
   if(is.na(control$dnvmix.reltol)) {
      ## Use absolute error
      tol <- control$dnvmix.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error (=> default with tol = 0.01) 
      tol <- control$dnvmix.reltol
      do.reltol <- TRUE
   }
   ## Store and create seeds if 'sobol' is used to get the same shifts later
   if(control$method == "sobol") {
      seeds_ <- sample(1:(1e3*B), B) # B seeds for 'sobol()'
   }
   ## Additional variables needed if the increment chosen is "dblng"
   if(dblng) {
      if(control$method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Matrix to store RQMC estimates
   rqmc.estimates <- matrix(-Inf, ncol = n, nrow = B)
   ## Define trafo-function that maps u to (u.left, u.right) 
   trafo <- function(u) u.left + (u.right - u.left)*u
   ## Initialize 'max.error' to > tol so that we can enter the while loop
   max.error <- tol + 42
   ## Matrix to store (u, W_1(u),...,W_k(u)) 
   if(return.all) {
      max.nrow <- if(dblng) current.n*B*2^(max.iter.rqmc-1) else 
         max.iter.rqmc*B*current.n
      UsWs <- matrix(NA, ncol = numgroups + 1, nrow = max.nrow)
      curr.lastrow <- 0 # counter row-index additional points are being inserted after
   }
   ## Some more constants needed multiple times
   lconst <- -d/2 * log(2*pi) - lrdet
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ZERO <- .Machine$double.neg.eps # avoid evaluation at 0 < ZERO 
   
   ## 2. Main loop #############################################################
   
   ## while() runs until precision abstol is reached or the number of function
   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
   ## the desired log-densities are calculated.
   while(max.error > tol & numiter < max.iter.rqmc &
         total.fun.evals < control$fun.eval[2])
   {
      ## In each randomization ...
      for(b in 1:B) {
         ## Get the point set (*not* sorted here!)
         U <- switch(
            control$method,
            "sobol" = {
               if(dblng) {
                  qrng::sobol(n = current.n, d = 1, randomize = "digital.shift", 
                              seed = seeds_[b], skip = (useskip * current.n))
               } else {
                  qrng::sobol(n = current.n, d = 1, randomize = "digital.shift", 
                              seed = seeds_[b], skip = (numiter * current.n))
               }
            },
            "ghalton" = {
               qrng::ghalton(n = current.n, d = 1, method = "generalized")
            },
            "PRNG" = {
               runif(current.n)
            }) 
         ## Obtain realizations of 1 / sqrt(W_j), j = 1,..,numgroups 
         mixings <- qW((U <- trafo(U)))
         if(!is.matrix(mixings)) mixings <- cbind(mixings) # can happen in ungrouped case
         ## Store if needed
         if(return.all) {
            UsWs[(curr.lastrow + 1) : (curr.lastrow + current.n), ] <- 
               cbind(U, mixings)
            curr.lastrow <- curr.lastrow + current.n
         }
         ## Transform to  1/sqrt(W_j), j=1,..,*d* 
         mixings <- mixings[, groupings, drop = FALSE] # reorder => now d columns
         rt.mix.i <- 1/sqrt(mixings) # 1/sqrt(W_j), j=1,..,d
         ## Compute mahalanobis distances for each row in 'x' 
         mahasq <- sapply(1:n, function(i){
            Dix <- t(rt.mix.i * matrix(x[i, ], ncol = d, nrow = current.n, byrow = TRUE))
            colSums(Dix * (scale.inv %*% Dix))
         }) # (current.n, n) matrix 
         ## Compute values of log h(u) (i.e., the logarithmic integrand)
         lhvals <- lconst - rowSums(log(mixings))/2 - mahasq/2 # (current.n, n) matrix 
         ## Compute next estimates via LogSumExp trick
         next.estimate <- -log(current.n) + logsumexp(lhvals) # n-vector
         ## Update RQMC estimates
         rqmc.estimates[b,] <-
            if(dblng) {
               ## In this case both, rqmc.estimates[b,] and
               ## next.estimate depend on n.current points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b, ]),
                     b = as.double(next.estimate),
                     n = as.integer(n)) - log(denom)
            } else {
               ## In this case, rqmc.estimates[b,] depends on
               ## numiter * current.n points whereas next.estimate
               ## depends on current.n points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b,] + log(numiter)),
                     b = as.double(next.estimate),
                     n = as.integer(n)) - log(numiter + 1)
            }
      } # end for(b in 1:B)
      ## Update of various variables
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(dblng) {
         ## Change denom and useksip (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Total number of function evaluations
      total.fun.evals <- total.fun.evals + B * current.n
      numiter <- numiter + 1
      ## Update error. The following is slightly faster than 'apply(..., 2, var)'
      ldens <- logsumexp(rqmc.estimates) - log(B) # performs better than .colMeans
      vars <- .colMeans((rqmc.estimates - rep(ldens, each = B))^2, B, n, 0)
      error <- if(!do.reltol) { # absolute error
         sqrt(vars)*CI.factor.sqrt.B
      } else { # relative error
         sqrt(vars)/abs(ldens)*CI.factor.sqrt.B
      }
      max.error <- max(error)
   }
   ## 3. Return ################################################################
   
   if(return.all) {
      list(ldensities = ldens, numiter = numiter, error = error,
           UsWs = UsWs[1:curr.lastrow,])
   } else {
      list(ldensities = ldens, numiter = numiter, error = error)
   }
   
}

##' @title Density of a Grouped Normal Variance Mixture
##' @param x (n, d)-matrix of evaluation points
##' @param qmix see ?pgnvmix()
##' @param loc see ?pgnvmix()
##' @param scale see ?pgnvmix()
##' @param factor Cholesky factor (lower triangular matrix) of 'scale';
##'        important here so that det(scale) is computed correctly!
##' @param control list; see ?get_set_param()
##' @param log logical indicating whether the logarithmic density is to be computed
##' @param verbose logical indicating whether warnings shall be thrown.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return n-vector with computed density values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert
dgnvmix <- function(x, groupings = 1:d, qmix, loc = rep(0, d), scale = diag(d),
                    factor = NULL, # needs to be lower triangular!
                    control = list(), log = FALSE, verbose = TRUE, ...)
{
   
   ## 1 Setup ##################################################################
   
   if(!is.matrix(x)) x <- rbind(x)
   d <- ncol(x) # dimension
   if(!is.matrix(scale)) scale <- as.matrix(scale)
   stopifnot(length(loc) == d, dim(scale) == c(d, d))
   verbose <- as.logical(verbose)
   numgroups <- length(unique(groupings))
   ## Deal with algorithm parameters
   control <- get_set_param(control)
   ## If 'factor' provided, compute 'scale' ('factor' not used in grouped case)
   if(!is.null(factor)) scale <- tcrossprod(factor)
   ## Deal with 'qmix' 
   mix_list      <- get_mix_(qmix = qmix, groupings = groupings, 
                             callingfun = "dnvmix", ... ) 
   qW            <- mix_list[[1]] # function(u)
   special.mix   <- mix_list[[2]] # string or NA
   ## Build result object (log-density)
   lres <- rep(-Inf, (n <- nrow(x))) # n-vector of results
   abserror <- rep(NA, n)
   relerror <- rep(NA, n)
   notNA <- rowSums(is.na(x)) == 0
   lres[!notNA] <- NA
   x <- x[notNA,, drop = FALSE] # non-missing data (rows)
   n <- nrow(x) # update 
   numiter <- 0 # initialize counter 
   scale.inv <- solve(scale)
   lrdet     <- log(det(scale))/2
   ## Absolte/relative precision?
   tol <- if(is.na(control$dnvmix.reltol)) { # if 'reltol = NA' use absolute precision
      do.reltol <- FALSE
      control$dnvmix.abstol
   } else { # otherwise use relative precision (default)
      do.reltol <- TRUE
      control$dnvmix.reltol
   }
   ## 2 Actual computation ####################################################
   
   numiter <- 0 # initialize counter 
   ## Deal with the different distributions
   if(numgroups == 1 & !is.na(special.mix)) { # case of a classical NVM dist'n 
      maha2 <- mahalanobis(x, loc, scale) # squared mahalanobis distances 
      lres[notNA] <- switch(
         special.mix,
         "inverse.gamma" = {
            df <- mix_list$param
            lgamma((df + d) / 2) - lgamma(df/2) - (d/2) * log(df*pi) - lrdet - 
               (df+d)/2 * log1p(maha2/df)
         },
         "constant" = {
            -(d/2) * log(2 * pi) - lrdet - maha2/2
         },
         "pareto" = {
            alpha <- mix_list$param
            log(alpha) - d/2*log(2*pi) - lrdet - (alpha+d/2)*log(maha2/2) +
               pgamma(maha2/2, scale = 1, shape = alpha+d/2, log.p = TRUE) + 
               lgamma(alpha+d/2)
         })
      if(!log) lres <- exp(lres) # already exponentiate
      relerror <- rep(0, length(maha2))
      abserror <- rep(0, length(maha2))
   } else {
      ## General case of a grouped or ungrouped normal variance mixture (=> RQMC)
      ## Apply non-adaptive RQMC on all inputs first
      rqmc.obj <- densgmix_rqmc(qW, x = x, scale.inv = scale.inv, u.left = 0,
                                u.right = 1, lrdet = lrdet, groupings = groupings, 
                                max.iter.rqmc = control$dnvmix.max.iter.rqmc.pilot,
                                return.all = TRUE, control = control)
      ## Extract results
      ldens   <- rqmc.obj$ldensities
      numiter <- rep(rqmc.obj$numiter, n)
      error   <- rqmc.obj$error
      if(any(error > tol)){
         ## Call adaptive procedure here
         ## Accuracy not reached for some inputs => use adaptive method there 
         if(control$dnvmix.doAdapt){
            ZERO <- .Machine$double.neg.eps # avoid evaluation at 0 < ZERO 
            ONE <- 1-.Machine$double.neg.eps # avoid evaluation at 1 > ONE  
            notRchd <- which(error > tol)
            x. <- x[notRchd, , drop = FALSE]
            n.notRchd <- nrow(x.) 
            ## Prepare realizations of (u, F_W^\i(u)) needed for all 'x[notRchd, ]'
            UsWs <- rqmc.obj$UsWs
            Us <- UsWs[, 1] # vector!
            ord <- order(Us)
            Us <- c(ZERO, Us[ord], ONE) # 'Us' are sorted and 'ZERO'/'ONE' included
            Ws <- rbind(qW(ZERO), UsWs[ord, -1, drop = FALSE], qW(ONE)) # 'Ws' are sorted 
            stopifnot( (n.UsWs <- nrow(Ws)) == nrow(Us), ncol(Ws) == numgroups) # check
            Ws <- Ws[, groupings, drop = FALSE] # now has exactly 'd' columns 
            n.UsWs <- nrow(Ws) 
            rt.mix.i <- 1/sqrt(Ws) # 1/sqrt(W_j), j=1,..,d
            ## Compute realizations of log h(u) for *all* 'x[notRchd, ]' at once
            mahasq <- sapply(1:n.notRchd, function(i){
               Dix <- t(rt.mix.i * matrix(x.[i, ], ncol = d, nrow = n.UsWs, byrow = TRUE))
               colSums(Dix * (scale.inv %*% Dix))
            }) # (n.UsWs, n.notRchd) matrix of mahalanobis distances for each row in 'x.' 
            lhvals <- # (n.UsWs, n.notRchd) matrix 
               -d/2 * log(2*pi) - lrdet - rowSums(log(Ws))/2 - mahasq/2 
            ## Result objects
            ldens_adapt    <- rep(NA, n.notRchd)
            error_adapt   <- rep(NA, n.notRchd)
            numiters_adapt <- rep(NA, n.notRchd)
            ## Go through all rows in 'x.' (= all columns of 'lhvals') to find 
            ## limits for the adaptive procedure
            for(i in 1:n.notRchd){
               ## Initialize 
               u.left <- NA
               u.right <- NA
               dhvals <- diff(lhvals[, i]) # log h(u_{i+1}) - log h(u_{i})
               if(all(dhvals > 0)) u.right <- 1 # maximum at the right endpoint (or close to it)
               if(all(dhvals < 0)) u.left <- 0 # maximum at left endpoint
               l.max <- lhvals[ (ind.max <- which.max(lhvals[, i])), i] # *observed* maximum
               u.max <- Us[ind.max] 
               ## Tolerance above which RQMC is used 
               l.tol.int.lower <- min(log(control$dnvmix.tol.int.lower), 
                                      l.max - control$dnvmix.order.lower * log(10))
               if(any(lhvals[, i] > l.tol.int.lower)){
                  ## Indices of 'u's so corresponding log h(u) >  l.tol.int.lower
                  ind.gr <- which(lhvals[, i] > l.tol.int.lower)
                  if(is.na(u.left)) u.left <- Us[ind.gr[1]]
                  if(is.na(u.right))
                     u.right <- if(ind.gr[length(ind.gr)] == n.UsWs) 1 else Us[ind.gr[length(ind.gr)]+1]
               } else {
                  ## No obs > threshold 
                  if(!is.na(u.right)){
                     ## 'u.right' was set above => set 'u.left' if not set yet 
                     if(is.na(u.left)){
                        u.left <- if(u.right >= ONE) 1 - 1e-10 else u.max - control$dnvmix.tol.stratlength
                     }
                  } else {
                     ## 'u.right' was not set
                     if(!is.na(u.left)){
                        ## But 'u.left' was 
                        u.right <- if(u.left == 0) 1e-10 else u.max + control$dnvmix.tol.stratlength
                     } else {
                        ## Neither 'u.left' nor 'u.right' was set
                        u.left <- u.max - control$dnvmix.tol.stratlength
                        u.right <- u.max + control$dnvmix.tol.stratlength
                     } 
                  }
               }
               ## Integrate the two regions outside (u.left, u.right) via trapezoidal rules
               ## ... (0, u.left):
               ldens.left <- if(u.left == 0) -Inf else {
                  ## 0 < u.left < 1 => Find obs in (0, u.left)
                  u_sml <- (Us <= u.left)
                  sum_u_sml <- sum(u_sml)
                  if(sum_u_sml > 1) {
                     ## Case 1: We have >1 observations in (0, u.left)
                     last_sml <- which(u_sml)[sum_u_sml]
                     weights <- abs(c(Us[1], diff(Us[1:last_sml])))
                     logsumexp(
                        as.matrix(log(weights) + (c(-Inf, lhvals[1:(last_sml-1), i]) + 
                                                     lhvals[1:last_sml, i])/2, ncol = 1))
                     
                  } else {
                     ## Case 2: No or only one observations in (0, u.left) 
                     log(u.left) + l.tol.int.lower - log(2)
                  }
               }
               ## ... (u.right, 1):
               ldens.right <- if(u.right == 1) -Inf else {
                  ## 0 < u.right < 1 => Find obs in (u.right, 1)
                  u_gtr <- (Us >= u.right)
                  sum_u_gtr <- sum(u_gtr)
                  if(sum_u_gtr > 1) {
                     ## Case 1: We have >1 observations in (u.right, 1)
                     first_gtr <- which(u_gtr)[1]
                     weights <- abs(c(Us[1], diff(Us[first_gtr:n.UsWs])))
                     logsumexp(
                        as.matrix(log(weights) + (c(lhvals[(first_gtr+1):n.UsWs, i], -Inf) + 
                                                     lhvals[first_gtr:n.UsWs, i])/2, ncol = 1))
                  } else {
                     ## Case 2: No or only one observations in (u.right, 1)
                     log1p(-u.right) + l.tol.int.lower - log(2)
                  }
               }
               ## Integrate from 'u.left' to 'u.right' via RQMC  
               ldens.stratum.obj <- 
                  densgmix_rqmc(qW, x = x.[i, , drop = FALSE], lrdet = lrdet,
                                scale.inv = scale.inv,
                                u.left = u.left, u.right = u.right, groupings = groupings,
                                max.iter.rqmc = control$max.iter.rqmc - 
                                   control$dnvmix.max.iter.rqmc.pilot,
                                control = control)
               numiters_adapt[i] <- ldens.stratum.obj$numiter
               error_adapt[i] <- ldens.stratum.obj$error
               ldens_adapt[i] <- logsumexp(rbind(
                  ldens.left, ldens.right, ldens.stratum.obj$ldensities + 
                     log(u.right - u.left), deparse.level = 0))
            }
            ## End of adaptive procedure. Store results
            ldens[notRchd] <- ldens_adapt
            error[notRchd] <- error_adapt
            numiter[notRchd] <- numiter[notRchd] + numiters_adapt
            ## Handle warnings
            if(as.logical(verbose)) { 
               if(any(is.na(error))) {
                  ## At least one error is NA
                  warning("Estimation unreliable, corresponding error estimate NA")
               }
               whichNA <- which(is.na(error))
               if(any(error[setdiff(1:length(error), whichNA)] > tol)) # 'setdiff' needed if 'whichNA' is empty
                  warning("Tolerance not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
            }
         } else if (as.logical(verbose)){
            ## Adaptive method *not* used, print warning 
            warning("Tolerance not reached for all inputs; consider using the adaptive method by setting 'control$dnvmix.doAdapt' to 'TRUE'.")
         }
      }
      ## Compute absolute and relative error on log mu_hat ('ldens' has length(notNA))
      abserror[notNA] <- if(do.reltol){
         relerror[notNA] <- error
         relerror[notNA] * abs(ldens) 
      } else { # error is absolute error
         relerror[notNA] <- error / abs(ldens) 
         error 
      }
      ## Correct results and error if 'log = FALSE'
      if(!log){
         ldens <- exp(ldens)
         ## CI for mu: exp(logmu_hat +/- abserr(logmu_hat))) = (lower, upper)
         ## => compute max. error on mu_hat as max( (upper - mu), (mu - lower) ) 
         relerror[notNA] <- max( (exp(abserror[notNA]) - 1), (1 - exp(-abserror[notNA])) )
         abserror[notNA] <- ldens * relerror[notNA] # ldens already exponentiated 
      }
      ## Grab results, correct 'error' and 'lres' if 'log = FALSE'
      lres[notNA] <- ldens
   }
   
   ## 3. Return ################################################################
   
   ## Note that 'lres' was exponentiated already if necessary.
   attr(lres, "abs. error") <- abserror
   attr(lres, "rel. error") <- relerror
   attr(lres, "numiter") <- numiter
   lres
}