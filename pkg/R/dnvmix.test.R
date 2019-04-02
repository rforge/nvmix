# # ### dnvmix() ###################################################################
# 
# ##' @title Exp - log trick for log(exp(a) + exp(b)) [not exported]
# ##' @param a numeric vector 
# ##' @param b numeric vector, same length as a
# ##' @return numeric vector log(exp(a) + exp(b)) 
# ##' @author Erik Hintz
# ##' @note NO checking is done for efficiency reasons
# logsumexp <- function(a, b){
#   logsumexp.matrix(rbind(a, b, deparse.level = 0))
# }
# 
# ##' @title Exp - log trick for log(sum_i exp(a_i)) [not exported]
# ##' @param M (n1, n2) matrix
# ##' @return n2-vector log(colSums(exp(M))) 
# ##' @author Erik Hintz
# ##' @note NO checking is done for efficiency reasons
# logsumexp.matrix <- function(M){
#   cmax <- apply(M, 2, max)
#   cmax + log(colSums(exp( M - rep(cmax, each = dim(M)[1]))))
# }
# 
# ##' @title Conditional density of a Multivariate Normal Variance Mixture 
# ##' @param qW function of one variable specifying the quantile function of W.
# ##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
# ##' @param lrdet log(sqrt(det(scale))) where 'scale' is the scale matrix of 
# ##'        the normal variance mixture distribution. 
# ##' @param U0 vector of first input uniforms. Length has to be multiple of B. 
# ##' @param d dimension of the Normal Variance Mixture        
# ##' @param control see ?dnvmix
# ##' @param tol desired tolerance (dnvmix.reltol / dnvmix.abstol)
# ##' @param do.reltol logical; if TRUE, 'tol' is interpreted as relative tolerance
# ##' @param lower.q numeric in (0,1)
# ##' @param upper.q numeric in (0,1), > lower.q. Density will be estimated 
# ##'         conditional on W being between its lower.q and upper.q quantile. 
# ##' @param return.all logical; if true, matrix (U, qW(U)) also returned.       
# ##' @return List of three:
# ##'         $ldensities n-vector with computed log-density values 
# ##'         $error n-vector of error estimates for log-densities; either relative
# ##'         error or absolte error depending on is.na(control$dnvmix.reltol)
# ##'         $numiter number of iterations needed 
# ##'         $rqmc.estimates (B, n) matrix of rqmc estimates for the log-density
# ##' @author Erik Hintz and Marius Hofert
# dnvmix.int.t <- function(qW, maha2.2, lrdet, d, control, lower.q, upper.q, 
#                          max.iter.rqmc, return.all)
# {
#   ## 1 Basics ##################################################################
#   ## Define various quantites:
#   dblng           <- (control$increment == "doubling")
#   B               <- control$B # number of randomizations
#   n               <- length(maha2.2) # sample size 
#   current.n       <- control$fun.eval[1] #initial sample size 
#   numiter         <- 0 # counter for the number of iterations 
#   total.fun.evals <- 0
#   ## Absolte/relative precision?
#   if(is.na(control$dnvmix.reltol)){
#     ## Use absolute error
#     tol <- control$dnvmix.abstol
#     do.reltol <- FALSE
#   } else {
#     ## Use relative error
#     tol <- control$dnvmix.reltol
#     do.reltol <- TRUE
#   }
#   ## Store seed if 'sobol' is used to get the same shifts later:
#   if(control$method == "sobol") {
#     if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
#     seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
#   }
#   ## Additional variables needed if the increment chosen is "doubling"
#   if(dblng) {
#     if(control$method == "sobol") useskip <- 0
#     denom <- 1
#   }
#   ## Matrix to store RQMC estimates
#   rqmc.estimates <- matrix(-Inf, ncol = n, nrow = B) 
#   ## Will be needed a lot:
#   CI.factor.sqrt.B <- control$CI.factor / sqrt(B) 
#   ## Define trafo-function that maps u to (q,1) or (1,q) depending on 'up'
#   trafo <- function(u) lower.q + (upper.q - lower.q)*u
#   ## Initialize 'max.error' to > tol so that we can enter the while loop:
#   max.error <- tol + 42 
#   ## Matrix to store U, W values => nrows = maximal number of funevals 
#   if(return.all){
#     max.nrow <- if(doubling) current.n*B*2^(max.iter.rqmc-1) else max.iter.rqmc*B*current.n
#     UsWs <- matrix(NA, ncol = 2, nrow = max.nrow)
#     curr.lastrow <- 0 # will count row-index additional points are being inserted after
#   } 
#   ## 2 Main loop ###############################################################
#   
#   ## while() runs until precision abstol is reached or the number of function
#   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
#   ## the desired log-densities are calculated.
#   while(max.error > tol && numiter < max.iter.rqmc && 
#         total.fun.evals < control$fun.eval[2]) 
#     {
#     ## Reset seed to have the same shifts in sobol( ... )
#     if(control$method == "sobol" && numiter > 0)
#       .Random.seed <<- seed # reset seed to have the same shifts in sobol( ... )
#     
#     
#     for(b in 1:B){
#       ## 2.1 Get the point set ###########################################
#       U <- sort(switch(control$method,
#                   "sobol" = {
#                     if(dblng) {
#                       qrng::sobol(n = current.n, d = 1,
#                                   randomize = TRUE,
#                                   skip = (useskip * current.n))
#                     } else {
#                       qrng::sobol(n = current.n, d = 1,
#                                   randomize = TRUE,
#                                   skip = (numiter * current.n))
#                     }
#                   },
#                   "ghalton" = {
#                     qrng::ghalton(n = current.n, d = 1,
#                                   method = "generalized")
#                   },
#                   "PRNG" = {
#                     runif(current.n)
#                   })) # sorted for later!
#       
#       ## 2.2 Evaluate the integrand at the (next) point set #############
#       
#       W <- qW(U <- trafo(U)) # realizations of the mixing variable; sorted!
#       if(return.all){
#         UsWs[(curr.lastrow + 1) : (curr.lastrow + current.n + 1)] <- cbind(U, W)
#         curr.lastrow <- curr.lastrow + current.n 
#       }
#       next.estimate <- .Call("eval_dnvmix_integrand",
#                              W          = as.double(sort(W)),
#                              maha2_2    = as.double(maha2.2),
#                              current_n  = as.integer(current.n),
#                              n          = as.integer(n),
#                              d          = as.integer(d),
#                              k          = as.integer(d),
#                              lrdet      = as.double(lrdet))
#       ## 2.3 Update RQMC estimates #######################################
#       rqmc.estimates[b,] <-
#         if(dblng) {
#           ## In this case both, rqmc.estimates[b,] and
#           ## next.estimate depend on n.current points
#           logsumexp(rqmc.estimates[b,], next.estimate) - log(denom)
#           #b(rqmc.estimates[b,] + next.estimate) / denom
#         } else {
#           ## In this case, rqmc.estimates[b,] depends on
#           ## numiter * n.current points whereas next.estimate
#           ## depends on n.current points
#           logsumexp(log(numiter) + rqmc.estimates[b, ], next.estimate) - 
#               log(numiter + 1)
#           # (numiter * rqmc.estimates[b,] + next.estimate) / (numiter + 1)
#         }
#       
#     } # end for(b in 1:B)
#     ## Update of various variables
#     ## Double sample size and adjust denominator in averaging as well as useskip
#     if(dblng) {
#       ## Change denom and useksip (exactly once, in the first iteration)
#       if(numiter == 0){
#         denom <- 2
#         useskip <- 1
#       } else {
#         ## Increase sample size n. This is done in all iterations
#         ## except for the first two
#         current.n <- 2 * current.n
#       }
#     }
#     ## Total number of function evaluations:
#     total.fun.evals <- total.fun.evals + B * current.n
#     numiter <- numiter + 1
#     ## Update error. 
#     vars <- apply(rqmc.estimates, 2, var) # variances
#     errors <- if(!do.reltol){
#       sqrt(vars)*CI.factor.sqrt.B
#     } else {
#       sqrt(vars)/abs(.colMeans(rqmc.estimates, B, n, 0))*CI.factor.sqrt.B
#     }
#     max.error <- max(errors)
#   } # while()
#   
#   ## Finalize 
#   ldensities <- log(.colMeans(exp(rqmc.estimates), B, n, 0))
#   
#   ## 4 Return ##################################################################
#   ret.obj <- if(return.all){
#     list(ldensities = ldensities, numiter = numiter, error = errors,
#          UsWs = UsWs)
#   } else {
#     list(ldensities = ldensities, numiter = numiter, error = errors)
#   }
#   ret.obj
# }
# 
# ##' @title Conditional density of a Multivariate Normal Variance Mixture 
# ##' @param qW function of one variable specifying the quantile function of W.
# ##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
# ##' @param lrdet log(sqrt(det(scale))) where 'scale' is the scale matrix of 
# ##'        the normal variance mixture distribution. 
# ##' @param d dimension of the Normal Variance Mixture        
# ##' @param control see ?dnvmix
# ##' @param verbose see ?dnvmix
# ##' @return List of three:
# ##'         $ldensities n-vector with computed log-density values 
# ##'         $error n-vector of error estimates for log-densities; either relative
# ##'         error or absolte error depending on is.na(control$dnvmix.reltol)
# ##'         $numiter number of iterations needed 
# ##'         $rqmc.estimates (B, n) matrix of rqmc estimates for the log-density
# ##' @author Erik Hintz and Marius Hofert
# dnvmix.internal(qW, maha2.2 = maha2.2, lrdet = lrdet, d = d, control = control, 
#                 verbose = verbose)
# {
#   ## 1 Basics ##################################################################
#   
#   ## Absolte/relative precision?
#   if(is.na(control$dnvmix.reltol)){
#     tol <- control$dnvmix.abstol
#     do.reltol <- FALSE
#   } else {
#     ## Use relative error
#     tol <- control$dnvmix.reltol
#     do.reltol <- TRUE
#   }
#   ## Call RQMC procedure without any stratification 
#   rqmc.obj <- dnvmix.int.t(qW, maha2.2 = maha2.2, lrdet = lrdet, d = d, 
#                          control = control, 
#                          lower.q = 0, upper.q = 1, max.iter.rqmc = 4,
#                          return.all = TRUE)
#   ## Extract results
#   ldens   <- rqmc.obj$ldensities
#   numiter <- rqmc.obj$numiter
#   error   <- rqmc.obj$error
#   if(any(error > tol)){
#     ## Accuracy not reached for at least one 'maha2.2' value
#     ## => Use adaptive approach for those
#     notRchd <- which(error > tol)
#     rqmc.strat.obj <- dnvmix.int.strat(qW, maha2.2[notRchd], lrdet = lrdet, d = d,
#                                        UsWs = rqmc.obj$UsWs, control = contol)
#     ldens[notRchd] <- rqmc.strat.obj$ldensities
#     numiter[notRchd] <- numiter[notRchd] + rqmc.strat.obj$numiter
#     error[notRchd] <- rqmc.strat.obj$error
#   }
#   list(ldensities = ldens, numiter = numiter, error = error)
# }
# 
# 
# 
# ##' @title Density of a Multivariate Normal Variance Mixture
# ##' @param x (n, d)-matrix of evaluation points
# ##' @param qmix specification of the (mixture) distribution of W. This can be:
# ##'        1) a character string specifying a supported distribution (additional
# ##'           arguments of this distribution are passed via '...').
# ##'        2) a list of length at least one; the first argument specifies
# ##'           the base name of an existing distribution which can be sampled
# ##'           with prefix "q", the other elements denote additional parameters
# ##'           passed to this "rmix" random number generator.
# ##'        3) a function being interpreted as the quantile function F_W^-.
# ##' @param loc d-vector (location vector)
# ##' @param scale (d, d)-covariance matrix (scale matrix)
# ##' @param factor Cholesky factor (lower triangular matrix) of 'scale';
# ##'        important here so that det(scale) is computed correctly!
# ##' @param method character string indicating the method to be used:
# ##'         - "sobol":   Sobol sequence
# ##'         - "ghalton": generalized Halton sequence
# ##'         - "PRNG":    pure Monte Carlo
# ##' @param abstol numeric >= 0 providing the absolute precision required.
# ##'        If abstol = 0, algorithm will run until total number of function
# ##'        evaluations exceeds fun.eval[2].
# ##' @param CI.factor Monte Carlo confidence interval multiplier. Algorithm runs
# ##'        CI.factor * (estimated standard error) < abstol. If CI.factor = 3.3
# ##'        (default), one can expect the actual absolute error to be less than
# ##'        abstol in 99.9% of the cases
# ##' @param fun.eval 2-vector giving the initial function evaluations (in the
# ##'        first loop; typically powers of 2) and the maximal number of
# ##'        function evaluations
# ##' @param max.iter.rqmc maximum number of iterations in the RQMC approach        
# ##' @param B number of randomizations to get error estimates.
# ##' @param log logical indicating whether the logarithmic density is to be computed
# ##' @param verbose logical indicating whether a warning is given if the required
# ##'        precision 'abstol' has not been reached.
# ##' @param ... additional arguments passed to the underlying mixing distribution
# ##' @return n-vector with computed density values and attributes 'error'
# ##'         (error estimate) and 'numiter' (number of while-loop iterations)
# ##' @author Erik Hintz and Marius Hofert
# dnvmix.t <- function(x, qmix, loc = rep(0, d), scale = diag(d),
#                      factor = NULL, # needs to be lower triangular!
#                      control = list(), log = FALSE, verbose = TRUE,...)
# {
#   ## Checks 
#   if(!is.matrix(x)) x <- rbind(x)
#   d <- ncol(x) # dimension
#   if(!is.matrix(scale)) scale <- as.matrix(scale)
#   stopifnot(length(loc) == d, dim(scale) == c(d, d))
#   ## Deal with algorithm parameters, see also get.set.parameters():
#   ## get.set.parameters() also does argument checking, so not needed here. 
#   control <- nvmix:::get.set.parameters(control)
#   ## If factor is not provided, determine it here as a *lower* triangular matrix
#   if(is.null(factor)) factor <- t(chol(scale)) # lower triangular
#   
#   ## 1 Define the quantile function of the mixing variable ###################
#   ## If 'mix' is "constant" or "inverse.gamma", we use the analytical formulas
#   is.const.mix <- FALSE # logical indicating whether we have a multivariate normal
#   inv.gam <- FALSE # logical indicating whether we have a multivariate t
#   qW <- if(is.character(qmix)) { # 'qmix' is a character vector
#     qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
#     switch(qmix,
#            "constant" = {
#              is.const.mix <- TRUE
#              function(u) 1
#            },
#            "inverse.gamma" = {
#              if(hasArg(df)) {
#                df <- list(...)$df
#              } else {
#                stop("'qmix = \"inverse.gamma\"' requires 'df' to be provided.")
#              }
#              ## Still allow df = Inf (normal distribution)
#              stopifnot(is.numeric(df), length(df) == 1, df > 0)
#              if(is.finite(df)) {
#                inv.gam <- TRUE
#                df2 <- df / 2
#                mean.sqrt.mix <- sqrt(df) * gamma(df2) / (sqrt(2) * gamma((df+1)/2)) # used for preconditioning
#                function(u) 1 / qgamma(1 - u, shape = df2, rate = df2)
#              } else {
#                is.const.mix <- TRUE
#                mean.sqrt.mix <- 1 # used for preconditioning
#                function(u) 1
#              }
#            },
#            stop("Currently unsupported 'qmix'"))
#   } else if(is.list(qmix)) { # 'mix' is a list of the form (<character string>, <parameters>)
#     stopifnot(length(qmix) >= 1, is.character(distr <- qmix[[1]]))
#     qmix. <- paste0("q", distr)
#     if(!existsFunction(qmix.))
#       stop("No function named '", qmix., "'.")
#     function(u)
#       do.call(qmix., append(list(u), qmix[-1]))
#   } else if(is.function(qmix)) { # 'mix' is the quantile function F_W^- of F_W
#     function(u)
#       qmix(u, ...)
#   } else stop("'qmix' must be a character string, list or quantile function.")
#   ## Build result object (log-density)
#   n <- nrow(x)
#   lres <- rep(-Inf, n) # n-vector of results
#   notNA <- rowSums(is.na(x)) == 0
#   lres[!notNA] <- NA
#   x <- x[notNA,, drop = FALSE] # non-missing data (rows)
#   
#   ## 2 Actual computation ######################################################
#   ## Recall that 'scale' is *lower triangular*. For short, let 'scale' = L
#   ## Solve L * z = x_i - mu for z, so z = L^{-1} * (x_i - mu)   (d vector)
#   ## => z^2 (=> componentwise) = z^T z = (x_i - mu)^T * (L^{-1})^T L^{-1} (x_i - mu)
#   ##                           = z^T z = (x_i - mu)^T * (L L^T )^{-1} (x_i - mu)
#   ##                           = (x_i - mu)^T * scale^{-1} * (x_i - mu)
#   ##                           = quadratic form
#   ## Now do this for *all* x_i simultaneously using that L is lower triangular:
#   ## Forwardsolve: "right-hand sides" of equation must be in the *columns*, thus t(x)
#   z <- forwardsolve(factor, t(x) - loc, transpose = FALSE)
#   maha2 <- colSums(z^2) # = sum(z^T z); n-vector of squared Mahalanobis distances from x to mu w.r.t. scale
#   ## Note: could probably be done with mahalanobis() but unclear how we would
#   ##       get det(scale) then.
#   ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
#   ##                       = log(prod(diag(R))) = sum(log(diag(R)))
#   lrdet <- sum(log(diag(factor)))
#   if(!is.finite(lrdet)) stop(paste("Density not defined for singular 'scale' "))
#   ## Counter
#   numiter <- 0 # initialize counter (0 for 'inv.gam' and 'is.const.mix')
#   ## Deal with the different distributions
#   if(inv.gam) { # multivariate t
#     df.d.2 <- (df + d) / 2
#     lres[notNA] <- lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) -
#       lrdet - df.d.2 * log1p(maha2 / df)
#     if(!log) lres <- exp(lres) # already exponentiate
#     error <- rep(0, length(maha2))
#   } else if(is.const.mix) { # multivariate normal
#     lres[notNA] <- -(d/2) * log(2 * pi) - lrdet - maha2/2
#     if(!log) lres <- exp(lres) # already exponentiate
#     error <- rep(0, length(maha2))
#   } else { 
#     ## General case of a multivariate normal variance mixture (RQMC)
#     ## Prepare inputs for dnvmix.int.t
#     ## Sort maha-distance and divide by 2; store ordering to recover original
#     ## ordering later:
#     ordering.maha <- order(maha2)
#     maha2.2 <- maha2[ordering.maha]/2
#     ## Call internal dnvix (which itself calls C-Code)
#     ests <- nvmix:::dnvmix.internal(qW, maha2.2 = maha2.2, lrdet = lrdet, d = d,
#                          control = control, verbose = verbose)
#     ## Grab results, correct 'error' and 'lres' if 'log = FALSE'
#     lres[notNA] <- ests$ldensities[order(ordering.maha)]
#     error <- if(log){
#       ests$error[order(ordering.maha)]
#     } else {
#       lres <- exp(lres)
#       ests$error[order(ordering.maha)]*pmax(lres[notNA], 1)
#     }
#     numiter <- ests$numiter 
#   }
#   ## Return
#   ## Note that 'lres' was exponentiated already if necessary. 
#   attr(lres, "error")   <- error
#   attr(lres, "numiter") <- numiter
#   lres 
# }
# 





