# require(QRM)
# require(nvmix)
# require(qrng)
# require(copula)
# 
# ## Comments:
# ##
# ## Variant *2* leads to reasonably good estimates of nu (almost comparable
# ## with EM algorithm) and is very fast (only one 2 dim optimization)
# ## Q1:  Can we find a good procedure to estimate 'scale' given that nu is known?
# ## Q2:  If the answer is yes and the algorithm is efficient, maybe can use
# ##      a few iterations (estimate nu => estimate scale given new => estimate
# ##      nu again with new estimate of scale etc)
# ##
# ## Q1 + Q2: YES <3
# ##
# ## Q3:  Can we get the Fisher information or at least the observed Fisher
# ##      information? => numerically approximate it straight from dnmvix?
# ##      Our algorithm is not a 'true' EM algorithm; hence, methods that
# ##      approximate I from EM may not work here.
# ## Q4:  fitnvmix() and dnvmix.() are "large functions" (>1MB). Does that affect
# ##      run-time? if we make it smaller, is it gonna be faster?
# 
# 
# ### fitnvmix() #################################################################
# 
# ##' @title Fitting Multivariate Normal Variance Mixtures
# ##'        (only Student-t right now)
# ##' @param X (n,d) data matrix
# ##' @param qmix character string ("constant", "inverse.gamma") or function. If
# ##'        function, it *has* to be qmix(u, nu) and the length of nu has
# ##'        to be provided in mix.param.length.
# ##' @param mix.param.length
# ##' @param mix.param.bounds either NA or (mix.param.length, 2) matrix, where
# ##'         1st/2nd column corresponds to lower/upper limits for i'th element
# ##'         of nu
# ##'         (eg if W~exp(nu), then mix.param.length = 1 and mix.param.bounds =
# ##'         c(0, Inf))
# ##'         Can have NAs for unrestricted parameters
# ##'         TODO: Re-think if to use NA or Inf etc
# ##'         TODO: Replace '0' by 'zero' internally?
# ##' @param fullstep
# ##' @param ECME.step
# ##' @param ECME.maxiter
# ##' @param conv.tol
# ##' @param zero
# ##' @return list of three (if qmix = "constant"), otherwise four:
# ##'         $nu: estimate for nu (omitted if qmix = "constant")
# ##'         $loc: estimate for the location (here mean) vector
# ##'         $scale: estimate for scale matrix
# ##'         $num.iter: # of EM iterations
# ##' TODO    include option to give names to parameters etc
# ##' TODO    maybe hide some arguments (CI.factor etc) or have them passed via ...
# ##' @author Erik Hintz
# 
# fitnvmix <- function(X, qmix,
#                      mix.param.length = 1, mix.param.bounds = NA,
#                      ECMEstep = TRUE, ECMEstep.do.nu = TRUE,
#                      fitting.method = c("cScov", "CopX", "CopU"),
#                      ECME.maxiter = 25, conv.tol = 1e-3,
#                      control.optim = list(), control = list(),
#                      conv.tol.nu = 1e-2, do.nu = FALSE)
#   {
# 
#   ## 0: Initialize various quantities: #########################################
# 
#   fitting.method <- match.arg(fitting.method)
#   ## Get algorithm specific parameters that are needed often: (=> readability)
#   names.control <- names(control)
#   if(!any(names.control == "fun.eval")){
#     ## 'fun.eval' was *not* provided:
#     control <- get.set.parameters(control)
#     control$fun.eval <- c(2^11, 1e8)
#   }
#   control       <- get.set.parameters(control)
#   method        <- control$method
#   dnvmix.abstol <- control$dnvmix.abstol
#   CI.factor     <- control$CI.factor
#   B             <- control$B
#   fun.eval      <- control$fun.eval
# 
# 
#   ## Get quantile function:
#   ## If 'mix' is "constant" or "inverse.gamma", we use the analytical formulas
#   ## TODO: Do that for inverse.gamma (not yet done for testing)
#   is.const.mix  <- FALSE # logical indicating whether we have a multivariate normal
#   inv.gam       <- FALSE # logical indicating whether we have a multivariate t
#   ## set up qW as function(u, nu)
#   qW <- if(is.character(qmix)) { # 'qmix' is a character vector
#     qmix <- match.arg(qmix, choices = c("constant", "inverse.gamma"))
#     switch(qmix,
#            "constant" = {
#              is.const.mix <- TRUE
#              function(u) 1
#            },
#            "inverse.gamma" = {
#              mix.param.length <- 1
#              mix.param.bounds <- matrix(c(1,NA), ncol = 2)
#              function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
#            },
#            stop("Currently unsupported 'qmix'"))
#   } else if(is.list(qmix)) { # 'mix' is a list of the form (<character string>, <parameters>)
#     ## TODO: Do this!
#     ## stopifnot(length(qmix) >= 1, is.character(distr <- qmix[[1]]))
#     ## qmix. <- paste0("q", distr)
#     ## if(!existsFunction(qmix.))
#     ##   stop("No function named '", qmix., "'.")
#     ## function(u)
#     ##  do.call(qmix., append(list(u), qmix[-1]))
#   } else if(is.function(qmix)) { # 'mix' is the quantile function F_W^- of F_W
#     function(u, nu)
#       qmix(u, nu)
#   } else stop("'qmix' must be a character string, list or quantile function.")
# 
#   ## Case of MVN: MLEs are sample mean and sample cov matrix
#   if(is.const.mix){
#     loc.est <- colMeans(X)
#     ## TODO: Do this better (as.matrix to remove attributes, can be done better)
#     scale.est <- as.matrix(nearPD(cov(X))$mat) # sample covariance matrix
#     return(list(loc = loc.est, scale = scale.est, iter = 0))
#   }
# 
#   ## Check inputs, get dimensions
#   ## TODO: More checking (INFs, ...)
#   if(!is.matrix(X)) X <- rbind(X)
#   notNA <- rowSums(is.na(X)) == 0
#   X <- X[notNA,] # non-missing data (rows)
#   tX <- t(X)
#   n <- nrow(X)
#   d <- ncol(X)
# 
# 
#   ## Grab parameter bounds on nu
#   ## TODO: improve (is.na does not work with matrix)
#   if(!any(!is.na(mix.param.bounds))){
#     mix.param.bounds.lower <- rep(-Inf, mix.param.length)
#     mix.param.bounds.upper <- rep(Inf, mix.param.length)
#   } else{
#     stopifnot(all.equal( dim(mix.param.bounds), c(mix.param.length, 2)))
#     mix.param.bounds.lower <- mix.param.bounds[, 1]
#     mix.param.bounds.upper <- mix.param.bounds[, 2]
#   }
# 
#   ## Get Sobol pointset that is being reused again and again:
#   if(!exists(".Random.seed")) runif(1)
#   seed <- .Random.seed
#   U0 <- switch(method,
#                "sobol"   = {
#                  as.vector(sapply(1:B, function(i)
#                    sobol(control$fun.eval[1], d = 1, randomize = TRUE)))
#                },
#                "gHalton" = {
#                  as.vector(sapply(1:B, function(i)
#                    ghalton(control$fun.eval[1], d = 1, method = "generalized")))
#                },
#                "prng"    = {
#                  runif(control$fun.eval[1]*B)
#                })
# 
#   ## 1: Initial estimates for nu, loc, scale: ##################################
#   loc.est <- colMeans(X) # unbiased estimator for loc
# 
#   ## Sample covariance matrix
#   SCov <- as.matrix(nearPD(cov(X))$mat) # TODO do this smarter
#   ## Determine maha distance and determinant of scale factor:
#   chol.SCov <- t(chol(SCov)) # cholesky factor
#   z <- forwardsolve(chol.SCov, tX - loc.est, transpose = FALSE)
#   maha2.2 <- sort(colSums(z^2))/2 # = (x-mu)^T sig^{-1} (x-mu)/2; sorted for dnvmix.int
#   lrdet <- sum(log(diag(chol.SCov)))
# 
#   if(fitting.method != "cScov"){
#     ## Copula based method => work with pseudo-obs
#     U <- copula::pobs(X)
#     P <- switch (fitting.method,
#                  "CopX" = {as.matrix(nearPD(cor(X))$mat)},
#                  "CopU" = {as.matrix(nearPD(cor(U))$mat)})
#     chol.P <- t(chol(P))
#     ## Set up -log likelihood as a function of nu (given P), based on copula
#     neg.log.likelihood.init.nu <- function(nu){
#       qmix. <- function(u) qW(u, nu = nu) # function of u
#       - sum(dnvmixcop(U, qmix = qmix., factor = chol.P, control = control,
#                   verbose = verbose, log = TRUE))
#     }
#     ## Optimize neg.log.likelihood over nu
#     nu.est <- optim(1, fn = neg.log.likelihood.init.nu,
#                    lower = mix.param.bounds.lower,
#                    upper = mix.param.bounds.upper,
#                    method = "L-BFGS-B", control = control.optim)$par
#     ## Estimate 'scale' as multiple of SCov; choose the factor to max -loglike
#     neg.log.likelihood.c <- function(c){
#       ## Define a qmix function that can be passed to dnvmix()
#       qmix. <- function(u) qW(u, nu = nu.est) # function of u
#       ## Call dnvmix.int which by default returns the log-density
#       ldens.obj <- nvmix:::dnvmix.int(qmix = qmix., maha2.2 = maha2.2/c, lrdet = (lrdet + d/2*log(c)),
#                                       U0 = U0, d = d, method = method, abstol = dnvmix.abstol,
#                                       CI.factor = CI.factor, fun.eval = fun.eval, max.iter.rqmc = max.iter.rqmc,
#                                       B = B, seed = seed)
#       ## Check if error tolerance reached
#       if(ldens.obj$error > abstol && verbose)
#         warning("'dnvmix.abstol' not reached when estimating log-likelihood; consider increasing 'maxiter.rqmc'")
#       ## Return -log-density
#       -sum(ldens.obj$ldensities)
#     }
#     ## Optimize neg.log.likelihood over c
#     ## Starting value is E(W) where nu = nu.est
#     start.c <- mean(qW(sobol(n = 100, d = 1, randomize = TRUE)))
#     c.est <- optim(start.c, fn = neg.log.likelihood.c,
#                    lower = 0.1, # for stability
#                    upper = NA,
#                    method = "L-BFGS-B", control = control.optim)$par
#     scale.est <- c.est*SCov
#   } else {
#     ## Sample covariance matrix
#     SCov <- as.matrix(nearPD(cov(X))$mat) # TODO do this smarter
#     ## Determine maha distance and determinant of scale factor:
#     chol.SCov <- t(chol(SCov)) # cholesky factor
#     z <- forwardsolve(chol.SCov, tX - loc.est, transpose = FALSE)
#     maha2.2 <- sort(colSums(z^2))/2 # = (x-mu)^T sig^{-1} (x-mu)/2; sorted for dnvmix.int
#     lrdet <- sum(log(diag(chol.SCov)))
#     ## -log likelihood as function of param=(nu,c) of length mix.param.length+1
#     ## Will take scale = c * Scov, loc = loc.est
#     neg.log.likelihood.init <- function(param){
#       ## Define a qmix function that can be passed to dnvmix()
#       qmix. <- function(u) qW(u, nu = param[1:mix.param.length]) # function of u
#       c <- param[mix.param.length + 1]
#       ## Call dnvmix.int which by default returns the log-density
#       ldens.obj <- nvmix:::dnvmix.int(qmix = qmix., maha2.2 = maha2.2/c, lrdet = (lrdet + d/2*log(c)),
#                                       U0 = U0, d = d, method = method, abstol = dnvmix.abstol,
#                                       CI.factor = CI.factor, fun.eval = fun.eval, max.iter.rqmc = max.iter.rqmc,
#                                       B = B, seed = seed)
#       ## Check if error tolerance reached
#       if(ldens.obj$error > abstol && verbose)
#         warning("'dnvmix.abstol' not reached when estimating log-likelihood; consider increasing 'maxiter.rqmc'")
#       ## Return -log-density
#       -sum(ldens.obj$ldensities)
#     }
#     ## Optimize neg.log.likelihood over (nu, c)
#     param <- optim(rep(1, mix.param.length + 1), fn = neg.log.likelihood.init,
#                    lower = c(mix.param.bounds.lower, 0.1),
#                    upper = c(mix.param.bounds.upper, NA),
#                    method = "L-BFGS-B", control = control.optim)$par
# 
#     ## Grab estimate for nu as well as for the scale matrix
#     nu.est <- param[1:mix.param.length]
#     scale.est <- param[mix.param.length + 1] * SCov
#   }
# 
# 
#   ## 2: ECME step: #############################################################
#   if(ECMEstep){
#     ## Initialize variables for the estimation of the weights
#     CI.factor <- CI.factor / sqrt(B)
# 
#     ## Function estimating the weights given current values of (nu, loc, scale)
#     ## U is the initial point-set. This function is defined in this local
#     ## environment so that we don't have to pass non-changing inputs (notably
#     ## qW, tX, n, B, d) everytime we call it.
#     get.weights <- function(nu, loc, scale, U, abstol = control$weights.abstol){
#       ## Get (x-loc)^T scale^{-1} (x-loc)
#       factor <- t(chol(scale))
#       z <- forwardsolve(factor, tX - loc, transpose = FALSE)
#       maha2 <- colSums(z^2) # = sum(z^T z); n-vector of squared Mahalanobis distances from x to mu w.r.t. scale
#       ## log(sqrt(det(scale))):
#       lrdet <- sum(log(diag(factor)))
#       ## Initialize RQMC procedure to estimate the weights
#       ## Matrices to store RQMC estimates for weights and density
#       rqmc.estimates.density <- matrix(0, ncol = n, nrow = B)
#       rqmc.estimates.condexp <- matrix(0, ncol = n, nrow = B)
#       error <- abstol + 42 # initialize error to > abstol to enter while loop
#       total.fun.evals <- 0
#       current.n <- length(U)/B
#       ## First pointset that was passed as vector: U has length B * current.n
#       ## Realizations of l'th shift are elements (l-1)*current.n + (1:current.n)
#       current.n <- length(U)/B
#       W <- qW(U, nu = nu)
#       ## as.vector(1/W) %*% t(as.vector(maha2/2)) is what outer(1/W, maha2/2)
#       ## does after checking, so shd be faster than using outer()
#       b <- - (d/2) * log(2 * pi * W) - lrdet - as.vector(1/W) %*% t(as.vector(maha2/2))
#       b. <- -log(W) + b
#       for(l in 1:B){
#         ## TODO maybe call C here.
#         ## Grab realizations corresponding to l'th shift and use exp-log trick
#         bmax  <- apply(b[ (l-1)*current.n + (1:current.n),], 2, max)
#         rqmc.estimates.density[l,] <- -log(current.n) + bmax +
#           log(.Internal(colSums(exp(b[ (l-1)*current.n + (1:current.n), ] - rep(bmax, each = current.n)),
#                                 current.n, n, 0)))
#         bmax. <- apply(b.[(l-1)*current.n + (1:current.n),], 2, max)
#         rqmc.estimates.condexp[l,] <- -log(current.n) + bmax. +
#           log(.Internal(colSums(exp(b.[ (l-1)*current.n + (1:current.n), ] - rep(bmax., each = current.n)),
#                                 current.n, n, 0)))
#       }
# 
#       error <- CI.factor * max( sd( rqmc.estimates.density[, min.maha.index]),
#                                 sd( rqmc.estimates.condexp[, min.maha.index]))
#       total.fun.evals <- B * current.n
# 
#       ## Main loop
#       while(error > abstol && total.fun.evals < max.fun.evals){
#         .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
#         for(l in 1:B){
#           ## Get realizations of W
#           W <- qW(qrng::sobol(current.n, d = 1, randomize = TRUE,
#                               skip = current.n), nu = nu)
# 
#           b <- - (d/2) * log(2 * pi * W) - lrdet - as.vector(1/W) %*% t(as.vector(maha2/2)) # (current.n, n)-matrix, each column corresponds to "one x"
#           b. <- -log(W) + b
#           ## Exp-log trick
#           bmax  <- apply(b, 2, max) # n-vector of maximal b's
#           bmax. <- apply(b., 2, max) # n-vector of maximal b's
# 
#           rqmc.estimates.density[l,] <- (rqmc.estimates.density[l,] - log(current.n) + bmax +
#                                            log(.Internal(colSums(exp(b - rep(bmax, each = current.n)), current.n, n, 0)))) / 2
#           rqmc.estimates.condexp[l,] <- (rqmc.estimates.condexp[l,] - log(current.n) + bmax. +
#                                            log(.Internal(colSums(exp(b. - rep(bmax., each = current.n)), current.n, n, 0)))) / 2
#         }
# 
#         total.fun.evals <- total.fun.evals + B * current.n
# 
#         current.n <- 2 * current.n
# 
#         ## Compute error measures
#         error <- CI.factor * max( sd( rqmc.estimates.density[, min.maha.index]),
#                                   sd( rqmc.estimates.condexp[, min.maha.index]))
#       } # while()
# 
#       ## Return
#       exp(.Internal(colMeans(rqmc.estimates.condexp, B, n, 0))
#           - .Internal(colMeans(rqmc.estimates.density, B, n, 0)))
#     }
# 
#     num.iter <- 0
#     converged <- FALSE
#     converged.nu <- FALSE
# 
#     while(num.iter < ECME.maxiter && !converged){
#       ## Get new weights
#       weights <- get.weights(nu.est, loc.est, scale.est, U0)
# 
#       ## Get new scale.est: 1/n * sum_{i=1}^n weights_i (X_i-mu)(X_i-mu)^T
#       diff.scale <- scale.est -
#         (scale.est <- as.matrix(nearPD(cov.wt(X, wt = weights, center = loc.est,
#                                               method = "ML")$cov)$mat))
# 
#       ## Get new loc.est: sum_{i=1}^n weights_i X_i / sum weights
#       ## as.vector because we need loc.est as a vector, not (d,1) matrix
#       diff.loc <- loc.est - (loc.est <- as.vector(crossprod(X, weights)/sum(weights)))
# 
# 
#       ## Update nu, if desired/necessary:
#       if(!converged.nu && do.nu){
#         ## Determine new maha distance and determinant of scale factor:
#         factor.est <- t(chol(scale.est)) # cholesky factor
#         z <- forwardsolve(factor.est, tX - loc.est, transpose = FALSE)
#         maha2 <- colSums(z^2) # = (x-mu)^T sig^{-1} (x-mu)
#         min.maha.index <- which.min(maha2)
#         lrdet <- sum(log(diag(factor.est)))
# 
# 
#         neg.log.likelihood.em <- function(nu){
#           .Random.seed <- seed
#           qmix. <- function(u) qW(u, nu = nu) # function of u only
#           -sum( dnvmix.(qmix = qmix., maha2 = maha2, lrdet = lrdet, n = n,
#                         U = U0, d = d, min.maha.index = min.maha.index,
#                         log = TRUE, abstol = abstol,
#                         fun.eval = c(2^6, max.fun.evals), B = B))
#         }
# 
#         ## Optimize neg.log.likelihood over nu
#         diff.nu <- nu.est -
#           (nu.est <- optim(nu.est, fn = neg.log.likelihood.em,
#                            lower = mix.param.bounds.lower,
#                            upper = mix.param.bounds.upper, method = "L-BFGS-B",
#                            control = list(factr = factr))$par)
# 
#         converged.nu <- (abs(diff.nu) < conv.tol.nu)
# 
#       } else {
#         diff.nu <- 0
#       }
# 
#       num.iter <- num.iter + 1
# 
#       #TODO think of something smarter here
#       converged <- (abs(diff.nu) < conv.tol) &&
#         (sqrt(sum(diff.loc^2)) < conv.tol) && (norm(diff.scale, "F") < conv.tol)
#     }
# 
# 
#     # One last nu update
#     factor.est <- t(chol(scale.est)) # cholesky factor
#     z <- forwardsolve(factor.est, tX - loc.est, transpose = FALSE)
#     maha2 <- colSums(z^2) # = (x-mu)^T sig^{-1} (x-mu)
#     min.maha.index <- which.min(maha2)
#     lrdet <- sum(log(diag(factor.est)))
# 
# 
#     #neg.log.likelihood.em <- function(nu){
#     #  .Random.seed <- seed
#      # qmix. <- function(u) qW(u, nu = nu) # function of u only
#      # -sum( dnvmix.(qmix = qmix., maha2 = maha2, lrdet = lrdet, n = n, U = U0, d = d,
#      #               min.maha.index = min.maha.index, log = TRUE, abstol = abstol,
#      #               fun.eval = c(2^6, max.fun.evals), B = B))
#     #}
# 
#     #nu.est <- optim(nu.est, fn = neg.log.likelihood.em, lower = mix.param.bounds.lower,
#     #                   upper = mix.param.bounds.upper, method = "L-BFGS-B", control = list(factr = factr))$par
# 
# 
#     neg.log.likelihood.2 <- function(param){
#       .Random.seed <- seed
#       ## Define a qmix function that can be passed to dnvmix()
#       qmix. <- function(u) qW(u, nu = param[1:mix.param.length]) # function of u
#       c <- param[mix.param.length + 1]
#       -sum( dnvmix.(qmix = qmix., maha2 = maha2/c, lrdet = (lrdet + d/2*log(c)),
#                     n = n, U = U0, d = d, min.maha.index = min.maha.index,
#                     log = TRUE, fun.eval = c(2^6, max.fun.evals), verbose = TRUE,
#                     abstol = abstol, B = B))
#     }
# 
#     ## Optimize neg.log.likelihood over (nu, c)
#     param <- optim(c(nu.est, 1), fn = neg.log.likelihood.2,
#                    lower = c(mix.param.bounds.lower, 0.1),
#                    upper = c(mix.param.bounds.upper, NA),
#                    method = "L-BFGS-B", control = list(factr = factr))$par
# 
#     ## Grab estimate for nu als well as for the full scale matrix
#     nu.est <- param[1:mix.param.length]
#     scale.est <- param[mix.param.length + 1] * scale.est
#   }
# 
#   ## 3: Return #################################################################
# 
#   list(nu = nu.est, loc = loc.est, scale = scale.est, iter = num.iter)
# }
# 
# 
# 
# df. <- 2.2
# d. <- 10
# n. <- 500
# loc. <- rep(0,d.)
# A <- matrix(runif(d. * d.), ncol = d.)
# diag_vars <- diag(runif(d., min = 2, max = 5))
# P <- diag_vars %*% cov2cor(A %*% t(A)) %*% diag_vars
# X <- rStudent(n., loc = loc., scale = P, df = df.)
# system.time(My <- fitnvmix(X, qmix = "inverse.gamma", abstol = 1e-3, ECME.maxiter = 20,
#                            max.fun.evals = 98304 ))
# 
# 
# 
# My$nu
# My <- fitnvmix(X, qmix = "inverse.gamma", max.fun.evals = 98304, abstol = 5e-3, ECME.maxiter = 20)
# 
# 
# EM.df <- fit.mst(X, method = "BFGS")$df
# EM$df
# EM$mu
# mean((P  - My$scale)/P)
# mean((P  - EM$Sigma)/P)
# 
# 
# My$nu - EM$df
# norm(My$scale - P, "F")
# norm(EM$Sigma - P , "F")
# 
# 
# mean( (P-scale.est)/P)
# 

