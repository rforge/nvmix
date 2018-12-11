##### fitnvmix() #################################################################



#' Estimate nu using copula density (internal function)
#'
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
#' @return nu.est (scalar of vector of length init.nu); MLE estimate of nu
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
            .Random.seed <- seed 
            qmix. <- function(u) qW(u, nu = nu) # function of u
            - sum(dnvmixcop(U, qmix = qmix., factor = factor, control = control,
                            verbose = verbose, log = TRUE))}
    }
    ## Optimize neg.log.likelihood over nu
    nu.est <- optim(init.nu, fn = neg.log.likelihood.nu,
                    lower = mix.param.bounds[, 1],
                    upper = mix.param.bounds[, 2],
                    method = "L-BFGS-B", control = control.optim)$par
    nu.est
}


#' Estimate nu given loc, scale by maximizing log-likelihood (internal function)
#'
#' @param tX t(X) where X is as in ?fitnmvix()
#' @param qW quantile function of W; must be function(u, nu)
#' @param init.nu initial estimate of nu
#' @param factor cholesky factor of the scale matrix
#' @param control see ?fitnvmix()
#' @param control.optim passed to optim; see ?optim 
#' @param mix.param.bounds see ?fitnvmix()
#' @param verbose see ?fitnvmix()
#' @param inv.gam logical indicating if W is inv.gamma (special case)
#' @param U0 vector of uniforms, see also dnvmix.int
#' @param seed seed used to produce U0
#' @return nu.est (scalar of vector of length init.nu); MLE estimate of nu
#' @author Erik Hintz

estim.nu <- function(tX, qW, init.nu, loc, scale, control, control.optim,
                     mix.param.bounds, verbose, inv.gam, U0, seed){
    factor <- t(chol(scale))
    if(inv.gam){ ## in this case, dnvmix() uses analytical formula for density
        neg.log.likelihood.nu <- function(nu){
            -sum(dnvmix(tX, qmix = "inverse.gamma", loc = loc, scale = scale,
                        df = nu, log = TRUE, verbose = verbose))
        }
    } else {
        ## Get various quantitites passed to dnvmix.int
        z <- forwardsolve(factor, tX - loc, transpose = FALSE)
        maha2.2 <- sort(colSums(z^2)/2)
        lrdet <- sum(log(diag(factor)))
        d <- ncol(factor)
        .Random.seed <- seed # reset seed => monotonicity
        neg.log.likelihood.nu <- function(nu){
            qmix. <- function(u) qW(u, nu = nu) # function of u only
            ## Call dnvmix.int which by default returns the log-density
            ldens.obj <- dnvmix.int(qW = qmix., maha2.2 = maha2.2, lrdet = lrdet,
                                    U0 = U0, d = d, method = control$method, 
                                    abstol = control$dnvmix.abstol,
                                    CI.factor = control$CI.factor, 
                                    fun.eval = control$fun.eval, 
                                    max.iter.rqmc = control$max.iter.rqmc,
                                    B = control$B, seed = seed)
            # Check if error tolerance reached
            if(ldens.obj$error > control$dnvmix.abstol && verbose)
                warning("'dnvmix.abstol' not reached when estimating log-likelihood; consider increasing 'maxiter.rqmc'")
            ## Return -log-density
            -sum(ldens.obj$ldensities)
        }
    }
    ## Optimize neg.log.likelihood over nu
    nu.est <- optim(init.nu, fn = neg.log.likelihood.nu,
                    lower = mix.param.bounds[, 1],
                    upper = mix.param.bounds[, 2],
                    method = "L-BFGS-B", control = control.optim)$par
    nu.est
}

##' @title Fitting Multivariate Normal Variance Mixtures
##'        (only Student-t right now)
##' @param X (n,d) data matrix
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
##' @param fullstep
##' @param ECME.step
##' @param ECME.maxiter
##' @param conv.tol
##' @param zero
##' @return list of three (if qmix = "constant"), otherwise four:
##'         $nu: estimate for nu (omitted if qmix = "constant")
##'         $loc: estimate for the location vector
##'         $scale: estimate for scale matrix
##'         $num.iter: of EM iterations
##' TODO    include option to give names to parameters etc
##' TODO    maybe hide some arguments (CI.factor etc) or have them passed via ...
##' @author Erik Hintz
fitnvmix <- function(X, qmix,
                     mix.param.length = 1, mix.param.bounds = NA,
                     fitting.method = c("cScov", "CopX", "CopU"),
                     ECMEstep = TRUE, control.optim = list(),  control = list(),
                     verbose = TRUE)
{
    ## 0: Initialize various quantities: #######################################
    fitting.method <- match.arg(fitting.method)
    ## Get algorithm specific parameters that are needed often: (=> readability)
    names.control <- names(control)
    if(!any(names.control == "fun.eval")){
        ## 'fun.eval' was *not* provided:  
        control <- get.set.parameters(control)
        control$fun.eval <- c(2^11, 1e8)
    } else {
        control <- get.set.parameters(control)
    }
    method        <- control$method
    dnvmix.abstol <- control$dnvmix.abstol
    CI.factor     <- control$CI.factor
    B             <- control$B
    fun.eval      <- control$fun.eval
    max.iter.rqmc <- control$max.iter.rqmc
    
    ## Get quantile function:
    ## If 'mix' is "constant" or "inverse.gamma", we use the analytical formulas
    ## TODO: Do that for inverse.gamma (not yet done for testing)
    is.const.mix  <- FALSE # logical indicating whether we have a multivariate normal
    inv.gam       <- FALSE # logical indicating whether we have a multivariate t
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
        loc.est <- colMeans(X)
        ## TODO: Do this better (as.matrix to remove attributes, can be done better)
        scale.est <- as.matrix(nearPD(cov(X))$mat) # sample covariance matrix
        return(list(loc = loc.est, scale = scale.est, iter = 0))
    }
    
    ## Check inputs, get dimensions
    ## TODO: More checking (INFs, ...)
    if(!is.matrix(X)) X <- rbind(X)
    notNA <- rowSums(is.na(X)) == 0
    X <- X[notNA,] # non-missing data (rows)
    tX <- t(X)
    n <- nrow(X)
    d <- ncol(X)
    
    ## Check/define parameter bounds on nu
    if(!any(!is.na(mix.param.bounds))){
        mix.param.bounds <- cbind(rep(-Inf, mix.param.length), rep(Inf, mix.param.length))
    } else{
        stopifnot(all.equal( dim(mix.param.bounds), c(mix.param.length, 2)))
    }
    
    ## Get initial pointset that is being reused again and again:
    if(!exists(".Random.seed")) runif(1)
    seed <- .Random.seed
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
    
    ## 1: Initial estimates for nu, loc, scale: ##################################
    
    ## Unbiased estimator for 'loc":
    loc.est <- colMeans(X) 
    ## Sample covariance matrix
    SCov <- as.matrix(nearPD(cov(X))$mat) # TODO do this smarter
    ## Determine maha distance and determinant of 'scale' 
    chol.SCov <- t(chol(SCov)) 
    z <- forwardsolve(chol.SCov, tX - loc.est, transpose = FALSE)
    maha2.2 <- sort(colSums(z^2))/2 # sorted for dnvmix.int
    lrdet <- sum(log(diag(chol.SCov)))
    
    if(fitting.method != "cScov"){
        ## Copula based method => work with pseudo-obs
        pObs <- copula::pobs(X)
        P   <- switch (fitting.method,
                       "CopX" = {as.matrix(nearPD(cor(X))$mat)},
                       "CopU" = {as.matrix(nearPD(cor(pObs))$mat)})
        chol.P <- t(chol(P))
        ## Estimate nu:
        nu.est <- estim.nu.cop(U = pObs, qW = qW, init.nu = rep(1, mix.param.length), 
                               factor = chol.P, control = control,
                               control.optim = control.optim, 
                               mix.param.bounds = mix.param.bounds, 
                               verbose = verbose, inv.gam = inv.gam)
        ## Estimate 'scale' as multiple of SCov; choose 'c' to max loglikelihood
        neg.log.likelihood.c <- if(inv.gam){
            function(c){
                -sum(dnvmix(tX, qmix = "inverse.gamma", loc = loc.est, 
                            factor = sqrt(c)*chol.SCov,
                            df = nu.est, log = TRUE, verbose = verbose))}
        } else {
            qmix. <- function(u) qW(u, nu = nu.est) # function of u
            function(c){
                ldens.obj <- nvmix:::dnvmix.int(qW = qmix., maha2.2 = maha2.2/c, 
                                                lrdet = (lrdet + d/2*log(c)),
                                                U0 = U0, d = d, method = method, 
                                                abstol = dnvmix.abstol,
                                                CI.factor = CI.factor, 
                                                fun.eval = fun.eval, 
                                                max.iter.rqmc = max.iter.rqmc,
                                                B = B, seed = seed)
                ## Check if error tolerance reached
                if(ldens.obj$error > dnvmix.abstol && verbose)
                    warning("'dnvmix.abstol' not reached when estimating log-likelihood; consider increasing 'maxiter.rqmc'")
                ## Return -log-density
                -sum(ldens.obj$ldensities)
            }
        }
        ## Optimize neg.log.likelihood over c
        ## Starting value is E(W) with nu = nu.est
        start.c <- mean(qW(sobol(n = 250, d = 1, randomize = TRUE), nu = nu.est))
        c.est   <- optim(start.c, fn = neg.log.likelihood.c,
                         lower = 0.1, # for stability
                         upper = NA,
                         method = "L-BFGS-B", control = control.optim)$par
        scale.est <- c.est*SCov
    } else {
        ## -loglikelihood as function of param=(nu,c) of length mix.param.length+1
        neg.log.likelihood.init <- function(param){
            if(inv.gam){
                ## In case of inv.gam, a closed formula for the density exists:
                return(-sum(dnvmix(X, qmix = "inverse.gamma", factor = sqrt(param[2])*chol.SCov, 
                                   loc = loc.est, df = param[1], log = TRUE)))
            } else {
                ## Define a qmix function that can be passed to dnvmix()
                qmix. <- function(u) qW(u, nu = param[1:mix.param.length]) # function of u
                c <- param[mix.param.length + 1]
                ## Call dnvmix.int which by default returns the log-density
                ldens.obj <- nvmix:::dnvmix.int(qW = qmix., maha2.2 = maha2.2/c, 
                                                lrdet = (lrdet + d/2*log(c)),
                                                U0 = U0, d = d, method = method, 
                                                abstol = dnvmix.abstol,
                                                CI.factor = CI.factor, 
                                                fun.eval = fun.eval, 
                                                max.iter.rqmc = max.iter.rqmc,
                                                B = B, seed = seed)
                ## Check if error tolerance reached
                if(ldens.obj$error > dnvmix.abstol && verbose)
                    warning("'dnvmix.abstol' not reached when estimating log-likelihood; consider increasing 'maxiter.rqmc'")
                ## Return -log-density
                -sum(ldens.obj$ldensities)
            }
        }
        ## Optimize neg.log.likelihood over (nu, c)
        param <- optim(rep(1, mix.param.length + 1), fn = neg.log.likelihood.init,
                       lower = c(mix.param.bounds[, 1], 0.1),
                       upper = c(mix.param.bounds[, 2], NA),
                       method = "L-BFGS-B", control = control.optim)$par
        ## Grab estimate for nu as well as for the scale matrix
        nu.est    <- param[1:mix.param.length]
        scale.est <- param[mix.param.length + 1] * SCov
    }
    
    ## 2: ECME step: #############################################################
    if(ECMEstep){
        ## Initialize variables for the estimation of the weights
        CI.factor.sqrt.B <- CI.factor / sqrt(B)
        
        ## Function estimating the weights given current values of (nu, loc, scale)
        ## U is the initial point-set. This function is defined in this local
        ## environment so that we don't have to pass non-changing inputs (notably
        ## qW, tX, n, B, d) everytime we call it.
        get.weights <- function(nu, loc, scale, U, abstol = control$weights.abstol){
            ## Get (x-loc)^T scale^{-1} (x-loc)/2, sorted for eval_dnvmix_integrand
            factor <- t(chol(scale))
            z <- forwardsolve(factor, tX - loc, transpose = FALSE)
            maha2.2 <- colSums(z^2)/2
            # need maha2.2 ordered for eval_dnvmix_int:
            ordering.maha2.2 <- order(maha2.2)
            maha2.2 <- maha2.2[ordering.maha2.2]
            ## log(sqrt(det(scale))):
            lrdet <- sum(log(diag(factor)))
            ## Initialize RQMC procedure to estimate the weights
            ## Matrices to store RQMC estimates for weights and density
            rqmc.estimates.density <- matrix(0, ncol = n, nrow = B)
            rqmc.estimates.condexp <- matrix(0, ncol = n, nrow = B)
            error <- abstol + 42 # initialize error to > abstol to enter while loop
            total.fun.evals <- 0
            current.n <- length(U)/B
            ## First pointset that was passed as vector: U has length B * current.n
            ## Realizations of l'th shift are elements (l-1)*current.n + (1:current.n)
            W <- qW(U, nu = nu)
            for(l in 1:B){
                ## Grab realizations corresponding to l'th shift and use exp-log trick
                ## Note: The C function eval_dnvmix_integrand performs the following
                ##  b <- - (d/2) * log(2 * pi) - k/2 * log(W) - lrdet - outer(1/W, maha2.2)
                ##  bmax <- apply(b, 2, max) 
                ##  log(current_n) + bmax + log(colSums(exp(b - rep(bmax, each = current_n))))
                W.current.sorted <- sort(W[(l-1)*current.n + (1:current.n)])
                rqmc.estimates.density[l,] <- .Call("eval_dnvmix_integrand",
                                                    W          = as.double(W.current.sorted),
                                                    maha2_2    = as.double(maha2.2),
                                                    current_n  = as.integer(current.n),
                                                    n          = as.integer(n),
                                                    d          = as.integer(d),
                                                    k          = as.integer(d),
                                                    lrdet      = as.double(lrdet))
                rqmc.estimates.condexp[l,] <- .Call("eval_dnvmix_integrand",
                                                    W          = as.double(W.current.sorted),
                                                    maha2_2    = as.double(maha2.2),
                                                    current_n  = as.integer(current.n),
                                                    n          = as.integer(n),
                                                    d          = as.integer(d),
                                                    k          = as.integer(d+2), 
                                                    ##  Note k = d+2 here!
                                                    lrdet      = as.double(lrdet))
            }
            error <- CI.factor.sqrt.B * max( sd( rqmc.estimates.density[, 1]), # smallest maha index in column 1
                                             sd( rqmc.estimates.condexp[, 1]))
            total.fun.evals <- B * current.n
            ## Main loop
            while(error > abstol && total.fun.evals < fun.eval[2]){
                if(method == "sobol") .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
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
                for(l in 1:B){
                    ## Get realizations of W
                    W.current.sorted <- sort(W[(l-1)*current.n + (1:current.n)])
                    ## See above for details on eval_dnvmix_integrand
                    rqmc.estimates.density[l,] <- (rqmc.estimates.density[l,] +
                                                       .Call("eval_dnvmix_integrand",
                                                             W          = as.double(W.current.sorted),
                                                             maha2_2    = as.double(maha2.2),
                                                             current_n  = as.integer(current.n),
                                                             n          = as.integer(n),
                                                             d          = as.integer(d),
                                                             k          = as.integer(d),
                                                             lrdet      = as.double(lrdet)))/2
                    rqmc.estimates.condexp[l,] <- (rqmc.estimates.condexp[l,] +
                                                       .Call("eval_dnvmix_integrand",
                                                             W          = as.double(W.current.sorted),
                                                             maha2_2    = as.double(maha2.2),
                                                             current_n  = as.integer(current.n),
                                                             n          = as.integer(n),
                                                             d          = as.integer(d),
                                                             k          = as.integer(d+2),
                                                             lrdet      = as.double(lrdet)))/2
                }
                total.fun.evals <- total.fun.evals + B * current.n
                current.n <- 2 * current.n
                ## Compute error measures
                error <- CI.factor.sqrt.B * max( sd( rqmc.estimates.density[, 1]),
                                                 sd( rqmc.estimates.condexp[, 1]))
            }
            ## Return
            weights <- exp(colMeans(rqmc.estimates.condexp) - colMeans(rqmc.estimates.density))
            weights[order(ordering.maha2.2)] ## Recover original ordering
        }
        
        num.iter <- 0
        converged <- FALSE
        converged.nu <- FALSE
        while(num.iter < control$ECME.maxiter && !converged){
            ## Get new weights
            weights <- get.weights(nu.est, loc.est, scale.est, U0)
            ## Get new scale.est: 1/n * sum_{i=1}^n weights_i (X_i-mu)(X_i-mu)^T
            diff.scale <- scale.est -
                (scale.est <- as.matrix(nearPD(cov.wt(X, wt = weights, 
                                                      center = loc.est,
                                                      method = "ML")$cov)$mat))
            ## Get new loc.est: sum_{i=1}^n weights_i X_i / sum weights
            ## as.vector because we need loc.est as a vector, not (d,1) matrix
            diff.loc <- loc.est - 
                (loc.est <- as.vector(crossprod(X, weights)/sum(weights)))
            ## Update nu, if desired/necessary:
            if(!converged.nu && control$ECMEstep.do.nu){
                if(method != "cScov"){
                    ## Correlation matrix:
                    chol.P <- t(chol(cov2cor(scale.est)))
                    nu.est <- estim.nu.cop(pObs, qW = qW, init.nu = nu.est, 
                                           factor = chol.P, control = control, 
                                           control.optim = control.optim,
                                           mix.param.bounds = mix.param.bounds,
                                           verbose = verbose, inv.gam = inv.gam)
                } else {
                    ## Optimize neg.log.likelihood over nu
                    diff.nu <- nu.est -
                        (nu.est <- estim.nu(tX, qW = qW, init.nu = nu.est,
                                            loc = loc.est, scale = scale.est,
                                            control = control,
                                            control.optim = control.optim,
                                            mix.param.bounds = mix.param.bounds,
                                            verbose = verbose, inv.gam = inv.gam,
                                            U0 = U0, seed = seed))
                    converged.nu <- (abs(diff.nu) < control$ECME.conv.tol[3])
                }
            } else {
                diff.nu <- 0
            }
            num.iter <- num.iter + 1
            ## TODO think of something smarter here
            converged <- converged.nu &&
                (sqrt(sum(diff.loc^2)) < control$ECME.conv.tol[1]) && 
                (norm(diff.scale, "F") < control$ECME.conv.tol[2])
        }
        
        ## Another last nu update?
        if(control$laststep.do.nu){
            ## One last nu update
            if(method != "cScov"){
                ## Correlation matrix:
                chol.P <- t(chol(cov2cor(scale.est)))
                nu.est <- estim.nu.cop(pObs, qW, nu.est, chol.P, control, control.optim,
                                       mix.param.bounds, verbose, inv.gam)
            } else {
                nu.est <- estim.nu(tX = tX, qW = qW, loc = loc.est, scale = scale.est,
                                   control = control, control.optim = control.optim,
                                   mix.param.bounds = mix.param.bounds, 
                                   verbose = verbose, inv.gam = inv.gam, U0 = U0,
                                   seed = seed)
            }
        }
    }
    ## 3: Return #################################################################
    list(nu = nu.est, loc = loc.est, scale = scale.est, iter = num.iter)
}







