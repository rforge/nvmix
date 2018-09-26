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
    temp_i <- as.matrix(scale[wo.ij,i])
    temp_j <- as.matrix(scale[wo.ij,j])
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
##' @param mean.sqrt.mix E(sqrt(W))
##' @return list with reordered integration limits, scale matrix and Cholesky factor
##' @author Erik Hintz and Marius Hofert
##' @note See Genz and Bretz (2002, p. 957)
precond <- function(lower, upper, scale, cholScale, mean.sqrt.mix)
{
    d <- length(lower)
    y <- rep(0, d - 1)

    ## Find i = argmin_j { <expected length of interval> }
    i <- which.min(apply(pnorm(cbind(upper, lower) / (mean.sqrt.mix * sqrt(diag(scale)))), 1, diff))
    if(i != 1) {
        ## Swap 1 and i in 'lower', 'upper' and 'scale'
        tmp <- swap(i = 1, j = i, lower = lower, upper = upper, scale = scale)
        lower <- tmp$lower
        upper <- tmp$upper
        scale <- tmp$scale
    }

    ## Store y1
    y[1] <- -(dnorm(upper[1]/mean.sqrt.mix) - dnorm(lower[1]/mean.sqrt.mix)) /
        (pnorm(upper[1]/mean.sqrt.mix) - pnorm(lower[1]/mean.sqrt.mix))

    ## Update the Cholesky factor
    cholScale[1, 1] <- sqrt(scale[1, 1])
    cholScale[2:d, 1] <- as.matrix(scale[2:d, 1] / cholScale[1, 1])
    for(j in 2:(d-1)) {
        denom <- sqrt(diag(scale)[j:d] - rowSums(as.matrix(cholScale[j:d, 1:(j-1)])^2))
        c <- as.matrix(cholScale[j:d, 1:j-1]) %*% y[1:(j-1)]
        ## Find i = argmin { <expected length of interval j> }
        i <- which.min(pnorm((upper[j:d] / mean.sqrt.mix - c) / denom) -
                       pnorm((lower[j:d] / mean.sqrt.mix - c) / denom)) + j - 1
        if(i != j){ # swap i and j
            tmp <- swap(i = i, j = j, lower = lower, upper = upper, scale = scale)
            lower <- tmp$lower
            upper <- tmp$upper
            scale <- tmp$scale
            cholScale[c(i,j),]   <- as.matrix(cholScale[c(j,i),])
            cholScale[j,(j+1):i] <- as.matrix(0, ncol = i - j, nrow = 1)
        }
        ## Update Cholesky factor
        cholScale[j,j] <- sqrt(scale[j,j] - sum(cholScale[j,1:(j-1)]^2))
        if(j < d-1)
            cholScale[(j+1):d, j] <- (scale[(j+1):d, j] - as.matrix(cholScale[(j+1):d, 1:(j-1)]) %*%
                                       cholScale[j, 1:(j-1)]) / cholScale[j, j]
        else cholScale[(j+1):d, j] <- (scale[(j+1):d, j] - cholScale[(j+1):d, 1:(j-1)] %*%
                                        cholScale[j, 1:(j-1)]) / cholScale[j, j]
        ## Get yj
        low.j.up.j <- c(lower[j] / mean.sqrt.mix - cholScale[j, 1:(j-1)] %*% y[1:(j-1)],
                        upper[j] / mean.sqrt.mix - cholScale[j, 1:(j-1)] %*% y[1:(j-1)]) / cholScale[j, j]
        y[j] <- (dnorm(low.j.up.j[1]) - dnorm(low.j.up.j[2])) / (pnorm(low.j.up.j[2]) - pnorm(low.j.up.j[1]))
    }
    cholScale[d, d] <- sqrt(scale[d, d] - sum(cholScale[d, 1:(d-1)]^2))

    ## Return
    list(lower = lower, upper = upper, scale = scale, cholScale = cholScale)
}

##' @title Distribution Function of the Multivariate t Distribution
##' @param upper d-vector of upper evaluation limits
##' @param lower d-vector of lower evaluation limits
##' @param mix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param mean.sqrt.mix expectation of sqrt(W)
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
##'        as described in Genz and Bretz (2002, pp. 955--956) is performed.
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
##' @param B number of randomizations to get error estimates.
##' @param ... additional arguments passed to the underlying mixing distributions
##' @return list with the
##'         - computed probability
##'         - total number of function evaluations
##'         - number of iterations in the while loop
##'         - error estimate
##'         - variance estimate
##' @author Erik Hintz
pnvmix <- function(upper, lower = rep(-Inf, d), mix, mean.sqrt.mix = NULL,
                   loc = rep(0, d), scale = diag(d), standardized = FALSE,
                   method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
                   abstol = 1e-3, CI.factor = 3.3, fun.eval = c(2^6, 1e8), B = 12, ...)
{
    ## Checks
    d <- length(upper) # dimension of the problem
    if(!is.matrix(scale)) scale <- as.matrix(scale)
    stopifnot(length(lower) == d, lower < upper, length(loc) == d, # note: mean.sqrt.mix tested later
              dim(scale) == c(d, d), is.logical(standardized), is.logical(precond),
              abstol >= 0, CI.factor >= 0, length(fun.eval) == 2, fun.eval >= 0, B >= 1)
    method <- match.arg(method)

    ## Deal with infinite limits
    lowFin <- is.finite(lower)
    upFin  <- is.finite(upper)
    lowupFin <- lowFin | upFin # at least one finite limit
    if(any(!lowupFin)) {
        lower <- lower[lowupFin]
        upper <- upper[lowupFin]
        scale <- scale[lowupFin, lowupFin, drop = FALSE]
        d <- ncol(scale)
    }

    ## Standardize if necessary
    ## Shift
    if(any(loc != 0)) {
        lower <- lower - loc
        upper <- upper - loc
    }
    ## Scale
    if(!standardized) {
        Dinv <- diag(1/sqrt(diag(scale))) # TODO: why not work with cov2cor()?; these 4 lines can most likely be improved
        scale <- Dinv %*% scale %*% Dinv
        lower[lowFin] <- as.vector(Dinv[lowFin, lowFin] %*% lower[lowFin]) # only works for those values which are not +/- Inf
        upper[upFin]  <- as.vector(Dinv[upFin,  upFin]  %*% upper[upFin])
    }

    ## Define the quantile function of the mixing variable
    const <- FALSE # logical indicating whether we have a multivariate normal
    inv.gam <- FALSE # logical indicating whether we have a multivariate t
    W <- if(is.character(mix)) { # 'mix' is a character vector specifying supported mixture distributions (utilizing '...')
             mix <- match.arg(mix, choices = c("constant", "inverse.gamma"))
             switch(mix,
                    "constant" = {
                        const <- TRUE
                        function(u) 1
                    },
                    "inverse.gamma" = {
                        if(hasArg(df)) df <- list(...)$df else
                            stop("'mix = \"inverse.gamma\"' requires 'df' to be provided.")
                        ## Still allow df = Inf (normal distribution)
                        stopifnot(is.numeric(df), length(df) == 1, df > 0)
                        if(is.finite(df)) {
                            inv.gam <- TRUE
                            df2 <- df / 2
                            mean.sqrt.mix <- sqrt(df) * gamma(df2) / ( sqrt(2) * gamma( (df+1) / 2 ) ) # used for preconditioning
                            function(u) 1 / qgamma(u, shape = df2, rate = df2)
                        } else {
                            const <- TRUE
                            mean.sqrt.mix <- 1 # used for preconditioning
                            function(u) 1
                        }
                    },
                    stop("Currently unsupported 'mix'"))
         } else if(is.list(mix)) { # 'mix' is a list of the form (<character string>, <parameters>)
             stopifnot(length(mix) >= 1, is.character(distr <- mix[[1]]))
             qmix <- paste0("q", distr)
             if(!existsFunction(qmix))
                 stop("No function named '", qmix, "'.")
             function(u){
                 return(  do.call(qmix, c(u, mix[-1])))
             }
         } else if(is.function(mix)) { # 'mix' is interpreted as the quantile function F_W^- of the mixture distribution F_W of W
             function(u){
                 return(mix(u, ...))
             }
         } else stop("'mix' must be a character string, list or quantile function.")

    ## If d = 1, deal with multivariate normal or t via pnorm() and pt()
    if(d == 1){
        if(const){
            Prob <- pnorm(upper) - pnorm(lower)
            return(list(Prob = Prob, N = 0, i = 0, ErrEst = 0, Var = 0))
        }
        if(inv.gam){
            Prob <- pt(upper, df = df) - pt(lower, df = df)
            return(list(Prob = Prob, N = 0, i = 0, ErrEst = 0, Var = 0))
        }
    }

    ## Preconditioning (resorting the limits; only for d > 2)
    cholScale <- t(chol(scale)) # get Cholesky factor (lower triangular)
    if(precond && d > 2) {
        if(is.null(mean.sqrt.mix)) # approximate E(sqrt(W))
            mean.sqrt.mix <- mean(sqrt(W(qrng::sobol(n = 5000, d = 1, randomize = TRUE))))
        if(any(mean.sqrt.mix <= 0))
            stop("'mean.sqrt.mix' has to be positive (possibly after being generated in pnvmix())")
        temp <- precond(lower = lower, upper = upper, scale = scale,
                        cholScale = cholScale, mean.sqrt.mix = mean.sqrt.mix)
        lower <- temp$lower
        upper <- temp$upper
        cholScale <- temp$cholScale
    }

    ## TODO: put meaningful header here (and in what follows)
    CI.factor <- CI.factor / sqrt(B) # instead of dividing sigma by sqrt(B) each time
    n. <- fun.eval[1] # initial n
    T. <- rep(0, B) # vector to store RQMC estimates

    ZERO <- .Machine$double.eps # TODO: shouldn't this be double.xmin? see ?.Machine
    ONE <- 1-.Machine$double.neg.eps

    if(method == "sobol") {
        if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
        seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used.
    }

    err <- abstol + 42 # initialize err to something bigger than abstol so that we can enter the while loop
    N. <- 0 # N. will count the total number of function evaluations
    i. <- 0 # initialize counter; this will count the number of iterations in the while loop

    ## useskip <- 0 # will need that because the first iteration is a little different from all the others
    ## denom <- 1
    while(err > abstol && N. < fun.eval[2])
    {
        if(method == "sobol")
            .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )

        ## Get N RQCM estimates
        for(l in 1:B)
        {
            ## Get the point set
            ## If const = TRUE, we only need d - 1 (quasi) random numbers
            ## (the case const = TRUE and d = 1 has already been dealt with)
            U <- if(const){
                     U. <- switch(method,
                                  "sobol"   = {
                                      ## qrng::sobol(n = n., d = d - 1, randomize = TRUE, skip = (useskip * n.))
                                      qrng::sobol(n = n., d = d - 1, randomize = TRUE, skip = (i. * n.) )
                                  },
                                  "ghalton" = {
                                      qrng::ghalton(n = n., d = d - 1, method = "generalized")
                                  },
                                  "prng"    = {
                                      matrix(runif( n. * (d - 1)), ncol = d - 1)
                                  })

                     ## First and last column contain 1s corresponding to simulated values from sqrt(mix)
                     cbind( rep(1, n.), U., rep(1, n.))
                 } else {
                     U. <- switch(method,
                                  "sobol"   = {
                                      ## qrng::sobol(n = n., d = d, randomize = TRUE, skip = (useskip * n.))
                                      qrng::sobol(n = n., d = d, randomize = TRUE, skip = (i. * n.))
                                  },
                                  "ghalton" = {
                                      qrng::ghalton(n = n., d = d, method = "generalized")
                                  },
                                  "prng"    = {
                                      matrix(runif( n. * d), ncol = d)
                                  })

                     ## Case d = 1 somewhat special again:
                     if(d == 1){
                         cbind(sqrt(W(U.)), sqrt(W(1 - U.)))
                     } else {
                         ## Column 1:sqrt(mix), columns 2--d: unchanged (still uniforms),
                         ## column d + 1: antithetic realization of sqrt(mix)
                         cbind(sqrt(W(U.[, 1])), U.[, 2:d], sqrt(W(1 - U.[, 1])))
                     }
                 }

            ## Evaluate the integrand at the point set and save it in T[]; calls C.
            ## Both T.[l] and the new estimate are based on n. evaluations, so we
            ## can just average them unless we are in the first iteration in which
            ## case denom = 1 and T.[l] = 0.
            if(d == 1) {
                ## Case of dimension 1: Don't need to approximate the multivariate
                ##                      normal df and can just use pnorm()
                ## Note that d = 1 for a pure normal or t df has already been addressed

                ## TODO: once 'converged', remove the commented parts (or put in a switch and keep both)
                ## T.[l] <- (T.[l] + mean((pnorm(upper/U[,1])   - pnorm(lower/U[,1]) +
                ##                         pnorm(upper/U[,d+1]) - pnorm(lower/U[,d+1]))/2)) / denom
                T.[l] <- (i. * T.[l] + mean((pnorm(upper/U[,1])   - pnorm(lower/U[,1]) +
                                             pnorm(upper/U[,d+1]) - pnorm(lower/U[,d+1])) / 2)) / (i. + 1)
            } else {
                ## TODO: why keep this code? give a reason
                ## T.[l] <- (T.[l] + .Call("eval_nvmix_integral",
                ##                lower     = as.double(lower),
                ##                upper     = as.double(upper),
                ##                U         = as.double(U),
                ##                n         = as.integer(n.),
                ##                d         = as.integer(d),
                ##                cholScale = as.double(cholScale),
                ##                ZERO      = as.double(ZERO),
                ##                ONE       = as.double(ONE)) )/denom
                T.[l] <- (i. * T.[l] + .Call("eval_nvmix_integral",
                                             lower     = as.double(lower),
                                             upper     = as.double(upper),
                                             U         = as.double(U),
                                             n         = as.integer(n.),
                                             d         = as.integer(d),
                                             cholScale = as.double(cholScale),
                                             ZERO      = as.double(ZERO),
                                             ONE       = as.double(ONE)) ) / (i. + 1)
            }
        } # end for(l in 1:B)

        ## Update various variables
        N. <- N. + 2 * B * n. # number of function evaluations; (* 2 since antithetic variates are used in eval_nvmix_integral())
        ## TODO why keep this code? give a reason
        ## ## Change denom and useksip. This is done exactly once, namely in the first iteration.
        ## if(i. == 0){
        ##
        ##   denom <- 2
        ##   useskip <- 1
        ##
        ## } else {
        ##
        ##   ## Increase sample size n. This is done in all iterations except for the first two.
        ##   n. <- 2 * n.
        ##
        ## }
        sig <- sd(T.) # get standard deviation of the estimator
        err <- CI.factor * sig # update error; note that this CI.factor is actually CI.factor/sqrt(N) (see above)
        i. <- i. + 1 # update counter
    }

    ## Finalize
    T <- mean(T.) # calculate the RQMC estimator
    var <- (sig / sqrt(B))^2 # its variance
    if(err > abstol)
        warning("Precision level 'abstol' not reached; consider increasing the second component of 'fun.eval'")

    ## Return
    list(Prob = T, N = N., i = i., ErrEst = err, Var = var)
}
