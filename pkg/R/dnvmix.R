### dnvmix() ###################################################################

##' @title Density of a Multivariate Normal Variance Mixture
##' @param x (n, d)-matrix of evaluation points
##' @param mix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param loc d-vector (location vector)
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param factor *upper triangular* factor R of the covariance matrix 'scale'
##'        such that R^T R = 'scale' here (otherwise det(scale) not computed
##'        correctly!)
##' @param method character string indicating the method to be used:
##'         - "sobol":   Sobol sequence
##'         - "ghalton": generalized Halton sequence
##'         - "prng":    pure Monte Carlo
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
##' @param log logical indicating whether the logarithmic density is to be computed
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return list with the
##'         - computed probability
##'         - total number of function evaluations
##'         - number of iterations in the while loop
##'         - error estimate
##'         - variance estimate
##' @author Erik Hintz and Marius Hofert
dnvmix <- function(x, mix, loc = rep(0, d), scale = diag(d), # TODO: do we need a 'standardized = FALSE' here? (see pnvmix())
                   factor = factorize(scale), # needs to be triangular!
                   method = c("sobol", "ghalton", "PRNG"),
                   abstol = 0.001, CI.factor = 3.3, fun.eval = c(2^6, 1e8), B = 12,
                   log = FALSE, ...)
{
    ## Checks
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x) # dimension
    if(!is.matrix(scale)) scale <- as.matrix(scale)
    stopifnot(length(loc) == d, dim(scale) == c(d, d), # note: 'mix' is tested later
              abstol >= 0, CI.factor >= 0, length(fun.eval) == 2, fun.eval >= 0, B >= 1,
              is.logical(log))
    method <- match.arg(method)

    ## Define the quantile function of the mixing variable
    ## If 'mix' is "constant" or "inverse.gamma", we use the analytical formulas
    const <- FALSE # logical indicating whether we have a multivariate normal
    inv.gam <- FALSE # logical indicating whether we have a multivariate t
    qW <- if(is.character(mix)) { # 'mix' is a character vector specifying supported mixture distributions (utilizing '...')
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
              function(u)
                  return(do.call(qmix, c(u, mix[-1])))
          } else if(is.function(mix)) { # 'mix' is interpreted as the quantile function F_W^- of the mixture distribution F_W of W
              function(u)
                  return(mix(u, ...))
          } else stop("'mix' must be a character string, list or quantile function.")

    ## Build result object (log-density)
    n <- nrow(x)
    lres <- rep(-Inf, n) # n-vector of results
    notNA <- rowSums(is.na(x)) == 0
    lres[!notNA] <- NA
    x <- x[notNA,] # non-missing data (rows)

    ## Actual computation
    if(inherits(factor, "error") || is.null(factor)) { # if 'factor' has a problem
        lres[notNA & (rowSums(x == loc) == d)] <- Inf # if a row of x is equal to the location vector; TODO: why??? density there should be finite, right?
    } else {
        ## Solve R^T * z = x - mu for z, so z = (R^T)^{-1} * (x - mu) (a (d, d)-matrix)
        ## => z^2 (=> componentwise) = z^T z = (x - mu)^T * ((R^T)^{-1})^T (R^T)^{-1} (x - mu)
        ##                           = z^T z = (x - mu)^T * R^{-1} (R^T)^{-1} (x - mu)
        ##                           = (x - mu)^T * (R^T R)^{-1} * (x - mu)
        ##                           = (x - mu)^T * scale^{-1} * (x - mu) = quadratic form
        z <- backsolve(factor, t(x) - loc, transpose = TRUE) # TODO: can't we work with transpose = FALSE and avoid transposing 'x'?
        maha2 <- colSums(z^2) # = sum(z^T z); n-vector of squared Mahalanobis distances from x to mu w.r.t. scale
        ## TODO: can't we do this with mahalanobis()?
        ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
        ##                       = log(prod(diag(R))) = sum(log(diag(R)))
        lrdet <- sum(log(diag(factor)))

        ## Counters TODO: better names (as in pnvmix())
        N. <- 0 # N. will count the total number of function evaluations
        i. <- 0 # initialize counter; this will count the number of iterations in the while loop

        ## Deal with the different distributions
        if(inv.gam) { # multivariate t
            df.d.2 <- (df + d) / 2
            lres[notNA] <- lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) - lrdet - df.d.2 * log1p(maha2 / df)
            err <- 0
            var <- 0
        } else if(const) { # multivariate normal
            lres[notNA] <- -(d/2) * log(2 * pi) - lrdet - maha2/2
            err <- 0
            var <- 0
        } else { # general case of a multivariate normal variance mixture (RQMC)

            ## Basics
            CI.factor <- CI.factor / sqrt(B) # instead of dividing sigma by sqrt(B) each time
            n. <- fun.eval[1] # initial n
            T. <- matrix(0, ncol = n, nrow = B) # matrix to store RQMC estimates
            err <- abstol + 42 # initialize err to something bigger than abstol so that we can enter the while loop
            useskip <- 0 # will need that because the first iteration is a little different from all the others
            denom <- 1

            ## Make sure seed exists for 'method' being "sobol"
            if(method == "sobol") {
                if(!exists(".Random.seed")) runif(1) # dummy to generate .Random.seed
                seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used
            }

            ## Main loop
            while(err > abstol && N. < fun.eval[2]) {
                if(method == "sobol") .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )

                ## Get B RQCM estimates
                for(l in 1:B) {
                    ## Get the point set
                    U <- switch(method,
                                "sobol"   = {
                                    qrng::sobol(n., d = 1, randomize = TRUE, skip = (useskip * n.))
                                },
                                "gHalton" = {
                                    qrng::ghalton(n., d = 1, method = "generalized")
                                },
                                "prng"    = {
                                    cbind(runif(n.)) # 1-column matrix
                                })
                    W <- qW(U) # n.-vector of W's

                    ## exp-log trick
                    b <- - (d/2) * log(2 * pi * W) - lrdet - outer(1/W, maha2 / 2) # (n., n)-matrix, each column corresponds to "one x"
                    bmax <- apply(b, 2, max) # n-vector of maximal b's
                    T.[l,] <- (T.[l,] - log(n.) + bmax +
                               log(colSums(exp(b - rep(bmax, each = n.))))) / denom
                }

                ## Update various variables
                N. <- N. + B * n.

                ## Change denom and useskip; this is done exactly once, namely in the first iteration.
                if(i. == 0){
                    denom <- 2
                    useskip <- 1

                } else {
                    ## Increase the sample size n; this is done in all iterations except the first
                    n. <- 2 * n.
                }

                ## Compute error measures and update counter
                sig <- max(apply(T., 2, sd)) # get standard deviation of the column with the largest standard deviation
                err <- CI.factor * sig # update error. Note that this CI.factor is actually CI.factor/sqrt(N)
                i. <- i. + 1 # update counter
            } # while()

            ## Finalize
            var <- (sig / sqrt(B))^2
            if(err > abstol)
                warning("Precision level 'abstol' not reached; consider increasing the second component of 'fun.eval'")
            lres[notNA] <- colMeans(T.)
        }
    }

    ## Return
    if(log) {
        list(Density = lres,      N = N., i = i., ErrEst = err, Var = var)
    } else {
        list(Density = exp(lres), N = N., i = i., ErrEst = err, Var = var)
    }
}
