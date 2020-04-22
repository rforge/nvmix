### d/p/r/fitStudent() #########################################################

##' @title Density of the Multivariate Student t Distribution
##' @param x (n, d)-matrix of evaluation points
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix, positive definite (scale != covariance
##'        matrix here)
##' @param factor *lower triangular* factor R of the covariance matrix 'scale'
##'        such that R^T R = 'scale' here (otherwise det(scale) not computed
##'        correctly!)
##' @param log logical indicating whether the logarithmic density is computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @param ... additional arguments passed to the underlying dnvmix()
##' @return n-vector of t_nu(loc, scale) density values
##' @author Erik Hintz and Marius Hofert
dStudent <- function(x, df, loc = rep(0, d), scale = diag(d),
                     factor = NULL, # needs to be triangular!
                     log = FALSE, verbose = TRUE, ...)
{
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x) # for 'loc', 'scale'
    dnvmix(x, qmix = "inverse.gamma", loc = loc, scale = scale,
           factor = factor, log = log, verbose = verbose, df = df, ...)
}

##' @title Density of the t copula 
##' @param u (n, d)-matrix of evaluation points
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param scale (d, d)-covariance matrix, positive definite (scale != covariance
##'        matrix here)
##' @param log logical indicating whether the logarithmic density is computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @param ... additional arguments passed to the underlying dnvmix()
##' @return n-vector of t_nu(loc, scale) density values
##' @author Erik Hintz and Marius Hofert
dStudentcop <- function(u, df, scale = diag(d), log = FALSE, verbose = TRUE)
{
   if(!is.matrix(u)) u <- rbind(u)
   d <- ncol(u) # for 'loc', 'scale'
   ## Call 'dnvmixcop'
   dnvmixcop(u, qmix = "inverse.gamma", scale = scale, log = log, verbose = verbose,
             df = df)
}

##' @title Density of the grouped t distribution
##' @param x (n, d)-matrix of evaluation points
##' @param groupings see ?pgnvmix() 
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix, positive definite (scale != covariance
##'        matrix here)
##' @param control list; see ?get_set_param() 
##' @param log logical indicating whether the logarithmic density is computed
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @return n-vector of t_nu(loc, scale) density values
##' @author Erik Hintz and Marius Hofert
dgStudent <- function(x, groupings = 1:d, df, loc = rep(0, d), scale = diag(d), 
                      control = list(), log = FALSE, verbose = TRUE)
{
   if(!is.matrix(x)) x <- rbind(x)
   d <- ncol(x) # for 'loc', 'scale'
   ## Call 'dgnvmix()'
   dgnvmix(x, groupings = groupings, qmix = "inverse.gamma", loc = loc, 
           scale = scale, df = df, factor = NULL, control = control, log = log,
           verbose = verbose)
}


##' @title Distribution Function of the Multivariate Student t Distribution
##' @param upper d-vector of upper evaluation limits
##' @param lower d-vector of lower evaluation limits
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param control ?get_set_param() 
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pStudent <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d),
                     df, loc = rep(0, d), scale = diag(d), standardized = FALSE,
                     control = list(), verbose = TRUE)
{
    ## Checks (needed to get the default for 'lower' correctly)
    if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
    n <- nrow(upper) # number of evaluation points
    d <- ncol(upper) # dimension
    pnvmix(upper, lower = lower, qmix = "inverse.gamma", loc = loc, scale = scale,
           standardized = standardized, control = control,
           verbose = verbose, df = df)
}

##' @title Distribution Function of the generalized/grouped Multivariate t Distribution
##' @param upper d-vector of upper evaluation limits
##' @param lower d-vector of lower evaluation limits
##' @param groupings ?pgnvmix() 
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param control ?get_set_param() 
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pgStudent <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d),
                     groupings = 1:d, 
                     df, loc = rep(0, d), scale = diag(d), standardized = FALSE,
                     control = list(), verbose = TRUE)
{
   ## Checks (needed to get the default for 'lower' correctly)
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   ## Call 'pgnvmix()' 
   pgnvmix(upper, lower = lower, groupings = groupings, qmix = "inverse.gamma", 
           loc = loc, scale = scale, standardized = standardized, control = control,
          verbose = verbose, df = df)
}

##' @title Distribution Function of the t copula
##' @param upper d-vector of upper evaluation points in [0,1]^d
##' @param lower d-vector of lower evaluation limits in [0,1]^d
##' @param df dof parameter 
##' @param scale (d, d)-covariance matrix 
##' @param control ?get_set_param() 
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pStudentcop <- function(upper, lower = matrix(0, nrow = n, ncol = d), df, 
                        scale = diag(d), control = list(), verbose = TRUE)
{
   ## Checks 
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   ## Call more general pgStudentcop() 
   pgStudentcop(upper, lower = lower, groupings = rep(1, d), df = df, scale = scale,
                control = control, verbose = verbose)
}

##' @title Distribution Function of the generalized/grouped Multivariate t Distribution
##' @param upper d-vector of upper evaluation points in [0,1]^d
##' @param lower d-vector of lower evaluation limits in [0,1]^d
##' @param groupings ?pgnvmix() 
##' @param df see ?pgStudent()
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param control ?get_set_param() 
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pgStudentcop <- function(upper, lower = matrix(0, nrow = n, ncol = d),
                      groupings = 1:d, df, scale = diag(d), control = list(), 
                      verbose = TRUE)
{
   ## Checks 
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   ## Transform limits via qt(..., df)
   upper_ <- sapply(1:d, function(i) qt(upper[, i], df = df[groupings[i]]))
   lower_ <- sapply(1:d, function(i) qt(lower[, i], df = df[groupings[i]]))
   ## Call 'pgnvmix()' 
   pgnvmix(upper_, lower = lower_, groupings = groupings, qmix = "inverse.gamma", 
           scale = scale, control = control, verbose = verbose, df = df)
}

##' @title Random Number Generator for the Multivariate Student t Distribution
##' @param n sample size
##' @param df degrees of freedom > 0; if df = Inf, sample from a Normal distribution
##'        is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param factor factor R of the covariance matrix 'scale' with d rows
##'        such that R R^T = 'scale'.
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rStudent <- function(n, df, loc = rep(0, d), scale = diag(2),
                     factor = NULL, method = c("PRNG", "sobol", "ghalton"), 
                     skip = 0)
{
   method <- match.arg(method) 
   d <- if(!is.null(factor)) { # for 'loc', 'scale'
             nrow(factor <- as.matrix(factor))
         } else {
             nrow(scale <- as.matrix(scale))
         }
    
    if(method == "PRNG"){
       ## Provide 'rmix' and no 'qmix' => typically faster
       rnvmix(n, rmix = "inverse.gamma", 
              loc = loc, scale = scale, factor = factor, df = df,
              method = method, skip = skip)
    } else {
       ## Provide 'qmix' for inversion based methods
       rnvmix(n, qmix = "inverse.gamma", 
              loc = loc, scale = scale, factor = factor, df = df,
              method = method, skip = skip)
    }
}

##' @title Random Number Generator for the generalzied Multivariate t Distribution
##' @param n sample size
##' @param groupings see ?pgnvmix() 
##' @param df degrees of freedom > 0; if df = Inf, sample from a Normal distribution
##'        is returned
##' @param loc d-vector (location != mean vector here)
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param factor factor R of the covariance matrix 'scale' with d rows
##'        such that R R^T = 'scale'.
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rgStudent <- function(n, groupings = 1:d, df, loc = rep(0, d), scale = diag(2),
                      factor = NULL, method = c("PRNG", "sobol", "ghalton"), 
                      skip = 0)
{
   method <- match.arg(method) 
   d <- if(!is.null(factor)) { # for 'loc', 'scale'
      nrow(factor <- as.matrix(factor))
   } else {
      nrow(scale <- as.matrix(scale))
   }
   ## Provide 'qmix' (=> 'rmix' not allowed for gNVM() distributions)
   rgnvmix(n, qmix = "inverse.gamma", groupings = groupings,
           loc = loc, scale = scale, factor = factor, df = df,
           method = method, skip = skip)
}

##' @title Random Number Generator for grouped t copla 
##' @param n sample size
##' @param groupings see ?pgnvmix() 
##' @param df degrees of freedom > 0; if df = Inf, sample from a Normal distribution
##'        is returned
##' @param scale (d, d)- correlation matrix
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rgStudentcop <- function(n, groupings = 1:d, df, scale = diag(2),
                         method = c("PRNG", "sobol", "ghalton"), skip = 0)
{
   method <- match.arg(method) 
   d <- nrow(scale <- as.matrix(scale))
   ## Sample from the grouped t distribution 
   t_sample <- rgnvmix(n, qmix = "inverse.gamma", groupings = groupings, 
                       scale = scale, df = df, method = method, skip = skip)
   ## Apply the correct pt(, df) columnwise and return
   sapply(1:d, function(i) pt(t_sample[, i], df = df[groupings[i]]))
}

##' @title Random Number Generator for the t copla 
##' @param n sample size
##' @param df degrees of freedom > 0; if df = Inf, sample from a Normal distribution
##'        is returned
##' @param scale (d, d)- correlation matrix
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rStudentcop <- function(n, df, scale = diag(2), 
                        method = c("PRNG", "sobol", "ghalton"), skip = 0)
{
   d <- nrow(scale <- as.matrix(scale))
   method <- match.arg(method)
   ## Call more general 'rgStudentcop' without grouping 
   rgStudentcop(n, groupings = rep(1, d), df = df, scale = scale, 
                method = method, skip = skip)
}


##' @title Fitting the Parameters of a Multivariate Student t Distribution
##' @param x (n,d) data matrix
##' @param mix.param.bounds see ?fitnvmix
##' @param ... additional arguments passed to the underlying fitnvmix()
##' @return see ?fitnvmix
##' @author Marius Hofert
fitStudent <- function(x, mix.param.bounds = c(1e-3, 1e2), ...)
{
    fit <- fitnvmix(x, qmix = "inverse.gamma", mix.param.bounds = mix.param.bounds, ...)
    ## Consistency with other *Student() functions
    nms <- names(fit)
    nms[nms == "nu"] <- "df"
    names(fit) <- nms
    ## Return
    fit
}
