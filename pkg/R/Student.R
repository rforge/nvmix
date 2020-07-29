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
dStudentcopula <- function(u, df, scale = diag(d), factor = NULL, log = FALSE, 
                           verbose = TRUE)
{
   ## Checks 
   if(!is.matrix(u)) u <- rbind(u)
   d <- ncol(u) 
   n <- nrow(u)
   stopifnot(all(u <= 1), all(u >= 0)) 
   ## Result object
   res <- rep(-Inf, n)
   notNA <- rowSums(is.na(u)) == 0 
   not01 <- rowSums( u <= 0 | u >= 1 ) == 0 # rows where no component is <= 0 or >=1
   ## Fill in NAs where needed
   res[!notNA] <- NA
   ## Density is zero outside (0,1)^d 
   res[!not01 & notNA] <- if(log) -Inf else 0 
   u <- u[notNA & not01,, drop = FALSE] # non-missing data inside (0,1)^d
   ## Call 'dnvmixcopula()' with non-missing, (0,1)^d rows 
   res[notNA & not01] <- 
      dnvmixcopula(u, qmix = "inverse.gamma", scale = scale, factor = factor, 
                   verbose = verbose, df = df, log = log)
   res
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
                      factor = NULL, factor.inv = NULL, control = list(), 
                      log = FALSE, verbose = TRUE)
{
   if(!is.matrix(x)) x <- rbind(x)
   d <- ncol(x) # for 'loc', 'scale'
   ## Call 'dgnvmix()'
   dgnvmix(x, groupings = groupings, qmix = "inverse.gamma", loc = loc, 
           scale = scale, factor = factor, factor.inv = factor.inv, df = df, 
           control = control, log = log, verbose = verbose)
}

##' @title Density Function of the grouped t copula  
##' @param u (n, d) matrix of evaluation points 
##' @param groupings ?pgnvmix() 
##' @param df see ?pgStudent()
##' @param scale (d, d)-covariance matrix (scale != covariance matrix here)
##' @param control ?get_set_param() 
##' @param verbose logical indicating whether a warning is given if the required
##'        precision 'abstol' (see dnvmix()) has not been reached.
##' @param log logical if log-density required         
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
dgStudentcopula <- function(u, groupings = 1:d, df, scale = diag(d), factor = NULL,
                            factor.inv = NULL, control = list(), verbose = TRUE, 
                            log = FALSE)
{
   ## Checks 
   if(!is.matrix(u)) u <- rbind(u)
   d <- ncol(u) 
   n <- nrow(u)
   stopifnot(all(u <= 1), all(u >= 0)) 
   ## Result object
   res <- rep(-Inf, n)
   notNA <- rowSums(is.na(u)) == 0 
   not01 <- rowSums( u <= 0 | u >= 1 ) == 0 # rows where no component is <= 0 or >=1
   ## Fill in NAs where needed
   res[!notNA] <- NA
   ## Density is zero outside (0,1)^d 
   res[!not01 & notNA] <- if(log) -Inf else 0 
   u <- u[notNA & not01,, drop = FALSE] # non-missing data inside (0,1)^d
   ## Compute quantiles
   qu <- sapply(1:d, function(i) qt(u[, i], df = df[groupings[i]]))
   if(!is.matrix(qu)) qu <- rbind(qu) # otherwise dimension not correct in dgnvmix()
   num <- dgnvmix(qu, qmix = "inverse.gamma", scale = scale, factor = factor, 
                  factor.inv = factor.inv, df = df, groupings = groupings, 
                  verbose = verbose, 
                  control = control, log = TRUE) # vector 
   ## Matrix with marginal density applied on the columns of 'qu' 
   temp <- sapply(1:d, function(i) dt(qu[, i], df = df[groupings[i]], log = TRUE))
   if(!is.matrix(temp)) temp <- rbind(temp)
   denom <- rowSums(temp)
   ## Store results and return
   res[notNA & not01] <- if(!log) exp(num - denom) else num - denom
   ## Return 
   res
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
   if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
   pnvmix(upper, lower = lower, qmix = "inverse.gamma", loc = loc, scale = scale,
          standardized = standardized, control = control,
          verbose = verbose, df = df)
}

##' @title Distribution Function of the grouped Multivariate t Distribution
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
   if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
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
pStudentcopula <- function(upper, lower = matrix(0, nrow = n, ncol = d), df, 
                           scale = diag(d), control = list(), verbose = TRUE)
{
   ## Checks 
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
   ## Call more general pgStudentcopula() 
   pgStudentcopula(upper, lower = lower, groupings = rep(1, d), df = df, scale = scale,
                   control = control, verbose = verbose)
}

##' @title Distribution Function of the grouped t copula 
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
pgStudentcopula <- function(upper, lower = matrix(0, nrow = n, ncol = d),
                            groupings = 1:d, df, scale = diag(d), control = list(), 
                            verbose = TRUE)
{
   ## Checks 
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
   upper <- pmax( pmin(upper, 1), 0) 
   lower <- pmax( pmin(lower, 1), 0) 
   ## Transform limits via qt(..., df)
   upper_ <- sapply(1:d, function(i) qt(upper[, i], df = df[groupings[i]]))
   lower_ <- if(all(lower == 0)){
      matrix(-Inf, nrow = n, ncol = d) # avoid estimation of the quantile 
   } else {
      sapply(1:d, function(i) qt(lower[, i], df = df[groupings[i]]))
   }
   ## Call 'pgnvmix()' (which handles NA correctly) 
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
##' @param factor factor R of the covariance matrix 'scale' with d rows
##'        such that R R^T = 'scale'.
##' @return (n, d)-matrix with t_nu(loc, scale) samples
##' @author Erik Hintz and Marius Hofert
rgStudentcopula <- function(n, groupings = 1:d, df, scale = diag(2), factor = NULL, 
                            method = c("PRNG", "sobol", "ghalton"), skip = 0)
{
   method <- match.arg(method) 
   d <- if(!is.null(factor)) { 
      nrow(factor <- as.matrix(factor))
   } else {
      nrow(scale <- as.matrix(scale))
   }
   ## Sample from the grouped t distribution 
   t_sample <- 
      rgnvmix(n, qmix = "inverse.gamma", groupings = groupings, scale = scale,
              factor = factor, df = df, method = method, skip = skip)
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
rStudentcopula <- function(n, df, scale = diag(2), 
                           method = c("PRNG", "sobol", "ghalton"), skip = 0)
{
   d <- nrow(scale <- as.matrix(scale))
   method <- match.arg(method)
   ## Call more general 'rgStudentcop' without grouping 
   rgStudentcopula(n, groupings = rep(1, d), df = df, scale = scale, 
                   method = method, skip = skip)
}


##' @title Fitting the Parameters of a Multivariate Student t Distribution
##' @param x (n,d) data matrix
##' @param loc location vector; estimated if not supplied
##' @param scale (d,d) scale matrix; estimated if not supplied 
##' @param mix.param.bounds see ?fitnvmix
##' @param ... additional arguments passed to the underlying fitnvmix()
##' @return see ?fitnvmix
##' @author Marius Hofert
fitStudent <- function(x, loc = NULL, scale = NULL, mix.param.bounds = c(1e-3, 1e2), ...)
{
   fit <- fitnvmix(x, qmix = "inverse.gamma", 
                   loc = loc, scale = scale, mix.param.bounds = mix.param.bounds, ...)
   ## Consistency with other *Student() functions
   nms <- names(fit)
   nms[nms == "nu"] <- "df"
   names(fit) <- nms
   ## Return
   fit
}



#' Fitting grouped t-copulas
#' @param x (n, d) matrix of data the underlying copula of which is to be estimated
#' @param u (n, d) matrix of copula observations in (0,1) 
#' @param df.init NULL or vector with initial estimates for 'df'; can contain NAs 
#' @param scale NULL or known 'scale' matrix (estimated via p.w. Kendall's tau if not provided)
#' @param groupings see ?pgnvmix()
#' @param df.bounds 2-vector giving bounds on the dof parameter
#' @param control see ?get_set_param()
#' @param verbose logical if warnings shall be returned 
#' @return List of three:
#'         'df'      vector of estimated dof parameters
#'         'scale'   estimated scale matrix 
#'         'll'      estimated log-likelihood at optimum
#'         'df.init' vector of initial estimates for the dof parameters 
#' @author Erik Hintz  
#' @note Either 'x' or 'u' or both can be provided   
fitgStudentcopula <- function(x, u, df.init = NULL, scale = NULL, 
                           groupings = rep(1, d), df.bounds = c(0.5, 30),
                           control = list(), verbose = TRUE){
   
   ## 0 Setup ##################################################################
   ## Both 'x' and 'u' can be provided
   x.provided <- FALSE 
   if(hasArg(x)){
      if(!is.matrix(x)) x <- cbind(x) 
      x.provided <- TRUE
      notNA <- rowSums(is.na(x)) == 0
      x <- x[notNA,, drop = FALSE] # non-missing data (rows)
      n <- nrow(x) # sample size
      d <- ncol(x) # dimension 
      if(!hasArg(u)){ # pseudo-observations *not* provided 
         u <- copula::pobs(x) 
      } else { # pseudo-observations provided; remove NA and check dimension
         if(!is.matrix(u)) u <- cbind(u) 
         notNA <- rowSums(is.na(u)) == 0
         u <- u[notNA,, drop = FALSE] # non-missing data (rows)
         ## Check
         if(!all.equal(dim(u), c(n, d)))
            stop("Dimensions of 'u' and 'x' do not match.")
         if(any(u >= 1 | u <= 0))
            stop("Elements in 'u' must be in (0,1).")
      }
   } else {
      if(!hasArg(u))
         stop("Either 'u' or 'x' or both must be provided.")
      if(!is.matrix(u)) u <- rbind(u) 
      notNA <- rowSums(is.na(u)) == 0
      u <- u[notNA,, drop = FALSE] # non-missing data (rows)
      n <- nrow(u) # sample size
      d <- ncol(u) # dimension 
   }
   ## At least two data points must be provided
   if(n <= 1)
      stop("Data-set must have at least two rows.")
   
   ## Initialize various quantities
   control <- get_set_param(control)
   numgroups <- length(unique(groupings)) # number of groups 
   stopifnot(all(groupings %in% 1:numgroups)) 
   
   ## 1 Estimation of 'scale' ##################################################   
   
   if(is.null(scale)){ # 'scale' not provided => estimate it
      scale <- sin(pcaPP::cor.fk(u) * pi/2)
      ## Ensure positive-definitness
      scale <- as.matrix(Matrix::nearPD(scale)$mat)
   } else stopifnot(all.equal(dim(scale), c(d, d)))
   factor.inv <- solve(t(chol(scale))) # for repeated calls of 'dgStudentcopula()' => faster 
   
   ## 2 Estimation of 'df' #####################################################   
   
   ## 2.1 Find starting values #################################################  
   
   if(!is.null(df.init)){
      stopifnot(length(df.init) == numgroups)
      initNA <- which(is.na(df.init))
      if(length(initNA) < numgroups) 
         stopifnot(all(df.init[!initNA] >= df.bounds[1]), all(df.init[!initNA] <= df.bounds[2]))
   } else {
      df.init <- rep(NA, numgroups)
      initNA <- 1:numgroups
   }
   
   if(length(initNA) > 0){
      ## -log-likelihood (for faster evaluation)
      nLLt <- function(nu, P, u) {
         x <- qt(u, df = nu)
         -sum(dStudent(x, sigma = P, df = nu, log = TRUE) - rowSums(dt(x, df = nu, log = TRUE)))
      }
      ## Estimate dof in each group where 'df.init' was not provided 
      for(k in initNA){
         ind.sub <- which(groupings == k)
         d.sub <- length(ind.sub) # dimension of the group 
         if(d.sub > 1){
            ## Group at least bivariate => Estimate 'df' of t-copula 
            df.init[k] <- optimize(nLLt, interval = df.bounds, u = u[, ind.sub],
                                   P = scale[ind.sub, ind.sub])$minimum
         } else {
            ## Group 'univariate', so margin is merely uniform 
            if(x.provided){
               ## If 'x' provided, assume x[, ind.sub] ~ t_{nu} => estimate 'nu' as MLE
               df.init[k] <- fitStudent(cbind(x[, ind.sub]))$df
            } else {
               df.init[k] <- 5
            }
         }
      }
   }
   
   ## 2.2 Joint estimation of 'df' #############################################   
   
   if(numgroups == 1){
      ## One group => classical t copula => 'df.init' is MLE
      df <- df.init 
   } else {
      ## -loglikelihood as a function of 'df' 
      seed <- sample(1:1e3, 1) # for reproducibility 
      nLLgt <- function(df){
         set.seed(seed) # => monotonicity 
         if(any(df < df.bounds[1]) | any(df > df.bounds[2])) return(Inf)
         -sum(dgStudentcopula(u, groupings = groupings, df = df, control = control,
                              factor.inv = factor.inv, log = TRUE))
      }
      ll.init <- -nLLgt(df.init) # likelihood of initial parameter 'df.init' 
      ## Call 'optim' and grab 'df' along with likelihood 
      opt.obj <- optim(df.init, nLLgt, control = control$control.optim)
      df <- opt.obj$par
      ## Check 'convergence' returned by optim()
      if(verbose){
         if(opt.obj$convergence == 1)
            warning("Maximum number of iterations exhaused in optim(); consider increasing 'optim.maxit' in the control argument.")
         if(opt.obj$convergence == 10)
            warning("optim() detected degeneracy of the Nelder-Mead simplex.")
      }
      ll.mle <- -opt.obj$value 
      ## Check if likelihood increased
      if(verbose & (ll.mle < ll.init) )
         warning("'df.init' yields larger likelihood than 'df' returned from 'optim()'.")
   }
   
   ## 3. Return ################################################################ 
   
   list(df = df, scale = scale, ll = ll.mle, df.init = df.init)
   
}