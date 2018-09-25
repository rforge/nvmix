### dStudent() #################################################################

##' @title Density of the Multivariate t Distribution
##' @param x (n, d)-matrix of evaluation points
##' @param df degrees of freedom > 0; if df = Inf, the normal density is returned
##' @param loc location vector of dimension d
##' @param scale (d, d)-covariance matrix, positive definite
##' @param factor *upper triangular* factor of the covariance matrix 'scale' such
##'        that factor^T factor = scale here (otherwise det(scale) not computed
##'        correctly!)
##' @param log logical indicating whether the logarithmic density is computed
##' @return n-vector with t_nu(loc, scale) density values
##' @author Marius Hofert and Erik Hintz
dStudent <- function(x, df, loc = rep(0, d), scale = diag(d),
                     factor = tryCatch(factorize(scale), error = function(e) e), # needs to be triangular!
                     log = FALSE)
{
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x)
    dnvmix(x = x, loc = loc, scale = scale, mix = "inverse.gamma",
           factor = factor, log = log, df = df)$Density
}

