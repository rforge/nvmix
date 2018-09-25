### dNorm() ####################################################################

##' @title Density of the Multivariate Normal Distribution
##' @param x (n, d)-matrix of evaluation points
##' @param loc location vector of dimension d (= mean vector here)
##' @param scale (d, d)-covariance matrix, positive definite
##' @param factor *upper triangular* factor of the covariance matrix 'scale' such
##'        that factor^T factor = scale here (otherwise det(scale) not computed
##'        correctly!)
##' @param log logical indicating whether the logarithmic density is computed
##' @return n-vector with N(loc, scale) density values
##' @author Marius Hofert and Erik Hintz
dmultinorm <- function(x, loc = rep(0, d), scale = diag(d),
                       factor = tryCatch(factorize(scale), error = function(e) e), # needs to be triangular!
                       log = FALSE)
{
    if(!is.matrix(x)) x <- rbind(x)
    d <- ncol(x)
    dnvmix(x = x, loc = loc, scale = scale, mix = "constant",
           factor = factor, log = log)$Density
}


