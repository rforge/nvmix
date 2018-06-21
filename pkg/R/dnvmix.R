### dnvmix() ###################################################################

##' @title Density of Multivariate Normal Variance Mixtures
##' @param x (n, d)-matrix of evaluation points
##' @param df degrees of freedom (positive real or Inf in which case the density
##'        of a N(loc, scale) is evaluated)
##' @param loc location vector of dimension d
##' @param scale covariance matrix of dimension (d, d)
##' @param factor factorization matrix of the covariance matrix scale;
##'        caution: this has to be an *upper triangular* matrix R
##'        such that R^T R = scale here (otherwise det(scale) not computed correctly)
##' @return n-vector with t_nu(loc, scale) density values
##' @author Marius Hofert
dnvmix <- function(x, df, loc = rep(0, d), scale,
                   factor = tryCatch(factorize(scale), error = function(e) e), # needs to be triangular!
                   log = FALSE)
{
    if(!is.matrix(x)) x <- rbind(x)
    n <- nrow(x)
    d <- ncol(x)
    stopifnot(df > 0, length(loc) == d)
    notNA <- apply(!is.na(x), 1, all)
    lres <- rep(-Inf, n)
    lres[!notNA] <- NA
    x <- x[notNA,] # available points
    tx <- t(x) # (d, n)-matrix
    if(inherits(factor, "error") || is.null(factor)) {
        lres[notNA & (colSums(tx == loc) == d)] <- Inf
    } else {
        ## Solve R^T * z = x - mu for z, so z = (R^T)^{-1} * (x - mu) (a (d, d)-matrix)
        ## => z^2 (=> componentwise) = z^T z = (x - mu)^T * ((R^T)^{-1})^T (R^T)^{-1} (x - mu)
        ##                           = z^T z = (x - mu)^T * R^{-1} (R^T)^{-1} (x - mu)
        ##                           = (x - mu)^T * (R^T R)^{-1} * (x - mu)
        ##                           = (x - mu)^T * scale^{-1} * (x - mu) = quadratic form
        z <- backsolve(factor, tx - loc, transpose = TRUE)
        maha2 <- colSums(z^2) # = sum(z^T z); squared Mahalanobis distance from x to mu w.r.t. scale
        ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
        ## = log(prod(diag(R))) = sum(log(diag(R)))
        lrdet <- sum(log(diag(factor)))
        lres[notNA] <- if(is.finite(df)) {
                           df.d.2 <- (df + d) / 2
                           lgamma(df.d.2) - lgamma(df/2) - (d/2) * log(df * pi) - lrdet - df.d.2 * log1p(maha2 / df)
                       } else {
                           - (d/2) * log(2 * pi) - lrdet - maha2/2
                       }
    }
    if(log) lres else exp(lres) # also works with NA, -Inf, Inf
}
