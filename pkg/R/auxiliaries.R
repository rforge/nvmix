### Auxiliary tools ############################################################

##' @title Matrix Factorizations
##' @param x matrix
##' @param method factorization method (from 'specific' to more 'general')
##'        "chol": Cholesky decomposition (here: upper triangular matrix) for
##'                positive definite matrices
##'        "chol.pivot": Cholesky decomposition with pivoting (see rmvnorm())
##'        "eigen": eigendecomposition (rmvt() -> rmvnorm() => default for rmvt())
##'        "svd": singular-value decomposition (see rmvnorm())
##' @param ... additional arguments passed to the underlying functions
##' @return factor (factorized matrix)
##' @author Marius Hofert
##' @note Could be called with tryCatch(factorize(x), error = function(e) e)
##'       but probably not necessary
factorize <- function(x, method = c("chol", "chol.pivot", "eigen", "svd"),
                      ...)
{
    method <- match.arg(method)
    switch(method,
           "chol" = { # for positive definite matrices; typically fastest
               chol(x, ...)
           },
           "chol.pivot" = { # for positive semidefinite matrices
               ## ... but can result in non-upper triangular factor, see:
               ## set.seed(271)
               ## A <- matrix(runif(d * d), ncol = d)
               ## P <- cov2cor(A %*% t(A))
               ## chol(P) # upper triangular
               ## (R <- chol(P, pivot = TRUE)) # also upper triangular
               ## R[, order(attr(R, "pivot"))] # not upper triangular anymore
               R <- chol(x, pivot = TRUE, ...)
               R[, order(attr(R, "pivot"))] # t(L) for L the Cholesky factor; upper triangular
           },
           "eigen" = { # eigendecomposition; in general not upper triangular
               ev <- eigen(x, ...) # uses 'isSymmetric()' to determine whether symmetric
               ## Note: eigen() *normalizes* the eigenvectors to unit length
               ##       (only then Q^{-1} = Q^T and thus AA\T = Sigma for A = Q\Lambda^{1/2})
               t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
           },
           "svd" = { # singular-value decomposition; in general not upper triangular
               sv <- svd(x, ...)
               sv$u %*% sqrt(diag(pmax(sv$d, 0))) %*% t(sv$v) # or t(sv$v %*% (t(sv$u) * sqrt(pmax(sv$d, 0))))
           },
           stop("Wrong 'method'"))
}

