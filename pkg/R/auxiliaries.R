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

##' @title Swap Variables i and j in a, b and R
##' @param i variable to be switched with j
##' @param j variable to be switched with i
##' @param a vector
##' @param b vector
##' @param Sigma matrix
##' @return list a, b, P after components/rows/columns i and j have been switched
##' @author Erik Hintz
swap <- function(i, j, a, b, Sigma)
{
    ## Build vector
    ij <- c(i,j)
    ji <- c(j,i)
    ## Reorder a, b
    a[ij] <- a[ji]
    b[ij] <- b[ji]
    ## Reorder Sigma
    wo.ij <- setdiff(seq_len(nrow(Sigma)), ij)
    temp_i <- as.matrix(Sigma[wo.ij,i])
    temp_j <- as.matrix(Sigma[wo.ij,j])
    temp_ii <- Sigma[i,i]
    Sigma[wo.ij,i] <- temp_j
    Sigma[wo.ij,j] <- temp_i
    Sigma[i,wo.ij] <- temp_j
    Sigma[j,wo.ij] <- temp_i
    Sigma[i,i] <- Sigma[j,j]
    Sigma[j,j] <- temp_ii
    ## Return
    list(a = a, b = b, Sigma = Sigma)
}
