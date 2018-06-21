### Auxiliary tools ############################################################

##' @title Matrix Factorizations
##' @param x matrix
##' @param method factorization method (from 'specific' to more 'general')
##'        "chol": Cholesky decomposition (here: upper triangular matrix) for positive definite matrices
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
        ev <- eigen(x, symmetric = TRUE, ...)
        t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
    },
    "svd" = { # singular-value decomposition; in general not upper triangular
        sv <- svd(x, ...)
        t(sv$v %*% (t(sv$u) * sqrt(pmax(sv$d, 0))))
    },
    stop("Wrong 'method'"))
}

##' @title Swap Variables i and j in a, b and R
##' @param a vector
##' @param b vector
##' @param R covariance matrix
##' @param i which variables shall be switched?
##' @param j which variables shall be switched?
##' @return list a, b, R after variables i and j have been switched.
##' @author Erik Hintz
swap <- function(a, b, R, i, j)
{
    ## Reorder a, b
    a[c(i,j)] <- a[c(j,i)]
    b[c(i,j)] <- b[c(j,i)]
    ## Reorder R
    woij <- setdiff(seq_len(nrow(R)), c(i,j))
    temp_i <- as.matrix(R[woij,i])
    temp_j <- as.matrix(R[woij,j])
    temp_ii <- R[i,i]
    R[woij,i] <- temp_j
    R[woij,j] <- temp_i
    R[i,woij] <- temp_j
    R[j,woij] <- temp_i
    R[i,i] <- R[j,j]
    R[j,j] <- temp_ii
    ## Return
    list(a = a, b = b, R = R)
}
