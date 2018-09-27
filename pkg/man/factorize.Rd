\name{factorize}
\alias{factorize}
\title{Matrix Factorization}
\description{
  Factorizing matrices with various method.
}
\usage{
factorize(x, method = c("chol", "chol.pivot", "eigen", "svd"), ...)
}
\arguments{
  \item{x}{\code{\link{matrix}}.}
  \item{method}{\code{\link{character}} string specifying the
    factorization method. Available are:
    \describe{
      \item{\code{"chol"}}{Cholesky decomposition of a positive definite
	matrix based on \code{\link{chol}()}.}
      \item{\code{"chol.pivot"}}{Cholesky decomposition of a positive
	semidefinite matrix based on \code{\link{chol}(, pivot = TRUE)};
        see \code{?chol} for more details.}
      \item{\code{"eigen"}}{factor from an eigendecomposition.}
      \item{\code{"svd"}}{factor from a singular-value decomposition.}
    }
  }
  \item{\dots}{additional arguments passed to the underlying main
    functions; see the source code for details.}
}
\value{
  The factor, that is, the factorized \code{\link{matrix}}. This is a
  matrix \eqn{R} such that such that \eqn{R^T R} equals \code{x}.
}
\details{
  When the factor \eqn{R} is multiplied from the right to an
  \eqn{(n, d)}-matrix of independent standard normals, one obtains
  a sample from a multivariate normal with zero mean vector and
  covariance matrix \code{x}.
}
\author{Marius Hofert}
\seealso{
  \code{\link{rnvmix}()} where \code{factorize()} is used.
}
\examples{
## Generate a random correlation matrix in dimension 5
d <- 5
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Factorizations
A.chol       <- factorize(P) # Cholesky factor
A.chol.pivot <- factorize(P, method = "chol.pivot") # with pivoting
A.eigen      <- factorize(P, method = "eigen") # factor from eigendecomposition
A.svd        <- factorize(P, method = "svd") # from singular-value decomposition
## Check
P.chol       <- t(A.chol) \%*\% A.chol
P.chol.pivot <- t(A.chol.pivot) \%*\% A.chol.pivot
P.eigen      <- t(A.eigen) \%*\% A.eigen
P.svd        <- t(A.svd) \%*\% A.svd
stopifnot(all.equal(P.chol,  P), all.equal(P.chol.pivot, P),
          all.equal(P.eigen, P), all.equal(P.svd,        P))
}
\keyword{array}