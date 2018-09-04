\name{dmultinorm}
\alias{dmultinorm}
\title{Density of the Multivariate Normal Distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t}).
}
\usage{
dmultinorm(x, loc = rep(0, d), scale, 
       factor = tryCatch(factorize(scale), error = function(e) e),
       log = FALSE)
}
\arguments{
  \item{x}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points.}
  \item{loc}{location vector of dimension \eqn{d} (the number of columns
    of \code{x}).}
  \item{scale}{covariance matrix of dimension \eqn{(d, d)}.}
  \item{factor}{factorization matrix of the covariance matrix
    \code{scale}. A \code{\link{matrix}} \eqn{R} with \eqn{d} rows
    such that \eqn{R^T R} equals \code{scale}. Although not tested
    (for efficiency reasons), \eqn{R} (thus \code{factor}) has
    to be \strong{upper triangular} (otherwise the determinant
        of \code{scale} involved in the density is not computed correctly).
    Also note that \code{factor} can be \code{NULL} or inherit from
    \code{"error"} in which case the degenerate density is computed
    as \code{Inf} at \code{loc} and 0 everywhere else.}
  \item{log}{\code{\link{logical}} indicating whether the logarithm
    of the density is to be computed.}
}
\value{
  \code{dnvmix()} returns an \eqn{n}-\code{\link{vector}} with the
  density values (default) or log-density values (if \code{log})
  of the mutlivariate normal distribution with
  location vector \code{loc} and scale matrix \code{scale}
  (a covariance matrix).
}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used is the Cholesky decomposition.
  To this end, \code{scale} needs to have full rank.
  
  This is a user-friendly wrapper for \code{dnvmix(x, mix = "constant",...)}.
}
\author{Marius Hofert and Erik Hintz}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{pnvmix}()}, \code{\link{rnvmix}()}, \code{\link{dnvmix}()}
}
\examples{
## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
x <- matrix(1:12/12, ncol = d) # evaluation points

## Evaluate normal density
dn <- dmultinorm(x, scale = P)
stopifnot(all.equal(dn, c(0.013083858, 0.011141923, 0.009389987, 0.007831596),
                    tol = 1e-7))
}
\keyword{distribution}