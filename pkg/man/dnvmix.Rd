\name{dnvmix}
\alias{dnvmix}
\title{Density of Multivariate Normal Variance Mixtures}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t}).
}
\usage{
dnvmix(x, df, loc = rep(0, d), scale,
       factor = tryCatch(factorize(scale), error = function(e) e),
       log = FALSE)
}
\arguments{
  \item{x}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points.}
  \item{df}{degrees of freedom (positive real number or \code{Inf} in
    which case the density of a multivariate normal distribution with
    mean vector \code{loc} and covariance matrix \code{scale} is
    evaluated).}
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
  of the specified multivariate normal variance mixture with
  location vector \code{loc} and scale matrix \code{scale}
  (a covariance matrix).
}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used is the Cholesky decomposition.
  To this end, \code{scale} needs to have full rank.
}
\author{Marius Hofert}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{pnvmix}()}, \code{\link{rnvmix}()}
}
\examples{
## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Evaluate t_{3.5} density
df <- 3.5
x <- matrix(1:12/12, ncol = d) # evaluation points
dt <- dnvmix(x, df = df, scale = P)
stopifnot(all.equal(dt, c(0.013266542, 0.011967156, 0.010760575, 0.009648682),
                    tol = 1e-7))
## Evaluate normal density
dn <- dnvmix(x, df = Inf, scale = P)
stopifnot(all.equal(dn, c(0.013083858, 0.011141923, 0.009389987, 0.007831596),
                    tol = 1e-7))

## Missing data
x[3,2] <- NA
x[4,3] <- NA
dt <- dnvmix(x, df = df, scale = P)
stopifnot(is.na(dt) == rep(c(FALSE, TRUE), each = 2))

## Univariate case
x <- matrix(1:10/10, ncol = 1)
dt <- dnvmix(x, df = df, factor = 1)
dt. <- dt(as.vector(x), df = df)
stopifnot(all.equal(dt, dt.))
}
\keyword{distribution}