\name{dStudent}
\alias{dStudent}
\title{Density of the Multivariate Student t Distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t}).
}
\usage{
dStudent(x, df, loc = rep(0, d), scale, 
       factor = tryCatch(factorize(scale), error = function(e) e), 
       log = FALSE)
}
\arguments{
  \item{x}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points.}
  \item{df}{degrees of freedom, positive. If \eqn{df=Inf}, the density of a Multivariate Normal Distribution is returned.}
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
  of the multivariate t distribution with
  location vector \code{loc}, scale matrix \code{scale}
  (a covariance matrix) and degrees of freedom \code{df}.
}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used is the Cholesky decomposition.
  To this end, \code{scale} needs to have full rank.
  
  This is a user-friendly wrapper for \code{dnvmix(x, mix = "inverse.gamma",...)}.
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
## Evaluate t_{3.5} density
df <- 3.5
x <- matrix(1:12/12, ncol = d) # evaluation points
dt <- dStudent(x, df = df, scale = P)
stopifnot(all.equal(dt, c(0.013266542, 0.011967156, 0.010760575, 0.009648682),
                    tol = 1e-7))
}
\keyword{distribution}