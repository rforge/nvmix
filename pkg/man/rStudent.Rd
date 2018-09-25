\name{rStudent}
\alias{rStudent}
\title{Random Number Generator for the Multivariate Student t Distribution}
\description{
  Generate vectors of random variates from the multivariate Student \emph{t} distribution.
}
\usage{
rStudent(n, loc = rep(0, d), scale, df, factor = factorize(scale))
}
\arguments{
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{loc}{location vector of dimension \eqn{d}.}
  \item{scale}{covariance matrix of dimension \eqn{(d, d)}.}
  \item{df}{degress of freedom (positive real or Inf in which case a sample from N(loc, scale) is returned.}
  \item{factor}{factorization matrix of the covariance matrix
    \code{scale}. A \code{\link{matrix}} \eqn{R} with \eqn{d}
    rows such that \eqn{R^T R} equals \code{scale}. \eqn{R} is internally
    multiplied from the right to an \eqn{(n, d)}-matrix of independent
    standard normals to obtain a sample from a multivariate normal
    with zero mean vector and covariance matrix \code{scale}.}
}
\value{
  \code{rStudent()} returns an \eqn{(n, d)}-\code{\link{matrix}}
  containing \eqn{n} samples of the
  \eqn{d}-dimensional multivariate \emph{t} distribution with location
  vector \code{loc} and scale matrix \code{scale} (a covariance matrix).
}
\details{
  This is a user friendly wrapper and calls \code{rnvmix(...,mix = 
  "inverse.gamma", df = df)}, see \code{\link{rnvmix}()}.
  
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used to obtain \code{factor} is the Cholesky
  decomposition via \code{\link{chol}()}. To this end, \code{scale}
  needs to have full rank. See \code{\link{factorize}()} for other options.

  The number of rows of \code{factor} are the dimension \eqn{d} of
  the sample. Typically (but not necessarily), \code{factor} is square.
  
}
\author{Marius Hofert, Erik Hintz}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{dStudent}()}, \code{\link{pStudent}()}, \code{\link{rnvmix}()}, \code{\link{factorize}()}
}
\examples{
## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))

## Draw random variates and compare
df <- 3.5
n <- 1000
set.seed(271)
X <- rStudent(n, df = df, scale = P) # providing scale
set.seed(271)
X. <- rStudent(n, df = df, factor = chol(P)) # providing the factor
stopifnot(all.equal(X, X.))

## Univariate case (dimension = number of rows of 'factor' = 1 here)
set.seed(271)
X.1d <- rStudent(n, df = df, factor = 1/2)
set.seed(271)
X.1d. <- rStudent(n, df = df, factor = 1)/2 # manual scaling
stopifnot(all.equal(X.1d, X.1d.))
}
\keyword{distribution}