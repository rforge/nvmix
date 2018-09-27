\name{rStudent}
\alias{rStudent}
\title{Random Number Generator for the Multivariate Student t Distribution}
\description{
  Generate vectors of random variates from the multivariate Student
  \emph{t} distribution.
}
\usage{
rStudent(n, df, loc = rep(0, d), scale = diag(2), factor = factorize(scale))
}
\arguments{
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{df}{degress of freedom (positive real or Inf in which case a
    sample from a multivariate normal with mean vector \code{loc}
    and covariance matrix \code{scale} is generated).}
  \item{loc}{location vector of dimension \eqn{d}; this equals the mean
    vector of a random vector following the specified Student \emph{t}
    distribution if (and only if) the latter exists.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)}
    (\eqn{d = 2} per default); this equals the covariance matrix of a
    random vector following the specified Student \emph{t} distribution divided
    by  the mean of the mixing variable if (and only if) the former exists.}
  \item{factor}{\code{\link{matrix}} \eqn{R} with \eqn{d}
    rows such that \eqn{R^T R} equals \code{scale}. \eqn{R} is internally
    multiplied from the right to an \eqn{(n, d)}-matrix of independent
    standard normals to obtain a sample from a multivariate normal
    with zero mean vector and covariance matrix \code{scale}.}
}
\value{
  \eqn{(n, d)}-\code{\link{matrix}} containing \eqn{n} samples of the
  \eqn{d}-dimensional multivariate Student \emph{t} distribution with
  location vector \code{loc} and scale matrix \code{scale}.
}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used to obtain \code{factor} is the Cholesky
  decomposition via \code{\link{chol}()}. To this end, \code{scale}
  needs to have full rank. See \code{\link{factorize}()} for other options.

  The number of rows of \code{factor} are the dimension \eqn{d} of
  the sample. Typically (but not necessarily), \code{factor} is square.

  \code{rStudent()} and \code{rNorm()} are wrappers of
  \code{rnvmix(, mix = "inverse.gamma", df = df)} and
  \code{rnvmix(, mix = "constant", df = df)}, respectively.
}
\author{Erik Hintz and Marius Hofert}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{rnvmix}()} (of which \code{rStudent()} is a wrapper),
  \code{\link{dStudent}()}, \code{\link{pStudent}()}, \code{\link{factorize}()}
}
\examples{
## Generate a random correlation matrix in d dimensions
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
X.1d  <- rStudent(n, df = df, factor = 1/2)
set.seed(271)
X.1d. <- rStudent(n, df = df, factor = 1)/2 # manual scaling
stopifnot(all.equal(X.1d, X.1d.))
}
\keyword{distribution}