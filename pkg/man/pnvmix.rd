\name{pnvmix}
\alias{pnvmix}
\title{Distribution Function of Multivariate Normal Variance Mixtures}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t} for non-integer degrees of freedom).
}
\usage{
pnvmix(a, b, R, nu, gam = 3.3, eps = 0.001, Nmax = 1e8, N = 10,
       n_init = 2^5, precond = TRUE)
}
\arguments{
  \item{a}{vector of length \eqn{d}.}
  \item{b}{vector of length \eqn{d}.}
  \item{R}{positive definite \eqn{(d,d)}-covariance matrix.}
  \item{nu}{degress of freedom (any positive value).}
  \item{eps}{error tolerance.}
  \item{gam}{determines the stopping criterion of the algorithm; it will
    run until \eqn{err < gam * eps}.}
  \item{Nmax}{maximum number of function evaluations, can be used to
    control run time.}
  \item{N}{Number of repetitions to get an error estimate in the
    randomized quasi-Monte Carlo approach.}
  \item{n_init}{size of the first point set being used to estimate
    the probability.}
  \item{precond}{\code{\link{logical}} indicating if preconditioning
    is applied, that is, reordering the variables.}
}
\value{
  \code{pnvmix()} returns a list of length four, containing the
  the estimated probabilities, the number of iterations, the total
  number of function evaluations and an error estimate.
}
\details{
  Note that this procedure calls underlying C code. Currently, the
  dimensions \eqn{d\ge 16510}{d >= 16510} are not supported.
}
\author{Marius Hofert, Erik Hintz and Christiane Lemieux}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Evaluate t_{3.5} distribution function
a <- runif(d) * sqrt(d) * (-3) # random lower limit
b <- runif(d) * sqrt(d) * 3 # random upper limit
pt <- pnvmix(a = a, b = b, R = P, nu = 3.5)
stopifnot(all.equal(pt$Prob, 0.8061, tol = 5e-4))
}
\keyword{distribution}