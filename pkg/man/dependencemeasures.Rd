\name{dependencemeasures}
\alias{lambda_gStudent}
\title{Dependence Measures for grouped normal variance mixture copulas}
\description{
  Computation of rank correlation coefficients Spearman's rho and Kendall's tau
  for grouped normal variance mixture copulas as well as computation of the 
  (lower and upper) tail dependence coefficient of a grouped \emph{t} copula.
}
\usage{
lambda_gStudent(df, rho, control = list(), verbose = TRUE)
}
\arguments{
  \item{df}{either scalar or \eqn{2}-\code{\link{vector}} giving the degrees-of-
  freedoms for the t copula; if provided as scalar, the copula is an (ungrouped) 
  t copula and \code{lambda_gStudent()} uses a closed formula.}
  \item{rho}{\eqn{n}-\code{\link{vector}} giving the \eqn{\rho} parameters of the copula.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{verbose}{
    \code{\link{logical}}
    indicating whether a warning is given if the required precision
    has not been reached.}
}
\value{
  \code{lambda_gStudent()} returns
  a \code{\link{numeric}} \eqn{n}-vector with the computed
  dependence measure with corresponding attributes 
  \code{"abs. error"} and \code{"rel. error"}(error estimates of the RQMC estimator) 
  and \code{"numiter"} (number of iterations).
}
\details{
  In the ungrouped case, \code{lambda_gStudent()} computes the tail dependence
  coefficient \eqn{lambda} based on the knwon formula 
  \code{2 * pt( -sqrt( (df + 1)*(1 - rho) / (1 + rho)), df = df + 1)}.
  
  In the grouped case, RQMC methods are used to efficiently approximate the 
  integral given in Eq. (26) of Luo and Shevchenko (2010). 
  
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2019),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \url{https://arxiv.org/abs/1911.03017}.
  
  Luo, X. and Shevchenko, P. (2010).
  The \emph{t} copula with multiple parameters of degrees of freedom: bivariate 
  characteristics and application to risk management.
  \emph{Quantitative Finance} 10(9),
  1039-1054.
}
\seealso{
  \code{\link{dgStudentcopula}()}, \code{\link{pgStudentcopula}()},  
  \code{\link{rgStudentcopula}()}
}
\examples{
## Create a plot displaying 'lambda' as a function of the copula parameter
## for various choices of the degrees-of-freedom 
df <- c(3, 6, 9)
df_ <- list( rep(df[1], 2), rep(df[2], 2), rep(df[3], 2), # ungrouped 
             c(df[1], df[2]), c(df[1], df[3]), c(df[2], df[3])) # grouped  
l.df_ <- length(df_)
rho <- seq(from = -0.99, to = 0.99, length.out = 112) # rho parameters
set.seed(1) # for reproducibilty
lambdas <-
   sapply(seq_len(l.df_), function(i) nvmix:::lambda_gStudent(df_[[i]], rho = rho))
lgnd <- character(length(df_))
plot(NA, xlim = range(rho), ylim = range(lambdas), xlab = expression(rho), 
     ylab = expression(lambda))
for(i in seq_len(l.df_)){
   lines(rho, lambdas[, i], col = i, lty = i)
   lgnd[i] <- if(df_[[i]][1] == df_[[i]][2]) paste0("df = ", df_[[i]][1]) else
      paste0("df1 = ", df_[[i]][1], ", df2 = ", df_[[i]][2])
}
legend("topleft", lgnd, col = seq_len(l.df_), lty = seq_len(l.df_),
       bty = 'n')
## If called with 'df' a 1-vector, closed formula for lambda is used => check
lambda.true <- sapply(1:3, function(i) nvmix:::lambda_gStudent(df_[[i]][1], rho = rho))
stopifnot(max(abs( lambda.true - lambdas[, 1:3])) < 4e-4)
}
\keyword{distribution}