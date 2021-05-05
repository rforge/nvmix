### QQ-plot of Mahalanobis distances for visual GOF ############################

##' @title QQ Plot of Mahalanobis distances versus their theoretical quantiles
##' @param x (n,d) data matrix
##' @param qmix see ?pnvmix()
##' @param loc see ?pnvmix()
##' @param scale see ?pnvmix()
##' @param control see ?pnvmix()
##' @param plot.diag logical; if TRUE the curve f(x) = x is plotted additionally
##' @param plot.band logical if CI bands shall be plotted as well 
##' @param verbose logical if warnings from underlying 'qgammamix()' shall be
##'        thrown
##' @param control list of algorithm specific parameters, see ?get_set_param()
##'        and ?fitnvmix
##' @param ... see ?pnvmix()
##' @return invisibly returns a list of three: 'maha2' (sorted squared mahalanobis
##'         distances obtained from 'x', sorted); 'q' (theoretical quantiles);
##'         'ldens' (log-density at 'q')
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
qqplot_maha <- function(x, qmix, loc, scale, plot.diag = TRUE, plot.band = TRUE,
                        verbose = TRUE, control = list(), ...)
{
   ## Initialize and check inputs
   control <- get_set_param(control)
   control$newton.df.reltol <- control$qqplot.df.reltol
   if(!is.matrix(x)) x <- rbind(x)
   notNA <- rowSums(is.na(x)) == 0
   x     <- x[notNA,, drop = FALSE] # non-missing data (rows)
   n     <- nrow(x)
   d     <- ncol(x)
   ## Obtain sorted Mahalanobis distances (X-loc)^T scale^{-1} (X-loc)
   maha2 <- sort(mahalanobis(x, center = loc, cov = scale))
   pp <- ppoints(n)
   ## Obtain theoretical quantiles
   theo_quant <-
      qgammamix(pp, qmix = qmix, d = d, control = control,
                verbose = verbose, q.only = FALSE, stored.values = NULL, ...)
   ## Plot
   plot(theo_quant$q, maha2, xlab = "Theoretical quantiles",
        ylab = "Sample quantiles", main = "")
   ## Add a diagonal
   if(plot.diag) lines(theo_quant$q, theo_quant$q, lty = 2)
   ## Add the CI bands
   if(plot.band){
      ## Compute SE of the empirical quantiles (see Fox (2008), pp 35-36)
      logSE <- (log(pp) + log(1-pp) - log(n))/2 - 
         theo_quant$log.density # logarithmic SE 
      SE <- exp(logSE) 
      ## Add the band
      lines(theo_quant$q, maha2 + 2 * SE, lty = 3)
      lines(theo_quant$q, maha2 - 2 * SE, lty = 3)
   }
   invisible(list(maha2 = maha2, q = theo_quant$q, 
                  ldens = theo_quant$log.density))
}
