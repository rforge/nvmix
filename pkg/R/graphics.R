### QQ-plot of Mahalanobis distances for visual GOF ############################

##' @title QQ Plot of Mahalanobis distances versus their theoretical quantiles
##' @param x (n,d) data matrix
##' @param qmix see ?pnvmix()
##' @param loc see ?pnvmix()
##' @param scale see ?pnvmix()
##' @param fitnvmix_object object of class 'fitnvmix'; if provided, x, qmix, 
##'        loc and scale are ignored.
##' @param trafo_to_normal logical, if TRUE the 
##'        underlying Mahalanobis distances are mapped to normals by a probability-
##'        quantile-transform so that the resulting QQ plot is essentially a normal
##'        QQ plot.
##' @param test character; specifying if (and which) GoF test shall be performed.
##' @param bootstrap_pars  list with elements 'B' (Bootstrap sample size for computing
##'        CIs) and 'level' (confidence level). 
##' @param plot logical if the result should be plotted.       
##' @param verbose see ?pnvmix() 
##' @param control see ?pnvmix()
##' @param digits number of digits for the test statistic and the p-value
##' @param plotpars list; see ?get_set_qqplot_param()
##' @param ... see ?pnvmix()
##' @return invisibly returns an object of class 'qqplot_maha'. 
##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
qqplot_maha <- function(x, qmix, loc, scale, fitnvmix_object, 
                        trafo_to_normal = FALSE, test = c("KS", "AD", "none"),
                        bootstrap_pars = list(B = 500, level = 0.95),
                        plot = TRUE, verbose = TRUE, control = list(), 
                        digits = max(3, getOption("digits") - 3),
                        plotpars = list(), ...){
   ## Initialize and check inputs
   control <- get_set_param(control)
   control$newton.df.reltol <- control$qqplot.df.reltol
   test <- match.arg(test)
   plotpars <- get_set_qqplot_param(plotpars)
   call <- match.call() # for return
   if(hasArg(fitnvmix_object)){
      stopifnot(class(fitnvmix_object) == "fitnvmix")
      ## Grab elements from the fitnvmix_object and then call recursively 
      return(qqplot_maha(
         x = fitnvmix_object$data, qmix = fitnvmix_object$qmix, 
         loc = fitnvmix_object$loc, scale = fitnvmix_object$scale,
         trafo_to_normal = trafo_to_normal, test = test, 
         bootstrap_pars = bootstrap_pars, plot = plot, verbose = verbose, 
         control = control, digits = digits, plotpars = plotpars, 
         nu = fitnvmix_object[[1]]))
      ## Note: 'fitnvmix_object$qmix' is either a function(u, nu) or a string
      ## ("pareto", "inverse.gamma", "constant"). In all cases, 'nu' can
      ## be as the name of the argument supplied via the ellipsis argument 
   } else {
      if(!is.matrix(x)) x <- rbind(x)
      notNA <- rowSums(is.na(x)) == 0
      x     <- x[notNA,, drop = FALSE] # non-missing data (rows)
      n     <- nrow(x)
      d     <- ncol(x)
      ## Obtain sorted Mahalanobis distances (X-loc)^T scale^{-1} (X-loc)
      maha2 <- sort(mahalanobis(x, center = loc, cov = scale))
   }
   pp <- ppoints(n)
   ## Compute theoretical quantiles 
   theo_quant <- if(trafo_to_normal){
      maha2 <- qnorm(pgammamix(maha2, qmix = qmix, d = d, control = control, 
                               verbose = verbose, ...)) # maha2 mapped to N(0,1)
      x_ <- qnorm(c(0.25, 0.75))
      theodf <- function(x) pnorm(x) 
      list(q = qnorm(pp), log.density = dnorm(pp, log = TRUE))
   } else {
      x_ <- qgammamix(c(0.25, 0.75), qmix = qmix, d = d, control = control,
                     verbose = verbose, q.only = TRUE, ...)
      theodf <- function(x) pgammamix(x, qmix = qmix, d = d, control = control,
                                      ...) # hypothesized cdf for below
      qgammamix(pp, qmix = qmix, d = d, control = control,
                verbose = verbose, q.only = FALSE, stored.values = NULL, ...)
   }
   ## Already compute intercept and slope as in 'qqline()'
   y_ <- quantile(maha2, c(0.25, 0.75))
   slope <- diff(y_) / diff(x_) 
   int <- x_[1L] - slope * y_[1L]
   ## Obtain asymptotic CI (see Fox (2008), pp 35-36)
   logSE <- (log(pp) + log(1-pp) - log(n))/2 - 
      theo_quant$log.density # logarithmic SE 
   asymptSE <- exp(logSE) 
   ## Obtain bootstrap CI
   boot_CI <- if(bootstrap_pars$B > 1){
      alpha.2 <- (1 - bootstrap_pars$level) / 2
      bootsample <- sapply(1:bootstrap_pars$B, function(i) maha2[sort(
         sample(1:n, size = n, replace = TRUE))]) 
      apply(bootsample, 1, quantile, probs = c(alpha.2, 1 - alpha.2), 
            names = FALSE) # (2, n) matrix 
   } else NULL
   ## Compute test statistic
   testout <- switch(test, 
                     "KS" = {
                        tmp <- ks.test(maha2, y = theodf)
                        tmp$method <- "Kolmogorov-Smirnov GoF test" # shorter
                        tmp
                     }, 
                     "AD" = {
                        ad.test(maha2, distr.fun = theodf)
                     },
                     NULL) # for "none" or anything else
   ## Create S3-class object of type 'qqplot_maha'
   out <- class_qqplot_maha(
      maha2 = maha2, theo_quant = theo_quant$q, boot_CI = boot_CI, 
      trafo_to_normal = trafo_to_normal, asymptSE = asymptSE, test = test, 
      testout = testout, int = int, slope = slope, B = bootstrap_pars$B,
      call = call)
   ## Create plot
   if(plot) {
      plot(out, digits = digits, plotpars = plotpars)
   } else {
      invisible(out)
   }
}



# OLD (for now just commented out, to be removed soon)
# ##' @title QQ Plot of Mahalanobis distances versus their theoretical quantiles
# ##' @param x (n,d) data matrix
# ##' @param qmix see ?pnvmix()
# ##' @param loc see ?pnvmix()
# ##' @param scale see ?pnvmix()
# ##' @param control see ?pnvmix()
# ##' @param plot.diag logical; if TRUE the curve f(x) = x is plotted additionally
# ##' @param plot.band logical if CI bands shall be plotted as well 
# ##' @param verbose logical if warnings from underlying 'qgammamix()' shall be
# ##'        thrown
# ##' @param plot logical if QQ plot shall be plotted     
# ##' @param control list of algorithm specific parameters, see ?get_set_param()
# ##'        and ?fitnvmix
# ##' @param ... see ?pnvmix()
# ##' @return invisibly returns a list of three: 'maha2' (sorted squared mahalanobis
# ##'         distances obtained from 'x', sorted); 'q' (theoretical quantiles);
# ##'         'ldens' (log-density at 'q')
# ##' @author Erik Hintz, Marius Hofert, Christiane Lemieux
# qqplot_maha <- function(x, qmix, loc, scale, plot.diag = TRUE, plot.band = TRUE,
#                         verbose = TRUE, plot = TRUE, control = list(), ...)
# {
#    ## Initialize and check inputs
#    control <- get_set_param(control)
#    control$newton.df.reltol <- control$qqplot.df.reltol
#    if(!is.matrix(x)) x <- rbind(x)
#    notNA <- rowSums(is.na(x)) == 0
#    x     <- x[notNA,, drop = FALSE] # non-missing data (rows)
#    n     <- nrow(x)
#    d     <- ncol(x)
#    ## Obtain sorted Mahalanobis distances (X-loc)^T scale^{-1} (X-loc)
#    maha2 <- sort(mahalanobis(x, center = loc, cov = scale))
#    pp <- ppoints(n)
#    ## Obtain theoretical quantiles
#    theo_quant <-
#       qgammamix(pp, qmix = qmix, d = d, control = control,
#                 verbose = verbose, q.only = FALSE, stored.values = NULL, ...)
#    ## Plot
#    if(plot){
#       plot(theo_quant$q, maha2, xlab = "Theoretical quantiles",
#            ylab = "Sample quantiles", main = "")
#       ## Add a diagonal
#       if(plot.diag) lines(theo_quant$q, theo_quant$q, lty = 2)
#    }
#    ## Add the CI bands
#    if(plot.band & plot){
#       ## Compute SE of the empirical quantiles (see Fox (2008), pp 35-36)
#       logSE <- (log(pp) + log(1-pp) - log(n))/2 - 
#          theo_quant$log.density # logarithmic SE 
#       SE <- exp(logSE) 
#       ## Add the band
#       lines(theo_quant$q, maha2 + 2 * SE, lty = 3)
#       lines(theo_quant$q, maha2 - 2 * SE, lty = 3)
#    }
#    invisible(list(maha2 = maha2, q = theo_quant$q, 
#                   ldens = theo_quant$log.density))
# }


### S3 class functions and methods #############################################

#' Function to define S3 class 'qqplot_maha'
#'
#' @param maha2 observed Mahalanobis distances 
#' @param theo_quant theoretical quantiles 
#' @param method string, "bootstrap" or "normaltrafo"
#' @param asymptSE asymptotic standard errors
#' @param test string ("KS", "AD" or "none")
#' @param testout list returned by ks.test() or ad.test() or NULL if test = "none"
#' @param int intercept of the fitted line
#' @param slope slope of the fitted line
#' @param call function call
#' @return S3 object of class 'qqplot_maha'
#' @author Erik Hintz
class_qqplot_maha <- function(maha2, theo_quant, boot_CI, trafo_to_normal, asymptSE, 
                              test, testout, int, slope, B, call){
   res <- list(maha2 = maha2, theo_quant = theo_quant, boot_CI = boot_CI,
               trafo_to_normal = trafo_to_normal, asymptSE = asymptSE, test = test, 
               testout = testout, int = int, slope = slope, B = B, call = call)
   ## Return object of class "qqplot_maha"
   structure(res, class = "qqplot_maha")
}

## Method 'print' for S3 class 'qqplot_maha'
print.qqplot_maha <- function(x, ..., digits = max(3, getOption("digits") - 3)){
   ## Print function call to qqplot_maha()
   cat("Call: ", deparse(x$call), "\n", sep = "")
   cat("\n")
   ## Provide sample size
   cat(paste0("Input: ", length(x$maha2), " squared Mahalanobis distances."), "\n")
   ## Check 'trafo_to_normal'
   if(x$trafo_to_normal)
      cat("Squared Mahalanobis distances were transformed to N(0, 1).", "\n", sep = "")
   cat("\n")
   ## Check if a test was performed
   test_text <- if(is.null(x$testout)) "No GoF test performed." else {
      paste0(x$testout$method, ": D = ", round(x$testout$statistic, digits), ", p = ", 
             format(x$testout$p.value, scientific = TRUE, digits = digits), ".")
   }
   cat(test_text, "\n")
   cat("\n")
   cat("Computed results stored in the object:", "\n")
   cat("- theoretical quantiles in $theo_quant;", "\n")
   cat("- sorted, squared Mahalanobis distances in $maha2;", "\n")
   cat("- estimated, asymptotic standard errors in $asymptSE;", "\n")
   if(!is.null(x$boot_CI))
      cat(paste0("- Bootstrap CIs (estimated from ", x$B, " resamples) in $boot_CI;"), "\n")
   if(!is.null(x$testout))
      cat("- GoF test results in $testout;", "\n")
   invisible(x)
}


## Method 'plot' for S3 class 'qqplot_maha'
plot.qqplot_maha <- 
   function(x, ..., digits = max(3, getOption("digits") - 3), plotpars = list())
{
   ## Set parameters
   plotpars <- get_set_qqplot_param(plotpars)
   ## Construct the text for the axis
   axistext <- if(is.null(x$testout)) "" else {
      paste0(x$testout$method, ": D = ", round(x$testout$statistic, digits), 
             ", p = ", format(x$testout$p.value, scientific = TRUE, 
                               digits = digits), ".")}
   ## Plot
   plot(x$theo_quant, x$maha2, xlab = plotpars$xlab, ylab = plotpars$ylab,
        xlim = plotpars$xlim, ylim = plotpars$ylim, main = plotpars$main,
        sub = plotpars$sub, log = plotpars$log, pch = plotpars$pch,
        col = plotpars$col[1])
   ## Add line
   if(plotpars$plot_line)
      abline(x$int, x$slope, lty = plotpars$lty[1], col = plotpars$col[2])
   ## Add asymptotic CI 
   for(i in c(-1, 1))
      lines(x$theo_quant, x$maha2 + i * x$asymptSE, lty = plotpars$lty[2], 
            col = plotpars$col[3])
   lgnd <- "Asymptotic CI"
   numlegend <- 1 # number of elements in the legend
   ## Add bootstrap CI
   if(!is.null(x$boot_CI)){ 
      lgnd <- c(lgnd, "Bootstrap CI")
      numlegend <- 2
      for(i in 1:2) 
         lines(x$theo_quant, x$boot_CI[i, ], lty = plotpars$lty[3], 
               col = plotpars$col[4])
      
   }
   if(plotpars$plot_legend) legend("topleft", lgnd, 
                          lty = plotpars$lty[2:(1+numlegend)], 
                          col = plotpars$col[3:(2+numlegend)], 
                          bty = "n")
   if(plotpars$plot_test) mtext(axistext, side = 4) # empty if no test performed
   ## Invisibly return input object
   invisible(x)
}

