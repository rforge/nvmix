lambda_gStudent <- function(df, rho, control = list(), verbose = TRUE)
{

   ## 1. Checks and set-up  ####################################################
   stopifnot(all(df > 0), all(rho > -1), all(rho < 1), is.logical(verbose))
   ## 'df' can have length 1 (ungrouped) or length 2 (grouped)
   l.df <- length( df <- as.vector(df) )
   stopifnot(1 <= l.df, l.df <= 2)
   l.rho <- length(rho <- as.vector(rho)) # 'rho' can be a vector

   ## In ungrouped case, closed formula for lambda is available
   if(l.df == 1){
      res <- 2 * pt( -sqrt( (df + 1)*(1 - rho) / (1 + rho)), df = df + 1)
      attr(res, "abs. error") <- rep(0, l.rho)
      attr(res, "rel. error") <- rep(0, l.rho)
      return(res)
   }
   ## Grab and declare variables for RQMC estimation below
   #control <- get_set_param(control)
   control = list(lambda.abstol = 1e-3)
   ## Absolute or relative tolerance?
   tol <- if(is.na(control$lambda.abstol)){
      do.reltol <- TRUE # use relative tolerance
      stopifnot(control$lambda.reltol > 0)
      control$lambda.reltol
   } else {
      do.reltol <- FALSE
      control$lambda.abstol
   }
   control <- get_set_param(list())
   ## Grab method and prepare variables for RQMC
   method <- control$method
   B      <- control$B
   increment <- control$increment
   if(method == "sobol") seeds_ <- sample(1:(1e5*B), B) # B seeds for 'sobol()'
   rqmc.estimates   <- matrix(0, ncol = l.rho, nrow = B)
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   error <- rep(NA, l.rho) # store error for all 'rho'
   ## Additional variables needed if the increment chosen is "doubling"
   if(increment == "doubling") {
      if(method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Initialize max.error to > tol to enter while loop
   max.error <- tol + 42
   total.fun.evals <- 0
   numiter <- 0
   current.n <- control$fun.eval[1]
   ## Constant needed repeatedly
   B_df <-
      c((2^(df[2]/2-df[1]/2) * gamma((1+df[2])/2) / gamma((1+df[1])/2))^(1/df[2]),
        (2^(df[1]/2-df[2]/2) * gamma((1+df[1])/2) / gamma((1+df[2])/2))^(1/df[1]))

   ## 2. Actual computation ####################################################
   ## while() runs until precision 'tol' is reached or the number of function
   ## evaluations exceed fun.eval[2] or the number of iterations exceed
   ## control$max.iter.rqmc. In each iteration, B RQMC estimates of
   ## of the tail dependence coefficient lambda are computed; if 'rho' is a vector
   ## the same mixing realizations are used for all elements in 'rho'

   while(max.error > tol & total.fun.evals < control$fun.eval[2] &
         numiter < control$max.iter.rqmc)
   {

      ## Get B RQCM estimates
      for(b in 1:B)
      {
         ## 2.1 Get the point set ###########################################
         U <- switch(method,
                     "sobol" = {
                        if(increment == "doubling") {
                           qrng::sobol(n = current.n, d = 1,
                                       randomize = "digital.shift",
                                       seed = seeds_[b],
                                       skip = (useskip * current.n))
                        } else {
                           qrng::sobol(n = current.n, d = 1,
                                       randomize = "digital.shift",
                                       seed = seeds_[b],
                                       skip = (numiter * current.n))
                        }
                     },
                     "ghalton" = {
                        qrng::ghalton(n = current.n, d = 1,
                                      method = "generalized")
                     },
                     "PRNG" = {
                        matrix(runif( current.n ), ncol = 1)
                     })
         ## 2.2 Evaluate the integrand at the (next) point set #################
         ## Realizations chi^2(df[1]+1) and chi^2(df[2]+1) via qgamma():
         chisq_sample <-
            cbind(qgamma(U, shape = (df[1] + 1)/2, scale = 2),
                  qgamma(U, shape = (df[2] + 1)/2, scale = 2))
         ## Compute next estimate based on 'chisq_sample'
         next.estimate <-
            colMeans(sapply(1:l.rho, function(i){
               pnorm(- (B_df[1] * chisq_sample[, 1]^(df[1]/(2*df[2])) -
                           rho[i]*sqrt(chisq_sample[, 1]))*(1-rho[i]^2)^(-1/2)) +
                  pnorm(- (B_df[2] * chisq_sample[, 2]^(df[2]/(2*df[1])) -
                              rho[i]*sqrt(chisq_sample[, 2])) *(1-rho[i]^2)^(-1/2))
            })) # l.rho - vector

         ## 2.3 Update RQMC estimates ##########################################
         rqmc.estimates[b, ] <-
            if(increment == "doubling") {
               ## In this case both, rqmc.estimates[b] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b, ] + next.estimate) / denom
            } else {
               ## In this case, rqmc.estimates[b] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b, ] + next.estimate) / (numiter + 1)
            }
      } # end for(b in 1:B)

      ## Update of various variables
      ## Number of function evaluations
      total.fun.evals <- total.fun.evals + B * current.n
      if(increment == "doubling") {
         ## Change 'denom' and 'useksip' (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Update error depending on 'do.reltol'
      error <- if(!do.reltol) { # absolute error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)
      } else { # relative error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)/
            .colMeans(rqmc.estimates, m = B, n = l.rho)
      }
      max.error <- max(error)
      numiter <- numiter + 1 # update counter
   } # while ()

   ## 3. Finalize and return ###################################################
   res <- .colMeans(rqmc.estimates, m = B, n = l.rho)
   ## Handle warnings
   reached <- (error <= tol)
   if(any(!reached) & verbose > 0) {
      ii <- which(!reached)
      if(verbose == 1) {
         strng <- if(length(ii) > 6) {
            paste0(paste(head(ii), collapse = ", "), ",...")
         } else {
            paste(ii, collapse = ", ")
         }
         warning("Tolerance not reached for entries ",strng," of 'rho'; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument.")
      } else {
         for(i in 1:length(ii)) {
            warning(sprintf("Tolerance not reached for entries %d of 'rho'; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument", ii[i]))
         }
      }
   }
   ## Compute absolute and relative error for attributes
   abserror <- if(do.reltol){
      relerror <- error
      error * res
   } else {
      relerror <- error / res # 'error' is absolute error
      error
   }
   ## Return
   attr(res, "abs. error") <- abserror
   attr(res, "rel. error") <- relerror
   attr(res, "numiter") <- numiter
   res
}

