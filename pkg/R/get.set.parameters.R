### get.set.parameters() #######################################################

##' @title  Retrieve algorithm specific default parameters and overwrite them
##' @return list with default values for all functions in the 'nvmix' package:
##              $method 
##              $mean.sqrt.mix 
##              $precond 
##              $pnvmix.abstol 
##              $increment 
##              $dnvmix.abstol 
##              $dnvmix.abstol.log 
##              $max.iter.newton 
##              $newton.conv.abstol 
##              $newton.df.abstol 
##              $newton.logdens.abstol 
##              $max.iter.rqmc 
##              $CI.factor 
##              $fun.eval 
##              $B 
##' @note newton.logdens.abstol not very small as it is only needed to get the
##' the next iteration in the Newton procedure *unless* it is being called 
##' from dnvmixcop() in which case this tolerance is changed there. 
##' @author Erik Hintz and Marius Hofert

get.set.parameters <- function(control = list()){
  ## Set up default controls
  control.int <- list(
    ## For pnvmix(): 
    method = "sobol", 
    mean.sqrt.mix = NULL, 
    precond = TRUE, 
    pnvmix.abstol = 1e-3, 
    increment = "doubling", 
    ## For dnvmix():
    dnvmix.abstol = 1e-3, 
    dnvmix.abstol.log = 1e-3, # not used (yet)
    ## For qnvmix():
    max.iter.newton = 40, 
    newton.conv.abstol = 1e-4,
    newton.df.abstol = 1e-4,
    newton.logdens.abstol = 1e-2, 
    ## For all (randomized) algorithms:
    max.iter.rqmc = 500, 
    CI.factor = 3.3,
    fun.eval = c(2^7, 1e8), 
    B = 12)
  
  if(length(control) > 0){
    ## If input provided, grab input controls and overwrite:
    names.control <- names(control.int)
    control.int[(names.provided <- names(control))] <- control
    
    ## Did the user provide something that is not used?
    if (length(unmatched <- names.provided[!names.provided %in% names.control])) 
      warning("unknown names in control: ", paste(unmatched, collapse = ", "))
    
    ## Now some checkings:
    stopifnot(is.logical(control.int$precond),
              control.int$pnvmix.abstol >= 0,
              control.int$dnvmix.abstol >= 0,
              control.int$dnvmix.abstol.log >= 0,
              control.int$max.iter.newton >= 0,
              control.int$newton.conv.abstol >= 0,
              control.int$newton.df.abstol >= 0,
              control.int$newton.logdens.abstol >= 0,
              control.int$max.iter.rqmc >= 2,
              control.int$CI.factor >= 0,
              length(control.int$fun.eval) == 2, control.int$fun.eval >= 0,
              control.int$B > 1) # If B=1 error estimates are NA => need B>1
    
    ## Check if 'method' and 'increment' were provided correctly
    control.int$method <- match.arg(control.int$method, 
                                    choices = c("sobol", "ghalton", "PRNG"))
    control.int$increment <- match.arg(control.int$increment, 
                                       choices = c("doubling", "num.init"))
  }
  control.int
}


