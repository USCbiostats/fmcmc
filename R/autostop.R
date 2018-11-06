


with_autostop <- function(expr, has_converged) {
  
  # Getting the parent environment
  parenv <- parent.frame()
  
  # Retrieving parameters from the MCMC call
  nbatch   <- parenv$nbatch
  nchains  <- parenv$nchains
  autostop <- parenv$autostop

  # Checking frequency to measure batch
  if (!is.numeric(autostop))
    stop("The `autostop` parameter must be a number. autostop=", autostop,
         call. = FALSE)
  else if (length(autostop) != 1L)
    stop("The `autostop` parameter must be of length 1. autostop=", autostop,
         call. = FALSE)
  
  # Capturing the expression
  expr <- sys.call()[[2]]
  
  # Calculating lengths. The bulk vector sets what will be
  # nbatch in each call. This excludes burnin
  if (autostop > 0L) {
    bulks <- rep(autostop, (nbatch - parenv$burnin) %/% autostop)
    if ((nbatch - parenv$burnin) %% autostop)
      bulks <- c(bulks, (nbatch - parenv$burnin) - sum(bulks))
    
    # We need to add the burnin to the first
    bulks[1] <- bulks[1] + parenv$burnin
    
  } else
    # If there's no autostop, then no need to split
    bulks <- nbatch
  
    # Do while no convergence
  converged <- FALSE
  i <- 0L
  ans <- NULL
  while (!converged) {
    
    # Updating the nbatch argument
    parenv$nbatch <- bulks[i <- i + 1]
    
    # Running the MCMC
    ans <- c(ans, eval(expr, envir = parenv))
    
    if (converged <- has_converged(ans)) {
      
      message(
        "Convergence has been reached with ", coda::niter(ans), " iterations."
        )
      
    }
      
    
  }
  
  # Returning
  return(ans)
  
}

# autostop(1 + 1)