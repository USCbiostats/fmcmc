#' Automatic stop
#' 
#' @param threshold Numeric value. A Gelman statistic below the threshold
#' will return `TRUE`.
#' @param check_invariant Logical. When `TRUE` the function only computes
#' the gelman diagnostic using variables with greater than `1e-10` variance
#' @return A function passed to [MCMC] to check automatic convergence.
#' @export
gelman_convergence <- function(threshold = 1.10, check_invariant=TRUE) {

  function(x) {
    
    if (coda::nchain(x) > 1L) {
      
      # Checking invariant
      if (check_invariant) {
        
        variances <- which(stats::sd(do.call(rbind, x))^2 < 1e-10)
        
        if (length(variances) > 0) {
        
          if (length(variances) == coda::nvar(x))
            return(FALSE)
            
          for (i in 1:coda::nchain(x))
            x[[i]] <- x[[i]][, -variances,drop=FALSE]
        
        }
        
      }
      
      # Computing gelman test
      d <- tryCatch(coda::gelman.diag(x), error = function(e) e)
      
      if (inherits(d, "error")) {
        
        warning("At ", coda::niter(x), " `gelman.diag` failed to be computed.",
                " Will skip and try with the next batch.", call. = FALSE,
                immediate. = TRUE)
        
        return(FALSE)
        
      }
      
      # Depending on multivariable or not
      val <- ifelse(coda::nvar(x) > 1L, d$mpsrf, d$psrf[1,"Point est."])
      
      # Is it lower than the treshold?
      return(val < threshold)
        
      
    } else {
      
      stop("Convergence test with the Gelman is only available when `nchains` > 1L.",
           call. = FALSE)
      
    }
    
  }
  
}

with_autostop <- function(expr, conv_checker) {
  
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
  
  # Correcting the autostop
  if (autostop*2 > nbatch) {
    autostop <- 0L
    parenv$autostop <- 0L
  }

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
  i         <- 0L
  ans       <- NULL
  for (i in seq_along(bulks)) {
    
    # Updating the nbatch argument
    parenv$nbatch <- bulks[i]
    
    if (i > 1) {
      parenv$burnin  <- 0L
      parenv$initial <- do.call(rbind, ans[coda::niter(ans),])
    }
      
        # Running the MCMC and adding it to the tail
    tmp <- eval(expr, envir = parenv)
    if (is.list(tmp) & !coda::is.mcmc.list(tmp))
      tmp <- coda::as.mcmc.list(tmp)
    
    # Appending retults
    ans <- append_chains(ans, tmp)
    
    if ((autostop > 0L) && (converged <- conv_checker(ans))) {
      message(
        "Convergence has been reached with ", coda::niter(ans), " iterations (",
        sum(bulks[1:i]), " steps)."
        )
      break
    }
    
  }
  
  # Did it converged?
  if (autostop && (i == length(bulks) & !converged))
    warning("No convergence reached.", call. = FALSE)
  
  # Returning
  return(ans)
  
}

# autostop(1 + 1)