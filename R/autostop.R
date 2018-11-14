#' Convergence Monitoring
#' 
#' Built-in set of functions to be used in companion with the argument 
#' `conv_checker` in [MCMC]. These functions are not intended to be used
#' in a context other than the `MCMC` function.
#'  
#' @param threshold Numeric value. A Gelman statistic below the threshold
#' will return `TRUE`.
#' @param check_invariant Logical. When `TRUE` the function only computes
#' the gelman diagnostic using variables with greater than `1e-10` variance.
#' @param ... Further arguments passed to the method.
#' @return A function passed to [MCMC] to check automatic convergence.
#' @name convergence-checker
#' @aliases automatic-stop
NULL

#' Removes invariant columns
#' @noRd
rm_invariant <- function(x) {
  
  variances <- which(stats::sd(do.call(rbind, x))^2 < 1e-10)
  
  if (length(variances) > 0) {
    
    if (length(variances) == coda::nvar(x))
      return(FALSE)
    
    for (i in 1:coda::nchain(x))
      x[[i]] <- x[[i]][, -variances,drop=FALSE]
    
  }
  
  x
}

#' @export
#' @param autoburnin Logical. Passed to [coda::gelman.diag]. By default this
#' option is set to `FALSE` since the [MCMC] function already a first chunk of
#' steps via the argument `burnin`.
#' 
#' @rdname convergence-checker
gelman_convergence <- function(threshold = 1.10, check_invariant=TRUE, autoburnin=FALSE,...) {

  function(x) {
    
    if (coda::nchain(x) > 1L) {
      
      # Checking invariant
      if (check_invariant) 
        x <- rm_invariant(x)
        
      # Computing gelman test
      d <- tryCatch(coda::gelman.diag(x, autoburnin = autoburnin, ...), error = function(e) e)
      
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

#' @rdname convergence-checker
#' @details
#' In the case of `geweke_convergence`, `threshold` sets the p-value 
#' for the null \eqn{H0: Z = 0}, i.e. equal means between the first and last
#' chunks of the chain. See [coda::geweke.diag]. This implies that the higher
#' the threshold, the lower the probability of stopping the chain.
#' 
#' In the case that the chain has more than one parameter, the algorithm will
#' return true if and only if the test fails to reject the null for all the
#' parameters.
#' @export
geweke_convergence <- function(threshold=.025, check_invariant=TRUE,...) {
  
  
  function(x) {
    
    if (coda::nchain(x) > 1L) 
      stop("The `geweke` convergence check is only available with runs of a single chain.",
           call. = FALSE)
    
    # Checking invariant
    if (check_invariant) 
      x <- rm_invariant(x)
    
    d <- tryCatch(coda::geweke.diag(x, ...)[[1]]$z, error = function(e) e)
    
    if (inherits(d, "error")) {
      
      warning("At ", coda::niter(x), " `geweke.diag` failed to be computed.",
              " Will skip and try with the next batch.", call. = FALSE,
              immediate. = TRUE)
      
      return(FALSE)
      
    }
    
    d <- stats::pnorm(d)
    d <- ifelse(d > .5, 1 - d, d)*2
    
    if (any(!is.finite(d)))
      return(FALSE)
    
    all(d > threshold)
  }
  
}

#' @rdname convergence-checker
#' @details The `auto_convergence` function is the default and is just a wrapper
#' of `gelman_convergence` and `geweke_convergence`. This function returns a 
#' convergence checker that will be either of the other two depending on wether
#' `nchains` in `MCMC` is greater than one--in which case it will use the Gelman
#' test--or not--in which case it will use the Geweke test.
#' @export
auto_convergence <- function() {
  
  gelman_conv <- gelman_convergence()
  geweke_conv <- geweke_convergence()
  
  function(x) {
    
    if (as.character(sys.call(-1L)[[1]]) != "with_autostop") {
      warning("This function should not be used in a context other than ",
              "the argument `conv_checker` in `MCMC`.")
      
      return(FALSE)
    }
    
    env <- parent.frame()
    
    if (env$nchains > 1L)
      gelman_conv(x)
    else
      geweke_conv(x)
    
  }
  
}

with_autostop <- function(expr, conv_checker) {
  
  # Getting the parent environment
  parenv <- parent.frame()
  
  # Retrieving parameters from the MCMC call
  nsteps   <- parenv$nsteps
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
  if (autostop*2 > nsteps) {
    autostop <- 0L
    parenv$autostop <- 0L
  }

  # Capturing the expression
  expr <- sys.call()[[2]]
  
  # Calculating lengths. The bulk vector sets what will be
  # nsteps in each call. This excludes burnin
  if (autostop > 0L) {
    bulks <- rep(autostop, (nsteps - parenv$burnin) %/% autostop)
    if ((nsteps - parenv$burnin) %% autostop)
      bulks <- c(bulks, (nsteps - parenv$burnin) - sum(bulks))
    
    # We need to add the burnin to the first
    bulks[1] <- bulks[1] + parenv$burnin
    
  } else
    # If there's no autostop, then no need to split
    bulks <- nsteps
  
    # Do while no convergence
  converged <- FALSE
  i         <- 0L
  ans       <- NULL
  for (i in seq_along(bulks)) {
    
    # Updating the nsteps argument
    parenv$nsteps <- bulks[i]
    
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
        "Convergence has been reached with ", sum(bulks[1:i]), " steps (",
        coda::niter(ans), " final count of observations)."
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