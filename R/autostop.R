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
#' @param freq Integer scalar. Frequency of checking.
#' @param ... Further arguments passed to the method.
#' @return A function passed to [MCMC] to check automatic convergence.
#' @name convergence-checker
#' @aliases automatic-stop
NULL

#' Removes invariant columns
#' @noRd
rm_invariant <- function(x) {
  
  variances <- which(stats::sd(
    if (is.list(x)) do.call(rbind, x) else x
    )^2 < 1e-10)
  
  if (length(variances) > 0) {
    
    if (length(variances) == coda::nvar(x))
      return(FALSE)
    
    for (i in 1:coda::nchain(x))
      x[[i]] <- x[[i]][, -variances,drop=FALSE]
    
  }
  
  x
}

#' @export 
#' @rdname convergence-checker
convergence_gelman <- function(
  freq            = 1000L,
  threshold       = 1.10,
  check_invariant = TRUE,
  ...
  ) {

  structure(function(x) {
    
    if (coda::nchain(x) > 1L) {
      
      # Checking invariant
      if (check_invariant) 
        x <- rm_invariant(x)
        
      # Computing gelman test
      d <- tryCatch(coda::gelman.diag(x, ...), error = function(e) e)
      
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
    
  }, freq = freq)
  
}

#' @rdname convergence-checker
#' @details
#' In the case of `convergence_geweke`, `threshold` sets the p-value 
#' for the null \eqn{H_0: Z = 0}, i.e. equal means between the first and last
#' chunks of the chain. See [coda::geweke.diag]. This implies that the higher
#' the threshold, the lower the probability of stopping the chain.
#' 
#' In the case that the chain has more than one parameter, the algorithm will
#' return true if and only if the test fails to reject the null for all the
#' parameters.
#' @export
convergence_geweke <- function(freq = 1000L, threshold=.025, check_invariant=TRUE,...) {
  
  
  structure(function(x) {
    
    if (coda::nchain(x) > 1L) 
      stop("The `geweke` convergence check is only available with runs of a single chain.",
           call. = FALSE)
    
    # Checking invariant
    if (check_invariant) 
      x <- rm_invariant(x)
    
    d <- tryCatch(coda::geweke.diag(x, ...)$z, error = function(e) e)
    
    if (inherits(d, "error")) {
      
      warning("At ", coda::niter(x), " `geweke.diag` failed to be computed.",
              " Will skip and try with the next batch.", call. = FALSE,
              immediate. = TRUE)
      
      return(FALSE)
      
    }
    
    d <- 1 - stats::pnorm(-abs(d))*2.0
    
    if (any(!is.finite(d)))
      return(FALSE)
    
    all(d > threshold)
  }, freq=freq)
  
}


#' @rdname convergence-checker
#' @details
#' For the `convergence_heildel`, see [coda::heidel.diag] for details.
#' @export
convergence_heildel <- function(freq = 1000L, ..., check_invariant=TRUE) {
  
  structure(function(x) {
    
    if (coda::nchain(x) > 1L) 
      stop("The `heidel` convergence check is only available with runs of a single chain.",
           call. = FALSE)
    
    # Checking invariant
    if (check_invariant) 
      x <- rm_invariant(x)
    
    d <- tryCatch(coda::heidel.diag(x, ...)[, c("stest", "htest")], error = function(e) e)
    
    if (inherits(d, "error")) {
      
      warning("At ", coda::niter(x), " `heidel.diag` failed to be computed.",
              " Will skip and try with the next batch.", call. = FALSE,
              immediate. = TRUE)
      
      return(FALSE)
      
    }
    
    if (any(!is.finite(d)))
      return(FALSE)
    
    # The tests return a 1
    all(d == 1)
  }, freq = freq)
  
}

#' @rdname convergence-checker
#' @details The `convergence_auto` function is the default and is just a wrapper
#' for `convergence_gelman` and `convergence_geweke`. This function returns a 
#' convergence checker that will be either of the other two depending on whether
#' `nchains` in `MCMC` is greater than one--in which case it will use the Gelman
#' test--or not--in which case it will use the Geweke test.
#' @export
convergence_auto <- function(freq = 1000L) {
  
  gelman_conv <- convergence_gelman(freq)
  geweke_conv <- convergence_geweke(freq)
  
  structure(function(x) {
    
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
    
  }, freq = freq)
  
}

#' Run MCMC with convergence checker
#' @noRd
#' @param expr The expression to parse
#' @param conv_checker A function to be used as a convergence checker.
with_autostop <- function(expr, conv_checker) {
  
  # Getting the parent environment
  freq   <- attr(conv_checker, "freq")
  parenv <- parent.frame()
  
  # Retrieving parameters from the MCMC call
  nsteps    <- parenv$nsteps
  nchains   <- parenv$nchains

  # Correcting the freq
  if (freq*2 > nsteps) 
    freq <- 0L

  # Capturing the expression
  expr <- sys.call()[[2]]
  
  # Calculating lengths. The bulk vector sets what will be
  # nsteps in each call. This excludes burnin
  bulks <- rep(freq, (nsteps - parenv$burnin) %/% freq)
  if ((nsteps - parenv$burnin) %% freq)
    bulks <- c(bulks, (nsteps - parenv$burnin) - sum(bulks))
  
  # We need to add the burnin to the first
  bulks[1] <- bulks[1] + parenv$burnin
  
  # Do while no convergence
  converged <- FALSE
  i         <- 0L
  ans       <- NULL
  for (i in seq_along(bulks)) {
    
    # Updating the nsteps argument
    parenv$nsteps <- bulks[i]

    if (i > 1) {
      parenv$burnin  <- 0L
      parenv$initial <- ans[coda::niter(ans),]
      if (coda::is.mcmc.list(ans))
        parenv$initial <- do.call(rbind, parenv$initial)
    }
      
    # Running the MCMC and adding it to the tail
    tmp <- eval(expr, envir = parenv)
    if (is.list(tmp) & !coda::is.mcmc.list(tmp))
      tmp <- coda::as.mcmc.list(tmp)
    
    # Appending retults
    ans <- append_chains(ans, tmp)
    
    if ((converged <- conv_checker(ans))) {
      message(
        "Convergence has been reached with ", sum(bulks[1:i]), " steps (",
        coda::niter(ans), " final count of samples)."
        )
      break
    }
    
  }
  
  # Did it converged?
  if (!is.null(conv_checker) && (i == length(bulks) & !converged))
    message("No convergence reached after ", sum(bulks[1:i]), " steps (",
            coda::niter(ans), " final count of samples).")
  
  # Returning
  return(ans)
  
}

