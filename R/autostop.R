#' Convergence Monitoring
#' 
#' Built-in set of functions to be used in companion with the argument 
#' `conv_checker` in [MCMC]. These functions are not intended to be used
#' in a context other than the `MCMC` function.
#'  
#' @param threshold Numeric value. A Gelman statistic below the threshold
#' will return `TRUE`.
#' @param check_invariant Logical. When `TRUE` the function only computes
#' the Gelman diagnostic using variables with greater than `1e-10` variance.
#' @param freq Integer scalar. Frequency of checking.
#' @param ... Further arguments passed to the method.
#' @return A function passed to [MCMC] to check automatic convergence.
#' @name convergence-checker
#' @aliases automatic-stop
NULL

#' @export
#' @rdname convergence-checker
#' @description The object `LAST_CONV_CHECK` is an environment that holds
#' information regarding the convergence checker used. This information can be
#' updated every time that the `conv_checker` function is called by `MCMC` using
#' the functions `convergence_data_set` and `convergence_msg_set`. The function
#' `convergence_data_get` is just a wrapper of [get()].
#' 
#' The `msg` member of `LAST_CONV_CHECK` is resetted before `conv_checker` is
#' called.
LAST_CONV_CHECK <- structure(new.env(), class = c("fmcmc_run_conv_check", "environment"))
assign("msg", NA_character_, envir = LAST_CONV_CHECK)


#' @export
print.fmcmc_run_conv_check <- function(x, ...) {
  
  cat("LAST_CONV_CHECK holds the following information:\n")
  print(utils::ls.str(x))
  
  invisible(x)
}


convergence_data_flush <- function() {
  rm(list = ls(envir = LAST_CONV_CHECK, all.names = TRUE), envir = LAST_CONV_CHECK)
  assign("msg", NA_character_, envir = LAST_CONV_CHECK)
}

#' @rdname convergence-checker
#' @section Building a convergence checker:
#' Convergence checkers are simply a function that receives as argument a matrix
#' (or list of them) with sampled values, and returns a logical scalar with the
#' value `TRUE` if the chain converged. An example of a personalized convergence
#' checker is provided below. The frequency with which the check is performed is
#' retrieved from the attribute `"freq"` from the convergence checker function,
#' i.e., `attr(..., "freq")`. If missing, convergence will be checked halfway
#' the number of steps in the chain, i.e., `floor(nsteps/2)`.
#' 
#' @examples 
#' # Example 1: Presonalized conv checker --------------------------------------
#' # Dummy rule, if acceptance rate is near between .2 and .3.
#' convergence_example <- function(x) {
#'   arate <- 1 - coda::rejectionRate(x)
#'   all(
#'     abs(arate - .25) < .05
#'   )
#' }
#' 
#' # Tell fmcmc::MCMC what is the frequency
#' attr(convergence_example, "freq") <- 2e3
#' 
#' set.seed(223)
#' x <- rnorm(1000)
#' y <- x * 2 + rnorm(1000)
#' logpost <- function(p) {
#'   sum(dnorm(y, mean = x * p, log = TRUE))
#' }
#' 
#' ans <- MCMC(
#'   initial = 0, fun = logpost, nsteps = 5e4,
#'   kernel= kernel_ram(),
#'   conv_checker = convergence_example
#' )
#' 
#' # Example 2: Adding information ---------------------------------------------
#' # Here we do two things: Save a value and set a message for the user
#' convergence_example_with_info <- structure(function(x) {
#'   arate <- 1 - coda::rejectionRate(x)
#'   
#'   # Saving a value
#'   if (!exists("arates", envir = LAST_CONV_CHECK)) {
#'     convergence_data_set(list(arates = arate))
#'   } else {
#'     convergence_data_set(list(
#'       arates = rbind(convergence_data_get("arates"), arate)
#'     ))
#'   }
#'   
#'   # Setting up the message
#'   convergence_msg_set(
#'     sprintf("Current Avg. Accept. Rate: %.2f", mean(arate))
#'   )
#'   
#'   all(
#'     abs(arate - .25) < .05
#'   )
#' }, freq = 2000)
#' 
#' 
#' ans <- MCMC(
#'   initial = 0, fun = logpost, nsteps = 5e4,
#'   kernel= kernel_ram(),
#'   conv_checker = convergence_example_with_info,
#'   seed = 112,
#'   progress = FALSE
#' )

#' @export
#' @param x In the case of `convergence_data_set`, a named list. For 
#' `convergence_data_get`, a character vector.
convergence_data_set <- function(x) {
  
  if (!is.list(x))
    stop("-x- must be a named list.")
  if (length(names(x)) != length(x))
    stop("Not all the elements in -x- are named.")
  
  if ("msg" %in% names(x))
    stop("'msg' should be set using -convergence_msg_set-")
  
  for (i in names(x))
    assign(i, x[[i]], envir = LAST_CONV_CHECK)
  
  invisible(NULL)
  
}

#' @export
#' @rdname convergence-checker
convergence_data_get <- function(x) {
  if (length(x) > 1L) {
    mget(x, envir = LAST_CONV_CHECK)
  } else
    get(x, envir = LAST_CONV_CHECK)
}

#' @export
#' @param msg Character scalar. Personalized message to print.
#' @rdname convergence-checker
convergence_msg_set <- function(msg = NA_character_) {
  
  if (!is.character(msg))
    stop("-msg- must be a character.")
  else if (length(msg) != 1L)
    stop("-msg- must be of length 1.")
  
  assign("msg", msg, envir = LAST_CONV_CHECK)
  
  invisible(NULL)
  
}

#' @export
#' @rdname convergence-checker
convergence_msg_get <- function() {
  get("msg", envir = LAST_CONV_CHECK)
}

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
#' @details `convergence_gelman` is a wrapper of [coda::gelman.diag()].
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
      
      # Updating the convergence checker
      if (exists("dat", envir = LAST_CONV_CHECK)) {
        convergence_data_set(list(
          dat = c(convergence_data_get("dat"), stats::setNames(list(d), stats::end(x)))
          ))
      } else {
        convergence_data_set(list(dat = stats::setNames(list(d), stats::end(x))))
      }
      
      # Depending on multivariate or not
      val <- ifelse(coda::nvar(x) > 1L, d$mpsrf, d$psrf[1,"Point est."])
      convergence_msg_set(sprintf("Gelman-Rubin's R: %.4f.", val))
      
      
      # Is it lower than the threshold?
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
    
    # Updating the convergence checker
    if (exists("dat", envir = LAST_CONV_CHECK)) {
      convergence_data_set(list(
        dat = c(convergence_data_get("dat"), stats::setNames(list(d), stats::end(x)))
      ))
    } else {
      convergence_data_set(list(dat = stats::setNames(list(d), stats::end(x))))
    }
    
    convergence_msg_set(sprintf("avg Geweke's Z: %.4f.", mean(d[is.finite(d)])))
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
      stop(
        "The -heidel- convergence check is only available with runs of a single chain.",
        call. = FALSE
        )
    
    # Checking invariant
    if (check_invariant) 
      x <- rm_invariant(x)
    
    d <- tryCatch(coda::heidel.diag(x, ...), error = function(e) e)
    
    if (inherits(d, "error")) {
      
      warning(
        "At ", coda::niter(x), " -coda::heidel.diag- failed to be computed.",
        " Will skip and try with the next batch.",
        call.      = FALSE,
        immediate. = TRUE
        )
      
      return(FALSE)
      
    }
    
    tests <- d[, c("stest", "htest")]

    # Updating the convergence checker
    if (exists("dat", envir = LAST_CONV_CHECK)) {
      convergence_data_set(list(
        dat = c(convergence_data_get("dat"), stats::setNames(list(d), stats::end(x)))
      ))
    } else {
      convergence_data_set(list(dat = stats::setNames(list(d), stats::end(x))))
    }
    
    convergence_msg_set(sprintf("Heidel's Avg. pval: %.2f", mean(d[,"pvalue"])))
    
    
    if (any(!is.finite(tests)))
      return(FALSE)
    
    # The tests return a 1
    all(tests == 1)
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
    
    who <- as.character(sys.call(-1L)[[1]])
    if (!length(who) || (who != "with_autostop"))
      warning("This function should not be used in a context other than ",
              "the argument `conv_checker` in `MCMC`.")
    
    env <- parent.frame()
    
    if (env$nchains > 1L)
      gelman_conv(x)
    else
      geweke_conv(x)
    
  }, freq = freq)
  
}
