#' Checks the dimensions of the parameter according to k.
#' @noRd
check_dimensions <- function(x, k) {
  
  # Getting the name of the variable
  name <- deparse(match.call()$x)
  
  if (length(x) > 1L &&  length(x) != k)
    stop("Incorrect length of -", name, "-.", call. = FALSE)
  
  if (length(x) == 1L && k > 1L)
    return(rep(x, k))
  
  return(x)
  
  
}

#' Paramaters' update sequence
#' @param k Integer. Number of parameters
#' @param nsteps Integer. Number of steps.
#' @param fixed Logical scalar or vector of length `k`. Indicates which parameters
#' will be treated as fixed or not. Single values are recycled.
#' @param scheme Scheme in which the proposals are made (see details).
#' 
#' @details 
#' The parameter `scheme` present on the currently available kernels sets the way
#' in which proposals are made. By default, `scheme = "joint"`, proposals are done
#' jointly, this is, at each step of the chain we are proposing new states for
#' each parameter of the model. When `scheme = "ordered"`, a sequential update schema
#' is followed, in which, at each step of the chain, proposals are made one
#' variable at a time, If `scheme = "random"`, proposals are also made one
#' variable at a time but in a random scheme.
#' 
#' Finally, users can specify their own sequence of proposals for the variables
#' by passing a numeric vector to `scheme`, for example, if the user wants to make
#' sequential proposals following the scheme 2, 1, 3, then scheme must be set to
#' be `scheme = c(2, 1, 3)`.
#' 
#' @return A logical vector of size `nsteps` x `k`.
#' @export
plan_update_sequence <- function(k, nsteps, fixed, scheme) {
  
  # Setting the scheme in which the variables will be updated
  if (length(scheme) > 1L && is.numeric(scheme)) {
    
    # Is it the right length?
    if (length(scheme) != sum(!fixed))
      stop(
        "When setting the update scheme, it should have the same length ",
        "as the number of variables that will not be fixed. ",
        "Right now length(scheme) = ", length(scheme), " while sum(!fixed) = ",
        sum(!fixed),".", call. = FALSE
        )
    
    # Is the full sequence included?
    test <- which(!(which(!fixed) %in% scheme))
    if (length(test))
      stop(
        "One or more variables was not included in the ordering sequence. ",
        "The full list follows: ", paste(test, collapse=", "), ". ",
        "Only variables that are not fixed can be included in this list.",
        call. = FALSE
        )
    
    update_sequence <- matrix(FALSE, nrow = nsteps, ncol = k)
    suppressWarnings(update_sequence[cbind(1:nsteps, scheme)] <- TRUE)
    
    
  } else if (scheme == "joint") {
    
    update_sequence <- matrix(TRUE, nrow = nsteps, ncol = k)
    
    for (j in which(fixed))
      update_sequence[, j] <- FALSE
  
  } else if (scheme == "ordered") {
    
    update_sequence <- matrix(FALSE, nrow = nsteps, ncol = k)
    suppressWarnings(update_sequence[cbind(1:nsteps, which(!fixed))] <- TRUE)
    
  } else if (scheme == "random") {
    
    update_sequence <- matrix(FALSE, nrow = nsteps, ncol = k)
    update_sequence[cbind(
      1:nsteps,
      sample(which(!fixed), nsteps, TRUE))
      ] <- TRUE
    update_sequence <- update_sequence
    
  } else {
    
    stop(
      "-scheme- update must be either an integer sequence, 'joint', 'ordered', ",
      "or 'random'.",
      call. = FALSE
      )
    
  }
  
  if (sum(update_sequence[1,]) == 0L)
    stop("The number of parameters to update, i.e. not fixed, cannot be ",
         "zero. Check the value -fixed- in the kernel initialization.", 
         call. = FALSE)
    
  
  return(update_sequence)
  
}

#' Create Personalized Transition Kernels for MCMC
#' 
#' The function `kernel_new` is a helper function that allows creating
#' `fmcmc_kernel` objects which are passed to the [MCMC()] function.
#' 
#' @param proposal,logratio Functions. Both receive a single argument, an environment.
#' This functions are called later within [MCMC] (see details).
#' @param kernel_env Environment. This will be used as the main container of the
#' kernel's components. It is returned as an object of class `c("environment", "fmcmc_kernel")`.
#' @param ... In the case of `kernel_new`, further arguments to be stored with
#' the kernel.
#' 
#' @details
#' 
#' The objects `fmcmc_kernels` are environments that in general contain the 
#' following objects:
#' 
#' -  `proposal`: The function used to propose changes in the chain based
#'    on the current state. The function must return a vector of length equal
#'    to the number of parameters in the model.
#'    
#' -  `logratio`: This function is called after a new state has been proposed,
#'    and is used to compute the log of the Hastings ratio.
#'    
#'    In the case that the `logratio` function is not specified, then it is assumed
#'    that the transition kernel is symmetric, this is, log-ratio is then implemented
#'    as `function(env) {env$f1 - env$f0}`
#'    
#' -  `...`: Further objects that are used within those functions.
#' 
#' Both functions, `proposal` and `logratio`, receive a single argument, an
#' environment, which is passed by the [MCMC] function during each step using
#' the function [environment()].
#' 
#' The passed environment is actually the environment in which the `MCMC`
#' function is running, in particular, this environment contains the following
#' objects:
#' 
#' \tabular{lcl}{
#' \strong{Object} \tab \tab \strong{Description} \cr
#' `i` \tab \tab Integer. The current iteration. \cr
#' `chain_id` \tab \tab Integer. The current chain (relevent for some kernels). \cr
#' `theta1` \tab \tab Numeric vector. The last proposed state. \cr
#' `theta0` \tab \tab Numeric vector. The current state \cr
#' `f`\tab \tab The log-unnormalized posterior function (a wrapper of `fun` passed 
#' to [MCMC]). \cr
#' `f1` \tab \tab The last value of `f(theta1)` \cr
#' `f0` \tab \tab The last value of `f(theta0)` \cr
#' `kernel` \tab \tab The actual `fmcmc_kernel` object. \cr
#' `ans` \tab \tab The matrix of samples defined up to `i - 1`.
#' }
#' 
#' These are the core component of the `MCMC` function. The following block
#' of code is how this is actually implemented in the package:
#' 
#' ```
#' for (i in 1L:nsteps) {
#'   # Step 1. Propose
#'   theta1[] <- kernel$proposal(environment())
#'   f1       <- f(theta1)
#'   
#'   # Checking f(theta1) (it must be a number, can be Inf)
#'   if (is.nan(f1) | is.na(f1) | is.null(f1)) 
#'     stop(
#'       "fun(par) is undefined (", f1, ")",
#'       "Check either -fun- or the -lb- and -ub- parameters.",
#'       call. = FALSE
#'     )
#'   
#'   # Step 2. Hastings ratio
#'   if (R[i] < kernel$logratio(environment())) {
#'     theta0 <- theta1
#'     f0     <- f1
#'   }
#'   
#'   # Step 3. Saving the state
#'   ans[i,] <- theta0
#'   
#' }
#' ```
#' 
#' For more details see the vignette `vignette("user-defined-kernels", "fmcmc")`.
#' 
#' @family kernels
#' @aliases fmcmc_kernel kernels
#' @examples 
#' 
#' # Example creating a multivariate normal kernel using the mvtnorm R package
#' # for a bivariate normal distribution
#' library(mvtnorm)
#' 
#' # Define your Sigma
#' sigma <- matrix(c(1, .2, .2, 1), ncol = 2)
#' 
#' # How does it looks like?
#' sigma
#' #      [,1] [,2]
#' # [1,]  1.0  0.2
#' # [2,]  0.2  1.0
#' 
#' # Create the kernel
#' kernel_mvn <- kernel_new(
#'   proposal = function(env) {
#'   env$theta0 + as.vector(mvtnorm::rmvnorm(1, mean = 0, sigma = sigma.))
#'   },
#'   sigma. = sigma
#' )
#' 
#' # As you can see, in the previous call we passed sigma as it will be used by
#' # the proposal function
#' # The logaratio function was not necesary to be passed since this kernel is
#' # symmetric.
#' 
#' @return An environment of class `fmcmc_kernel` which contains the following:
#' 
#' - `proposal` A function that receives a single argument, an environment. This
#'   is the proposal function used within [MCMC()].
#'   
#' - `logratio` A function to compute log ratios of the current vs the proposed
#'   step of the chain. Also used within [MCMC()].
#' 
#' - `...` Further arguments passed to `kernel_new`.
#' 
#' @references 
#' Brooks, S., Gelman, A., Jones, G. L., & Meng, X. L. (2011). Handbook of
#' Markov Chain Monte Carlo. Handbook of Markov Chain Monte Carlo.
#' @export
kernel_new <- function(
  proposal,
  ... ,
  logratio   = NULL,
  kernel_env = new.env(hash = TRUE)
  ) {
  
  # Checks
  if (length(formals(proposal)) != 1L)
    stop(
      "The `proposal` function should receive a single arguments (an environment).",
      call. = FALSE
      )
  if (!is.null(logratio) && length(formals(logratio)) != 1L)
    stop(
      "The `logratio` function should receive a single argument (an environment) ",
      call. = FALSE
      )
  
  if (is.null(logratio))
    logratio <- function(env) {env$f1 - env$f0}
  
  environment(proposal) <- kernel_env
  environment(logratio) <- kernel_env
  
  kernel_env$proposal <- proposal
  kernel_env$logratio <- logratio
  
  dots <- list(...)
  for (n in names(dots))
    kernel_env[[n]] <- dots[[n]]
  
  structure(kernel_env, class = c("environment", "fmcmc_kernel"))
  
}

#' @export
#' @rdname kernel_new
#' @param x An object of class `fmcmc_kernel`.
print.fmcmc_kernel <- function(x, ...) {
  
  cat("\nAn environment of class fmcmc_kernel:\n\n")
  print(utils::ls.str(x))
  cat("\n")
  invisible(x)
  
}

#' Reflective Boundaries
#' 
#' Adjust a proposal according to its support by reflecting it. This is the workhorse
#' of [kernel_normal_reflective] and [kernel_unif_reflective]. It is intended
#' for internal use only. 
#' 
#' @param x A numeric vector. The proposal
#' @param lb,ub Numeric vectors of length `length(x)`. Lower and upper bounds.
#' @param which Integer vector. Index of variables to be updated.
#' 
#' @return An adjusted proposal vector.
#' @export
reflect_on_boundaries <- function(
  x,
  lb,
  ub,
  which
) {
  
  test_above <- which(x[which] > ub[which])
  test_below <- which(x[which] < lb[which])
  
  d <- ub - lb
  
  if (length(test_above)) {
    
    # Putting in context
    test_above <- which[test_above]
    
    # Direction of update
    d_above <- x[test_above] - ub[test_above]
    odd     <- (d_above %/% d[test_above]) %% 2
    d_above <- d_above %% d[test_above]
    
    x[test_above] <- (lb[test_above] + d_above)*odd + 
      (ub[test_above] - d_above)*(1-odd)
  } 
  
  if (length(test_below)) {
    
    # Putting in context
    test_below <- which[test_below]
    
    # Direction of update
    d_below <- lb[test_below] - x[test_below]
    odd     <- (d_below %/% d[test_below]) %% 2
    d_below <- d_below %% d[test_below]
    
    x[test_below] <- (ub[test_below] - d_below)*odd + 
      (lb[test_below] + d_below)*(1-odd)
    
  } 
  
  x
  
}
