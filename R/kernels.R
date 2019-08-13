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

#' Generates a update sequence accordignly
#' @noRd
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

#' Various kernel functions for MCMC
#' 
#' @param mu,scale Either a numeric vector or a scalar. Proposal mean and scale.
#' @param lb,ub Either a numeric vector or a scalar. Lower and upper bounds for
#' bounded kernels.
#' @param fixed Logical scalar or vector. When `TRUE` fixes the corresponding
#' parameter, avoiding new proposals.
#' @param scheme scheme in which proposals are made (see details).
#' @details
#' The objects `fmcmc_kernels` are environments that in general contain the 
#' following objects:
#' 
#' -  `proposal`: The function used to propose changes in the chain based
#'    on the current state. The function must return a vector of length equal
#'    to the number of parameters in the model.
#'    
#' -  `logration`: This function is called after a new state has been proposed,
#'    and is used to compute the log of the hastings ratio.
#'    
#'    In the case that the `logratio` function is not specified, then it is assumed
#'    that the transition kernel is symmetric, this is, logratio is then implemented
#'    as `function(env) {env$f1 - env$f0}`
#'    
#' -  `...`: Further objects that are used within those functions.
#' 
#' Both functions, `proposal` and `logratio`, receive a single argument, an
#' environment, which is passed by the [MCMC] function during each step using
#' the function [environment]. The passed
#' environment is actually the environemnt in which the `MCMC` function is running,
#' in particular, this environment contains the following objects:
#' 
#' \tabular{lcl}{
#' \strong{Object} \tab \tab \strong{Description} \cr
#' `i` \tab \tab Integer. The current iteration. \cr
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
#' For more defails see the vignnete `vignette("user-defined-kernels", "fmcmc")`.
#' 
#' @section Proposal scheme:
#'  
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
#' @name kernels
#' @aliases fmcmc_kernel fmcmc-kernel
#' @examples 
#' \dontrun{
#' # Example creating a multivariate normal kernel using the mvtnorm R package
#' # for a bivariate normal distribution
#' library(mvtnorm)
#' 
#' # Define your Sigma
#' sigma <- matrix(c(1, .2, .2, 1), ncol = 2)
#' 
#' # How does it looks like?
#' sigma
#'      [,1] [,2]
#' [1,]  1.0  0.2
#' [2,]  0.2  1.0
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
#' }
NULL

#' @export
#' @rdname kernels
#' @param proposal,logratio Functions. The function receives a single argument, an environment.
#' This functions are called later within [MCMC] (see details).
#' @param ... In the case of `kernel_new`, further arguments to be stored with
#' the kernel.
#' @section Creating your own kernels:
#' The function `kernel_new` is a helper function that allows creating
#' `fmcmc_kernel` which is used with the `MCMC` function. The `fmcmc_kernel`
#' are the backbone of the [MCMC] function.
kernel_new <- function(proposal, logratio = NULL, ...) {
  
  # Checks
  if (length(formals(proposal)) != 1L)
    stop("The `proposal` function should receive a single argument.", call. = FALSE)
  if (!is.null(logratio) && length(formals(logratio)) != 1L)
    stop("The `logratio` function should receive a single argument.", call. = FALSE)
  
  if (is.null(logratio))
    logratio <- function(env) {env$f1 - env$f0}
  
  env <- new.env()
  environment(proposal) <- env
  environment(logratio) <- env
  
  env$proposal <- proposal
  env$logratio <- logratio
  
  dots <- list(...)
  for (n in names(dots))
    env[[n]] <- dots[[n]]
  
  structure(env, class = c("environment", "fmcmc_kernel"))
  
}

#' @export
#' @rdname kernels
#' @param x An object of class `fmcmc_kernel`.
print.fmcmc_kernel <- function(x, ...) {
  
  cat("\nAn object of class fmcmc_kernel. The following two functions are available:\n\n")
  print(utils::ls.str(x))
  cat("\n")
  invisible(x)
  
}

#' @export
#' @rdname kernels
#' @param min.,max. Passed to [stats::runif].
#' @section Kernels:
#' The `kernel_unif` function provides a uniform transition kernel. This (symmetric)
#' kernel function by default adds the current status values between \[-1,1\].
kernel_unif <- function(
  min.  = -1.0,
  max.  = 1.0,
  fixed = FALSE,
  scheme = "joint"
  ) {
  
  k               <- NULL
  update_sequence <- NULL
  
  kernel_new(
    proposal = function(env) {
      
      if (env$i == 1L | is.null(k)) {
        
        k <<- length(env$theta0)
        
        # Checking boundaries
        min.  <<- check_dimensions(min., k)
        max.  <<- check_dimensions(max., k)
        fixed <<- check_dimensions(fixed, k)
        
        if (any(max. <= min.))
          stop("-max.- cannot be <= than -min.-.", call. = FALSE)
        
        # Setting the scheme in which the variables will be updated
        update_sequence <<- plan_update_sequence(
          k      = k,
          nsteps = env$nsteps,
          fixed  = fixed,
          scheme  = scheme
          )
        
        # It is easier to do the updates accordignly
        k <<- sum(update_sequence[1,])
        
      }
      
      theta1 <- env$theta0
      which. <- which(update_sequence[env$i,])
      theta1[which.] <- theta1[which.] + 
        stats::runif(k, min.[which.], max.[which.])
      theta1
    },
    min.  = min.,
    max.  = max.,
    fixed = fixed,
    k     = k,
    scheme = scheme,
    update_sequence = update_sequence
  )
}

#' @export
#' @rdname kernels
#' @section Kernels:
#' The `kernel_unif_reflective` is similar to `kernel_unif` with the
#' main difference that proposals are bounded to be within `[lb, ub]`.
kernel_unif_reflective <- function(
  min.  = -1.0,
  max.  = 1.0,
  lb    = min.,
  ub    = max.,
  fixed = FALSE,
  scheme = "joint"
) {
  
  k               <- NULL
  update_sequence <- NULL
  
  kernel_new(
    proposal = function(env) {
      
      if (env$i == 1L | is.null(k)) {
        
        k <<- length(env$theta0)
        
        # Checking boundaries
        min.  <<- check_dimensions(min., k)
        max.  <<- check_dimensions(max., k)
        ub    <<- check_dimensions(ub, k)
        lb    <<- check_dimensions(lb, k)
        fixed <<- check_dimensions(fixed, k)
        
        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        if (any(max. <= min.))
          stop("-max.- cannot be <= than -min.-.", call. = FALSE)
        
        # Setting the scheme in which the variables will be updated
        update_sequence <<- plan_update_sequence(
          k      = k,
          nsteps = env$nsteps,
          fixed  = fixed,
          scheme  = scheme
        )
        
        # It is easier to do the updates accordignly
        k <<- sum(update_sequence[1,])
        
      }
      
      # normal_reflective(env$theta0, lb, ub, mu, scale, fixed)
      which. <- which(update_sequence[env$i, ])
      any_reflective(
        randfun = stats::runif,
        # Parameters for rnorm
        n       = k,
        min     = min.[which.],
        max     = max.[which.],
        x       = env$theta0,
        # Parameters for the function
        lb      = lb,
        ub      = ub,
        which   = which.
      )
      
    },
    lb    = lb,
    ub    = ub,
    min.  = min.,
    max.  = max.,
    fixed = fixed,
    k     = k,
    scheme = scheme,
    update_sequence = update_sequence
  )
}

#' @export
#' @rdname kernels
#' @section Kernels:
#' The `kernel_normal` function provides the cannonical normal kernel
#' with symmetric transition probabilities.
kernel_normal <- function(
  mu    = 0,
  scale = 1,
  fixed = FALSE,
  scheme = "joint"
  ) {
  
  k               <- NULL
  update_sequence <- NULL
  
  kernel_new(
    proposal = function(env) {
      
      if (env$i == 1L | is.null(k)) {
        
        k <<- length(env$theta0)
        
        # Checking boundaries
        mu    <<- check_dimensions(mu, k)
        scale <<- check_dimensions(scale, k)
        fixed <<- check_dimensions(fixed, k)
        
        # Setting the scheme in which the variables will be updated
        update_sequence <<- plan_update_sequence(
          k      = k,
          nsteps = env$nsteps,
          fixed  = fixed,
          scheme  = scheme
        )
        
        if (sum(!fixed) == 0L)
          stop("The number of parameters to update, i.e. not fixed, cannot be ",
               "zero. Check the value -fixed- in the kernel initialization.", 
               call. = FALSE)
        
        k <<- sum(update_sequence[1,])
        
      }
      
      # Which to update (only whichever is to be updated in the sequence)
      # of updates.
      which. <- which(update_sequence[env$i, ])
      
      theta1 <- env$theta0
      theta1[which.] <- theta1[which.] +
        stats::rnorm(k, mean = mu[which.], sd = scale[which.])
      theta1
    },
    logratio = function(env) env$f1 - env$f0,
    mu       = mu,
    scale    = scale,
    k        = k,
    fixed    = fixed,
    update_sequence = update_sequence,
    scheme    = scheme
    )
}

#' @export 
#' @rdname kernels
#' @section Kernels:
#' The `kernel_normal_reflective` implements the normal kernel with reflective
#' boundaries. Lower and upper bounds are treated using reflecting boundaries, this is, 
#' if the proposed \eqn{\theta'} is greater than the \code{ub}, then \eqn{\theta' - ub}
#' is substracted from \eqn{ub}. At the same time, if it is less than \code{lb}, then
#' \eqn{lb - \theta'} is added to \code{lb} iterating until \eqn{\theta} is within
#' \code{[lb, ub]}.
#' 
#' In this case, the transition probability is symmetric (just like the normal
#' kernel).
kernel_normal_reflective <- function(
  mu    = 0,
  scale = 1,
  lb    = -.Machine$double.xmax,
  ub    = .Machine$double.xmax,
  fixed = FALSE,
  scheme = "joint"
) {
  
  k               <- NULL
  update_sequence <- NULL
  
  kernel_new(
    proposal = function(env) {
      
      # Checking whether k exists. We should do this only once. When doing this
      # we have to make sure that the length of lb and ub are according to the
      # length of model parameter.
      #
      # We also want to restart k in the first run. Since it is based on
      # environments, the user may move this around... so it is better to just
      # restart this every time that the MCMC function starts from scratch.
      if (env$i == 1L | is.null(k)) {
        
        k <<- length(env$theta0)
        
        # Checking boundaries
        ub    <<- check_dimensions(ub, k)
        lb    <<- check_dimensions(lb, k)
        mu    <<- check_dimensions(mu, k)
        scale <<- check_dimensions(scale, k)
        fixed <<- check_dimensions(fixed, k)

        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        
        # Setting the scheme in which the variables will be updated
        update_sequence <<- plan_update_sequence(
          k      = k,
          nsteps = env$nsteps,
          fixed  = fixed,
          scheme  = scheme
        )
        
        k <<- sum(update_sequence[1, ])
        
      }
      
      # Which to update (only whichever is to be updated in the sequence)
      # of updates.
      which. <- which(update_sequence[env$i, ])
      any_reflective(
        randfun = stats::rnorm,
        # Parameters for rnorm
        n       = k,
        mean    = mu[which.],
        sd      = scale[which.],
        x       = env$theta0,
        # Parameters for the function
        lb      = lb,
        ub      = ub,
        which   = which.
        )
    },
    logratio = function(env) env$f1 - env$f0,
    mu       = mu,
    scale    = scale, 
    ub       = ub, 
    lb       = lb, 
    fixed    = fixed,
    k        = k,
    scheme    = scheme,
    update_sequence = update_sequence
  )
  
}

any_reflective <- function(
  randfun,
  ...,
  x,
  lb,
  ub,
  which
) {
  
  x[which] <- x[which] + randfun(...)
  
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


