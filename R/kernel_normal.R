#' Gaussian Transition Kernel
#' @template mu-scale
#' @template scheme
#' @details
#' The `kernel_normal` function provides the canonical normal kernel
#' with symmetric transition probabilities.
#' @export
#' @return An object of class [fmcmc_kernel].
#' @family kernels
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
#' @rdname kernel_normal
#' @family kernels
#' @details
#' The `kernel_normal_reflective` implements the normal kernel with reflective
#' boundaries. Lower and upper bounds are treated using reflecting boundaries, this is, 
#' if the proposed \eqn{\theta'} is greater than the \code{ub}, then \eqn{\theta' - ub}
#' is subtracted from \eqn{ub}. At the same time, if it is less than \code{lb}, then
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
      theta1 <- env$theta0
      
      # Proposal
      theta1[which.] <- theta1[which.] + 
        stats::rnorm(k, mean = mu[which.], sd = scale[which.])
      
      # Reflecting
      reflect_on_boundaries(
        x       = theta1,
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
