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
      theta1 <- env$theta0
      
      theta1[which.] <- theta1[which.] +
        stats::runif(k, min = min.[which.], max = max.[which.])
      
      reflect_on_boundaries(
        x       = theta1,
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
