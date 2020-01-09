#' @export
#' @rdname kernels
#' @param eta A function that receives the MCMC environment. This is to calculate
#' the scalig factor for the adaptation.
#' @param arate Numeric scalar. Objective acceptance rate.
#' @section Kernels: 
#' 
#' The `kernel_ram` Implements Vihola (2012)'s Robust Adaptive Metropolis.
kernel_ram <- function(
  mu     = 0,
  eta    = function(i, k) min(c(1.0, i^(-2.0/3.0) * k)),
  arate  = 0.234,
  freq   = 1L,
  warmup = 0L,
  Sigma  = NULL,
  eps    = 1e-4,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  fixed  = FALSE
) {
  
  k          <- NULL
  which.     <- NULL
  Ik         <- NULL
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k     <<- length(env$theta0)
        
        mu     <<- check_dimensions(mu, k)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        which. <<- which(!fixed)
        
        k      <<- sum(!fixed)
        Ik     <<- diag(k)
        
        if (is.null(Sigma))
          Sigma <<- Ik * eps
        
        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        
      }
      
      # Making proposal
      U      <- stats::rnorm(k)
      theta1 <- env$theta1
      theta1[which.] <- env$theta0[which.] + (Sigma %*% U)[,1L]
      
      # Updating the scheme
      if (env$i > warmup && !(env$i %% freq)) {
        
        # Computing
        a_n <- min(1, exp(env$f(theta1) - env$f0))
        if (!is.finite(a_n))
          a_n <- 0.0
        
        Sigma <<- t(chol(Sigma %*% (
          Ik + eta(env$i, k)*(a_n - arate) * U %*% t(U) /
            norm(rbind(U), "2")
        ) %*% t(Sigma)))
        
        
      }
      
      # Reflecting
      reflect_on_boundaries(theta1, lb, ub, which.)
      
    },
    mu         = mu,
    eta        = eta, 
    arate      = arate,
    freq       = freq,
    warmup     = warmup,
    Sigma      = Sigma,
    eps        = eps,
    lb         = lb,
    ub         = ub,
    fixed      = fixed,
    k          = k,
    which.     = which.,
    Ik         = Ik
  )
  
}
