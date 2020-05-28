#' Robust Adaptive Metropolis (RAM) Transition Kernel
#' 
#' Implementation of Vihola (2012)'s Robust Adaptive Metropolis.
#' 
#' @export
#' @template lb-ub
#' @template mu-Sigma
#' @template until
#' @param eta A function that receives the MCMC environment. This is to calculate
#' the scaling factor for the adaptation.
#' @param arate Numeric scalar. Objective acceptance rate.
#' @param qfun Function. As described in Vihola (2012)'s, the `qfun` function is
#' a symmetric function used to generate random numbers.
#' @param fixed Logical scalar or vector of length `k`. Indicates which parameters
#' will be treated as fixed or not. Single values are recycled.
#' @param freq Integer scalar. Frequency of updates. How often the
#' variance-covariance matrix is updated.
#' @param constr Logical lower-diagonal square matrix of size `k`. **Not** in the
#' original paper, but rather a tweak that imposes a constraint on the `S_n`
#' matrix. If different from `NULL`, the kernel multiplates `S_n` by this
#' constraint so that zero elements are pre-imposed.
#' 
#' @details 
#' 
#' The idea is similar to that of the Adaptive Metropolis algorithm (AM implemented
#' as [kernel_adapt()] here) with the difference that it takes into account a
#' target acceptance rate.
#' 
#' The `eta` function regulates the rate of adaptation. The default implementation
#' will decrease the rate of adaptation exponentially as a function of the iteration
#' number.
#' 
#' @return An object of class [fmcmc_kernel].
#' 
#' @references 
#' Vihola, M. (2012). Robust adaptive Metropolis algorithm with coerced acceptance
#' rate. Statistics and Computing, 22(5), 997–1008.
#' \url{https://doi.org/10.1007/s11222-011-9269-5}
#' @family kernels
#' @examples 
#' # Setting the acceptance rate to 30 % and deferring the updates until
#' # after 1000 steps
#' kern <- kernel_ram(arate = .3, warmup = 1000)
kernel_ram <- function(
  mu     = 0,
  eta    = function(i, k) min(c(1.0, i^(-2.0/3.0) * k)),
  qfun   = function(k) stats::rt(k, k),
  arate  = 0.234,
  freq   = 1L,
  warmup = 0L,
  Sigma  = NULL,
  eps    = 1e-4,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  fixed  = FALSE,
  until  = Inf,
  constr = NULL
) {
  
  k       <- NULL
  which.  <- NULL
  Ik      <- NULL
  nerrors <- 0L
  
  # We create copies of this b/c each chain has its own values
  abs_iter <- 0L
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k     <<- length(env$theta0)
        
        # mu     <<- check_dimensions(mu, k)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        which. <<- which(!fixed)
        
        k      <<- sum(!fixed)
        Ik     <<- diag(k)
        
        # # We wont be reusing mu, so we need to adjust it
        # mu <<- mu[which.]
        
        # Initializing Sigma (for all chains)
        if (is.null(Sigma))
          Sigma <<- Ik * eps
        
        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        
      }
      
      # Making proposal
      U      <- qfun(k)
      theta1 <- env$theta1
      theta1[which.] <- env$theta0[which.] + (Sigma %*% U)[, 1L]
      
      # Updating the scheme
      if ((until > abs_iter) && (abs_iter > warmup) && !(env$i %% freq)) {
        
        # Computing
        a_n <- min(1, exp(env$f(theta1) - env$f0))
        if (!is.finite(a_n))
          a_n <- 0.0
        
        Sigma <<- Sigma %*% (
          Ik + eta(env$i, k) * (a_n - arate) * tcrossprod(U) /
            norm(rbind(U), "2") ^ 2.0
        ) %*% t(Sigma)
        
        Sigma_temp <<- tryCatch(t(chol(Sigma)), error = function(e) e)
        if (inherits(Sigma_temp, "error")) {
          nerrors <<- nerrors + 1L
          Sigma <<- t(chol(Matrix::nearPD(Sigma)$mat))
        } else
          Sigma <<- Sigma_temp
        
        # Applying a constraint on sigma
        if (!is.null(constr))
          Sigma <<- constr[which., , drop = FALSE][, which., drop = FALSE] * Sigma

      }
      
      # Increasing the absolute number of iteration
      abs_iter <<- abs_iter + 1L
      
      # Reflecting
      reflect_on_boundaries(theta1, lb, ub, which.)
      
    },
    mu         = mu,
    eta        = eta, 
    qfun       = qfun,
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
    Ik         = Ik,
    abs_iter   = abs_iter,
    until      = until,
    constr     = constr,
    nerrors    = nerrors
  )
  
}
