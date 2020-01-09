#' @export 
#' @rdname kernels
#' @param bw Integer scalar. The bandwidth, is the number of observations to
#' include in the computation of the variance-covariance matrix.
#' @param freq Integer scalar. Frequency of updates. How often the
#' variance-covariance matrix is updated.
#' @param warmup Integer scalar. The number of iterations that the algorithm has
#' to wait before starting to do the updates.
#' @param Sigma The variance-covariance matrix. By default this will be an
#' identity matrix during the warmup period.
#' @param eps Double scalar. Default size of the initial step (see details).
#' @param Sd Overall scale for the algorithm. By default, the variance-covariance
#' is scaled to \eqn{2.4^2/d}, with \eqn{d} the number of dimensions.
#' @section Kernels:
#' `kernel_adapt` Implements the adaptive Metropolis (AM) algorithm of Haario
#' et al. (2001). If the value of bw is greater than zero, then the algorithm
#' folds back AP, a  previous version which is known to have ergodicity problems.
#' 
#' The parameter `eps` has two functions. The first one is to set the initial
#' scale for the multivariate normal kernel, which is replaced after `warmup`
#' steps with the actual variance-covariance computed by the main algorithm.
#' The second usage is in the equation that ensures that the variance-covariance
#' is greater than zero, this is, the \eqn{\varepsilon}{epsilon} parameter in the
#' original paper.
#' 
#' 
kernel_adapt <- function(
  mu     = 0,
  bw     = 0L,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  freq   = 50L,
  warmup = 500L,
  Sigma  = NULL,
  Sd     = NULL,
  eps    = 1e-4,
  fixed  = FALSE
) {
  
  k      <- NULL
  sd     <- NULL
  Ik     <- NULL
  which. <- NULL
  
  if (bw > 0L && bw > warmup)
    stop("The `warmup` parameter must be greater than `bw`.", call. = FALSE)
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k     <<- length(env$theta0)
        Ik    <<- diag(k) * eps
        
        if (is.null(Sigma))
          Sigma <<- Ik
        
        mu     <<- check_dimensions(mu, k)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        which. <<- which(!fixed)
        
        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        
        if (is.null(Sd))
          Sd <<- 5.76 / length(which.) # 2.38^2
      }
      
      # Updating the scheme
      if (env$i > warmup && env$i > 2 && !(env$i %% freq)) {
        
        ran <- if (bw <= 0L) 1L:(env$i - 1L) 
        else (env$i - bw + 1L):(env$i - 1L)
        
        Sigma <<- Sd * (stats::cov(env$ans[ran, , drop = FALSE]) + Ik)
      }
      
      # Making the proposal
      theta1 <- env$theta0
      theta1[which.] <- env$theta0[which.] +
        MASS::mvrnorm(mu = mu[which.], Sigma = Sigma[which., , drop = FALSE][,which. , drop = FALSE])
      
      reflect_on_boundaries(theta1, lb, ub, which.)
      
    },
    mu     = mu,
    bw     = bw,
    lb     = lb,
    ub     = ub,
    freq   = freq,
    warmup = warmup,
    Sigma  = Sigma,
    Sd     = Sd,
    eps    = eps,
    fixed  = fixed,
    k      = k,
    Ik     = Ik,
    which. = which.
  )
  
}
