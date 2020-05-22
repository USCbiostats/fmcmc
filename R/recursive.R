#' Recursive algorithms for computing variance and mean
#' 
#' These algorithms are used in [kernel_adapt()] to simplify variance-covariance
#' recalculation at every step of the algorithm.
#' 
#' @param X_t Last value of the sample
#' @param Cov_t Covariance in t
#' @param Mean_t,Mean_t_prev Vectors of averages in time `t` and `t-1` respectively.
#' @param t. Sample size up to `t-1`.
#' @param Sd,eps,Ik See [kernel_adapt()].
#' @export
#' @details The variance covariance algorithm was described in Haario, Saksman and
#' Tamminen (2002).
#' 
#' @references 
#' Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis algorithm.
#' Bernoulli, 7(2), 223–242.
#' \url{https://projecteuclid.org/euclid.bj/1080222083}
#' @examples 
#' # Generating random data (only four points to see the difference)
#' set.seed(1231)
#' n <- 3
#' X <- matrix(rnorm(n*4), ncol = 4)
#' 
#' # These two should be equal
#' mean_recursive(
#'   X_t         = X[1,],
#'   Mean_t_prev = colMeans(X[-1,]),
#'   t.          = n - 1
#' )
#' colMeans(X)
#' 
#' # These two should be equal
#' cov_recursive(
#'   X_t         = X[1, ], 
#'   Cov_t       = cov(X[-1,]), 
#'   Mean_t      = colMeans(X),
#'   Mean_t_prev = colMeans(X[-1, ]),
#'   t           = n-1
#' )
#' cov(X)
#' 
#' # Speed example -------------------------------------------------------------
#' set.seed(13155511)
#' X <- matrix(rnorm(1e3*100), ncol = 100)
#' 
#' ans0 <- cov(X[-1,])
#' t0 <- system.time({
#'   ans1 <- cov(X)
#' })
#' 
#' t1 <- system.time(ans2 <- cov_recursive(
#'   X[1, ], ans0,
#'   Mean_t      = colMeans(X),
#'   Mean_t_prev = colMeans(X[-1,]),
#'   t. = 1e3 - 1
#' ))
#' 
#' # Comparing accuracy and speed
#' range(ans1 - ans2)
#' t0/t1
#' 
cov_recursive <- function(
  X_t,
  Cov_t,
  Mean_t_prev,
  t.,
  Mean_t = NULL,
  eps    = 0,
  Sd     = 1,
  Ik     = diag(ncol(Cov_t))
) {
  
  # Computing the current average
  if (is.null(Mean_t))
    Mean_t <- mean_recursive(X_t = X_t, Mean_t_prev = Mean_t_prev, t. = t.)
  
  if (is.matrix(X_t)) {
    
    # First the mean, which is easier
    ans <- array(dim = c(ncol(X_t), rev(dim(X_t))))
    for (i in 1:nrow(X_t)) {
      if (i == 1L)
        ans[, , i] <- cov_recursive(
          X_t         = X_t[i, ],
          Cov_t       = Cov_t,
          Mean_t      = Mean_t[i, ],
          Mean_t_prev = Mean_t_prev,
          t.          = t.,
          eps         = eps,
          Sd          = Sd,
          Ik          = Ik
        )
      else 
        ans[, , i] <- cov_recursive(
          X_t         = X_t[i, ],
          Cov_t       = ans[, , i - 1L],
          Mean_t      = Mean_t[i, ],
          Mean_t_prev = Mean_t[i - 1L, ],
          t.          = t. + i - 1,
          eps         = eps,
          Sd          = Sd,
          Ik          = Ik
        )
        
    }

    return(ans)
    
  }
  
  (t. - 1)/t. * Cov_t + 
    Sd/t. * (
      t. * tcrossprod(Mean_t_prev) -
        (t. + 1) * tcrossprod(Mean_t) +
        tcrossprod(X_t) + 
        eps * Ik
    )
  
}

#' @export
#' @rdname cov_recursive
mean_recursive <- function(X_t, Mean_t_prev, t.) {
  
  if (!is.matrix(X_t))
    return((Mean_t_prev * t. + X_t)/ (t. + 1))
  else {
    ans <- matrix(nrow = nrow(X_t), ncol = ncol(X_t))
    for (i in 1:nrow(ans)) {
      if (i == 1)
        ans[i,] <- (Mean_t_prev * t. + X_t[i, ])/ (t. + 1)
      else
        ans[i,] <- (ans[i - 1L, ] * (t. + i - 1L) + X_t[i, ])/ (t. + i)
    }
    return(ans)
  }
  
}

