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
#' @details The variance covariance algorithm was described in Haario, Sksman and
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
  Mean_t,
  Mean_t_prev,
  t.,
  eps = 0, Sd = 1, Ik = diag(length(X_t))
) {
  
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
  
  (Mean_t_prev * t. + X_t)/ (t. + 1)
  
}

