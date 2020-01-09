# plan_update_sequence ---------------------------------------------------------

set.seed(1313)
x <- rgamma(5e2, shape = .1, rate = 2)
f <- function(p) sum(dgamma(x, shape = p[1], rate = p[2], log = TRUE))

# Wild
lb <- -1
ub <- .5
k  <- kernel_unif(min. = lb, max. = ub)
N  <- 5e3

set.seed(1)
expect_error({
  suppressWarnings(MCMC(
    c(.5, .5), f, nsteps = 2000,
    kernel = kernel_normal_reflective(
      lb = 0, fixed = c(TRUE, FALSE), scheme = c(1, 2)), burnin = 1000,
    nchains = 2L
    ))},
  "not be fixed")

set.seed(1)
expect_error({
  suppressWarnings(MCMC(
    c(.5, .5), f, nsteps = 2000,
    kernel = kernel_normal_reflective(
      lb = 0, fixed = c(TRUE, TRUE)), burnin = 1000,
    nchains = 2L
  ))},
  "cannot be zero")

ans0 <- fmcmc:::plan_update_sequence(2, 20000, c(FALSE, FALSE), "random")
expect_equal(colMeans(ans0), c(.5,.5), tol = .05)

# Testing errors ---------------------------------------------------------------

# Wrong length of parameters
k <- kernel_unif(fixed = c(TRUE, TRUE))
set.seed(1)
expect_error({
  ans <- local({
    theta0 <- c(0, 0, 0)
    theta1 <- theta0
    nsteps <- N
    for (i in 1:N) {
      ans[i,] <- k$proposal(environment())
    }
    
    ans
  })
}, "length of")

# Wrong length of scheme
k <- kernel_unif(fixed = c(TRUE, TRUE), scheme = c(1,2,3))
set.seed(1)
expect_error({
  ans <- local({
    theta0 <- c(0, 0)
    theta1 <- theta0
    nsteps <- N
    for (i in 1:N) {
      ans[i,] <- k$proposal(environment())
    }
    
    ans
  })
}, "same length as")


