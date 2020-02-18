# kernel_adapt -----------------------------------------------------------------

# Better than regular kernel
set.seed(1313)
x <- rgamma(5e2, shape = .1, rate = 2)
f <- function(p) sum(dgamma(x, shape = p[1], rate = p[2], log = TRUE))

set.seed(1)
k0 <-  kernel_normal_reflective(lb = 0)
ans0 <- suppressWarnings(MCMC(
  c(.5, .5), f, nsteps = 5000, kernel = k0, burnin = 1000,
  nchains = 2L
))
set.seed(1)
k1 <- kernel_adapt(lb = 0)
ans1 <- suppressWarnings(MCMC(
  c(.5, .5), f, nsteps = 5000, kernel = k1, burnin = 1000,
  nchains = 2L
))

mean_ans0 <- summary(ans0)$statistics[,"Mean"]
mean_ans1 <- summary(ans1)$statistics[,"Mean"]

expect_true(dist(rbind(mean_ans0, c(.1, 2))) < 1)
expect_true(dist(rbind(mean_ans1, c(.1, 2))) < 1)

expect_lt(
  coda::gelman.diag(ans0)$mpsrf,
  coda::gelman.diag(ans1)$mpsrf
)

# Tests for recursive variance, co-variance ------------------------------------
set.seed(1231)
n <- 3
X <- matrix(rnorm(n*4), ncol = 4)

# These two should be equal
ans0 <- mean_recursive(
  X_t         = X[1,],
  Mean_t_prev = colMeans(X[-1,]),
  t.          = n - 1
)
ans1 <- colMeans(X)
expect_equivalent(ans0, ans1)

# These two should be equal
ans0 <- cov_recursive(
  X_t         = X[1, ], 
  Cov_t       = cov(X[-1,]), 
  Mean_t      = colMeans(X),
  Mean_t_prev = colMeans(X[-1, ]),
  t           = n-1
)
ans1 <- cov(X)
expect_equivalent(ans0, ans1, tol = 1e-10)
