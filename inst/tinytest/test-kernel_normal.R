# kernel_normal ----------------------------------------------------------------

# Wild
lb <- -1
ub <- .5
k  <- kernel_normal(mu = 10, scale = 1.5)
N  <- 1e4

ans <- numeric(N)

set.seed(1)
ans <- local({
  theta0 <- 0
  theta1 <- 0
  nsteps <- N
  f <- function(p) dunif(p, log=TRUE, min = lb, max = ub)
  for (i in 1:N) {
    ans[i] <- k$proposal(environment())
  }
  
  ans
})

expect_equal(mean(ans), 10, .025)
expect_equal(sd(ans), 1.5, .025)

# Reflective 1
lb <- -1
ub <- .5
k  <- kernel_normal_reflective(mu = 10, scale = 5, lb = -.5, ub = 2)
N  <- 2e3

ans <- numeric(N)

set.seed(1)
ans <- local({
  theta0 <- 0
  theta1 <- 0
  nsteps <- N
  for (i in 1:N) {
    ans[i] <- k$proposal(environment())
  }
  
  ans
})

expect_equal(range(ans), c(-.5, 2), .025)

# Reflective 2
lb <- c(-1, 0)
ub <- c(10, 1)
k  <- kernel_normal_reflective(mu = 0, scale = 5, lb = lb, ub = ub)
N  <- 2e3

ans <- matrix(nrow = N, ncol = 2)

set.seed(1)
ans <- local({
  theta0 <- c(0,0)
  theta1 <- theta0
  nsteps <- N
  for (i in 1:N) {
    ans[i,] <- k$proposal(environment())
  }
  
  ans
})

expect_equal(range(ans[,2]), c(lb[2], ub[2]), tol = .025)
expect_equal(range(ans[,1]), c(lb[1], ub[1]), tol = .025)

# Reflective 1by1
lb <- c(-1, 0)
ub <- c(10, 1)
k  <- kernel_normal_reflective(mu = 0, scale = 5, lb = lb, ub = ub, scheme = "ordered")
N  <- 2e3

ans <- matrix(nrow = N, ncol = 2)

set.seed(1)
ans <- local({
  theta0 <- c(0,0)
  theta1 <- theta0
  nsteps <- N
  for (i in 1:N) {
    ans[i,] <- k$proposal(environment())
  }
  
  ans
})

even <- which((1:N %% 2) == 0)
odd  <- setdiff(1:N, even)

expect_true(all(ans[even,1] == 0))
expect_true(all(ans[odd,2] == 0))
expect_equal(range(ans[,2]), c(lb[2], ub[2]), tol = .025)
expect_equal(range(ans[,1]), c(lb[1], ub[1]), tol = .025)

# Expecting output
expect_output(k, "object of class fmcmc_kernel")
