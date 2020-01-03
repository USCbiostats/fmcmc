
# kernel_unif baseline ---------------------------------------------------------

# Wild
lb <- -1
ub <- .5
k  <- kernel_unif(min. = lb, max. = ub)
N  <- 5e3

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

expect_equal(mean(ans), (lb + ub)/2, .025)

# Fixed parameters
k  <- kernel_unif(min. = lb, max. = ub, fixed = c(FALSE, TRUE, FALSE))
N  <- 5e3

ans <- matrix(nrow = N, ncol = 3)

set.seed(1)
ans <- local({
  theta0 <- c(0, 0, 0)
  theta1 <- c(0, 0, 0)
  nsteps <- N
  f <- function(p) sum(dunif(p, log=TRUE, min = lb, max = ub))
  for (i in 1:N) {
    ans[i,] <- k$proposal(environment())
  }
  
  ans
})

expect_true(all(ans[,2] == 0))
expect_equal(colMeans(ans), c((lb + ub)/2, 0, (lb + ub)/2), .025)

# kernel_unif_reflective -------------------------------------------------------

# Wild
lb <- -1
ub <- .5
k  <- kernel_unif_reflective(min. = lb, max. = ub, lb = lb / 2, ub = ub / 2)
N  <- 5e3

ans <- numeric(N)

set.seed(1)
ans <- local({
  theta0 <- 0
  theta1 <- 0
  nsteps <- N
  f <- function(p) dunif(p, log=TRUE, min = lb / 2, max = ub / 2)
  for (i in 1:N) {
    ans[i] <- k$proposal(environment())
  }
  
  ans
})

expect_equal(mean(ans), (lb + ub)/4, .025)

# Fixed parameters
k  <- kernel_unif_reflective(
  min. = lb, max. = ub, lb = lb / 2, ub = ub / 2,
  fixed = c(TRUE, FALSE, TRUE))
N  <- 5e3

ans <- matrix(nrow = N, ncol = 3)

set.seed(1)
ans <- local({
  theta0 <- c(0, 0, 0)
  theta1 <- c(0, 0, 0)
  nsteps <- N
  f <- function(p) sum(dunif(p, log=TRUE, min = lb / 2, max = ub / 2))
  for (i in 1:N) {
    ans[i,] <- k$proposal(environment())
  }
  
  ans
})

expect_true(all(ans[,-2] == 0))
expect_equal(colMeans(ans), c(0, (lb + ub)/4, 0), .025)

# kernel_unif scheme -----------------------------------------------------------

# Wild
lb <- -1
ub <- .5
k  <- kernel_unif(min. = lb, max. = ub, scheme = "ordered")
N  <- 5e3

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

expect_equal(mean(ans), (lb + ub)/2, .025)

# Fixed parameters
k  <- kernel_unif(min. = lb, max. = ub, fixed = c(FALSE, TRUE, FALSE),
                       scheme = "ordered")

ans <- matrix(nrow = N, ncol = 3)

set.seed(1)
ans <- local({
  theta0 <- c(0, 0, 0)
  theta1 <- c(0, 0, 0)
  nsteps <- N
  f <- function(p) sum(dunif(p, log=TRUE, min = lb, max = ub))
  for (i in 1:N) {
    ans[i,] <- k$proposal(environment())
  }
  
  ans
})

expect_true(all(ans[,2] == 0))
even <- which((1:N %% 2) == 0)
odd  <- setdiff(1:N, even)
expect_equal(mean(ans[odd,1]), (lb + ub)/2, .025)
expect_equal(mean(ans[even,3]), (lb + ub)/2, .025)

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

# kernel_adapt -----------------------------------------------------------------

# Better than regular kernel
set.seed(1313)
x <- rgamma(5e2, shape = .1, rate = 2)
f <- function(p) sum(dgamma(x, shape = p[1], rate = p[2], log = TRUE))

set.seed(1)
ans0 <- suppressWarnings(MCMC(
  c(.5, .5), f, nsteps = 2000, kernel = kernel_normal_reflective(lb = 0), burnin = 1000,
  nchains = 2L
  ))
set.seed(1)
ans1 <- suppressWarnings(MCMC(
  c(.5, .5), f, nsteps = 2000, kernel = kernel_adapt(lb = 0), burnin = 1000,
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


