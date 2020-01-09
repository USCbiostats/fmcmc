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
