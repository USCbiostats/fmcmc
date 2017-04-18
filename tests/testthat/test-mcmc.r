context("MCMC")

# ------------------------------------------------------------------------------
test_that("Reasonable values", {
  
  # Simulating data
  set.seed(981)
  
  D <- rnorm(1000, 0)
  
  # Preparing function
  fun <- function(x) {
    res <- log(dnorm(D, x))
    if (any(is.infinite(res) | is.nan(res)))
      return(.Machine$double.xmin)
    sum(res)
  }
  
  # Running the algorithm and checking expectation
  ans <- suppressWarnings(
    MCMC(fun, 1, 5e3, burnin = 500, ub = 3, lb = -3, scale = 1)
  )
  expect_equal(mean(ans), mean(D), tolerance = 0.1, scale = 1)
})

# ------------------------------------------------------------------------------
test_that("Multiple chains", {
  # Simulating data
  set.seed(981)
  
  D <- rnorm(1000, 0)
  
  # Preparing function
  fun <- function(x, D) {
    res <- log(dnorm(D, x))
    if (any(is.infinite(res) | is.nan(res)))
      return(.Machine$double.xmin)
    sum(res)
  }
  
  # Running the algorithm and checking expectation
  ans <- suppressWarnings(
    MCMC(fun, 1, 5e3, burnin = 500, ub = 3, lb = -3, scale = 1, D=D, nchains=2)
  )
  expect_equal(sapply(ans, mean), rep(mean(D),2), tolerance = 0.1, scale = 1)
})

# ------------------------------------------------------------------------------
test_that("Error passing arguents", {
  fun0 <- function(x) {}
  fun1 <- function(x, a) {}
  fun2 <- function(x, a, b) {}
  
  expect_error(MCMC(fun0, 1, a=1), "not present in -fun")
  expect_error(MCMC(fun1, 1), "has extra arguments")
  expect_error(MCMC(fun2, 1, a=1), "requires more arguments")
})

# ------------------------------------------------------------------------------
test_that("Repeating the chains in parallel", {
  # Simulating data
  set.seed(981)
  
  D <- rnorm(500, 0)
  
  # Preparing function
  fun <- function(x, D) {
    res <- log(dnorm(D, x))
    if (any(is.infinite(res) | is.nan(res)))
      return(.Machine$double.xmin)
    sum(res)
  }
  
  # Running the algorithm and checking expectation
  set.seed(1)
  ans0 <- suppressWarnings(
    MCMC(fun, 1, 200, burnin = 0, ub = 3, lb = -3, scale = 1, nchains=2, D=D)
  )
  
  set.seed(1)
  ans1 <- suppressWarnings(
    MCMC(fun, 1, 200, burnin = 0, ub = 3, lb = -3, scale = 1, nchains=2, D=D)
  )
  
  expect_equal(ans0, ans1)
  
})

