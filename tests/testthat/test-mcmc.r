context("MCMC")

test_that("Error checks", {
  
  f <- function(i) {}
  expect_error(
    MCMC(f, initial = 1, nsteps = 100, burnin = 100),
    "burnin"
  )
  
  expect_error(
    MCMC(f, initial = 1, nsteps = 100, burnin = 0, thin = -1),
    "thin"
  )
  
  expect_error(
    MCMC(f, initial = 1, nsteps = 100, burnin = 0, thin = 100),
    "thin"
  )
  
  f <- function(i) {NaN}
  expect_error(
    MCMC(f, initial = 1, nsteps = 100, burnin = 10),
    "undefined"
  )
  
})

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
  set.seed(111)
  ans0 <- suppressWarnings(
    MCMC(fun, initial = 1, nsteps = 5e3, burnin = 500,
         kernel = kernel_normal(scale = 1))
  )
  expect_equal(mean(ans0), mean(D), tolerance = 0.05, scale = 1)
  
})

# ------------------------------------------------------------------------------
test_that("Reasonable values after changing the scale", {
  
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
  set.seed(111)
  ans0 <- suppressWarnings(
    MCMC(fun, initial = 1, nsteps = 5e3, burnin = 500,
         kernel = kernel_normal(scale = 2))
  )
  expect_equal(mean(ans0), mean(D), tolerance = 0.05, scale = 1)
  
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
  ans0 <- suppressWarnings(
    MCMC(fun, initial = 1, nsteps = 5e3, burnin = 500,
         kernel = kernel_normal(scale = 1),
         D=D, nchains=2)
  )
  
  ans1 <- suppressWarnings(
    MCMC(fun, initial = 1, nsteps = 5e3, burnin = 500,
         kernel = kernel_normal(scale = 1),
         D=D, nchains=2, multicore = TRUE)
  )
  
  expect_equal(sapply(ans0, mean), rep(mean(D),2), tolerance = 0.1, scale = 1)
  expect_equal(sapply(ans1, mean), rep(mean(D),2), tolerance = 0.1, scale = 1)
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
    MCMC(fun, initial = 1, nsteps = 200, burnin = 0,
         kernel = kernel_normal(scale = 1),
         nchains=2, D=D)
  )
  
  set.seed(1)
  ans1 <- suppressWarnings(
    MCMC(fun, initial = 1, nsteps = 200, burnin = 0,
         kernel = kernel_normal(scale = 1),
         nchains=2, D=D)
  )
  
  expect_equal(ans0, ans1)
  
})


# ------------------------------------------------------------------------------
test_that("Fixed parameters", {
  
  # Simulating data
  set.seed(981)
  
  D <- rnorm(1000, 0, 2)
  
  # Preparing function
  fun <- function(x) {
    res <- log(dnorm(D, x[1], x[2]))
    if (any(is.infinite(res) | is.nan(res)))
      return(.Machine$double.xmin)
    sum(res)
  }
  
  # Running the algorithm and checking expectation
  ans <- suppressWarnings(
    MCMC(fun, initial = c(1, 2), nsteps = 5e3, burnin = 500,
         kernel = kernel_reflective(
           ub = 3, lb = -3, scale = 1, fixed = c(FALSE, TRUE)
         ))
  )
  expect_true(all(ans[,2] == 2))
  expect_equivalent(colMeans(ans), c(0, 2), tolerance = 0.1, scale = 1)
  
})
