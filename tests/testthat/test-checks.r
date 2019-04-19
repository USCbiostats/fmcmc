context("Checks and errors")

# ------------------------------------------------------------------------------
test_that("Error passing arguents", {
  fun0 <- function(x) {}
  fun1 <- function(x, a) {}
  fun2 <- function(x, a, b) {}
  
  expect_error(MCMC(fun0, initial = 1, nsteps = 10000, a=1), "not present in -fun")
  expect_error(MCMC(fun1, initial = 1, nsteps = 10000), "has extra arguments")
  expect_error(MCMC(fun2, initial = 1, nsteps = 10000, a=1), "requires more arguments")
  
})

test_that("Checking errors", {
  
  ll <- function(p) {
    ans <- sum(log(dnorm(y - p[1] - x*p[2], sd = p[3])))
    if (!is.finite(ans))
      return(-.Machine$double.xmax*1e-10)
    ans
  }
  
  expect_error(MCMC(ll, initial = 1, nsteps = 2000, autostop = "1"), "must be a number")
  expect_error(MCMC(ll, initial = 1, nsteps = 2000, autostop = c(1,1)), "must be of length")
})

# ------------------------------------------------------------------------------
test_that("Check of initial values", {
  
  init <- c(.4, .1)
  expect_equivalent(check_initial(init, 1), matrix(init, nrow=1))
  expect_equivalent(
    suppressWarnings(check_initial(init, 2)), matrix(init, nrow=2, ncol=2, byrow=TRUE)
    )
  
  init <- matrix(1:9, ncol=3)
  expect_equivalent(check_initial(init, 3), init)
  ## Not run: 
  
  expect_error(check_initial(init, 2), "number of rows") # Returns an error
  
})
