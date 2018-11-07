context("Checks and errors")

# ------------------------------------------------------------------------------
test_that("Error passing arguents", {
  fun0 <- function(x) {}
  fun1 <- function(x, a) {}
  fun2 <- function(x, a, b) {}
  
  expect_error(MCMC(fun0, initial = 1, nbatch = 10000, a=1), "not present in -fun")
  expect_error(MCMC(fun1, initial = 1, nbatch = 10000), "has extra arguments")
  expect_error(MCMC(fun2, initial = 1, nbatch = 10000, a=1), "requires more arguments")
})

# ------------------------------------------------------------------------------
test_that("Check of initial values", {
  
  init <- c(.4, .1)
  expect_equivalent(check_initial(init, 1), matrix(init, nrow=1))
  expect_equivalent(
    suppressWarnings(check_initial(init, 2)), matrix(init, nrow=2, ncol=2)
    )
  
  init <- matrix(1:9, ncol=3)
  expect_equivalent(check_initial(init, 3), init)
  ## Not run: 
  
  expect_error(check_initial(init, 2), "number of rows") # Returns an error
  
})
