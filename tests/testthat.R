library(testthat)
library(amcmc)

Sys.setenv("R_TESTS" = "")
test_check("amcmc")

