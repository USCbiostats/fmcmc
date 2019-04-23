library(testthat)
library(fmcmc)

Sys.setenv("R_TESTS" = "")
test_check("fmcmc")

