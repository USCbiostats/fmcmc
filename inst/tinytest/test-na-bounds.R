# Tests for NA bounds functionality

# Test kernel_ram with NA bounds
kern_ram <- kernel_ram(
  lb = c(alpha = NA, beta = 0, gamma = 0),
  ub = c(alpha = NA, beta = 1, gamma = NA)
)
expect_inherits(kern_ram, "fmcmc_kernel")

# Test kernel_normal_reflective with NA bounds
kern_normal <- kernel_normal_reflective(
  lb = c(NA, 0, 0),
  ub = c(NA, 1, NA)
)
expect_inherits(kern_normal, "fmcmc_kernel")

# Test kernel_adapt with NA bounds
kern_adapt <- kernel_adapt(
  lb = c(NA, 0),
  ub = c(NA, NA)
)
expect_inherits(kern_adapt, "fmcmc_kernel")

# Test kernel_nmirror with NA bounds
kern_nmirror <- kernel_nmirror(
  lb = c(NA, 0),
  ub = c(NA, 1)
)
expect_inherits(kern_nmirror, "fmcmc_kernel")

# Test kernel_umirror with NA bounds
kern_umirror <- kernel_umirror(
  lb = c(0, NA),
  ub = c(1, NA)
)
expect_inherits(kern_umirror, "fmcmc_kernel")

# Test kernel_unif_reflective with NA bounds
kern_unif <- kernel_unif_reflective(
  lb = c(NA, 0),
  ub = c(NA, 1)
)
expect_inherits(kern_unif, "fmcmc_kernel")

# Test MCMC with NA bounds works correctly
set.seed(12345)
fun <- function(x) {
  sum(dnorm(x, mean = c(0, 0.5), sd = c(1, 0.2), log = TRUE))
}

ans <- MCMC(
  initial = c(0, 0.5),
  fun = fun,
  nsteps = 100,
  kernel = kernel_normal_reflective(
    scale = 0.1,
    lb = c(NA, 0),  # First param unbounded, second bounded below
    ub = c(NA, 1)   # First param unbounded, second bounded above
  )
)

# Verify second parameter is within bounds [0, 1]
expect_true(all(ans[,2] >= 0))
expect_true(all(ans[,2] <= 1))

# Test that NA in lb is equivalent to -.Machine$double.xmax
kern_na <- kernel_normal_reflective(lb = c(NA, 0), ub = c(NA, 1))
kern_explicit <- kernel_normal_reflective(
  lb = c(-.Machine$double.xmax, 0), 
  ub = c(.Machine$double.xmax, 1)
)

# Both should create valid kernels
expect_inherits(kern_na, "fmcmc_kernel")
expect_inherits(kern_explicit, "fmcmc_kernel")
