set.seed(12314)
dat <- rnorm(1e3, mean = 2, sd = 1.5)
logpost <- function(x) {
  sum(dnorm(dat, mean = x[1], sd = x[2], log = TRUE))
}

nsteps <- 1e3

# Simple chains ----------------------------------------------------------------
ans0 <- MCMC(
  c(1,1), logpost, nsteps = nsteps, nchains = 1,
  kernel = kernel_normal_reflective(lb = c(-1, 0))
  )
ans1 <- MCMC(
  c(1,1), logpost, nsteps = nsteps, nchains = 1,
  kernel = kernel_normal_reflective(lb = c(-1, 0))
  )

ans <- append_chains(ans0, ans1)
expect_true(all(rownames(ans) == 1:(nsteps*2)))

# MCMC lists -------------------------------------------------------------------
ans0 <- MCMC(
  matrix(runif(4), nrow = 2), logpost, nsteps = nsteps, nchains = 2,
  kernel = kernel_normal_reflective(lb = c(-1, 0))
)
ans1 <- MCMC(
  matrix(runif(4), nrow = 2), logpost, nsteps = nsteps, nchains = 2,
  kernel = kernel_normal_reflective(lb = c(-1, 0))
)

ans <- append_chains(ans0, ans1)
ans <- sapply(lapply(ans, rownames), function(i) all(i == 1:(nsteps*2)))
expect_equal(ans, c(TRUE, TRUE))

# Different number of chains
expect_error(append_chains(ans0[1], ans1), "same number of chains")

# Error in thin parameter ------------------------------------------------------
ans0 <- MCMC(
  c(1,1), logpost, nsteps = nsteps, nchains = 1,
  kernel = kernel_normal_reflective(lb = c(-1, 0)), thin = 10
)
ans1 <- MCMC(
  c(1,1), logpost, nsteps = nsteps, nchains = 1,
  kernel = kernel_normal_reflective(lb = c(-1, 0)), thin = 15
)

expect_error(append_chains(ans0, ans1), "have the same `thin")
expect_error(append_chains(1, ans0), "method available")

logpost2 <- function(x) logpost(c(x, 1.5))

ans2 <- MCMC(
  c(1), logpost2, nsteps = nsteps, nchains = 1,
  kernel = kernel_normal_reflective(lb = c(-1)), thin = 15
)

expect_error(append_chains(ans1, ans2), "same number of parameters")
