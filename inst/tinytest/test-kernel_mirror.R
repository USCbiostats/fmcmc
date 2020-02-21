# kernel_normal ----------------------------------------------------------------

set.seed(1313)
x <- rgamma(5e2, shape = 1, rate = 2)
f <- function(p) {
  sum(dgamma(x, shape = p[1], rate = p[2], log = TRUE))
}

# Better than regular kernel
set.seed(1)
ans0 <- suppressWarnings(MCMC(
  c(.5, .5), f, nsteps = 2000,
  kernel = kernel_normal_reflective(lb = 0), burnin = 1000,
  nchains = 2L
))

set.seed(1)
ans1 <- suppressWarnings(MCMC(
  c(.5, .5), f, nsteps = 2000,
  kernel = kernel_nmirror(lb = 0, warmup = 1000),
  burnin = 1000,
  nchains = 2L
))

set.seed(1)
ans2 <- suppressWarnings(MCMC(
  c(.5, .5), f, nsteps = 2000,
  kernel = kernel_umirror(lb = 0, warmup = 1000),
  burnin = 1000,
  nchains = 2L
))

mean_ans0 <- summary(ans0)$statistics[,"Mean"]
mean_ans1 <- summary(ans1)$statistics[,"Mean"]
mean_ans2 <- summary(ans2)$statistics[,"Mean"]

expect_true(dist(rbind(mean_ans0, c(.1, 2))) < 1.5)
expect_true(dist(rbind(mean_ans1, c(.1, 2))) < 1.5)
expect_true(dist(rbind(mean_ans2, c(.1, 2))) < 1.5)

expect_lt(
  coda::gelman.diag(ans1)$mpsrf[1],
  coda::gelman.diag(ans0)$mpsrf[1]
)

expect_lt(
  coda::gelman.diag(ans2)$mpsrf[1],
  coda::gelman.diag(ans0)$mpsrf[1]
)
