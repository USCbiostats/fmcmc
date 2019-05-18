context("Convergence checkers work") 

test_that("Multiple chains mcmc behaves", {
  
  set.seed(71)
  n <- 1000
  x <- rnorm(n)
  s <- 2
  y <- 1 + .2*x + rnorm(n, sd = s)
  
  ll <- function(p) {
    ans <- sum(log(dnorm(y - p[1] - x*p[2], sd = p[3])))
    if (!is.finite(ans))
      return(-.Machine$double.xmax*1e-10)
    ans
  }
  
  # Gwedeke autostop convergence
  ans0 <- suppressWarnings(MCMC(
    ll, initial = rep(5, 3), nsteps = 2000,
    kernel = kernel_reflective(
      lb = 0, ub = 10, mu = 0, scale = .1
      ), conv_checker = convergence_auto(500L),
    thin = 10
    ))
  
  ans0 <- MCMC(
    ll, initial = tail(ans0, 0L), nsteps = 2000,
    kernel = kernel_reflective(
      lb = 0, ub = 10, mu = 0, scale = .1
    ), conv_checker = NULL,
    thin = 10, burnin = 0L
  )
  
  # Gelman autostop convergence
  ans1 <- suppressWarnings(MCMC(
    ll, initial = rep(5, 3), nsteps = 2000,
    kernel = kernel_reflective(
      lb = 0, ub = 10, mu = 0, scale = .1
    ),
    thin         = 10,
    conv_checker = convergence_gelman(500L),
    nchains      = 2
  ))
  
  ans1 <- MCMC(
    ll,
    initial = colMeans(do.call(rbind, tail(ans1, 0))),
    nsteps = 2000,
    kernel = kernel_reflective(
      lb = 0, ub = 10, mu = 0, scale = .1
    ), conv_checker = NULL,
    thin = 10, burnin = 0L
  )
  
  # Heidel autostop convergence
  ans2 <- suppressWarnings(MCMC(
    ll, initial = rep(5, 3), nsteps = 2000,
    kernel = kernel_reflective(
      lb = 0, ub = 10, mu = 0, scale = .1
    ), 
    thin = 10, conv_checker = convergence_heildel(500L)
  ))
  
  ans2 <- MCMC(
    ll, initial = tail(ans2, 0), nsteps = 2000,
    kernel = kernel_reflective(
      lb = 0, ub = 10, mu = 0, scale = .1
    ), conv_checker = NULL,
    thin = 10, burnin = 0L
  )
  
  # Should reach an estimates close to the model parameters
  sol <- c(1, .2, 2)
  expect_equivalent(colMeans(ans0), sol, tol = 0.2)
  expect_equivalent(colMeans(ans1), sol, tol = 0.2)
  expect_equivalent(colMeans(ans2), sol, tol = 0.2)
  
})

