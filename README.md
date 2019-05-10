
# fmcmc: A friendly MCMC framework <img src="man/figures/logo.png" align="right" height="140"/>

[![Travis-CI Build
Status](https://travis-ci.org/USCbiostats/fmcmc.svg?branch=master)](https://travis-ci.org/USCbiostats/fmcmc)
[![Build
status](https://ci.appveyor.com/api/projects/status/3x9qj7imvoijb1vf?svg=true)](https://ci.appveyor.com/project/gvegayon/fmcmc)
[![Coverage
Status](https://img.shields.io/codecov/c/github/USCbiostats/fmcmc/master.svg)](https://codecov.io/github/USCbiostats/fmcmc?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/fmcmc)](https://cran.r-project.org/package=fmcmc)
[![status](http://joss.theoj.org/papers/2e86b709451443990c1c6776ebb7f756/status.svg)](http://joss.theoj.org/papers/2e86b709451443990c1c6776ebb7f756)

The `fmcmc` R package implements a flexible

Current features:

1.  Automatic stop using convergence checker.

2.  Parallel chains using `parallel`

3.  Flexible framework to specify different transition kernels.

4.  Implements the Normal (random walk) with reflective boundaries
    kernel.

# Installing

From github:

``` r
devtools::install_github("USCbiostats/fmcmc")
```

# Example

Linear regression model

``` r
library(fmcmc)

# Simulating data
set.seed(78845)
n <- 1000
X <- rnorm(n)
y <- 3 + 2*X + rnorm(n, sd = 2)

# Loglikelihood function
ll <- function(par, X., y.) {
  
  ans <- sum(log(dnorm((y. - (par[1] + X.*par[2]))/par[3])/par[3]))
  
  if (!is.finite(ans))
    return(-Inf)
  
  ans
}
```

``` r
# Running the MCMC
ans <- MCMC(
  ll, X. = X, y. = y,
  initial = c(1, 1, 1),
  # It allows to specify different kernels
  kernel   = kernel_reflective(
    lb    = .00001,
    ub    = 100,
    scale = .1
  ),
  # As well as convergence checkers
  conv_checker = convergence_gelman(threshold = 1.05),
  nsteps   = 2e5,
  thin     = 10,
  autostop = 2e3,
  burnin   = 1e4,
  # As well as parallel chains
  nchains  = 4L,
  multicore = TRUE
  )
```

    ## Warning: While using multiple chains, a single initial point has been
    ## passed via `initial`: c(1, 1, 1). The values will be recycled. Ideally you
    ## would want to start each chain from different locations.

    ## Convergence has been reached with 12000 steps. (200 final count of observations).

``` r
library(coda)

summary(ans)
```

    ## 
    ## Iterations = 10010:12000
    ## Thinning interval = 10 
    ## Number of chains = 4 
    ## Sample size per chain = 200 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##       Mean      SD Naive SE Time-series SE
    ## par1 3.069 0.06742 0.002384       0.002910
    ## par2 1.990 0.06542 0.002313       0.002350
    ## par3 2.049 0.04779 0.001690       0.002064
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##       2.5%   25%   50%   75% 97.5%
    ## par1 2.930 3.025 3.069 3.116 3.197
    ## par2 1.862 1.942 1.990 2.037 2.110
    ## par3 1.963 2.018 2.048 2.081 2.141

``` r
plot(ans)
```

![](man/figures/summary-and-plot1-1.png)<!-- -->

``` r
gelman.diag(ans)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## par1       1.00       1.01
    ## par2       1.01       1.02
    ## par3       1.00       1.00
    ## 
    ## Multivariate psrf
    ## 
    ## 1.01

# Other tools

  - <https://cran.r-project.org/web/packages/mcmc/mcmc.pdf>

  - <https://cran.r-project.org/web/packages/HybridMC/HybridMC.pdf>

  - <https://cran.r-project.org/web/packages/adaptMCMC/adaptMCMC.pdf>

  - <https://cran.r-project.org/web/packages/elhmc/elhmc.pdf>

# Contributing to `fmcmc`

Please note that the ‘fmcmc’ project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.

# Funding

Supported by National Cancer Institute Grant \#1P01CA196596. cd
