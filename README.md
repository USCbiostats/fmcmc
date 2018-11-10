amcmc: A flexible MCMC estimation framework
================

[![Travis-CI Build
Status](https://travis-ci.org/USCbiostats/amcmc.svg?branch=master)](https://travis-ci.org/USCbiostats/amcmc)
[![Build
status](https://ci.appveyor.com/api/projects/status/3x9qj7imvoijb1vf?svg=true)](https://ci.appveyor.com/project/gvegayon/amcmc)
[![Coverage
Status](https://img.shields.io/codecov/c/github/USCbiostats/amcmc/master.svg)](https://codecov.io/github/USCbiostats/amcmc?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/amcmc)](https://cran.r-project.org/package=amcmc)

Current features:

1.  Automatic stop using convergence checker.

2.  Parallel chains using `parallel::lapply`

3.  Normal Random Walk with Reflective Boundaries kernel.

# Installing

From github:

``` r
devtools::install_github("USCbiostats/amcmc")
```

# Example

Linear regression model

``` r
library(amcmc)

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

# Running the MCMC
ans <- MCMC(
  ll, X. = X, y. = y,
  initial = c(1, 1, 1),
  nsteps   = 2e4,
  nchains  = 4L,
  autostop = 1e3,
  burnin   = 1e4,
  scale    = .1,
  multicore = TRUE
  )
```

    ## Warning: A single initial point has been passed via `initial`: c(1, 1, 1).
    ## The values will be recycled.

    ## Convergence has been reached with 11000 steps (1000 final count of observations).

``` r
library(coda)

summary(ans)
```

    ## 
    ## Iterations = 10001:11000
    ## Thinning interval = 1 
    ## Number of chains = 4 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##       Mean      SD  Naive SE Time-series SE
    ## par1 3.072 0.06856 0.0010840       0.004018
    ## par2 1.991 0.06352 0.0010044       0.003638
    ## par3 2.050 0.04591 0.0007259       0.002161
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##       2.5%   25%   50%   75% 97.5%
    ## par1 2.939 3.028 3.073 3.121 3.201
    ## par2 1.870 1.948 1.992 2.033 2.118
    ## par3 1.960 2.020 2.051 2.082 2.139

``` r
plot(ans)
```

![](README_files/figure-gfm/summary-and-plot1-1.png)<!-- -->

``` r
gelman.diag(ans)
```

    ## Potential scale reduction factors:
    ## 
    ##      Point est. Upper C.I.
    ## par1       1.01       1.03
    ## par2       1.01       1.03
    ## par3       1.01       1.02
    ## 
    ## Multivariate psrf
    ## 
    ## 1.01

# Other tools

  - <https://cran.r-project.org/web/packages/mcmc/mcmc.pdf>

  - <https://cran.r-project.org/web/packages/HybridMC/HybridMC.pdf>

  - <https://cran.r-project.org/web/packages/adaptMCMC/adaptMCMC.pdf>

  - <https://cran.r-project.org/web/packages/elhmc/elhmc.pdf>

# Contributing to `amcmc`

Please note that the ‘amcmc’ project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.

# Funding

Supported by National Cancer Institute Grant \#1P01CA196596.
