---
title: 'amcmc: A flexible MCMC framework'
authors:
- affiliation: 1
  name: George G Vega Yon
  orcid: 0000-0002-3171-0844
- affiliation: 1
  name: Paul Marjoram
  orcid: 0000-0003-0824-7449
date: "22 April 2019"
output:
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- metropolis-hastings
- mcmc
- markov chain monte carlo
- transition kernel
- automatic convergence
affiliations:
- index: 1
  name: Department of Preventive Medicine, University of Southern California
---

# Summary

Markov Chain Monte Carlo (MCMC) is used in a variaty of statistical and computational venues such as: statistical inference, Markov quadrature (also known as Monte Carlo integration), stochastic optimization, among others. The **amcmc** R [@R] package provides a flexible framework for implementing MCMC methods that use the Metropolis-Hastings algorithm [@hastings1970; @metropolis1953]. Notwithstanding the simplicity of its implementation, **amcmc** provides the following out-of-the-box features that can be valued for both practitioners of MCMC and educators: seemensly multiple-chains sampling using parallel computing, user-defined transition kernel, and automatic stop using convergence monitoring.

In the case of transition kernels, users can either use one of the transition kernels shipped with the package (e.g. the gaussian kernel and its bounded version, the gaussian kernel with reflective boundaries), or create their own without having to write the MCMC algorithm from scracth.

Likewise, users can use either use one of the convergence monitoring checking functions that are part of the package, for example: The Gelman and Rubin's [@Gelman1992], Geweke's [@Geweke1991], etc. Or build their own to be used within the framework.

While there are several other R packages that either implement MCMC algorithms or provide wrappers for implementations in other languages are several [see for example: @Sturtz2005; @Morey2009; @Stan2018; @Geyer2019; @Scheidegger2019], to the authors knowledge, the `amcmc` package is the first one on its kind providing a framework as flexible as the one described here implemented fully within R.

# Funding and Support

This work is supported by the National Cancer Institute (NCI), Grant Number 5P01CA196569-02

# References
