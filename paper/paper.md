---
title: 'amcmc: A flexible MCMC framework'
authors:
- affiliation: 1
  name: George G Vega Yon
  orcid: 0000-0002-3171-0844
- affiliation: 1
  name: Paul Marjoram
  orcid: 0000-0003-0824-7449
date: "4 October 2018"
output:
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- text mining
- natural language processing
affiliations:
- index: 1
  name: Department of Preventive Medicine, University of Southern California
---

# Summary

**amcmc** is an R package

The package is designed for R users needing to apply .

# Feature 1

**amcmc** blah

# Feature 2

 
# Package design

**quanteda** has been carefully designed with several key aims in mind.

_Consistency_.  **quanteda** functions and objects are named systematically such
that `corpus()`, `tokens()` and `dfm()` construct those object types, and that
`corpus_*()`, `tokens_*()` and `dfm_*()` functions return a modified version of these
objects. Naming consistency applies also to the extensive built-in data objects in the
package, whose names always start with `data_*` followed by object types. This
not only gives the users a clear overview of the package, but also makes the
package more reliable for other packages that depend on it.

_Accessibility_.  **quanteda** contains extensive manual pages structured around
the naming rules. Furthermore, there are references, package vignettes, examples,
and tutorials on the website at https://quanteda.io. 
These materials help beginner users understand how to use these functions for
basic operations and expert users how to combine the functions for advanced text
processing and analysis.

_Performance_.  **quanteda**'s performance is enhanced by token hashing and
parallel computation implemented in C++, permitting large and fast text analysis
even on computers with relatively limited resources (such as laptop computers).
Built to use sparse data structures, **quanteda** can efficiently performs
complex textual data analyses, such as computing distances, calculating feature
discrimination statistics (keyness), or model fitting, even for large
document-feature matrices.

_Transparency and reproducibility_.  **quanteda** is designed to facilitate rigorous,
transparent, and reproducible scientific analysis of text. Being open-source
software, its source code can be scrutinized and corrected by other experts. Its
functions are designed to encourage a reproducible workflow by linking successive 
processing tasks in a clear, readable manner.

_Compatibility with other packages_.  For analysis not provided by built-in
functions, users can move **quanteda** objects seamlessly to other packages,
such as the
**stm** package for structural topic models [@STM] or word embedding packages
like **text2vec** [@text2vec].  **quanteda** also works well with companion
packages such as **spacyr** [@spacyr], an R wrapper to the spaCy Python library
[@spacy2], and **readtext** [@readtext], a package for converting and importing
text files into R.

# Funding and Support

This work is supported by the National Cancer Institute (NCI), Grant Number 5P01CA196569-02

# References
