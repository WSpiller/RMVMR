# RMVMR

<!-- badges: start -->
[![R-CMD-check](https://github.com/WSpiller/RMVMR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WSpiller/RMVMR/actions/workflows/R-CMD-check.yaml)
[![r-universe](https://mrcieu.r-universe.dev/badges/RMVMR)](https://mrcieu.r-universe.dev/RMVMR)
<!-- badges: end -->

## Installation

RMVMR can be installed from the [MRCIEU R-Universe](https://mrcieu.r-universe.dev/) with

```r
install.packages("RMVMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
```

To install `RMVMR` directly from the GitHub repository, first make sure you have the `remotes` package installed:

```r
install.packages("remotes")
```

Then the `RMVMR` package can be installed using:

```r
remotes::install_github("WSpiller/RMVMR")
```

To update the package just run the following command again.

```r
remotes::install_github("WSpiller/RMVMR")
``` 

## Description

We have written the `RMVMR` R package to perform radial multivariable Mendelian randomization analyses, including heterogeneity
statistics for assessing instrument strength and validity. The package accommodates any number of exposures less than 6,
and is currently includes a range of functions for estimating causal effects, as well as assessing conditional instrument strength and pleiotropic bias.

## Citation

Spiller, W. Bowden, J and Sanderson E. "Estimating and visualising multivariable Mendelian randomization analyses within a radial framework" MedRxiv, 2023, doi: https://doi.org/10.1101/2023.04.04.23288134

## License

This project is licensed under GNU GPL v2.
