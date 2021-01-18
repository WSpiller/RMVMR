# RMVMR

## Installation

To install `RMVMR` directly from the GitHub repository, first make sure you have the `remotes` package installed:

```r
install.packages("remotes")
```

Then the `RMVMR` package can be installed using:
```r
remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```

To update the package just run the following command again.
```r
remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
``` 

## Description

We have written the `RMVMR` R package to perform radial multivariable Mendelian randomization analyses, including heterogeneity
statistics for assessing instrument strength and validity. The package accommodates any number of exposures less than 6,
and is currently includes a range of functions for estimating causal effects, as well as assessing conditional instrument strength and pleiotropic bias.



## Citation

The corresponding citation for the paper is forthcoming:

## License




This project is licensed under GNU GPL v2.
