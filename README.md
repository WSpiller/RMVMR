# MVMR

## Installation

To install `RMVMR` directly from the GitHub repository, first make sure you have the `remotes` package installed:

    install.packages("remotes")

Then the `RMVMR` package can be installed using:

    library(remotes)
    install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

    install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
    
    install_github("WSpiller/RadialMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
    
    
    
To update the package just run the `remotes::install_github("WSpiller/RMVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)` command again.

## Description

We have written the `RMVMR` R package to perform radial multivariable Mendelian randomization analyses, including heterogeneity
statistics for assessing instrument strength and validity. The package accommodates any number of exposures less than 6,
and is currently includes a range of functions for estimating causal effects, as well as assessing conditional instrument strength and pleiotropic bias.



## Citation

The corresponding citation for the paper is forthcoming:

## License

This project is licensed under GNU GPL v3.



