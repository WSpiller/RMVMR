# ivw_rmvmr

Fits a radial IVW multivariable Mendelian randomization model using
first order weights.

## Usage

``` r
ivw_rmvmr(r_input, summary = TRUE)
```

## Arguments

- r_input:

  A formatted data frame using the
  [`format_rmvmr`](https://wspiller.github.io/RMVMR/reference/format_rmvmr.md)
  function or an object of class `MRMVInput` from
  [`MendelianRandomization::mr_mvinput`](https://rdrr.io/pkg/MendelianRandomization/man/mr_mvinput.html)

- summary:

  A logical argument (`TRUE` or `FALSE`) indicating whether a summary of
  results should be presented (default= `TRUE`).

## Value

An dataframe containing MVMR results, including estimated coefficients,
their standard errors, t-statistics, and corresponding (two-sided)
p-values.

## References

Spiller, W., et al., Estimating and visualising multivariable Mendelian
randomization analyses within a radial framework. Forthcoming.

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
# Example using format_rmvmr formatted data
f.data <- format_rmvmr(
    BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
    BYG = rawdat_rmvmr$sbp_beta,
    seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
    seBYG = rawdat_rmvmr$sbp_se,
    RSID = rawdat_rmvmr$snp)
ivw_rmvmr(f.data, TRUE)
#> 
#> Radial Multivariable MR
#> 
#>               Estimate Std. Error    t value  Pr(>|t|)
#> exposure1 -0.021845534 0.01417258 -1.5413941 0.1254465
#> exposure2  0.003735376 0.01033780  0.3613319 0.7183884
#> exposure3  0.025570592 0.01601916  1.5962501 0.1126561
#> 
#> Residual standard error: 2.197 on 142 degrees of freedom
#> 
#> 

# Example using MRMVInput formatted data from the
#  MendelianRandomization package
if (require("MendelianRandomization", quietly = TRUE)) {
bx <- as.matrix(rawdat_rmvmr[,c("ldl_beta", "hdl_beta", "tg_beta")])
bxse <- as.matrix(rawdat_rmvmr[,c("ldl_se", "hdl_se", "tg_se")])
dat <- MendelianRandomization::mr_mvinput(bx = bx,
                                          bxse = bxse,
                                          by = rawdat_rmvmr$sbp_beta,
                                          byse = rawdat_rmvmr$sbp_se,
                                          snps = rawdat_rmvmr$snp)
ivw_rmvmr(r_input = dat, summary = TRUE)
}
#> 
#> Radial Multivariable MR
#> 
#>               Estimate Std. Error    t value  Pr(>|t|)
#> exposure1 -0.021845534 0.01417258 -1.5413941 0.1254465
#> exposure2  0.003735376 0.01033780  0.3613319 0.7183884
#> exposure3  0.025570592 0.01601916  1.5962501 0.1126561
#> 
#> Residual standard error: 2.197 on 142 degrees of freedom
#> 
#> 
```
