# format_rmvmr

Reads in summary data. Checks and organises columns for use in
calculating multivariable Mendelian Randomization analyses. Where
variant IDs are not provided, a vector is generated for variant
identification.

## Usage

``` r
format_rmvmr(BXGs, BYG, seBXGs, seBYG, RSID)
```

## Arguments

- BXGs:

  A matrix containing beta-coefficient values for genetic associations
  with the each exposure. Columns should indicate exposure number, with
  rows representing estimates for a given genetic variant.

- BYG:

  A numeric vector of beta-coefficient values for genetic associations
  with the outcome.

- seBXGs:

  A matrix containing standard errors corresponding to the matrix of
  beta-coefficients `BXGs`.

- seBYG:

  A numeric vector of standard errors corresponding to the
  beta-coefficients `BYG`.

- RSID:

  A vector of names for genetic variants included in the analysis. If
  variant IDs are not provided (`RSID="NULL"`), a vector of ID numbers
  will be generated.

## Value

A formatted data frame with additional classes `rmvmr_format` and
`mvmr_format`

## References

Spiller, W., et al., Estimating and visualising multivariable Mendelian
randomization analyses within a radial framework. Forthcoming.

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
f.data <- format_rmvmr(
 BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
 BYG = rawdat_rmvmr$sbp_beta,
 seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
 seBYG = rawdat_rmvmr$sbp_se,
 RSID = rawdat_rmvmr$snp)
names(f.data)
#> [1] "SNP"      "betaYG"   "sebetaYG" "betaX1"   "betaX2"   "betaX3"   "sebetaX1"
#> [8] "sebetaX2" "sebetaX3"
class(f.data)
#> [1] "data.frame"   "rmvmr_format" "mvmr_format" 
```
