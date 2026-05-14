# plot_rmvmr

Generates two radial multivariable Mendelian randomization (MVMR) plots.
The first plot shows the estimated direct effect for each exposure
obtained by fitting a radial MVMR model. Each data point shows the
square root weighting for each SNP on the x-axis, and product of the
ratio estimate and square root weighting for each SNP on the y-axis.
These values are obtained by performing a univariate radial MR analysis
for each exposure using the SNPs displayed, specifically through use of
the
[`RadialMR::ivw_radial`](https://wspiller.github.io/RadialMR/reference/ivw_radial.html)
function. Only SNPs strongly associated with the corresponding exposure
are used, such that their first stage F-statistic is greater than 10.
The second plot applies a correction to each ratio estimate. In both
plots, the distance of each observation from the corresponding
regression line is proportional to the contribution of that SNP towards
global heterogeneity.

## Usage

``` r
plot_rmvmr(r_input, rmvmr, cordat = NULL)
```

## Arguments

- r_input:

  A formatted data frame using the
  [`format_rmvmr`](https://wspiller.github.io/RMVMR/reference/format_rmvmr.md)
  function or an object of class `MRMVInput` from
  [`MendelianRandomization::mr_mvinput`](https://rdrr.io/pkg/MendelianRandomization/man/mr_mvinput.html)

- rmvmr:

  An object containing the output from the
  [`ivw_rmvmr`](https://wspiller.github.io/RMVMR/reference/ivw_rmvmr.md)
  function of class `IVW_RMVMR`.

- cordat:

  Optional. A pre-computed object from
  [`pleiotropy_rmvmr`](https://wspiller.github.io/RMVMR/reference/pleiotropy_rmvmr.md).
  If `NULL` (default),
  [`pleiotropy_rmvmr`](https://wspiller.github.io/RMVMR/reference/pleiotropy_rmvmr.md)
  is called internally.

## Value

An object of class `"RMVMR_plot"` containing the following components:

- `p1`:

  A radial MVMR plot without correction

- `p2`:

  A radial MVMR plot with correction

## References

Spiller, W., et al., Estimating and visualising multivariable Mendelian
randomization analyses within a radial framework. Forthcoming.

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
# \donttest{
f.data <- format_rmvmr(
    BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
    BYG = rawdat_rmvmr$sbp_beta,
    seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
    seBYG = rawdat_rmvmr$sbp_se,
    RSID = rawdat_rmvmr$snp)
rmvmr_output <- ivw_rmvmr(f.data, FALSE)
plot_object <- plot_rmvmr(f.data, rmvmr_output)
plot_object$p1

plot_object$p2

# }
```
