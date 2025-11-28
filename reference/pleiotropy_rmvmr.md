# pleiotropy_rmvmr

Generates Q-statistics quantifying the degree of heterogeneity in
univariate Radial MR analyses applying a correction using the output
from
[`ivw_rmvmr`](https://wspiller.github.io/RMVMR/reference/ivw_rmvmr.md).
The function returns two data frames. The first data frame includes the
global Q-statistic for each exposure after applying a correction, as
well as a corresponding p-value. The second data frame contains the
individual Q-statistic for each SNP in the corrected univariate
analyses, relative to the exposure given in column `exposure`.

## Usage

``` r
pleiotropy_rmvmr(r_input, rmvmr)
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

## Value

An object of class `"RMVMR_Q"` containing the following components:

- `gq`:

  A data frame containing the global Q-statistic and p-value after
  applying a correction for each exposure

- `qdat`:

  A data frame containing the individual Q-statistic and p-value for
  each SNP after applying a correction for each exposure

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
rmvmr_output <- ivw_rmvmr(f.data, FALSE)
q_object <- pleiotropy_rmvmr(f.data, rmvmr_output)
q_object$gq
#>            q_statistic   p_value
#> Exposure_1    39.00391 0.9999987
#> Exposure_2    32.72136 0.9999999
#> Exposure_3    33.84794 0.9998315
head(q_object$qdat)
#>          snp        wj corrected_beta         qj       qj_p ref_exposure
#> 1 rs10019888  9.638623    0.182234939 0.40143746 0.52634782   Exposure_1
#> 2 rs10468017 51.780930    0.001079676 0.02721426 0.86896961   Exposure_1
#> 3  rs1047891 12.076698    0.490835137 3.17425712 0.07480722   Exposure_1
#> 4 rs10761762  9.195911   -0.413504104 1.41062004 0.23495349   Exposure_1
#> 5 rs11045163 10.340424   -0.176298188 0.24667725 0.61942451   Exposure_1
#> 6 rs11065987 10.561270   -0.802588890 6.43772951 0.01117214   Exposure_1
```
