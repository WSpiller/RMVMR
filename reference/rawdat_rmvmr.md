# Raw multivariable MR summary data using lipid fractions as exposures and systolic blood pressure as an outcome.

A dataset containing summary data on 145 genetic variants associated
with either low-density lipoprotein (LDL), high-density lipoprotein
(HDL), or triglycerides. Data includes variant rsid numbers,
associations with each lipid fraction, the associations between genetic
variants and systolic blood pressure (SBP), and corresponding standard
errors.

## Usage

``` r
rawdat_rmvmr
```

## Format

A data frame with 145 rows and 9 variables. Specifically this includes
the following information:

- `snp`:

  The identification number for each variant

- `ldl_beta`:

  The association estimate for the genetic variant obtained by
  regressing LDL-C upon the genetic variant

- `hdl_beta`:

  The association estimate obtained by regressing HDL-C upon the genetic
  variant

- `tg_beta`:

  The association estimate obtained by regressing triglycerides upon the
  genetic variant

- `sbp_beta`:

  The association estimate for SBP obtained by regressing SBP upon the
  genetic variant

- `ldl_se`:

  The standard error corresponding to association estimate `ldl_beta`

- `hdl_se`:

  The standard error corresponding to association estimate `hdl_beta`

- `tg_se`:

  The standard error corresponding to association estimate `tg_beta`

- `sbp_se`:

  The standard error corresponding to association estimate `sbp_beta`

## Source

- <http://www.mrbase.org/>

- <https://www.nature.com/articles/ng.2797>

- <https://www.nature.com/articles/ng.3768>

## Details

rawdat_rmvmr

## Author

Wes Spiller; Eleanor Sanderson; Jack Bowden.

## Examples

``` r
head(rawdat_rmvmr)
#>          snp ldl_beta hdl_beta tg_beta    sbp_beta ldl_se hdl_se  tg_se
#> 1 rs10019888  -0.0270   0.0182  0.0228 -0.00426935 0.0046 0.0050 0.0045
#> 2 rs10468017   0.1179   0.0020  0.0379  0.00110389 0.0038 0.0042 0.0039
#> 3  rs1047891  -0.0269   0.0079  0.0000 -0.01317370 0.0039 0.0042 0.0038
#> 4 rs10490626   0.0081  -0.0508  0.0085 -0.00111303 0.0064 0.0069 0.0062
#> 5 rs10761762   0.0191   0.0103 -0.0270 -0.00854986 0.0034 0.0036 0.0033
#> 6 rs10832962   0.0043   0.0320  0.0109  0.00472509 0.0038 0.0040 0.0037
#>       sbp_se
#> 1 0.00280123
#> 2 0.00227690
#> 3 0.00222743
#> 4 0.00374264
#> 5 0.00207701
#> 6 0.00237133
```
