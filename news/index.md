# Changelog

## RMVMR 0.4.5

- Fixed a bug in
  [`strength_rmvmr()`](https://wspiller.github.io/RMVMR/reference/strength_rmvmr.md)
  where supplying `gencov` as a matrix or vector (rather than the scalar
  default) errored before
  [`MVMR::strength_mvmr()`](https://wspiller.github.io/MVMR/reference/strength_mvmr.html)
  was called.
- [`format_rmvmr()`](https://wspiller.github.io/RMVMR/reference/format_rmvmr.md)
  now generates placeholder variant IDs when `RSID = NULL`, matching the
  documentation. The `RSID` argument description previously showed a
  typo (`RSID="NULL"`).
- [`plot_rmvmr()`](https://wspiller.github.io/RMVMR/reference/plot_rmvmr.md)
  no longer recomputes the univariate radial analyses when called with
  the default `cordat = NULL`, roughly halving the number of
  [`RadialMR::ivw_radial()`](https://wspiller.github.io/RadialMR/reference/ivw_radial.html)
  calls. Results are unchanged.
- [`ivw_rmvmr()`](https://wspiller.github.io/RMVMR/reference/ivw_rmvmr.md)
  now fits the radial IVW model explicitly rather than via variables
  left over from the orientation loop. The estimate is invariant to the
  choice of orientation, so this is fitted directly as the canonical
  model; the returned coefficients (estimates, standard errors, t-values
  and p-values), residual standard error and degrees of freedom are
  unchanged.
- [`ivw_rmvmr()`](https://wspiller.github.io/RMVMR/reference/ivw_rmvmr.md)
  no longer returns the undocumented `data` element. This element held
  intermediate orientation data frames that were not used elsewhere in
  the package. The estimated coefficients remain unchanged and are still
  available in the `coef` element.

## RMVMR 0.4.4

- Bump roxygen2 to 8.0.0 and add a package level helpfile.

## RMVMR 0.4.3

- Implement some code optimizations.

## RMVMR 0.4.2

- RMVMR now requires MVMR version 0.4.3 or later.

## RMVMR 0.4.1

- Fixed a bug in
  [`plot_rmvmr()`](https://wspiller.github.io/RMVMR/reference/plot_rmvmr.md)
  in the
  [`ggplot2::scale_x_continuous()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  calls. RMVMR should now run under version 4 of ggplot2, which is due
  to be released soon.

## RMVMR 0.4

- RMVMR now required R 4.1.0 or newer. This is because of the
  requirement of a dependency package of ggplot2, the scales package,
  now having this requirement (which is the new tidyverse requirement on
  release of R 4.5.0).
