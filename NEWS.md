# RMVMR 0.4.7

* `pleiotropy_rmvmr()`, `plot_rmvmr()` and `strength_rmvmr()` no longer pass on the warning from RadialMR 1.2.5 and later that the fixed-effect exact confidence interval of `RadialMR::ivw_radial()` is empty. RMVMR does not use that interval, so the warning, which was repeated for each exposure, was not relevant. Results are unchanged.
* Now requires RadialMR 1.2.5 or later.

# RMVMR 0.4.6

* Fixed a bug in the corrected Q-statistic calculation used by `pleiotropy_rmvmr()` and `plot_rmvmr()`: the square-root radial weight `Wj` was not squared, so SNP contributions were weighted by the square-root weight rather than the full inverse-variance weight. Q-statistics, their p-values and outlier detection now differ from previous versions (thanks @nvitt, #8).

# RMVMR 0.4.5

* Fixed a bug in `strength_rmvmr()` where supplying `gencov` as a matrix or vector (rather than the scalar default) errored before `MVMR::strength_mvmr()` was called.
* `format_rmvmr()` now generates placeholder variant IDs when `RSID = NULL`, matching the documentation. The `RSID` argument description previously showed a typo (`RSID="NULL"`).
* `plot_rmvmr()` no longer recomputes the univariate radial analyses when called with the default `cordat = NULL`, roughly halving the number of `RadialMR::ivw_radial()` calls. Results are unchanged.
* `ivw_rmvmr()` now fits the radial IVW model explicitly rather than via variables left over from the orientation loop. The estimate is invariant to the choice of orientation, so this is fitted directly as the canonical model; the returned coefficients (estimates, standard errors, t-values and p-values), residual standard error and degrees of freedom are unchanged.
* `ivw_rmvmr()` no longer returns the undocumented `data` element. This element held intermediate orientation data frames that were not used elsewhere in the package. The estimated coefficients remain unchanged and are still available in the `coef` element.

# RMVMR 0.4.4

* Bump roxygen2 to 8.0.0 and add a package level helpfile.

# RMVMR 0.4.3

* Implement some code optimizations.

# RMVMR 0.4.2

* RMVMR now requires MVMR version 0.4.3 or later.

# RMVMR 0.4.1

* Fixed a bug in `plot_rmvmr()` in the `ggplot2::scale_x_continuous()` calls. RMVMR should now run under version 4 of ggplot2, which is due to be released soon.

# RMVMR 0.4

* RMVMR now required R 4.1.0 or newer. This is because of the requirement of a dependency package of ggplot2, the scales package, now having this requirement (which is the new tidyverse requirement on release of R 4.5.0).
