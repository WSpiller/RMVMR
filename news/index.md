# Changelog

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
