# Internal helpers shared by pleiotropy_rmvmr() and plot_rmvmr(). These are not
# exported and take an already validated rmvmr_format data frame.

# Perform a univariate radial MR analysis for each exposure using only the SNPs
# whose first-stage F-statistic exceeds 10, and assemble the per-SNP results
# into a single data frame with a Group factor labelling the reference exposure.
rmvmr_univariate_radial <- function(r_input, exp.number) {
  f.vec <- matrix(0L, nrow = length(r_input[, 1]), ncol = exp.number)

  for (i in 1:exp.number) {
    f.vec[, i] <- as.integer(r_input[, 3 + i]^2 / r_input[, 3 + exp.number + i]^2 >= 10)
  }

  p.list <- vector("list", exp.number)

  for (i in 1:exp.number) {
    #Format data for univariate MR using significant SNPs for each exposure
    Xsub <- r_input[f.vec[, i] == 1, ]
    Xrad.dat <- RadialMR::format_radial(
      Xsub[, 3 + i],
      Xsub[, 2],
      Xsub[, 3 + exp.number + i],
      Xsub[, 3],
      Xsub[, 1]
    )
    Xfit <- RadialMR::ivw_radial(
      Xrad.dat,
      0.05 / nrow(Xrad.dat),
      1,
      0.0001,
      FALSE
    )
    Xdat <- data.frame(Xfit[5])
    Xdat$Group <- i
    names(Xdat) <- c("SNP", "Wj", "BetaWj", "Qj", "Qj_Chi", "Outliers", "Group")
    p.list[[i]] <- Xdat
  }

  p.dat <- do.call(rbind, p.list)

  p.dat[, 7] <- as.factor(p.dat[, 7])
  for (i in 1:exp.number) {
    levels(p.dat[, 7])[i] <- paste0("Exposure_", i, collapse = "")
  }

  p.dat
}

# Apply the MVMR correction to the univariate radial results in p.dat using the
# fitted IVW coefficients (rmvmr, the coef matrix from ivw_rmvmr). Returns the
# global corrected Q-statistics and the per-SNP corrected data frame.
rmvmr_correction <- function(r_input, rmvmr, exp.number, p.dat) {
  Ratio_list <- vector("list", exp.number)

  for (i in 1:exp.number) {
    tdat <- r_input[
      r_input$SNP %in% p.dat[p.dat$Group == levels(p.dat$Group)[i], ]$SNP,
    ]
    Ratio_temp <- tdat[, 2] / tdat[, (3 + i)]

    for (j in 1:exp.number) {
      if (i != j) {
        Ratio_temp <- Ratio_temp -
          ((tdat[, (3 + j)] * rmvmr[j, 1]) / tdat[, (3 + i)])
      }
    }
    Ratio_list[[i]] <- Ratio_temp
  }
  Ratios <- unlist(Ratio_list)

  #QJ calculations
  p.dat$coratios <- Ratios

  Qj_list <- vector("list", exp.number)

  for (i in 1:exp.number) {
    tdat <- p.dat[p.dat$Group == levels(p.dat$Group)[i], ]
    Qj_list[[i]] <- tdat$Wj^2 * (tdat$coratios - rmvmr[i, 1])^2
  }
  Qjvec <- unlist(Qj_list)

  p.dat$Qjcor <- Qjvec

  #Define matrix for recording total Q statistics.
  Qj_out <- matrix(0L, nrow = exp.number, ncol = 2)

  for (i in 1:exp.number) {
    Qj_out[i, 1] <- sum(p.dat[p.dat$Group == levels(p.dat$Group)[i], ]$Qjcor)
    Qj_out[i, 2] <- stats::pchisq(
      Qj_out[i, 1],
      nrow(p.dat[p.dat$Group == levels(p.dat$Group)[i], ]) - exp.number,
      lower.tail = FALSE
    )
  }

  TotalQs <- data.frame(Qj_out)
  names(TotalQs) <- c("q_statistic", "p_value")
  row.names(TotalQs) <- levels(p.dat$Group)

  indqj <- stats::pchisq(p.dat$Qjcor, 1, lower.tail = FALSE)

  p.dat$corQjchi <- indqj

  out_data <- p.dat[, c(1, 2, 8, 9, 10, 7)]

  names(out_data) <- c(
    "snp",
    "wj",
    "corrected_beta",
    "qj",
    "qj_p",
    "ref_exposure"
  )

  list("gq" = TotalQs, "qdat" = out_data)
}
