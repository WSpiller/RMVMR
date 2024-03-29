% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ivw_rmvmr.R
\name{ivw_rmvmr}
\alias{ivw_rmvmr}
\title{ivw_rmvmr}
\usage{
ivw_rmvmr(r_input, summary = TRUE)
}
\arguments{
\item{r_input}{A formatted data frame using the \code{\link{format_rmvmr}} function or an object of class \code{MRMVInput} from \code{\link[MendelianRandomization:mr_mvinput]{MendelianRandomization::mr_mvinput}}}

\item{summary}{A logical argument (\code{TRUE} or \code{FALSE}) indicating whether a summary of results should be presented (default= \code{TRUE}).}
}
\value{
An dataframe containing MVMR results, including estimated coefficients, their standard errors, t-statistics, and corresponding (two-sided) p-values.
}
\description{
Fits a radial IVW multivariable Mendelian randomization model using first order weights.
}
\examples{
# Example using format_rmvmr formatted data
f.data <- format_rmvmr(
    BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
    BYG = rawdat_rmvmr$sbp_beta,
    seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
    seBYG = rawdat_rmvmr$sbp_se,
    RSID = rawdat_rmvmr$snp)
ivw_rmvmr(f.data, TRUE)

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
}
\references{
Spiller, W., et al., Estimating and visualising multivariable Mendelian randomization analyses within a radial framework. Forthcoming.
}
\author{
Wes Spiller; Eleanor Sanderson; Jack Bowden.
}
