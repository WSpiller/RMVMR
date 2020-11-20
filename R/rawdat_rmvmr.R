#' Raw multivariable MR summary data using lipid fractions as exposures and systolic blood pressure as an outcome.
#'
#' A dataset containing summary data on 145 genetic variants associated with either
#' low-density lipoprotein (LDL), high-density lipoprotein (HDL), or triglycerides. Data includes variant rsid numbers,
#' associations with each lipid fraction, the associations between genetic variants and systolic blood pressure (SBP),
#' and corresponding standard errors.
#' 
#' rawdat_rmvmr
#'
#' @format A data frame with 145 rows and 9 variables. Specifically this includes the following information:
#' \describe{
#'  \item{\code{snp}} {The identification number for each variant}
#'  \item{\code{ldl_beta}} {The association estimate for the genetic variant obtained by regressing LDL-C upon the genetic variant}
#'  \item{\code{hdl_beta}} {The association estimate obtained by regressing HDL-C upon the genetic variant}
#'  \item{\code{tg_beta}} {The association estimate obtained by regressing triglycerides upon the genetic variant}
#'  \item{\code{sbp_beta}} {The association estimate for SBP obtained by regressing SBP upon the genetic variant}
#'  \item{\code{ldl_se}} {The standard error corresponding to association estimate `ldl_beta`}
#'  \item{\code{hdl_se}} {The standard error corresponding to association estimate `hdl_beta`}
#'  \item{\code{tg_se}} {The standard error corresponding to association estimate `tg_beta`}
#'  \item{\code{sbp_se}} {The standard error corresponding to association estimate `sbp_beta`}
#'  
#' }
#' 
#' @source 
#' \itemize{
#' \item \url{http://www.mrbase.org/}
#' \item \url{https://www.nature.com/articles/ng.2797}
#' \item \url{https://www.nature.com/articles/ng.3768}
#' }
#'@author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#'
#' @examples
#'
#' head(rawdat_rmvmr)

"rawdat_rmvmr"
