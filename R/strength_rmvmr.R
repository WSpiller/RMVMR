#' strength_rmvmr
#'
#' Calculates Q-statistics quantifying instrument strength. Each exposure is treated as an outcome sequentially, fitting the remaining
#' exposures within a radial MVMR model. High Q-statistics indicate a high instrument strength, comparable to the Q_x statistic in conventional
#' MVMR analyses. The function outputs a list of plots, global Q-statistics, and individual Q-contributions indexed by the exposure number ordered
#' using the [`format_rmvmr`] function. Named exposures in each list refer to the remaining exposures in the strength RMVMR model.
#'
#' @param r_input A formatted data frame using the [`format_rmvmr`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput`]
#' @param gencov Calculating heterogeneity statistics using the \code{MVMR} package requires the covariance between the
#'  effect of the genetic variants on each exposure to be known. This can either be estimated from individual level data,
#'  be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. A value of 0 is used by default.
#'
#' @return An object of class \code{"S_RMVMR"} containing the following components:\describe{
#' \item{\code{plot}}{A list containing plots for RMVMR analyses regressing each exposure sequentially upon remaining exposures in the \code{r_input} object. Plots are indexed by the exposure number serving as the outcome for the RMVMR analysis}
#' \item{\code{qstat}}{A list containing global Q-statistics for RMVMR analyses regressing each exposure sequentially upon remaining exposures in the \code{r_input} object. Indexing follows that of \code{plots} and p-values for global heterogeneity are provided}
#' \item{\code{qall}}{A list containing the individual Q-statistics and data for RMVMR analyses regressing each exposure sequentially upon remaining exposures in the \code{r_input} object. Indexing follows that of \code{plots}}
#' }
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Spiller, W., et al., Estimating and visualising multivariable Mendelian randomization analyses within a radial framework. Forthcoming.
#' @export
#' @examples
#' f.data <- format_rmvmr(
#'     BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
#'     BYG = rawdat_rmvmr$sbp_beta,
#'     seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
#'     seBYG = rawdat_rmvmr$sbp_se,
#'     RSID = rawdat_rmvmr$snp)
#' output <- strength_rmvmr(f.data)
#'
#' # The following shows the strength plot and Q statistics for exposure 2,
#' # regressing exposure 2 upon exposures 1 and 3 (which are labeled exposure 1
#' # and exposure 2 based on ordering in the RMVMR model).
#'
#' output$plot[[2]]
#' output$qstat[[2]]
strength_rmvmr <- function(r_input, gencov=0){

  # convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_rmvmr_format(r_input)
  }

  # Perform check that r_input has been formatted using format_rmvmr function
  if(!("rmvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "rmvmr_format", please resave the object with the output of format_rmvmr().')
  }

  if(!is.list(gencov) && gencov == 0) {
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }

  invisible(utils::capture.output(MVMR_S <- MVMR::strength_mvmr(r_input,gencov)))

  exp.number<-length(names(r_input)[-c(1,2,3)])/2

  plots <- vector('list', exp.number)
  Qs <- vector('list', exp.number)
  Qall <- vector('list', exp.number)

  for(i in 1:exp.number){

    if(exp.number == 2){

      tdat<-RadialMR::format_radial(r_input[,(4):(3+exp.number)][-i],r_input[,(3+i)],r_input[,(4+exp.number):(3+exp.number+exp.number)][-i],
                         r_input[,(3+exp.number+i)],r_input[,1])

      A <- RadialMR::ivw_radial(tdat,0.05/nrow(tdat),1,0.0001,FALSE)

      plots[[i]] <- local({
        i <- i
        p1 <- RadialMR::plot_radial(A)
      })

      Qs[[i]] <- local({
        i <- i
        qst <- A$qstatistic
      })

      Qall[[i]] <- local({
        i <- i
        qll <- A$data
      })

    }else{

      tdat<-format_rmvmr(r_input[,(4):(3+exp.number)][-i],r_input[,(3+i)],r_input[,(4+exp.number):(3+exp.number+exp.number)][-i],
                         r_input[,(3+exp.number+i)],r_input[,1])

      A<-ivw_rmvmr(tdat, FALSE)

      G<-pleiotropy_rmvmr(tdat,A)

      plots[[i]] <- local({
        i <- i
        p1 <- plot_rmvmr(tdat,A)$p2
      })

      Qs[[i]] <- local({
        i <- i
        qst <- G$gq
      })

      Qall[[i]] <- local({
        i <- i
        qll <- G$qdat
      })

    }

    }

  multi_return <- function() {
    Out_list <- list("plot" = plots,"qstat"= Qs,"qall"= Qall, "f"=MVMR_S)
    class(Out_list)<-"S_RMVMR"

    return(Out_list)
  }

  OUT<-multi_return()

}
