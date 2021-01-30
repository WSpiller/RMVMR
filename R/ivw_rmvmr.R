#' ivw_rmvmr
#'
#' Fits a radial IVW multivariable Mendelian randomization model using first order weights.
#'
#' @param r_input A formatted data frame using the \code{format_rmvmr} function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput`]
#' @param summary A logical argument (\code{TRUE} or \code{FALSE}) indicating whether a summary of results should be presented (default= \code{TRUE}).
#'
#' @return An dataframe containing MVMR results, including estimated coefficients, their standard errors, t-statistics, and corresponding (two-sided) p-values.
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Spiller, W., et al., Estimating and visualising multivariable Mendelian randomization analyses within a radial framework. Forthcoming.
#' @importFrom stats lm
#' @export
#' @examples
#' # Example using format_rmvmr formatted data
#' f.data <- format_rmvmr(
#'     BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
#'     BYG = rawdat_rmvmr$sbp_beta,
#'     seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
#'     seBYG = rawdat_rmvmr$sbp_se,
#'     RSID = rawdat_rmvmr$snp)
#' ivw_rmvmr(f.data, TRUE)
#'
#' # Example using MRMVInput formatted data from the
#' #  MendelianRandomization package
#' bx <- as.matrix(rawdat_rmvmr[,c("ldl_beta", "hdl_beta", "tg_beta")])
#' bxse <- as.matrix(rawdat_rmvmr[,c("ldl_se", "hdl_se", "tg_se")])
#' dat <- MendelianRandomization::mr_mvinput(bx = bx,
#'                                           bxse = bxse,
#'                                           by = rawdat_rmvmr$sbp_beta,
#'                                           byse = rawdat_rmvmr$sbp_se,
#'                                           snps = rawdat_rmvmr$snp)
#' ivw_rmvmr(r_input = dat, summary = TRUE)

# Define IVW Radial Multivariable MR function: This takes the formatted dataframe from
# the format_MVMR function as an input, and outputs a summary of effect estimates as well as formatted radial data frames
# for downstream plotting

ivw_rmvmr<-function(r_input,summary){

  # convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_rmvmr_format(r_input)
  }

  # Perform check that r_input has been formatted using format_rmvmr function
  if(!("rmvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "rmvmr_format", please resave the object with the output of format_rmvmr().')
  }

  #Determine the number of exposures included in the model
  exp.number<-length(names(r_input)[-c(1,2,3)])/2

  #Create zero matrix for weight calculations
  tm.weights<-matrix(0L, nrow = length(r_input[,1]), ncol = exp.number)

  #Create subset of exposure summary data
  exp.dat<-r_input[,4:(3+exp.number)]


  #Calculate square root weights wj
  for(i in 1:exp.number){
    tm.weights[,i] = sqrt((exp.dat[,i]^2)/r_input[,3]^2)
  }

  #Create zero matrix for ratio calculations
  tm.ratios<-matrix(0L, nrow = length(r_input[,1]), ncol = exp.number)

  #Calculate ratio estimates
  for(i in 1:exp.number){
    tm.ratios[,i] = r_input[,2]/exp.dat[,i]
  }

  #Create zero matrix for weight times ratio calculations
  tm.wr<-matrix(0L, nrow = length(r_input[,1]), ncol = exp.number)

  #Multiply weighting by ratio estimates
  for(i in 1:exp.number){
    tm.wr[,i]<-tm.weights[,i]*tm.ratios[,i]
  }

  #Create empty list for plotting data frames

  t.list <- vector(mode = "list", length = 3)

  # orientate data for X1

  for(j in 1:exp.number){

  expvec<- 1:exp.number


  tm.oriented<-matrix(0L, nrow = length(r_input[,1]), ncol = exp.number)


  for(i in 1:(exp.number)){

    if(j == i){
      tm.oriented[,i]<-tm.weights[,i]
    }else{
      tm.oriented[,i]<-exp.dat[,expvec[(i)]]*sign(exp.dat[,j])
      tm.oriented[,i]<-tm.oriented[,i] / r_input[,3]

    }

  }

  tempdat<-data.frame(tm.oriented)

  t.list[[j]] <- cbind(tm.wr[,j],tempdat)

  names(t.list[[j]])[1]<-paste0("Bwj_",j,collapse="")

  #Rename columns for ease of interpretation
  for(i in 1:exp.number){
    names(t.list[[j]])[i+1]<- paste0("wj_",i,collapse="")
    }

  }

  A_sum<-summary(lm(tm.wr[,j]~ -1 + ., tempdat))

  A<-summary(lm(tm.wr[,j]~ -1 + ., tempdat))$coef

  #Rename the regressors for ease of interpretation
  for(i in 1:exp.number){
    dimnames(A)[[1]][i]<- paste0("exposure",i,collapse="")
  }

  if(summary == T){

    # Print a few summary elements that are common to both lm and plm model summary objects
    cat("\n")
    cat("Radial Multivariable MR\n")
    cat("\n")
    print(A)
    cat("\nResidual standard error:", round(A_sum$sigma,3), "on", A_sum$df[2], "degrees of freedom")
    cat("\n")
    cat("\n")
    cat("\n")
  }

  multi_return <- function() {
    Out_list <- list("coef" = A,"data"= t.list)
    class(Out_list)<-"IVW_RMVMR"

    return(Out_list)
}

OUT<-multi_return()


}
