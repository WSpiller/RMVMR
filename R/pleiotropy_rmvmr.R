#' pleiotropy_rmvmr
#'
#' Generates Q-statistics quantifying the degree of heterogeneity in univariate Radial MR analyses applying a correction using the
#' output from \code{ivw_rmvmr}. The function returns two data frames. The first data frame includes the global Q-statistic for each exposure after applying
#' a correction, as well as a corresponding p-value. The second data frame contains the individual Q-statistic for each SNP in the corrected univariate
#' analyses, relative to the exposure given in column \code{exposure}.
#'
#' @param r_input A formatted data frame using the \code{format_rmvmr} function.
#' @param rmvmr An object containing the output from the \code{ivw_rmvmr} function of class \code{IVW_RMVMR}.
#'
#' @return An object of class \code{"RMVMR_Q"} containing the following components:\describe{
#' \item{\code{gq}}{A data frame containing the global Q-statistic and p-value after applying a correction for each exposure}
#' \item{\code{qdat}}{A data frame containing the individual Q-statistic and p-value for each SNP after applying a correction for each exposure}
#'}
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#'@references Spiller, W., et al., Estimating and visualising multivariable Mendelian randomization analyses within a radial framework. Forthcoming.
#' @importFrom stats pchisq
#' @export
#' @examples
#'
#' f.data <- format_rmvmr(
#'     BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
#'     BYG = rawdat_rmvmr$sbp_beta,
#'     seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
#'     seBYG = rawdat_rmvmr$sbp_se,
#'     RSID = rawdat_rmvmr$snp)
#'     
#' rmvmr_output<-ivw_rmvmr(f.data,F)
#'
#' q_object<-pleiotropy_rmvmr(f.data,rmvmr_output)
#' 
#' q_object$gq
#' q_object$qdat


pleiotropy_rmvmr<-function(r_input,rmvmr){
  
  #Load RadialMR package
  library(RadialMR)
  
  # Extract MVMR estimates
  rmvmr<-rmvmr$coef
  
  #Define number of exposures included in MVMR model
  exp.number<-length(names(r_input)[-c(1,2,3)])/2
  
  #Define matrix for identifying IV1 satisfying variants using F>10.
  f.vec<-matrix(0L, nrow = length(r_input[,1]), ncol = exp.number)
  
  for(i in 1:exp.number){
    f.vec[,i]<- r_input[,3+i]^2/r_input[,3 + exp.number + i]^2
    for(j in 1:length(r_input[,1])){
      if(f.vec[j,i] < 10){
        f.vec[j,i]<-0
      }else{
        f.vec[j,i]<-1
      }
    }
  }
  
  #Define null variable for univariate MR data
  Xlist<-NULL
  
  #Obtain univariate MR data for each exposure
  for(i in 1:exp.number){
    
    Xsub<-r_input[f.vec[,i] == 1,]
    Xrad.dat<-format_radial(Xsub[,3+i],Xsub[,2],Xsub[,3 + exp.number + i],Xsub[,3],Xsub[,1])
    X.res<-ivw_radial(Xrad.dat,0.05/nrow(Xrad.dat),1,0.0001,F)
    if(is.null(Xlist)){
      Xlist<-X.res
    }else{
      Xlist<-append(Xlist,X.res)
    }
  }
  
  #Create combined data frame of univariate values
  
  p.dat<-NULL
  for(i in 1:exp.number){
    Xdat<-data.frame(Xlist[5 + ((i-1)* 13)])
    Xdat$Group<-rep(i,nrow(Xdat))
    names(Xdat)<-c("SNP","Wj","BetaWj","Qj","Qj_Chi","Outliers","Group")
    
    if(is.null(p.dat)){
      p.dat<-Xdat
    }else{
      p.dat<-rbind(p.dat,Xdat)
    }
    
  }
  
  p.dat[,7]<-as.factor(p.dat[,7])
  for(i in 1:exp.number){
    levels(p.dat[,7])[i] <- paste0("Exposure_",i,collapse="")
  }
  
  ##AL
  
  
  Ratios<-NULL
  
  for(i in 1:exp.number){
    
    tdat<-r_input[r_input$SNP %in% p.dat[p.dat$Group==levels(p.dat$Group)[i],]$SNP,]
    Ratio_temp<-tdat[,2] / tdat[,(3+i)]
    
    for(j in 1:exp.number){
      if(i == j){
        Ratio_temp<-Ratio_temp
      }else{
        Ratio_temp<- Ratio_temp - ((tdat[,(3+j)]*rmvmr[j,1]) / tdat[,(3+i)])
      }
    }
    if(is.null(Ratios)){
      Ratios<-Ratio_temp
    }else{
      Ratios<-c(Ratios,Ratio_temp)
    }
  }
  
  #QJ calculations
  
  p.dat$coratios<-Ratios
  
  Qjvec<-NULL
  
  for(i in 1:exp.number){
    
    tdat<-p.dat[p.dat$Group==levels(p.dat$Group)[i],]
    
    Qj<-tdat$Wj * (tdat$coratios - rmvmr[i,1])^2
    
    if(is.null(Qjvec)){
      Qjvec<-Qj
    }else{
      Qjvec<-c(Qjvec,Qj)
    }
  }
  
  p.dat$Qjcor<-Qjvec
  
  #Define matrix for recording total Q statistics.
  Qj_out<-matrix(0L, nrow = exp.number, ncol = 2)
  indqj<-rep(0,length(p.dat[,1]))
  
  for(i in 1:exp.number){
    Qj_out[i,1]<-sum(p.dat[p.dat$Group==levels(p.dat$Group)[i],]$Qjcor)
    Qj_out[i,2]<-pchisq(Qj_out[i,1],nrow(p.dat[p.dat$Group==levels(p.dat$Group)[i],])-exp.number,lower.tail = FALSE)
  }
  
  TotalQs<-data.frame(Qj_out)
  names(TotalQs)<- c("q_statistic","p_value")
  row.names(TotalQs)<-levels(p.dat$Group)
  
  
  for(i in 1:length(p.dat[,1])){
    indqj[i]<- pchisq(p.dat$Qjcor[i],1,lower.tail = FALSE)
  }
  
  p.dat$corQjchi<-indqj
  
  out_data<-p.dat[,c(1,2,8,9,10,7)]
  
  names(out_data)<-c("snp","wj","corrected_beta","qj","qj_p","ref_exposure")
  
  multi_return <- function() {
    Out_list <- list("gq" = TotalQs, "qdat" = out_data)
    class(Out_list)<-"RMVMR_Q"
    
    return(Out_list)
  }
  
  OUT<-multi_return()
  
}