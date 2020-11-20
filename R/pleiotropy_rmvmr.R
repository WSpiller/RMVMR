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
#' \item{\code{qstat}}{A data frame containing the global Q-statistic and p-value after applying a correction for each exposure}
#' \item{\code{qall}}{A data frame containing the individual Q-statistic and p-value for each SNP after applying a correction for each exposure}
#'}
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Sanderson, E., et al., An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology, 2019, 48, 3, 713-727. <https://dx.doi.org/10.1093/ije/dyy262>
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
#' q_object$qstat
#' q_object$qall

pleiotropy_rmvmr<-function(r_input,rmvmr){
  
  library(RadialMR)
  
  rmvmr<-rmvmr$coef
  
  exp.number<-length(names(r_input)[-c(1,2,3)])/2
  
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
  
  Xlist<-NULL
  
  for(i in 1:exp.number){
    
    #Format data for univariate MR using significant SNPs for each exposure
    Xsub<-r_input[f.vec[,i] == 1,]
    Xrad.dat<-format_radial(Xsub[,3+i],Xsub[,2],Xsub[,3 + exp.number + i],Xsub[,3],Xsub[,1])
    
    X.res<-ivw_radial(Xrad.dat,0.05/nrow(Xrad.dat),1,0.0001,F)
    
    
    if(is.null(Xlist)){
      Xlist<-X.res
    }else{
      Xlist<-append(Xlist,X.res)
    }
    
  }
  
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
  
  ###############
  ###############
  ###############
  
  #### Correction Plot
  
  correction.vec<-NULL
  
  for(i in 1:exp.number){
    
    if(is.null(correction.vec)){
      Tempdat<-r_input[r_input$SNP %in% p.dat$SNP[p.dat$Group== levels(p.dat$Group)[i]],]
      correction.vec<-p.dat[p.dat$Group== levels(p.dat$Group)[i],]$BetaWj /
        p.dat[p.dat$Group== levels(p.dat$Group)[i],]$Wj
      
      for(j in 1:exp.number){
        
        if(i==j){
          correction.vec<-correction.vec
        }else{
          correction.vec<- correction.vec - ((Tempdat[,3+j]*rmvmr[j,1])/Tempdat[,3+i])
        }
        
      }
    }else{
      Tempdat<-r_input[r_input$SNP %in% p.dat$SNP[p.dat$Group== levels(p.dat$Group)[i]],]
      correction.vec2<-p.dat[p.dat$Group== levels(p.dat$Group)[i],]$BetaWj /
        p.dat[p.dat$Group== levels(p.dat$Group)[i],]$Wj
      
      for(j in 1:exp.number){
        
        if(i==j){
          correction.vec2<-correction.vec2
        }else{
          correction.vec2<- correction.vec2 - ((Tempdat[,3+j]*rmvmr[j,1])/Tempdat[,3+i])
        }
        
      }
      
      correction.vec<-c(correction.vec,correction.vec2)
      
    }
  }
  
  indiq<-NULL
  Tempdat<-NULL
  
  for(i in 1:exp.number){
    
    Xdat<-data.frame(Xlist[5+((i-1)*13)])
    
    pdat.t<-p.dat
    
    #Define Q statistic for each individual variant
    Tempdat<-r_input[r_input$SNP %in% p.dat$SNP[p.dat$Group== levels(p.dat$Group)[i]],]
    Qj1<- (Xdat[,3]/Xdat[,2])
    
    for(j in 1:exp.number){
      if(i == j){
        Qj1<-Qj1
      }else{
        Qj1<- Qj1 - ((Tempdat[,3+j]*rmvmr[j,1])/Tempdat[,3+i])
      }
    }
    
    if(is.null(indiq)){
      indiq<-Qj1
      
    }else{
      indiq<-c(indiq,Qj1)
    }
  }
  
  indiqchi<-NULL
  
  for(i in 1:length(indiq)){
    
    
    if(is.null(indiqchi)){
      indiqchi<-pchisq(indiq[i],1,lower.tail = FALSE)
    }else{
      indiqchi<-c(indiqchi,pchisq(indiq[i],1,lower.tail = FALSE))
    }
    
  }
  
  pdat.t<-cbind(p.dat,indiq)
  
  TQ<-NULL
  TQchi<-NULL
  
  for(i in 1:exp.number){
    
    Xt<-pdat.t[pdat.t$Group == levels(pdat.t$Group)[i],]
    
    indiq2<-Xt$Wj*(Xt$indiq-rmvmr[i,1])^2
    
    if(is.null(TQ)){
      TQ<-sum(indiq2)
      TQchi<-pchisq(sum(indiq2),(nrow(Xt)-1),lower.tail = FALSE)
    }else{
      TQ<-c(TQ,sum(indiq2))
      TQchi<-c(TQchi,pchisq(sum(indiq2),(nrow(Xt)-1),lower.tail = FALSE))
    }
    
  }
  
  
  IQs<-data.frame(p.dat$SNP,indiq,indiqchi,p.dat$Group)
  names(IQs)<-c("SNP","qstat","q_chi","exposure")
  
  Total_Q<-data.frame(TQ,TQchi)
  names(Total_Q)<-c("global_qstat","gq_chi")
  
  
  for(i in 1:exp.number){
    row.names(Total_Q)[i] <- paste0("Exposure_",i,collapse="")
    
  }
  
  ###############
  ###############
  ###############
  
  multi_return <- function() {
    Out_list <- list("qstat" = Total_Q, "qall" = IQs)
    class(Out_list)<-"RMVMR_Q"
    
    return(Out_list)
  }
  
  OUT<-multi_return()
  
}
