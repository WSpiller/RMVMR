#' plot_rmvmr
#'
#' Generates two radial multivariable Mendelian randomization (MVMR) plots. The first plot shows shows the estimated direct effect for each exposure obtained by fitting a radial MVMR model.
#' Each data point shows the square root weighting for each SNP on the x-axis, and product of the ratio estimate and square root weighting for each SNP on the y-axis. These values are obtained
#' exposure are used, such that their first stage F-statistic is greater than 10. The second plot applies a correction to each ratio estimate. In both plots, the distance of each observation from the 
#' by performing a univariate radial MR analysis for each exposure using the SNPs displayed, specifically through use of the \code{RadialMR::ivw_radial} function. Only SNPs strongly associated with the corresponding
#' corresponding regression line is proportional to the contribution of that SNP towards global heterogeneity.
#'
#' @param r_input A formatted data frame using the \code{format_rmvmr} function.
#' @param rmvmr An object containing the output from the \code{ivw_rmvmr} function of class \code{IVW_RMVMR}.
#'
#' @return An object of class \code{"RMVMR_plot"} containing the following components:\describe{
#' \item{\code{p1}}{A radial MVMR plot without correction}
#' \item{\code{p2}}{A radial MVMR plot with correction}
#'}
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Spiller, W., et al., Estimating and visualising multivariable Mendelian randomization analyses within a radial framework. Forthcoming.
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
#'
#' 
#' rmvmr_output <- ivw_rmvmr(f.data, FALSE)
#' plot_object <- plot_rmvmr(f.data, rmvmr_output)
#' plot_object$p1
#' plot_object$p2

  
  
plot_rmvmr <- function(r_input, rmvmr){

  # to suppress the R CMD check note about: no visible binding for global variable
  BetaWj <- Group <- Wj <- wj <- ref_exposure <- corrected_beta <- NULL

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
    Xrad.dat <- RadialMR::format_radial(Xsub[,3+i],Xsub[,2],Xsub[,3 + exp.number + i],Xsub[,3],Xsub[,1])

    X.res <- RadialMR::ivw_radial(Xrad.dat,0.05/nrow(Xrad.dat),1,0.0001,FALSE)


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
  
  cpalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#D55E00", "#0072B2", "#CC79A7")

  B <- ggplot2::ggplot(p.dat, ggplot2::aes(x = Wj, y = BetaWj)) +
    ggplot2::labs(title="Radial MVMR without correction") +
    ggplot2::geom_point(ggplot2::aes(colour = Group)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::ylab(expression(hat(beta)[j]~sqrt(W[j]))) +
    ggplot2::xlab(expression(sqrt(W[j]))) +
    ggplot2::scale_x_continuous(limits = c(0,max(p.dat$Wj + 5), expand = c(0,0))) +
    ggplot2::scale_y_continuous(limits = c(min(p.dat$BetaWj - 5), max(p.dat$BetaWj + 5)))

  for(i in 1:exp.number){

    B <- B + ggplot2::geom_segment(x = 0,
                                   xend = max(p.dat$Wj + 5),
                                   y = 0,
                                   yend = rmvmr$coef[i,1]*max(p.dat$Wj + 5),
                                   color = cpalette[i])
  }

  B <- B +
    ggplot2::scale_color_manual(name = "Estimates",
                                breaks = levels(p.dat[,7]),
                                values = cpalette[1:(exp.number)])

  #### Correction Plot

  cordat <- pleiotropy_rmvmr(r_input, rmvmr)

  cpalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#D55E00", "#0072B2", "#CC79A7")
  
  p.dat<-cordat$qdat
  
  C<-ggplot(p.dat,aes(x=wj,y=wj*corrected_beta))+labs(title="Radial MVMR with correction")+ geom_point(aes(colour=ref_exposure))+
    theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
    scale_x_continuous(limits = c(0,max(p.dat$wj+5),expand=c(0,0)))+scale_y_continuous(limits = c(min((p.dat$wj*p.dat$corrected_beta)-5),max((p.dat$wj*p.dat$corrected_beta)+5)))
  
  for(i in 1:exp.number){
    C<- C + geom_segment(x = 0, xend = max(p.dat$wj+5), y = 0, yend = rmvmr$coef[i,1]*max(p.dat$wj+5),color=cpalette[i])
    
  }
  
  C<- C + 
    scale_color_manual(name="Estimates",breaks=levels(p.dat[,6]),values=cpalette[1:(exp.number)])
  
  multi_return <- function() {
    Out_list <- list("p1" = B,"p2"= C)
    class(Out_list)<-"RMVMR_plot"
    
    return(Out_list)
  }
  
  OUT<-multi_return()
  
}

