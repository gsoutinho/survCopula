#' Random generator for multivariate survival data via copulas.
#' 
#' @description This function can be used to generate multivariate 
#' survival data in a variety of scenarios including competing risks, recurrent event and multi-state models.
#' 
#' @param typeCopula Type of copula. Possible options are \code{"clayton"}, \code{"frank"} \code{"FGM"}, 
#' \code{"AMH"}, \code{"gumbel-hougaard"} and \code{"joe"}. Defaults to \code{"clayton"}.
#' @param theta A numeric value for the space parameter.
#' @param typeX Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"} \code{"Unif"} 
#' and \code{"Gamma"}. Defaults to \code{"Exp"}.
#' @param num1_X A numeric value for the first parameter of the first marginal distribution.
#' @param num2_X A numeric value for the second parameter of the first marginal distribution. 
#' Only required for two parameter distributions.
#' @param typeY Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"},
#'  \code{"Unif"} and \code{"Gamma"}. Defaults to \code{"Exp"}.
#' @param num1_Y A numeric value for the first parameter of the second marginal distribution.
#' @param num2_Y A numeric value for the second parameter of the second marginal distribution.
#' Only required for two parameter distributions.
#' @param typeCens Type of censuring distribution. Possible options are \code{"None"}, \code{"Unif"},
#'  \code{"Exp"} and \code{"Wei"}. Defaults to \code{"None"}.
#' @param num1_Cens A numeric value for the first parameter of the censoring distribution.
#' @param num2_Cens A numeric value for the second parameter of the censoring distribution.
#' Only required for two parameter distributions.
#' @param typeSurvData Type of survival data. Possible options are \code{"time-to-event"}, 
#' \code{"recurrent"}, #' \code{"competing-risks"} and \code{"illness-death"}. Defaults to 
#' \code{"illness-death"}.
#' @param state2.prob Probability of a individual move to the intermediate state in illness-death model. 
#' Only required if typeSurvData =”illness-death”. Default to 0.7.
#' @param nsim Number of observations to be generated.
#' @return 
#' A numeric vector with the random multivariate survival data. 
#' Meaning of the colums for the type of survival data:
#' Time-to-event: 
#' \item{T}{Surival time. T = min(Y,Z).}
#' \item{Z}{Censoring time variable.}
#' \item{Delta}{Indicator status for censoring. Delta takes 1 when Y <= Z.}
#' 
#' Recurrent:
#' \item{T1}{First gap time.}
#' \item{Delta1}{Censoring indicator variable for the first gap time.}
#' \item{T2}{Second gap time.}
#' \item{Delta2}{Censoring indicator variable for the second gap time.}
#' \item{Z}{Censoring time variable.}
#' 
#' Competing risks:
#' \item{T}{Survival time.}
#' \item{Z}{Censoring time variable.}
#' \item{Delta}{Indicator status for censoring. Delta takes 0 if the competing risk process 
#' does not move from the initial state at the survival time T, or the value 1 and 2 for the 
#' possible causes of death 1 and 2.} 
#' 
#' Illness-death:
#' \item{T1}{Sojourn time in the initial state.}
#' \item{Delta1}{Indicator status. Delta1 takes 1 when T1<T and T1<Z.}
#' \item{T}{Total survival time.}
#' \item{Delta}{Indicator status for censoring. Delta takes 1 when T<Z.} 
#' \item{Z}{Censoring time variable.}
#' 
#' @examples
#' 
#' sim.data<-dgCopula(typeCopula ='clayton', theta=1, typeX='Unif', num1_X=0, num2_X=5, 
#'                   typeY='Unif',  num1_Y=0, num2_Y=5,  typeCens='Unif', num1_Cens=0, 
#'                   num2_Cens=7,  nsim=250,typeSurvData='time-to-event')
#'
#' head(sim.data)
#'             
#'sim.data2<-dgCopula(typeCopula ='frank', theta=10, typeX='Exp', num1_X=0.5, 
#'                    typeY='Exp', num1_Y=1.5, nsim=250,typeSurvData='illness-death', 
#'                    typeCens='Unif', num1_Cens=0, 
#'                    num2_Cens=4, state2.prob=0.6)
#'             
#' @references 
#' Meira-Machado, L.; de Una-Alvarez, J.; Cadarso-Suarez, C. and Andersen, P.K. Multi-state models for 
#' the analysis of time to event data. Statistical Methods in Medical Research, 2009, 18, 195-222.
#' 
#' Meira-Machado, L., Sestelo, M.; Gonlcalves, A. Nonparametric estimation of the survival function 
#' for ordered multivariate failure time data: A comparative study. Biometrical Journal, 2016, 58, 623-634.
#' 
#' L Meira-Machado, S Faria, A simulation study comparing modeling approaches in an illness-death multi-state model,
#' Communications in Statistics-Simulation and Computation, 2014, 43 (5), 929-946
#' 
#' A Moreira, J de Una-Alvarez, L Machado Presmoothing the Aalen-Johansen estimator in the illness-death model, 
#' Electronic Journal of Statistics, 2013, 7, 1491-1516
#' 
#' A Moreira, L Meira-Machado, survivalBIV: Estimation of the bivariate distribution function for sequentially
#' ordered events under univariate censoring, J Stat Softw, 2012, 46 (13), 1-16
#' 
#' @author Gustavo Soutinho, Luis Meira-Machado
#'  
#' @importFrom "stats" "qgamma" "qnorm" "qunif" "rbinom" "rexp" "runif" "rweibull" "uniroot"

#' @export invF
#' @export copula
#' @export rcopula
#' @export dgCopula


dgCopula<- function(typeCopula = 'clayton', theta = 1, typeX='Exp', num1_X=1, num2_X=NULL, 
                    typeY='Exp', num1_Y=1, num2_Y=NULL, typeCens='None', num1_Cens=NULL, num2_Cens=NULL,
                    typeSurvData='illness-death', 
                    state2.prob=0.7, nsim=250){
  
  TAB<-NULL
  
  for(i in 1:nsim){
    
    #i<-1
    
    v1<-runif(1,0,1)
    
    v2<-runif(1,0,1)
    
    res<-copula(v1, v2, theta=theta,type=typeCopula,typeX=typeX, num1_X=num1_X, num2_X=num2_X,
                typeY=typeY, num1_Y=num1_Y, num2_Y=num2_Y)
    
    x<-res[1]
    
    y<-res[2]
    
    if(typeCens=='None'){
      
      z<-10000000000000000000000000000
      
    }else{
      
      if(typeCens=='Unif'){
        
        z<-runif(1,num1_Cens,num2_Cens)
        
      }
      
      if(typeCens=='Exp'){
        
        z<-rexp(1,num1_Cens)
        
      }
      
      if(typeCens=='Wei'){
        
        z<-rweibull(1,num1_Cens,num2_Cens )
        
      }
      
    }
    
    if(typeSurvData=='time-to-event'){
      
      t<-min(y, z)
      
      delta<-ifelse(y<z, 1,0)
      
      TAB<-rbind(TAB,cbind(t,z, delta))
      
    }
    
    if(typeSurvData=='recurrent'){
      
      t1<-min(x, z)
      
      delta1<-ifelse(x<z, 1,0)
      
      if(x>z){
        
        t2<-0
        
        delta2<-0
        
      }else{
        
        t2<-min(y, z-x)
        
        #delta2<-ifelse(t1+t2<=z, 1,0) estava isto
        
        delta2<-ifelse(x+y<=z,1,0)
      }
      
      
      TAB<-rbind(TAB,cbind(t1, delta1, t2, delta2, z))
      
    }
    
    
    if(typeSurvData=='competing-risks'){
      
      D<-ifelse(x<=y, 1, 2)
      
      t<-min(x,y,z)
      
      delta<-(ifelse(min(x,y)<=z, 1, 0)) * D
      
      TAB<-rbind(TAB,cbind(t,z, delta))
      
    }
    
    if(typeSurvData=='illness-death'){
      
      prob <- rbinom(1, 1, state2.prob)
      
      if(prob==1){
        
        t1<-min(x, z)
        
        delta1<-ifelse(x<z, 1, 0)
        
        t<-min(x+y,z)
        
        delta<-ifelse(x+y<z, 1, 0)
        
        TAB<-rbind(TAB,cbind(t1, delta1, t, delta, z))
        
        
      }else{
        
        w<-runif(1,0,2) #exemplo de distribuicao
        
        t1<-min(w, z)
        
        delta1<-ifelse(w<z, 1, 0)
        
        t<-t1
        
        delta<-delta1
        
        TAB<-rbind(TAB,cbind(t1, delta1, t, delta, z))
      }
      
    }
    
  }#fim simulacao
  
  if(typeSurvData=='time-to-event'){
    
    TAB<-as.data.frame(TAB)
    
    colnames(TAB)<-c('T','Z', 'Delta')
    
    #head(TAB)
    
  }
  
  if(typeSurvData=='recurrent'){
    
    TAB<-as.data.frame(TAB)
    
    colnames(TAB)<-c('T1', 'Delta1', 'T2', 'Delta2', 'Z')
    
    #head(TAB)
    
  }
  
  if(typeSurvData=='competing-risks'){
    
    TAB<-as.data.frame(TAB)
    
    colnames(TAB)<-c('T','Z', 'Delta')
    
    #head(TAB)
  }
  
  if(typeSurvData=='illness-death'){
    
    TAB<-as.data.frame(TAB)
    
    colnames(TAB)<-c('T1','Delta1', 'T','Delta', 'Z')
    
    #head(TAB)
    
    #TAB[1:10,]
    
  }
  
  #res<-list(tab=TAB, typeCopula=typeCopula,teta=teta,nsim=nsim)
  
  res<-TAB
  
  return(res)
  
}