#' Random Number Generation for Copula Functions.
#' @description Random Number Generation for Bivariate Copula Functions. Returns a number of pairs of random values 
#' from a bivariate copula with marginal distributions X and Y. 
#' @param typeCopula Type of copula. Possible options are \code{"clayton"}, \code{"frank"} \code{"FGM"}, \code{"AMH"},
#' \code{"gumbel-hougaard"} and \code{"joe"}. Defaults to \code{"clayton"}.
#' @param theta A numeric value for the space parameter.
#' @param typeX Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"} \code{"Unif"} and
#' \code{"Gamma"}. Defaults to \code{"Exp"}.
#' @param num1_X A numeric value for the first parameter of the first marginal distribution.
#' @param num2_X A numeric value for the second parameter of the first marginal distribution.
#' Only required for two parameter distributions.
#' @param typeY Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"} \code{"Unif"} and \code{"Gamma"}.
#' Defaults to \code{"Exp"}.
#' @param num1_Y A numeric value for the first parameter of the second marginal distribution.
#' @param num2_Y A numeric value for the second parameter of the second marginal distribution. 
#' Only required for two parameter distributions. 
#' @param nsim Number of observations to be generated.
#' 
#' @return
#' 2-dimensional random vector with the results of the simulation.
#' @examples 
#' 
#' res<-rcopula(typeCopula = 'clayton', theta = 2, typeX='Exp', num1_X=0.9, 
#'              typeY='Exp', num1_Y=0.3, nsim=1000)
#' 
#' res
#' 
#' res2<-rcopula(typeCopula = 'AMH', theta = 2, typeX='Norm', num1_X=0.9, num2_X=0.3, 
#'               typeY='Gamma', num1_Y=3, num2_Y=2, nsim=1000)
#'               
#' res2[,2]
#' 
#' @author Gustavo Soutinho, Luis Meira-Machado

rcopula<- function(typeCopula = 'clayton', theta = 1, typeX='Exp', num1_X=1, num2_X=NULL, 
                   typeY='Exp', num1_Y=1, num2_Y=NULL, nsim=500){
  
  TAB<-NULL
  
  
  for(i in 1:nsim){
    
    #i<-1
    
    v1<-runif(1,0,1)
    
    v2<-runif(1,0,1)
    
    res<-copula(v1, v2, theta=theta,type=typeCopula,typeX=typeX, num1_X=num1_X, num2_X=num2_X,
                typeY=typeY, num1_Y=num1_Y, num2_Y=num2_Y)
    
    x<-res[1]
    
    y<-res[2]
    
    TAB<-rbind(TAB,cbind(i, x, y))
    
  }
  
  TAB<-as.data.frame(TAB)
  
  colnames(TAB)<-c('ID','X','Y')
  
  #res<-list(tab=TAB, typeCopula=typeCopula,teta=teta,typeX=typeX, num1_X=num1_X, num2_X=num2_X, 
  #          typeY=typeY, num1_Y=num1_Y, num2_Y=num2_Y,
  #          nsim=nsim)
  
  
  res<-TAB[,2:3]
  
  return(res)
  
}