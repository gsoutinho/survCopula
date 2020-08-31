#' Inverse of cumulative function of a marginal distribution.
#' 
#' @description The value of the inverse cumulative distribution.
#' @param u A numeric value belong to the interval [0,1], corresponding to the cumulative density.
#' @param type Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"}, 
#' \code{"Unif"} and \code{"Gamma"}. Defaults to \code{"Exp"}.
#' @param num1 A numeric value for the first parameter of the marginal distribution.
#' @param num2 A numeric value for the second parameter of the marginal distribution.
#' @return
#' A numeric value corresponding to the distribution of the marginal function.
#' @examples 
#' invF(0.2)
#' invF(0.2,type = 'Norm', num1 = 0.2,num2 = 0.1)

#' @author Gustavo Soutinho, Luis Meira-Machado

invF<-function(u, type='Exp',num1=1, num2=NULL){
  
  if(type=='Exp'){
    
    res<-(-1/num1)*log(1-u)  #qexp(0.5, par1) #u=0.5
    
  }
  
  if(type=='Norm'){
    
    
    res<-qnorm(u, mean=num1, sd=num2)
    
  }
  
  if(type=='Unif'){
    
    
    res<-qunif(u, min=num1, max=num2)
    
  }
  
  if(type=='Gamma'){
    
    res<-qgamma(u, shape= num1, scale = num2)
    
  }
  
  return(res)
}