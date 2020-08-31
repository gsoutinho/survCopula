#' Bivariate copula random generation.
#' 
#' @description Random Number Generation for Bivariate Copula Functions. Only returns a
#' single pair of random values from  a bivariate copula with marginal distributions X and Y.
#' 
#' 
#' @param v1 A numeric value belong to the interval [0,1], corresponding to the cumulative 
#' density of the first marginal distribution.
#' @param v2 A numeric value belong to the interval [0,1], corresponding to the cumulative
#' density of the first marginal distribution.
#' @param theta A numeric value for the space parameter.
#' @param type Type of copula. Possible options are \code{"clayton"}, \code{"frank"} \code{"FGM"},
#' \code{"AMH"}, \code{"gumbel-hougaard"} and \code{"joe"}. Defaults to \code{"clayton"}.
#' @param typeX Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"},
#' \code{"Unif"} and \code{"Gamma"}. Defaults to \code{"Exp"}.
#' @param num1_X A numeric value for the first parameter of the first marginal distribution. Defaults to \code{"Exp"}.
#' @param num2_X A numeric value for the second parameter of the first marginal distribution. 
#' Only required for two parameter distributions.
#' @param typeY Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"},
#' \code{"Unif"} and \code{"Gamma"}.Defaults to \code{"Exp"}.
#' @param num1_Y A numeric value for the first parameter of the second marginal distribution.
#' @param num2_Y A numeric value for the second parameter of the second marginal distribution.
#' Only required for two parameter distributions.
#' @return 
#' 2-dimensional vector for the random variables.
#' See also
#' dgCopula, rcopula
#' @examples
#' clay<-copula(0.6, 0.4, theta=2, type='clayton', typeX='Exp', num1_X=0.56, typeY='Exp', num1_Y=0.90) 
#' clay[1]
#' AMH<-copula(0.4, 0.4, theta=0.59, type='AMH', typeX='Norm', num1_X=0.56, num2_X=0.3, typeY='Gamma', 
#' num1_Y=0.90, num2_Y=0.30)
#' AMH[2]
#' 
#' @references
#' Nelsen R.B. (2006). An Introduction to Copulas (2nd ed.), Springer-Verlag.
#' 

copula<-function(v1, v2, theta=theta, type='clayton', typeX='Exp', num1_X=1, num2_X=NULL, 
                 typeY='Exp', num1_Y=1, num2_Y=NULL){
  
  
  if(!is.numeric(v1))
    stop("v1 must be a numeric value")
  
  if(!is.numeric(v2))
    stop("v2 must be a numeric value")
  
  
  #algoritmos distribuicoes condicionais:
  
  if(type=='clayton'){
    
    u1<-v1
    
    u2<-((v1^(-theta))*(v2^(-theta/(1+theta))-1)+1)^(-1/theta)
    
  }
  
  
  if(type=='frank'){
    
    u1<-v1
    
    u2<-(-1/theta)*log(1+(v2*(1-exp(-theta)))/(v2*(exp(-theta*u1)-1)-exp(-theta*u1)))
    
  }
  
  if(type=='FGM'){
    
    u1<-v1
    
    a<-1+theta*(1-2*u1)
    
    b<-sqrt(a^2-4*(a-1)*v2)
    
    u2<-2*v2/(a+b)
    
  }
  
  if(type=='AMH'){
    
    t<-runif(1,0,1)
    
    u1<-v1
    
    a<-(1-v1)
    
    b<-1-theta*(1+2*a*v2)+2*(theta^2)*(a^2)*v2
    
    c<-1+theta*(2-4*a+4*a*v2)+(theta^2)*(1-4*a*v2+4*(a^2)*v2)
    
    u2<-(2*t*(a*theta-1)^2)/(b+sqrt(c))
    
  }
  
  #metodos para distribuicao bivariada:
  
  if(type=='gumbel-hougaard'){
    
    fun<-function (x) x*((1-log(x))/theta)-v2
    
    t <- uniroot(fun, c(0.00001, 1))$root
    
    u1<-exp(v1^(1/theta)*log(t))
    
    u2<-exp((1-v1)^(1/theta)*log(t))
    
    
  }
  
  
  if(type=='joe'){
    
    fun<-function (x) x-((log(1-(1-x)^theta))*(1-(1-x)^theta))/(theta*(1-x)^(theta-1))-v2
    
    t <- uniroot(fun, c(0.00001, 1))$root
    
    u1<-1-(1-(1-(1-t)^theta)^v1)^(1/theta)
    
    u2<-1-(1-(1-(1-t)^theta)^(1-v1))^(1/theta)
    
  }
  

  x<-invF(u1, type=typeX, num1=num1_X, num2=num2_X)
  
  y<-invF(u2, type=typeY, num1=num1_Y, num2=num2_Y)
  
  #res<-list(x=x, y=y, theta=theta, typeX=typeX,  num1=num1_X, num2=num2_X,typeY=typeY, num1=num1_Y, num2=num2_Y)
  
  res<-c(x,y)
  
  return(res)
  
}