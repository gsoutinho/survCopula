% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcopula.R
\name{rcopula}
\alias{rcopula}
\title{Random Number Generation for Copula Functions.}
\usage{
rcopula(
  typeCopula = "clayton",
  theta = 1,
  typeX = "Exp",
  num1_X = 1,
  num2_X = NULL,
  typeY = "Exp",
  num1_Y = 1,
  num2_Y = NULL,
  nsim = 500
)
}
\arguments{
\item{typeCopula}{Type of copula. Possible options are \code{"clayton"}, \code{"frank"} \code{"FGM"}, \code{"AMH"},
\code{"gumbel-hougaard"} and \code{"joe"}. Defaults to \code{"clayton"}.}

\item{theta}{A numeric value for the space parameter.}

\item{typeX}{Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"} \code{"Unif"} and
\code{"Gamma"}. Defaults to \code{"Exp"}.}

\item{num1_X}{A numeric value for the first parameter of the first marginal distribution.}

\item{num2_X}{A numeric value for the second parameter of the first marginal distribution.
Only required for two parameter distributions.}

\item{typeY}{Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"} \code{"Unif"} and \code{"Gamma"}.
Defaults to \code{"Exp"}.}

\item{num1_Y}{A numeric value for the first parameter of the second marginal distribution.}

\item{num2_Y}{A numeric value for the second parameter of the second marginal distribution. 
Only required for two parameter distributions.}

\item{nsim}{Number of observations to be generated.}
}
\value{
2-dimensional random vector with the results of the simulation.
}
\description{
Random Number Generation for Bivariate Copula Functions. Returns a number of pairs of random values 
from a bivariate copula with marginal distributions X and Y.
}
\examples{

res<-rcopula(typeCopula = 'clayton', theta = 2, typeX='Exp', num1_X=0.9, 
             typeY='Exp', num1_Y=0.3, nsim=1000)

res

res2<-rcopula(typeCopula = 'AMH', theta = 2, typeX='Norm', num1_X=0.9, num2_X=0.3, 
              typeY='Gamma', num1_Y=3, num2_Y=2, nsim=1000)
              
res2[,2]

}
\author{
Gustavo Soutinho, Luis Meira-Machado
}
