% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copula.R
\name{copula}
\alias{copula}
\title{Bivariate copula random generation.}
\usage{
copula(
  v1,
  v2,
  theta = theta,
  type = "clayton",
  typeX = "Exp",
  num1_X = 1,
  num2_X = NULL,
  typeY = "Exp",
  num1_Y = 1,
  num2_Y = NULL
)
}
\arguments{
\item{v1}{A numeric value belong to the interval [0,1], corresponding to the cumulative 
density of the first marginal distribution.}

\item{v2}{A numeric value belong to the interval [0,1], corresponding to the cumulative
density of the first marginal distribution.}

\item{theta}{A numeric value for the space parameter.}

\item{type}{Type of copula. Possible options are \code{"clayton"}, \code{"frank"} \code{"FGM"},
\code{"AMH"}, \code{"gumbel-hougaard"} and \code{"joe"}. Defaults to \code{"clayton"}.}

\item{typeX}{Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"},
\code{"Unif"} and \code{"Gamma"}. Defaults to \code{"Exp"}.}

\item{num1_X}{A numeric value for the first parameter of the first marginal distribution. Defaults to \code{"Exp"}.}

\item{num2_X}{A numeric value for the second parameter of the first marginal distribution. 
Only required for two parameter distributions.}

\item{typeY}{Type of marginal distribution. Possible options are \code{"Exp"}, \code{"Norm"},
\code{"Unif"} and \code{"Gamma"}.Defaults to \code{"Exp"}.}

\item{num1_Y}{A numeric value for the first parameter of the second marginal distribution.}

\item{num2_Y}{A numeric value for the second parameter of the second marginal distribution.
Only required for two parameter distributions.}
}
\value{
2-dimensional vector for the random variables.
See also
dgCopula, rcopula
}
\description{
Random Number Generation for Bivariate Copula Functions. Only returns a
single pair of random values from  a bivariate copula with marginal distributions X and Y.
}
\examples{
clay<-copula(0.6, 0.4, theta=2, type='clayton', typeX='Exp', num1_X=0.56, typeY='Exp', num1_Y=0.90) 
clay[1]
AMH<-copula(0.4, 0.4, theta=0.59, type='AMH', typeX='Norm', num1_X=0.56, num2_X=0.3, typeY='Gamma', 
num1_Y=0.90, num2_Y=0.30)
AMH[2]

}
\references{
Nelsen R.B. (2006). An Introduction to Copulas (2nd ed.), Springer-Verlag.
}
