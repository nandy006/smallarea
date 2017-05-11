\name{prasadraomult}
\alias{prasadraomult}
\title{
Estimate of the matrix of variance component in Fay Herriot Model using Prasad Rao Method
}
\description{
This function returns a matrix of the variance component in the multivariate Fay Herriot Model. 
It is an unbiased estimator based on method of moments and is analogous to the Univariate Fay herriot estimator. Note that our function does not accept any missing values. Also this function should not be use if you have only one small area attribute. Use the function prasadrao instead.

}
\usage{
prasadraomult(directest.mat, samplingvar.mat, samplingcov.mat, design.mat)
}
\arguments{
  \item{directest.mat}{
  An nxq matrix of direct survey based estimates of q small area attributes on n small areas
}
  \item{samplingvar.mat}{
  An nxq matrix of known sampling variances of the q small area attributes on n small areas. Note that the i-th row of this matrix represents the sampling variances of the q small area attributes associated with the i-th small area.
}
  \item{samplingcov.mat}{
  An nx(q(q-1)/2) matrix of known sampling covariances of the q small area attributes on n small areas. Note that the i-th row of this matrix represents the sampling covariances of the q small area attributes associated with the i-th small area. The covariance terms in the matrix should occur in natural order, for example with 4 attributes the 6 columns of this matrix are cov(y1,y2), cov(y1,y3), cov(y1,y4), cov(y2,y3),  cov(y2,y4),cov(y3,y4) respectively.
}
  \item{design.mat}{
A numeric matrix. The first column is a column of ones(also called the intercept). 
The other columns consist of observations of each of the covariates or the explanatory 
variable in Fay Herriot Model.

}
}
\details{

}
\value{
\item{estimate}{estimate of the matrix of variance component}
}
\references{

}
\author{
Abhishek Nandy
}
\note{

}

\seealso{
\code{\link{prasadraoest}}
}
\examples{

set.seed(55)
n <- 500                                          # Number of small areas
p <- 1                                            # Number of covariates (1 implies intercept only)
q <- 5                                            # Number of small area attributes for each small area
rho <- 0.3
psi <- matrix( c( 1, rho, rho^2, rho^3, rho^4,    # Matrix of variance component
                  rho, 1, rho,rho^2, rho^3,
                  rho^2, rho, 1 , rho, rho^2,
                  rho^3, rho^2, rho, 1, rho^3,
                  rho^4, rho^3,rho^2, rho, 1), q, q )
require( MASS )
Var <- matrix( rep( 1 , times = n*q ), n, q )     # Matrix of Sampling variances
samplingvar.mat <- Var

rho <- c( runif( n/2 , 0.8, 0.9 ), runif( n/2, -0.9, -0.8 ) )
rho2 <- rho^2
rho3 <- rho^3
rho4 <- rho^4

Cov <- as.matrix( cbind( rho, rho2, rho3, rho4, rho, rho2, rho3, rho, rho2, rho ) ) # Matrix of sampling covariances
samplingcov.mat <- Cov

theta <- rep( 0, times = n*q )                    
directest <- rep( 0, times = n*q )




for( i in 1:n ){                              # Simulating unknown small area means
  theta[ ( ( i -1 )*q + 1 ):( i*q ) ]  <- mvrnorm( 1, mu =  rep(0, q) , 
                                                   Sigma = psi )
}




  directest <- rep( 0, times = n*q )
  for( i in 1:n ){                            # Simulating direct survey based estimators
    directest[ ( ( i -1 )*q + 1 ):( i*q ) ]  <- mvrnorm( 1, mu = theta[ ( ( i -1 )*q + 1 ):( i*q ) ], 
                                                         Sigma = matrix( c( 1, rho[i], rho[i]^2, rho[i]^3, rho[i]^4,
                                                                            rho[i], 1, rho[i],rho[i]^2, rho[i]^3,
                                                                            rho[i]^2, rho[i], 1 , rho[i], rho[i]^2,
                                                                            rho[i]^3, rho[i]^2, rho[i], 1, rho[i]^3,
                                                                            rho[i]^4, rho[i]^3,rho[i]^2, rho[i], 1)
                                                                         , q, q ))
  }
  
  
  directest.mat <- matrix( directest, n, q, byrow = TRUE )
  prasadraomult( directest.mat,  samplingvar.mat, samplingcov.mat )

  }



\keyword{ prasad rao }
\keyword{ fay Herriot }
\keyword{ small area estimation }
\keyword{ variance component }
\keyword{multivariate prasad rao}
\keyword{ multivariate fay Herriot }
