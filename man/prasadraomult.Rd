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
  }



\keyword{ prasad rao }
\keyword{ fay Herriot }
\keyword{ small area estimation }
\keyword{ variance component }
\keyword{multivariate prasad rao}
\keyword{ multivariate fay Herriot }
