\name{prasadraomult}
\alias{prasadraomult}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prasadraomult(directest.mat, samplingvar.mat, samplingcov.mat, design.mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{directest.mat}{
%%     ~~Describe \code{directest.mat} here~~
}
  \item{samplingvar.mat}{
%%     ~~Describe \code{samplingvar.mat} here~~
}
  \item{samplingcov.mat}{
%%     ~~Describe \code{samplingcov.mat} here~~
}
  \item{design.mat}{
%%     ~~Describe \code{design.mat} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (directest.mat, samplingvar.mat, samplingcov.mat, design.mat) 
{
    require(Matrix)
    require("Rmpfr")
    options(digits = 16)
    prasadraomult.nocov <- function(directest.mat, samplingvar.mat, 
        samplingcov.mat) {
        require(Matrix)
        require("Rmpfr")
        options(digits = 16)
        n <- nrow(directest.mat)
        design.mat <- matrix(1, n, 1)
        p <- 1
        q <- ncol(directest.mat)
        areas <- nrow(directest.mat)
        residual.transpose <- t(directest.mat) - rowMeans(t(directest.mat))
        leverage <- rep(1 - 1/areas, areas)
        adjust.varmat <- replicate(ncol(samplingvar.mat), leverage) * 
            samplingvar.mat
        adjust.var <- colSums(adjust.varmat)
        adjust.covmat <- replicate(ncol(samplingcov.mat), leverage) * 
            samplingcov.mat
        adjust.cov <- colSums(adjust.covmat)
        if (length(adjust.cov) == 1) {
            adjust.cov = as.numeric(adjust.cov)
        }
        adjustment.mat <- diag(adjust.var)
        adjustment.mat[lower.tri(adjustment.mat, diag = FALSE)] <- adjust.cov
        adjustment.mat <- forceSymmetric(adjustment.mat, uplo = "L")
        estimate <- (residual.transpose \%*\% t(residual.transpose) - 
            adjustment.mat)/(areas - 1)
        if (det(estimate) <= 0) {
            return(1e-04 * diag(q))
        }
        else {
            return(estimate)
        }
    }
    if (missing(design.mat)) {
        return(prasadraomult.nocov(directest.mat, samplingvar.mat, 
            samplingcov.mat))
    }
    else {
        projection <- function(mat) {
            mat \%*\% solve(t(mat) \%*\% mat) \%*\% t(mat)
        }
        areas <- nrow(directest.mat)
        proj <- projection(design.mat)
        fit <- lm(directest.mat ~ design.mat)
        residual.transpose <- t(fit$residuals)
        leverage <- as.vector(diag(proj))
        adjust.varmat <- replicate(ncol(samplingvar.mat), 1 - 
            leverage) * samplingvar.mat
        adjust.var <- colSums(adjust.varmat)
        adjust.covmat <- replicate(ncol(samplingcov.mat), 1 - 
            leverage) * samplingcov.mat
        adjust.cov <- colSums(adjust.covmat)
        if (length(adjust.cov) == 1) {
            adjust.cov = as.numeric(adjust.cov)
        }
        adjustment.mat <- diag(adjust.var)
        adjustment.mat[upper.tri(adjustment.mat, diag = FALSE)] <- adjust.cov
        adjustment.mat[lower.tri(adjustment.mat, diag = FALSE)] <- adjust.cov
        estimate <- (residual.transpose \%*\% t(residual.transpose) - 
            adjustment.mat)/(areas - ncol(design.mat) - 1)
        return(estimate)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
