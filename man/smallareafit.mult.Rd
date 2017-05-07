\name{smallareafit.mult}
\alias{smallareafit.mult}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
smallareafit.mult(directest.mat, samplingvar.mat, samplingcov.mat, design.mat)
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
    require("Rmpfr")
    options(digits = 16)
    if (missing(design.mat)) {
        n <- nrow(directest.mat)
        design.mat <- matrix(1, n, 1)
        p <- 1
        q <- ncol(directest.mat)
        prasadraomult.nocov <- function(directest.mat, samplingvar.mat, 
            samplingcov.mat) {
            areas <- nrow(directest.mat)
            residual.transpose <- t(directest.mat) - rowMeans(t(directest.mat))
            leverage <- rep(1/areas, areas)
            adjust.varmat <- replicate(ncol(samplingvar.mat), 
                1 - leverage) * samplingvar.mat
            adjust.var <- colSums(adjust.varmat)
            adjust.covmat <- replicate(ncol(samplingcov.mat), 
                1 - leverage) * samplingcov.mat
            adjust.cov <- colSums(adjust.covmat)
            if (length(adjust.cov) == 1) {
                adjust.cov = as.numeric(adjust.cov)
            }
            adjustment.mat <- diag(adjust.var)
            adjustment.mat[lower.tri(adjustment.mat, diag = FALSE)] <- adjust.cov
            adjustment.mat <- forceSymmetric(adjustment.mat, 
                uplo = "L")
            estimate <- (residual.transpose \%*\% t(residual.transpose) - 
                adjustment.mat)/(areas - 1)
            if (det(estimate) <= 0) {
                return(1e-04 * diag(q))
            }
            else {
                return(estimate)
            }
        }
        require(Matrix)
        x <- 1:n
        psi.hat <- prasadraomult.nocov(directest.mat, samplingvar.mat, 
            samplingcov.mat)
        options(digits = 16)
        Sigmainv.list <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            D <- D + psi.hat
            D <- forceSymmetric(D, uplo = "L")
            return(solve(D))
        })
        inv.Sigma <- bdiag(Sigmainv.list)
        designmatexpanded.list <- lapply(x, function(x) {
            as.matrix(bdiag(replicate(q, t(as.matrix(design.mat[x, 
                ])), simplify = F)))
        })
        designmat.expanded <- do.call(rbind, designmatexpanded.list)
        designmat.expanded <- Matrix(designmat.expanded, sparse = TRUE)
        directest.vec <- as.vector(t(directest.mat))
        beta.estimate <- solve(t(designmat.expanded) \%*\% inv.Sigma \%*\% 
            designmat.expanded, sparse = TRUE) \%*\% (t(designmat.expanded) \%*\% 
            directest.vec)
        options(digits = 16)
        B.est <- matrix(beta.estimate, q, p, byrow = TRUE)
        if (det(psi.hat) <= 0.01) {
            thetaest.list <- lapply(x, function(x) {
                D <- diag(samplingvar.mat[x, ])
                D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                  ]
                D <- forceSymmetric(D, uplo = "L")
                B.est \%*\% design.mat[x, ]
            })
            theta.est <- t(do.call(cbind, thetaest.list))
        }
        else {
            thetaest.list <- lapply(x, function(x) {
                D <- diag(samplingvar.mat[x, ])
                D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                  ]
                D <- forceSymmetric(D, uplo = "L")
                psi.hat <- prasadraomult.nocov(directest.mat, 
                  samplingvar.mat, samplingcov.mat)
                solve((solve(D) + solve(psi.hat))) \%*\% (solve(D) \%*\% 
                  directest.mat[x, ] + solve(psi.hat) \%*\% B.est \%*\% 
                  design.mat[x, ])
            })
            theta.est <- t(do.call(cbind, thetaest.list))
        }
        g1 <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            g1 <- sum(diag(psi.hat \%*\% solve(D + psi.hat) \%*\% 
                D))
            return(g1)
        })
        g1 <- unlist(g1)
        mymat <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            return(solve(psi.hat + D))
        })
        mymat <- solve(Reduce("+", mymat))
        g2 <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            mymat1 <- psi.hat \%*\% solve(D + psi.hat) \%*\% D
            g2 <- sum(diag(solve(psi.hat) \%*\% mymat1 \%*\% mymat1 \%*\% 
                solve(psi.hat) \%*\% mymat))
            return(g2)
        })
        g2 <- unlist(g2)
        z <- 1:n
        fun5 <- function(z) {
            u <- matrix(NA, q, q)
            for (index1 in 1:q) {
                for (index2 in 1:q) {
                  fun4 <- function(x) {
                    D <- diag(samplingvar.mat[x, ])
                    D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                      ]
                    D <- forceSymmetric(D, uplo = "L")
                    mat <- psi.hat + D
                    fun1 <- function(index1, index2, index3, 
                      index4) {
                      mat[index1, index2] * mat[index3, index4] + 
                        mat[index3, index2] * mat[index1, index4]
                    }
                    fun3 <- function(index3) {
                      index4 <- 1:q
                      fun2 <- function(index4) {
                        fun1(index1, index2, index3, index4)
                      }
                      u <- lapply(index4, fun2)
                      u <- unlist(u)
                      return(u)
                    }
                    index3 <- 1:q
                    v <- lapply(index3, fun3)
                    return(matrix(unlist(v), nrow = q, byrow = TRUE))
                  }
                  v <- Reduce("+", lapply(x, fun4))
                  D <- diag(samplingvar.mat[z, ])
                  D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[z, 
                    ]
                  D <- forceSymmetric(D, uplo = "L")
                  psi.hat <- prasadraomult.nocov(directest.mat, 
                    samplingvar.mat, samplingcov.mat)
                  dummy <- solve(D + psi.hat) \%*\% v
                  u[index1, index2] <- sum(diag(dummy))/(nrow(design.mat) - 
                    ncol(design.mat))^2
                }
            }
            D <- diag(samplingvar.mat[z, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[z, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            g3 <- solve(psi.hat + D) \%*\% D \%*\% solve(psi.hat + 
                D) \%*\% u
            g3 <- sum(diag(g3))
            return(g3)
        }
        g3 <- lapply(z, fun5)
        g3 <- unlist(g3)
        mspe <- g1 + 2 * g2 + g3
        return(list(theta.est = theta.est, B.estimate = B.est, 
            mspe = mspe, psi.hat = psi.hat))
    }
    else {
        n <- nrow(directest.mat)
        q <- ncol(directest.mat)
        p <- ncol(design.mat)
        prasadraomult <- function(directest.mat, samplingvar.mat, 
            samplingcov.mat, design.mat) {
            require(Matrix)
            require("Rmpfr")
            projection <- function(mat) {
                mat \%*\% solve(t(mat) \%*\% mat) \%*\% t(mat)
            }
            areas <- nrow(directest.mat)
            proj <- projection(design.mat)
            fit <- lm(directest.mat ~ design.mat)
            residual.transpose <- t(fit$residuals)
            leverage <- as.vector(diag(proj))
            adjust.varmat <- replicate(ncol(samplingvar.mat), 
                1 - leverage) * samplingvar.mat
            adjust.var <- colSums(adjust.varmat)
            adjust.covmat <- replicate(ncol(samplingcov.mat), 
                1 - leverage) * samplingcov.mat
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
        require(Matrix)
        x <- 1:n
        psi.hat <- prasadraomult(directest.mat, samplingvar.mat, 
            samplingcov.mat, design.mat)
        options(digits = 16)
        Sigmainv.list <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            D <- D + psi.hat
            D <- forceSymmetric(D, uplo = "L")
            return(solve(D))
        })
        inv.Sigma <- bdiag(Sigmainv.list)
        designmatexpanded.list <- lapply(x, function(x) {
            as.matrix(bdiag(replicate(q, t(as.matrix(design.mat[x, 
                ])), simplify = F)))
        })
        designmat.expanded <- do.call(rbind, designmatexpanded.list)
        designmat.expanded <- Matrix(designmat.expanded, sparse = TRUE)
        directest.vec <- as.vector(t(directest.mat))
        beta.estimate <- solve(t(designmat.expanded) \%*\% inv.Sigma \%*\% 
            designmat.expanded, sparse = TRUE) \%*\% (t(designmat.expanded) \%*\% 
            directest.vec)
        options(digits = 16)
        B.est <- matrix(beta.estimate, q, p, byrow = TRUE)
        if (det(psi.hat) <= 0.01) {
            thetaest.list <- lapply(x, function(x) {
                D <- diag(samplingvar.mat[x, ])
                D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                  ]
                D <- forceSymmetric(D, uplo = "L")
                psi.hat <- prasadraomult(directest.mat, samplingvar.mat, 
                  samplingcov.mat, design.mat)
                B.est \%*\% design.mat[x, ]
            })
            theta.est <- t(do.call(cbind, thetaest.list))
        }
        else {
            thetaest.list <- lapply(x, function(x) {
                D <- diag(samplingvar.mat[x, ])
                D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                  ]
                D <- forceSymmetric(D, uplo = "L")
                psi.hat <- prasadraomult(directest.mat, samplingvar.mat, 
                  samplingcov.mat, design.mat)
                solve((solve(D) + solve(psi.hat))) \%*\% (solve(D) \%*\% 
                  directest.mat[x, ] + solve(psi.hat) \%*\% B.est \%*\% 
                  design.mat[x, ])
            })
            theta.est <- t(do.call(cbind, thetaest.list))
        }
        g1 <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            g1 <- sum(diag(psi.hat \%*\% solve(D + psi.hat) \%*\% 
                D))
            return(g1)
        })
        g1 <- unlist(g1)
        mymat <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            design <- as.vector(design.mat[x, ])
            design <- design \%o\% design
            return(kronecker(solve(psi.hat + D), design))
        })
        mymat <- solve(Reduce("+", mymat))
        g2 <- lapply(x, function(x) {
            D <- diag(samplingvar.mat[x, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            mymat1 <- psi.hat \%*\% solve(D + psi.hat) \%*\% D
            design <- as.vector(design.mat[x, ])
            mymat2 <- kronecker(solve(psi.hat), t(design))
            g2 <- sum(diag(t(mymat2) \%*\% mymat1 \%*\% mymat1 \%*\% 
                mymat2 \%*\% mymat))
            return(g2)
        })
        g2 <- unlist(g2)
        z <- 1:n
        fun5 <- function(z) {
            u <- matrix(NA, q, q)
            for (index1 in 1:q) {
                for (index2 in 1:q) {
                  fun4 <- function(x) {
                    D <- diag(samplingvar.mat[x, ])
                    D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[x, 
                      ]
                    D <- forceSymmetric(D, uplo = "L")
                    mat <- psi.hat + D
                    fun1 <- function(index1, index2, index3, 
                      index4) {
                      mat[index1, index2] * mat[index3, index4] + 
                        mat[index3, index2] * mat[index1, index4]
                    }
                    fun3 <- function(index3) {
                      index4 <- 1:q
                      fun2 <- function(index4) {
                        fun1(index1, index2, index3, index4)
                      }
                      u <- lapply(index4, fun2)
                      u <- unlist(u)
                      return(u)
                    }
                    index3 <- 1:q
                    v <- lapply(index3, fun3)
                    return(matrix(unlist(v), nrow = q, byrow = TRUE))
                  }
                  v <- Reduce("+", lapply(x, fun4))
                  D <- diag(samplingvar.mat[z, ])
                  D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[z, 
                    ]
                  D <- forceSymmetric(D, uplo = "L")
                  dummy <- solve(D + psi.hat) \%*\% v
                  u[index1, index2] <- sum(diag(dummy))/(nrow(design.mat) - 
                    ncol(design.mat))^2
                }
            }
            D <- diag(samplingvar.mat[z, ])
            D[lower.tri(D, diag = FALSE)] <- samplingcov.mat[z, 
                ]
            D <- forceSymmetric(D, uplo = "L")
            psi.hat <- prasadraomult(directest.mat, samplingvar.mat, 
                samplingcov.mat, design.mat)
            g3 <- solve(psi.hat + D) \%*\% D \%*\% solve(psi.hat + 
                D) \%*\% u
            g3 <- sum(diag(g3))
            return(g3)
        }
        g3 <- lapply(z, fun5)
        g3 <- unlist(g3)
        mspe <- g1 + 2 * g2 + g3
        return(list(theta.est = theta.est, B.estimate = B.est, 
            mspe = mspe, psi.hat = psi.hat))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
