########################################################################################################
################################### Prasad Rao Function ################################################
########################################################################################################

prasadraomult <- function( directest.mat, samplingvar.mat, samplingcov.mat, design.mat ){
  # the data input format is for each small area you have a row with columns giving
  # the values of the direct estimators, followed by colmns giving values of the 
  # covariates, followed by variance of each of the small area attributes followed by 
  # covariance terms in natural order i.e. for example with 4 attributes 
  # cov(y1,y2), cov(y1,y3), cov(y1,y4), cov(y2,y3), 
  # cov(y2,y4),cov(y3,y4) 
  
  require(Matrix)
  require("Rmpfr")
  #options( digits = 16)
  
  prasadraomult.nocov <- function( directest.mat, samplingvar.mat, samplingcov.mat ){
    # the data input format is for each small area you have a row with columns giving
    # the values of the direct estimators, followed by colmns giving values of the 
    # covariates, followed by variance of each of the small area attributes followed by 
    # covariance terms in natural order i.e. for example with 4 attributes 
    # cov(y1,y2), cov(y1,y3), cov(y1,y4), cov(y2,y3), 
    # cov(y2,y4),cov(y3,y4) 
    
    require(Matrix)
    require("Rmpfr")
    options( digits = 16)
    n <- nrow( directest.mat )
    design.mat <- matrix(1, n, 1 )
    p <- 1 
    q <- ncol( directest.mat )
    
    # the projection matrix function
    
    
    # total number of areas
    areas <- nrow( directest.mat )
    
    
    
    #residual.transpose <- t( directest.mat )%*%( diag( areas ) - proj )
    
    residual.transpose <- t( directest.mat ) - rowMeans( t( directest.mat ) )
    
    #residual.transpose <- sqrt(1 + c*log(n)/n)*residual.transpose
    
    # p_ii's  used to compute the prasad rao estimator
    leverage <- rep( 1 - 1/areas, areas )
    
    # this is the adjustment for the variance covariance matrix in the 
    # formula P.R.M-1 in page 67
    
    adjust.varmat <- replicate( ncol( samplingvar.mat ) , leverage )*samplingvar.mat
    adjust.var <- colSums( adjust.varmat )
    
    
    adjust.covmat <- replicate( ncol( samplingcov.mat ) , leverage )*samplingcov.mat
    adjust.cov <- colSums( adjust.covmat )
    
    
    # This condition below takes care of the bivariate fay herriot
    # This is a technical condition that makes sure we can construct
    # the sampling variance covariance matrix without error for
    # the bivariate case.
    
    if( length( adjust.cov ) ==  1 ){
      
      adjust.cov = as.numeric( adjust.cov)
      
    }
    
    
    # constructing the adjusted variance covariance matrix term in the formula
    # P.R.M-1 
    adjustment.mat <- diag( adjust.var )
    adjustment.mat[ lower.tri( adjustment.mat, diag = FALSE ) ] <- adjust.cov
    adjustment.mat <- forceSymmetric(adjustment.mat, uplo = "L")
    
    estimate <- ( residual.transpose%*%t( residual.transpose ) - 
                    adjustment.mat )/( areas -  1 )
    
    if (det(estimate) <= 0){
      return(0.0001*diag(q))
    }else{
      return( estimate )
    }
    
  }
  
  
  
  
  if(missing(design.mat)){
    return(prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat ))
  }
  else{
    projection <- function( mat ){
      
      mat%*%solve( t( mat )%*%mat )%*%t( mat )
      
    }
    
    # total number of areas
    areas <- nrow( directest.mat )
    
    proj <- projection( design.mat )
    
    fit <- lm(directest.mat ~ design.mat)
    residual.transpose <- t( fit$residuals )
    #residual.transpose <- t( directest.mat )%*%( diag( areas ) - proj )
    
    # p_ii's  used to compute the prasad rao estimator
    leverage <- as.vector( diag( proj ) )
    
    # this is the adjustment for the variance covariance matrix in the 
    # formula P.R.M-1 in page 67
    
    adjust.varmat <- replicate( ncol( samplingvar.mat ) , 1 - leverage )*samplingvar.mat
    adjust.var <- colSums( adjust.varmat )
    
    
    adjust.covmat <- replicate( ncol( samplingcov.mat ) , 1 - leverage )*samplingcov.mat
    adjust.cov <- colSums( adjust.covmat )
    
    
    # This condition below takes care of the bivariate fay herriot
    # This is a technical condition that makes sure we can construct
    # the sampling variance covariance matrix without error for
    # the bivariate case.
    
    if( length( adjust.cov ) ==  1 ){
      
      adjust.cov = as.numeric( adjust.cov)
      
    }
    
    
    # constructing the adjusted variance covariance matrix term in the formula
    # P.R.M-1 
    adjustment.mat <- diag( adjust.var )
    adjustment.mat[ upper.tri( adjustment.mat, diag = FALSE ) ] <- adjust.cov
    adjustment.mat[ lower.tri( adjustment.mat, diag = FALSE ) ] <- adjust.cov
    
    estimate <- ( residual.transpose%*%t( residual.transpose ) - 
                    adjustment.mat )/( areas - ncol ( design.mat ) - 1 )
    
    
    
    return(estimate)
  }
  
}
