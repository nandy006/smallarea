smallareafit.mult <- function( directest.mat,  samplingvar.mat, samplingcov.mat, design.mat ){
  
  require("Rmpfr")
  options( digits = 16)
  
  
  if(missing(design.mat)){
    
    
    n <- nrow( directest.mat )
    design.mat <- matrix(1, n, 1 )
    p <- 1 
    q <- ncol( directest.mat )
    
    ########################################################################################################
    ################################### Prasad Rao Function ################################################
    ########################################################################################################
    
    
    prasadraomult.nocov <- function( directest.mat, samplingvar.mat, samplingcov.mat ){
      # the data input format is for each small area you have a row with columns giving
      # the values of the direct estimators, followed by colmns giving values of the 
      # covariates, followed by variance of each of the small area attributes followed by 
      # covariance terms in natural order i.e. for example with 4 attributes 
      # cov(y1,y2), cov(y1,y3), cov(y1,y4), cov(y2,y3), 
      # cov(y2,y4),cov(y3,y4) 
      
      
      
      # the projection matrix function
      
      
      # total number of areas
      areas <- nrow( directest.mat )
      
      
      
      #residual.transpose <- t( directest.mat )%*%( diag( areas ) - proj )
      
      residual.transpose <- t( directest.mat ) - rowMeans( t( directest.mat ) )
      
      
      # p_ii's  used to compute the prasad rao estimator
      leverage <- rep( 1/areas, areas )
      
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
    
    
    ################################## end of prasad rao function ##########################################
    
    
    
    
    
    
    
    
    require( Matrix)
    
    x <- 1:n
    
    # Creating the list of matrices psi + Di for each small area
    
    # Sigma.list <- lapply(x, function( x ) {
    #  D <- diag( samplingvar.mat[ x, ] )
    #  D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
    # D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
    #  D <- D + prasadraomult( directest.mat, design.mat, samplingvar.mat, samplingcov.mat )
    # return( D )
    #} )
    
    psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
    
    options(digits = 16)
    Sigmainv.list <- lapply(x, function( x ) {
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      #D <- D + prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
      D <- D + psi.hat
      D <- forceSymmetric(D, uplo = "L")
      return( solve( D ) )
    } )
    
    # The big variance covariance matrix for the entire y Vector stacked verically
    # i.e y1
    #     y2
    #     .
    #     .
    #     yn
    #Sigma <- bdiag( Sigma.list )
    inv.Sigma <- bdiag( Sigmainv.list )
    
    # Make the bigger design matrix 
    # [x11....x1p...........................]
    # [x11....x1p...........................]
    # [.....................................]
    # [.....................................]
    # [x11....x1p...........................]
    # [..........x21.....x2p................]
    # [.....................................]
    # [.....................xn1..........xnp]
    #require( data.table )
    designmatexpanded.list <- lapply(x, function( x ) {
      as.matrix( bdiag( replicate( q, t( as.matrix( design.mat[ x, ] ) ) , simplify=F ) ) )
      
    } )
    designmat.expanded <- do.call( rbind, designmatexpanded.list )
    designmat.expanded <- Matrix( designmat.expanded, sparse = TRUE )
    
    #  replicate( q, t( as.matrix( design.mat[ 1, ] ) ) , simplify=F )
    
    
    directest.vec <- as.vector( t( directest.mat ) )
    
    #beta.estimate <- solve( t( designmat.expanded )%*%solve( Sigma )%*%designmat.expanded )%*%( t( 
    #designmat.expanded )%*%solve( Sigma )%*%directest.vec )
    
    #beta.estimate <- solve( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded )%*%( t( 
    #designmat.expanded )%*%inv.Sigma%*%directest.vec )
    
    #correction.mat <- ( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded + 
    #t( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded ) )/2
    
    beta.estimate <- solve( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded, sparse = TRUE )%*%( t( 
      designmat.expanded )%*%directest.vec )
    
    
    
    options(digits = 16)
    # beta.estimate <- chol2inv( chol( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded ) )%*%( t(
    #  designmat.expanded )%*%inv.Sigma%*%directest.vec )
    
    
    B.est <- matrix( beta.estimate, q, p , byrow = TRUE)
    if( det( psi.hat ) <= 0.01 ){
      thetaest.list <- lapply(x, function( x ) {
        D <- diag( samplingvar.mat[ x, ] )
        #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
        D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
        D <- forceSymmetric(D, uplo = "L")
        #psi.hat <- prasadraomult.nocov( directest.mat, 
        #samplingvar.mat, samplingcov.mat )
        B.est%*%design.mat[ x, ]  
      } )
      
      theta.est <- t( do.call( cbind, thetaest.list ) )
      
      
    }  
    else{  thetaest.list <- lapply(x, function( x ) {
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      psi.hat <- prasadraomult.nocov( directest.mat, 
                                      samplingvar.mat, samplingcov.mat )
      solve( ( solve( D ) + solve( psi.hat ) ) )%*%( solve( D )%*%directest.mat[ x, ]  +
                                                       solve( psi.hat )%*%B.est%*%design.mat[ x, ]  )
    } )
    
    theta.est <- t( do.call( cbind, thetaest.list ) )
    }
    
    #psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
    
    g1 <- lapply(x, function( x ){
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      #psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
      g1 <- sum(diag( psi.hat%*%solve( D + psi.hat )%*%D ))
      
      return( g1 )
    } )
    
    g1 <- unlist(g1)
    
    mymat <- lapply(x, function(x){
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      #psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
      #design <- as.vector(design.mat[x,])
      #design <- design%o%design
      #return(kronecker(solve(psi.hat + D), design))
      return(solve(psi.hat + D))
    })
    
    mymat <- solve( Reduce('+',mymat) )
    
    
    
    
    g2 <- lapply(x, function( x ){
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      #psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
      mymat1 <- psi.hat%*%solve( D + psi.hat )%*%D 
      #design <- as.vector(design.mat[x,])
      #mymat2 <- kronecker(solve(psi.hat),t(design))
      #g2 <- sum(diag(t(mymat2)%*%mymat1%*%mymat1%*%mymat2%*%mymat))
      g2 <- sum(diag( solve(psi.hat)%*%mymat1%*%mymat1%*%solve(psi.hat)%*%mymat))
      return( g2 )
    } )
    
    g2 <- unlist(g2)
    
    
    z <- 1:n
    
    fun5 <- function(z){
      
      u <- matrix(NA, q, q)
      
      for(index1 in 1:q){
        for(index2 in 1:q){
          fun4 <- function(x){
            
            #q <- 
            D <- diag( samplingvar.mat[ x, ] )
            #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
            D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
            D <- forceSymmetric(D, uplo = "L")
            #psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
            mat <- psi.hat + D
            
            fun1 <- function(index1,index2,index3,index4){
              mat[index1,index2]*mat[index3,index4] + mat[index3,index2]*mat[index1,index4]
              
            }
            
            
            fun3 <- function(index3){
              index4 <- 1:q
              #index3 <- 1
              fun2 <- function(index4){
                fun1(index1, index2, index3, index4)
              }
              u <- lapply(index4,fun2)
              u <- unlist(u)
              return(u)
            }
            
            index3 <- 1:q
            
            v <- lapply(index3, fun3)
            
            return(matrix(unlist(v), nrow = q,byrow = TRUE))
            
          }
          
          v <- Reduce('+',lapply(x, fun4))
          
          D <- diag( samplingvar.mat[ z, ] )
          #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
          D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ z, ]
          D <- forceSymmetric(D, uplo = "L")
          psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
          
          dummy <- solve(D + psi.hat)%*%v
          
          u[index1,index2] <- sum(diag(dummy))/(nrow(design.mat)- ncol(design.mat))^2
          
        }
      }
      
      D <- diag( samplingvar.mat[ z, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ z, ]
      D <- forceSymmetric(D, uplo = "L")
      #psi.hat <- prasadraomult.nocov( directest.mat, samplingvar.mat, samplingcov.mat )
      
      g3 <- solve(psi.hat + D)%*%D%*%solve(psi.hat + D)%*%u
      g3 <- sum(diag(g3))
      return(g3)
      
    }
    
    g3 <- lapply(z, fun5)
    
    g3<- unlist(g3)
    
    mspe <- g1 + 2*g2 + g3
    
    
    return(list( theta.est = theta.est, B.estimate = B.est, mspe = mspe,psi.hat = psi.hat))
  }
  
  
  else{
    
    n <- nrow( directest.mat )
    q <- ncol( directest.mat )
    p <- ncol(design.mat)
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
    
    require( Matrix)
    
    x <- 1:n
    
    # Creating the list of matrices psi + Di for each small area
    
    # Sigma.list <- lapply(x, function( x ) {
    #  D <- diag( samplingvar.mat[ x, ] )
    #  D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
    # D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
    #  D <- D + prasadraomult( directest.mat, design.mat, samplingvar.mat, samplingcov.mat )
    # return( D )
    #} )
    
    psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
    
    options(digits = 16)
    Sigmainv.list <- lapply(x, function( x ) {
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      D <- D + psi.hat
      D <- forceSymmetric(D, uplo = "L")
      return( solve( D ) )
    } )
    
    # The big variance covariance matrix for the entire y Vector stacked verically
    # i.e y1
    #     y2
    #     .
    #     .
    #     yn
    #Sigma <- bdiag( Sigma.list )
    inv.Sigma <- bdiag( Sigmainv.list )
    
    # Make the bigger design matrix 
    # [x11....x1p...........................]
    # [x11....x1p...........................]
    # [.....................................]
    # [.....................................]
    # [x11....x1p...........................]
    # [..........x21.....x2p................]
    # [.....................................]
    # [.....................xn1..........xnp]
    #require( data.table )
    designmatexpanded.list <- lapply(x, function( x ) {
      as.matrix( bdiag( replicate( q, t( as.matrix( design.mat[ x, ] ) ) , simplify=F ) ) )
      
    } )
    designmat.expanded <- do.call( rbind, designmatexpanded.list )
    designmat.expanded <- Matrix( designmat.expanded, sparse = TRUE )
    
    #  replicate( q, t( as.matrix( design.mat[ 1, ] ) ) , simplify=F )
    
    
    directest.vec <- as.vector( t( directest.mat ) )
    
    #beta.estimate <- solve( t( designmat.expanded )%*%solve( Sigma )%*%designmat.expanded )%*%( t( 
    #designmat.expanded )%*%solve( Sigma )%*%directest.vec )
    
    #beta.estimate <- solve( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded )%*%( t( 
    #designmat.expanded )%*%inv.Sigma%*%directest.vec )
    
    #correction.mat <- ( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded + 
    #t( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded ) )/2
    
    beta.estimate <- solve( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded, sparse = TRUE )%*%( t( 
      designmat.expanded )%*%inv.Sigma%*%directest.vec )
    
    
    
    options(digits = 16)
    # beta.estimate <- chol2inv( chol( t( designmat.expanded )%*%inv.Sigma%*%designmat.expanded ) )%*%( t(
    #  designmat.expanded )%*%inv.Sigma%*%directest.vec )
    
    
    B.est <- matrix( beta.estimate, q, p , byrow = TRUE)
    if( det( psi.hat ) <= 0.01 ){
      thetaest.list <- lapply(x, function( x ) {
        D <- diag( samplingvar.mat[ x, ] )
        #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
        D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
        D <- forceSymmetric(D, uplo = "L")
        psi.hat <- prasadraomult( directest.mat, 
                                  samplingvar.mat, samplingcov.mat, design.mat )
        B.est%*%design.mat[ x, ]  
      } )
      
      theta.est <- t( do.call( cbind, thetaest.list ) )
      
      
    }  
    else{  thetaest.list <- lapply(x, function( x ) {
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      psi.hat <- prasadraomult( directest.mat, 
                                samplingvar.mat, samplingcov.mat, design.mat )
      psi.hat%*%solve( D + psi.hat )%*%D%*%( solve( D )%*%directest.mat[ x, ]  +
                                               solve( psi.hat )%*%B.est%*%design.mat[ x, ]  )
    } )
    
    theta.est <- t( do.call( cbind, thetaest.list ) )
    }
    
    #psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
    
    
    g1 <- lapply(x, function( x ){
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      #psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
      g1 <- sum(diag( psi.hat%*%solve( D + psi.hat )%*%D ))
      
      return( g1 )
    } )
    
    g1 <- unlist(g1)
    
    mymat <- lapply(x, function(x){
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      #psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
      design <- as.vector(design.mat[x,])
      design <- design%o%design
      return(kronecker(solve(psi.hat + D), design))
    })
    
    mymat <- solve( Reduce('+',mymat) )
    
    
    
    
    g2 <- lapply(x, function( x ){
      D <- diag( samplingvar.mat[ x, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D <- forceSymmetric(D, uplo = "L")
      #psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
      mymat1 <- psi.hat%*%solve( D + psi.hat )%*%D 
      design <- as.vector(design.mat[x,])
      mymat2 <- kronecker(solve(psi.hat),t(design))
      g2 <- sum(diag(t(mymat2)%*%mymat1%*%mymat1%*%mymat2%*%mymat))
      return( g2 )
    } )
    
    g2 <- unlist(g2)
    
    
    z <- 1:n
    
    fun5 <- function(z){
      
      u <- matrix(NA, q, q)
      
      for(index1 in 1:q){
        for(index2 in 1:q){
          fun4 <- function(x){
            
            #q <- 
            D <- diag( samplingvar.mat[ x, ] )
            #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
            D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
            D <- forceSymmetric(D, uplo = "L")
            #psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
            mat <- psi.hat + D
            
            fun1 <- function(index1,index2,index3,index4){
              mat[index1,index2]*mat[index3,index4] + mat[index3,index2]*mat[index1,index4]
              
            }
            
            
            fun3 <- function(index3){
              index4 <- 1:q
              #index3 <- 1
              fun2 <- function(index4){
                fun1(index1, index2, index3, index4)
              }
              u <- lapply(index4,fun2)
              u <- unlist(u)
              return(u)
            }
            
            index3 <- 1:q
            
            v <- lapply(index3, fun3)
            
            return(matrix(unlist(v), nrow = q,byrow = TRUE))
            
          }
          
          v <- Reduce('+',lapply(x, fun4))
          
          D <- diag( samplingvar.mat[ z, ] )
          #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
          D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ z, ]
          D <- forceSymmetric(D, uplo = "L")
          #psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
          
          dummy <- solve(D + psi.hat)%*%v
          
          u[index1,index2] <- sum(diag(dummy))/(nrow(design.mat)- ncol(design.mat))^2
          
        }
      }
      
      D <- diag( samplingvar.mat[ z, ] )
      #D[ upper.tri( D, diag = FALSE ) ] <- samplingcov.mat[ x, ]
      D[ lower.tri( D, diag = FALSE ) ] <- samplingcov.mat[ z, ]
      D <- forceSymmetric(D, uplo = "L")
      psi.hat <- prasadraomult( directest.mat, samplingvar.mat, samplingcov.mat, design.mat )
      
      g3 <- solve(psi.hat + D)%*%D%*%solve(psi.hat + D)%*%u
      g3 <- sum(diag(g3))
      return(g3)
      
    }
    
    g3 <- lapply(z, fun5)
    
    g3<- unlist(g3)
    
    mspe <- g1 + 2*g2 + g3
    
    
    return(list( theta.est = theta.est, B.estimate = B.est, mspe = mspe, psi.hat = psi.hat))
    
  }
  
  
  
}





