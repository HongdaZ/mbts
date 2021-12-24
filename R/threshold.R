## Find threshold for tumor
threshold <- function( res, shift ) {
  nc <- ncol( res )
  m <- res[ 1, nc ]
  sigma2 <- res[ 2, nc ]
  m + sqrt( sigma2 ) * shift[ 3 ]
} 