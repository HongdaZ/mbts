## Get shift for tumor
## Find threshold for tumor
getShift <- function( res, m ) {
  nc <- ncol( res )
  m1 <- res[ 1, nc ]
  sigma2 <- res[ 2, nc ]
  ( m - m1 ) / sqrt( sigma2 ) 
} 