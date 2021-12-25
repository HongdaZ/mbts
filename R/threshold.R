#' Get Current Threshold
#' 
#' \code{threshold} gets current threshold for bright signals.
#' @author Hongda Zhang
#' @param res Matrix of estimated mean, variance and interaction coefficients.
#' @param shift Current shift for unary potentials.
#' @export
## Find threshold for tumor
threshold <- function( res, shift ) {
  nc <- ncol( res )
  m <- res[ 1, nc ]
  sigma2 <- res[ 2, nc ]
  m + sqrt( sigma2 ) * shift[ 3 ]
} 