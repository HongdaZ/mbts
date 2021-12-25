#' Get Shift for Unary Potentials
#' 
#' \code{getShift} get the shift for unary potentials.
#' @author Hongda Zhang
#' @param res Matrix of estimated mean, variance and interaction coefficients.
#' @param m Threshold for bright signals.
#' @export
## Get shift for tumor
## Find threshold for tumor
getShift <- function( res, m ) {
  nc <- ncol( res )
  m1 <- res[ 1, nc ]
  sigma2 <- res[ 2, nc ]
  ( m - m1 ) / sqrt( sigma2 ) 
} 