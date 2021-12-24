# Initialize data for estimation
initEst <- function( label, intst ) {
  info <- indexMat( intst, label )
  seg <- label[ info$idx ]
  pad <- vector( mode = "integer", length( seg ) )
  seg_pad <- rbind( seg, pad )
  
  res <- list( info = info, seg = seg_pad, dim = dim( label ) )
  return( res )
} 