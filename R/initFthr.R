## Initialize data for further segmentation of edema
initFther <- function( label, intst ) {
  info <- infoMat( intst, label )
  seg <- label[ info$idx ]
  pad <- vector( mode = "integer", length( seg ) )
  seg_pad <- rbind( seg, pad )
  
  res <- list( info = info, seg = seg_pad )
  return( res )
}