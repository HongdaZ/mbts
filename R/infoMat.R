## Get info list for further segmentation
infoMat <- function( img, label ) {
  .Call( "infoMat", img, label )
}